#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import time
import shutil
import argparse
import subprocess as sp
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

# ---------------------------
# 你的两个函数（保持签名，内部改为在脚本目录下执行）
# ---------------------------

def submit_sbatch_script(script_path):
    """Submit a job with sbatch and return the job ID.
    为了让 slurm-<jobid>.out 落在脚本所在目录，这里在脚本目录下执行 sbatch。
    """
    script_path = os.path.abspath(script_path)
    cwd = os.path.dirname(script_path)

    # 提交脚本并捕获输出（工作目录=脚本所在目录）
    result = sp.run(["sbatch", os.path.basename(script_path)],
                    stdout=sp.PIPE, stderr=sp.PIPE, text=True, cwd=cwd)

    if result.returncode != 0:
        print(f"[{script_path}] submit unsuccessfully:", result.stderr.strip())
        return None

    # sbatch 输出通常为: "Submitted batch job <job_id>"
    out = result.stdout.strip()
    toks = out.split()
    job_id = toks[-1] if toks else None
    print(f"[{script_path}] submit successfully, job ID: {job_id}")
    return job_id


def monitor_job(job_id):
    """Monitor the status of a job using squeue and display runtime."""
    start_time = time.time()  # 记录作业提交的开始时间

    # 用更稳定的格式输出，只取 JOBID 和 状态（%i %t）
    while True:
        result = sp.run(["squeue", "--noheader", "-o", "%i %t", "-j", str(job_id)],
                        stdout=sp.PIPE, stderr=sp.PIPE, text=True)

        line = result.stdout.strip()
        if result.returncode != 0:
            # squeue 异常（例如短暂的调度器抖动），稍后重试
            print(f"[{job_id}] squeue error: {result.stderr.strip()}")
            time.sleep(5)
            continue

        if not line:
            # squeue 不再显示该作业：说明作业已退出（成功或失败）
            print(f"job {job_id} completed.")
            break

        # 解析状态（PD/R/RU/CG 等）
        try:
            jid, st = line.split()[0], line.split()[1]
        except Exception:
            jid, st = job_id, "UNK"

        elapsed_time = time.time() - start_time
        elapsed_minutes = int(elapsed_time // 60)
        elapsed_seconds = int(elapsed_time % 60)

        print(f"job {jid} is running... elapsed: {elapsed_minutes}分 {elapsed_seconds}秒  status: {st}")
        time.sleep(5)

# ---------------------------
# 帮助函数
# ---------------------------

def find_subdirs(work_dir):
    """返回工作目录下一层的所有子目录（按名称排序）"""
    p = Path(work_dir)
    return sorted([str(d) for d in p.iterdir() if d.is_dir()])

def copy_hm_sh(src_hm_sh, dst_dir):
    """把 HM.sh 拷贝到 dst_dir。若已有则覆盖。返回目标脚本路径。"""
    dst = Path(dst_dir) / "HM.sh"
    shutil.copy2(src_hm_sh, dst)
    # 确保可执行
    try:
        os.chmod(dst, 0o755)
    except Exception:
        pass
    return str(dst)

def run_gnuplot_if_exists(work_dir, plot_file="Nodes_dis.gnu"):
    """如果工作目录下存在 plot_file，则调用 gnuplot 绘图"""
    plot_path = Path(work_dir) / plot_file
    if not plot_path.exists():
        print(f"[{work_dir}] {plot_file} not found. skip gnuplot.")
        return True

    print(f"[{work_dir}] running: gnuplot {plot_file}")
    proc = sp.run(["gnuplot", plot_file], cwd=str(work_dir),
                  stdout=sp.PIPE, stderr=sp.PIPE, text=True)
    if proc.returncode != 0:
        print(f"[{work_dir}] gnuplot error:\n{proc.stderr}")
        return False
    print(f"[{work_dir}] gnuplot done.")
    return True

def job_flow_for_subdir(subdir, hm_sh_src, dry_run=False):
    """
    进入 subdir：
      1) 拷贝 HM.sh
      2) sbatch 提交并监控直到结束（在 subdir 内执行，保证 slurm 输出落在 subdir）
      3) 结束后运行 gnuplot Nodes_dis.gnu
    """
    try:
        print(f"=== [{subdir}] start ===")

        # 1) 拷贝 HM.sh
        script_path = copy_hm_sh(hm_sh_src, subdir)
        print(f"[{subdir}] HM.sh copied to {script_path}")

        if dry_run:
            print(f"[{subdir}] dry-run: skip submit & monitor.")
            run_gnuplot_if_exists(subdir)
            print(f"=== [{subdir}] done (dry-run) ===")
            return True

        # 2) 提交作业并监控（submit 内部在脚本目录执行 sbatch）
        job_id = submit_sbatch_script(script_path)
        if job_id is None:
            print(f"[{subdir}] submit failed. skip gnuplot.")
            return False

        monitor_job(job_id)

        # 3) gnuplot（在子目录内执行）
        ok = run_gnuplot_if_exists(subdir)
        print(f"=== [{subdir}] done ===")
        return ok

    except KeyboardInterrupt:
        print(f"[{subdir}] interrupted by user.")
        return False
    except Exception as e:
        print(f"[{subdir}] error: {e}")
        return False

# ---------------------------
# 主逻辑：并发提交/监控
# ---------------------------

def main():
    ap = argparse.ArgumentParser(
        description="批量并发提交 SLURM 任务：每个子目录复制 HM.sh → sbatch 提交 → 监控 → gnuplot")
    ap.add_argument("work_dir", help="顶层工作目录（将扫描其下一层子目录）")
    ap.add_argument("--hm-sh", required=True, help="HM.sh 模板的绝对或相对路径")
    ap.add_argument("--max-parallel", type=int, default=5, help="并发上限（默认 4）")
    ap.add_argument("--include", nargs="*", default=None,
                    help="只处理这些子目录名（空格分隔）。默认处理全部子目录")
    ap.add_argument("--exclude", nargs="*", default=None,
                    help="排除这些子目录名（空格分隔）")
    ap.add_argument("--dry-run", action="store_true", help="不提交，仅演示拷贝与流程")
    args = ap.parse_args()

    work_dir = os.path.abspath(args.work_dir)
    hm_sh_src = os.path.abspath(args.hm_sh)

    if not os.path.isdir(work_dir):
        print(f"work_dir 不存在或不是目录: {work_dir}")
        return 1
    if not os.path.isfile(hm_sh_src):
        print(f"HM.sh 模板不存在: {hm_sh_src}")
        return 1

    all_subdirs = find_subdirs(work_dir)

    # include / exclude 过滤
    if args.include:
        inc = set(args.include)
        all_subdirs = [d for d in all_subdirs if Path(d).name in inc]
    if args.exclude:
        exc = set(args.exclude)
        all_subdirs = [d for d in all_subdirs if Path(d).name not in exc]

    if not all_subdirs:
        print("没有待处理的子目录。")
        return 0

    print(f"将处理 {len(all_subdirs)} 个子目录；并发上限={args.max_parallel}")

    # 线程池做并发
    results = {}
    with ThreadPoolExecutor(max_workers=max(1, args.max_parallel)) as ex:
        fut2dir = {
            ex.submit(job_flow_for_subdir, d, hm_sh_src, args.dry_run): d
            for d in all_subdirs
        }
        for fut in as_completed(fut2dir):
            d = fut2dir[fut]
            ok = False
            try:
                ok = fut.result()
            except Exception as e:
                print(f"[{d}] raised exception: {e}")
                ok = False
            results[d] = ok

    # 汇总
    success = [Path(d).name for d, ok in results.items() if ok]
    failed  = [Path(d).name for d, ok in results.items() if not ok]
    print("\n=== 汇总 ===")
    print("成功：", success)
    print("失败：", failed)

    return 0 if len(failed) == 0 else 2


if __name__ == "__main__":
    raise SystemExit(main())