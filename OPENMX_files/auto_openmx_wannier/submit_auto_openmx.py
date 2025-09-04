#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import subprocess as sp
import time
import re
from pathlib import Path
import shutil
import sys, os

# ======== 必改配置 ========
INPUT_DIR = Path("/home/users/shenyc/openmx_wannier/test_POSCARS/")   # 需要读取文件名的目录
WORKFILE_DIR = Path("/home/users/shenyc/Myscripts/OPENMX_files/auto_openmx_wannier/")       # 模板 HM.sh 与 auto_openmx.py 路径
WORK_DIR = Path("./res")                 # 生成临时 HM 脚本的目录
FILE_SUFFIX = ".POSCAR"                           # 可选：只匹配某后缀，如 ".vasp"；留空=不过滤
POLL_INTERVAL = 10                         # 监控轮询间隔（秒）
MAX_INFLIGHT = 4
# =========================




HM_TEMPLATE = WORKFILE_DIR / "HM.sh"
AUTO_OPENMX_PY =  WORKFILE_DIR / "auto_openmx.py"

# 匹配并替换 HM.sh 中的目标行：python auto_openmx.py xxxxx
# - 容忍前后空格、python/python3、额外参数
RE_LINE = re.compile(
    r"^(?P<prefix>\s*python(?:3)?\s+\./auto_openmx\.py\s+)(?P<arg>\S+)(?P<suffix>.*)$"
)

FINAL_STATES = {
    "COMPLETED", "FAILED", "CANCELLED", "TIMEOUT", "OUT_OF_MEMORY",
    "PREEMPTED", "BOOT_FAIL", "NODE_FAIL"
}

def run(cmd, check=True, text=True):
    return sp.run(cmd, check=check, text=text, stdout=sp.PIPE, stderr=sp.PIPE)

def list_files(input_dir: Path, suffix: str):
    files = []
    for p in sorted(input_dir.iterdir()):
        if p.is_file():
            if suffix:
                if p.name.endswith(suffix):
                    files.append(p)
            else:
                files.append(p)
    return files

def make_job_script(template_path: Path, outfile: Path, filename_to_fill: str):
    content = template_path.read_text(encoding="utf-8", errors="ignore").splitlines()
    new_lines = []
    replaced = False

    # 提取编号部分，例如 mp-12345.vasp -> 12345
    m = re.search(r'^([0-9]+(?:\.[0-9]+)?)_', filename_to_fill)
    if m:
        number = m.group(1)
    else:
        number = filename_to_fill   # 如果没找到数字，就用完整文件名

    for line in content:
        m = RE_LINE.match(line)
        if m and not replaced:
            # 在替换 auto_openmx.py 之前，插入 mkdir
            new_lines.append(f"mkdir mat-{number}")
            # 替换 xxxxx 为实际文件名
            newline = f"{m.group('prefix')}{filename_to_fill}{m.group('suffix')}"
            new_lines.append(newline)
            replaced = True
        else:
            new_lines.append(line)

    if not replaced:
        raise RuntimeError(
            f"未在 {template_path} 中找到形如 `python auto_openmx.py xxxxx` 的行，请检查模板。"
        )

    outfile.write_text("\n".join(new_lines) + "\n", encoding="utf-8")
    # 拷贝可执行权限
    st = template_path.stat()
    outfile.chmod(st.st_mode)

def sbatch_submit(job_script: Path):
    # 调用 sbatch 提交，解析 jobid
    res = run(["sbatch", str(job_script)])
    out = res.stdout.strip()
    # 典型输出：Submitted batch job 123456
    m = re.search(r"Submitted batch job (\d+)", out)
    if not m:
        raise RuntimeError(f"无法解析 sbatch 输出：{out}\nSTDERR: {res.stderr}")
    jobid = m.group(1)
    print(f"[提交] {job_script.name} -> JobID {jobid}")
    return jobid

def query_squeue_state(jobid: str):
    # 返回当前 squeue 状态（如 R/PD/CG...）；不在队列则返回 None
    res = run(["squeue", "-h", "-j", jobid, "-o", "%T"], check=False)
    txt = res.stdout.strip()
    if txt == "":
        return None
    return txt.splitlines()[0].strip()

def query_squeue_states(jobids):
    """批量获取 squeue 状态；返回 {jobid: StateText}，不在队列的 jobid 不包含在结果里。"""
    if not jobids:
        return {}
    try:
        res = run(["squeue", "-h", "-j", ",".join(jobids), "-o", "%i|%T"], check=False)
        states = {}
        for line in res.stdout.splitlines():
            line = line.strip()
            if not line:
                continue
            jid, st = line.split("|", 1)
            states[jid.strip()] = st.strip()
        return states
    except Exception:
        return {}

def query_sacct_state(jobid: str):
    # 通过 sacct 获取最终状态（对已结束作业更可靠）
    # 使用 -P -n 便于解析；State 字段可能含子状态（例如 COMPLETED）
    res = run(["sacct", "-j", jobid, "--format=State", "-n", "-P"], check=False)
    lines = [x.strip() for x in res.stdout.splitlines() if x.strip()]
    if not lines:
        return None
    # 通常第一行是主作业状态，后续是 step；取第一行
    state = lines[0].split("|")[0].strip()
    # 去掉可能的后缀（如 COMPLETED,COMPLETED）
    state = state.split()[0].strip()
    return state

def wait_until_done(jobid: str, poll: int = 30):
    last_print = ""
    while True:
        sq = query_squeue_state(jobid)
        if sq is not None:
            msg = f"[监控] Job {jobid} 当前状态（squeue）: {sq}"
            if msg != last_print:
                print(msg)
                last_print = msg
            time.sleep(poll)
            continue

        # 不在 squeue，转 sacct 查最终状态
        st = query_sacct_state(jobid)
        if st:
            print(f"[完成] Job {jobid} 最终状态（sacct）: {st}")
            return st
        else:
            # sacct 可能有轻微延迟，稍等再查
            print(f"[监控] Job {jobid} 不在队列中，等待 sacct 更新…")
            time.sleep(poll)

def main():
    if not HM_TEMPLATE.exists():
        print(f"模板不存在: {HM_TEMPLATE}", file=sys.stderr)
        sys.exit(1)
    WORK_DIR.mkdir(parents=True, exist_ok=True)

    files = list_files(INPUT_DIR, FILE_SUFFIX)
    if not files:
        print(f"在 {INPUT_DIR} 未找到匹配文件（后缀过滤：{FILE_SUFFIX or '无'}）")
        sys.exit(0)

    print(f"共发现 {len(files)} 个文件，将依次提交：")
    for i, f in enumerate(files, 1):
        print(f"  [{i:02d}] {f.name}")

    # 先为所有文件生成脚本（与你原来一样）
    job_scripts = []
    for f in files:
        job_script = WORK_DIR / f"HM_{f.stem}.sh"
        try:
            make_job_script(HM_TEMPLATE, job_script, f.name+" "+str(INPUT_DIR)+"/")  # 或换成 str(f.resolve()) 看你原来写法
            job_scripts.append((f, Path(f"HM_{f.stem}.sh")))
        except Exception as e:
            print(f"[跳过] 生成 {job_script.name} 失败：{e}")
    
    os.chdir(WORK_DIR)
    sp.run(["cp", str(AUTO_OPENMX_PY), "./"])


    next_index = 0
    active = {}        # {jobid: (file_path, job_script)}
    finished = []      # [(jobid, file_path, job_script, final_state)]

    while next_index < len(job_scripts) or active:
        # 补齐窗口
        while next_index < len(job_scripts) and len(active) < MAX_INFLIGHT:
            f, js = job_scripts[next_index]
            try:
                # 与你原来一致：cwd 设为 INPUT_DIR，便于相对路径
                res = sp.run(["sbatch", str(js)],
                            text=True, stdout=sp.PIPE, stderr=sp.PIPE, check=True)
                out = res.stdout.strip()
                m = re.search(r"Submitted batch job (\d+)", out)
                if not m:
                    raise RuntimeError(f"无法解析 sbatch 输出：{out}\nSTDERR: {res.stderr}")
                jobid = m.group(1)
                active[jobid] = (f, js)
                print(f"[提交] {js.name} (文件: {f.name}) -> JobID {jobid}  ｜ 当前并发: {len(active)}/{MAX_INFLIGHT}")
            except sp.CalledProcessError as e:
                print(f"[错误] 提交 {js.name} 失败：{e.stdout}\n{e.stderr}")
            except Exception as e:
                print(f"[错误] {js.name} 发生异常：{e}")
            finally:
                next_index += 1

        if not active:
            break

        # 轮询一次所有激活作业
        time.sleep(POLL_INTERVAL)
        jobids = list(active.keys())
        sq_states = query_squeue_states(jobids)

        to_remove = []
        for jid in jobids:
            if jid in sq_states:
                # 仍在队列/运行，打印一次状态即可
                print(f"[监控] Job {jid} （squeue）: {sq_states[jid]}")
                continue

            # 不在 squeue，用 sacct 查最终态
            st = query_sacct_state(jid)
            if st:
                f, js = active[jid]
                print(f"[完成] Job {jid} 最终状态（sacct）: {st}  <- {js.name}（文件: {f.name}）")
                finished.append((jid, f, js, st))
                to_remove.append(jid)
            else:
                print(f"[监控] Job {jid} 不在队列中，等待 sacct 更新…")

        for jid in to_remove:
            active.pop(jid, None)

    # 汇总
    ok = sum(1 for _,_,_,st in finished if st == "COMPLETED")
    print("\n========== 汇总 ==========")
    print(f"总计：{len(files)}，已提交：{len(finished)}，成功：{ok}，失败/其他：{len(finished)-ok}")
    for jid, f, js, st in finished:
        tag = "OK" if st == "COMPLETED" else "!!"
        print(f"[{tag}] {jid}  {st:>12}  {js.name}  ({f.name})")

if __name__ == "__main__":
    main()