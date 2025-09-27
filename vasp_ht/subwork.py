import os
import sys
import subprocess
import threading
import time
from queue import Queue

# ================= 配置 =================
work_root = "/data/home/ycshen/openMX_wannier/multi_test/work_soc"
max_concurrent = 5
refresh_interval = 10
progress_log = os.path.join(work_root, "progress.log")
failed_log = os.path.join(work_root, "failed_tasks.log")

STEPS = ["SCF", "BD", "AUTOWR", "WRSCF", "WR"]
STEP_MARKERS = [
    "@@@ SCF completed @@@",
    "@@@ BD completed @@@",
    "@@@ AUTOWR completed @@@",
    "@@@ WRSCF completed @@@",
    "@@@ WR completed @@@"
]

lock = threading.Lock()
active_progress = {}  # 材料名 -> 当前完成步骤索引
task_queue = Queue()

# ================= 日志刷新线程 =================
def log_progress(stop_event):
    while not stop_event.is_set():
        with lock:
            lines = []
            for mat, step_idx in active_progress.items():
                bar = "█" * step_idx + "-" * (len(STEPS) - step_idx)
                lines.append(f"{mat}: [{bar}] ({step_idx}/{len(STEPS)})")
            with open(progress_log, "w") as f:
                f.write("\n".join(lines) + "\n")
        time.sleep(refresh_interval)

# ================= 材料任务 =================
def run_material(mat_dir):
    mat_name = os.path.basename(mat_dir)
    log_file = os.path.join(mat_dir, "auto_log")
    step_idx = 0

    with lock:
        active_progress[mat_name] = step_idx

    script_path = os.path.join(mat_dir, "semiauto_dft_wannier.py")
    if not os.path.isfile(script_path):
        print(f"[{mat_name}] 没有 semiauto_dft_wannier.py，跳过")
        with lock:
            active_progress.pop(mat_name, None)
        return

    with open(log_file, "w") as f:
        process = None
        try:
            process = subprocess.Popen(
                [sys.executable, "-u", "semiauto_dft_wannier.py"],
                cwd=mat_dir,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                bufsize=1,
                universal_newlines=True
            )

            last_refresh = time.time()
            while True:
                line = process.stdout.readline()
                if line == "" and process.poll() is not None:
                    break
                if line:
                    f.write(line)
                    f.flush()
                    print(f"[{mat_name}] {line.strip()}")  # 输出到控制台，便于排查

                    while step_idx < len(STEP_MARKERS) and STEP_MARKERS[step_idx] in line:
                        step_idx += 1
                        with lock:
                            active_progress[mat_name] = step_idx

                if time.time() - last_refresh >= refresh_interval:
                    with lock:
                        active_progress[mat_name] = step_idx
                    last_refresh = time.time()

        except Exception as e:
            print(f"[{mat_name}] 执行出错：{e}")
            f.write(f"\n[ERROR] 执行异常：{e}\n")

        finally:
            try:
                if process:
                    if process.poll() is None:
                        process.terminate()
                        time.sleep(2)
                    process.wait()
            except Exception as e:
                print(f"[{mat_name}] 子进程清理失败：{e}")

            # 记录失败任务
            if process and process.returncode != 0:
                with open(failed_log, "a") as ff:
                    ff.write(f"{mat_name} failed with code {process.returncode}\n")

            with lock:
                active_progress.pop(mat_name, None)

# ================= 工作者线程 =================
def worker():
    while True:
        mat_dir = task_queue.get()
        if mat_dir is None:
            break
        try:
            run_material(mat_dir)
        except Exception as e:
            print(f"[ERROR] 材料任务失败：{mat_dir}, 原因：{e}")
        task_queue.task_done()

# ================= 主程序 =================
def main():
    materials_dirs = [os.path.join(work_root, d) for d in os.listdir(work_root)
                      if os.path.isdir(os.path.join(work_root, d))]

    stop_event = threading.Event()
    refresher_thread = threading.Thread(target=log_progress, args=(stop_event,))
    refresher_thread.start()

    # 启动 worker 线程
    threads = []
    for _ in range(max_concurrent):
        t = threading.Thread(target=worker)
        t.start()
        threads.append(t)

    # 填充任务队列
    for mat_dir in materials_dirs:
        task_queue.put(mat_dir)

    task_queue.join()  # 等待所有任务完成

    # 停止 worker
    for _ in range(max_concurrent):
        task_queue.put(None)
    for t in threads:
        t.join()

    stop_event.set()
    refresher_thread.join()
    print("所有材料任务完成。")

if __name__ == "__main__":
    main()