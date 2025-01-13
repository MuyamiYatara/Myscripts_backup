#!/bin/bash

# 定义一些变量
hrdat="Bismuth_hr.dat"
input_file="wt.in"
job_script="NANO.sh"
output_dir_base="theta"  # 输出文件夹的基础名称
param_line_number=15  # 需要修改的行号
initial_value=0  # 初始值
increment=10  # 每次递增的值
total_jobs=36  # 总共需要提交的任务数
max_concurrent_jobs=4  # 最大并发任务数

#清除上次计算遗留的文件和文件夹
#./clear.sh
# 创建批量任务提交
for i in $(seq 1 $total_jobs); do
    # 检查当前正在运行的任务数量，如果达到最大并发任务数，则等待
    while [ $(squeue -u $USER | grep $USER | wc -l) -ge $max_concurrent_jobs ]; do
        sleep 10  # 等待10秒后再次检查
    done
    
    # 创建输出目录
    new_value=$(echo "$initial_value + ($i - 1) * $increment" | bc)
    output_dir="${output_dir_base}_${new_value}"
    mkdir -p "$output_dir"
    
    # 修改输入文件
    sed "${param_line_number}s/.*/Btheta= ${new_value}, Bphi= 90/" "$input_file" > "${output_dir}/${input_file}"
    
    # 拷贝作业脚本到输出目录
    cp "$job_script" "$output_dir/"
    cp "$hrdat" "$output_dir/"
    
    # 提交作业
    cd "$output_dir"
    sbatch "$job_script"
    cd ..
    
    echo "Submitted job $i with parameter $new_value, waiting for available slot..."
done

echo "All jobs submitted!"

while [ $(squeue -u $USER | grep $USER | wc -l) -ge 1 ]; do
        sleep 10  # 等待10秒后再次检查
done

echo "All jobs completed!"

