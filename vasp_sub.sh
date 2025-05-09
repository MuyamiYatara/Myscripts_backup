#!/bin/bash

# 定义一些变量
hrdat="INCAR"
kpoints="KPOINTS"
potcar="POTCAR"
input_file="POSCAR"
job_script="NANO.sh"
output_dir_base="van"  # 输出文件夹的基础名称
param_line_number=0  # 需要修改的行号
initial_value=0  # 初始值
increment=1  # 每次递增的值
total_jobs=27  # 总共需要提交的任务数
max_concurrent_jobs=9  # 最大并发任务数

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
    param_line_number=$((i + 89))
    echo "Deleting line: $param_line_number"
    sed "${param_line_number}d" "$input_file" > "${output_dir}/${input_file}"
    
    # 拷贝作业脚本到输出目录
    cp "$job_script" "$output_dir/"
    cp "$hrdat" "$output_dir/"
    cp "$kpoints" "$output_dir/"
    cp "$potcar" "$output_dir/"
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

