#!/bin/bash

# 定义一些变量
output_dir_base="theta"  # 输出文件夹的基础名称
total_jobs=36  # 总共需要处理的任务数
increment=10
initial_value=0
target_line=10  # 需要提取的行号
output_file="combined_output.txt"  # 输出文件
every_output_file="rho_total_mu_0.00000eV.dat"

# 清空或创建输出文件
> "$output_file"

# 遍历所有输出目录
for i in $(seq 1 $total_jobs); do
    # 生成当前输出目录名称
    new_value=$(echo "$initial_value + ($i - 1) * $increment" | bc)
    output_dir="${output_dir_base}_${new_value}"
    
    # 检查rho文件是否存在
    if [ -f "${output_dir}/${every_output_file}" ]; then
        # 提取文件中的指定行
        line_content=$(sed -n "${target_line}p" "${output_dir}/${every_output_file}")
        
        # 从wt.in中获取提交任务时修改的参数值
        input_value=${new_value}
        
        # 将参数值和提取的行内容写入输出文件
        echo "$input_value $line_content" >> "$output_file"
    else
        echo "Warning: ${output_dir}/${every_output_file} not found"
    fi
done

echo "All data combined into $output_file"
