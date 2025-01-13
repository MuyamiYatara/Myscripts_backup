#!/bin/bash

# 定义一些变量
output_dir_base="theta"  # 输出文件夹的基础名称
total_jobs=31 # 总共的温度数目
increment=9
initial_value=30
block_line_number=10 #每一块有多少行
first_target_line=4  # 需要提取的第一行行号
output_file="Tdata_output.txt"  # 输出文件
every_output_file="seebeck_total_mu_0.000eV.dat"

# 清空或创建输出文件
> "$output_file"

# 遍历所有输出目录
for i in $(seq 1 $total_jobs); do
    # 生成当前行（温度）的数值
    new_value=$(echo "$initial_value + ($i - 1) * $increment" | bc)
    
    # 生成需要读取的行
    target_line=$(echo "${first_target_line}+($i-1)*(${block_line_number}+3)" | bc)

    # 提取文件中的指定行
    line_content=$(sed -n "${target_line}p" "${every_output_file}")
                
    # 将参数值和提取的行内容写入输出文件
    echo "$new_value $line_content" >> "$output_file"

done

echo "All data combined into $output_file"
