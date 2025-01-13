#先计算左旋光
#cp wt.in_L wt.in
job_id1=$(sbatch W90.sh | awk '{print $4}')

echo $job_id1
# 等待第一个任务完成

while squeue -j $job_id1 | grep -q "^ *$job_id1 "; do
    sleep 10
done
# 第一个任务完成后执行一系列操作
# 这里可以是任意你需要的操作，比如数据处理、文件移动等
echo "第一个任务已完成，开始执行后续操作..."
# (在此处添加你的操作代码)
# 移动dat文件

#cp spectrum_unfold_kpath.dat L.dat
#再计算右旋光
#cp wt.in_R wt.in

# 提交第二个任务，并获取任务ID
#job_id2=$(sbatch NANO.sh | awk '{print $4}')
#
#echo $job_id2
#
#while squeue -j $job_id2 | grep -q "^ *$job_id2 "; do
#    sleep 10
#done
#
#echo "第二个任务已完成，开始执行后续操作..."
#
#cp spectrum_unfold_kpath.dat R.dat
# 绘图
gnuplot w9.gnu
