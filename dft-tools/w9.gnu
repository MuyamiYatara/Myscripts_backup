set terminal pngcairo enhanced color font ",60" size 1920,1080
set output 'spectrum_unfold_kpath.png'
set style data dots
#set nokey
#set xrange [0: 6]
set xtics font ",30"
set ytics font ",30"
set xlabel font ",35"
set ylabel font ",35"
upper_bound = 2
lower_bound = -2

set yrange [lower_bound : upper_bound]
### arrows and xsticks copied from wannier90_band.gnu file ###
set arrow from  0.76509,  lower_bound to  0.76509,  upper_bound nohead
set arrow from  1.53018,  lower_bound to  1.53018,  upper_bound nohead
set arrow from  2.61218,  lower_bound to  2.61218,  upper_bound nohead
set arrow from  3.69419,  lower_bound to  3.69419,  upper_bound nohead
set arrow from  4.45928,  lower_bound to  4.45928,  upper_bound nohead
set arrow from  5.22437,  lower_bound to  5.22437,  upper_bound nohead
set xtics (" � "  0.00000," X "  0.76509," M "  1.53018," � "  2.61218," Z "  3.69419," R "  4.45928," A "  5.22437," Z "  6.30637)

set arrow from graph 0, first 0 to graph 1, first 0 nohead lw 4 lc rgb "green" lt 1 
set key outside spacing 1 font ",25" sample 3
set ylabel "E-{E}_{F}/ev" offset 2.0,0
set xlabel "KPATH" offset 0,0.5
plot "../BD/BAND.dat" w l lw 3 lc rgb "blue" title "dft" ,\
"wannier90_band.dat" u 1:($2 - 6.3) with line linetype 1 linewidth 3 linecolor rgb "red" title "wannier"
