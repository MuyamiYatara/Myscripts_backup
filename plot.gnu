#set terminal  postscript enhanced color font ",30"
#set output 'spectrum_unfold.eps'
set terminal pngcairo enhanced color font ",60" size 1920,1680
set palette defined ( 0 "white", 1  "#D72F01" )
set output 'spectrum_unfold_kpath.png'
#set style data linespoints
#set size 0.9, 1
#set origin 0.05,0
#unset key
#set border lw 3
#set view map
##set xtics font ",24"
##set ytics font ",24"
##set ylabel font ",24"
##set ylabel offset 1.5,0
#emin=   -0.500000
#emax=    0.500000
#set xrange [0:    2.5]
#set ylabel "Energy (eV)"
set yrange [ -0.2 : 0.2 ]
#set xtics ("H  "    0.00000,"G  "    0.46391,"F  "    0.92781,"G  "    1.39172,"L  "    1.85562)
#set arrow from    0.46391, emin to    0.46391, emax nohead front lw 3
#set arrow from    0.92781, emin to    0.92781, emax nohead front lw 3
#set arrow from    1.39172, emin to    1.39172, emax nohead front lw 3
#set colorbox
#set pm3d interpolate 2,2
plot 'BAND.dat' u 1:2 w l
