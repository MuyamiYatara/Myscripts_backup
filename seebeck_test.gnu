set encoding iso_8859_1
set terminal pdfcairo enhanced color font ",30" size 10, 6
set output 'sigma.pdf'
set border lw 2
set autoscale fix
set ylabel '{/Symbol s}_{xx}/{/Symbol t} (({/Symbol W}*m*s)^{-1})'
set xlabel 'B{/Symbol t} (T.ps)'
set key outside
set palette defined (0 'red', 1 'green')
unset colorbox
set ylabel offset 0.0,0
Tmin =  30.00
Tmax = 300.00
NumT =   31
plot_NumT = 5
OmegaMin =  -0.01
OmegaMax =   0.01
OmegaNum =    3
lw =    4

#plot conductivity/tau
set ylabel '{S}_{xx}/{/Symbol t} (({/Symbol W}*m*s)^{-1})'
set xlabel 'B{/Symbol t} (T.ps)'
plot for [i=0:plot_NumT] 'seebeck_total_mu_0.000eV.dat' every :::i::i+1 u 1:3 w l lw lw lt palette frac i/(plot_NumT*1.0)title sprintf('T=%.0f K', Tmin + (Tmax-Tmin)/(NumT*1.0-1.0)*i)
 
set ylabel '{S}/{/Symbol t}'
set xlabel 'T (K)'
set xrange [0:300]
set yrange [0:-0.00012]
plot 'Tdata_output.txt' u 1:4 w l lw lw lt palette frac 0 title 'B=0,S_{xx}',\
'Tdata_output.txt' u 1:12 w l lw lw lt palette frac 1 title 'B=0,S_{zz}'
 
#plot resistivity*tau
 
#set ylabel '{/Symbol r}_{xx}*{/Symbol t} ({/Symbol W}*m*s)'
#plot for [i=0:NumT-1] 'rho_total_mu_0.00eV.dat' every :::i::i+1 u 1:2 w l lw lw lt palette frac i/(NumT*1.0)title sprintf('T=%.0f K', Tmin + (Tmax-Tmin)/(NumT*1.0-1.0)*i)
 
#set ylabel '{/Symbol r}_{yx}*{/Symbol t} ({/Symbol W}*m*s)'
#plot for [i=0:NumT-1] 'rho_total_mu_0.00eV.dat' every :::i::i+1 u 1:5 w l lw lw lt palette frac i/(NumT*1.0)title sprintf('T=%.0f K', Tmin + (Tmax-Tmin)/(NumT*1.0-1.0)*i)
