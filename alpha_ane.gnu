set terminal pdfcairo enhanced color font ",30" size 13, 6
set output 'ane.pdf'
set key outside
set palette defined (0 'red', 1 'green')
unset colorbox
Tmin =  20.00
Tmax = 320.00
NumT =   21
OmegaMin =  -1.00
OmegaMax =   1.00
OmegaNum =  501
#control the line we start for counting chemical potential
i_offset = 250

set arrow from graph 0, first 0 to graph 1, first 0 nohead lw 4 lc rgb "blue" lt 1

#plot the temperature-dependent alpha_yx for the first six chemical potentials
set xlabel "T (K)"
set ylabel "{/Symbol a}_{yx} (A/(mK))"
set ylabel offset 0.0,0
plot for [i=(0+i_offset):(5+i_offset)] 'alpha_ane_eta10.00meV.txt' every ::(i*NumT)::i*NumT+(NumT-1)u 2:(-$3) w l lt palette frac (i-i_offset)/5. title sprintf('{/Symbol m}=%.3f eV', OmegaMin+(OmegaMax-OmegaMin)/(OmegaNum*1.0-1.0)*i )

#plot the chemical potential dependent alpha_yx
set xlabel "E-E_f (eV)"
set ylabel "{/Symbol a}_{yx} (A/(mK))"
set ylabel offset 0.0,0
plot for [i=0:NumT-1] 'alpha_ane_eta10.00meV.txt' every NumT::i u 1:(-$3) w l lt palette frac i/(NumT*1.0-1.0) title sprintf('T=%.3f K',Tmin+(Tmax-Tmin)/(NumT*1.0-1.0)*i)
