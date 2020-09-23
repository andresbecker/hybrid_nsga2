set terminal png nocrop enhanced size 800,800 font "arial,10" 
set output 'FRENTES.png'
set samples 800, 800 
set title "{/=15 FRENTES}"
set ylabel "{/=13 f_2}" 
set xlabel "{/=13 f_1}"

#FIJA EL ANCHO DE LAS CAJAS
#set bar 10.000000

#set yrange [ -0.8 : -0.4 ] noreverse nowriteback
#set xrange [ 0.3 : 0.5 ] noreverse nowriteback

plot 'pareto.dat' title 'FRENTE DE PARETO' with points lc rgb "black", 'g.dat' with points lc rgb "red" lt 5 title 'NSGA-II HIBRIDO', 'o.dat' with points lc rgb "cyan" lt 5 title 'NSGA-II ORIGINAL'
