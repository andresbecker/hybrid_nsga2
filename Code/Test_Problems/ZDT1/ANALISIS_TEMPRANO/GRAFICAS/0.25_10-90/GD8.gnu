set terminal png nocrop enhanced size 900,500 font "arial,10" 
set output 'GD_gen_8.png'
set samples 800, 800 
set title "{/=15 GD con p=2 y 30 Ejecuciones para ZDT1}"
set ylabel "{/=13 GD}" 
set xlabel "{/=13 GENERACIONES DONDE SE APLICO EL MDPM}"

#FIJA EL ANCHO DE LAS CAJAS
set bar 10.000000

#set yrange [ 0 : * ] noreverse nowriteback
set xrange [ 60 : 82 ] noreverse nowriteback

#GRAFICA DE PUNTOS
plot 'gd_g.dat' using 1:7:xticlabel(8) with points lc "red" lw 2 title 'NSGA-II HIBRIDO', 'gd_o.dat' using 1:7 with points lc "black" lw 2 title 'NSGA-II ORIGINAL'

#GRAFICA DE CAJA
#plot 'gd_g.dat' using 1:4:2:3:6:xticlabels(8) with candlesticks lw 2 title 'Rango Intercuartil' whiskerbars, 'gd_g.dat' using 1:7:7:7:7 with candlesticks lc "red" lw 2 title 'Promedio', 'gd_g.dat' using 1:5:5:5:5 with candlesticks lw 2 title 'Mediana', 'gd_o.dat' using 1:4:2:3:6:xticlabels(8) with candlesticks lw 2 title 'Rango Intercuartil' whiskerbars, 'gd_o.dat' using 1:7:7:7:7 with candlesticks lc "red" lw 2 title 'Promedio', 'gd_o.dat' using 1:5:5:5:5 with candlesticks lw 2 title 'Mediana'
