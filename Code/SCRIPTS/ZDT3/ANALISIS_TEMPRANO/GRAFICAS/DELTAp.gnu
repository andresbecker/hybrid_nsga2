set terminal png nocrop enhanced size 900,500 font "arial,10" 
set output 'DELTAp.png'
set samples 800, 800 
set title "{/=15 DELTAp con p=2 y 30 Ejecuciones para ZDT6}"
set ylabel "{/=13 DELTAp}" 
set xlabel "{/=13 GENERACIONES DONDE SE APLICO EL MDPM}"

#FIJA EL ANCHO DE LAS CAJAS
set bar 10.000000

#set yrange [ 0 : * ] noreverse nowriteback
set xrange [ 0 : 50 ] noreverse nowriteback

#GRAFICA DE PUNTOS
plot 'deltap_g.dat' using 1:7 with points lc "red" lw 2 title 'NSGA-II HIBRIDO', 'deltap_o.dat' using 1:7 with points lc "black" lw 2 title 'NSGA-II ORIGINAL'

#GRAFICA DE CAJA
#plot 'deltap_g.dat' using 1:4:2:3:6:xticlabels(8) with candlesticks lw 2 title 'Rango Intercuartil' whiskerbars, 'deltap_g.dat' using 1:7:7:7:7 with candlesticks lc "red" lw 2 title 'Promedio', 'deltap_g.dat' using 1:5:5:5:5 with candlesticks lw 2 title 'Mediana', 'deltap_o.dat' using 1:4:2:3:6:xticlabels(8) with candlesticks lw 2 title 'Rango Intercuartil' whiskerbars, 'deltap_o.dat' using 1:7:7:7:7 with candlesticks lc "red" lw 2 title 'Promedio', 'deltap_o.dat' using 1:5:5:5:5 with candlesticks lw 2 title 'Mediana'
