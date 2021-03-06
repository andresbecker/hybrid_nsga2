set terminal png nocrop enhanced size 900,500 font "arial,10" 
set output 'DELTAp.png'
set samples 800, 800 
set title "{/=15 {/=18 {/Symbol D}_p}}"
set ylabel "{/=13 {/Symbol D}_p}" 
set xlabel "{/=13 ALGORITMO}"

#FIJA EL ANCHO DE LAS CAJAS
set bar 10.000000

#set yrange [ 0 : * ] noreverse nowriteback
set xrange [ 0 : 3 ] noreverse nowriteback

plot 'deltap.dat' using 1:4:2:3:6:xticlabels(8) with candlesticks lw 2 title 'Rango Intercuartil' whiskerbars, '' using 1:7:7:7:7 with candlesticks lc "red" lw 2 title 'Promedio', '' using 1:5:5:5:5 with candlesticks lw 2 title 'Mediana'
