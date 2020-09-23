set terminal png nocrop enhanced size 500,500 font "arial,12" 
set output 'comparativo.png'
#set key bmargin left horizontal Right noreverse enhanced autotitle box lt black linewidth 1.000 dashtype solid
set grid
#set style data linespoints
#set samples 800, 800 
set title "NSGA-II vs. HIBRIDO"
set ylabel "# de evaluaciones necesarias de F" 
set xlabel "{/Symbol e}"
set yrange [ 8000 : 30000 ] noreverse nowriteback
#set xrange [ 0 : 7 ] noreverse nowriteback 
x = 0.0
set key left

plot 'epsilon_hibrido.out' using 1:3:xticlabels(2) pt 7 lc "red" title 'HIBRIDO', 'epsilon_nsga2.out' using 1:3:xticlabels(2) pt 9 lc "blue" title 'NSGA-II'


