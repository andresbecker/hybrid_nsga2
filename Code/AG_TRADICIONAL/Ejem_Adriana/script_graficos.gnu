do for [i=0: 100] {

set terminal png nocrop enhanced size 500,355 font "arial,8" 
set output 'n'.i.'.png'
#set key bmargin left horizontal Right noreverse enhanced autotitle box lt black linewidth 1.000 dashtype solid
set samples 800, 800 
set title "{/=15 f(x)=(5{/Symbol p}/2x)^2sin(x)  y  P_{".i."}}"
set ylabel "f(x)" 
set xlabel "Eje X"
set yrange [ -.6 : 1.1 ] noreverse nowriteback
set xrange [2*pi:9*pi] noreverse nowriteback
x = 0.0

plot (25*pi*pi/(4*x*x))*sin(x) title 'f(x)', "<awk '{if($1==".i."){print $3, $5}}' report.out" pt 5 title 'P_'.i.'

}
