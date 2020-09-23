#!/bin/bash
# -*- ENCODING: UTF-8 -*-

rm -f medidas_g.out medidas_o.out alfa.out GRAFICAS/*.dat

b=0.123
n_eje=30
f_prueba="zdt1"
final=$(awk '{if(NR==5) print (5+$0+14)}' INPUT_DATA/$f_prueba.in)
t_frente=$(awk '{if(NR==1) pop=$0; if(NR==2){print $0*pop}}' INPUT_DATA/$f_prueba.in)
comp=225

for((i=1; i<=n_eje; i++))
do
    semilla=$(bc <<< 'scale=8; '$b' + (1-'$b')*('$i'-1)/'$n_eje'')
    echo "REALIZANDO EJECUCION: "$i;
    #Ejecuta el NSGA-II+MDPM
    ./nsga2r $semilla < INPUT_DATA/$f_prueba.in > borrar.txt
    echo "Numero de evaluaciones NSGA-II+MDPM:"
    awk '{print $1+$2}' eval.out
    #Extrae la ultima poblacion para sacar GD e IGD
    awk '{if($0!~/#/)print $1, $2;}' final_pop.out > frente_aprox.dat
    #Calculamos GD e IGD para NSGA-II+MDPM
    ./h < INPUT_DATA/h.in >> medidas_g.out

    #Ejecuta el NSGA-II ORIGINAL
    awk '{if(NR==2){gen='$comp'; print gen;}else{if(NR=='$final'){print gen+1}else{print $0}}}' INPUT_DATA/$f_prueba.in | ./nsga2r $semilla > borrar.txt
    echo "Numero de evaluaciones NSGA-II original:"
    awk '{print $1+$2}' eval.out
    awk '{if($0!~/#/)print $1, $2;}' final_pop.out > frente_aprox.dat
    #Calculamos GD e IGD para NSGA-II
    ./h < INPUT_DATA/h.in >> medidas_o.out
done


#Generamos archivos para graficas de caja del NSGA-II+MDPM
echo "medidas_g.out" | ./qs $n_eje
cp gd_box.out gd_g.out
cp igd_box.out igd_g.out
cp deltap_box.out deltap_g.out
rm gd_box.out igd_box.out deltap_box.out

#Generamos archivos para graficas de caja del NSGA-II
echo "medidas_o.out" | ./qs $n_eje
cp gd_box.out gd_o.out
cp igd_box.out igd_o.out
cp deltap_box.out deltap_o.out
rm gd_box.out igd_box.out deltap_box.out

#GENERAMOS GRAFICO DE CAJA
awk '{if($0!~/\#/){print 1, $0, "ORI"}}' gd_o.out >> GRAFICAS/gd.dat
awk '{if($0!~/\#/){print 2, $0, "HIB"}}' gd_g.out >> GRAFICAS/gd.dat
awk '{if($0!~/\#/){print 1, $0, "ORI"}}' igd_o.out >> GRAFICAS/igd.dat
awk '{if($0!~/\#/){print 2, $0, "HIB"}}' igd_g.out >> GRAFICAS/igd.dat
awk '{if($0!~/\#/){print 1, $0, "ORI"}}' deltap_o.out >> GRAFICAS/deltap.dat
awk '{if($0!~/\#/){print 2, $0, "HIB"}}' deltap_g.out >> GRAFICAS/deltap.dat

#GENERAMOS GRAFICO DEL FRENTE
#PARA GD
#awk '{if($0!~/#/){print $1}}' gd_g.out > borrar.txt
#i=$(awk 'BEGIN{getline ejec < "borrar.txt"}{if($1==ejec){print NR}}' medidas_g.out)
#PARA IGD
awk '{if($0!~/#/){print $1}}' igd_g.out > borrar.txt
i=$(awk 'BEGIN{getline ejec < "borrar.txt"}{if($2==ejec){print NR}}' medidas_g.out)
semilla=$(bc <<< 'scale=8; '$b' + (1-'$b')*('$i'-1)/'$n_eje'')
echo Mejor semilla IGD=$semilla > porc_eval.out
./nsga2r $semilla < INPUT_DATA/$f_prueba.in > borrar.txt
awk '{if($0!~/#/)print $1, $2;}' final_pop.out > GRAFICAS/g.dat
awk '{if(NR==2){gen=($0+'$comp'); print gen;}else{if(NR=='$final'){print gen+1}else{print $0}}}' INPUT_DATA/$f_prueba.in | ./nsga2r $semilla > borrar.txt
awk '{if($0!~/#/)print $1, $2;}' final_pop.out > GRAFICAS/o.dat
cp pareto.dat GRAFICAS/
cd GRAFICAS/
gnuplot *.gnu
eog FRENTE.png 
cd ..

rm -f all_pop.out best_pop.out borrar.txt final_pop.out frente_aprox.dat initial_pop.out params.out qs.out eval.out
