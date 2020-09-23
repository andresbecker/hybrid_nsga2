#!/bin/bash
# -*- ENCODING: UTF-8 -*-
rm -f alfa.out GRAFICAS/*.dat

n_iter=199
n_eje=30
b=0.123
f_prueba="zdt3"

final=$(awk '{if(NR==5) print (5+$0+14)}' INPUT_DATA/$f_prueba.in)

for ((j=0; j<=n_iter ; j++))
do
    rm -f medidas_o.out

    for ((i=1; i<=n_eje ; i++))
    do
	semilla=$(bc <<< 'scale=8; '$b' + (1-'$b')*('$i'-1)/'$n_eje'')
        echo "REALIZANDO EJECUCION: "$j.$i "de" $n_iter.$n_eje;

	awk '{if(NR==2){print $0+'$j'}else{if(NR=='$final'){print ($0+'$j'+1)}else{print $0}}}' INPUT_DATA/$f_prueba.in | ./nsga2r $semilla > borrar.txt	
        echo "Numero de evaluaciones NSGA-II original:"
        awk '{print $1}' eval.out
        awk '{if($0!~/#/)print $1, $2;}' final_pop.out > frente_aprox.dat
        #Calculamos GD e IGD para NSGA-II
        ./h < INPUT_DATA/h.in >> medidas_o.out
    done
    rm -f alfa.out
    #Generamos archivos para graficas de caja del NSGA-II
    echo "medidas_o.out" | ./qs $n_eje
    cp gd_box.out gd_o.out
    cp igd_box.out igd_o.out
    cp deltap_box.out deltap_o.out
    rm gd_box.out igd_box.out deltap_box.out

    #GENERAMOS ARCHIVOS PARA GRAFICO
    awk '{if(NR==2){print $0}}' INPUT_DATA/$f_prueba.in > borrar.txt
    awk 'BEGIN{getline gen < "borrar.txt"}{if($0!~/\#/){print '$j', $0, "G"(gen+'$j')}}' gd_o.out >> GRAFICAS/gd_o.dat
    awk 'BEGIN{getline gen < "borrar.txt"}{if($0!~/\#/){print '$j', $0, "G"(gen+'$j')}}' igd_o.out >> GRAFICAS/igd_o.dat
    awk 'BEGIN{getline gen < "borrar.txt"}{if($0!~/\#/){print '$j', $0, "G"(gen+'$j')}}' deltap_o.out >> GRAFICAS/deltap_o.dat
done

cd GRAFICAS/
gnuplot *.gnu
cd ..

rm -f all_pop.out best_pop.out borrar.txt final_pop.out frente_aprox.dat initial_pop.out params.out qs.out
