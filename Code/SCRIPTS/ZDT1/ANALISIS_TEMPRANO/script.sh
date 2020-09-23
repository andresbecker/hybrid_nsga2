#!/bin/bash
# -*- ENCODING: UTF-8 -*-

echo "Para ejecutar este Script solamente es necesario modificar el archivo INPUT_DATA/zdt1.in"

rm -f alfa.out GRAFICAS/gd_g.dat GRAFICAS/gd_o.dat GRAFICAS/igd_g.dat GRAFICAS/igd_o.dat

b=0.123
n_eje=30
t_frente=500

#DETERMINAMOS CUANTAS VECES SE LLEVARA A CABO EL PRIMER BUQLUE
n_gen=$(awk '{if(NR==2){print $0}}' INPUT_DATA/zdt1.in)
gen_grad=$(awk '{if(NR==5)final=5+$0+14; if(NR==final)print $0;}' INPUT_DATA/zdt1.in)
n_iter=$(expr $n_gen - $gen_grad + 1 )

for ((j=1; j<=n_iter ; j++))
do
    rm -f medidas_g.out medidas_o.out
    #LLEVA CONTEO DE LAS ITERACIONES
    echo $j > temp.out
    #GENERAMOS ARCHIVO DE ENTRADA PARA CADA CICLO
    awk -f AWK/zdt1_iter.awk INPUT_DATA/zdt1.in > INPUT_DATA/zdt1_iter.in

    for ((i=1; i<=n_eje ; i++))
    do
	semilla=$(bc <<< 'scale=8; '$b' + (1-'$b')*('$i'-1)/'$n_eje'')
        echo "REALIZANDO EJECUCION: "$j.$i "de" $n_iter.$n_eje" faltantes";
        #Ejecuta el NSGA-II+MDPM
        ./nsga2r $semilla < INPUT_DATA/zdt1_iter.in > borrar.txt
        echo "Numero de evaluaciones NSGA-II+MDPM:"
        awk '{print $0}' eval.out
        #Extrae la ultima poblacion para sacar GD e IGD
        awk '{if($0!~/#/)print $1, $2;}' final_pop.out > frente_aprox.dat
        #Calculamos GD e IGD para NSGA-II+MDPM
        ./h $t_frente < INPUT_DATA/h.in >> medidas_g.out
	#CREAMOS ARCHIVO PARA EJECUTAR NSGA-II ORIGINAL CON IGUAL NUMERO DE EVALUACIONES DE LA FUNCION OBJETIVO
	awk -f AWK/archivo_nsga2.awk INPUT_DATA/zdt1.in > INPUT_DATA/zdt1_o.in
        ./nsga2r $semilla < INPUT_DATA/zdt1_o.in > borrar.txt
        echo "Numero de evaluaciones NSGA-II original:"
        awk '{print $0}' eval.out
        awk '{if($0!~/#/)print $1, $2;}' final_pop.out > frente_aprox.dat
        #Calculamos GD e IGD para NSGA-II
        ./h $t_frente < INPUT_DATA/h.in >> medidas_o.out
    done
    rm -f alfa.out
    #Generamos archivos para graficas de caja del NSGA-II+MDPM
    echo "medidas_g.out" | ./qs1 $n_eje
    cp gd_box.out gd_g.out
    cp igd_box.out igd_g.out
    rm gd_box.out igd_box.out

    #Generamos archivos para graficas de caja del NSGA-II
    echo "medidas_o.out" | ./qs1 $n_eje
    cp gd_box.out gd_o.out
    cp igd_box.out igd_o.out
    rm gd_box.out igd_box.out

    #GENERAMOS ARCHIVOS PARA GRAFICO
    awk '{if($0!~/\#/){print '$j', $0, "G"('$gen_grad'+'$j'-1)}}' gd_o.out >> GRAFICAS/gd_o.dat
    awk '{if($0!~/\#/){print ('$j'), $0, "G"('$gen_grad'+'$j'-1)}}' gd_g.out >> GRAFICAS/gd_g.dat
    awk '{if($0!~/\#/){print '$j', $0, "G"('$gen_grad'+'$j'-1)}}' igd_o.out >> GRAFICAS/igd_o.dat
    awk '{if($0!~/\#/){print ('$j'), $0, "G"('$gen_grad'+'$j'-1)}}' igd_g.out >> GRAFICAS/igd_g.dat
done

cd GRAFICAS/
gnuplot *.gnu
cd ..

rm -f all_pop.out best_pop.out borrar.txt eval.out final_pop.out frente_aprox.dat initial_pop.out params.out qs.out temp.out
