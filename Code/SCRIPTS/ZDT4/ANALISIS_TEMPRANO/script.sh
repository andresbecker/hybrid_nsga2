#!/bin/bash
# -*- ENCODING: UTF-8 -*-

#DETENEMOS SINCRONIZACION DE DROPBOX
dropbox stop

echo "Para ejecutar este Script solamente es necesario modificar el archivo INPUT_DATA/$f_prueba.in"

rm -f alfa.out medidas_o.out GRAFICAS/*_g.dat GRAFICAS/*_o.dat
b=0.123
n_eje=30
f_prueba="zdt4"
min_gd=0
min_igd=0
min_deltap=0
temp=0
eval=0
error=0

#DETERMINAMOS CUANTAS VECES SE LLEVARA A CABO EL PRIMER BUQLUE
final=$(awk '{if(NR==5) print (5+$0+14)}' INPUT_DATA/$f_prueba.in)
pop_size=$(awk '{if(NR==1){print $0}}' INPUT_DATA/$f_prueba.in)
n_gen=$(awk '{if(NR==2){print $0}}' INPUT_DATA/$f_prueba.in)
gen_grad=$(awk '{if(NR=='$final')print $0;}' INPUT_DATA/$f_prueba.in)
n_iter=$(expr $n_gen - $gen_grad )

for ((j=0; j<n_iter && error==0 ; j++))
do
    rm -f medidas_g.out

    for ((i=1; i<=n_eje ; i++))
    do
	semilla=$(bc <<< 'scale=8; '$b' + (1-'$b')*('$i'-1)/'$n_eje'')
        echo "REALIZANDO EJECUCION: "$j.$i "de" $n_iter.$n_eje;
        #Ejecuta el NSGA-II+MDPM
        awk '{if(NR=='$final'){print $0+'$j'}else print $0}' INPUT_DATA/$f_prueba.in | ./nsga2r $semilla > borrar.txt
        echo "Numero de evaluaciones NSGA-II+MDPM:"
        awk '{printf("NSGA-II: %d,  MDPM: %d,  TOTAL: %d\n", $1, $2, $1+$2)}' eval.out
	#CORROBORA QUE NO 
	eval=$(awk '{print $1+$2}' eval.out)
	if [ "$eval" -gt "$(expr $n_gen \* $pop_size + $pop_size )" ]
	then
	    echo "ERROR EN LA GEN "$(expr $gen_grad + $j )
	    error=1
	fi
        #Extrae la ultima poblacion para sacar GD e IGD
        awk '{if($0!~/#/)print $1, $2;}' final_pop.out > frente_aprox.dat
        #Calculamos GD e IGD para NSGA-II+MDPM
        ./h < INPUT_DATA/h.in >> medidas_g.out
	if [ "$j" -eq 0 ]
    	then
	    awk '{if(NR==2){gen=$0; print $0;}else{if(NR=='$final'){print gen+1}else{print $0}}}' INPUT_DATA/$f_prueba.in | ./nsga2r $semilla > borrar.txt
	    echo "Numero de evaluaciones NSGA-II original:"
	    awk '{print $1+$2}' eval.out
	    awk '{if($0!~/#/)print $1, $2;}' final_pop.out > frente_aprox.dat
	    #Calculamos GD e IGD para NSGA-II
	    ./h < INPUT_DATA/h.in >> medidas_o.out
	fi
    done
    rm -f alfa.out
    #Generamos archivos para graficas de caja del NSGA-II+MDPM
    echo "medidas_g.out" | ./qs $n_eje
    cp gd_box.out gd_g.out
    cp igd_box.out igd_g.out
    cp deltap_box.out deltap_g.out
    rm gd_box.out igd_box.out deltap_box.out
    
    #BUSCAMOS LOS MEJORES VALORES
    if [ "$j" -eq 0 ]
    then
	min_gd=$(awk '{if($0!~/#/) print $6}' gd_g.out)
	gd_iter=$gen_grad
        min_igd=$(awk '{if($0!~/#/) print $6}' igd_g.out)
	igd_iter=$gen_grad
        min_deltap=$(awk '{if($0!~/#/) print $6}' deltap_g.out)
	deltap_iter=$gen_grad
    else
	temp=$(awk '{if($0!~/#/) print $6}' gd_g.out)
        if [ $(awk 'BEGIN{if('$temp'<'$min_gd'){print "1"}else {print 0}}') -eq 1 ]
	then
	    min_gd=$temp
	    gd_iter=$(bc <<< ''$gen_grad'+'$j'')
	fi
	temp=$(awk '{if($0!~/#/) print $6}' igd_g.out)
        if [ $(awk 'BEGIN{if('$temp'<'$min_igd'){print "1"}else {print 0}}') -eq 1 ]
	then
	    min_igd=$temp
	    igd_iter=$(bc <<< ''$gen_grad'+'$j'')
	fi
	temp=$(awk '{if($0!~/#/) print $6}' deltap_g.out)
        if [ $(awk 'BEGIN{if('$temp'<'$min_deltap'){print "1"}else {print 0}}') -eq 1 ]
	then
	    min_deltap=$temp
	    deltap_iter=$(bc <<< ''$gen_grad'+'$j'')
	fi
    fi
    #GENERAMOS ARCHIVOS PARA GRAFICO
    awk '{if($0!~/\#/){print '$j'+1, $0, "G"('$gen_grad'+'$j')}}' gd_g.out >> GRAFICAS/gd_g.dat
    awk '{if($0!~/\#/){print '$j'+1, $0, "G"('$gen_grad'+'$j')}}' igd_g.out >> GRAFICAS/igd_g.dat
    awk '{if($0!~/\#/){print '$j'+1, $0, "G"('$gen_grad'+'$j')}}' deltap_g.out >> GRAFICAS/deltap_g.dat
done

#GUARDAMOS MEJORES VALORES
echo "GD=$min_gd    IGD=$min_igd    DELTAp=$min_deltap" > mejores.out
echo "GD_gen=$gd_iter    IGD_gen=$igd_iter    DELTAp_gen=$deltap_iter" >> mejores.out

#Generamos archivos para graficas de caja del NSGA-II
echo "medidas_o.out" | ./qs $n_eje
cp gd_box.out gd_o.out
cp igd_box.out igd_o.out
cp deltap_box.out deltap_o.out
rm gd_box.out igd_box.out deltap_box.out
awk '{if($0!~/\#/){print '$j'+'$gen_grad'+1, $0, "NSGA-II"}}' gd_o.out > GRAFICAS/gd_o.dat
awk '{if($0!~/\#/){print '$j'+'$gen_grad'+1, $0, "NSGA-II"}}' igd_o.out > GRAFICAS/igd_o.dat
awk '{if($0!~/\#/){print '$j'+'$gen_grad'+1, $0, "NSGA-II"}}' deltap_o.out > GRAFICAS/deltap_o.dat

cd GRAFICAS/
gnuplot *.gnu
cd ..

rm -f all_pop.out best_pop.out borrar.txt eval.out final_pop.out frente_aprox.dat initial_pop.out params.out qs.out temp.out

#REANUDAMOS SINCRONIZACION DE DROPBOX
dropbox start
