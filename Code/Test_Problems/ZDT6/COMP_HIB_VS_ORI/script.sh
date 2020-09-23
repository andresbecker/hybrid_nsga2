#!/bin/bash
# -*- ENCODING: UTF-8 -*-

#DETENEMOS SINCRONIZACION DE DROPBOX
dropbox stop

echo "Para ejecutar este Script solamente es necesario modificar el archivo INPUT_DATA/$f_prueba.in"

rm -f alfa.out medidas_o.out medidas_g.out GRAFICAS/*_g.dat GRAFICAS/*_o.dat hibri_vs_ori.out
b=0.123
n_eje=30
f_prueba="zdt6"
prom_gd_H=0
prom_igd_H=0
prom_deltap_H=0
prom_gd_O=0
prom_igd_O=0
prom_deltap_O=0
flag_gd=0
flag_igd=0
flag_deltap=0
eval=0
j=0
aumento=1

#DETERMINAMOS CUANTAS VECES SE LLEVARA A CABO EL PRIMER BUQLUE
final=$(awk '{if(NR==5) print (5+$0+14)}' INPUT_DATA/$f_prueba.in)
n_gen=$(awk '{if(NR==2){print $0}}' INPUT_DATA/$f_prueba.in)

#EJECUTAMOS EL HIBRIDO
for ((i=1; i<=n_eje ; i++))
do
  echo "REALIZANDO EJECUCION DEL HIBRIDO: $i DE $n_eje ";
  semilla=$(bc <<< 'scale=8; '$b' + (1-'$b')*('$i'-1)/'$n_eje'')
  ./nsga2r $semilla < INPUT_DATA/$f_prueba.in > borrar.txt
  echo "Numero de evaluaciones NSGA-II Híbrido:"
  awk '{print $1+$2}' eval.out
  awk '{if($0!~/#/)print $1, $2;}' final_pop.out > frente_aprox.dat
  #Calculamos GD e IGD para NSGA-II
  ./h < INPUT_DATA/h.in >> medidas_g.out
done
#Generamos archivos para graficas de caja del NSGA-II Hibrido
echo "medidas_g.out" | ./qs $n_eje
cp gd_box.out gd_g.out
cp igd_box.out igd_g.out
cp deltap_box.out deltap_g.out
rm gd_box.out igd_box.out deltap_box.out

#EXTRAEMOS LOS PROMEDIOS DEL HIBRIDO
prom_gd_H=$(awk '{if($0!~/#/) print $6}' gd_g.out)
prom_igd_H=$(awk '{if($0!~/#/) print $6}' igd_g.out)
prom_deltap_H=$(awk '{if($0!~/#/) print $6}' deltap_g.out)


#Ejecuta el NSGA-II Original AUMENTANDO LAS GENERACIONES
while [ "$(bc <<< ''$flag_gd'+'$flag_igd'+'$flag_deltap'')" -lt 3 ];
do
    rm -f medidas_o.out

    for ((i=1; i<=n_eje ; i++))
    do
	semilla=$(bc <<< 'scale=8; '$b' + (1-'$b')*('$i'-1)/'$n_eje'')
        echo "REALIZANDO EJECUCION del NSGA-II: $j.$i ";

        awk '{if(NR==2){gen=($0+'$j'); print gen;}else{if(NR=='$final'){print gen+1}else{print $0}}}' INPUT_DATA/$f_prueba.in | ./nsga2r $semilla > borrar.txt
        echo "Numero de evaluaciones NSGA-II original:"
        awk '{print $1+$2}' eval.out
	#CORROBORA QUE NO 
	eval=$(awk '{print $1+$2}' eval.out)
        awk '{if($0!~/#/)print $1, $2;}' final_pop.out > frente_aprox.dat
        #Calculamos GD e IGD para NSGA-II+MDPM
        ./h < INPUT_DATA/h.in >> medidas_o.out

    done
    rm -f alfa.out
    #Generamos archivos para graficas de caja del NSGA-II+MDPM
    echo "medidas_o.out" | ./qs $n_eje
    cp gd_box.out gd_o.out
    cp igd_box.out igd_o.out
    cp deltap_box.out deltap_o.out
    rm gd_box.out igd_box.out deltap_box.out
    
    #EXTRAEMOS LOS PROMEDIOS DEL ORIGINAL
    prom_gd_O=$(awk '{if($0!~/#/) print $6}' gd_o.out)
    prom_igd_O=$(awk '{if($0!~/#/) print $6}' igd_o.out)
    prom_deltap_O=$(awk '{if($0!~/#/) print $6}' deltap_o.out)

    #CHECAMOS SI EL PROMEDIO DEL ORIGINAL YA A ALCANZADO AL DEL HIBRIDO
    if [ $(awk 'BEGIN{if('$prom_gd_O'<='$prom_gd_H'){print "1"}else {print 0}}') -eq 1 ] && [ "$flag_gd" -eq 0 ];
    then
	flag_gd=1;
        echo "Generacion GD=$(bc <<< ''$n_gen'+'$j'')" >> hibri_vs_ori.out
	awk '{if($0!~/\#/){print '$j'+1, $0, "HIBRIDO"}}' gd_g.out > GRAFICAS/gd_g.dat
    fi
    if [ $(awk 'BEGIN{if('$prom_igd_O'<='$prom_igd_H'){print "1"}else {print 0}}') -eq 1 ] && [ "$flag_igd" -eq 0 ];
    then
	flag_igd=1;
        echo "Generacion IGD=$(bc <<< ''$n_gen'+'$j'')" >> hibri_vs_ori.out
	awk '{if($0!~/\#/){print '$j'+1, $0, "HIBRIDO"}}' igd_g.out > GRAFICAS/igd_g.dat
    fi
    if [ $(awk 'BEGIN{if('$prom_deltap_O'<='$prom_deltap_H'){print "1"}else {print 0}}') -eq 1 ] && [ "$flag_deltap" -eq 0 ];
    then
	flag_deltap=1;
        echo "Generacion Delta_p=$(bc <<< ''$n_gen'+'$j'')" >> hibri_vs_ori.out
	awk '{if($0!~/\#/){print '$j'+1, $0, "HIBRIDO"}}' deltap_g.out > GRAFICAS/deltap_g.dat
    fi

    #GENERAMOS ARCHIVOS PARA GRAFICO
    awk '{if($0!~/\#/){print '$j'+1, $0, "G"('$n_gen'+'$j')}}' gd_o.out >> GRAFICAS/gd_o.dat
    awk '{if($0!~/\#/){print '$j'+1, $0, "G"('$n_gen'+'$j')}}' igd_o.out >> GRAFICAS/igd_o.dat
    awk '{if($0!~/\#/){print '$j'+1, $0, "G"('$n_gen'+'$j')}}' deltap_o.out >> GRAFICAS/deltap_o.dat

j=$(bc <<< ''$j'+'$aumento'')
done

cd GRAFICAS/
gnuplot *.gnu
cd ..

rm -f all_pop.out best_pop.out borrar.txt eval.out final_pop.out frente_aprox.dat initial_pop.out params.out qs.out

#REANUDAMOS SINCRONIZACION DE DROPBOX
dropbox start
