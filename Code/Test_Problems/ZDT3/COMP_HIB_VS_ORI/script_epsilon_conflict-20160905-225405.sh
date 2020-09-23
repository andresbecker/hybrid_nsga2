#!/bin/bash
# -*- ENCODING: UTF-8 -*-

#DETENEMOS SINCRONIZACION DE DROPBOX
dropbox stop

rm -f alfa.out medidas_o.out medidas_g.out GRAFICAS/*_g.dat GRAFICAS/*_o.dat

b=0.123
n_eje=30
f_prueba="zdt3"
prom_deltap_O=0
prom_deltap_H=0
epsilon1=0.1
epsilon2=0.0816
epsilon3=0.0632
epsilon4=0.0448
epsilon5=0.0264
epsilon6=0.008
epsilons=6
flag_o_e1=0
flag_h_e1=0
flag_o_e2=0
flag_h_e2=0
flag_o_e3=0
flag_h_e3=0
flag_o_e4=0
flag_h_e4=0
flag_o_e5=0
flag_h_e5=0
flag_o_e6=0
flag_h_e6=0
eval=0
j=0
aumento=1
cota_superior=70000

#DETERMINAMOS CUANTAS VECES SE LLEVARA A CABO EL PRIMER BUQLUE
final=$(awk '{if(NR==5) print (5+$0+14)}' INPUT_DATA/$f_prueba.in)
n_gen=$(awk '{if(NR==2){print $0}}' INPUT_DATA/$f_prueba.in)

#ARCHIVO PARA GRAFICO NSGA-II
echo "#NSGA-II" > GRAFICAS/epsilon_nsga2.out
echo "#No. EPSILON NO.EVAL" >> GRAFICAS/epsilon_nsga2.out
#ARCHIVO PARA GRAFICO HIBRIDO
echo "#HIBRIDO" > GRAFICAS/epsilon_hibrido.out
echo "#No. EPSILON NO.EVAL" >> GRAFICAS/epsilon_hibrido.out
#GENERAMOS ARCHIVO PARA GUARDAD PROMEDIO DE EVAL
echo "#No. EPSILON NSGA-II MDPM" > porc_eval.out

flag_nsga2=$(bc <<< ''$flag_o_e1'+'$flag_o_e2'+'$flag_o_e3'+'$flag_o_e4'+'$flag_o_e5'+'$flag_o_e6'')
flag_hibrido=$(bc <<< ''$flag_h_e1'+'$flag_h_e2'+'$flag_h_e3'+'$flag_h_e4'+'$flag_h_e5'+'$flag_h_e6'')

#Ejecuta el NSGA-II Original AUMENTANDO LAS GENERACIONES
while [ "$flag_nsga2" -lt "$epsilons" ] || [ "$flag_hibrido" -lt "$epsilons" ] && [ "$eval" -lt "$cota_superior" ];
do
    rm -f medidas_o.out medidas_g.out alfa.out
    flag_nsga2=$(bc <<< ''$flag_o_e1'+'$flag_o_e2'+'$flag_o_e3'+'$flag_o_e4'+'$flag_o_e5'+'$flag_o_e6'')
    flag_hibrido=$(bc <<< ''$flag_h_e1'+'$flag_h_e2'+'$flag_h_e3'+'$flag_h_e4'+'$flag_h_e5'+'$flag_h_e6'')
    t_frente=$(awk '{if(NR==1) pop=$0; if(NR==2){print (($0+'$j')*pop)}}' INPUT_DATA/$f_prueba.in)
    eval_nsga2=0
    eval_grad=0
    mean_nsga2=0
    mean_grad=0

    for ((i=1; i<=n_eje ; i++))
    do
        echo "NSGA-II no. $j.$i"
	semilla=$(bc <<< 'scale=8; '$b' + (1-'$b')*('$i'-1)/'$n_eje'')

	#HIBRIDO
	if [ "$flag_hibrido" -lt "$epsilons" ]
	then
            awk '{if(NR==2){print $0+'$j';}else{print $0}}' INPUT_DATA/$f_prueba.in | ./nsga2r $semilla > borrar.txt
            echo "Numero de evaluaciones NSGA-II Hibrido:"
            awk '{print $1+$2}' eval.out
	    eval=$(awk '{print $1+$2}' eval.out)
            awk '{if($0!~/#/)print $1, $2;}' final_pop.out > frente_aprox.dat
            #Calculamos GD e IGD para NSGA-II+MDPM
            ./h < INPUT_DATA/h.in >> medidas_g.out

    	    eval_nsga2=$(awk '{print $1}' eval.out)
    	    mean_nsga2=$(bc <<< 'scale=4; '$mean_nsga2' + '$eval_nsga2'')
    	    eval_grad=$(awk '{print $2}' eval.out)
    	    mean_grad=$(bc <<< 'scale=4; '$mean_grad' + '$eval_grad'')
	fi
	#NSGA-II
	if [ "$flag_nsga2" -lt "$epsilons" ]
	then
            awk '{if(NR==2){gen=($0+'$j'); print gen;}else{if(NR=='$final'){print gen+1}else{print $0}}}' INPUT_DATA/$f_prueba.in | ./nsga2r $semilla > borrar.txt
            echo "Numero de evaluaciones NSGA-II original:"
            awk '{print $1+$2}' eval.out
	    eval=$(awk '{print $1+$2}' eval.out)
            awk '{if($0!~/#/)print $1, $2;}' final_pop.out > frente_aprox.dat
            #Calculamos GD e IGD para NSGA-II+MDPM
            ./h < INPUT_DATA/h.in >> medidas_o.out
	fi
    done

    #Generamos archivos para graficas de caja del NSGA-II y extraemos promedio del deltap
    if [ "$flag_nsga2" -lt "$epsilons" ]
    then
	echo "medidas_o.out" | ./qs $n_eje
    	cp deltap_box.out deltap_o.out
    	rm gd_box.out igd_box.out deltap_box.out
	prom_deltap_O=$(awk '{if($0!~/#/) print $6}' deltap_o.out)
    fi

    #Generamos archivos para graficas de caja del NSGA-II HIBRIDO y extraemos promedio del deltap
    if [ "$flag_hibrido" -lt "$epsilons" ]
    then
    	echo "medidas_g.out" | ./qs $n_eje
    	cp deltap_box.out deltap_g.out
    	rm gd_box.out igd_box.out deltap_box.out
   	prom_deltap_H=$(awk '{if($0!~/#/) print $6}' deltap_g.out)
	#CALCULAMOS ARCHIVO CON PROMEDIO DE EVALUACIONES
	mean_nsga2=$(bc <<< 'scale=2; 100*'$mean_nsga2'/('$t_frente'*'$n_eje')')
	mean_grad=$(bc <<< 'scale=2; 100*'$mean_grad'/('$t_frente'*'$n_eje')')

    fi

    #COMPROBAMOS SI EL PROMEDIO ES MENOR QUE EPSILON_1 PARA NSGA-II
    if [ $(awk 'BEGIN{if('$prom_deltap_O'<='$epsilon1'){print "1"}else {print 0}}') -eq 1 ] && [ "$flag_o_e1" -eq 0 ];
    then
	flag_o_e1=1;
        echo "1 $epsilon1 $eval" >> GRAFICAS/epsilon_nsga2.out
    fi
    #COMPROBAMOS SI EL PROMEDIO ES MENOR QUE EPSILON_1 PARA HIBRIDO
    if [ $(awk 'BEGIN{if('$prom_deltap_H'<='$epsilon1'){print "1"}else {print 0}}') -eq 1 ] && [ "$flag_h_e1" -eq 0 ];
    then
	flag_h_e1=1;
        echo "1 $epsilon1 $eval" >> GRAFICAS/epsilon_hibrido.out
	echo "1 $epsilon1 $mean_nsga2 $mean_grad" >> porc_eval.out
    fi

    #COMPROBAMOS SI EL PROMEDIO ES MENOR QUE EPSILON_2 PARA NSGA-II
    if [ $(awk 'BEGIN{if('$prom_deltap_O'<='$epsilon2'){print "1"}else {print 0}}') -eq 1 ] && [ "$flag_o_e2" -eq 0 ];
    then
	flag_o_e2=1;
        echo "2 $epsilon2 $eval" >> GRAFICAS/epsilon_nsga2.out
    fi
    #COMPROBAMOS SI EL PROMEDIO ES MENOR QUE EPSILON_2 PARA HIBRIDO
    if [ $(awk 'BEGIN{if('$prom_deltap_H'<='$epsilon2'){print "1"}else {print 0}}') -eq 1 ] && [ "$flag_h_e2" -eq 0 ];
    then
	flag_h_e2=1;
        echo "2 $epsilon2 $eval" >> GRAFICAS/epsilon_hibrido.out
	echo "2 $epsilon2 $mean_nsga2 $mean_grad" >> porc_eval.out
    fi

    #COMPROBAMOS SI EL PROMEDIO ES MENOR QUE EPSILON_3 PARA NSGA-II
    if [ $(awk 'BEGIN{if('$prom_deltap_O'<='$epsilon3'){print "1"}else {print 0}}') -eq 1 ] && [ "$flag_o_e3" -eq 0 ];
    then
	flag_o_e3=1;
        echo "3 $epsilon3 $eval" >> GRAFICAS/epsilon_nsga2.out
    fi
    #COMPROBAMOS SI EL PROMEDIO ES MENOR QUE EPSILON_3 PARA HIBRIDO
    if [ $(awk 'BEGIN{if('$prom_deltap_H'<='$epsilon3'){print "1"}else {print 0}}') -eq 1 ] && [ "$flag_h_e3" -eq 0 ];
    then
	flag_h_e3=1;
        echo "3 $epsilon3 $eval" >> GRAFICAS/epsilon_hibrido.out
	echo "3 $epsilon3 $mean_nsga2 $mean_grad" >> porc_eval.out
    fi

    #COMPROBAMOS SI EL PROMEDIO ES MENOR QUE EPSILON_4 PARA NSGA-II
    if [ $(awk 'BEGIN{if('$prom_deltap_O'<='$epsilon4'){print "1"}else {print 0}}') -eq 1 ] && [ "$flag_o_e4" -eq 0 ];
    then
	flag_o_e4=1;
        echo "4 $epsilon4 $eval" >> GRAFICAS/epsilon_nsga2.out
    fi
    #COMPROBAMOS SI EL PROMEDIO ES MENOR QUE EPSILON_4 PARA HIBRIDO
    if [ $(awk 'BEGIN{if('$prom_deltap_H'<='$epsilon4'){print "1"}else {print 0}}') -eq 1 ] && [ "$flag_h_e4" -eq 0 ];
    then
	flag_h_e4=1;
        echo "4 $epsilon4 $eval" >> GRAFICAS/epsilon_hibrido.out
	echo "4 $epsilon4 $mean_nsga2 $mean_grad" >> porc_eval.out
    fi

    #COMPROBAMOS SI EL PROMEDIO ES MENOR QUE EPSILON_5 PARA NSGA-II
    if [ $(awk 'BEGIN{if('$prom_deltap_O'<='$epsilon5'){print "1"}else {print 0}}') -eq 1 ] && [ "$flag_o_e5" -eq 0 ];
    then
	flag_o_e5=1;
        echo "5 $epsilon5 $eval" >> GRAFICAS/epsilon_nsga2.out
    fi
    #COMPROBAMOS SI EL PROMEDIO ES MENOR QUE EPSILON_5 PARA HIBRIDO
    if [ $(awk 'BEGIN{if('$prom_deltap_H'<='$epsilon5'){print "1"}else {print 0}}') -eq 1 ] && [ "$flag_h_e5" -eq 0 ];
    then
	flag_h_e5=1;
        echo "5 $epsilon5 $eval" >> GRAFICAS/epsilon_hibrido.out
	echo "5 $epsilon5 $mean_nsga2 $mean_grad" >> porc_eval.out
    fi

    #COMPROBAMOS SI EL PROMEDIO ES MENOR QUE EPSILON_6 PARA NSGA-II
    if [ $(awk 'BEGIN{if('$prom_deltap_O'<='$epsilon6'){print "1"}else {print 0}}') -eq 1 ] && [ "$flag_o_e6" -eq 0 ];
    then
	flag_o_e6=1;
        echo "6 $epsilon6 $eval" >> GRAFICAS/epsilon_nsga2.out
    fi
    #COMPROBAMOS SI EL PROMEDIO ES MENOR QUE EPSILON_6 PARA HIBRIDO
    if [ $(awk 'BEGIN{if('$prom_deltap_H'<='$epsilon6'){print "1"}else {print 0}}') -eq 1 ] && [ "$flag_h_e6" -eq 0 ];
    then
	flag_h_e6=1;
        echo "6 $epsilon6 $eval" >> GRAFICAS/epsilon_hibrido.out
	echo "6 $epsilon6 $mean_nsga2 $mean_grad" >> porc_eval.out
    fi

    j=$(bc <<< ''$j'+'$aumento'')
done

cd GRAFICAS/
gnuplot comparativo.gnu
cd ..

rm -f all_pop.out best_pop.out borrar.txt eval.out final_pop.out frente_aprox.dat initial_pop.out params.out qs.out

#REANUDAMOS SINCRONIZACION DE DROPBOX
dropbox start
