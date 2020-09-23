#include <stdio.h> 
#include <stdlib.h> 

/*VARIABLES GLOBALES*/
/*Arreglos bidimencionales*/
int no_eje;
double *GD, *IGD, *DELTAP, *temp;

/*FUNCIONES*/
void qs(int, int);

int main(int argc, char **argv)
{
    FILE *fp_gd_box, *fp_igd_box, *fp_deltap_box, *fp_qs, *fp_medidas;

    int i, j, k, q1, q2, q3;
    double sum_gd, sum_igd, sum_deltap;
    char medidas[30];
    
    if (argc<2)
    {
        printf("\n Escribe './qs no_eje' donde no_eje es # de ejecuciones del NSGA-II \n");
        exit(1);
    }
    no_eje = (int)atof(argv[1]);
    scanf("%s", medidas);

    fp_medidas=fopen(medidas, "r");
    if(fp_medidas==NULL)
    {
        fputs ("ERROR, no se pudo abrir el archivo de datos\n",stderr);
        return 0;
    }
    
    /*Reservamos memoria para arreglos bidimencionales*/
    GD = (double*)malloc(no_eje*sizeof(double));
    IGD = (double*)malloc(no_eje*sizeof(double));
    DELTAP = (double*)malloc(no_eje*sizeof(double));
    temp = (double*)malloc(no_eje*sizeof(double));

    for(i=0; i<no_eje ;i++)
    {
        GD[i]=0.0;
        IGD[i]=0.0;
	DELTAP[i]=0.0;
        temp[i]=0.0;
    }

    fp_gd_box=fopen("gd_box.out", "w");
    fp_igd_box=fopen("igd_box.out", "w");
    fp_deltap_box=fopen("deltap_box.out", "w");
    fp_qs=fopen("qs.out", "w");
    
    if(fp_gd_box==NULL)
    {
        fputs ("ERROR, no se pudo abrir el archivo 'gd_box.out'\n",stderr);
        return 0;
    }
    if(fp_igd_box==NULL)
    {
        fputs ("ERROR, no se pudo abrir el archivo 'igd_box.out'\n",stderr);
        return 0;
    }
    if(fp_deltap_box==NULL)
    {
	fputs ("ERROR, no se pudo abrir el archivo 'deltap_box.out'\n",stderr);
        return 0;
    }
    if(fp_qs==NULL)
    {
        fputs ("ERROR, no se pudo abrir el archivo 'qs.out'\n",stderr);
        return 0;
    }
    
    for(i=0; i<no_eje ;i++)
        fscanf(fp_medidas, "%lf %lf %lf", &GD[i], &IGD[i], &DELTAP[i]);
    
    /*ORDENAMOS GD*/
    for(i=0; i<no_eje ;i++)
        temp[i]=GD[i];
    qs(0, no_eje-1);
    for(i=0; i<no_eje ;i++)
        GD[i]=temp[i];

    /*ORDENAMOS IGD*/
    for(i=0; i<no_eje ;i++)
        temp[i]=IGD[i];
    qs(0, no_eje-1);
    for(i=0; i<no_eje ;i++)
        IGD[i]=temp[i];

    /*ORDENAMOS DELTAP*/
    for(i=0; i<no_eje ;i++)
        temp[i]=DELTAP[i];
    qs(0, no_eje-1);
    for(i=0; i<no_eje ;i++)
        DELTAP[i]=temp[i];

    /*Imprime las medidas ordenadas*/
    fprintf(fp_qs, "#GD IGD DELTAp\n");
    for(i=0; i<no_eje ;i++)
        fprintf(fp_qs, "%e %e %e\n", GD[i], IGD[i], DELTAP[i]);

    /*Creamos las variables para las graficas de caja*/
    fprintf(fp_gd_box, "#MIN MAX Q1 Q2 Q3 MU\n");
    fprintf(fp_igd_box, "#MIN MAX Q1 Q2 Q3 MU\n");
    
    q1=(no_eje/4)+1;
    q2=(no_eje/2)+1;
    q3=(no_eje*3/4)+1;
    
    sum_gd=0;
    sum_igd=0;
    sum_deltap=0;
    for(j=0; j<no_eje; j++)
    {
        sum_gd=sum_gd+GD[j];
        sum_igd=sum_igd+IGD[j];
	sum_deltap=sum_deltap+DELTAP[j];
    }
    sum_gd=sum_gd/no_eje;
    sum_igd=sum_igd/no_eje;
    sum_deltap=sum_deltap/no_eje;
        
    fprintf(fp_gd_box, "%e %e %e %e %e %e\n", GD[0], GD[no_eje-1], GD[q1], GD[q2], GD[q3], sum_gd);
    fprintf(fp_igd_box, "%e %e %e %e %e %e\n", IGD[0], IGD[no_eje-1], IGD[q1], IGD[q2], IGD[q3], sum_igd);
    fprintf(fp_deltap_box, "%e %e %e %e %e %e\n", DELTAP[0], DELTAP[no_eje-1], DELTAP[q1], DELTAP[q2], DELTAP[q3], sum_deltap);
    
    fclose (fp_qs); fclose (fp_gd_box); fclose(fp_deltap_box); fclose (fp_igd_box); fclose (fp_medidas);
    
    free (GD);
    free (IGD);
    free (DELTAP);
    
    return 0;
}

void qs(int limite_izq, int limite_der)
{
    int izq, der;
    double temporal, pivote;

    izq=limite_izq;
    der = limite_der;
    pivote = temp[(izq+der)/2];

    do{
        while(temp[izq]<pivote && izq<limite_der)
            izq++;
        while(pivote<temp[der] && der > limite_izq)
            der--;
        if(izq <=der)
        {
            temporal= temp[izq];
            temp[izq]=temp[der];
            temp[der]=temporal;
            izq++;
            der--;

        }

    }while(izq<=der);
    
    if(limite_izq<der)
        qs(limite_izq, der);

    if(limite_der>izq)
        qs(izq, limite_der);
}
