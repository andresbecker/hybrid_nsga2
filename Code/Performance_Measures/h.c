#include <stdio.h> 
#include <stdlib.h> 
#include <math.h>

# define P 2
# define R 0.5

/*FUNCIONES*/
double distancia_conjuntos(double **, double **, int, int, int);
double distancia_euclidiana(double **, double **, int, int);

int main()
{
    FILE *fp_pareto, *fp_datos;
    
    /*Arreglos bidimencionales*/
    double **pareto, **datos, gd, igd, gd_p, igd_p, delta_p;
    int i, j, k, pareto_size, pop_size, delta=0;
    char frente_pareto[30], frente_aprox[30];
    
    scanf("%s", frente_aprox);
    scanf("%d", &pop_size);
    scanf("%s", frente_pareto);
    scanf("%d", &pareto_size);

    fp_pareto = fopen (frente_pareto, "r");
    if (fp_pareto==NULL)
    {
	fputs ("ERROR, no se pudo abrir el archivo del frente de pareto\n",stderr);
	return 0;
    }
    fp_datos = fopen (frente_aprox, "r");
    if (fp_datos==NULL)
    {
	fputs ("ERROR, no se pudo abrir el archivo del frente aprox\n",stderr);
	return 0;
    }
    
    /*Reservamos memoria para arreglos bidimencionales*/
    pareto = (double **)malloc(pareto_size*sizeof(double *));
    for(i=0; i<pareto_size; i++)
        pareto[i]=(double *)malloc(2*sizeof(double));

    datos = (double **)malloc(pop_size*sizeof(double *));
    for(i=0; i<pop_size; i++)
        datos[i]=(double *)malloc(2*sizeof(double));
    
    //Inicializamos lor arreglos
    for(i=0; i<pareto_size; i++)
        for(j=0; j<2; j++)
            pareto[i][j]=0.0;

    for(i=0; i<pop_size; i++)
       for(j=0; j<2; j++)
            datos[i][j]=0.0;

    //COPIAMOS DEL ARCHIVO A LOS ARREGLOS
    for(i=0; i<pareto_size; i++)
        fscanf(fp_pareto, "%lf %lf", &pareto[i][0], &pareto[i][1]);

    for(i=0; i<pop_size; i++)
        fscanf(fp_datos, "%lf %lf", &datos[i][0], &datos[i][1]);

    //CALCULAMOS DISTANCIAS
    gd=distancia_conjuntos(datos, pareto, pop_size, pareto_size, delta);
    igd=distancia_conjuntos(pareto, datos, pareto_size, pop_size, delta);
    delta=1;
    gd_p=distancia_conjuntos(datos, pareto, pop_size, pareto_size, delta);
    igd_p=distancia_conjuntos(pareto, datos, pareto_size, pop_size, delta);

    //ENCONTRAMOS DELTA_P 
    delta_p=gd_p;
    if(delta_p<igd_p)
	delta_p=igd_p;
    
    //IMPRIMIMOS
    printf("%e %e %e\n", gd, igd, delta_p);
    
    
    fclose (fp_pareto); fclose (fp_datos);
    
    for(i=0; i<pareto_size; i++)
    {
        free(pareto[i]);
    }
    free(pareto);
    for(i=0; i<pop_size; i++)
    {
        free(datos[i]);
    }
    free(datos);
    
    return 0;
}

double distancia_conjuntos(double **x, double **y, int n, int m, int delta)
{
    int i, j;
    double d_puntos, inf, sum, d_conjunto;
    
    sum=0; 
    for(j=0; j<n; j++)
    {
        inf=distancia_euclidiana(x,y,0,j);
        for(i=1; i<m; i++)
        {
            d_puntos=distancia_euclidiana(x,y,i,j);
            if(d_puntos<inf)
                inf=d_puntos;
        }
        sum=sum+pow(inf,P);
    }
    if(delta==0)
        d_conjunto=pow(sum,R)/n;
    if(delta==1)
        d_conjunto=pow(sum/n,R);
    
    return d_conjunto;
}

double distancia_euclidiana(double **x, double **y, int ii, int jj)
{
    double d;
    int i;
    
    d=0;
    for(i=0; i<2; i++)
        d=d+pow(x[jj][i]-y[ii][i],2);
    d=sqrt(d);
    
    return d;
}
