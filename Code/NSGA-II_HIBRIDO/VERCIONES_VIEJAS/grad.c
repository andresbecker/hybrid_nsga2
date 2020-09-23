# include <stdio.h>
# include <stdlib.h>
# include <math.h>

# include "global.h"
# include "rand.h"


void metodo_grad(population *popu)
{
    individual *itemp;
    int i, j, k, kk, iter, pop;
    double *norm, alpha, c_parada, temp;
    FILE *fp_alfa;

    fp_alfa=fopen("alfa.out", "ab");
    if(fp_alfa==NULL)
    {
	fputs ("ERROR, no se pudo abrir el archivo alfa.out\n",stderr);
	return 0;
    }
    
    alpha=0;
    
    /*RESERVAMOS MEMORIA*/
    norm=(double *)malloc(nobj*sizeof(double));
    direccion=(double *)malloc(nreal*sizeof(double));
    gradiente = (double**)malloc(nobj*sizeof(double));
    for(i=0; i<nobj; i++)
    {
        gradiente[i]=(double*)malloc(nreal*sizeof(double));
    }
 
    /*COMIENZA PROCESO PARA CALCULAR NUEVOS PUNTOS*/
 
    for(pop=0; pop<popsize; pop++)
    {
        iter=0;
        itemp=&(popu->ind[pop]);
        if(itemp->rank==1)
        {
            printf("\nIterendo Individuo %d\n", pop+1);
	    while(iter<maxiter)
	    {
                /*INICIALIZAMOS ARREGLOS*/
                for(i=0; i<nobj; i++)
                {
                    for(j=0; j<nreal; j++)
                    {
                        gradiente[i][j]=0.0;
                    }
                }
                for(k=0; k<nreal;k++)
                {
                    direccion[k]=0;
                }
                /*CALCULAMOS GRADIENTE DE CADA Fi*/
                grad(&(popu->ind[pop]));
                /*CALCULAMOS LA NORMA DEL GRADIENTE DE CADA Fi*/
               	for(k=0; k<nobj; k++)
                {
                    norm[k]=norma(&(gradiente[k][0]));
                }
		/*VEMOS SI SE CUMPLE CRITERIO PARA MEJORA LOCAL*/
     		c_parada=0;
                for(k=0; k<nreal; k++)
                {
               	    temp=1;
                    for(kk=0; kk<nobj; kk++)
                    {
                        temp=temp*gradiente[kk][k]/norm[kk];
                    }
                    c_parada=c_parada+temp;
                }
		printf("\n%lf>=%f\n",c_parada,epsilon-1);
		if(c_parada >= (epsilon-1))
		{
                    /*CALCLAMOS DIRECCION
                    direccion=(graf1(x1)/|graf1(x1)|+graf2(x1)/|graf2(x1)|,...,graf1(xn)/|graf1(xn)|+graf2(xn)/|graf2(xn)|)*/
                    for(k=0; k<nreal; k++)
                    {
                        for(kk=0; kk<nobj; kk++)
                        {
                            direccion[k]=direccion[k]+(gradiente[kk][k]/norm[kk]);
                        }
                    }
                    /*CALCULAMOS ALPHA*/
                    alpha=armijo(&(popu->ind[pop]));
                    printf("\nalfa=%0.50lf", alpha);
  		    fprintf(fp_alfa,"%0.80lf\n", alpha);
                    if(alpha>0)
                    {
                        /*GENERAMOS EL NUEVO PUNTO*/
                        for(k=0; k<nreal; k++)
                        {
                            (itemp->xreal[k])=(itemp->xreal[k])-alpha*direccion[k];
                            /*CHECAMOS QUE EL X PROPUESTO ESTE DENTRO DE LOS LIMITES PARA CADA Xi*/
                            if((itemp->xreal[k])<min_realvar[k])
                                (itemp->xreal[k])=min_realvar[k]+EPS;
                            if((itemp->xreal[k])>max_realvar[k])
                                (itemp->xreal[k])=max_realvar[k]-EPS;
                        }
                        /*ACTULIZAMOS EL VALOR DE OBJ*/
                        test_problem (itemp->xreal, itemp->xbin, itemp->gene, itemp->obj, itemp->constr);
			iter++;
		    }
		    else
			iter=maxiter;
		}
		else
		    iter=maxiter;
	    }
	}
    }

    free(norm);
    free(direccion);
    for(i=0; i<nobj; i++)
    {
        free(gradiente[i]);
    }
    free(gradiente);
    fclose(fp_alfa);
}

double armijo(individual *ind)
{
    int i, j, cond_ob;
    double *der, *izq, *x, *o, *g_por_d, alfa;
    
    alfa=2*alfa_max;
            
    izq=(double *)malloc(nobj*sizeof(double));
    der=(double *)malloc(nobj*sizeof(double));
    g_por_d=(double *)malloc(nobj*sizeof(double));
    x=(double *)malloc(nreal*sizeof(double));
    o=(double *)malloc(nobj*sizeof(double));

    for(i=0; i<nobj; i++)
    {
        g_por_d[i]=0;
    }
    
    for(i=0; i<nobj; i++)
        for(j=0; j<nreal; j++)
            g_por_d[i]=g_por_d[i]+gradiente[i][j]*direccion[j];
    do
    {
        cond_ob=0;
        alfa=alfa/2;
        
        for(i=0; i<nreal; i++)
        {
            x[i]=ind->xreal[i];
            x[i]=x[i]-alfa*direccion[i];
            /*CHECAMOS QUE EL X PROPUESTO ESTE DENTRO DE LOS LIMITES PARA CADA Xi*/
            if(x[i]<min_realvar[i])
                x[i]=min_realvar[i]+EPS;
            if(x[i]>max_realvar[i])
                x[i]=max_realvar[i]-EPS; 
        }
        test_problem (x, ind->xbin, ind->gene, o, ind->constr);
        for(i=0; i<nobj; i++)
        {
            izq[i]=o[i];
        }
        for(i=0; i<nobj; i++)
        {
            der[i]=(ind->obj[i])-c1*alfa*g_por_d[i];
        }
        
        /*CORROBORA QUE ARMIJO SE CUMPLA PARA TODAS LAS Fi*/
        for(i=0; i<nobj; i++)
            if(izq[i]<=der[i])
                cond_ob++;
            
        /*CHECAMOS QUE ALFA NO HAYA DISMINUIDO DEMACIADO*/
        if(alfa==0)
        {
            printf("\n\n¡¡¡¡¡ALFA=0!!!!!!\n\n");
            return 0;
        }
            
    }while(cond_ob<nobj);
    
    return alfa;
}

void grad(individual *ind)
{
    int i, j;
    double *x, *o, original;
    
    x=(double *)malloc(nreal*sizeof(double));
    o=(double *)malloc(nobj*sizeof(double));

    for(i=0; i<nreal; i++)
    {
        x[i]=ind->xreal[i];
    }
    
    for(i=0; i<nreal; i++)
    {
        original=x[i];
        x[i]=x[i]+2*h;
        test_problem (x, ind->xbin, ind->gene, o, ind->constr);

        for(j=0; j<nobj; j++)
        {
            gradiente[j][i]=o[j];
        }
        x[i]=original;
        x[i]=x[i]+h;
        test_problem (x, ind->xbin, ind->gene, o, ind->constr);
        for(j=0; j<nobj; j++)
        {
            gradiente[j][i]=gradiente[j][i]-4*o[j];
            gradiente[j][i]=gradiente[j][i]+3*(ind->obj[j]);
            gradiente[j][i]=(-1)*gradiente[j][i]/(2*h);
        }
        x[i]=original;
    }
    return;
}

double norma(double *x)
{
    int ii;
    double n;
    n=0;
    
    for(ii=0; ii<nreal; ii++)
    {
        n=n+pow(x[ii],2);
    }
    n=sqrt(n);
    
    return n;
}
