///////////////Librerias//////////////////     
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "PSO.h"
///////////Ficheros Globales///////////////
FILE *Fichero_Grafica;
FILE *Fichero_Reporte;
FILE *Fichero_Population;
FILE *fp_aux;
FILE *Track_Particles;
///////////////////////////////////////////
int gpopsize=0;
double varibles_interval[21]={-0.5, -0.45, -0.4, -0.35, -0.3, -0.25, -0.2, -0.15, -0.1, -0.05, 0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5};
char nombreArchivo[100];
char comando[300];
#define Num_Variables 3 
///////////////////////////////////////////
int main(int argc , char *argv[])
{
	int popsize,max_generations,g=0,n=0;
	float g_pond,p_pond,max_inertia,min_inertia,inertia,seed;
	char Reporte[100];
	if(argc != 9)
	{
		printf("Faltan argumentos\n");
		printf("Sintaxis: ./prog , #Particles , Generations , Global_P , Personal_P , max_inertia , min_inertia, seed(0,1), output_file\n");
	}
	else
	{
		popsize=atoi(argv[1]);
		max_generations=atoi(argv[2]);
		g_pond=atof(argv[3]);
		p_pond=atof(argv[4]);
		max_inertia=atof(argv[5]);
		min_inertia=atof(argv[6]);
		seed=atof(argv[7]);
		strcpy(Reporte,argv[8]);		
		strcat(Reporte,".txt");
		Fichero_Grafica = fopen("grafica.txt", "w" );
		Fichero_Reporte = fopen(Reporte,"w");
		Allocate_Memory(popsize);
		Get_Previous_Individuals();
		printf("Initialize\n");			
		bestp.fitness = 0;
		bestp.frecuencia = 0;
 		bestp.amplitud = 0;
		bestp.generation = 0;
		/* Initialize the populations and Initial Report */		
		randomize(seed);			
		Initialize_Population(popsize);		
		Initial_Report(popsize,max_generations,p_pond,g_pond,max_inertia,min_inertia,seed);
		Track_Particles = fopen("First_Track.txt", "w" );
			Track_Particles_Function(oldpop,popsize);
		fclose(Track_Particles);
		show_pop(popsize);
		//PSO Algorithm
		while(g<max_generations)
		{
			statistics(oldpop,g,popsize,Num_Variables);
			//Update best personal and global position
			Update_Bp_Gg(popsize,g);
			//Linear decreasing Inertia
			inertia=(max_inertia-min_inertia)*((max_generations-g)/max_generations)+0.1;
			printf("Inertia:%f\n",inertia);
			Update_Position(popsize,p_pond,g_pond,inertia);
			g++;			
		}		
		//////////////////////////////////////////////////////
		statistics(oldpop,g,popsize,Num_Variables);	
		Track_Particles = fopen("Final_Track.txt", "w" );
			Track_Particles_Function(oldpop,popsize);
		fclose(Track_Particles);
		show_pop(popsize);					
		fclose(Fichero_Grafica);
		fclose(Fichero_Reporte);		
		Free_All();
	}
}

void Allocate_Memory(int popsize)
{
	printf("Allocating Memory\n");
	unsigned numbytes,numbytes2;	
	/* Allocate memory for old and new populations of particles */
	numbytes = popsize*sizeof(struct particle);
	numbytes2=2000*sizeof(struct Population);	
	if ((oldpop = (struct particle *) malloc(numbytes)) == NULL)
	{
		nomemory("Old Population");
	}	
	if ((globalpop = (struct Population *) malloc(numbytes2)) == NULL)
	{
		nomemory("Global Population");
	}
}

void nomemory(char *string)
{
	printf("ERROR!! --> malloc: out of memory making %s\n",string);
	exit(1);
}

void Free_All()
{
	printf("Free Memory\n");
	free(oldpop);	
	free(globalpop);
}

void Get_Previous_Individuals()
{
	FILE *file = fopen("Poblacion.dat", "r");
	struct Population aux;
	printf("Leyendo archivo de Poblacion Guardada\n");   
	while(fread(&aux, sizeof(aux), 1, file))
	{		
    	globalpop[gpopsize].x1=aux.x1;
		globalpop[gpopsize].x2=aux.x2;
		globalpop[gpopsize].x3=aux.x3;
		globalpop[gpopsize].fitness = aux.fitness;
		globalpop[gpopsize].frecuencia = aux.frecuencia;
		globalpop[gpopsize].amplitud = aux.amplitud;
		gpopsize++;
		printf("Global Pop size is now:%d\n",gpopsize );		
		/*printf("x1=%lf\n",ind.x1 );
		printf("x2=%lf\n",ind.x2 );
		printf("x3=%lf\n",ind.x3 );
		printf("VFuncion=%lf\n",ind.valorFuncion );
		printf("Frecuencia=%lf\n",ind.frecuencia );
		printf("Amplitud=%lf\n\n",ind.amplitud );*/
	}		
    fclose(file); 

}

void Initialize_Population(int popsize)
{
	int i,k;
	for(i=0;i<popsize;i++)
	{
		oldpop[i].x1=rndreal(-0.5 ,0.5);
		oldpop[i].x2=rndreal(-0.5 ,0.5);
		oldpop[i].x3=rndreal(-0.5 ,0.5);
		oldpop[i].v1=rndreal(0.1 ,0.2);
		oldpop[i].v2=rndreal(0.1 ,0.2);
		oldpop[i].v3=rndreal(0.1 ,0.2);;
		Evaluate_Individual(&(oldpop[i]));
		oldpop[i].best_x1=oldpop[i].x1;
		oldpop[i].best_x2=oldpop[i].x2;
		oldpop[i].best_x3=oldpop[i].x3;
		oldpop[i].best_p_f=oldpop[i].fitness;
	}
}

void Evaluate_Individual(struct particle *ind)
{
	int indice =-1;
	double x1,x2,x3;	
	ind->fitness=0.0;
	char c[10];
	x1=round_number(ind->x1);
	ind->x1=x1;
	x2=round_number(ind->x2);
	ind->x2=x2;
	x3=round_number(ind->x3);
	ind->x3=x3;
	if(gpopsize>0)
	{
		indice = buscar(ind);
	}
	/*Script que usa FDTD para el decaimiento de la energia*/
	if(indice==-1)
	{
		//meep shift1=0.1 shift2=0.4 shift3=0.5 funcion/l3defect.ctl | tee a.out
		//mpirun.openmpi -np 2 meep-mpi-default shift1=0.1 shift2=0.4 shift3=0.5 funcion/l3defect.ctl | tee a.out
		//strcpy(comando,"mpirun.openmpi -np 3 meep-mpi-default shift1=");
		strcpy(comando,"meep shift1=");
		sprintf(c, "%g", x1);strcat(comando,c);
		strcat(comando," shift2=");
		sprintf(c, "%g", x2);strcat(comando,c);
		strcat(comando," shift3=");
		sprintf(c, "%g", x3);strcat(comando,c);
		strcat(comando," funcion/l3defect.ctl > funcion/a.out");
		printf("\n\n%s\n",comando);
		system(comando);	

		//PARA EL FACTOR FRECUENCIA AL FINAL ES -f 1
		printf("Obteniendo Frecuencia\n");	
		strcpy(comando,"grep harminv funcion/a.out | cut -d "); 
		strcat(comando,",");
		strcat(comando," -f 2,4,5 | sort -r -g -k 3 -t , | sed -n 1,1p | cut -d ");
		strcat(comando,",");
		strcat(comando," -f 1");
		fp_aux = popen(comando, "r");
			fscanf(fp_aux, "%lf", &ind->frecuencia);
			printf("Frecuencia:%lf \n",ind->frecuencia );
		pclose(fp_aux);
	
		//PARA EL FACTOR Q AL FINAL ES -f 2  lo guardaremos en el ind->valorfunción 
		printf("Obteniendo Q\n");	
		strcpy(comando,"grep harminv funcion/a.out | cut -d ");
		strcat(comando,",");
		strcat(comando," -f 2,4,5 | sort -r -g -k 3 -t , | sed -n 1,1p | cut -d ");
		strcat(comando,",");
		strcat(comando," -f 2");
		fp_aux = popen(comando, "r");
			fscanf(fp_aux, "%lf", &ind->fitness);
			printf("Fitness:%lf \n",ind->fitness );
		pclose(fp_aux);
     
		//PARA LA AMPLITUD AL FINAL ES -f 3
		//grep harminv a.out | cut -d "," -f 2,4,5 | sort -r -g -k 3 -t , | sed -n 1,1p | cut -d "," -f 2
		printf("Obteniendo Amplitud\n");	
		strcpy(comando,"grep harminv funcion/a.out | cut -d "); 
		strcat(comando,",");
		strcat(comando," -f 2,4,5 | sort -r -g -k 3 -t , | sed -n 1,1p | cut -d ");
		strcat(comando,",");
		strcat(comando," -f 3");
		fp_aux = popen(comando, "r"); 
			fscanf(fp_aux, "%lf", &ind->amplitud);
			printf("Amplitud:%lf \n\n",ind->amplitud );
		pclose(fp_aux);
		//En este caso el valor de la función el fitness es igual y es dado por el factor Q		
		add_individual(ind);
	}
	else
	{
		printf("Individuo ya evaluado\n");		
		ind->fitness = globalpop[indice].fitness;
		ind->frecuencia = globalpop[indice].frecuencia;
		ind->amplitud = globalpop[indice].amplitud;
		printf("Frecuencia:%lf \n",ind->frecuencia );
		printf("Fitness:%lf \n",ind->fitness );
		printf("Amplitud:%lf \n\n",ind->amplitud );		
	}
}

void randomize(float seed)
{
    int j1;
    for(j1=0; j1<=54; j1++)
    {
      oldrand[j1] = 0.0;
    }
    jrand=0;
    warmup_random(seed);
}

void warmup_random(float random_seed)
{
    int j1, ii;
    double new_random, prev_random;
    oldrand[54] = random_seed;
    new_random = 0.000000001;
    prev_random = random_seed;
    for(j1 = 1 ; j1 <= 54; j1++)
    {
		ii = (21*j1)%54;
		oldrand[ii] = new_random;
		new_random = prev_random-new_random;
		if(new_random<0.0) new_random = new_random + 1.0;
		prev_random = oldrand[ii];
    }
    advance_random();
    advance_random();
    advance_random();
    jrand = 0;
}

void advance_random()
{
    int j1;
    double new_random;
    for(j1 = 0; j1 < 24; j1++)
    {
		new_random = oldrand[j1] - oldrand[j1+31];
		if(new_random < 0.0) new_random = new_random + 1.0;
			oldrand[j1] = new_random;		
    }
    for(j1 = 24; j1 < 55; j1++)
    {
		new_random = oldrand [j1] - oldrand [j1-24];
		if(new_random < 0.0) new_random = new_random + 1.0;
			oldrand[j1] = new_random;
    }
}

float rndreal(float lo ,float hi)
{
    return((randomperc() * (hi - lo)) + lo);
}

float randomperc()
{
    jrand++;
    if(jrand >= 55)
    {
		jrand = 1;
		advance_random();
    }
    return((float) oldrand[jrand]);
}

void add_individual(struct particle * ind)
{
	if(gpopsize<(2000))
	{
		//Agregar Individuo a Arreglo en Memoria
		globalpop[gpopsize].x1=ind->x1;
		globalpop[gpopsize].x2=ind->x2;
		globalpop[gpopsize].x3=ind->x3;		
		globalpop[gpopsize].fitness = ind->fitness;
		globalpop[gpopsize].frecuencia = ind->frecuencia;
		globalpop[gpopsize].amplitud = ind->amplitud;		
		//Agregar Individuo a Archivo
		printf("Agregare:\n");
		printf("x1:%f,x2:%f,x3:%f\n",ind->x1,ind->x2,ind->x3 );
		FILE *file = fopen("Poblacion.dat", "a");
		struct Population individual;
		individual.x1=ind->x1;
		individual.x2=ind->x2;
		individual.x3=ind->x3;
		individual.frecuencia=ind->frecuencia;
		individual.fitness=ind->fitness;
		individual.amplitud=ind->amplitud;    
			fwrite(&individual, sizeof(individual), 1, file);	
		printf("Lo agregue\n");		
    	fclose(file);  
    	gpopsize++;
	}
}

int buscar (struct particle * ind)
{
	for (int i=0; i < gpopsize; i++) 
	{
		if ((globalpop[i].x1==ind->x1)&&(globalpop[i].x2==ind->x2)&&(globalpop[i].x3==ind->x3))
			{
				return i;
			}
	}
	return -1;
}

void show_pop(int popsize)
{
	int i;
	for(i=0;i<popsize;i++)
	{
		printf("I%i\n",i );
		printf("[X1:%f] [X2:%f] [X3:%f]	\n",oldpop[i].x1,oldpop[i].x2,oldpop[i].x3);
		printf("[V1:%f] [V2:%f] [V3:%f]	\n",oldpop[i].v1,oldpop[i].v2,oldpop[i].v3);
		printf("Q:%f\n",oldpop[i].fitness );
		printf("Amplitud:%f\n",oldpop[i].amplitud );
		printf("Frecuencia:%f\n",oldpop[i].frecuencia );
		printf("\n");
	}
}

void Initial_Report(int popsize,int gmax,float personal,float global,float i_max,float i_min,float seed)
{
	printf("\n >>>>> Parameters used with PSO <<<<<\n");
 	fprintf( Fichero_Reporte ,">>>>> Parameters used with PSO <<<<< \n");
	printf("Total population size\t=\t%d\n",popsize);
	fprintf( Fichero_Reporte ,"Total population size\t=\t%d\n",popsize);
	printf("Max Generations\t\t=\t%d\n",gmax);
	fprintf( Fichero_Reporte ," Max Generations\t\t=\t%d\n",gmax);
	printf("Personal Ponderation\t=\t%f\n",personal);
	fprintf( Fichero_Reporte ,"Personal Ponderation\t=\t%f\n",personal);
	printf("Global Ponderation\t=\t%f\n",global);
	fprintf( Fichero_Reporte ,"Global Ponderation\t=\t%f\n",global);
	printf("Maximum Inertia\t\t=\t%f\n",i_max);
	fprintf( Fichero_Reporte ,"Maximum Inertia\t\t=\t%f\n",i_max);
	printf("Minimum Inertia\t=\t%f\n",i_min);
	fprintf( Fichero_Reporte ,"Minimum Inertia\t\t\t=\t%f\n",i_min);
	printf("Seed\t=\t%f\n",seed);
	fprintf( Fichero_Reporte ,"Seed\t\t\t=\t%f\n",seed);	
	printf("\n\n");
	fprintf( Fichero_Reporte,"\n\n");	
}

void Update_Bp_Gg(int popsize,int gen)
{
	int i;
	bestp_neighborhood.fitness=oldpop[0].fitness;
	bestp_neighborhood.generation=gen;
	for(i=0;i<popsize;i++)
	{
		//Update Best Personal
		if(oldpop[i].fitness>oldpop[i].best_p_f)
		{
			oldpop[i].best_p_f=oldpop[i].fitness;
			oldpop[i].best_x1=oldpop[i].x1;
			oldpop[i].best_x2=oldpop[i].x2;
			oldpop[i].best_x3=oldpop[i].x3;	
		}
		//Update Best Global of Neighborhood
		if(oldpop[i].fitness>bestp_neighborhood.fitness)
		{
			bestp_neighborhood.fitness=oldpop[i].fitness;
			bestp_neighborhood.x1=oldpop[i].x1;
			bestp_neighborhood.x2=oldpop[i].x2;
			bestp_neighborhood.x3=oldpop[i].x3;
			bestp_neighborhood.amplitud=oldpop[i].amplitud;	
			bestp_neighborhood.frecuencia=oldpop[i].frecuencia;
		}
	}
}

void Update_Position(int popsize,float p_pond,float g_pond,float inertia)
{
	int i,n;
	for(i=0;i<popsize;i++)
	{
		oldpop[i].v1=(inertia * oldpop[i].v1) + 
					 (p_pond)*(oldpop[i].best_x1-oldpop[i].x1) +
					 (g_pond)*(bestp_neighborhood.x1-oldpop[i].x1);

		oldpop[i].x1=oldpop[i].x1+oldpop[i].v1;

		oldpop[i].v2=(inertia * oldpop[i].v2) + 
					 (p_pond)*(oldpop[i].best_x2-oldpop[i].x2) +
					 (g_pond)*(bestp_neighborhood.x2-oldpop[i].x2);

		oldpop[i].x2=oldpop[i].x2+oldpop[i].v2;

		oldpop[i].v3=(inertia * oldpop[i].v3) + 
					 (p_pond)*(oldpop[i].best_x3-oldpop[i].x3) +
					 (g_pond)*(bestp_neighborhood.x3-oldpop[i].x3);

		oldpop[i].x3=oldpop[i].x3+oldpop[i].v3;

		Evaluate_Individual(&(oldpop[i]));	
	}
}

double round_number(double number)
{
	int i=0;
	for (i=0;i<21;i++)
		if (number<=varibles_interval[i])
			break;
	return varibles_interval[i];
}

void statistics(struct particle *pop,int gen,int popsize,int variables)
{
	int i,j;
	int mejor=0,peor=0;
	float sumfitness = 0.0;	
	float min = pop[0].fitness;
	float max = pop[0].fitness;
	float mean_gen=0.0;
	float ecm=0.0;
	/* Loop for max, min, sumfitness */
	for (j=0; j < popsize; j++)
	{	
		sumfitness = sumfitness + pop[j].fitness;            
		if (pop[j].fitness > max)
		{ 
			/* New maximum  fitness*/
			max = pop[j].fitness;			
			/* New maximum  fitness index*/
		    mejor=j;
		}		
		if (pop[j].fitness < min)
		{
			/* New minimim fitness */
			min = pop[j].fitness;
			/* New minimum index*/
			peor=j;
		}				
	}
	for(j=0;j<popsize;j++)
	{

	}
	//Update Best Ever
	if(pop[mejor].fitness>bestp.fitness)
	{
		bestp.fitness=pop[mejor].fitness;
		bestp.x1=pop[mejor].x1;
		bestp.x2=pop[mejor].x2;
		bestp.x3=pop[mejor].x3;
		bestp.amplitud=pop[mejor].amplitud;	
		bestp.frecuencia=pop[mejor].frecuencia;
		bestp.generation=gen;
	}
	mean_gen=sumfitness/popsize;
	fprintf(Fichero_Grafica,"%d\t%.10f\n",gen,bestp.fitness);
	ecm=ECM(oldpop,popsize,variables);
	report(min,peor,max,mejor,mean_gen,gen,ecm);	
}

void report(float min,int peor,float max,int mejor,float mean_gen,int gen,float ecm)
{	
	printf("\nEstadisticas de generacion #%d\n\n",gen);
	printf("Q_max:%f\n",max);
	printf("Q_min:%f\n",min);
	printf("Q_mean:%f\n",mean_gen);
	printf("Best Individual [x1,x2,x3]=[ %f , %f , %f ]\n",oldpop[mejor].x1,oldpop[mejor].x2,oldpop[mejor].x3);	
	printf("Amplitud:%f\n",oldpop[mejor].amplitud );
	printf("Frecuencia:%f\n",oldpop[mejor].frecuencia );
	printf("ECM:%f\n",ecm);
	printf("Best Individual Ever in gen #%d [x1,x2,x3]=[ %f , %f , %f ]\n",bestp.generation,bestp.x1,bestp.x2,bestp.x3);
	printf("\n");
	fprintf(Fichero_Reporte,"\nEstadisticas de generacion #%d\n\n",gen);
	fprintf(Fichero_Reporte,"Q_max:%f\n",max);
	fprintf(Fichero_Reporte,"Q_min:%f\n",min);
	fprintf(Fichero_Reporte,"Q_mean:%f\n",mean_gen);
	fprintf(Fichero_Reporte,"Best Individual [x1,x2,x3]=[ %f , %f , %f ]\n",oldpop[mejor].x1,oldpop[mejor].x2,oldpop[mejor].x3);	
	fprintf(Fichero_Reporte,"Amplitud:%f\n",oldpop[mejor].amplitud );
	fprintf(Fichero_Reporte,"Frecuencia:%f\n",oldpop[mejor].frecuencia );
	fprintf(Fichero_Reporte,"ECM:%f\n",ecm);
	fprintf(Fichero_Reporte,"Best Individual Ever in gen #%d [x1,x2,x3]=[ %f , %f , %f ]\n",bestp.generation,bestp.x1,bestp.x2,bestp.x3);
	fprintf(Fichero_Reporte,"\n");
}

void Track_Particles_Function(struct particle *pop,int popsize)
{
	int i;
	for(i=0;i<popsize;i++)
	{
		fprintf(Track_Particles,"%f\t%f\t%f\n",oldpop[i].x1,oldpop[i].x2,oldpop[i].x3);
	}
}

double ECM(struct particle *pop,int popsize,int variables)
{
	int j,k;
	double ecm=0.0,ecm_sub=0.0;
	for(j=0;j<popsize;j++)
	{
		ecm_sub=pow((pop[j].x1-bestp.x1),2);			
		ecm_sub=ecm_sub+pow((pop[j].x2-bestp.x2),2);
		ecm_sub=ecm_sub+pow((pop[j].x2-bestp.x2),2);			
		ecm_sub=sqrt(ecm_sub);		
		ecm=ecm+ecm_sub;
		ecm_sub=0.0;
	}
	ecm=ecm/variables;
	return ecm;
}
