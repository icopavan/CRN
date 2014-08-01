#include "mex.h"
/*#include <stdio.h>*/
#include <stdlib.h>
#include <math.h>
#include "global.h"
#include "rand.h"
#include "rand2.h"
#include "list.h"
#include "allocate.h"
#include "initialize.h"
#include "decode.h"
#include "eval.h"
#include "problemdef.h"
#include "dominance.h"
#include "crowddist.h"
#include "sort.h"
#include "rank.h"
#include "auxiliary.h"
#include "report.h"
#include "display.h"
#include "mutation.h"
#include "tourselect.h"
#include "crossover.h"
#include "merge.h"
#include "fillnds.h"

int nreal;
int nbin;
int nobj;
int ncon;
int popsize;
double pcross_real;
double pcross_bin;
double pmut_real;
double pmut_bin;
double eta_c;
double eta_m;
int ngen;
int nbinmut;
int nrealmut;
int nbincross;
int nrealcross;
int *nbits;
double *min_realvar;
double *max_realvar;
double *min_binvar;
double *max_binvar;
int bitlength;
int choice;
int obj1; /*objectives to be printed in case choice is 1 and angle of view*/
int obj2;
int obj3;
int angle1;
int angle2;
/* MATLAB imported variables*/
double **pu;
double **su;
double **alpha_su,**alpha_pu;
double N,B,C,pwr_max,f,pwr_pu;
int N_f,M;
double **d_pu,**d_su,**d_supu;
double **L_pu,**L_supu;
double ***L_su;
double SINRSULT,SINRPULT;


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	int i,j,k,freq;
    double *tempsu;
    double *temppu;
    double *tempalpha_su,*tempalpha_pu;
    double *input_values;
    mxArray *rhs[3];
    /* MATLAB importated variable assignment*/
    rhs[0] = mxDuplicateArray(prhs[0]);
    rhs[1] = mxDuplicateArray(prhs[1]);
    rhs[2] = mxDuplicateArray(prhs[2]);
    tempsu = mxGetPr(rhs[0]);
    temppu = mxGetPr(rhs[1]);
    input_values = mxGetPr(rhs[2]);
    N = input_values[0];
    N_f = input_values[1];
    M = input_values[2];
    B = input_values[3];
    C = input_values[4];
    pwr_max = input_values[5];
    f = input_values[6];
    pwr_pu = input_values[7];
    SINRPULT = input_values[8];
    SINRSULT = input_values[9];
    
    
   /* printf("here");
    mexPrintf("\n N_f = %lf\t%lf\t%lf\t",tempsu[0],tempsu[1],tempsu[2]);*/
    
    /* allocation for su*/
    su = (double **) malloc((2*M)*sizeof(double*));
    for(i=0;i<2*M;i++)
    {
        su[i] = (double*)malloc((2)*sizeof(double));
    }
    
    /* allocation for pu*/
    pu = (double **) malloc((2*N_f)*sizeof(double*));
    for(i=0;i<2*N_f;i++)
    {
        pu[i] = (double*)malloc((2)*sizeof(double));
    }
 
    /* reading su values*/
    for(i=0;i<2*M;i++)
    {
        for(j=0;j<2;j++)
        {
           su[i][j] = tempsu[(j*2*M)+i];
           mexPrintf("su[%d][%d] = %lf\n",i,j,su[i][j]);
        }
    }
    
    /* Reading pu values*/
    for(i=0;i<2*N_f;i++)
    {
        for(j=0;j<2;j++)
        {
           pu[i][j] = temppu[(j*2*N_f)+i];
           mexPrintf("pu[%d][%d] = %lf\n",i,j,pu[i][j]);
        }
    }

    /*Define pu*/
    d_pu = (double **)malloc(2*N_f*sizeof(double*));
    L_pu = (double **) malloc(2*N_f*sizeof(double*));
    for(i=0;i<2*N_f;i++)
    {
        d_pu[i] = (double*) malloc(2*N_f*sizeof(double));
        L_pu[i] = (double*) malloc(2*N_f*sizeof(double));
    }
    
    /* Pu distance calculation*/
    for(i=0;i<N_f;i++)
    {
        d_pu[i][N_f+i] = sqrt((pu[i][0]-pu[N_f+i][0])*(pu[i][0]-pu[N_f+i][0])+(pu[i][1]-pu[N_f+i][1])*(pu[i][1]-pu[N_f+i][1]));
        mexPrintf("d_pu[%d][%d] = %lf\n",i,N_f+i,d_pu[i][N_f+i]);
        L_pu[i][N_f+i] = (C*C)/((4*pi*(f+B*(i)))*(4*pi*(f+B*(i)))*pow((d_pu[i][N_f+i]),4));
        mexPrintf("L_pu[%d][%d] = %lf\n",i,N_f+i,L_pu[i][N_f+i]);
    }
    
   /* Define su*/
   d_su = (double**) malloc(2*M*sizeof(double*));
   L_su = (double***) malloc(2*M*sizeof(double**));
   for(i = 0;i<2*M;i++)
   {
        d_su[i] = (double*) malloc(2*M*sizeof(double));
        L_su[i] = (double**)malloc(2*M*sizeof(double*));
        
   }
    
    for(i=0;i<2*M;i++)
    {
        for(j=0;j<(2*M);j++)
        {
            L_su[i][j] = (double*) malloc(N_f*sizeof(double));
           
        }
    }
    
    
   /* Define L_supu*/
    d_supu = (double**) malloc(2*M*sizeof(double*));
    L_supu = (double**) malloc(2*M*sizeof(double*));
    for(i = 0;i<2*M;i++)
    {
        d_supu[i] = (double*) malloc(2*N_f*sizeof(double));
        L_supu[i] = (double*) malloc(2*N_f*sizeof(double));
    }
    
    
    for(i =0;i<M;i++)
    {
        /* Distance between transmitter reciever*/
        d_su[i][M+i]= sqrt((su[i][0]-su[M+i][0])*(su[i][0]-su[M+i][0])+(su[i][1]-su[M+i][1])*(su[i][1]-su[M+i][1]));
        L_su[i][M+i][N_f-1] = (C*C)/(pow(4*pi*(f+B*(N_f-1)),2) *pow(d_su[i][M+i],4));
        for(j=0;j<N_f;j++)
        {
            d_supu[M+i][j] = sqrt((su[M+i][0]-pu[j][0])*(su[M+i][0]-pu[j][0])+(su[M+i][1]-pu[j][1])*(su[M+i][1]-pu[j][1]));
        }
    }

    /* L_su calculations*/
    for(i=0;i<M;i++)
    {
        for(j=0;j<M;j++)
        {
            /* Distance from su i to su M+j (reciever)*/
            d_su[i][M+j] = sqrt((su[i][0]-su[M+j][0])*(su[i][0]-su[M+j][0])+(su[i][1]-su[M+j][1])*(su[i][1]-su[M+j][1]));
            for(k = 0;k<N_f;k++) 
            {
                /* Loss due to secondary user i at reciever M+j at frequency k*/
                L_su[i][M+j][k] = (C*C)/(pow((4*pi*(f+B*(k))),2)*pow(d_su[i][M+j],4));
                 mexPrintf("L_su[%d][%d][%d] = %e\n",i,M+j,k,L_su[i][M+j][k]);
            }
        }
    }
    for(i=0;i<M;i++)
    {
        for(j=0;j<N_f;j++)
        {
            d_supu[i][N_f+j] = sqrt((su[i][0]-pu[N_f+j][0])*(su[i][0]-pu[N_f+j][0])+(su[i][1]-pu[N_f+j][1])*(su[i][1]-pu[N_f+j][1]));
            /*loss due to secondary user i at primary user N_f+j(only freq.
            corresponding to primary user N_f+j)*/
            L_supu[i][N_f+j] = (C*C)/(pow((4*pi*(f+B*(j))),2)*pow(d_supu[i][N_f+j],4));
        }
    }
    
    for(i=0;i<M;i++)
    {
        for(j=0;j<N_f;j++)
        {
            freq = f+B*(j);
            L_supu[M+i][j] = (C*C)/(pow((4*pi*freq),2)*pow(d_supu[M+i][j],4));
        }
    }
    
    mexPrintf("\nM = %d, N = %lf, N_f = %d, B = %lf, f = %lf \n" ,M,N,N_f,B,f);
    
    /*display both su and pu */
   display_initialpts();
    
    FILE *fpt1;
    FILE *fpt2;
    FILE *fpt3;
    FILE *fpt4;
    FILE *fpt5;
    FILE *gp;
    FILE *fp_test;
    population *parent_pop;
    population *child_pop;
    population *mixed_pop;
    seed = .123;
    if (seed<=0.0 || seed>=1.0)
    {
        printf("\n Entered seed value is wrong, seed value must be in (0,1) \n");
        exit(1);
    }
    
    
    
    fpt1 = fopen("initial_pop.out","w");
    fpt2 = fopen("final_pop.out","w");
    fpt3 = fopen("best_pop.out","w");
    fpt4 = fopen("all_pop.out","w");
    fpt5 = fopen("params.out","w");
    fp_test = fopen("input.in","r");
    fprintf(fpt1,"# This file contains the data of initial population\n");
    fprintf(fpt2,"# This file contains the data of final population\n");
    fprintf(fpt3,"# This file contains the data of final feasible population (if found)\n");
    fprintf(fpt4,"# This file contains the data of all generations\n");
    fprintf(fpt5,"# This file contains information about inputs as read by the program\n");
    printf("\n Enter the problem relevant and algorithm relevant parameters ... ");
    printf("\n Enter the population size (a multiple of 4) : ");
    fscanf(fp_test,"%d",&popsize);
    printf("%d\n",popsize);
    if (popsize<4 || (popsize%4)!= 0)
    {
        printf("\n population size read is : %d",popsize);
        printf("\n Wrong population size entered, hence exiting \n");
        exit (1);
    }
    printf("\n Enter the number of generations : ");
    fscanf(fp_test,"%d",&ngen);
    printf("%d\n",ngen);
    if (ngen<1)
    {
        printf("\n number of generations read is : %d",ngen);
        printf("\n Wrong nuber of generations entered, hence exiting \n");
        exit (1);
    }
    printf("\n Enter the number of objectives : ");
    fscanf(fp_test,"%d",&nobj);
    printf("%d\n",nobj);
    if (nobj<1)
    {
        printf("\n number of objectives entered is : %d",nobj);
        printf("\n Wrong number of objectives entered, hence exiting \n");
        exit (1);
    }
    printf("\n Enter the number of constraints : ");
    fscanf(fp_test,"%d",&ncon);
    printf("%d\n",ncon);
    if (ncon<0)
    {
        printf("\n number of constraints entered is : %d",ncon);
        printf("\n Wrong number of constraints enetered, hence exiting \n");
        exit (1);
    }
    printf("\n Enter the number of real variables : ");
    fscanf(fp_test,"%d",&nreal);
    printf("%d\n",nreal);
    if (nreal<0)
    {
        printf("\n number of real variables entered is : %d",nreal);
        printf("\n Wrong number of variables entered, hence exiting \n");
        exit (1);
    }
    if (nreal != 0)
    {
        min_realvar = (double *)malloc(nreal*sizeof(double));
        max_realvar = (double *)malloc(nreal*sizeof(double));
        for (i=0; i<nreal; i++)
        {
            printf ("\n Enter the lower limit of real variable %d : ",i+1);
            fscanf (fp_test,"%lf",&min_realvar[i]);
            printf("%lf\n",min_realvar[i]);
            printf ("\n Enter the upper limit of real variable %d : ",i+1);
            fscanf (fp_test,"%lf",&max_realvar[i]);
            printf("%lf\n",max_realvar[i]);
            if (max_realvar[i] <= min_realvar[i])
            {
                printf("\n Wrong limits entered for the min and max bounds of real variable, hence exiting \n");
                exit(1);
            }
        }
        printf ("\n Enter the probability of crossover of real variable (0.6-1.0) : ");
        fscanf (fp_test,"%lf",&pcross_real);
        printf("%lf\n",pcross_real);
        if (pcross_real<0.0 || pcross_real>1.0)
        {
            printf("\n Probability of crossover entered is : %e",pcross_real);
            printf("\n Entered value of probability of crossover of real variables is out of bounds, hence exiting \n");
            exit (1);
        }
        printf ("\n Enter the probablity of mutation of real variables (1/nreal) : ");
        fscanf (fp_test,"%lf",&pmut_real);
        printf("%lf\n",pmut_real);
        if (pmut_real<0.0 || pmut_real>1.0)
        {
            printf("\n Probability of mutation entered is : %e",pmut_real);
            printf("\n Entered value of probability of mutation of real variables is out of bounds, hence exiting \n");
            exit (1);
        }
        printf ("\n Enter the value of distribution index for crossover (5-20): ");
        fscanf (fp_test,"%lf",&eta_c);
        printf("%lf\n",eta_c);
        if (eta_c<=0)
        {
            printf("\n The value entered is : %e",eta_c);
            printf("\n Wrong value of distribution index for crossover entered, hence exiting \n");
            exit (1);
        }
        printf ("\n Enter the value of distribution index for mutation (5-50): ");
        fscanf (fp_test,"%lf",&eta_m);
        printf("%lf\n",eta_m);
        if (eta_m<=0)
        {
            printf("\n The value entered is : %e",eta_m);
            printf("\n Wrong value of distribution index for mutation entered, hence exiting \n");
            exit (1);
        }
    }
    printf("\n Enter the number of binary variables : ");
    fscanf(fp_test,"%d",&nbin);
    printf("%d\n",nbin);
    if (nbin<0)
    {
        printf ("\n number of binary variables entered is : %d",nbin);
        printf ("\n Wrong number of binary variables entered, hence exiting \n");
        exit(1);
    }
    if (nbin != 0)
    {
        nbits = (int *)malloc(nbin*sizeof(int));
        min_binvar = (double *)malloc(nbin*sizeof(double));
        max_binvar = (double *)malloc(nbin*sizeof(double));
        for (i=0; i<nbin; i++)
        {
            printf ("\n Enter the number of bits for binary variable %d : ",i+1);
            fscanf (fp_test,"%d",&nbits[i]);
            if (nbits[i] < 1)
            {
                printf("\n Wrong number of bits for binary variable entered, hence exiting");
                exit(1);
            }
            printf ("\n Enter the lower limit of binary variable %d : ",i+1);
            fscanf (fp_test,"%lf",&min_binvar[i]);
            printf ("\n Enter the upper limit of binary variable %d : ",i+1);
            fscanf (fp_test,"%lf",&max_binvar[i]);
            if (max_binvar[i] <= min_binvar[i])
            {
                printf("\n Wrong limits entered for the min and max bounds of binary variable entered, hence exiting \n");
                exit(1);
            }
        }
        printf ("\n Enter the probability of crossover of binary variable (0.6-1.0): ");
        fscanf (fp_test,"%lf",&pcross_bin);
        if (pcross_bin<0.0 || pcross_bin>1.0)
        {
            printf("\n Probability of crossover entered is : %e",pcross_bin);
            printf("\n Entered value of probability of crossover of binary variables is out of bounds, hence exiting \n");
            exit (1);
        }
        printf ("\n Enter the probability of mutation of binary variables (1/nbits): ");
        fscanf (fp_test,"%lf",&pmut_bin);
        if (pmut_bin<0.0 || pmut_bin>1.0)
        {
            printf("\n Probability of mutation entered is : %e",pmut_bin);
            printf("\n Entered value of probability  of mutation of binary variables is out of bounds, hence exiting \n");
            exit (1);
        }
    }
    if (nreal==0 && nbin==0)
    {
        printf("\n Number of real as well as binary variables, both are zero, hence exiting \n");
        exit(1);
    }
    choice=0;
    printf("\n Do you want to use gnuplot to display the results realtime (0 for NO) (1 for yes) : ");
    fscanf(fp_test,"%d",&choice);
    printf("%d\n",choice);
    if (choice!=0 && choice!=1)
    {
        printf("\n Entered the wrong choice, hence exiting, choice entered was %d\n",choice);
        exit(1);
    }
    if (choice==1)
    {
        gp = popen(GNUPLOT_COMMAND,"w");
        if (gp==NULL)
        {
            printf("\n Could not open a pipe to gnuplot, check the definition of GNUPLOT_COMMAND in file global.h\n");
            printf("\n Edit the string to suit your system configuration and rerun the program\n");
            exit(1);
        }
        if (nobj==2)
        {
            printf("\n Enter the objective for X axis display : ");
            fscanf(fp_test,"%d",&obj1);
            printf("%d\n",obj1);
            if (obj1<1 || obj1>nobj)
            {
                printf("\n Wrong value of X objective entered, value entered was %d\n",obj1);
                exit(1);
            }
            printf("\n Enter the objective for Y axis display : ");
            fscanf(fp_test,"%d",&obj2);
            printf("%d\n",obj2);
            if (obj2<1 || obj2>nobj)
            {
                printf("\n Wrong value of Y objective entered, value entered was %d\n",obj2);
                exit(1);
            }
            obj3 = -1;
        }
        else 
        {
            printf("\n #obj > 2, 2D display or a 3D display ?, enter 2 for 2D and 3 for 3D :");
            fscanf(fp_test,"%d",&choice);
            printf("%d\n",choice);
            if (choice!=2 && choice!=3)
            {
                printf("\n Entered the wrong choice, hence exiting, choice entered was %d\n",choice);
                exit(1);
            }
            if (choice==2)
            {
                printf("\n Enter the objective for X axis display : ");
                fscanf(fp_test,"%d",&obj1);
                printf("%d\n",obj1);
                if (obj1<1 || obj1>nobj)
                {
                    printf("\n Wrong value of X objective entered, value entered was %d\n",obj1);
                    exit(1);
                }
                printf("\n Enter the objective for Y axis display : ");
                fscanf(fp_test,"%d",&obj2);
                printf("%d\n",obj2);
                if (obj2<1 || obj2>nobj)
                {
                    printf("\n Wrong value of Y objective entered, value entered was %d\n",obj2);
                    exit(1);
                }
                obj3 = -1;
            }
            else
            {
                printf("\n Enter the objective for X axis display : ");
                fscanf(fp_test,"%d",&obj1);
                printf("%d\n",obj1);
                if (obj1<1 || obj1>nobj)
                {
                    printf("\n Wrong value of X objective entered, value entered was %d\n",obj1);
                    exit(1);
                }
                printf("\n Enter the objective for Y axis display : ");
                fscanf(fp_test,"%d",&obj2);
                printf("%d\n",obj2);
                if (obj2<1 || obj2>nobj)
                {
                    printf("\n Wrong value of Y objective entered, value entered was %d\n",obj2);
                    exit(1);
                }
                printf("\n Enter the objective for Z axis display : ");
                fscanf(fp_test,"%d",&obj3);
                printf("%d\n",obj3);
                if (obj3<1 || obj3>nobj)
                {
                    printf("\n Wrong value of Z objective entered, value entered was %d\n",obj3);
                    exit(1);
                }
                printf("\n You have chosen 3D display, hence location of eye required \n");
                printf("\n Enter the first angle (an integer in the range 0-180) (if not known, enter 60) :");
                fscanf(fp_test,"%d",&angle1);
                printf("%d\n",angle1);
                if (angle1<0 || angle1>180)
                {
                    printf("\n Wrong value for first angle entered, hence exiting \n");
                    exit(1);
                }
                printf("\n Enter the second angle (an integer in the range 0-360) (if not known, enter 30) :");
                fscanf(fp_test,"%d",&angle2);
                printf("%d\n",angle2);
                if (angle2<0 || angle2>360)
                {
                    printf("\n Wrong value for second angle entered, hence exiting \n");
                    exit(1);
                }
            }
        }
    }
    printf("\n Input data successfully entered, now performing initialization \n");
    fprintf(fpt5,"\n Population size = %d",popsize);
    fprintf(fpt5,"\n Number of generations = %d",ngen);
    fprintf(fpt5,"\n Number of objective functions = %d",nobj);
    fprintf(fpt5,"\n Number of constraints = %d",ncon);
    fprintf(fpt5,"\n Number of real variables = %d",nreal);
    if (nreal!=0)
    {
        for (i=0; i<nreal; i++)
        {
            fprintf(fpt5,"\n Lower limit of real variable %d = %e",i+1,min_realvar[i]);
            fprintf(fpt5,"\n Upper limit of real variable %d = %e",i+1,max_realvar[i]);
        }
        fprintf(fpt5,"\n Probability of crossover of real variable = %e",pcross_real);
        fprintf(fpt5,"\n Probability of mutation of real variable = %e",pmut_real);
        fprintf(fpt5,"\n Distribution index for crossover = %e",eta_c);
        fprintf(fpt5,"\n Distribution index for mutation = %e",eta_m);
    }
    fprintf(fpt5,"\n Number of binary variables = %d",nbin);
    if (nbin!=0)
    {
        for (i=0; i<nbin; i++)
        {
            fprintf(fpt5,"\n Number of bits for binary variable %d = %d",i+1,nbits[i]);
            fprintf(fpt5,"\n Lower limit of binary variable %d = %e",i+1,min_binvar[i]);
            fprintf(fpt5,"\n Upper limit of binary variable %d = %e",i+1,max_binvar[i]);
        }
        fprintf(fpt5,"\n Probability of crossover of binary variable = %e",pcross_bin);
        fprintf(fpt5,"\n Probability of mutation of binary variable = %e",pmut_bin);
    }
    fprintf(fpt5,"\n Seed for random number generator = %e",seed);
    bitlength = 0;
    if (nbin!=0)
    {
        for (i=0; i<nbin; i++)
        {
            bitlength += nbits[i];
        }
    }
    fprintf(fpt1,"# of objectives = %d, # of constraints = %d, # of real_var = %d, # of bits of bin_var = %d, constr_violation, rank, crowding_distance\n",nobj,ncon,nreal,bitlength);
    fprintf(fpt2,"# of objectives = %d, # of constraints = %d, # of real_var = %d, # of bits of bin_var = %d, constr_violation, rank, crowding_distance\n",nobj,ncon,nreal,bitlength);
    fprintf(fpt3,"# of objectives = %d, # of constraints = %d, # of real_var = %d, # of bits of bin_var = %d, constr_violation, rank, crowding_distance\n",nobj,ncon,nreal,bitlength);
    fprintf(fpt4,"# of objectives = %d, # of constraints = %d, # of real_var = %d, # of bits of bin_var = %d, constr_violation, rank, crowding_distance\n",nobj,ncon,nreal,bitlength);
    nbinmut = 0;
    nrealmut = 0;
    nbincross = 0;
    nrealcross = 0;
    parent_pop = (population *)malloc(sizeof(population));
    child_pop = (population *)malloc(sizeof(population));
    mixed_pop = (population *)malloc(sizeof(population));
    
    allocate_memory_pop (parent_pop, popsize);
    allocate_memory_pop (child_pop, popsize);
    allocate_memory_pop (mixed_pop, 2*popsize);
    
    randomize();
    initialize_pop (parent_pop);
    printf("\n Initialization done, now performing first generation");
    decode_pop(parent_pop);
    evaluate_pop (parent_pop);
    assign_rank_and_crowding_distance (parent_pop);
  
    report_pop (parent_pop, fpt1);
    fprintf(fpt4,"# gen = 1\n");
    report_pop(parent_pop,fpt4);
    printf("\n gen = 1");
    fflush(stdout);
    if (choice!=0)    onthefly_display (parent_pop,gp,1);
     fflush(fpt1);
    fflush(fpt2);
    fflush(fpt3);
    fflush(fpt4);
    fflush(fpt5);
    sleep(1);
    
    for (i=2; i<=ngen; i++)
    {
        selection (parent_pop, child_pop);
        mutation_pop (child_pop);
        decode_pop(child_pop);
        evaluate_pop(child_pop);
        merge (parent_pop, child_pop, mixed_pop);
        fill_nondominated_sort (mixed_pop, parent_pop);
        /* Comment following four lines if information for all
        generations is not desired, it will speed up the execution */
        fprintf(fpt4,"# gen = %d\n",i);
        report_pop(parent_pop,fpt4);
        fflush(fpt4);
        if (choice!=0) onthefly_display (parent_pop,gp,i);
        printf("\n gen = %d",i);
    }
    
     printf("\n Generations finished, now reporting solutions");
    report_pop(parent_pop,fpt2);
    report_feasible(parent_pop,fpt3);
    /* returning solution */
    plhs[0] = mxCreateDoubleMatrix(nobj,popsize,mxREAL);
    plhs[1] = mxCreateDoubleMatrix(ncon,1,mxREAL);
    plhs[2] = mxCreateDoubleMatrix(nreal,1,mxREAL);
    plhs[3] = mxCreateDoubleMatrix(1,1,mxREAL);
    plhs[4] = mxCreateDoubleMatrix(1,1,mxREAL);
    double *p1,*p2,*p3,*p4,*p5;
    p1 = mxGetPr(plhs[0]);
    p2 = mxGetPr(plhs[1]);
    p3 = mxGetPr(plhs[2]);
    p4 = mxGetPr(plhs[3]);
    p5 = mxGetPr(plhs[4]);
    int popnum = 0; /*Rank 1 solutions*/
    for (i=0; i<popsize; i++)
    {
        if (parent_pop->ind[i].rank==1)
        {
            p1[popnum] = parent_pop->ind[i].obj[0];
            p1[popnum+1] = parent_pop->ind[i].obj[1];
            popnum+=2;
            if (ncon!=0)
            {
                for (j=0; j<ncon; j++)
                {
                    p2[j] = parent_pop->ind[i].constr[j];
                }
            }
            if (nreal!=0)
            {
                for (j=0; j<nreal; j++)
                {
                    p3[j] = parent_pop->ind[i].xreal[j];
                }
            }
            p4[0]  = parent_pop->ind[i].constr_violation;
       }
    }
    p5[0] = popnum/2;
    for(i=0;i<2*M;i++)
    {
        for(j=0;j<2*N_f;j++)
        {
            printf("L_supu[%d][%d] = %e \n",i,j,L_supu[i][j]);
        }
    }
    
    if (nreal!=0)
    {
        fprintf(fpt5,"\n Number of crossover of real variable = %d",nrealcross);
        fprintf(fpt5,"\n Number of mutation of real variable = %d",nrealmut);
    }
    if (nbin!=0)
    {
        fprintf(fpt5,"\n Number of crossover of binary variable = %d",nbincross);
        fprintf(fpt5,"\n Number of mutation of binary variable = %d",nbinmut);
    }
    fflush(stdout);
    fflush(fpt1);
    fflush(fpt2);
    fflush(fpt3);
    fflush(fpt4);
    fflush(fpt5);
  
    if (choice!=0)
    {
        pclose(gp);
    }
    if (nreal!=0)
    {
        free (min_realvar);
        free (max_realvar);
    }
    if (nbin!=0)
    {
        free (min_binvar);
        free (max_binvar);
        free (nbits);
    }
    
    deallocate_memory_pop (parent_pop, popsize);
    deallocate_memory_pop (child_pop, popsize);
    deallocate_memory_pop (mixed_pop, 2*popsize);
    free (parent_pop);
    free (child_pop);
    free (mixed_pop);
   
    
    
    /* Free d_pu and L_pu*/
    for(i=0;i<2*N_f;i++)
    {
        free(d_pu[i]);
        free(L_pu[i]);
    }
    free(d_pu);
    free(L_pu);
    
    
    for(i=0;i<2*M;i++)
    {
        free(su[i]);
    }
    free(su);
    
    for(i=0;i<2*N_f;i++)
    {
        free(pu[i]);
    }

    free(pu);
    for(i=0;i<2*M;i++)
        free(d_su[i]);
    free(d_su);
     for(i=0;i<2*M;i++)
    {
        for(j=0;j<(2*M);j++)
        {
            free(L_su[i][j]);
           
        }
    }
    for(i=0;i<2*M;i++)
        free(L_su[i]);
    free(L_su);
     for(i = 0;i<2*M;i++)
    {
        free(d_supu[i]);
        free(L_supu[i]);
    }
    free(d_supu);
    free(L_supu);
    
    
    printf("\n Routine successfully exited \n");
   
        fclose(fp_test);
        fclose(fpt1);
        fclose(fpt2);
        fclose(fpt3);
        fclose(fpt4);
        fclose(fpt5);

	/*mxArrayDestroy(mx);
	mxArrayDestroy(mf);*/
}  
	 
