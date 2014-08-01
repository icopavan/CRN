/* Test problem definitions */

# define thru 
#ifdef thru
void test_problem(double *xreal,double *xbin, int **gene, double *obj,double *constr)
{
  double *powerlevels,*I,*SINR,*P;
  int *power_level_no;
  int *trans,*freq_trans;
  int i,puno,suno,j = 0;
  double pwr_sum=0;
  /*double SINRSULT = 10;*/
  /*double SINRPULT = 100;*/
  double N_fm = (double)N_f;
  double objective = 0.0;
  double total_through = 0.0, variance = 0.0, temp_through = 0.0;
  double mean_through;
  int no_active;
  int no_active_users = 0;
  for(i=0;i<M*N_f;i++)
  {   
      /*printf("\nxreal[%d] = %lf",i,xreal[i]);*/
      if(xreal[i] > EPS)
      {
          j++;
	  }
   }
  no_active = j;
  int bands;
  for(i=0;i<M;++i)
  {
      for(bands = 0;bands<N_f;++bands)
          if(xreal[i*N_f+bands]>EPS)
          {
            no_active_users++;
            break;
          }
  }
  /*printf("\nno_active = %d",no_active);*/
  j = 0;
  powerlevels = (double*)malloc(no_active*sizeof(double));
  P = (double*)malloc(no_active*sizeof(double));
  
  for(i=0;i<M*N_f;i++)
  {
      if(xreal[i] >EPS)
      {
         powerlevels[j] = i;
         P[j] =xreal[i];
         j++;
      }
  }
  double temp1,temp2;    
  trans = (int*) malloc(no_active*sizeof(int));
  freq_trans = (int*) malloc(no_active*sizeof(int));
  /*printf("N_fm = %lf",N_fm);*/
  for(i=0;i<no_active;i++)
  {
      /*printf("\npowerlevels[%d] = %lf",i,powerlevels[i]); */
      /*printf("%lf=\n",floor((powerlevels[i])/(N_fm)));*/
      temp1 =floor(powerlevels[i]/N_fm);
      temp2 =fmod(powerlevels[i],N_fm);
      trans[i] = temp1;
      freq_trans[i] = temp2;
      /*printf("\ntrans[%d] = %d, freq_trans[%d] = %d",i,trans[i],i,freq_trans[i]);*/
   }
  
  power_level_no = (int*) malloc(no_active*sizeof(int));
  
  for(i = 0;i<no_active;i++)
  {
     power_level_no[i] = i;
  }
  
  I = (double *)malloc(no_active*sizeof(double));
  SINR = (double*) malloc(no_active*sizeof(double));
  
  /* Comment this*/
 /* for(i=0;i<no_active;i++)   
        printf("\nP[%d] = %e",i, P[power_level_no[i]]);*/
  
  for(i = 0;i<no_active;i++)
  {
     I[i]=0;
     for(j =0;j<no_active;j++)
     {
         if(j==i)
         {
             continue;
         }
         if(freq_trans[j] == freq_trans[i])
         {   
             if(xreal[M*N_f+trans[j]]>=alpha_su[trans[j]][M+trans[i]]-(theta/2) &&  xreal[M*N_f+trans[j]]<=alpha_su[trans[j]][M+trans[i]]+(theta/2))
             {
                 G = G_M;
             }
             else
                 G = G_S;
             I[i] = I[i]+ P[power_level_no[j]]*G*L_su[trans[j]][M+trans[i]][freq_trans[i]];
         }
     }
     /*Considering PU are omindirectional*/
     I[i] = I[i]+N+pwr_pu*G_M*L_supu[M+trans[i]][freq_trans[i]];
     if(xreal[M*N_f+trans[i]]>=alpha_su[trans[i]][M+trans[i]]-(theta/2) &&  xreal[M*N_f+trans[i]]<=alpha_su[trans[i]][M+trans[i]]+(theta/2))
     {
         G = G_M;
     }
     else
     {
         G = G_S;
     }
     SINR[i] = P[power_level_no[i]]*G*(L_su[trans[i]][M+trans[i]][freq_trans[i]])/I[i];
     /*printf("\nI[%d] = %e",i,I[i]);
     printf("\nSINR[%d] = %lf",i,SINR[i]);*/
  }
 
 objective = 0; 
 double *through_per_su;
 through_per_su = (double*)malloc(M*sizeof(double));
 for(i = 0;i<M;++i)
     through_per_su[i] = 0.0;
 for(i = 0;i<no_active;i++)
 {
     objective = objective+ B*log2(1+SINR[i]);
     through_per_su[trans[i]]+=B*log2(1+SINR[i]);
 }
 mean_through = objective/(double)M;
 variance = 0.0;
 /*Variance calculation*/
 for(i=0;i<M;++i)
     variance+=(through_per_su[i]-mean_through)*(through_per_su[i]-mean_through);
 
 obj[0] = -objective;
 obj[1] = variance;
 
/* printf("\nobjective = %e",obj[0]);*/

 /* Constraints */
 /*su constraints*/
 j=0;
 for(i=0;i<M*N_f;i++)
 {
     if(powerlevels[j] ==i)
     {
         constr[i] = (SINR[j]-SINRSULT)/SINRSULT;
         j++;
     }
     else
         constr[i] = 0;
 }
 
 double *SINRPU,*I_pu;
 I_pu = (double*) malloc(N_f*sizeof(double));
 SINRPU = (double*) malloc(N_f*sizeof(double));
/*pu constraints*/
 for(i = M*N_f;i<M*N_f+N_f;i++)
 {
     puno = i-M*N_f;
     I_pu[puno] = 0;
     for(j=0;j<no_active;j++)
     {
         if(freq_trans[j] == puno)
         {
             if(xreal[M*N_f+trans[j]]>=alpha_pu[trans[j]][N_f+puno]-(theta/2) &&  xreal[M*N_f+trans[j]]<=alpha_pu[trans[j]][N_f+puno]+(theta/2))
                 G = G_M;
             else
                 G = G_S;
             I_pu[puno] = I_pu[puno]+ P[j]*G*L_supu[trans[j]][N_f+puno];
         }
     }
     SINRPU[puno] = (pwr_pu*G_M*L_pu[puno][N_f+puno])/(N+I_pu[puno]);
     /*printf("\nSINRPU[%d] = %lf",puno,SINRPU[puno]);*/
     constr[i] = (SINRPU[puno]/SINRPULT)-1;
 }
 
 /*powerbudget constraints*/
 for(i = M*N_f+N_f;i<(M*N_f+N_f+M);i++)
 {
     suno = i-(M*N_f+N_f);
     pwr_sum = 0;
     for(j=0;j<no_active;j++)
     {
         if(trans[j] == suno)
             pwr_sum = pwr_sum+P[j];
     }
     constr[i] = (pwr_max-pwr_sum)/pwr_max;
 }
/*print constraints*/
/* 
for(i = 0;i<M*N_f+M+N_f;i++)
{
     printf("constr[%d] = %lf \n", i, constr[i]);
}*/
free(powerlevels);
free(P);
free(trans);
free(freq_trans);
free(I);
free(SINR);
free(I_pu);
free(SINRPU);
free(through_per_su);
}
#endif
