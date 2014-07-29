/* Routines to display the population information using gnuplot */

/* Function to display the current population for the subsequent generation */
void onthefly_display (population *pop, FILE *gp, int ii)
{
    int i,j;
    int flag;
    FILE *fpt;
    fpt = fopen("plot.out","w");
    flag = 0;
    double *p;
    mxArray *rhs[1];
    rhs[0] = mxCreateDoubleMatrix(popsize,2,mxREAL);
    p = mxGetPr(rhs[0]);
    j = 0;
    for (i=0; i<popsize; i++)
    {
        if (pop->ind[i].constr_violation==0)
        {
            if (choice!=3)
            {
                /*fprintf(fpt,"%e\t%e\n",pop->ind[i].obj[obj1-1],pop->ind[i].obj[obj2-1]);*/
                p[j] = pop->ind[i].obj[obj1-1];
                p[j+popsize] = pop->ind[i].obj[obj2-1];
                j = j+1;
            }
            else
                fprintf(fpt,"%e\t%e\t%e\n",pop->ind[i].obj[obj1-1],pop->ind[i].obj[obj2-1],pop->ind[i].obj[obj3-1]);
            fflush(fpt);
            flag = 1;
        }
    }
    mexCallMATLAB(0,NULL,1,rhs,"plotter");
    if (flag==0)
    {
        printf("\n No feasible soln in this pop, hence no display");
    }
    /*else
    {
        if (choice!=3){
            fprintf(gp,"set title 'Generation #%d'\n unset key\n plot 'plot.out' w points pointtype 6 pointsize 1\n",ii);
         }
        else
            fprintf(gp,"set title 'Generation #%d'\n set view %d,%d\n unset key\n splot 'plot.out' w points pointtype 6 pointsize 1\n",ii,angle1,angle2);
        fflush(gp);
    }*/
    fclose(fpt);
    mxDestroyArray(rhs[0]);
    return;
}



/* function to display the initial population*/
void display_initialpts()
{
    int i,puno=0;
    double *p;
    mxArray *rhs[1];
    rhs[0] = mxCreateDoubleMatrix(2*(M+N_f),2,mxREAL);
    p = mxGetPr(rhs[0]);
    for(i = 0;i<2*M;i++)
    {
        p[i] = su[i][0];
        p[2*(M+N_f)+i] = su[i][1];
    }
    for(i =2*M;i<2*(M+N_f);i++)
    {   puno = i-2*M; 
        p[i] = pu[puno][0];
        p[2*(M+N_f)+i] = pu[puno][1];
    }
    mexCallMATLAB(0,NULL,1,rhs,"plotter");
    mxDestroyArray(rhs[0]);
    return;
    
}
    
        
    
    