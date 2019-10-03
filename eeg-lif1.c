#include <stdio.h>
#include <math.h>
#include <stdlib.h>
/* #include <mpi.h>  for parallel */

/* to run the program use */
/* gcc -std=c99 -o a.out eeg-lif1.c -lm */
/* ulimit -s unlimited (Type this at the Linux command line for large memory NEEDS) */
/* to run in mpi-parallel use */
/* mpicc p-lif1d.c -o a.out */
/* mpirun -np 4 a.out > out.out */

/* Modified LIF: du/dt = mu -alpha *u */

double myrand()
{
	return ((double)rand()/(double)RAND_MAX);
}

int main(int argc, char **argv)
{                
//Main starts

/* FILE *ifp,*/
FILE *opf1, *opf2, *opf3, *opf4, *ipf;
ipf = fopen("00166_40Hz_matrix.txt", "r");  //for reading
opf1 = fopen("profile-00166_40Hz_matrix.trial2", "w");
opf2 = fopen("spacetime-00166_40Hz_matrix.trial2", "w");
opf3 = fopen("wmega-00166_40Hz_matrix.trial2","w");
opf4 = fopen("time-00166_40Hz_matrix.trial2","w"); 

/* **********************    WORKING PARAMETERS  *************************** */ 
     
/* mu=1.0, alpha=1.0, uth=0.0248, sigm=100.0  (40Hz)  files ending with  "....40Hz_matrix.trial2" */

/* ************************************************************************* */

  /* variable definition: */
int  n=80;
int i=1, j, k, kk, jj;
int it;
double pi=3.14149;
char octo = 'z'; /* PROBLEM WITH CHARACTERS */ 
double dt=0.0001 ;
double tstep=1.0/dt;
int ntstep=(int) tstep;
int itime = 2001*ntstep;
int ttransient=500*ntstep ;
double  alpha=1.0;
double mu=1.0;
double uth=0.0248*mu;
double sigm=30.00;      /* ???????? !!!!!!!!!!!!SOS for positive/excitatory coupling */
/*  double sigm=-1.0;      !!!!!!!!!!!!SOS for negative/inhibitory coupling */
double u[n],uplus[n],divide[n];
double sigma[n][n];
double wmega1[n],wmega[n];
double sum,fj;
double time;
double ss;
double pr=0.0;
double tot=0.0; 

int seed= 395849566;


/* read connectivity matrix sigma[n][n] from file ifp */

    

for (i = 0; i < n; i++) {
     wmega[i]=0.0;
     wmega1[i]=0.0;
     divide[i]=0.0;
   for (j=0; j < n; j++){
     sigma[i][j]=0.0;
     }
     }
     
for (i = 0; i < n; i++){
       for (j=0; j < n; j++)
    {
     fscanf(ipf,"%d  %d %lf \n", &k, &kk, &ss);
     sigma[k-1][kk-1]=ss;
     }
     }

    for(i=0; i<6; i++)
     {
     kk=1;
     printf("%d \t %d \t %.18f \n",kk,i,sigma[kk][i]);
     } 
for (i = 0; i < n; i++) {
   for (j=0; j < n; j++){
     divide[i]=divide[i]+sigma[j][i];
     }
/*     printf("%d divide  %.18f \n",i,divide[i]); */
    divide[i]=1.0; 
     }

/*     for(i=0; i<n; i++)
     {
     tot = tot + sigma[i,j];
     }
     printf("%f",tot); */

  /* initialize elements of arrays to 0 */         
   for ( i = 0; i < n; i++ ) {
      u[i] = 0.0; /* set all elements equal to 0 */
      uplus[i] =0.0;
   }

/* initialize elements of array u[i] to random */         
   for ( i = 0; i < n; i++ ) {
      u[i] =uth*myrand() ; /* set all elements equal to random */
/*     printf(" %d \t %f \n",i,u[i]); */
    }
printf("Out1 \n");

/* TEMPORAL ITERATION--------------------------------*/
   for ( it = 0; it < itime; it++ ) {
          time=(double)it;
          time=time * dt;
 
      for (i=0; i<n; i++) {               /* iteration over elements */
           uplus[i] = u[i]+dt*(mu - alpha* u[i]) ;
           sum=0.0;
           for(j=0; j<n; j++)
              {
               kk=i+j;
               kk=kk % n;
               sum=sum+sigma[j][i]*(u[j]-u[i]);
              }
           uplus[i]=uplus[i]+dt*sigm*sum/divide[i];
           u[i]=uplus[i];
           if(u[i] > uth)
           {
           u[i]=0.0;
           } 
           if(it > ttransient && u[i] == 0) // calculate wmega's
           {
             wmega1[i]=wmega1[i]+1.0;
            }      
/*                  print-out                   */    
         if (it%(1000*ntstep)==0 )                  
           {
             fprintf(opf2,"%d \t %d \t %.15f \n",i,it,u[i]);
/*             printf("%d \t %d \t %f \n",it,i,u[i]); */

            }
          if (it%(1001)==0)
            {
            fprintf(opf1," %d \t %.15f \t %d \n",i,u[i],it); 
            }
         }                       /* end of iteration over elements */
         if (it%(1000)==0 )                  
           {
/*             printf(" time is %d \n",it); */
           fprintf(opf4,"%d %.15f %.15f %.15f \n",it,u[27],u[28],u[29]); 
           }
          if (it > ttransient && (it%(200*ntstep))==0  )              
            {
             printf(" time is %d sum is %.15f \n",it,sum); 
             fprintf(opf1," \n");
/*             fprintf(opf4,"%d %.15f %.15f %.15f \n",it,u[10],u[300],u[1111]); */
                 for (j=0; j<n; j++)
                     {
                      jj=j+1;
                     wmega[j]=wmega1[j]/(time-ttransient*dt);
                      fj=wmega[j];
                     wmega[j]=2.0*pi*wmega[j];
                     fprintf(opf3,"%d %.15f %.15f %d \n",jj,fj,wmega[j],it);
                     }
                               fprintf(opf3," \n");
              }
             if (it%(1001)==0)
             {  
                 fprintf(opf1," \n");
             }      
 /*                 end of print-out                   */
      
                                    

     }             
/* END OF TEMPORAL ITERATION-----------------------------*/


printf("Fin");

}           //End-of-main
