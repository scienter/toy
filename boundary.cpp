#include <cmath>
#include <iostream>
#include <mpi.h>
#include <complex>
#include <vector>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_matrix_complex_double.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_sf_bessel.h>

#include "mesh.h"
#include "constants.h"

std::vector<cplx> complexMemoryFlat(int harmony, int nz, int nx, int ny);

void boundary(Domain *D)
{

   int myrank=0, nTasks=1;
   //MPI_Status status;

   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);     
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);     


   D->subSliceN = D->sliceN;

   // Field memory setting
   D->Ux=complexMemoryFlat(D->numHarmony,D->subSliceN+2,D->nx,D->ny);
   D->Uy=complexMemoryFlat(D->numHarmony,D->subSliceN+2,D->nx,D->ny);
   D->ScUx=complexMemoryFlat(D->numHarmony,D->subSliceN+2,D->nx,D->ny);
   D->ScUy=complexMemoryFlat(D->numHarmony,D->subSliceN+2,D->nx,D->ny);
/*
   D->totalEnergyX=(double **)malloc(D->maxStep*sizeof(double *));
   D->totalEnergyY=(double **)malloc(D->maxStep*sizeof(double *));
   for(i=0; i<D->maxStep; i++) {
      D->totalEnergyX[i]=(double *)malloc(D->numHarmony*sizeof(double ));
      D->totalEnergyY[i]=(double *)malloc(D->numHarmony*sizeof(double ));
   }

   // space charge
   nz = D->subSliceN+2;
   D->Ez=(double complex ****)malloc(nz*sizeof(double complex ***));
   for(i=0; i<nz; i++) {
      D->Ez[i]=(double complex ***)malloc(D->nr*sizeof(double complex **));
      for(j=0; j<D->nr; j++) {
         D->Ez[i][j]=(double complex **)malloc(D->SCLmode*sizeof(double complex *));
         for(l=0; l<D->SCLmode; l++)
            D->Ez[i][j][l]=(double complex *)malloc(D->SCFmode*sizeof(double complex ));
     }
   }
   for(i=0; i<nz; i++)
      for(j=0; j<D->nr; j++)
         for(l=0; l<D->SCLmode; l++)
            for(m=0; m<D->SCFmode; m++)
               D->Ez[i][j][l][m]=0.0+I*0.0;
*/


   // setting up particle's pointer
   //LL=D->loadList;
   //s=0; 
   //while(LL->next) {
   //  totalCnt=LL->numBeamlet*LL->numInBeamlet;
   //  LL->totalCnt=totalCnt;
   //  LL=LL->next;
   //  s++;
   //}

   D->particle.resize(D->subSliceN+2);
   for(auto& p : D->particle) {
      p.head.resize(D->nSpecies);
   }
/*
   D->avePx=0.0;
   D->avePy=0.0;
   D->aveGam=D->gamma0;
   D->totalCnt=1.0;

   // initialize prev B field
   D->prevK=0.0;

   // twiss parameter
   D->twsBX = (double *)malloc((D->maxStep+1)*sizeof(double ));
   D->twsGX = (double *)malloc((D->maxStep+1)*sizeof(double ));
   D->twsAX = (double *)malloc((D->maxStep+1)*sizeof(double ));
   D->twsEmitX = (double *)malloc((D->maxStep+1)*sizeof(double ));
   D->twsBY = (double *)malloc((D->maxStep+1)*sizeof(double ));
   D->twsGY = (double *)malloc((D->maxStep+1)*sizeof(double ));
   D->twsAY = (double *)malloc((D->maxStep+1)*sizeof(double ));
   D->twsEmitY = (double *)malloc((D->maxStep+1)*sizeof(double ));
   D->twsG = (double *)malloc((D->maxStep+1)*sizeof(double ));

   // Bessel Table
   D->BesselMax = D->harmony[D->numHarmony-1]*2*M_PI;
   D->dBessel = D->BesselMax/(1.0*D->bn);
   D->BesselMaxOrder = (D->harmony[D->numHarmony-1]+1)/2+1;
   D->BesselJ = (double **)malloc(D->bn*sizeof(double *));
   for(i=0; i<D->bn; i++) 
      D->BesselJ[i] = (double *)malloc(D->BesselMaxOrder*sizeof(double ));
   for(i=0; i<D->bn; i++) {
      xi = D->dBessel*i;
      for(j=0; j<D->BesselMaxOrder; j++) 
         D->BesselJ[i][j]=gsl_sf_bessel_Jn(j,xi);
   }

   // define Matrix Ma
   gsl_complex g;
   int stat;
   nx=D->nx;   ny=D->ny;

   // wake field
   D->den=(double *)malloc(sliceN*sizeof(double ));   
   D->wakeF=(double *)malloc(sliceN*sizeof(double ));   
   D->wakeE=(double *)malloc(sliceN*sizeof(double ));   
   for(i=0; i<sliceN; i++) {
      D->den[i]=0.0;
      D->wakeF[i]=0.0;
      D->wakeE[i]=0.0;
   }

   D->currentFlag=OFF;
   D->ue=0;
*/
}

std::vector<cplx> complexMemoryFlat(int harmony, int nz, int nx, int ny)
{
   int total = harmony * nz * nx * ny;

   // 가장 추천: 값 초기화하면서 할당
   //return std::vector<cplx>(total);  // C++11 이후 default-init (complex는 0)
   // 또는 명시적으로 0 초기화
   return std::vector<cplx>(total, cplx{0.0, 0.0});
}




/*
double complex **complexMemoryAsign(int harmony,int nx,int ny)
{
   int h,i,N;
   double complex **field;

   N=nx*ny;
   field = (double complex **)malloc(harmony*sizeof(double complex *));
   for(h=0; h<harmony; h++) 
     field[h] = (double complex *)malloc(N*sizeof(double complex ));
   
   for(h=0; h<harmony; h++)
     for(i=0; i<N; i++)
       field[h][i]=0.0+I*0.0;

   return field;
}

double complex ***complexMemory3Asign(int harmony,int nz,int nx,int ny)
{
   int n,h,i,j,N;
   double complex ***field;

   N=nx*ny;
   field = (double complex ***)malloc(harmony*sizeof(double complex **));
   for(h=0; h<harmony; h++) {
      field[h] = (double complex **)malloc(nz*sizeof(double complex *));
      for(i=0; i<nz; i++) 
         field[h][i] = (double complex *)malloc(N*sizeof(double complex ));
   }
   
   for(h=0; h<harmony; h++)
     for(i=0; i<nz; i++)
       for(j=0; j<N; j++)
         field[h][i][j]=0.0+I*0.0;

   return field;
}

double ***doubleMemory3Asign(int harmony,int nz,int nx,int ny)
{
   int n,h,i,j,N;
   double ***field;

   N=nx*ny;
   field = (double ***)malloc(harmony*sizeof(double **));
   for(h=0; h<harmony; h++) {
      field[h] = (double **)malloc(nz*sizeof(double *));
      for(i=0; i<nz; i++) 
         field[h][i] = (double *)malloc(N*sizeof(double ));
   }
   
   for(h=0; h<harmony; h++)
     for(i=0; i<nz; i++)
       for(j=0; j<N; j++)
         field[h][i][j]=0.0;

   return field;
}
*/
