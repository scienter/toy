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

std::vector<std::vector<cplx>> complexMemoryFlat(int harmony, int size);
std::vector<std::vector<double>> doubleMemoryFlat(int size1, int size2);

void boundary(Domain *D)
{

   int myrank=0, nTasks=1;
   //MPI_Status status;

   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);     
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);     

   //Static
   D->subSliceN = D->sliceN;
   D->minI = 0;
   D->maxI = 1;

   // Field memory setting
   std::cout << "Here. nx=" << D->nx
             << ", ny=" << D->ny
             << std::endl;
   D->Ux=complexMemoryFlat(D->numHarmony,(D->subSliceN+2)*D->nx*D->ny);
   D->Uy=complexMemoryFlat(D->numHarmony,(D->subSliceN+2)*D->nx*D->ny);
   D->ScUx=complexMemoryFlat(D->numHarmony,(D->subSliceN+2)*D->nx*D->ny);
   D->ScUy=complexMemoryFlat(D->numHarmony,(D->subSliceN+2)*D->nx*D->ny);
   D->Ez=complexMemoryFlat(D->SCLmode,D->SCFmode*(D->subSliceN+2)*D->nr);
   
   D->totalEnergyX=doubleMemoryFlat(D->maxStep,D->numHarmony);
   D->totalEnergyY=doubleMemoryFlat(D->maxStep,D->numHarmony);

   D->particle.resize(D->subSliceN+2);
   for(auto& p : D->particle) {
      p.head.resize(D->nSpecies);
      for(size_t sp = 0; sp < p.head.size(); ++sp) {
        if (p.head[sp] == nullptr) {
            p.head[sp] = new ptclHead{};
            p.head[sp]->pt = nullptr;
        }
      }

   }

   // Bessel Table
   D->BesselMax = D->harmony[D->numHarmony - 1] * 2.0 * M_PI;
   D->dBessel = D->BesselMax / static_cast<double>(D->bn);
   D->BesselMaxOrder = (D->harmony[D->numHarmony - 1] + 1) / 2 + 1;

   // 2차원 vector로 변경 (메모리 관리 자동)
   D->BesselJ.resize(D->bn);                    // bn 개의 행
   for (int i = 0; i < D->bn; ++i)
      D->BesselJ[i].resize(D->BesselMaxOrder); // 각 행에 BesselMaxOrder 개의 열

   for (int i = 0; i < D->bn; ++i) {
      double xi = D->dBessel * i;
      for (int j = 0; j < D->BesselMaxOrder; ++j) {
         D->BesselJ[i][j] = gsl_sf_bessel_Jn(j, xi);
      }
   }

}

std::vector<std::vector<double>> doubleMemoryFlat(int size1, int size2)
{
   std::vector<std::vector<double>> field(size1);
   for(int h=0; h<size1; ++h) 
      field[h].resize(size2);

   return field;
}


std::vector<std::vector<cplx>> complexMemoryFlat(int harmony, int size)
{
   std::vector<std::vector<cplx>> field(harmony);
   for(int h=0; h<harmony; ++h) 
      field[h].resize(size, cplx{0.0,0.0});

   return field;
}

