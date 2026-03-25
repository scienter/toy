#ifndef MESH_H
#define MESH_H

#include <complex>
#include <vector>
#include "particle.h"

using cplx = std::complex<double>;

enum class OperationMode {
    Unknown    = 0,
    Static     = 1,
    Time_Dependent = 2,
    Twiss      = 3
};

typedef struct _Domain
{
   int dimension;
   double ks,lambda0;
   OperationMode mode;
   int sliceN,subSliceN,nx,ny;   

   int numHarmony, *harmony;
   std::vector<cplx> Ux,Uy,ScUx,ScUy;


   //Electron beam
   int nSpecies=0;
   std::vector<LoadList> loadList;
   std::vector<Particle> particle;


}  Domain; 

void parameterSetting(Domain *D,const char *input);
int FindParameters (const char *block, int rank, const char *options, const char *input, char *ret);
void boundary(Domain *D);


#endif   //MESH_H
