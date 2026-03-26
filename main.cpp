#include <cstdio>        // printf, fprintf 등
#include <cstdlib>       // exit, atoi, atof 등
#include <mpi.h>
#include <cmath>         // 수학 함수
#include <ctime>         // time
#include <complex>       // std::complex
#include <string>
#include <iostream>
#include <hdf5.h>
#include <hdf5_hl.h>

#include "mesh.h"        // Domain, parameterSetting 등이 정의된 헤더 가정

int main(int argc, char *argv[])
{
   int iteration,s=0;
 
   MPI_Init(&argc,&argv);
    
   int myrank=0, nTasks=1;

   //MPI_Status status; 
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   Domain D{};     // zero-initialization
   LoadList LL{};

   // 입력 파일이 제대로 주어졌는지 확인
   if (argc < 2)  {
      if (myrank == 0) {
         std::fprintf(stderr, "Usage: %s <input_file>\n", argv[0]);
      }
      MPI_Finalize();
      return 1;
   }

   parameterSetting(&D,argv[1]);

   boundary(&D);

   if(argc>=3) {
      iteration=atoi(argv[2]);


   } else {

      //loading  beam
      iteration=0;
      s=0;
      for (auto& LL : D.loadList) {
         loadBeam(&D,LL,s,iteration);
         s++;
      }


   }

   for (size_t s=0; s<D.loadList.size(); ++s) {
      saveParticlesToTxt(D, s, "test");
   }


   MPI_Finalize();

   return 0;
}
