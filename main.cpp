#include <cstdio>        // printf, fprintf 등
#include <cstdlib>       // exit, atoi, atof 등
#include <mpi.h>
#include <cmath>         // 수학 함수
#include <ctime>         // time
#include <complex>       // std::complex
#include <string>
#include <hdf5.h>
#include <hdf5_hl.h>

#include "mesh.h"        // Domain, parameterSetting 등이 정의된 헤더 가정

int main(int argc, char *argv[])
{
   MPI_Init(&argc,&argv);
    
   int myrank=0, nTasks=1;

   //MPI_Status status; 
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   Domain D{};     // zero-initialization

   // 입력 파일이 제대로 주어졌는지 확인
   if (argc < 2)  {
      if (myrank == 0) {
         std::fprintf(stderr, "Usage: %s <input_file>\n", argv[0]);
      }
      MPI_Finalize();
      return 1;
   }

   parameterSetting(&D,argv[1]);



   MPI_Finalize();

   return 0;
}
