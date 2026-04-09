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
      //Create "totalEnergy", "twissFile" file
      FILE *out1 = fopen("totalEnergy", "w");      
      fclose(out1);
      FILE *out2 = fopen("twissFile", "w");      
      fprintf(out2,"#%12s %12s %12s %12s %12s %12s %12s\n",
                  "z","emitX","betaX","alphaX","emitY","betaY","alphaY");
      fclose(out2);

      //loading  beam
      iteration=0;
      s=0;
      for (auto& LL : D.loadList) {
         loadBeam(&D,LL,s,iteration);
         s++;
      }

   }



   while(iteration<D.maxStep) 
   {
      if(iteration%20==0) {
         for (size_t s=0; s<D.loadList.size(); ++s) {
            std::string fileName = "Field" + std::to_string(iteration);
            saveFieldsToTxt(D, fileName);
            fileName = std::to_string(s) + "Particle" + std::to_string(iteration);
            saveParticlesToTxt(D, s, fileName);
         }
      }

      // Update Files
      updateTotalEnergy(&D,iteration);
      calculate_twiss(D,iteration);


      solveField(&D,iteration);
      //std::cout << "iteration=" << iteration << "after solveField" << std::endl;

      updateK_quadG(&D,iteration,0);
      //std::cout << "iteration=" << iteration << "after 1st updateK" << std::endl;

      if(D.dimension==3) 
         transversePush(&D,iteration);

      updateK_quadG(&D,iteration,0.5);
      //std::cout << "iteration=" << iteration << "after 2nd updateK" << std::endl;

      if(D.dimension==3) 
         transversePush(&D,iteration);

      if(D.driftFlag==false) push_theta_gamma(&D,iteration);
      else {
         std::cout << "iteration=" << iteration
                   << ", driftON"
                   << std::endl;
      }
      //std::cout << "iteration=" << iteration << "push_theta_gamma" << std::endl;

      if(iteration%10==0) {
         if(myrank==0) 
            printf("iteration=%d, z=%g\n",iteration,iteration*D.dz);
      }
      iteration++;
   }



   MPI_Finalize();

   return 0;
}
