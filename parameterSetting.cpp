#include <iostream>
#include <string>
#include <cstring>
#include <vector>
//#include <cctype>
#include <cstdlib>      // std::stoi, std::stod
#include <cmath>        // M_PI 등
#include <random>       // std::random_device, std::mt19937 등
#include <mpi.h>

#include "mesh.h"
#include "constants.h"
#include <gsl/gsl_sf_bessel.h>

using std::cout;


// 난수 생성기 (전역 또는 static으로 관리하는 것이 좋음)
static std::mt19937& get_random_engine()
{
    static std::random_device rd;
    static std::mt19937 gen(rd());
    return gen;
}

double randomV()
{
    static std::uniform_real_distribution<double> dist(0.0, 1.0);
    return dist(get_random_engine());
}

OperationMode whatOperMode(char *str);

void parameterSetting(Domain *D,const char *input)
{
   int myrank=0, nTasks=1;
   //MPI_Status status;

   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);     
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);     

   int i=0;
   double energy=0.0;
   bool fail = false;
   char str[100];

   //FILE *in=NULL;

   if(FindParameters("Domain",1,"dimension",input,str)) D->dimension=atoi(str);
   else  { 
      cout << "in [Domain], dimension=?  (1:1D, 2:2D, 3:3D) \n";
      fail=true;  
   }
   if(FindParameters("Domain",1,"mode",input,str)) D->mode=whatOperMode(str);
   else  { 
      cout << "in [Domain], mode=?  (static, time) \n";
      fail=true;  
   }
   if(FindParameters("Domain",1,"photon_energy",input,str)) energy=atof(str);
   else  { 
      cout << "in [Domain], photon_energy=?  [eV] \n";
      fail=true;  
   }

   if(FindParameters("Domain",1,"num_harmony",input,str)) D->numHarmony=atoi(str);
   else { D->numHarmony=1; }

	D->harmony = new int[D->numHarmony];
   for(i=0; i<D->numHarmony; i++) {
      std::string param_name = "harmony" + std::to_string(i);
      char buf[100] = {};
      if(FindParameters("Domain",1,param_name.c_str(),input,buf)) D->harmony[i] = std::atoi(buf);
      else  { 
         std::cerr << "In [Domain], {} should be defined.\n";
         fail=true; 
      }
   }

   D->ks=energy*2*M_PI/(plankH*velocityC);
   D->lambda0=2*M_PI/D->ks;

   //Setting
   if (D->mode == OperationMode::Static) D->sliceN=1;
   if (myrank==0) { 
      cout << "lambda0=" << D->lambda0 
           << ", numHarmony=" << D->numHarmony 
           << ", opermode=" << static_cast<int>(D->mode)
           << ", sliceN=" << D->sliceN
           << "\n"; 
   }


   if(fail) { 
      MPI_Finalize();
      std::exit(1);
   }
}

OperationMode whatOperMode(char *str)
{
   if (!str || !*str) {
      return OperationMode::Unknown;
   }

   if (strstr(str,"Static") != nullptr)
      return OperationMode::Static;
   else if (strstr(str,"Time_Dependent") != nullptr)
      return OperationMode::Time_Dependent;
   else if (strstr(str,"Twiss") != nullptr)
      return OperationMode::Twiss;

   return OperationMode::Unknown;
}
