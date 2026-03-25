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
BeamMode whatBeamMode(char *str);
bool findBeamLoadParameters(int rank,LoadList& LL,Domain *D,const char *input);

void parameterSetting(Domain *D,const char *input)
{
   LoadList LL{};

   int myrank=0, nTasks=1;
   //MPI_Status status;

   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);     
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);     

   int i=0, rank=1;
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



   // Beam parameter setting
   while(findBeamLoadParameters(rank, LL, D,input))
   {
      D->loadList.push_back(LL);
      std::cout << "loadType=" << static_cast<int>(LL.type) 
                << "rank=" << rank    
                << "\n";
      rank ++;
      LL = LoadList {};
   }
   D->nSpecies = rank-1;
   std::cout << "nSpecies=" << D->nSpecies << "\n";


   // extra setting
   D->nx=1;
   D->ny=1;

   D->ks=energy*2*M_PI/(plankH*velocityC);
   D->lambda0=2*M_PI/D->ks;

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

bool findBeamLoadParameters(int rank,LoadList& LL,Domain *D,const char *input)
{
   char str[100],name[100];
   bool fail = false;

   int myrank, nTasks;
   MPI_Status status;

   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   if(FindParameters("EBeam",rank,"load_type",input,str)) LL.type = whatBeamMode(str);
   else LL.type=BeamMode::Unknown;

   if (LL.type != BeamMode::Unknown) {
      if(FindParameters("EBeam",1,"beamlets_in_bucket",input,str)) LL.numBeamlet=atoi(str);
      else  { printf("In [EBeam], beamlets_in_bucket=? [ea/bucket].\n"); fail=true;  }
      if(FindParameters("EBeam",1,"number_in_beamlet",input,str)) LL.numInBeamlet=atoi(str);
      else  { printf("In [EBeam], number_in_beamlet=? [ea/beamlet].\n"); fail=true;  }


      switch (LL.type)  {
      case BeamMode::Unknown:
         std::cout << "In [EBeam], wrong load_type !!!\n";
         fail = true;
         break;
      case BeamMode::Polygon :
         if(FindParameters("EBeam",rank,"z_nodes",input,str)) LL.znodes=atoi(str);
         else  { printf("in [EBeam], z_nodes=?\n");  fail=true;   }
         if(LL.znodes>0)  {
             LL.zpoint.resize(LL.znodes);
             LL.zn.resize(LL.znodes);
             for(int i=0; i<LL.znodes; i++)  {
                sprintf(name,"z%d",i);
                if(FindParameters("EBeam",rank,name,input,str)) LL.zpoint[i] = atof(str);
                else  { printf("z%d should be defined.\n",i);  fail=1; }

                sprintf(name,"z_n%d",i);
                if(FindParameters("EBeam",rank,name,input,str)) LL.zn[i] = atof(str);
                else { printf("z_n%d should be defined.\n",i);  fail=1; }
             }
          }
         break;
      case BeamMode::Gaussian:
        // Gaussian 모드 처리 코드를 여기에 작성하세요
        // 예: LL.sigma = ... 등
         break;
      }
   }

   return (LL.type != BeamMode::Unknown);

   if(fail) {
      MPI_Finalize();
      std::exit(1);
   }

}


BeamMode whatBeamMode(char *str)
{
   if (strstr(str,"Polygon") != nullptr)
      return BeamMode::Polygon;
   else if (strstr(str,"Gaussian") != nullptr)
      return BeamMode::Gaussian;
   else 
      return BeamMode::Unknown;
}



OperationMode whatOperMode(char *str)
{
   //if (!str || !*str) {
   //   return OperationMode::Unknown;
   //}
   if (strstr(str,"Static") != nullptr)
      return OperationMode::Static;
   else if (strstr(str,"Time_Dependent") != nullptr)
      return OperationMode::Time_Dependent;
   else if (strstr(str,"Twiss") != nullptr)
      return OperationMode::Twiss;
   else 
      return OperationMode::Unknown;
}
