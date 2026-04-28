#include <iostream>
#include <vector>
#include <memory>
#include <string>
#include <complex>
#include <cstdio>
#include <mpi.h>
#include <hdf5.h>
#include <hdf5_hl.h>

#include "mesh.h"
#include "constants.h"


// MPI datatype을 C++ 타입에 따라 자동으로 반환
inline MPI_Datatype get_mpi_type(const int*)          { return MPI_INT; }
inline MPI_Datatype get_mpi_type(const long*)         { return MPI_LONG; }
inline MPI_Datatype get_mpi_type(const float*)        { return MPI_FLOAT; }
inline MPI_Datatype get_mpi_type(const double*)       { return MPI_DOUBLE; }
inline MPI_Datatype get_mpi_type(const unsigned*)     { return MPI_UNSIGNED; }
inline MPI_Datatype get_mpi_type(const unsigned long*){ return MPI_UNSIGNED_LONG; }
//inline MPI_Datatype get_mpi_type(const size_t*)       { return MPI_UNSIGNED_LONG; }



// std::vector용 (예: harmony)
template<typename T>
inline MPI_Datatype get_mpi_type(const std::vector<T>*) { 
    return get_mpi_type(static_cast<const T*>(nullptr)); 
}

template<typename T>
void bcast_meta(T* data, int count, int root = 0)
{
   MPI_Datatype mpi_type = get_mpi_type(data);
   MPI_Bcast(data, count, mpi_type, root, MPI_COMM_WORLD);
}

// 2. std::vector 편의 버전 (가장 자주 사용할 버전)
template<typename T>
void bcast_meta(std::vector<T>& vec, int root = 0)
{
    int myrank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    int size = static_cast<int>(vec.size());

    // 크기 broadcast
    MPI_Bcast(&size, 1, MPI_INT, root, MPI_COMM_WORLD);

    if (myrank != root) {
        vec.resize(size);
    }

    if (size > 0) {
        bcast_meta(vec.data(), size, root);   // 실제 데이터 broadcast
    }
}



template<typename T>
void restoreMeta(const std::string& fileName,
                    const std::string& dataName,
                    T *data,
                    int dataCnt);

template<typename T>
void restore_attr_HDF(const std::string& fileName,
                      const std::string& dataName,
                      const std::string& attrName,
                      T *data,
                      int dataCnt);

std::vector<double> readChunk(const std::string& fileName,
                              int species,
                              long long start,
                              long long row_counts,
                              int column,
                              int dataCnt);

//void restoreFieldComp(std::vector<std::vector<cplx>> &data,
//                      const std::string& fileName,
//                      const std::string& dataName,
//                      int h, int sliceN, int N, int subSliceN, int minI);


void restoreParticleHDF(Domain *D,int iteration)
{
   int myrank,nTasks;
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

   char fileNameBuf[100];
   sprintf(fileNameBuf, "particle%d.h5", iteration);
   std::string fileName = fileNameBuf;

   int startI=1;
   int endI=1+D->subSliceN;

   int nSpecies = 0;
   size_t totalCnt = 0;
   int dataCnt = 0;
   if(myrank==0) {
      restoreMeta(fileName,"nSpecies",&nSpecies,1);
      if(nSpecies!=D->nSpecies) { 
         printf("Setting nSpecies is wrong!! restored nSpecies = %d\n",nSpecies); 
         exit(0); 
      }
      restoreMeta(fileName,"numData",&dataCnt,1);
   }
   MPI_Barrier(MPI_COMM_WORLD);
   bcast_meta(&nSpecies, 1);

   const int sliceN = D->sliceN;
   const int minI = D->minI;
   const int subSliceN = D->subSliceN;

   for(int s=0; s<nSpecies; ++s) {
      if(myrank==0) {
         std::string dataName = std::to_string(s);
         restore_attr_HDF(fileName,dataName,"totalCnt",&totalCnt,1);
      }
      bcast_meta(&totalCnt, 1);

      // equally chunk division
      long long rows_per_rank = totalCnt / nTasks;
      long long extra = totalCnt % nTasks;
      long long my_start = myrank * rows_per_rank + std::min((long long)myrank, extra);
      long long my_count = rows_per_rank + (myrank < extra ? 1LL : 0);

      std::vector<double> local_slice = readChunk(fileName,s,my_start,my_count,6,1);
      std::vector<long long> local_minmax(nTasks+1, 0);
      bool reverseScan = true;
      for(long long idx=0; idx<my_count; ++idx) {
         if(D->minI == static_cast<int>(local_slice[idx])) {
            local_minmax[myrank] = idx+my_start;
            std::cout << "minI, myrank=" << myrank << ", idx=" << idx << std::endl;
            reverseScan = false;
            break;
         }
      }
      if(reverseScan == true) {
         for(long long idx=my_count-1; idx>=0 ; --idx) {
            if(D->maxI == static_cast<int>(local_slice[idx])) {
               local_minmax[myrank+1] = idx+my_start;
               std::cout << "maxI, myrank=" << myrank << ", idx=" << idx << std::endl;
               break;
            }
         }
      }
      std::vector<long long> global_minmax(nTasks+1,0);
      MPI_Allreduce(local_minmax.data(),
                   global_minmax.data(),
                   nTasks + 1, 
                   MPI_LONG_LONG,
                   MPI_SUM,
                   MPI_COMM_WORLD);
      global_minmax[nTasks] = static_cast<long long>(totalCnt);

     
      my_start = global_minmax[myrank]; 
      my_count = global_minmax[myrank+1]-my_start; 
      std::vector<double> chunkData = readChunk(fileName,s,my_start,my_count,0,8);
      std::vector<double> chunkSlice = readChunk(fileName,s,my_start,my_count,6,1);

      std::vector<long long> chunkMinMax(subSliceN+1, 0);

      long long startIdx = 0;
      for(int sliceI=0; sliceI<subSliceN; ++sliceI) {
         for(long long idx=startIdx; idx<my_count; ++idx) {
            if(sliceI+D->minI == static_cast<int>(chunkSlice[idx])) {
               chunkMinMax[sliceI+1] = idx+1;
               startIdx = idx;
            }
         }
      }

      for(int sliceI=startI; sliceI<endI; ++sliceI) {  
         auto New = std::make_unique<ptclList>();

         // head[s] is nullptr, the generate new.
         if (D->particle[sliceI].head[s] == nullptr) {
            D->particle[sliceI].head[s] = new ptclHead{};
            D->particle[sliceI].head[s]->pt = nullptr;
         }
         
         New->next = D->particle[sliceI].head[s]->pt;
         D->particle[sliceI].head[s]->pt = New.get();

         int numPtcl = static_cast<int>(chunkMinMax[sliceI-startI+1]-chunkMinMax[sliceI-startI]);
printf("numPtcl = %d\n",numPtcl);

         New->weight = 0.0;
         New->x.resize(numPtcl);
 //        New->y.resize(numPtcl);
 //        New->px.resize(numPtcl);
 //        New->py.resize(numPtcl);
 //        New->theta.resize(numPtcl);
 //        New->gamma.resize(numPtcl);

         New.release();

      }
      
      //if(myrank==3) {
      //   for(size_t i=0; i<subSliceN+1; ++i) {
       //     std::cout << "chunk_minmax[" << i << "]=" << chunkMinMax[i]
       //               << ", my_count=" << my_count 
       //               << std::endl;
       //  }
      //}
      
   }   // End of for(s)



/*
   // restoring Field data (parallel I/O)
   int N = D->nx * D->ny;
   int subSliceN = D->subSliceN;
   int minI = D->minI;
   for (int h = 0; h < D->numHarmony; ++h) {
      int H = D->harmony[h];
      
      std::string dataName = "Ux" + std::to_string(H);
      restoreFieldComp(D->Ux, fileName, dataName, h, D->sliceN, N, subSliceN, minI);

      dataName = "Uy" + std::to_string(H);
      restoreFieldComp(D->Uy, fileName, dataName, h, D->sliceN, N, subSliceN, minI);
   }
   MPI_Barrier(MPI_COMM_WORLD);
*/
   if(myrank==0)
      std::cout << fileName << " is restored." << std::endl;
}

/*
void restoreFieldComp(std::vector<std::vector<cplx>> &data,
                      const std::string& fileName,
                      const std::string& dataName,
                      int h, int sliceN, int N, int subSliceN, int minI)
{
   int myrank, nTasks;
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

   // Parallel file access
   hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
   H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);

   //open file
   hid_t file_id = H5Fopen(fileName.c_str(),H5F_ACC_RDONLY, plist_id);
   H5Pclose(plist_id);

   hid_t dset_id = H5Dopen2(file_id, dataName.c_str(), H5P_DEFAULT);

   // Memory space for this process
   hsize_t count[2] = {static_cast<hsize_t>(subSliceN), static_cast<hsize_t>(N * 2)};
   hid_t memspace = H5Screate_simple(2, count, nullptr);

   // File space hyperslab
   hsize_t offset[2] = {static_cast<hsize_t>(minI), 0};
   hid_t filespace = H5Dget_space(dset_id);
   H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, nullptr, count, nullptr);

   // Read buffer
   std::vector<double> field(static_cast<size_t>(subSliceN) * N * 2);

   // Collective read
   plist_id = H5Pcreate(H5P_DATASET_XFER);
   H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
   
   H5Dread(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, field.data());


   // Data saving
   int startI = 1;
   int endI   = subSliceN + 1;
   size_t idx = 0 ; //minI * N * 2;

   for (int i = startI; i < endI; ++i) {
      for (int j = 0; j < N; ++j) {
         double re = field[idx++];
         double im = field[idx++];
         data[h][i*N + j] = re + I*im;
      }
   }

   H5Pclose(plist_id);
   H5Sclose(memspace);
   H5Sclose(filespace);
   H5Dclose(dset_id);
   H5Fclose(file_id);
}
*/
/*
template<typename T>
void restoreMeta(const std::string& fileName,
                 const std::string& dataName,
                 T *data,
                 int dataCnt)
{
   hid_t file_id = H5Fopen(fileName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

   hid_t dset_id = H5Dopen2(file_id, dataName.c_str(), H5P_DEFAULT);

   // HDF5 type sellection in auto
   hid_t mem_type = H5T_NATIVE_INT;        //default
   if constexpr (std::is_same_v<T, double>) {
      mem_type = H5T_NATIVE_DOUBLE;
   }
   else if constexpr (std::is_same_v<T, float>) {
      mem_type = H5T_NATIVE_FLOAT;
   }
   else if constexpr (std::is_same_v<T, long>) {
      mem_type = H5T_NATIVE_LONG;
   }

   H5Dread(dset_id, mem_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
   
   H5Dclose(dset_id);
   H5Fclose(file_id);
}
*/

std::vector<double> readChunk(const std::string& fileName,
                           int species,
                           long long start,
                           long long row_counts,
                           int column,
                           int dataCnt)
{
   hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
   H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);

   hid_t file_id = H5Fopen(fileName.c_str(), H5F_ACC_RDONLY, plist_id);
   H5Pclose(plist_id);

   std::string dataName = std::to_string(species);
   hid_t dset_id = H5Dopen2(file_id, dataName.c_str(), H5P_DEFAULT); 

   hsize_t offset[2] = {static_cast<hsize_t>(start), column};
   hsize_t count[2]  = {static_cast<hsize_t>(row_counts), dataCnt};
   
   hid_t memspace = H5Screate_simple(2, count, nullptr);
   hid_t filespace = H5Dget_space(dset_id);
   H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, nullptr, count, nullptr);

   std::vector<double> chunk_data(row_counts*dataCnt);

   hid_t xfer_plist = H5Pcreate(H5P_DATASET_XFER);
   H5Pset_dxpl_mpio(xfer_plist, H5FD_MPIO_COLLECTIVE);

   H5Dread(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, xfer_plist, chunk_data.data());

   H5Pclose(xfer_plist);
   H5Sclose(memspace);
   H5Sclose(filespace);
   H5Dclose(dset_id);
   H5Fclose(file_id);

   return chunk_data;
}


template<typename T>
void restore_attr_HDF(const std::string& fileName,
                      const std::string& dataName,
                      const std::string& attrName,
                      T *data,
                      int dataCnt)
{
   int myrank;
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   // 파일 열기 (Read Only)
   hid_t file_id = H5Fopen(fileName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

   // Dataset 열기
   hid_t dataset_id = H5Dopen2(file_id, dataName.c_str(), H5P_DEFAULT);
   if (dataset_id < 0) {
      if (myrank == 0) {
         std::cerr << "Error: Cannot open dataset " << dataName << std::endl;
      }
      H5Fclose(file_id);
      return;
   }

   // Attribute 열기
   hid_t attribute_id = H5Aopen(dataset_id, attrName.c_str(), H5P_DEFAULT);
   if (attribute_id < 0) {
      if (myrank == 0) {
         std::cerr << "Error: Attribute '" << attrName 
                   << "' not found in dataset " << dataName << std::endl;
      }
      H5Dclose(dataset_id);
      H5Fclose(file_id);
      return;
   }

   // HDF5 타입 자동 선택 (save_attr_HDF와 동일)
   hid_t mem_type = H5T_NATIVE_INT; // default
   if constexpr (std::is_same_v<T, double>) {
      mem_type = H5T_NATIVE_DOUBLE;
   }
   else if constexpr (std::is_same_v<T, float>) {
      mem_type = H5T_NATIVE_FLOAT;
   }
   else if constexpr (std::is_same_v<T, long>) {
      mem_type = H5T_NATIVE_LONG;
   }
   else if constexpr (std::is_same_v<T, int>) {
      mem_type = H5T_NATIVE_INT;
   }

   // Attribute 읽기
   herr_t status = H5Aread(attribute_id, mem_type, data);

   if (status < 0 && myrank == 0) {
      std::cerr << "Warning: Failed to read attribute '" << attrName << "'" << std::endl;
   }

   // 정리
   H5Aclose(attribute_id);
   H5Dclose(dataset_id);
   H5Fclose(file_id);
}
