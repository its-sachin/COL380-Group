#include<bits/stdc++.h>
#include <mpi.h>

using namespace std;


int main(int argc, char* argv[]){

    int rank, size;
    const int INT = sizeof(int);
    const int DOUBLE = sizeof(double);

    //Starting MPI pipeline
    MPI_Init(NULL, NULL);
    
    // Extracting Rank and Processor Count
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // int sendarray[5];
    // for(int i=0; i<5; i++) sendarray[i] = rank; 
    // int *rbuf; 

    // if ( rank == 0) { 
    //    rbuf = new int[5*size];
    //    } 
    // MPI_Gather( sendarray, 5, MPI_INT, rbuf, 5, MPI_INT, 0, MPI_COMM_WORLD); 
    
    // if(rank ==0){
    //     for(int i=0; i<size; i++){
    //         for(int j=0; j<5; j++){
    //             cout<<rbuf[i*5+j]<<" ";
    //         }
    //         cout<<endl;
    //     }
    // }

    int d = 3;
    // MPI_Datatype row;
    // MPI_Type_contiguous(d, MPI_INT, &row);
    // MPI_Type_commit(&row);

    int sizee = sizeof(double);
    MPI_Datatype a = MPI_DOUBLE;


    MPI_File fs;
    MPI_File_open(MPI_COMM_WORLD, "vect.bin",MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fs);
    MPI_File_set_view(fs, rank*3*sizee, a, a, "native", MPI_INFO_NULL);

    double aa[3] = {1,2,4};
    MPI_File_write(fs, aa, sizee, a, MPI_STATUS_IGNORE);
    MPI_File_write(fs, aa, sizee, a, MPI_STATUS_IGNORE);
    MPI_File_write(fs, aa, sizee, a, MPI_STATUS_IGNORE);
    // vector<int> buff[3] = {1,2,3};
    // for(int i=0; i<3;i++)buff[i]+=rank;

    // if(rank==0){
    //     MPI_File_write(fs, buff, 3*sizee, a, MPI_STATUS_IGNORE);
    //     MPI_File_write(fs, buff, 3*sizee, a, MPI_STATUS_IGNORE);
        
    // }

    // MPI_Barrier(MPI_COMM_WORLD);

    // MPI_File_write(fs, buff, 3*sizee,a, MPI_STATUS_IGNORE);
    // int* arr = new int[6];

    // int* temp = new int[3];
    // MPI_File_read(fs, temp, 1, row, MPI_STATUS_IGNORE);
    MPI_File_close(&fs);

    // cout << endl;
    // for(int i=0; i<3; i++)cout << rank << " temp[" << i << "] : " << temp[i] << endl;

    // // MPI_Allgather(temp, 1, row, arr, 1, row, MPI_COMM_WORLD);
    // MPI_Allgather(temp, 1, row, arr, 1, row, MPI_COMM_WORLD);

    // for(int i=0; i<2; i++){
    //     for(int j=0; j<3; j++){
    //         cout << rank << " arr[" << i << "][" << j << "] : " << arr[2*i + j] << endl;
    //     }
    // }

    MPI_Finalize();
}