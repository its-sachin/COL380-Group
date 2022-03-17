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

    int sendarray[5];
    for(int i=0; i<5; i++) sendarray[i] = rank; 
    int *rbuf; 

    if ( rank == 0) { 
       rbuf = new int[5*size];
       } 
    MPI_Gather( sendarray, 5, MPI_INT, rbuf, 5, MPI_INT, 0, MPI_COMM_WORLD); 
    
    if(rank ==0){
        for(int i=0; i<size; i++){
            for(int j=0; j<5; j++){
                cout<<rbuf[i*5+j]<<" ";
            }
            cout<<endl;
        }
    }

    MPI_Finalize();
}