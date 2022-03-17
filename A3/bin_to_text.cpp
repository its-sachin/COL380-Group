#include<bits/stdc++.h>
#include <mpi.h>

using namespace std;

int getBuffSize(int len, int size, int rank){
    int buffsize = len/size;
    // if(rank == size-1){
    //     buffsize = len - (size-1)*buffsize;
    // }
    return buffsize;
}

void readMPI(const char *filename, int* arr, int len, int size, int rank, MPI_File &fs){
    
    int blockSize =  sizeof(int);
    int buffsize = getBuffSize(len, size, rank);
    int disp = rank*buffsize*blockSize;
    MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_RDONLY, MPI_INFO_NULL, &fs);
    MPI_File_set_view(fs, disp, MPI_INT, MPI_INT, "native", MPI_INFO_NULL);

    MPI_File_read(fs, arr, buffsize, MPI_INT, MPI_STATUS_IGNORE);
    MPI_File_close(&fs);
}


int main(int argc, char* argv[]){

    int rank, size;
    const int INT = sizeof(int);
    const int DOUBLE = sizeof(double);

    //Starting MPI pipeline
    MPI_Init(NULL, NULL);
    
    // Extracting Rank and Processor Count
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int read[7];
    
    MPI_File fs;

    readMPI("dummy_bin/sizes.bin", read, 7, 1, 0, fs);

    int n = read[0];
    int ep = read[1];
    int max_level = read[2];
    int s_level_offset = read[3];
    int s_index = read[4];
    int s_indptr = read[5];
    int d = read[6];

    int* indp = new int(s_indptr);
    readMPI("dummy_bin/indptr.bin", indp, s_indptr, size, rank, fs);

    int* indptr = new int(s_indptr);
    MPI_Gather(indp, s_indptr, MPI_INT, indptr, s_indptr/size, MPI_INT, 0, MPI_COMM_WORLD);
    
    // double vect[n][d];
    // for(int i=0; i<n; i++){
    //     for(int j=0; j<d; j++){
    //         f.read((char*)&vect[i][j], sizeof(double));
    //     }
    // }

    // f.close();

    cout << "ep: "<<ep<<endl;
    cout <<"max level: "<<max_level<<endl;
    cout <<"n: "<<n<<endl;
    cout <<"d: "<<d<<endl;
    cout <<"indptr: "<<endl;
    for(int i=0; i<s_indptr; i++){
        cout << indptr[i] << " ";
    }
    cout <<endl;

    // cout<<"vect: "<<endl;
    // for(int i=0; i<n; i++){
    //     for(int j=0; j<d; j++){
    //         cout<<vect[i][j]<<" ";
    //     }
    //     cout<<endl;
    // }

    MPI_Finalize();
}