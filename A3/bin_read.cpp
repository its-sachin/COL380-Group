#include<bits/stdc++.h>
#include <mpi.h>

using namespace std;

int getBuffSize(int len, int size, int rank){
    int buffsize = len/size;
    if(rank == size-1){
        buffsize = len - (size-1)*buffsize;
    }
    return buffsize;
}

void readMPI(const char *filename, int* arr, int len, int size, int rank, MPI_File &fs){
    
    int blockSize =  sizeof(int);
    int buffsize = getBuffSize(len, size, rank);
    int disp = rank*(len/size)*blockSize;

    MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_RDONLY, MPI_INFO_NULL, &fs);
    MPI_File_set_view(fs, disp, MPI_INT, MPI_INT, "native", MPI_INFO_NULL);

    MPI_File_read(fs, arr, buffsize, MPI_INT, MPI_STATUS_IGNORE);
    MPI_File_close(&fs);
}

void gatherMPI(const char *filename, int* arr, int len, int size, int rank, MPI_File &fs){
    int* temp = new int[getBuffSize(len, size, rank)];
    readMPI(filename, temp, len, size, rank, fs);

    MPI_Allgather(temp, len/size, MPI_INT, arr, len/size, MPI_INT, MPI_COMM_WORLD);

    int left = getBuffSize(len, size, size - 1) - len/size;

    if(left > 0){
        for(int i=0; i<left; i++){
            temp[i] = temp[len/size + i];
        }
    }

    MPI_Bcast(temp, left, MPI_INT, size-1, MPI_COMM_WORLD);
        for(int i=0; i<left; i++){
            arr[size*(len/size) + i] = temp[i];
        }

    delete (temp);
}

void gatherMPI(const char *filename, double* arr, int n, int d, int size, int rank, MPI_File &fs){

    MPI_Datatype row;
    MPI_Type_contiguous(d, MPI_DOUBLE, &row);
    MPI_Type_commit(&row);

    int blockSize =  sizeof(double)*d;
    int disp = rank*(n/size)*blockSize;

    MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_RDONLY, MPI_INFO_NULL, &fs);
    MPI_File_set_view(fs, disp, row, MPI_DOUBLE, "native", MPI_INFO_NULL);

    int start = rank*(n/size);
    int end = (rank+1)*(n/size);
    if(rank == size-1){
        end = n;
    }

    double *temp = new double[(n + end-start)*d];
    MPI_File_read(fs, temp, end-start, row, MPI_STATUS_IGNORE);

    MPI_File_close(&fs);
    MPI_Allgather(temp, n/size, row, arr, n/size, row, MPI_COMM_WORLD);

    int left = n - (size-1)*(n/size) - n/size;
    // cout << "left: " << left<<endl;/
    if(left > 0){
        for(int i=0; i<left; i++){
            for(int j=0; j<d; j++)
                temp[i*d + j] = temp[(n/size + i)*d + j];
        }

        MPI_Bcast(temp, left, row, size-1, MPI_COMM_WORLD);
        for(int i=0; i<left; i++){
            for(int j=0; j<d; j++)
                arr[(size*(n/size) + i)*d + j] = temp[i*d + j];
        }
    }

    delete (temp);
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

    readMPI("outFolder/sizes.bin", read, 7, 1, 0, fs);

    int n = read[0]; // size of levels 
    int ep = read[1]; // ep file value
    int max_level = read[2]; // max level value
    int s_level_offset = read[3]; // level offset size
    int s_index = read[4]; // index size 
    int s_indptr = read[5]; // indptr size
    int d = read[6]; // vect bin per line size

    cout << rank << ": "<<"n: "<<n<<endl;
    cout << rank << ": "<< "ep: "<<ep<<endl;
    cout << rank << ": "<<"max level: "<<max_level<<endl;
    cout << rank << ": "<<"s_level_offset: "<<s_level_offset<<endl;
    cout << rank << ": "<<"s_index: "<<s_index<<endl;
    cout << rank << ": "<<"s_indptr: "<<s_indptr<<endl;
    cout << rank << ": "<<"d: "<<d<<endl;

    int* indptr = new int[s_indptr];
    gatherMPI("outFolder/indptr.bin", indptr, s_indptr, size, rank, fs);

    cout << rank << ": "<<"indptr: "<<endl;
    cout<< rank << ": ";
    for(int i=0; i<s_indptr; i++){
        cout <<  indptr[i] << " ";
    }
    cout <<endl;

    int* index = new int[s_index];
    gatherMPI("outFolder/index.bin", index, s_index, size, rank, fs);

    cout<< rank << ": "<<"index: "<<endl;
    cout<< rank << ": ";
    for(int i=0; i<s_index; i++){
        cout <<  index[i] << " ";
    }
    cout <<endl;

    int* level_offset = new int[s_level_offset];
    gatherMPI("outFolder/level_offset.bin", level_offset, s_level_offset, size, rank, fs);

    cout<< rank << ": "<<"level_offset: "<<endl;
    cout<< rank << ": ";
    for(int i=0; i<s_level_offset; i++){
        cout <<  level_offset[i] << " ";
    }
    cout <<endl;

    // double* vect = new double[n*d];
    // gatherMPI(argv[1], vect, n, d, size, rank, fs);

    // cout<< rank << ": "<<"vect: "<<endl;
    // for(int i=0; i<n; i++){
    //     cout<< rank << ": ";
    //     for(int j=0; j<d; j++){
    //         cout<<vect[i*d + j]<<" ";
    //     }
    //     cout<<endl;
    // }

    double* q = new double[n*d];
    gatherMPI("outFolder/vect.bin", q, n, d, size, rank, fs);

    cout<< rank << ": "<<"q: "<<endl;
    for(int i=0; i<n; i++){
        cout<< rank << ": ";
        for(int j=0; j<d; j++){
            cout<<q[i*d + j]<<" ";
        }
        cout<<endl;
    }


    MPI_Finalize();
}