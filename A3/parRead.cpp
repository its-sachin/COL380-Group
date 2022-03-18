#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
    
using namespace std;
    
int main(int argc, char **argv) {
    
    MPI_File input;
    int rank, size;
    int maxTextSize = 5;
    
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    MPI_File_open(MPI_COMM_WORLD, argv[1], MPI_MODE_RDONLY, MPI_INFO_NULL, &input);

    
    MPI_Offset toStart,toEnd,filesize,toEndNoOffset;
    int curRankSize;
    char *buff;
    // buff = (char*)malloc( (2)*sizeof(char));
    // MPI_File_read_at_all(input, 0, buff, 1 , MPI_CHAR, MPI_STATUS_IGNORE);
    // cout<<buff<<endl;
    bool status = false;
    int start = 0;
    int realEnd;
    while (!status)
    {
        free(buff);
        status = false;
        //status = true;
        maxTextSize = maxTextSize*2;
        //cout<<maxTextSize<<endl;
        MPI_File_get_size(input, &filesize);
        //cout<<" filesize "<<filesize<<endl;
        curRankSize = filesize/size;

        toStart = rank*curRankSize;
        toEnd = (maxTextSize+((rank+1)*curRankSize)) - 1;
        toEndNoOffset =((rank+1)*curRankSize) - 1;;
        if(toEnd>=filesize-1){
            status = true;
            toEnd = filesize - 1;
        }
        

        buff = (char*)malloc( ((toEnd  -toStart + 1) + 1)*sizeof(char));

        
        MPI_File_read_at_all(input, toStart, buff, (toEnd  -toStart + 1) , MPI_CHAR, MPI_STATUS_IGNORE);

        
        buff[(toEnd  -toStart + 1)] = '\0';


        cout<<"At rank "<<rank<<endl;

        //cout<<"mybuff \n "<<buff<<"\n mybuff end"<<endl;

        //cout<<toStart<<" "<<(toEndNoOffset)<<" "<<toEnd<<endl;
        //cout<<toStart<<" "<<(toEndNoOffset-toStart)<<" "<<toEnd  -toStart<<endl;
        //cout<<buff<<endl;
        realEnd = (toEnd  -toStart + 1);
        for(int i = toEndNoOffset-toStart + 1;i<((toEnd  -toStart + 1));i++){
            if(buff[i]==' '||buff[i]=='\n'){
                //cout<<"is SPace "<<rank<<endl;
                realEnd = i;
                status = true;
                break;
            }
        }
        //cout<<realEnd<<endl;
        if(rank==size-1){
            status = true;
        }

        start = 0;
        if(rank!=0){
            for(int i = 0;i<((toEnd  -toStart + 1));i++){
                //cout<<i<<" "<<buff[i]<<endl;
                if(buff[i]==' '||buff[i]=='\n'){
                    start = i+1;
                    break;
                }
            }
        }
    }
    for(int i = start;i<realEnd;i++){
        cout<<buff[i];
    }

    cout<<endl;
    //cout<<"End at rank "<<rank<<endl;    
    // if (ierr) {
    //     if (rank == 0) fprintf(stderr, "%s: Couldn't open file %s\n", argv[0], argv[1]);
    //     MPI_Finalize();
    //     exit(2);
    // }
    
    // ierr = MPI_File_open(MPI_COMM_WORLD, argv[2], MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &out);
    // if (ierr) {
    //     if (rank == 0) fprintf(stderr, "%s: Couldn't open output file %s\n", argv[0], argv[2]);
    //     MPI_Finalize();
    //     exit(3);
    // }
        
    MPI_File_close(&input);
    // MPI_File_close(&out);
    
    MPI_Finalize();
    return 0;
}



// 10
// 10
// 10
// 10
// At rank 0
// -0.029495 0.030087 -0.013817 -0.017568 0.045015
// 0.004782 0.032958 -0.003865
// End at rank 0
// At rank 1
// 0.003602 0.022957
// 0.004957 0.030688 -0.036196 -0.010938 0.008610
// -0.028899
// End at rank 1
// At rank 2
// 0.061895 -0.045630 0.000977 0.053625
// 0.026941 -0.000656 -0.028331
// End at rank 2
// At rank 3
// 0.038961 0.003586
// -0.014396 0.033023 -0.012421 -0.027277 0.017062
// End at rank 3