/* This is an interactive version of cpi */
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


int main(int argc,char *argv[])
{

    int  namelen, numprocs, rank;
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&numprocs);    
    MPI_Get_processor_name(processor_name,&namelen);
    MPI_Status status;

    MPI_File fh;


    if (rank == 0) {
        MPI_File_open(MPI_COMM_SELF, "test.txt",MPI_MODE_CREATE | MPI_MODE_WRONLY,MPI_INFO_NULL,&fh);
        //FILE* f = fopen("test.txt","wb+");
        //if(f==NULL){
        //printf("failed to open file\n");exit(1);
        //}
        for (int i=0; i < 5; i++){
            char buf[42];
            //fprintf(f,"%d \n",i);
            snprintf(buf,42,"%d \n",i);
            MPI_File_write(fh,buf,strlen(buf), MPI_CHAR,&status);
        }
        //        fclose(f);
        MPI_File_close(&fh);
    }
    else {
        // do nothing
    }


    MPI_Finalize();
    return 0;
}
