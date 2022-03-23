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


    MPI_File_open(MPI_COMM_SELF, "test.txt",MPI_MODE_CREATE | MPI_MODE_WRONLY,MPI_INFO_NULL,&fh);
    int i = rank + 1;
    std::string s = std::to_string(i) + " hey\n";
    for(int j = 100*rank; j < 100*numprocs; j++);
    std::cout << "EXECUTING " << rank << std::endl;
    MPI_File_set_view(fh, (1+rank)*s.size()*sizeof(char),  MPI_CHAR, MPI_CHAR, "native", MPI_INFO_NULL);
    MPI_File_write(fh,s.c_str(),s.size()*sizeof(char), MPI_CHAR,&status);
    //        fclose(f);
    MPI_File_close(&fh);


    MPI_Finalize();
    return 0;
}