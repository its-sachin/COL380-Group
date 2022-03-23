/* This is an interactive version of cpi */
#include <mpi.h>
#include<bits/stdc++.h>

using namespace std;
int main(int argc,char *argv[])
{
    int rank, size;
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    string fullFilePath = string(argv[1]) + "/vect.txt";
    MPI_File input;
    MPI_File_open(MPI_COMM_WORLD, fullFilePath.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &input);
    MPI_Offset filesize;
    MPI_File_get_size(input, &filesize);
    MPI_File_close(&input);


    ifstream infile;
    infile.open(fullFilePath, std::ios::in);

    long long fsize = (long long)filesize;
    long long start = (rank*fsize)/size;
    long long end = ((rank+1)*fsize)/size;

    if(rank > 0)start += 1;
    infile.seekg(start, std::ios::beg);

    if(rank > 0){
        char curr = '\0';
        while(!(curr == '\n' or curr == ' ')){
            infile.read((char*)&curr, 1);
            start += 1;
        }
    }
    // infile.seekg(start, std::ios::beg);

    vector<double> v;
    bool toadd = (rank == 0);
    while(start <= end){
        double a;
        infile >> a;
        start += to_string(a).size() + 1;
        v.push_back(a);
    }

    infile.close();

    cout << rank << ": ";
    for(int i=0; i<v.size(); i++)cout << v[i] << " ";
    cout << endl;

    MPI_Finalize();
    return 0;
}