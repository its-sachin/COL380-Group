#include <stdio.h>
#include <mpi.h>
#include <omp.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include<bits/stdc++.h>

using namespace std;

void writeAll(int n, std::ofstream &out, int size){
    for(int i =0; i<n ;i++){
        out.write((char*)&i, size);
    }
}

    
int main(int argc, char **argv) {
    // cout<<argv[1]<<endl;
    

    
    //criticalData[6] = 0;
    int rank, size;
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    vector<string> files = {"level.txt","level_offset.txt","index.txt","indptr.txt","vect.txt"};
    int* sizes = new int[7];

    for (int f = 0; f < files.size(); f++)
    {       
        cout  << rank <<  ": Reading " << files[f] << endl;
        string fullFilePath = string(argv[1]) + "/" + files[f];

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

        // if(rank > 0)start += 1;
        infile.seekg(start, std::ios::beg);

        if(rank > 0){
            char curr = '\0';
            while(!(curr == '\n' or curr == ' ')){
                infile.read((char*)&curr, 1);
                start += 1;
            }
        }


        string myString = files[f];
        string output = string(argv[2])+"/" +(myString.substr(0, myString.size()-3))+"bin";
        std::ofstream outfile(output, std::ios::out | std::ios::binary);
        int writeOffset = 0;
        int *allSizes = new int[size];

        int sizeoftype = sizeof(int);

        if(files[f]=="vect.txt"){
        
            sizeoftype = sizeof(double);
            vector<double> nums;
            double a;
            while(start <= end && infile >> a){  
                start += to_string(a).size() + 1;
                nums.push_back(a);
            }

            allSizes[rank] = nums.size();

            MPI_Allgather(&allSizes[rank], 1, MPI_INT, allSizes,1, MPI_INT, MPI_COMM_WORLD);

            for (int i = 0; i < rank; i++)
            {
                writeOffset+=allSizes[i];
            }
            
            if(rank == size -1){
                writeAll(writeOffset, outfile, sizeoftype);
            }


            MPI_Barrier(MPI_COMM_WORLD);
            outfile.seekp(sizeoftype*writeOffset, std::ios::beg);

            for(int i=0; i<nums.size(); i++){
                outfile.write((char*)&nums[i], sizeoftype);
            }
            outfile.close();
            
        } else{ 

            vector<int> nums;
            int a;
            while(start <= end && infile >> a){  
                start += to_string(a).size() + 1;
                nums.push_back(a);
            }

            allSizes[rank] = nums.size();
            
            MPI_Allgather(&allSizes[rank], 1, MPI_INT, allSizes,1, MPI_INT, MPI_COMM_WORLD);
            
            for (int i = 0; i < rank; i++)
            {
                writeOffset+=allSizes[i];
            }

            if(rank == size -1 ){
                writeAll(writeOffset, outfile, sizeoftype);
            }

            MPI_Barrier(MPI_COMM_WORLD);
            outfile.seekp(sizeoftype*writeOffset, std::ios::beg);

            for(int i=0; i<nums.size(); i++){
                outfile.write((char*)&nums[i], sizeoftype);
            }

            outfile.close();
            // Val = nums[0];

        }
        sizes[f+2] = writeOffset + allSizes[rank];
        infile.close();
    }
    

    if(rank == size - 1){
        sizes[6] /= sizes[2];
        sizes[0] = sizes[2];

        cout << sizes[2] << " " << sizes[6] << endl;

        ifstream input;

        files = {"ep.txt","max_level.txt"};
        for(int f=0; f<2; f++){
            string fullFilePath = string(argv[1])+"/"+string(files[f]);
            input.open(fullFilePath);
            int k;
            input >>k;
            sizes[f+1] = k;
            input.close();
        }

        string output = string(argv[2])+"/sizes.bin";
        std::ofstream outfile(output, std::ios::out | std::ios::binary);
        
        for(int i=0; i<7; i++){
            outfile.write((char*)&sizes[i], sizeof(int));
        }

        outfile.close();


    }

    
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