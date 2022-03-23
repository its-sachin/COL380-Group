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
        MPI_File input;
        int maxTextSize = 5;
        
        string fullFilePath = string(argv[1])+"/"+string(files[f]);
        
        MPI_File_open(MPI_COMM_WORLD, fullFilePath.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &input);
        // MPI_File_open(MPI_COMM_WORLD, "vect.bin", MPI_MODE_CREATE |MPI_MODE_WRONLY, MPI_INFO_NULL, &outFile);
        

        
        MPI_Offset toStart,toEnd,filesize,toEndNoOffset;
        int curRankSize;
        char *buff;
        bool status = false;
        int start = 0;
        int realEnd;
        while (!status)
        {
            //free(buff);
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

        // cout<<"At rank "<<rank<<endl;
        string s = "";
        // #pragma omp parallel for num_threads(stoi(argv[3])) shared(s, buff) reduction(+: s)
        for(int i = start;i<realEnd;i++){
            s+=buff[i];
        }

        // char* buff1 = new char[60];
        // string s1 = "";
        // for (int i = 0; i < 60; i++)
        // {
        //     buff1[i] = '0'+i;
        //     // cout<<buff1[i]<<" ";
        // }

        // #pragma omp parallel for num_threads(4) shared(s1, buff1) reduction(+: s1)
        // for (int i = 0; i < 60; i++)
        // {
        //     s1+=buff1[i];
        // }

        MPI_File_close(&input);

        std::istringstream is( s);
        string myString = files[f];
        string output = string(argv[2])+"/" +(myString.substr(0, myString.size()-3))+"bin";
        std::ofstream outfile(output, std::ios::out | std::ios::binary);
        int writeOffset = 0;
        int *allSizes = new int[size];
        int Val;

        int sizeoftype = sizeof(int);

        if(files[f]=="vect.txt"){
        
            sizeoftype = sizeof(double);
            vector<double> nums;
            double n;
            while( is >> n ) {
                nums.push_back(n);
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
            int n;
            while( is >> n ) {
                nums.push_back(n);
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
    }

    if(rank == size - 1){
        sizes[6] /= sizes[2];
        sizes[0] = sizes[2];

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