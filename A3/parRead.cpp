#include <stdio.h>
#include <mpi.h>
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
    
    MPI_File input;
    MPI_File outFile;
    int rank, size;
    int maxTextSize = 5;
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    MPI_File_open(MPI_COMM_WORLD, argv[1], MPI_MODE_RDONLY, MPI_INFO_NULL, &input);
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

     cout<<"At rank "<<rank<<endl;
    string str = "";
    for(int i = start;i<realEnd;i++){
        str+=buff[i];
        //cout<<buff[i];
    }

    //cout<<str<<endl;

    //cout<<str<<endl;
    //cout<<"At rank "<<rank<<endl;
    std::istringstream is( str);

    
    vector<double> nums;
    double n;
    while( is >> n ) {
        nums.push_back(n);
        //cout<<n<<endl;
    }
    int *allSizes = (int *)malloc(sizeof(int) * size);

    //vector<int> allSizes(size);

    allSizes[rank] = nums.size();

    
    //cout<<"Real rank size on rank "<<rank<<" is "<<allSizes[rank]<<endl;
    MPI_Allgather(&allSizes[rank], 1, MPI_INT, allSizes,1, MPI_INT, MPI_COMM_WORLD);
    
    //MPI_File_set_view(outFile, 0, MPI_DOUBLE, MPI_DOUBLE,"native", MPI_INFO_NULL);
    
    int writeOffset = 0;
    int full = 0;
     
    for (int i = 0; i < rank; i++)
    {
        cout << allSizes[i] << " ";
        writeOffset+=allSizes[i];
        full += allSizes[i];
    }

    std::ofstream outfile("vect.bin", std::ios::out | std::ios::binary);
    if(rank == size -1 ){
        writeAll(full, outfile, sizeof(double));
    }


    MPI_Barrier(MPI_COMM_WORLD);
    outfile.seekp(sizeof(double)*writeOffset, std::ios::beg);

    for(int i=0; i<nums.size(); i++){
        outfile.write((char*)&nums[i], sizeof(double));
    }

    outfile.close();

    // MPI_File_seek(outFile, writeOffset*sizeof(double), MPI_SEEK_SET );
    // if(rank==0){
    //     double* a = (double*)malloc( full*sizeof(double));
    //     MPI_File_write(outFile, a, full*sizeof(double), MPI_DOUBLE, MPI_STATUS_IGNORE);
    //     MPI_File_seek(outFile, 0, MPI_SEEK_SET );
    // }

    // MPI_Barrier(MPI_COMM_WORLD);

    // double* t = new double[nums.size()];

    // for(int i=0; i<nums.size();i++){
    //     t[i] = nums[i];
    // }
    // MPI_File_write(outFile, t, nums.size()*sizeof(double), MPI_DOUBLE, MPI_STATUS_IGNORE);

    

    // MPI_Barrier(MPI_COMM_WORLD);

    // for (int i = 0; i < nums.size(); i++)
    // {
    //     //cout<<nums[i]<<endl;
    //     // cout << "rank: " <<rank << " offset: " <<   (writeOffset + i)*(__SIZEOF_DOUBLE__) << " val: " << nums[i] << endl; 
    //     MPI_File_write_at(outFile, (writeOffset + i)*(__SIZEOF_DOUBLE__), &t[i],sizeof(double), MPI_DOUBLE, MPI_STATUS_IGNORE);
    // }
    
    //double data = 1.53;
    // char* data;
    // data = (char*)malloc( ((1)*sizeof(char)));
    // data[0] = 'a';
    // if(rank==0){
    //     MPI_File_write_at(outFile, 0, &data,sizeof(char), MPI_CHAR, MPI_STATUS_IGNORE);
    //     // char* newData;
    //     // data = (char*)malloc( ((1)*sizeof(char)));
    //     // MPI_File_read_at(outFile, (sizeof(char))*writeOffset, &newData,sizeof(char), MPI_CHAR, MPI_STATUS_IGNORE);
    //     // cout<<newData<<endl;
    //     //MPI_File_write_at(MPI_File fh, MPI_Offset offset, void *buf,int count, MPI_Datatype datatype, MPI_Status *status)
    // }

    // MPI_File_write(outFile, &data, sizeof(double), MPI_DOUBLE, MPI_STATUS_IGNORE);
    // MPI_File_set_view(outFile, sizeof(double), MPI_DOUBLE, MPI_DOUBLE,"native", MPI_INFO_NULL);
    
    //  data = 1.54;
    // MPI_File_write(outFile, &data, sizeof(double), MPI_DOUBLE, MPI_STATUS_IGNORE);
    // double newData;
    // MPI_File_read(outFile, &newData, sizeof(double), MPI_DOUBLE, MPI_STATUS_IGNORE);
    // cout<<newData;
    // for (int i = 0; i < size; i++)
    // {
    //     cout<<i<<" Received "<<allSizes[i]<<" in rank "<<rank<< " ";
    // }
    // cout << endl;
    
    // for (auto &&i : allSizes)
    // {
    //     cout<<i<<endl;
    // }
    //cout<<endl;
    
    // if(rank==0){
    //     int BUFSIZE = 0;
    //     MPI_Offset writeAt = BUFSIZE * sizeof(double);
    //     // for (int i = 0; i < nums.size(); i++)
    //     // {
    //         //cout<<i<<endl;
            
    //         // cout<<ierr<<endl;
    //         // MPI_File_write(outFile, reinterpret_cast<char*>( &nums[i] ), BUFSIZE, MPI_DOUBLE, MPI_STATUS_IGNORE);
    //         BUFSIZE++;
    //     // }    
    // }

    //vector<int> vec;
    
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
    // MPI_File_close(&outFile);
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