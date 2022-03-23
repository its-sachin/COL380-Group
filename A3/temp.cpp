#include <mpi.h>
#include <omp.h>
#include<bits/stdc++.h>

using namespace std;

void writeAll(int n, std::ofstream &out, int size){
    for(int i =0; i<n ;i++){
        out.write((char*)&i, size);
    }
}

int main(int argc, char **argv) {
    int rank,size;
    
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    vector<string> files = {"level.txt","level_offset.txt","index.txt","indptr.txt","vect.txt"};
    int* sizes = new int[7];

    for (int i = 0; i < files.size(); i++)
    {
        MPI_File input;
        MPI_Offset filesize;
        int maxTextSize = 50,perRankReadSize,lastIndexToInclude,firstIndexToInclude;
        char* buffData;

        string inputFilePath = string(argv[1])+"/"+string(files[i]);
        MPI_File_open(MPI_COMM_WORLD, inputFilePath.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &input);
        MPI_File_get_size(input, &filesize);

        // 0..... perRankReadSize-1 --> rank 0
        // perRankReadSize ....  2*perRankReadSize - 1 ---> rank 1
        //.....
        // n*perRankReadSize ....  2*n*perRankReadSize - 1 (filesize-1)---> rank n
        perRankReadSize = filesize/size;
        
        while (true)
        {
            bool status1 = false;
            bool status2 = false;

            maxTextSize = maxTextSize*2;
            MPI_Offset startReadAt,endReadAt,extendedEnd;
            
            startReadAt = (perRankReadSize)*rank;
            endReadAt = min(((perRankReadSize)*(rank+1))-1,(int)(filesize-1));
            extendedEnd = min((((perRankReadSize)*(rank+1))+maxTextSize)-1,(int)(filesize-1));

            int numBytesToRead = (extendedEnd-startReadAt) + 1;

            buffData = new char[numBytesToRead];

            MPI_File_read_at(input, startReadAt, buffData, numBytesToRead , MPI_CHAR, MPI_STATUS_IGNORE);

            int indexEndReadAT = endReadAt - startReadAt;
            int indexStartReadAt = 0;
            int indexExtendedEnd = extendedEnd - startReadAt;
            
            if(indexEndReadAT!= numBytesToRead-1 && (buffData[indexEndReadAT]==' '||buffData[indexEndReadAT]=='\n')){
                indexEndReadAT++;
            }

            lastIndexToInclude = indexEndReadAT;
            for (int j = indexEndReadAT; j <= min(numBytesToRead-1,indexExtendedEnd); j++)
            {
                lastIndexToInclude = j;
                if(buffData[j]==' '||buffData[j]=='\n'){
                    status1 = true;
                    break;
                }
            }
            if(rank==size-1||indexEndReadAT== numBytesToRead-1){
                status1 = true;
            }
            if(rank!=0){
                for (int j = 0; j < numBytesToRead; j++)
                {
                    firstIndexToInclude = j;
                    if(buffData[j]==' '||buffData[j]=='\n'){
                        status2 = true;
                        break;
                    }
                }
            }else{
                firstIndexToInclude = 0;
                status2 = true;
            }

            if(status1&&status2){
                break;
            }
            
        }
        string s = "";
        // #pragma omp parallel for num_threads(stoi(argv[3])) shared(s, buff) reduction(+: s)
        for(int i = firstIndexToInclude;i<=lastIndexToInclude;i++){
            s+=buffData[i];
        }


        MPI_File_close(&input);

        std::istringstream is( s);
        string myString = files[i];
        string output = string(argv[2])+"/" +(myString.substr(0, myString.size()-3))+"bin";
        std::ofstream outfile(output, std::ios::out | std::ios::binary);
        int writeOffset = 0;
        int *allSizes = new int[size];
        int Val;

        int sizeoftype = sizeof(int);

        if(files[i]=="vect.txt"){
        
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
        sizes[i+2] = writeOffset + allSizes[rank];
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