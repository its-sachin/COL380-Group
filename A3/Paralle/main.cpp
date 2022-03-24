#include<bits/stdc++.h>
#include <mpi.h>
#include <omp.h>
#include <chrono>
// #include 

using namespace std;
using namespace std::chrono;


typedef pair<int, double> pid;

struct Less
{
   bool operator()(const pid& a, const pid& b){
       return a.second > b.second;
   }
};

struct More
{
   bool operator()(const pid& a, const pid& b){
       return a.second < b.second;
   }
};


typedef priority_queue<pid, vector<pid>, Less> min_pq;
typedef priority_queue<pid, vector<pid>, More> max_pq;



void writeAll(int n, std::ofstream &out, int size){
    for(int i =0; i<n ;i++){
        out.write((char*)&i, size);
    }
}

double dotProduct(double* a, int ai, double* b, int bi, int n)
{
    double sum = 0;
    for(int i = 0; i < n; i++)
        sum += a[i + ai*n] * b[i + bi*n];
    return sum;
}

double cosine_dist(double* a, int ai, double* b, int bi, int n)
{
    double dot = dotProduct(a, ai, b, bi, n);
    double mod_a = dotProduct(a, ai, a, ai, n);
    double mod_b = dotProduct(b, bi, b, bi, n);
    return 1.0 - dot / (sqrt(mod_a) * sqrt(mod_b));
}



// --------------------------------------------------

void searchLayer(double* q, int u, int d, max_pq& top_k, int* indptr, int* level_offset, int lc, 
                int* index, unordered_set<int>& visited, double* vect, int k)
{
    min_pq candidats;
    max_pq temp = top_k;
    while(!temp.empty()){
        candidats.push(temp.top());
        temp.pop();
    }

    while(candidats.size() > 0){
        int ep = candidats.top().first;
        candidats.pop();

        // cout << "ep = " << ep << " Size = " << candidats.getSize()<< endl;

        int start = indptr[ep] + level_offset[lc];
        int end = indptr[ep] + level_offset[lc+1];

        // cout<< "start = " << start << " end = " << end << endl;
        for(int i = start; i < end; i++){

            int px = index[i];
            // cout << "px = " << px << endl;

            if(visited.find(px) != visited.end() || px == -1) continue;
            visited.insert(px);
            double dist = cosine_dist(vect, px, q, u, d);

            if(dist > top_k.top().second && top_k.size() >= k) continue;
            
            top_k.push(pid(px,dist));
            if(top_k.size() > k){
                top_k.pop();
            }
            candidats.push(pid(px,dist));
        }
    }
} 

string queryHNSW(double* q, int u, int d, int k, int ep, int* indptr, int* index, int* level_offset, int max_level, double* vect)
{
    max_pq top_k;
    top_k.push(pid(ep,cosine_dist(vect, ep, q, u, d)));
    // cout << "pushing " << ep << " " << top_k.getMax() << endl;
    unordered_set<int> visited;
    visited.insert(ep);

    for (int level = max_level-1; level >= 0; level--){
        // cout << "----------------------level = " << level << endl;
        searchLayer(q, u, d, top_k, indptr, level_offset, level, index, visited, vect, k);
        // cout << "visited " ;
        // for(auto it = visited.begin(); it != visited.end(); it++){
        //     cout << *it << " ";
        // }
        // cout << endl;
        // cout << "Size: " << top_k.getSize() << endl;
    }
    // top_k.print();
    string ans = "";
    for(int i = 0; i < k; i++){
        ans = to_string(top_k.top().first) + " " + ans;
        top_k.pop();
    }
    return ans;
}

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

    auto startTime = high_resolution_clock::now();
    int rank, size;

    //Starting MPI pipeline
    MPI_Init(NULL, NULL);
    
    // Extracting Rank and Processor Count
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int read[7];
    
    MPI_File fs;

    string path = string(argv[1]);

    readMPI((path + "/sizes.bin").c_str(), read, 7, 1, 0, fs);

    int n = read[0]; // size of levels 
    int ep = read[1]; // ep file value
    int max_level = read[2]; // max level value
    int s_level_offset = read[3]; // level offset size
    int s_index = read[4]; // index size 
    int s_indptr = read[5]; // indptr size
    int d = read[6]; // vect bin per line size

    int usize[1];
    readMPI((path + "/numusers.bin").c_str(), usize, 1, 1, 0, fs);

    int* indptr = new int[s_indptr];
    gatherMPI((path + "/indptr.bin").c_str(), indptr, s_indptr, size, rank, fs);

    int* index = new int[s_index];
    gatherMPI((path + "/index.bin").c_str(), index, s_index, size, rank, fs);

    int* level_offset = new int[s_level_offset];
    gatherMPI((path + "/level_offset.bin").c_str(), level_offset, s_level_offset, size, rank, fs);

    double* vect = new double[n*d];
    gatherMPI((path + "/vect.bin").c_str(), vect, n, d, size, rank, fs);

    double* q = new double[n*d];
    gatherMPI((path + "/users.bin").c_str(), q, n, d, size, rank, fs);

    int numuser = usize[0]/d;

    cout << rank << ": " <<numuser << " " <<usize[0] << endl;

    int start = rank*(numuser/size);
    int end = (rank+1)*(numuser/size);
    if(rank == size-1){
        end = numuser;
    }

    int k = stoi(argv[2]);


    string* answers = new string[end-start];
    int totalLen = 0;
    #pragma omp parallel for
    for(int i=start; i<end; i++){
        string ans = queryHNSW(q, i, d, k, ep, indptr, index, level_offset, max_level, vect) + "\n";
        answers[i-start] = ans;

        #pragma omp critical
        totalLen+=ans.size();
    }

    int *allSizes = new int[size];

    allSizes[rank] = totalLen;

    MPI_Allgather(&allSizes[rank], 1, MPI_INT, allSizes,1, MPI_INT, MPI_COMM_WORLD);

    int writeOffset = 0;
    for (int i = 0; i < rank; i++)
    {
        writeOffset += allSizes[i];
    }
    int sizeoftype = sizeof(char);

    ofstream out(argv[3], ios::binary);

    if(rank == size - 1){
        writeAll(writeOffset + allSizes[rank], out, sizeoftype);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    
    out.seekp(writeOffset*sizeoftype, ios::beg);
    
    for(int i=0; i<end-start; i++){
        out.write(answers[i].c_str(), answers[i].size());
        answers[i].clear();
    }
    
    delete(indptr);
    delete(index);
    delete(level_offset);
    delete(vect);
    delete(q);
    delete(allSizes);

    MPI_Finalize();
    if(rank==0){
        auto stopTime = high_resolution_clock::now();
        auto duration = duration_cast<microseconds>(stopTime - startTime);
        cout << "TIME IN SEC: " <<duration.count()/1e6 << endl;
    }
}