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

void gatherMPI(const char *filename, double* arr, int len, int size, int rank, MPI_File &fs){

    int start = rank*(len/size);
    int end = (rank+1)*(len/size);
    if(rank == size-1){
        end = len;
    }

    double* temp = new double[len];

    ifstream infile;
    infile.open(filename, ios::binary);

    infile.seekg(start*sizeof(double), ios::beg);
    //cout<<" At rank "<<rank;
    for(int i=0; i<end-start; i++){
        double curr;
        infile.read((char*)&curr, sizeof(double));
        //cout<<" "<<curr<<" "<<infile.tellg();
        temp[i] = curr;
    }
    
    MPI_Allgather(temp, len/size, MPI_DOUBLE, arr, len/size,  MPI_DOUBLE, MPI_COMM_WORLD);

    // cout<<" At rank Before "<<rank<<" ";

    // for (int i = 0; i < end-start; i++)
    // {
    //     cout<<temp[i]<<" ";
    // }
    // cout<<endl;

    
    int left = len - (size-1)*(len/size) - len/size;
    // cout << "left: " << left<<endl;/
    if(left > 0){
        for(int i=0; i<left; i++){
            temp[i] = temp[len/size + i];
        }

        MPI_Bcast(temp, left, MPI_DOUBLE, size-1, MPI_COMM_WORLD);
        for(int i=0; i<left; i++){
            arr[size*(len/size) + i] = temp[i];
        }
    }
    // cout<<" At rank After "<<rank<<" ";

    // for (int i = 0; i < len; i++)
    // {
    //     cout<<arr[i]<<" ";
    // }
    // cout<<endl;

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
    int numuser = usize[0]/d;

    // cout << rank << ": "<<"n: "<<n<<endl;
    // cout << rank << ": "<< "ep: "<<ep<<endl;
    // cout << rank << ": "<<"max level: "<<max_level<<endl;
    // cout << rank << ": "<<"s_level_offset: "<<s_level_offset<<endl;
    // cout << rank << ": "<<"s_index: "<<s_index<<endl;
    // cout << rank << ": "<<"s_indptr: "<<s_indptr<<endl;
    // cout << rank << ": "<<"d: "<<d<<endl;

    cout << rank << ": reading indptr"<<endl;
    int* indptr = new int[s_indptr];
    gatherMPI((path+"/indptr.bin").c_str(), indptr, s_indptr, size, rank, fs);

    cout << rank << ": reading index"<<endl;
    int* index = new int[s_index];
    gatherMPI((path+"/index.bin").c_str(), index, s_index, size, rank, fs);


    cout << rank << ": reading level_offset"<<endl;
    int* level_offset = new int[s_level_offset];
    gatherMPI((path+"/level_offset.bin").c_str(), level_offset, s_level_offset, size, rank, fs);

    cout << rank << ": vect"<<endl;
    double* vect = new double[n*d];
    gatherMPI((path+"/vect.bin").c_str(), vect, n*d, size, rank, fs);

   
   
    
    cout << rank << ": reading user"<<endl;
    int start = (rank*numuser)/size;
    int end = ((rank+1)*numuser)/size;
    if(rank == size-1){
        end = numuser;
    }

    double* q = new double[(end-start)*d];
    ifstream users;
    users.open(path + "users.bin", ios::binary);
    users.seekg(start*d*sizeof(double), ios::beg);
    for(int i=0; i<end-start;i ++){
        for(int j=0; j<d; j++)
            users.read((char*)&q[i*d+j], sizeof(double));
    }
    users.close();

    int k = stoi(argv[2]);


    string* answers = new string[end-start];
    int totalLen = 0;
    #pragma omp parallel for
    for(int i=0; i<end-start; i++){
        string ans = queryHNSW(q, i, d, k, ep, indptr, index, level_offset, max_level, vect) + "\n";
        answers[i] = ans;

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

    MPI_File fh;
    MPI_File_open(MPI_COMM_SELF, string(argv[3]).c_str(),MPI_MODE_CREATE | MPI_MODE_WRONLY,MPI_INFO_NULL,&fh);
    MPI_File_set_view(fh, writeOffset*sizeof(char),  MPI_CHAR, MPI_CHAR, "native", MPI_INFO_NULL);
    MPI_Status status;

    for(int i=0; i<end-start; i++){
        MPI_File_write(fh,answers[i].c_str(),answers[i].size()*sizeof(char), MPI_CHAR,&status);
        answers[i].clear();
    }

    MPI_File_close(&fh);
    
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

// 