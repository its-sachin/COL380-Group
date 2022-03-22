#include<bits/stdc++.h>
#include <mpi.h>
// #include 

using namespace std;


typedef pair<int, double> pid;

struct Less
{
   bool operator()(const pid& a, const pid& b){
       return a.second > b.second;
   }
};


class Heap{

    private:
    void swap(pid& a, pid& b){
        pid temp = a;
        a = b;
        b = temp;
    }

    public:
    priority_queue<pid, vector<pid>, Less> heap;
    int size;

    Heap(int k){
        size=k;
    }

    void push(pid p){
        heap.push(p);
    }

    pid pop(){
        pid p = heap.top();
        heap.pop();
        return p;
    }

    bool has_space(){
        return heap.size()<size;
    }

    int getSize(){
        return heap.size();
    }
  
    double getMax(){
        return heap.top().second;
    }

    // WORKS IF heap.size = size + 1
    void trim(){
        if(heap.size()>size){
            heap.pop();
        }
    }

    void print(){
        priority_queue<pid, vector<pid>, Less> temp = heap;
        while(!temp.empty()){
            cout << "("<< temp.top().first << " " << temp.top().second << "), ";
            temp.pop();
        }
        cout << endl;
    }
};
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
    return dot / (sqrt(mod_a) * sqrt(mod_b));
}



// --------------------------------------------------

void searchLayer(double* q, int u, int d, Heap& top_k, int* indptr, int* level_offset, int lc, 
                int* index, unordered_set<int>& visited, double* vect)
{
    Heap candidats = Heap(top_k.size);
    candidats.heap = top_k.heap;

    while(candidats.getSize() > 0){
        int ep = candidats.pop().first;

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
            if(dist < top_k.getMax() && !top_k.has_space()) continue;
            top_k.push(pid(px,dist));
            top_k.trim();
            candidats.push(pid(px,dist));
        }
    }
} 

string queryHNSW(double* q, int u, int d, int k, int ep, int* indptr, int* index, int* level_offset, int max_level, double* vect)
{
    Heap top_k = Heap(k);
    top_k.push(pid(ep,cosine_dist(vect, ep, q, u, d)));
    // cout << "pushing " << ep << " " << top_k.getMax() << endl;
    unordered_set<int> visited;
    visited.insert(ep);

    for (int level = max_level-1; level >= 0; level--){
        // cout<< "level: " << level << endl;
        // cout << "----------------------level = " << level << endl;
        searchLayer(q, u, d, top_k, indptr, level_offset, level, index, visited, vect);
        // cout << "visited " ;
        // for(auto it = visited.begin(); it != visited.end(); it++){
        //     cout << *it << " ";
        // }
        // cout << endl;
        // cout << "Size: " << top_k.getSize() << endl;
    }
    // top_k.print();
    string ans = "";
    for(int i = k-1; i >= 0; i--){
        ans += to_string(top_k.pop().first) + " ";
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

        MPI_Bcast(temp, left, MPI_INT, size-1, MPI_COMM_WORLD);
        for(int i=0; i<left; i++){
            arr[size*(len/size) + i] = temp[i];
        }
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

    int rank, size;

    //Starting MPI pipeline
    MPI_Init(NULL, NULL);
    
    // Extracting Rank and Processor Count
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int read[7];
    
    MPI_File fs;

    readMPI("dummy_bin/sizes.bin", read, 7, 1, 0, fs);

    int n = read[0]; // level file size
    int ep = read[1]; // 
    int max_level = read[2];
    int s_level_offset = read[3];
    int s_index = read[4];
    int s_indptr = read[5];
    int d = read[6];

    int* indptr = new int[s_indptr];
    gatherMPI("dummy_bin/indptr.bin", indptr, s_indptr, size, rank, fs);


    int* index = new int[s_index];
    gatherMPI("dummy_bin/index.bin", index, s_index, size, rank, fs);

    int* level_offset = new int[s_level_offset];
    gatherMPI("dummy_bin/level_offset.bin", level_offset, s_level_offset, size, rank, fs);

    double* vect = new double[n*d];
    gatherMPI("dummy_bin/vect.bin", vect, n, d, size, rank, fs);

    double* q = new double[n*d];
    gatherMPI("dummy_bin/vect.bin", q, n, d, size, rank, fs);

    int start = rank*(n/size);
    int end = (rank+1)*(n/size);
    if(rank == size-1){
        end = n;
    }

    int k = 3;

    vector<string> answers(end-start);
    // #pragma omp parallel num_threads(4)
    for(int i=start; i<end; i++){
        //cout << rank << ": " <<"i: " << i << endl;
        string ans = queryHNSW(q, i, d, k, ep, indptr, index, level_offset, max_level, vect) + "\n";
        //cout<<ans;
        answers[i-start] = ans;
        // cout << "[" <<rank << ": top k for q[" << i << "] =";
        // for(int j=0; j<k; j++){
        //     cout << ans[j] << " ";
        // }
        // cout <<"]"<< endl;
    }

    int totalLen = 0;
    for (int i = 0; i < answers.size(); i++)
    {
        cout<<answers[0];
        totalLen+=answers.size();
    }

    int *allSizes = new int[size];

    allSizes[rank] = totalLen;

    MPI_Allgather(&allSizes[rank], 1, MPI_INT, allSizes,1, MPI_INT, MPI_COMM_WORLD);

    int writeOffset = 0;
    for (int i = 0; i < rank; i++)
    {
        writeOffset+=allSizes[i];
    }
    int sizeoftype = sizeof(char);
    

    string output =  "finaloutput";
    std::ofstream outfile(output, std::ios::out);
    

    if(rank == size -1){
        writeAll(writeOffset, outfile, sizeoftype);
    }

    
    MPI_Barrier(MPI_COMM_WORLD);
    outfile.seekp(writeOffset, std::ios::beg);

    for (int i = 0; i < answers.size(); i++)
    {
        
        outfile.write(answers[i].c_str(),answers[i].size());
        //outfile.write((char*)&nums[i], sizeoftype);
    }
    
    // for(int i=0; i<nums.size(); i++){
    //     outfile.write((char*)&nums[i], sizeoftype);
    // }


    outfile.close();
    
    MPI_Finalize();
}