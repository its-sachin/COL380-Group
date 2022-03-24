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
        priority_queue<pid, vector<pid>, Less> temp = heap;
        double maxx = temp.top().second;
        while(!temp.empty()){
            maxx = max(maxx, temp.top().second);
            temp.pop();
        }
        return maxx;
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

int getsize(string fullFilePath){
    ifstream infile;
    infile.open(fullFilePath, std::ios::in);
    int a;
    int size = 0;
    while(infile >> a){  
        size++;
    }
    infile.close();
    return size;
}

int readOne(string fullFilePath){
    ifstream infile;
    infile.open(fullFilePath, std::ios::in);
    int a;
    infile >> a;
    infile.close();
    return a;
}

int main(int argc, char* argv[]){

    auto startTime = high_resolution_clock::now();
    int rank, size;

    //Starting MPI pipeline
    MPI_Init(NULL, NULL);
    
    // Extracting Rank and Processor Count
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    string path = argv[1];
    int n = getsize(path + "/level.txt");
    int ep = readOne(path + "/ep.txt");
    int max_level = readOne(path + "/max_level.txt");

    cout << rank << ": reading indptr" << endl;

    ifstream infile;
    int a;
    vector<int> nums;

    infile.open(path + "/indptr.txt", std::ios::in);
    while(infile >> a){  
        nums.push_back(a);
    }
    infile.close();
    int* indptr = new int[nums.size()];
    for(int i=0; i<nums.size(); i++)
        indptr[i] = nums[i];
    nums.clear();

    cout << rank << ": reading index" << endl;

    infile.open(path + "/index.txt", std::ios::in);
    while(infile >> a){  
        nums.push_back(a);
    }
    infile.close();
    int* index = new int[nums.size()];
    for(int i=0; i<nums.size(); i++)
        index[i] = nums[i];
    nums.clear();

    cout << rank << ": reading level_offset" << endl;

    infile.open(path + "/level_offset.txt", std::ios::in);
    while(infile >> a){  
        nums.push_back(a);
    }
    infile.close();
    int* level_offset = new int[nums.size()];
    for(int i=0; i<nums.size(); i++)
        level_offset[i] = nums[i];
    nums.clear();

    cout << rank << ": reading vect" << endl;


    vector<double> nums2;
    double b;
    
    infile.open(path + "/vect.txt", std::ios::in);
    
    while(infile >> b){  
        nums2.push_back(b);
    }
    infile.close();
    double* vect = new double[nums2.size()];
    for(int i=0; i<nums2.size(); i++)
        vect[i] = nums2[i];

    int d = nums2.size()/n;

    nums2.clear();

    cout << rank << ": reading user" << endl;

    infile.open(argv[3], std::ios::in);
    
    while(infile >> b){  
        nums2.push_back(b);
    }
    infile.close();
    double* q = new double[nums2.size()];
    for(int i=0; i<nums2.size(); i++)
        q[i] = nums2[i];

    nums2.clear();

    int start = rank*(n/size);
    int end = (rank+1)*(n/size);
    if(rank == size-1){
        end = n;
    }

    int k = stoi(argv[2]);

    string* answers = new string[end-start];
    int totalLen = 0;
    #pragma omp parallel for
    for(int i=start; i<end; i++){
        cout << rank << ": " << i << endl;
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

    ofstream out(argv[4], ios::binary);

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
    answers->clear();
    delete(allSizes);

    MPI_Finalize();
    if(rank==0){
        auto stopTime = high_resolution_clock::now();
        auto duration = duration_cast<microseconds>(stopTime - startTime);
        cout << "TIME IN SEC: " <<duration.count()/1e6 << endl;
    }
}