#include<bits/stdc++.h>

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

int* queryHNSW(double* q, int u, int d, int k, int ep, int* indptr, int* index, int* level_offset, int max_level, double* vect)
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
    top_k.print();
    int* ans = new int[k];
    for(int i = k-1; i >= 0; i--){
        ans[i] = top_k.pop().first;
    }
    return ans;
}


int main(int argc, char const *argv[])
{
    int d = 5;
    int n = 6;
    int k = 3;
    int ep = 2;
    int max_level = 3;

    double a[6][5] = {
        {-0.029495, 0.030087, -0.013817, -0.017568, 0.045015},
        {0.004782, 0.032958, -0.003865, 0.003602 ,0.022957},
        {0.004957, 0.030688, -0.036196, -0.010938, 0.008610},
        {-0.028899, 0.061895, -0.045630, 0.000977 ,0.053625},
        {0.026941, -0.000656, -0.028331, 0.038961 ,0.003586},
        {-0.014396, 0.033023, -0.012421, -0.027277, 0.017062}
    };

    double* vect = new double[n*d];
    double* q = new double[n*d];
    for(int i = 0; i < n; i++){
        for(int j = 0; j < d; j++){
            vect[i*d+j] = a[i][j];
            q[i*d+j] = a[i][j];
        }
    }

    int indptr[7] = {0, 5, 8, 14, 17, 22, 25};
    int index[25] = {1, -1, -1, 4, -1, 3, -1, -1, 1, 5, -1, 0, -1, -1, 1, -1, -1, 3, 5, -1, 0, -1, -1, -1, -1};
    int level_offset[4] = {0, 3, 5, 6};

    for(int i=0; i<n; i++){
        cout << "finding top k for q[" << i << "] = \n";
        int* ans = queryHNSW(q, i, d, k, ep, indptr, index, level_offset, max_level, vect);
        for(int j=0; j<k; j++){
            cout << ans[j] << " ";
        }
        cout <<"\n"<< endl;

    }
    
    return 0;
}