#include<bits/stdc++.h>

using namespace std;
typedef pair<int, double> pid;


class Heap{

    public:
    priority_queue<pid, vector<pid>, greater<pid> > heap;
    // priority_queue<pid> heap;    
    int size;

    Heap(int k){
        size=k;
    }

    void push(pid p){
        heap.push({p.second, p.first});
    }

    pid pop(){
        pid p = heap.top();
        heap.pop();
        return p;
    }

    bool has_space(){
        return heap.size()<size;
    }
  
    // TODO:
    int getMax(){
        
    }

    void trim(){

    }
};

int dotProduct(int* a, int* b, int n)
{
    int sum = 0;
    for(int i = 0; i < n; i++)
        sum += a[i] * b[i];
    return sum;
}

double cosine_dist(int* a, int* b, int n)
{
    int dot = dotProduct(a, b, n);
    int mod_a = dotProduct(a, a, n);
    int mod_b = dotProduct(b, b, n);
    return (double) dot / (sqrt(mod_a) * sqrt(mod_b));
}



// --------------------------------------------------

void searchLayer(int* q, int d, Heap& top_k, int* indptr, int* level_offset, int lc, 
                unordered_set<int>& visited, int** vect)
{
    Heap candidats = Heap(top_k.size);
    candidats.heap = top_k.heap;

    while(candidats.size > 0){
        int ep = candidats.pop().first;

        int start = indptr[ep] + level_offset[lc];
        int end = indptr[ep] + level_offset[lc+1];

        for(int px = start; px < end; px++){
            if(visited.find(px) == visited.end() || px == -1) continue;
            visited.insert(px);
            int dist = cosine_dist(q, vect[px], d);
            if(dist > top_k.getMax() && !top_k.has_space()) continue;
            top_k.push(pid(px,dist));
            top_k.trim();
            candidats.push(pid(px,dist));
        }
    }
} 

int* queryHNSW(int* q, int d, int k, int ep, int* indptr, int* index, int* level_offset, int max_level, int** vect)
{
    Heap top_k = Heap(k);
    top_k.push(pid(ep,cosine_dist(vect[ep], q, d)));
    unordered_set<int> visited;
    visited.insert(ep);

    for (int level = 0; level < max_level; level++)
        searchLayer(q, d, top_k, indptr, level_offset, level, visited, vect);

    int ans[k];
    for(int i = 0; i < k; i++)
        ans[i] = top_k.pop().first;
    return ans;
}

int main(int argc, char const *argv[])
{
    return 0;
}