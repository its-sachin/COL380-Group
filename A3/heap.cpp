#include<bits/stdc++.h>


using namespace std;

// parent ->i then 
// left child -> 2*i+1
// right child -> 2*i+2 

// parent(i)-> i//2 -1

typedef pair<int, double> pid;

struct Less
{
   bool operator()(const pid& a, const pid& b){
       return a.second > b.second;
   }
};

struct Greater
{
   bool operator()(const pid& a, const pid& b){
       return a.second < b.second;
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
    vector<pid> heap;

    int size;

    Heap(int k){
        size=k;
        make_heap(heap.begin(), heap.end(), Less());
    }

    void push(pid p){
        heap.push_back(p);
        push_heap(heap.begin(), heap.end(), Less());
    }

    pid pop(){
        pid p = heap.front();
        pop_heap(heap.begin(), heap.end(), Less());
        heap.pop_back();

        return p;
    }

    bool has_space(){
        return heap.size()<size;
    }
  
    double getMax(bool isInd = false){
        int ind = 0;
        for(int i=0; i<heap.size(); i++){
            if(heap[i].second > heap[ind].second){
                ind = i;
            }
        }
        if(isInd){
            return ind;
        }
        return heap[ind].second;
    }

    void trim(){
        if(heap.size()>size){
            int m = getMax(true);
            swap(heap[m], heap[heap.size()-1]);
            heap.pop_back();
            make_heap(heap.begin(), heap.end(), Less());
        }
    }
};

// int main(int argc, char const *argv[])
// {
//     int size=10;
//     int array[size] = {9,3,8,-1,33,-43,9,11,99,0};

//     Heap* h = new Heap(8);
//     for(int i=0; i<size; i++){
//         cout << "pushing " << array[i] << endl;
//         h->push({i, array[i]});
        
//         for(int j=0; j<h->heap.size(); j++){
//             cout << "(" << h->heap[j].first << " " << h->heap[j].second << "), ";
//         }
//         cout << endl;
//     }

//     cout << "popping " << endl;
//     pid a = h->pop();
//     cout << "(" << a.first << " " << a.second << ")" << endl;
//     for(int j=0; j<h->heap.size(); j++){
//         cout << "(" << h->heap[j].first << " " << h->heap[j].second << "), ";
//     }

//     cout << endl;

//     cout << "getting max" << endl;
//     cout << h->getMax() << endl;

//     cout << "trimming" << endl;
//     h->trim();

//     for(int j=0; j<h->heap.size(); j++){
//         cout << "(" << h->heap[j].first << " " << h->heap[j].second << "), ";
//     }

//     return 0;
// }
