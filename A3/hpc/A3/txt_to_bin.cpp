#include<bits/stdc++.h>


int lTob(int n){
    // int words[4];
    // for(int i = 0; i < 4; i++){
    //     words[i] = n & 0xff;
    //     n >>= 8;
    // }

    // int ans = 0;
    // for(int i = 0; i < 4; i++){
    //     ans <<= 8;
    //     ans |= words[i];
    // }
    // return ans;
    return n;
}

void writeIt(int n, std::ofstream &out, int size){
    int conv = lTob(n);
    out.write((char*)&conv, size);
}


int main(int argc, char* argv[]){
    
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

    int indptr[7] = {0, 5, 8, 14, 17, 22, 25};
    int index[25] = {1, -1, -1, 4, -1, 3, -1, -1, 1, 5, -1, 0, -1, -1, 1, -1, -1, 3, 5, -1, 0, -1, -1, -1, -1};
    int level_offset[4] = {0, 3, 5, 6};

    std::string folder = "dummy_bin/";
    std::ofstream outfile;

    outfile.open(folder + "sizes.bin", std::ios::out | std::ios::binary);
    writeIt(6, outfile, sizeof(int));
    writeIt(ep, outfile, sizeof(int));
    writeIt(max_level, outfile, sizeof(int));
    writeIt(4, outfile, sizeof(int));
    writeIt(25, outfile, sizeof(int));
    writeIt(7, outfile, sizeof(int));
    writeIt(5, outfile, sizeof(int));

    outfile.close();

    outfile.open(folder + "index.bin", std::ios::out | std::ios::binary);
    for(int i = 0; i < 25; i++){
        writeIt(index[i], outfile, sizeof(int));
    }
    outfile.close();

    outfile.open(folder + "indptr.bin", std::ios::out | std::ios::binary);
    for(int i = 0; i < 7; i++){
        writeIt(indptr[i], outfile, sizeof(int));
    }
    outfile.close();

    outfile.open(folder + "level_offset.bin", std::ios::out | std::ios::binary);
    for(int i = 0; i < 4; i++){
        writeIt(level_offset[i], outfile, sizeof(int));
    }
    outfile.close();

    outfile.open(folder + "vect.bin", std::ios::out | std::ios::binary);
    outfile.write((char *) &a, sizeof a);
    outfile.close();

}