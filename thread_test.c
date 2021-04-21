#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <map>
#include <tuple>
#include <vector>
#include <set>
#include <time.h>
#include<thread>
#include<atomic>
#include<chrono>

using namespace std;
int main(){
    unsigned int nCores = thread::hardware_concurrency();
    int nThreads = nCores -2;
    cout << "Now we have nCores= " <<  nCores << "avaliable cores"<< endl;
    //cout << "nCores: " <<  nCores << endl;
    return 0;
}
