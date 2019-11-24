#include <iostream>
#include <omp.h>
#include <vector>
#include <queue>
#include <stdlib.h>
#include <time.h>
#include <algorithm>
#include <chrono>
#include <unistd.h>

using namespace std;

vector<vector<bool>> graph;
int threads;

vector<int> nbor(int vertex);

vector<int> nets(int vertex, vector<int> Va);

vector<int> vtxs(int vertex, vector<int> Vb);

void Algorithm_1(vector<int> Vcolor); ///Greedy Graph Coloring

void Algorithm_2(queue<int> W, vector<int> &c); ///ColorWorkQueue

void Algorithm_3(queue<int> &W, vector<int> c); ///RemoveConflicts

void Algorithm_4(queue<int> &W, vector<int> Va, vector<int> Vb, vector<int> &c); ///BGPC-ColorWorkQueue-Vertex

void Algorithm_5(queue<int> &W, vector<int> Va, vector<int> Vb, vector<int> c); ///BGPC-RemoveConflicts-Vertex

void Algorithm_6(); ///BGPC-ColorWorkQueue-Net-v1

void Algorithm_7(); ///BGPC-RemoveConflicts-Net

void Algorithm_8(); ///BGPC-ColorWorkQueue-Net

void Algorithm_9(); ///D2GC-ColorWorkQueue-Net

void Algorithm_10(); ///D2GC-RemoveConflicts-Net

void Algorithm_11(); ///ColorWorkQueue-B1

void Algorithm_12(); ///ColorWorkQueue-B2

void generate_bipartite_graph();

void print_graph();
void fill_Vcolor(vector<int> &Vcolor);
void init(int t);
