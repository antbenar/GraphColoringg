#include "color.h"

vector<int> getVa(vector<int> Vcolor){
    return Vcolor;
}

vector<int> getVb(vector<int> Vcolor){
    vector<int> Vb;
    sort(Vcolor.begin(), Vcolor.end());
    int pos = 0;
    for(int i = 0; i < graph.size(); ++i){
        if(i == Vcolor[pos]){
            Vb.push_back(i);
            pos++;
        }
        if(pos == Vcolor.size()){
            break;
        }
    }
    return Vb;
}

vector<int> nbor(int vertex){
    vector<int> result;
    for(int i = 0; i < graph.size(); ++i){
        if(graph[vertex][i]){
            result.push_back(i);
        }
    }
    return result;
}

vector<int> nets(int vertex, vector<int> Vb){
    vector<int> rpta;
    for(int i = 0; i < Vb.size(); ++i){
        if(graph[Vb[i]][vertex]){
            rpta.push_back(Vb[i]);
        }
    }
    return rpta;
}

vector<int> vtxs(int vertex, vector<int> Va){
    vector<int> rpta;
    for(int i = 0; i < Va.size(); ++i){
        if(graph[Va[i]][vertex]){
            rpta.push_back(Va[i]);
        }
    }
    return rpta;
}

void Algorithm_1(vector<int> Vcolor){ ///Greedy Graph Coloring
    auto start = chrono::steady_clock::now();

    queue<int> W;
    vector<int> c(graph.size());
    for(int i = 0; i < Vcolor.size(); ++i){
        W.push(Vcolor[i]);
    }
    fill(c.begin(), c.end(), -1);
    while (!W.empty())
    {
        Algorithm_2(W, c);
        #pragma omp barrier
        Algorithm_3(W, c);
    }
    auto end = chrono::steady_clock::now();
    /*for(int i = 0; i < c.size(); ++i){
        cout << c[i] << " ";
    }*/
    cout << "Elapsed time in milliseconds : "
        << chrono::duration_cast<chrono::milliseconds>(end - start).count()
        << " ms" << endl;

}

void Algorithm_1_VaVb(vector<int> Vcolor){ ///Greedy Graph Coloring
    auto start = chrono::steady_clock::now();

    queue<int> W;
    vector<int> c(graph.size());
    for(int i = 0; i < Vcolor.size(); ++i){
        W.push(Vcolor[i]);
    }
    fill(c.begin(), c.end(), -1);

    vector<int> Va = getVa(Vcolor);
    vector<int> Vb = getVb(Vcolor);

    while (!W.empty())
    {
        Algorithm_4(W, Va, Vb, c);
        #pragma omp barrier
        Algorithm_5(W, Va, Vb, c);
    }
    auto end = chrono::steady_clock::now();
    /*for(int i = 0; i < c.size(); ++i){
        cout << c[i] << " ";
    }*/
    cout << "Elapsed time in milliseconds : "
        << chrono::duration_cast<chrono::milliseconds>(end - start).count()
        << " ms" << endl;

}

void Algorithm_2(queue<int> W, vector<int> &c){ ///ColorWorkQueue
    while(!W.empty()){
        #pragma omp parallel num_threads(threads)
        {
            int w;
            #pragma omp critical
            {
                w = W.front();
                W.pop();
            }
            vector<int> F;
            vector<int> nb = nbor(w);
            for(int u = 0; u < nb.size(); ++u){
                if(c[nb[u]] != -1){
                    F.push_back(c[nb[u]]);
                }
            }
            int col = 0;
            while( find(F.begin(), F.end(), col) != F.end() ){
                col++;
            }
            c[w] = col;
        }
    }
}

void Algorithm_3(queue<int> &W, vector<int> c){ ///RemoveConflicts
    queue<int> Wnext;
    while(!W.empty()){
        #pragma omp parallel shared(Wnext) num_threads(threads)
        {
            int w;
            #pragma omp critical
            {
                w = W.front();
                W.pop();
            }
            vector<int> nb = nbor(w);
            for(int i = 0; i < nb.size(); ++i){
                int u = nb[i];
                if(c[u] == c[w] && w > u){
                    Wnext.push(w);
                }
            }
        }
    }
    W = Wnext;
}

void Algorithm_4(queue<int> &W,vector<int> Va, vector<int> Vb, vector<int> &c){ ///BGPC-ColorWorkQueue-Vertex
    while(!W.empty()){
        #pragma omp parallel num_threads(threads)
        {
            int w;
            #pragma omp critical
            {
                w = W.front();
                W.pop();
            }
            vector<int> F;
            vector<int> Nets = nets(w, Vb);
            for(int i = 0; Nets.size(); ++i){
                int v = Nets[i];
                vector<int> Vtxs = vtxs(v, Va);
                for(int j = 0; j < Vtxs.size(); ++j){
                    int u = Vtxs[j];
                    if(u == w) continue;
                    if(c[u] != -1){
                        F.push_back(c[u]);
                    }
                }
            }
            int col = 0;
            while( find(F.begin(), F.end(), col) != F.end() ){
                col++;
            }
            c[w] = col;
        }
    }
}

void Algorithm_5(queue<int> &W, vector<int> Va, vector<int> Vb, vector<int> c){ ///BGPC-RemoveConflicts-Vertex
    queue<int> Wnext;
    while(!W.empty()){
        #pragma omp parallel shared(Wnext) num_threads(threads)
        {
            int w;
            #pragma omp critical
            {
                w = W.front();
                W.pop();
            }
            vector<int> Nets = nets(w, Vb);
            for(int i = 0; Nets.size(); ++i){
                int v = Nets[i];
                vector<int> Vtxs = vtxs(v, Va);
                for(int j = 0; j < Vtxs.size(); ++j){
                    int u = Vtxs[j];
                    if(u == w) continue;
                    if(c[u] == c[w]){
                        Wnext.push(w);
                    }
                }
            }
        }
    }
    W = Wnext;
}

void Algorithm_6(){ ///BGPC-ColorWorkQueue-Net-v1

}

void Algorithm_7(){ ///BGPC-RemoveConflicts-Net

}

void Algorithm_8(){ ///BGPC-ColorWorkQueue-Net

}

void Algorithm_9(){ ///D2GC-ColorWorkQueue-Net

}

void Algorithm_10(){ ///D2GC-RemoveConflicts-Net

}

void Algorithm_11(){ ///ColorWorkQueue-B1

}

void Algorithm_12(){ ///ColorWorkQueue-B2

}

void generate_bipartite_graph(){
    int n = rand() % 10000 + 11000;
    int E = n*2;
    for(int i = 0; i < n; ++i){
        vector<bool> temp;
        for(int j = 0; j < n; ++j){
            temp.push_back(0);
        }
        graph.push_back(temp);
    }
    int t = 0;
    while(t < E){
        int a = rand() % (n/2);
        int b = rand() % (n/2) + (n/2);
        if(graph[a][b]){
            continue;
        }
        else{
            graph[a][b] = 1;
            graph[b][a] = 1;
            t++;
        }
    }
}

void print_graph(){
    for(int i = 0; i < graph.size(); ++i){
        for(int j = 0; j < graph.size(); ++j){
            cout << graph[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;
}

void fill_Vcolor(vector<int> &Vcolor){
    for(int i = 0; i < graph.size()/3; ++i){
        int temp = rand() % graph.size();
        if(find( Vcolor.begin(), Vcolor.end(), temp ) == Vcolor.end()){
            ///cout << temp << " ";
            Vcolor.push_back(temp);
        }
        else{
            i--;
        }
    }
    ///cout << endl;
}

void init(int t){
    threads = t;
}
