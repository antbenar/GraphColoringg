#include <iostream> 
#include <sstream>
#include <string>       
#include <map>
#include <fstream>
#include <algorithm>
#include <vector>
#include "lib/mmio.h"
using namespace std;



void ReadMatrix(string s_InputFile, vector<int>& m_vi_Vertices, vector<int>& m_vi_Edges)
{
	string m_s_InputFile = s_InputFile;
	
	//initialize local data
	int col=0, row=0, rowIndex=0, colIndex=0;
	int entry_counter = 0, num_of_entries = 0;

	bool b_symmetric;
	istringstream in2;
	string line = "";
	map<int,vector<int> > nodeList;
	map<int,vector<double> > valueList;
	
	//READ IN BANNER
	MM_typecode matcode;
	FILE *f;
	if ((f = fopen(m_s_InputFile.c_str(), "r")) == NULL)  {
		cout<<m_s_InputFile<<" not Found!"<<endl;
		return;
	}

	/*
	if (mm_read_banner(f, &matcode) != 0)
	{
		printf("Could not process Matrix Market banner.\n");
		return;
	}
	
	
	if( !mm_is_coordinate(matcode) ){
		printf("Sorry, %s is dense, this application only not support sparse.\n", m_s_InputFile.c_str());
		exit(-1);
	}
	*/
	if( mm_is_symmetric(matcode) || mm_is_skew(matcode) || mm_is_hermitian(matcode) ) 
		b_symmetric = true;
	else 
		b_symmetric = false;
	
	fclose(f);  //mm_read_mtx_crd_size(f, &row, &col, &num_of_entries);  //FILE sys is kind of old.
	
	ifstream in (m_s_InputFile.c_str());
	if(!in) { printf("cannot open %s",m_s_InputFile.c_str()); exit(1); }
	
	do{
		getline(in,line);
	}while(line.size()>0 && line[0]=='%');
	
	in2.str(line);
	in2>>row>>col>>num_of_entries;

	
	// DONE - FIND OUT THE SIZE OF THE MATRIX
	
	if(b_symmetric){
		do{
			getline(in,line);
			if(line.empty() || line[0]=='%') 
				continue;
			entry_counter++;
			in2.clear();
			in2.str(line);
			in2>>rowIndex>>colIndex;
			
			if(rowIndex==colIndex) continue;
			rowIndex--;  //to 0 base
			colIndex--;  //to 0 base
			if(rowIndex<colIndex) { 
				printf("Error find a entry in symmetric matrix %s upper part. row %d col %d"
					   ,m_s_InputFile.c_str(), rowIndex+1, colIndex+1); 
				exit(1); 
			} 
			nodeList[rowIndex].push_back(colIndex);
			nodeList[colIndex].push_back(rowIndex);
		}while(!in.eof() && entry_counter<num_of_entries);
	}//end of if b_symmetric
	else{
		// if the graph is non symmetric, this matrix represent a directed graph
		// We force directed graph to become un-directed graph by adding the 
		// corresponding edge. This may leads to duplicate edges problem.
		// we then later sort and unique the duplicated edges.
		do{
			getline(in,line);
			if(line.empty() || line[0]=='%') 
				continue;
			entry_counter++;
			in2.clear();
			in2.str(line);
			in2>>rowIndex>>colIndex;
			
			if(rowIndex==colIndex) continue;
			rowIndex--;  //to 0 base
			colIndex--;  //to 0 base
			nodeList[rowIndex].push_back(colIndex);
			nodeList[colIndex].push_back(rowIndex);
		}while(!in.eof() && entry_counter<num_of_entries);
		
		for(auto &piv : nodeList){
			vector<int>& vec = piv.second;
			sort(vec.begin(), vec.end());
			vec.erase( unique(vec.begin(), vec.end()), vec.end());
		}
		row=col=max(row,col);
		
	}
	
	if(entry_counter<num_of_entries) { //entry_counter should be == num_of_entries
		fprintf(stderr,"Error: GraphInputOutput::ReadMatrixMarketAdjacencyGraph()\n");
		fprintf(stderr,"       Tries to read matrix %s\n",m_s_InputFile.c_str());
		fprintf(stderr,"       entry_counter %d < expected_entries %d\n",entry_counter, num_of_entries);
		fprintf(stderr,"       This may caused by trancate of the matrix file\n");
		exit(1);
	}
	
	//now construct the graph
	m_vi_Vertices.push_back(m_vi_Edges.size()); //m_vi_Edges.size() == 0 at this point
	for(int i=0;i<row; i++) {
		m_vi_Edges.insert(m_vi_Edges.end(),nodeList[i].begin(),nodeList[i].end());
		m_vi_Vertices.push_back(m_vi_Edges.size());
	}
	
	//CalculateVertexDegrees();
}

void print(vector<int>& toPrint){
	for(int i=0; i<toPrint.size(); ++i){
		cout << toPrint[i] << ", ";
	}
	cout << endl;
}

void printGraph(vector<int>& vertices, vector<int>& edges){
	int cantidadAristas;
	int posEdges = 0;
	
	cout << "Vertices con sus respectivas aristas: " << endl << endl;
	
	for(int i=0; i<vertices.size()-1; ++i){
		cout << "Vertice " << i << ": ";
		
		cantidadAristas = vertices[i+1] - vertices[i]; //cantidad de aristas que tiene el vértice actual
		
		for(int j=0; j<cantidadAristas; ++j){
			cout << edges[posEdges] << ", ";
			++posEdges;
		}
		cout << endl;
	}
	cout << endl;
}


int main(int argc, char *argv[]) {
	vector<int> edges;
	vector<int> vertices;
	string s_InputFile = "graphs/mymatrix.mtx";
		
	ReadMatrix( s_InputFile, vertices, edges );
	printGraph(vertices,edges);
	return 0;
}
/*
///mymatrix.mtx contiene la siguiente informacion, de la cual no se obtienen los pesos ya que no son necesarios para el coloreo 
//y todos los indices se restan 1 para que empiece desde 0 y no desde 1
//la primera linea con % es comentario

%%MatrixMarket matrix coordinate real general
3 6 8
1 1 2
1 2 1
1 4 -0.9823
1 5 -0.0071
2 3 2.1972
2 4 2.7726
3 5 3.6428
3 6 3.2721
*/
