#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include "mmio.h"
#include <sstream>
#include <string>       
#include <map>
#include <fstream>
#include <algorithm>
#include <vector>

#define true 1
#define false 0

/* Modified the read_graph and write_file functions from github.com/sdasgup3/parallel-sudoku/tree/master/MpiVersion */

#define DEBUG_VERTEX_DISTRIBUTION 0
#define DEBUG_PARALLEL_COLORING 0

using namespace std;

char * OUTPUT_PATH = (char*)"ENS_CS249_FinalProject/";
int rango, num_processors, root = 0; 
int V, E; //number of vertices and edges
int chromaticity_upper = -1; //upper bound on the chromaticity of the graph, initialized to -1
//chromaticity = minimum number of colors required to color the graph
int *  adjmatrix; //the symmetric adjacency graph matrix.  
// The index of element at col,row is idx = row * V + col. 
// The element at index i is at row = i/V, col = i % V
int * colors; //the colors assigned to each vertex
int * boundary_matrix;

// Functions defined below
int compare (const void *a, const void *b);
int ReadMatrix(char * s_InputFile);
void write_colors (char * filename);
//void parallel_coloring(int* range,int * offsets,int * p_graph);
void printMatrix(int row, int col, int* matrix);
void printArray(int row, int* array);
void create_boundary_matrix (int* matrix);


void create_boundary_matrix(int* vertex_offsets) {
	int x = num_processors-1; 
	
	int i, j, remainder;
	for (i = 0; i < V*V ; i++){
		if (adjmatrix[i] == 1){
			remainder = i % V;
			for (j = 0; j < num_processors; j++){
				if (remainder >= vertex_offsets[j] && (remainder < vertex_offsets[j+1] || j ==  x )){
					boundary_matrix[i] = j;  
				}
			}
		}
	}  
}


void parallel_coloring(int* range,int * offsets,int * p_graph)
{
	int i,j,k;
	int num_v_per_p, remainder_v_per_p, first_v;
	int * vtx_colors, * neighbor_colors;
	int num_colors, min_color;
	int procesor_to_communicate;
	
	num_v_per_p = V/num_processors;
	vtx_colors = (int *) malloc(range[rango] * sizeof(int));   
	memset(vtx_colors,0,range[rango] * sizeof(int));
	
	
	// while the chromacity is less than the upper bound
	for (i = 0; i < chromaticity_upper; i++) { 
		remainder_v_per_p = (V + rango) % num_processors; 
		first_v =  num_v_per_p * rango +  remainder_v_per_p * (remainder_v_per_p < rango); //index of the first vertex for process i
		for (j=0; j < range[rango]; j++) {   //for each vertex of this process
			neighbor_colors = (int *) malloc(V * sizeof(int));
			memset(neighbor_colors, 0, V * sizeof(int));
			num_colors = 0;
			
			for (k = 0; k < V; k++) {              
				if (p_graph[j * V + k] == 1) {//if there is an edge between j vertex and neighbor k vertex
					if (colors[k] != 0) { //if neighbor is colored just add its color to the neighbor_colors
						neighbor_colors[num_colors] = colors[k];
						num_colors++;
					} // if neighbor is not colored, its colored when its turn comes, code below
					else {
						colors[k] = min_color + 1;
						min_color++;
					}
				} // end if
			}// end of for k loop  
			
			if (colors[first_v + j] == 0){ // the vertex hasn't been colored
				// color with the smallest color that is not in the neighbor colors
				
				// sort the neighbor_colors array
				qsort(neighbor_colors, num_colors, sizeof(int), compare);    
				
				if (num_colors == 0 || neighbor_colors[0] > 1){ //if none of the neighbors is colored or the smallest color of a neighbor is >1
					min_color = 1;
				}
				
				else {
					for (k = 0;k < V; k++) {
						// THIS COLORS AT LEAST 278
						//min_color = neighbor_colors[k] + 1;
						
						//In between a color in the array of neighbors colors if there is a gap between two of the (sorted) neighbors colors
						if (k<V-1 && (neighbor_colors[k+1]-neighbor_colors[k]>1)) {
							min_color = neighbor_colors[j-1] + 1;
							
							break;
						}
						else {
							min_color = neighbor_colors[num_colors-1] + 1;
						}
					}
				}
				//if vtx_colors[j-1]
				vtx_colors[j] = min_color; //color
				// printf("min_color %d\n", min_color);
				// printf("Process %d colored vertex %d with color %d\n", rango, j + first_v + 1, min_color);
				
#if DEBUG_PARALLEL_COLORING
				if (i==1) {
					printf("START DEBUG_PARALLEL_COLORING, its processor %d\n", rango);
					int m;
					if (num_colors == 0)
						printf("rango=%d j=%d color=%d\n",rango,j,min_color);
					for (m = 0; m < num_colors; m++) {
						printf("rango=%d j=%d color=%d neighbors colors %d\n",rango,j,min_color,neighbor_colors[m]);
					}
				}                  
				printf("END DEBUG_PARALLEL_COLORING, its processor %d\n", rango);
#endif
			} // end of if (colors[] == 0)

			free(neighbor_colors);
		} // end of j
		for (j=0; j < range[rango]; j++) {   //for each vertex of this process
			for (k = 0; k < V; k++) {   
				if (p_graph[j * V + k] == 1){
					if (vtx_colors[j] == vtx_colors[j+1]){
						//printf("HELLO ITS ME j %d k %d \n", j, k);
						vtx_colors[j]++; //color
					}}}}
		
		//each process sends the colors of its vertices to root
		MPI_Gatherv(vtx_colors,range[rango],MPI_INT,colors,range,offsets,MPI_INT,root,MPI_COMM_WORLD);
		//root synchronizes colors on all processes
		MPI_Bcast(colors,V,MPI_INT,root,MPI_COMM_WORLD);
#if DEBUG_PARALLEL_COLORING
		if (i==1) {
			int p;
			for (p=0;p < range[rango];p++){      
				printf("Checking copy from vtx_colors to colors:rango=%d vertex=%d j_color=%d colors=%d\n",rango,offsets[rango]+p,vtx_colors[offsets[rango]+p],colors[offsets[rango] + p]);
			}
		}
#endif
	}
	
	free(vtx_colors);
}



// be able to print for testing purposes
void printArray (int size, int* array){
	int i;
	printf("\n");
	for (i = 0; i < size; i++) {
		printf("%d ", array[i]);
	}
	printf("\n");
	printf("\n\n");
}

// be able to print for testing purposes
void printMatrix (int row, int col, int* matrix){
	int i,j;
	
	printf("\n");
	
	for (i = 0; i < row; i++) {
		for (j = 0; j < col; j++) {
			printf("%d ", matrix[i * col + j]);
		}
		printf("\n");
	}
	printf("\n\n");
}

int compare (const void *a, const void *b)
{
	int la = *(const int*) a;
	int lb = *(const int*) b;
	return (la > lb)-(la < lb);
}

//only root will call this function to read the file
//unless DEBUG_VERTEX_DISTRIBUTION is true in which case
//all processes call this function and build their own graph
int ReadMatrix(char * s_InputFile)
{
	vector<int> m_vi_Vertices;
	vector<int> m_vi_Edges;
	int cantidadAristas = 0, posEdges=0;
	
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
		return 0;
	}

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
	
	adjmatrix = (int *) malloc(row * row * sizeof(int));
	memset(adjmatrix, 0, row * row * sizeof(int));
	
	for(int i=0; i<m_vi_Vertices.size()-1; ++i){
		cantidadAristas = m_vi_Vertices[i+1] - m_vi_Vertices[i]; //cantidad de aristas que tiene el vértice actual
		
		for(int j=0; j<cantidadAristas; ++j){
			//cout << "row: " << row << ", col: " << col << " edge: " << m_vi_Edges[posEdges]<<endl;
			adjmatrix[i*row + m_vi_Edges[posEdges]] = 1;
			adjmatrix[i*m_vi_Edges[posEdges] + row] = 1;
			++posEdges;
		}
	}

	
	//some files dont have the max_degree. Set it to the max possible chromaticity
	//of a graph which is V if V is odd and V-1 if V is even and would only
	//happen if the graph were complete, i.e. every pair of distinct vertices is
	//connected by a unique edge. Source: Wikipedia page on complete graphs
	if (row % 2 == 0) 
	{
		chromaticity_upper = V - 1;
	}
	else {
		chromaticity_upper = V;
	}
	return row;
}


void write_colors(char * filename){
	FILE* file = fopen(filename,"w");
	if (!file) {
		printf("Unable to open file %s\n",filename);
		return;
	}
	int i;   
	for (i = 0; i < V; i++) {
		fprintf(file,"Vertex = %d has color = %d\n", (i+1), colors[i]);
	}
	//write largest color i.e. chromatic number
	qsort(colors, V, sizeof(int), compare);
	fprintf(file,"Largest color was %d\n", colors[V-1]);
	fclose(file);
}
