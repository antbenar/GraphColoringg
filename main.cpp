#include <iostream> 
#include "lib/lib.cpp"

using namespace std;

void printGraph(int*  adjmatrix, int V){
	for(int i=0; i<V; ++i){
		for(int j=0; j<V; ++j){
			cout << adjmatrix[i*V + j] << ", ";
		}
		cout << endl;
	}
	cout << endl;
}
/*
"Usage: mpirun/mpiexec -np xx ./graphcoloring\n"
*/

int main(int argc, char *argv[]) {
	char* s_InputFile = (char*)"graphs/mymatrix.mtx";
	char* output_filename = (char*)"output.mto";
	
	
	int num_v_per_p, remainder_v_per_p; 
	int first_v, last_v, *range; //id's of the first and last vertices in a given processor
	int maxColor; 
	int i, j, k;
	int *p_graph, *p_graph_size, *offsets;//, *boundary_graph;
	//graph with edges corresponding to the vertices in this process and their size. 
	//If process p has vertices 0 and 1, then p_graph will be a 2xV matrix, so it has all the edges between 1,2 
	// and all vertices offsets is the index in graph where to start copying to p_graph
	
	int *vertex_offsets;
	double start_time, end_time, runtime, largest_runtime;
	
	//Initialize MPI
	MPI_Init(&argc, &argv); 
	MPI_Comm_rank(MPI_COMM_WORLD, &rango);
	MPI_Comm_size(MPI_COMM_WORLD, &num_processors);
	
	
	//only root reads the file and loads the full graph
	if (rango == root) {
		int V = ReadMatrix( s_InputFile );
		//printGraph(adjmatrix, V);
		
		printf("V = %d E = %d Chromaticity Upper Bound = %d\n", V, E, chromaticity_upper);
		printf("Root finished reading the graph from file.\n");
		fflush(stdout);
	}
	
	//root broadcasts E,V and chromaticity_upper to all other processes
	MPI_Bcast(&V,1,MPI_INT,root,MPI_COMM_WORLD);
	MPI_Bcast(&E,1,MPI_INT,root,MPI_COMM_WORLD);
	MPI_Bcast(&chromaticity_upper,1,MPI_INT,root,MPI_COMM_WORLD);

	
	// Distribution: create an initial adjacency matrix
	num_v_per_p = V/num_processors;
	p_graph_size = (int *) malloc(num_processors * sizeof(int));
	offsets = (int *)malloc(num_processors * sizeof(int));
	vertex_offsets = (int *)malloc(num_processors * sizeof(int));
	range = (int *)malloc(num_processors * sizeof(int));
	
	
	// More distribution
	for (i = 0; i < num_processors; i++) {  
		remainder_v_per_p = (V + i) % num_processors;
		first_v =  num_v_per_p * i + remainder_v_per_p * (remainder_v_per_p < i);     //index of the first vertex for process i
		last_v = (i + 1) * num_v_per_p + (remainder_v_per_p+1) * (remainder_v_per_p < i) - 1; //index of the last vertex for process i
		range[i] = last_v - first_v + 1;
		p_graph_size[i] = range[i] * V;
		
		offsets[0] = 0;
		vertex_offsets[0] = 0;
		
		if (i > 0) {
			offsets[i] = offsets[i-1] + p_graph_size[i-1];
			vertex_offsets[i] = vertex_offsets[i-1] + range[i-1];
		}
		
#if DEBUG_VERTEX_DISTRIBUTION    
		printf("V = %d rango = %d first = %d last = %d range = %d ideal_range = %d\n",V,i,first_v,last_v,range[i],(V+i)/num_processors);
#endif
	} // done with the processors
	
	p_graph = (int *) malloc(p_graph_size[rango] * sizeof(int));
	//boundary_graph = (int *) malloc(p_graph_size[rango] * sizeof(int));
	
	MPI_Scatterv(adjmatrix, p_graph_size, offsets, MPI_INT, p_graph, p_graph_size[rango], MPI_INT, root, MPI_COMM_WORLD);
	
#if DEBUG_VERTEX_DISTRIBUTION
	//check whether p_graph and graph match in the corresponding positions
	for (i = 0; i < range[rango]; i++) {
		for (j = 0; j < V; j++) {
			if (p_graph[i*V + j] != graph[offsets[rango] + i * V + j])
				printf("Incorrect subgraph assignment in p_graph process %d at row %d col %d\n",rango, i, j);
		}
	}
#endif
	
	if (rango == root) {
		boundary_matrix = (int *) malloc(V * V * sizeof(int));
		memset(boundary_matrix, -1, V * V * sizeof(int));
		create_boundary_matrix(vertex_offsets);
		printf("---Boundary Matrix---\n");
		//printf("x == -1 where no edges, x >= 0 where the processor to communicate to is known\n");
		printMatrix(V,V,boundary_matrix);
		
	}
	
	start_time = MPI_Wtime();
	
	//initialize colors to 0
	colors = (int *) malloc(V * sizeof(int));
	memset(colors, 0, V * sizeof(int));
	
	
	printf("p_graph of process #%d\n", rango);
	printArray(p_graph_size[rango], p_graph);
	
	
	// INITIAL COLORING 
	parallel_coloring(range, vertex_offsets, p_graph);
	
	//printf("Printing p_graph\n");
	//printMatrix(V, num_v_per_p, p_graph);
	
	end_time = MPI_Wtime();
	runtime = end_time-start_time;
	
	//find the largest runtime (most likely it will be root's runtime since root
	//does a few extra things like generating weights, gathering colors and
	//synchronizing them)
	MPI_Allreduce(&runtime,&largest_runtime,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
	
	printf("Max runtime was %f\n",largest_runtime);
	
	int num_uncolored =0;
	for (i=0;i<V;i++) { //find out how many vertices are uncolered
		if (colors[i] == 0)
			num_uncolored++;
	}
	if (num_uncolored > 0)
		printf("Not all vertices have been colored\n");
	
	if (rango == root){//root writes colors to file
		write_colors(output_filename);    
	}

	free(p_graph);
	free(range);
	free(p_graph_size);
	free(offsets);
	free(adjmatrix);
	free(colors);
	free(boundary_matrix);
	MPI_Finalize();  //Finalize
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
