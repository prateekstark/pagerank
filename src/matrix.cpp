#include<iostream>
#include<vector>
#include<algorithm>
#include<memory.h>
#include<string>
#include<iostream>
#include<iomanip>
#include<fstream>
#include<sstream>
#include "MapReduce.cpp"
using namespace std;

void printMatrix(double* array, int n){
	for(int i=0;i<n;i++){
		cout << array[i] << " ";
	}
	cout << endl;
}

class Matrix{
public:

	double* H;
	int size;
	double alpha = 0.85;
	MapReduce mr;

	vector<int> dangling_nodes;

	void set_size(int size){
		this->size = size;
		this->H = (double*)malloc((size*size)*sizeof(double));
	}

	void update_dangling_nodes(int* num_links){
		dangling_nodes.clear();
		for(int i=0;i<this->size;i++){
			if(num_links[i] == 0){
				dangling_nodes.push_back(i);
			}
		}
	}

	void readMatrixFile(string filename){
		
		ifstream inFile;
	    inFile.open(filename);
	    double x;

	    int fromVertex, toVertex;
	    if(!inFile){
	        cout << "Unable to open file";
	        exit(1);
	    }
	    string line;
	    int presentVertex = -1;
	    vector<int> adjacencyList;
	    int tempSize;
	    
	    int size = this->size;
		int* num_links = (int*)malloc(size*sizeof(int));
	    
	    while(getline(inFile, line)){

	    	stringstream ss(line);
	    	ss >> fromVertex;
	    	ss >> toVertex;
	    	
	    	if(presentVertex != fromVertex){
	    		if(presentVertex != -1){
		    		num_links[presentVertex] = adjacencyList.size();
		    		tempSize = adjacencyList.size();

		    		for(int i=0;i<tempSize;i++){
		    			H[presentVertex + (size*adjacencyList[i])] = 1.0/tempSize;
		    		}
		    		adjacencyList.clear();
		    	}
	    		presentVertex = fromVertex;
	    	}

	    	adjacencyList.push_back(toVertex);

	    }

		num_links[presentVertex] = adjacencyList.size();
		tempSize = adjacencyList.size();

		for(int i=0;i<tempSize;i++){
			H[presentVertex + (size*adjacencyList[i])] = 1.0/tempSize;
		}
		adjacencyList.clear();

	    // printMatrix(H, 16);

	    update_dangling_nodes(num_links);   
	    free(num_links);
	    inFile.close();

    }

	// column wise is easy to load, but the dilemma lies later.
	// If we have it in row major form, easily parallelizable
	// For not we will assume row major form

	void multiply(double factor, double* H, double* I, double* product){
		double sum;
		int n = this->size;
		for(int i=0; i<n;i++){
			sum = 0.0;
			for(int j=0;j<n;j++){
				sum += H[j+(n*i)]*I[j];
			}
			product[i] = factor*sum;
		}
	}

	void multiplyMapReduce(double factor, double* H, double* I, double* product){
		mr.multiply_matrix(factor, H, I, this->size, product);
	}

	void pagerank(double* GI){
		
		int n = this->size;
		double factor = this->alpha;
		
		double I[n] = {0};
		I[0] = 1;
		
		double DI[n] = {0};
		
		bool converged = false;
		double dangling_add;
		double identity_add;
		
		int index = 0;

		while(index < 100){
			// this->multiplyMapReduce(factor, H, I, GI);
			this->multiply(factor, H, I, GI);
			for(int i=0;i<dangling_nodes.size();i++){
				DI[dangling_nodes[i]] = 1.0/n;
			}

			dangling_add = 0.0;
			identity_add = 0.0;
			
			for(int i=0;i<n;i++){
				dangling_add += DI[i]*I[i];
				identity_add += I[i];
			}

			dangling_add *= factor;
			identity_add *= ((1.0-factor)/n);

			for(int i=0;i<n;i++){
				GI[i] += dangling_add + identity_add;
			}

			memcpy(I, GI, n*sizeof(double));

			index++;
		}

	}
};