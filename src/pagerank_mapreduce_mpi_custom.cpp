#include<iostream>
#include<vector>
#include<algorithm>
#include<memory.h>
#include<string>
#include<iostream>
#include<iomanip>
#include<fstream>
#include<sstream>
#include<mpi.h>
using namespace std;

void printVectorMatrix(vector<vector<double> > matrix){
	for(int i=0;i<matrix.size();i++){
		for(int j=0;j<matrix[i].size();j++){
			cout << matrix[i][j] << " ";
		}
		cout << endl;
	}
}

void printMap(map<int, vector<double> > mymap){
	for(map<int, vector<double> >::iterator iter = mymap.begin(); iter != mymap.end(); ++iter){
		int index = iter->first;
		cout << index << ": ";
		vector<double> tempVector(iter->second);
		for(int i=0;i<tempVector.size();i++){
			cout << tempVector[i] << " ";
		}
		cout << endl;
	}
}

class MyMap{
public:
	map<int, vector<double> > mymap;
	void addElement(int key, double value){
		if(mymap.find(key) == mymap.end()){
			vector<double> tempVector;
			tempVector.push_back(value);
			this->mymap.insert({key, tempVector});
		}
		else{
			this->mymap[key].push_back(value);
		}
	}
};

void printMatrix(double* array, int n){
	for(int i=0;i<n;i++){
		cout << array[i] << " ";
	}
	cout << endl;
}

void printIntMatrix(int* array, int n){
	for(int i=0;i<n;i++){
		cout << array[i] << " ";
	}
	cout << endl;
}

void readMatrixFile(string filename, int size, vector<int> &dangling_nodes, double* H){
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
	int num_links[size] = {0};

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

	dangling_nodes.clear();
	for(int i=0;i<size;i++){
		if(num_links[i] == 0){
			dangling_nodes.push_back(i);
		}
	}

    inFile.close();
}


int main(int argc, char *argv[]){
	if(argc == 3){
  		int ierr, mype = 0, np;


  		double elapsed_time;
  		ierr = MPI_Init(&argc, &argv);
		MPI_Barrier(MPI_COMM_WORLD);

		ierr = MPI_Comm_size(MPI_COMM_WORLD, &np);
		ierr = MPI_Comm_rank(MPI_COMM_WORLD, &mype);
		// cout << np << endl;
		// cout << mype << endl;
		double randomDanglingSum = 0.0;
		int tag1 = 21;
		int tag2 = 22;
		int tag3 = 23;
		/*
		if(mype == 0){
			cout << np << endl;
		}
		*/

		int n = stoi(argv[1]);	
		int size = n;

		double* H = (double*)malloc((size*size)*sizeof(double));
		*H = {0};
		string filePath = "../test/" + (string)argv[2];
		
		double alpha = 0.85;
		
		vector<int> dangling_nodes;
		
		double DI[n] = {0};
		
		if(mype == 0){
			elapsed_time = -MPI_Wtime();
			readMatrixFile(filePath, n, dangling_nodes, H);
			for(int i=0;i<dangling_nodes.size();i++){
				DI[dangling_nodes[i]] = 1.0/n;
			}
		}
		
		double GI[n] = {0};
		double factor = alpha;
		double I[n] = {0};
		I[0] = 1;
		
		
		bool converged = false;
		double dangling_add;
		double identity_add;
		
		double localMap[2*n] = {0};
		int localMapSize = 0;
		int index = 0;
		int partner_rank;

		int personalIndex;
		double personalAnswer[n] = {0.0};
		double sum = 0;
		int pid;


		while(index < 100){
			// MyMap mymap;
			pid = 0;
			for(int i=0;i<size;i++){
				if(mype == 0){
					localMapSize = 0;
					for(int j=0;j<size;j++){
						if(H[j+(size)*i] != 0 && I[j]!=0){
							localMap[localMapSize++] = H[j+(size)*i];
							localMap[localMapSize++] = I[j];
							// mymap.addElement(i, H[j+(size)*i]);
							// mymap.addElement(i, I[j]);
						}
					}
					MPI_Send(&i, 1, MPI_INT, (i%np), tag1, MPI_COMM_WORLD);
					MPI_Send(&localMapSize, 1, MPI_INT, (i%np), tag2, MPI_COMM_WORLD);
					MPI_Send(&localMap, localMapSize, MPI_DOUBLE, (i%np), tag3, MPI_COMM_WORLD);
				}
				if(mype == (i%np)){
					MPI_Recv(&personalIndex, 1, MPI_INT, 0, tag1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
					MPI_Recv(&localMapSize, 1, MPI_INT, 0, tag2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
					MPI_Recv(&localMap, localMapSize, MPI_DOUBLE, 0, tag3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);	
					sum = 0;
					for(int tempIndex=0;tempIndex<localMapSize;tempIndex+=2){
						sum += localMap[tempIndex]*localMap[1+tempIndex];
					}
					personalAnswer[personalIndex] = factor*sum;
				}
			}
			
			// MPI_Barrier(MPI_COMM_WORLD);

			MPI_Reduce(&personalAnswer, &GI, n, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
			
			if(mype == 0){
				randomDanglingSum = 0;
				for(int i=0;i<n;i++){
					randomDanglingSum += factor*DI[i]*I[i] + ((1.0-factor)/n)*I[i];
				}

				for(int i=0;i<n;i++){
					GI[i] += randomDanglingSum;
				}

				memcpy(I, GI, n*sizeof(double));
			}
			
			index++;
			MPI_Barrier(MPI_COMM_WORLD);
		}
		if(mype == 0){
			elapsed_time += MPI_Wtime();
			printMatrix(GI, n);
			cout << "Time taken to compute pagerank is: " << elapsed_time << "s" << endl; 
		}
		MPI_Finalize();
	}
	
	else{
		cerr << "Invalid Arguments!" << endl;
	}

	return 0;
}
