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
// #include "MapReduce.cpp"

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

class MapReduce{
public:
	MyMap map_function(const double* M, const double* N, const int size){
		MyMap mymap;
		for(int i=0;i<size;i++){
			for(int j=0;j<size;j++){
				// for(int k = 0; k<1; k++){
				if(M[j+(size)*i] != 0 && N[j]!=0){
					mymap.addElement(i, M[j+(size)*i]);
					mymap.addElement(i, N[j]);
				}
				// }
			}
		}

		/*
		for(int j=0;j<size;j++){
			for(int i=0;i<size;i++){
				// if(mymap.mymap.find(i) != mymap.mymap.end()){
					// mymap.addElement(i, N[j]);	
				// }
			}
		}
		
		*/
		// printMap(mymap.mymap);
		// cout << mymap.mymap.size() << endl;
		return mymap;
	}

	void reduce_function(double factor, MyMap mymap, int size, double* answer){
		*answer = {0};
		for(map<int, vector<double> >::iterator iter = mymap.mymap.begin(); iter != mymap.mymap.end(); ++iter){
			int index = iter->first;
			double sum = 0;
			vector<double> tempVector(iter->second);
			for(int tempIndex=0;tempIndex<2*tempVector.size();tempIndex+=2){
				sum += tempVector[tempIndex]*tempVector[1+tempIndex];
			}
			
			answer[index] = factor*sum;
		}
	}

	void multiply_matrix(const double factor, const double* M, const double* N, const int size, double* answer){
		MyMap mymap = map_function(M, N, size);
		reduce_function(factor, mymap, size, answer);
	}

};

MapReduce mr;

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

    // update_dangling_nodes(num_links, dangling_nodes, size);
    inFile.close();
}



void multiplyMapReduce(double factor, double* H, double* I, double* product, int size){
	mr.multiply_matrix(factor, H, I, size, product);
}

void multiply(double factor, double* H, double* I, double* product, int size){
	double sum;
	int n = size;
	for(int i=0; i<n;i++){
		sum = 0.0;
		for(int j=0;j<n;j++){
			sum += H[j+(n*i)]*I[j];
		}
		product[i] = factor*sum;
	}
}



int main(int argc, char *argv[]){
	if(argc == 3){
  		int ierr, mype = 0, np;

  		ierr = MPI_Init(&argc, &argv);
		ierr = MPI_Comm_size(MPI_COMM_WORLD, &np);
		ierr = MPI_Comm_rank(MPI_COMM_WORLD, &mype);
		// cout << np << endl;
		// cout << mype << endl;
		
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
		
		if(mype == 0){
			readMatrixFile(filePath, n, dangling_nodes, H);	
		}
		
		double GI[n] = {0};
		double factor = alpha;
		double I[n] = {0};
		I[0] = 1;
		
		double DI[n] = {0};
		
		bool converged = false;
		double dangling_add;
		double identity_add;
		
		// double *localCalcA = (double *)malloc((n*n/np) * sizeof(double *));
		// double *localCRow = (double *)malloc((n/np) * sizeof(double *));
		double localMap[2*n] = {0};
		int localMapSize = 0;
		int index = 0;
		int partner_rank;

		int personalIndex;
		double personalAnswer[n] = {0.0};
		double sum = 0;
		int pid;

		while(index < 1){
			// cout << index << endl;
			// Matrix_Multiply(localCalcA, I, localCRow, N/np, 32, N);
			// mr.multiply_matrix(factor, localCalcA, I, n, localCRow);
			// mr.multiply_matrix(factor, H, I, n, GI);
			// multiplyMapReduce(factor, H, I, GI, n);
			// multiply(factor, H, I, GI, n);


			// MyMap mymap = mr.map_function(H, I, n);
			MyMap mymap;
			pid = 0;
			for(int i=0;i<size;i++){
				partner_rank = i%np; // error
				if(mype == 0){
					// partner_rank = (1+(pid++))%np; // error
					// cout << partner_rank << endl;
					// MPI_Bcast(&partner_rank, 1, MPI_INT, 0, MPI_COMM_WORLD);
					localMapSize = 0;
					for(int j=0;j<size;j++){
						if(H[j+(size)*i] != 0 && I[j]!=0){
							localMap[localMapSize++] = H[j+(size)*i];
							localMap[localMapSize++] = I[j];
							mymap.addElement(i, H[j+(size)*i]);
							mymap.addElement(i, I[j]);
						}
					}
					MPI_Send(&i, 1, MPI_INT, partner_rank, 0, MPI_COMM_WORLD);
					MPI_Send(&localMapSize, 1, MPI_INT, partner_rank, 0, MPI_COMM_WORLD);
					MPI_Send(&localMap, localMapSize, MPI_DOUBLE, partner_rank, 0, MPI_COMM_WORLD);
					// cout << "Sent Index: " << i << endl;
					// cout << "Sent Map Size: " << localMapSize << endl;
					// printMatrix(localMap, localMapSize);
				}
				// cout << "rank before: " << partner_rank << " " << mype << endl;
				// MPI_Barrier(MPI_COMM_WORLD);
				// cout << "rank after: " << partner_rank << " " << mype << endl;
				if(mype == partner_rank){
					// cout << "Is it entering?" << endl;
				// else{
					MPI_Recv(&personalIndex, 1, MPI_INT, 0, mype, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
					MPI_Recv(&localMapSize, 1, MPI_INT, 0, mype, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
					MPI_Recv(&localMap, localMapSize, MPI_DOUBLE, 0, mype, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
					
					// if(partner_rank == 1){
					// cout << "My ID: " << mype << endl;
					// cout << "Received Index: " << personalIndex << endl;
					// cout << "Received Map Size: " << localMapSize << endl;
					// printMatrix(localMap, localMapSize);

					// }
					sum = 0;
					for(int tempIndex=0;tempIndex<localMapSize;tempIndex+=2){
						sum += localMap[tempIndex]*localMap[1+tempIndex];
					}
					personalAnswer[personalIndex] = factor*sum;
				}
				// MPI_Barrier(MPI_COMM_WORLD);
			}
			cout << "Did it come out? " << mype << endl;

			/*for(int i=1;i<np;i++){
				for(int personalIndex=0;personalIndex<n;personalIndex++){
					if(personalAnswer[personalIndex] != 0){
						if(mype == i){
							MPI_Send(&personalIndex, 1, MPI_INT, 0, mype, MPI_COMM_WORLD);
							MPI_Send(&personalAnswer[personalIndex], 1, MPI_INT, 0, mype, MPI_COMM_WORLD);
						}
						if(mype == 0){
							MPI_Recv(&personalIndex, 1, MPI_INT, mype, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
							MPI_Recv(&localMapSize, 1, MPI_INT, mype, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
						}
					}
				}	
				
			}*/

			
			MPI_Barrier(MPI_COMM_WORLD);
			// cout << "+++++++++++++++++++++++++++++++++++++++++" << endl;
			// cout << mype << endl;
			// printMatrix(personalAnswer, n);
			// cout << "+++++++++++++++++++++++++++++++++++++++++" << endl;





			MPI_Reduce(&personalAnswer, &GI, n, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
			if(mype == 0){
				printMatrix(GI, n);
			}
			// mr.reduce_function(factor, mymap, n, GI);
			// if(mype == 0){
			// 	printMatrix(GI, n);
			// }
			if(mype == 0){
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
			}
			index++;
			MPI_Barrier(MPI_COMM_WORLD);
		}
		if(mype == 0){
			cout << "----------------------------------" << endl;
			cout << n << endl;
			printMatrix(GI, n);
			cout << "----------------------------------" << endl;
		}
		// cout << mype << endl;
		// cout << "It is out!" << endl;
		MPI_Finalize();
	}
	
	else{
		cerr << "Invalid Arguments!" << endl;
	}

	return 0;
}
