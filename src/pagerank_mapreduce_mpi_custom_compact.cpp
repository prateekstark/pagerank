#include<iostream>
#include<vector>
#include<algorithm>
#include<memory.h>
#include<iostream>
#include<iomanip>
#include<mpi.h>
#include "utils.cpp"
using namespace std;

int main(int argc, char *argv[]){
	if(argc == 2){
  		int ierr, mype = 0, np;

  		double elapsed_time;
  		ierr = MPI_Init(&argc, &argv);
		MPI_Barrier(MPI_COMM_WORLD);

		ierr = MPI_Comm_size(MPI_COMM_WORLD, &np);
		ierr = MPI_Comm_rank(MPI_COMM_WORLD, &mype);
		
		double randomDanglingSum = 0.0;
		int tag1 = 21;
		int tag2 = 22;
		int tag3 = 23;

		string filePath = "../test/" + (string)argv[1] + ".txt";

		int n = getMatrixSize(filePath);
		int size = n;
		
		vector<vector<pair<int, double> > > H_Compact;
		H_Compact.resize(n);
		
		double alpha = 0.85;
		
		vector<int> dangling_nodes;
		
		double DI[n] = {0};
		readCompactMatrixFile(filePath, n, dangling_nodes, H_Compact);
		
		if(mype == 0){
			cout << "Matrix Size: " << n << endl;
			elapsed_time = -MPI_Wtime();	
			for(int i=0;i<dangling_nodes.size();i++){
				DI[dangling_nodes[i]] = 1.0/n;
			}
		}
		
		double GI[n] = {0};
		double factor = alpha;
		double I[n] = {0};
		
		I[0] = 1;
		
		bool converged = false;
		
		int index = 0;

		double personalAnswer[n];
		
		while(index < 100){
			memset(personalAnswer, 0, n*sizeof(double));
			for(int i=0;i<n;i++){
				if(mype == (i%np)){
					for(int j=0;j<H_Compact[i].size();j++){
						int column = H_Compact[i][j].first;
						double value = H_Compact[i][j].second;
						personalAnswer[i] += alpha*value*I[column];
					}
				}
			}
			
			
			MPI_Reduce(&personalAnswer, &GI, n, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
			
			if(mype == 0){
				randomDanglingSum = 0;
				for(int i=0;i<n;i++){
					randomDanglingSum += alpha*DI[i]*I[i] + ((1.0-alpha)/n)*I[i];
				}
				for(int i=0;i<n;i++){
					GI[i] += randomDanglingSum;
				}
				memcpy(I, GI, n*sizeof(double));
			}
			
			index++;
			MPI_Bcast(I, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		}
		if(mype == 0){
			elapsed_time += MPI_Wtime();
			
			cout << "Time taken to compute pagerank is: " << elapsed_time << "s" << endl; 
			writeToFile((string)argv[1], 1, I, n);
		}
		MPI_Finalize();
	}
	
	else{
		cerr << "Invalid Arguments!" << endl;
	}

	return 0;
}
