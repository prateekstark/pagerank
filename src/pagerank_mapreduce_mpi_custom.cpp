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

		string filePath = "../test/" + (string)argv[1];
		

		int n = getMatrixSize(filePath);
		int size = n;
		

		double* H = (double*)malloc((size*size)*sizeof(double));
		// *H = {0};
		memset(H, 0, n*n*sizeof(double));
		
		double alpha = 0.85;
		
		vector<int> dangling_nodes;
		
		double DI[n] = {0};
		
		if(mype == 0){
			cout << "Matrix Size: " << n << endl;
			elapsed_time = -MPI_Wtime();
			readMatrixFile(filePath, n, dangling_nodes, H);
			for(int i=0;i<dangling_nodes.size();i++){
				DI[dangling_nodes[i]] = 1.0/n;
			}
		}
		// return 0;
		
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

		while(index < 2){
			printMatrix(I, n);
			for(int i=0;i<size;i++){
				if(mype == 0){
					localMapSize = 0;
					for(int j=0;j<size;j++){
						if(H[j+(size)*i] != 0 && I[j]!=0){
							localMap[localMapSize++] = H[j+(size)*i]*I[j];
							// localMap[localMapSize++] = I[j];
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
					for(int tempIndex=0;tempIndex<localMapSize;tempIndex+=1){
						// sum += localMap[tempIndex]*localMap[1+tempIndex];
						sum += localMap[tempIndex];
					}
					personalAnswer[personalIndex] = factor*sum;
				}
			}

			MPI_Reduce(&personalAnswer, &GI, n, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
			printMatrix(GI, n);
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
