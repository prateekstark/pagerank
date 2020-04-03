#include<iostream>
#include "matrix.cpp"

using namespace std;

int main(int argc, char *argv[]){

	if(argc == 3){

		int n = stoi(argv[1]);
		
		Matrix m;
		m.set_size(n);

		string filePath = "../test/" + (string)argv[2];

		m.readMatrixFile(filePath);

		double* pr = (double*)malloc(n*sizeof(double));
		m.pagerank(pr);

		printMatrix(pr, n);
	}
	
	else{
		cerr << "Invalid Arguments!" << endl;
	}

	return 0;
}

