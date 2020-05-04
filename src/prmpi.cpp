#include "../Thirdparty/mrmpi/src/mapreduce.h"
#include "../Thirdparty/mrmpi/src/keyvalue.h"
#include "utils.cpp"

using namespace std;
using namespace MAPREDUCE_NS;

double alpha = 0.85;

struct H_MapInfo{
	vector<int> dangling_nodes;
	int size;
};

struct UpdateInfo{
	double* I;
	double* GI;	
};

void readFile(int itask, char *fname, KeyValue *kv, void *ptr){
	ifstream inFile;
    inFile.open(fname);
    if(!inFile){
        cout << "Unable to open file";
        exit(1);
    }

    H_MapInfo* info = (H_MapInfo*)ptr;
    int size = info->size;
    
    double x;
    int fromVertex, toVertex;

    vector<vector<pair<double, double> > > H_Compact;
    H_Compact.resize(size);

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
	    			H_Compact[adjacencyList[i]].push_back(make_pair((double)presentVertex, 1.0/tempSize));
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
		H_Compact[adjacencyList[i]].push_back(make_pair((double)presentVertex, 1.0/tempSize));
	}

	for(int i=0;i<size;i++){
		if(H_Compact[i].size() != 0){
			for(int j=0; j<H_Compact[i].size(); j++){
				kv->add((char*)&i, sizeof(unsigned), (char*)&H_Compact[i][j].first, sizeof(double));
				kv->add((char*)&i, sizeof(unsigned), (char*)&H_Compact[i][j].second, sizeof(double));
			}
		}
	}

	adjacencyList.clear();
	H_Compact.clear();

	info->dangling_nodes.clear();
	for(int i=0;i<size;i++){
		if(num_links[i] == 0){
			info->dangling_nodes.push_back(i);

		}
	}

    inFile.close();
}


void myreduce(char *key, int keybytes, char *multivalue, int nvalues, int *valuebytes, KeyValue *kv, void *ptr){
	UpdateInfo *updateInfo = (UpdateInfo*)ptr;

	int index = *(int *)key;

	
	double* row = (double *)multivalue;

	double sum = 0;
	
	for(int i=0;i<nvalues;i+=2){
		int column = (int)row[i];
		double value = row[i+1];
		sum += updateInfo->I[column]*value;
	}
	sum *= alpha;
	
	updateInfo->GI[index] = sum;

}

int main(int argc, char *argv[]){
    
    if(argc < 2){
    	cerr << "Invalid Arguments!" << endl;
    	return 0;
    }

    string filename = "../test/" + (string)argv[1] + ".txt";
    char* fname = &filename[0];
    
	int mype, np;


    MPI_Init(&argc,&argv);

	MPI_Comm_rank(MPI_COMM_WORLD, &mype);
	MPI_Comm_size(MPI_COMM_WORLD, &np);

    double elapsed_time;
    elapsed_time = -MPI_Wtime();
    int size = getMatrixSize(filename);
    H_MapInfo info;
    info.size = size;
    
    MapReduce *mr = new MapReduce(MPI_COMM_WORLD);
    int num_keys = mr->map(argc-1, &fname, 0, 1, 0, readFile, &info);
    mr->collate(NULL);
    int n = size;
        
    if(mype == 0){
    	cout << "Matrix Size: " << size << endl;
    	cout << "Number of keys: " << num_keys << endl;
    	cout <<"\nPageRank analysis MapReduce..." << endl;	
    }

    double DI[n] = {0};
    if(mype == np-1){
    	for(int i=0;i<info.dangling_nodes.size();i++){
    		DI[info.dangling_nodes[i]] = 1.0/n;
		}
		cout << endl;
    }

    MPI_Bcast(DI, n, MPI_DOUBLE, np-1, MPI_COMM_WORLD);

	UpdateInfo updateInfo;
	updateInfo.I = (double*)malloc((size)*sizeof(double));
	updateInfo.GI = (double*)malloc((size)*sizeof(double));
	double *GI = (double*)malloc((size)*sizeof(double));
	
	memset(updateInfo.I, 0, n*sizeof(double));
	memset(updateInfo.GI, 0, n*sizeof(double));
	memset(GI, 0, n*sizeof(double));

	updateInfo.I[0] = 1;
	
	bool converged = false;
	int index = 0;
	double localAnswer[n];

	while(index < 100){

		MapReduce *mr2 = mr->copy();
		mr2->reduce(&myreduce, &updateInfo);

		MPI_Reduce(updateInfo.GI, GI, n, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

		if(mype == 0){
			double randomDanglingSum = 0;
			
			for(int i=0;i<n;i++){
				randomDanglingSum += alpha*DI[i]*updateInfo.I[i] + ((1.0-alpha)/n)*updateInfo.I[i];
			}

			for(int i=0;i<n;i++){
				GI[i] += randomDanglingSum;
			}

			memcpy(updateInfo.I, GI, n*sizeof(double));
			memset(GI, 0, n*sizeof(double));
		}
		
		memset(updateInfo.GI, 0, n*sizeof(double));
		index++;
		delete mr2;
		MPI_Bcast(updateInfo.I, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	}
	
	if(mype == 0){
		elapsed_time += MPI_Wtime();
		// printMatrix(updateInfo.I, n);
		cout << "Time taken to compute pagerank is: " << elapsed_time << "s" << endl; 
		writeToFile((string)argv[1], 2, updateInfo.I, n);
	}
	MPI_Finalize();

    return 0;
}