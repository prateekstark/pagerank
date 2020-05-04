#include<iostream>
#include<vector>
#include<fstream>
#include<sstream>
#include<string>

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

void readCompactMatrixFile(string filename, int size, vector<int> &dangling_nodes, vector<vector<pair<int, double> > > &H_Compact){
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
	    			H_Compact[adjacencyList[i]].push_back(make_pair(presentVertex, 1.0/tempSize));

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
		H_Compact[adjacencyList[i]].push_back(make_pair(presentVertex, 1.0/tempSize));
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

void readMapMatrixFile(string filename, int size, vector<int> &dangling_nodes, map<int, vector<pair<int, double> > > &H_Map, vector<int> &activeIndices){
	ifstream inFile;
    inFile.open(filename);
    double x;

    int fromVertex, toVertex;
    if(!inFile){
        cout << "Unable to open file";
        exit(1);
    }
    vector<vector<pair<int, double> > > H_Compact;
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
	    			H_Compact[adjacencyList[i]].push_back(make_pair(presentVertex, 1.0/tempSize));

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
		H_Compact[adjacencyList[i]].push_back(make_pair(presentVertex, 1.0/tempSize));
	}
	for(int i=0;i<size;i++){
		if(H_Compact[i].size() != 0){
			H_Map.insert(pair<int, vector<pair<int, double> > >(i, H_Compact[i]));
			activeIndices.push_back(i);
		}
	}

	adjacencyList.clear();
	H_Compact.clear();

	dangling_nodes.clear();
	for(int i=0;i<size;i++){
		if(num_links[i] == 0){
			dangling_nodes.push_back(i);
		}
	}

    inFile.close();
}


int getMatrixSize(string filename){
	int matSize = 0;
	ifstream inFile;
    inFile.open(filename);
    // double x;

    int fromVertex, toVertex;
    if(!inFile){
        cout << "Unable to open file";
        exit(1);
    }

    string line;
    
    while(getline(inFile, line)){	
    	stringstream ss(line);
    	ss >> fromVertex;
    	ss >> toVertex;
    	matSize = max(matSize, max(fromVertex, toVertex));
    }
    matSize++;
    inFile.close();
    return matSize;
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


void printMatrixMap(vector<vector<pair<int, double> > > matrix){
	for(int i=0;i<matrix.size();i++){
		cout << i << ": ";
		for(int j=0;j<matrix[i].size();j++){
			cout << matrix[i][j].first << " ";
		}
		cout << endl;
	}
}