#include<iostream>
#include<vector>
#include<map>
using namespace std;

void printVectorMatrix(vector<vector<double> > matrix){
	for(int i=0;i<matrix.size();i++){
		for(int j=0;j<matrix[i].size();j++){
			cout << matrix[i][j] << " ";
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
				mymap.addElement(i, M[j+(size)*i]);
			}
		}

		for(int j=0;j<size;j++){
			for(int i=0;i<size;i++){
				mymap.addElement(i, N[j]);
			}
		}
		return mymap;
	}

	void reduce_function(double factor, MyMap mymap, int size, double* answer){
		for(map<int, vector<double> >::iterator iter = mymap.mymap.begin(); iter != mymap.mymap.end(); ++iter){
			int index = iter->first;
			double sum = 0;
			vector<double> tempVector(iter->second);
			for(int tempIndex=0;tempIndex<size;tempIndex++){
				sum += tempVector[tempIndex]*tempVector[size+tempIndex];
			}
			answer[index] = factor*sum;
		}
	}

	void multiply_matrix(const double factor, const double* M, const double* N, const int size, double* answer){
		MyMap mymap = map_function(M, N, size);
		reduce_function(factor, mymap, size, answer);
	}

};


// int main(){
// 	double mat1[9] = {1, 2, 3, 4, 5, 6, 7, 8, 9};
// 	double mat2[3] = {1, 0, 2};
// 	MapReduce mr;
// 	double answer[3];
// 	mr.multipy_matrix(mat1, mat2, 3, answer);
// 	for(int i=0;i<3;i++){
// 		cout << answer[i] << endl;
// 	}

// }
