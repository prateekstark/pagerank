// Copyright (c) 2009-2016 Craig Henderson
// https://github.com/cdmh/mapreduce

#include <boost/config.hpp>
#if defined(BOOST_MSVC)
#   pragma warning(disable: 4127)

// turn off checked iterators to avoid performance hit
#   if !defined(__SGI_STL_PORT)  &&  !defined(_DEBUG)
#       define _SECURE_SCL 0
#       define _HAS_ITERATOR_DEBUGGING 0
#   endif
#endif

#include "../Thirdparty/mapreduce/include/mapreduce.hpp"
#include "utils.cpp"
#include<chrono>
using namespace std;
using namespace std::chrono;

double* I;
// vector<vector<pair<int, double> > > H_Compact;/
map<int, vector<pair<int, double> > > H_Map;
vector<int> activeIndices;
vector<int> dangling_nodes;
double alpha = 0.85;
int activeNodeSize;

template<typename MapTask>
class datasource : mapreduce::detail::noncopyable{
public:
    datasource(unsigned size_) : sequence_(0), size(size_){
    }

    bool const setup_key(typename MapTask::key_type &key){
    	int tempIndex = sequence_++;
        key = activeIndices[tempIndex];
        // return ((H_Compact[key].size() != 0) && (key < size));
        // return key < size;
        // return H_Map.find(key) != H_Map.end();
        return tempIndex < activeNodeSize;
    }

    bool const get_data(typename MapTask::key_type const &key, typename MapTask::value_type &value){
    	// value.clear();
    	for(int j=0;j<H_Map[key].size();j++){
			value.push_back(H_Map[key][j]);
		}
		// value = H_Map[key];

		// value.assign()
        return true;
    }

private:
    int sequence_;
    int size;
};

struct map_task : public mapreduce::map_task<int, vector<pair<int, double> > >{
    template<typename Runtime>
    void operator()(Runtime &runtime, key_type const &key, value_type const &value) const
    {
    	// if(value.size() != 0){
		typename Runtime::reduce_task_type::key_type const emit_key = key;
    	runtime.emit_intermediate(emit_key, value);	
    	// }
    }
};

struct reduce_task : public mapreduce::reduce_task<int, vector<pair<int, double> > >{
    template<typename Runtime, typename It>
    void operator()(Runtime &runtime, key_type const &key, It it, It ite) const
    {
		double sum = 0;
		vector<pair<int, double> > results;
		
		for(unsigned i=0;i<(*it).size();i++){
			sum += ((*it)[i].second)*(I[(*it)[i].first]);
		}
		sum *= alpha;
		if(sum != 0){
			results.push_back(make_pair(0, sum));
			runtime.emit(key, results);
		}
		return;
    }
};

typedef mapreduce::job<map_task, reduce_task, mapreduce::null_combiner, datasource<map_task> > job;

int main(int argc, char *argv[]){
    mapreduce::specification spec;
    
    if(argc < 2){
    	cerr << "Invalid Arguments!" << endl;
    	return 0;
    }

    string filename = "../test/" + (string)argv[1] + ".txt";

    if (argc > 2){
        spec.map_tasks = max(1, atoi(argv[1]));
    }

    if(argc > 3){
        spec.reduce_tasks = atoi(argv[2]);
    }

    else{
        spec.reduce_tasks = max(1U, thread::hardware_concurrency());
    }

    // Datasource maybe from a file!
    int size = getMatrixSize(filename);
    int n = size;
    cout << "Matrix Size: " << size << endl;

    // H_Compact.resize(n);
    // readCompactMatrixFile(filename, n, dangling_nodes, H_Compact);
    readMapMatrixFile(filename, n, dangling_nodes, H_Map, activeIndices);
    // printMap(H_Map);
   
    activeNodeSize = activeIndices.size();
    cout <<"\nPageRank analysis MapReduce..." << endl;
    auto start = high_resolution_clock::now(); 

    double DI[n] = {0};
    
    for(int i=0;i<dangling_nodes.size();i++){
    	DI[dangling_nodes[i]] = 1.0/n;
	}

	double GI[n] = {0};
	I = (double*)malloc((size)*sizeof(double));
	memset(I, 0, n*sizeof(double));
	I[0] = 1;
	
	bool converged = false;
	int index = 0;
	// return 0;
	while(index < 100){
		double randomDanglingSum = 0;
		
		for(int i=0;i<n;i++){
			randomDanglingSum += alpha*DI[i]*I[i] + ((1.0-alpha)/n)*I[i];
		}

		memset(GI, 0, n*sizeof(double));
		
		job::datasource_type datasource(size);
		job job1(datasource, spec);
	    mapreduce::results result;

		// job1.run<mapreduce::schedule_policy::sequential<job> >(result);
		job1.run<mapreduce::schedule_policy::cpu_parallel<job> >(result);
		
		for(auto it=job1.begin_results(); it!=job1.end_results(); ++it){
			GI[it->first] = it->second[0].second;
	    }
	    
		for(int i=0;i<n;i++){
			GI[i] += randomDanglingSum;
		}

		memcpy(I, GI, n*sizeof(double));
		
		index++;
	}

	// printMatrix(I, n);
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);
	cout << "Time taken by function: " << (double)duration.count()/1000000 << " seconds" << endl; 
	writeToFile((string)argv[1], 0, I, n);
    return 0;
}