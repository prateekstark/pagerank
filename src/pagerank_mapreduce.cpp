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
#include <bits/stdc++.h> 

using namespace std;

double* H;
double* I;
vector<int> dangling_nodes;
int size;

template<typename MapTask>
class datasource : mapreduce::detail::noncopyable{
public:
    datasource() : sequence_(0){
    }

    bool const setup_key(typename MapTask::key_type &key){
        key = sequence_++;
        return key < size;
    }

    bool const get_data(typename MapTask::key_type const &key, typename MapTask::value_type &value){
		for(int j = 0; j<size; j++){
			if(H[j+(size)*key] != 0 && I[j] != 0){
				value.push_back(H[j+(size)*key]*I[j]);
			}
		}
        return true;
    }

private:
    int sequence_;
};

struct map_task : public mapreduce::map_task<int, vector<double> >{
    template<typename Runtime>
    void operator()(Runtime &runtime, key_type const &key, value_type const &value) const
    {
    	typename Runtime::reduce_task_type::key_type const emit_key = key;
    	// runtime.emit_intermediate(key, value);
        runtime.emit_intermediate(emit_key, value);
    }
};

struct reduce_task : public mapreduce::reduce_task<int, vector<double> >{
    template<typename Runtime, typename It>
    void operator()(Runtime &runtime, key_type const &key, It it, It ite) const
    {
    	if(it == ite){
    		return;
    	}
    	else{
    		double sum = 0.0;
    		for(It it1=++it; it1!=ite; ++it1){
            	sum += *it1;
	        }
    		// runtime.emit(key, accumulate(*it, *ite, 0));
    		runtime.emit(key, sum);
    		return;
    	}
    }
};

typedef mapreduce::job<map_task, reduce_task, mapreduce::null_combiner, datasource<map_task> > job;

int main(int argc, char *argv[]){
    mapreduce::specification spec;
    if (argc > 1){
        spec.map_tasks = max(1, atoi(argv[1]));
    }
    cout << spec.map_tasks << endl;

    if(argc > 2){
        spec.reduce_tasks = atoi(argv[2]);
    }

    else{
        spec.reduce_tasks = max(1U, thread::hardware_concurrency());
    }

    // Datasource maybe from a file!
    // int size;
    cout << "Enter size of array!" << endl;
    cin >> size;

    H = (double*)malloc((size*size)*sizeof(double));
    *H = {0};

    string filename;
    cout << "Enter File Name: " << endl;
    cin >> filename;
    filename = "../test/" + filename;

    readMatrixFile(filename, size, dangling_nodes, H);
    printMatrix(H, size*size);

    
    cout <<"\nPageRank analysis MapReduce..." << endl;
    
    int n = size;
    double alpha = 0.85;
    double DI[n] = {0};
    
    for(int i=0;i<dangling_nodes.size();i++){
		DI[dangling_nodes[i]] = 1.0/n;
	}

	double GI[n] = {0};
	I = (double*)malloc((size)*sizeof(double));
	I[n] = {0};
	I[0] = 1;
	
	bool converged = false;
	int index = 0;

	job::datasource_type datasource;
	job job1(datasource, spec);
    mapreduce::results result;

	while(index < 1){
		job1.run<mapreduce::schedule_policy::sequential<job> >(result);
		// job1.run<mapreduce::schedule_policy::cpu_parallel<job> >(result);
		cout <<"\nMapReduce finished in " << result.job_runtime.count() << "s with " << distance(job1.begin_results(), job1.end_results()) << " results\n\n";

		for(auto it=job1.begin_results(); it!=job1.end_results(); ++it){
			GI[it->first] = it->second;
	    }

		int randomDanglingSum = 0;
		for(int i=0;i<n;i++){
			randomDanglingSum += alpha*DI[i]*I[i] + ((1.0-alpha)/n)*I[i];
		}

		for(int i=0;i<n;i++){
			GI[i] += randomDanglingSum;
		}
		memcpy(I, GI, n*sizeof(double));
		memset(GI, 0, n*sizeof(double));
		index++;
	}

    return 0;
}