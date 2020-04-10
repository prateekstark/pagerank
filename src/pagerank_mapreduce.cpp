#include<boost/config.hpp>

#if defined(BOOST_MSVC)
#   pragma warning(disable: 4127)

// turn off checked iterators to avoid performance hit
#   if !defined(__SGI_STL_PORT)  &&  !defined(_DEBUG)
#       define _SECURE_SCL 0
#       define _HAS_ITERATOR_DEBUGGING 0
#   endif
#endif

#include "../Thirdparty/mapreduce/include/mapreduce.hpp"
#include "matrix.cpp"
#include <iostream>
#include <vector>
namespace pagerank_calculator{
	bool const is_prime(long const number){
	    if (number > 2)
	    {
	        if (number % 2 == 0)
	            return false;

	        long const n = std::abs(number);
	        long const sqrt_number = static_cast<long>(std::sqrt(static_cast<double>(n)));

	        for (long i = 3; i < sqrt_number; i+=2)
	        {
	            if (n % i == 0)
	                return false;
	        }
	    }
	    else if (number == 0 || number == 1)
	        return false;
	    return true;
	}
	template<typename MapTask>
	class pagerank_source : mapreduce::detail::noncopyable
	{
	  public:
		pagerank_source(vector<double> H, vector<double> I, int step, int size, int sequence, int first, int last)
	   	: H_(H), I_(I), step_(step), size_(size), sequence_(sequence), first_(first), last_(last)
	   	{
	    }

	    bool const setup_key(typename MapTask::key_type &key){
	        key = sequence_++;
	        return (key * step_ < last_);
	        return true;
	    }

	    bool const get_data(typename MapTask::key_type const &key, typename MapTask::value_type &value){
	        typename MapTask::value_type val;
	        

	        val.first  = first_ + (key * step_);
	        val.second = std::min(val.first + step_ - 1, last_);
	        std::swap(val, value);
	        return true;
	    }

	  private:
	    int sequence_;
	    int step_;
	    int size_;
	    int last_;
	    int first_;
	    vector<double> H_;
	    vector<double> I_;
	};


	struct map_task : public mapreduce::map_task<pair<char a, pair<int, int> >, std::pair<vector<double>, pair<int, int> > >
	{
	    template<typename Runtime>
	    void operator()(Runtime &runtime, key_type const &/*key*/, value_type const &value) const
	    {
	        for (key_type loop=value.first; loop<=value.second; ++loop)
	            runtime.emit_intermediate(is_prime(loop), loop);
	    }
	};

	struct reduce_task : public mapreduce::reduce_task<bool, long>
	{
	    template<typename Runtime, typename It>
	    void operator()(Runtime &runtime, key_type const &key, It it, It ite) const
	    {
	        if (key)
	            std::for_each(it, ite, std::bind(&Runtime::emit, &runtime, true, std::placeholders::_1));
	    }
	};

	typedef
	mapreduce::job<pagerank_calculator::map_task,
	               pagerank_calculator::reduce_task,
	               mapreduce::null_combiner,
	               pagerank_calculator::pagerank_source<pagerank_calculator::map_task>
	> job;

} // namespace pagerank
using namespace std;

int main(int argc, char *argv[]){

	if(argc == 3){

		int n = std::stoi(argv[1]);
		
		Matrix m;
		m.set_size(n);

		string filePath = "../test/" + (string)argv[2];

		m.readMatrixFile(filePath);
		mapreduce::specification spec;
		spec.map_tasks = 1;
		int reduce_tasks = std::max(1U, std::thread::hardware_concurrency());
		spec.reduce_tasks = reduce_tasks;
		// double* H = m.H;
		// double* I = (double*)malloc(n*sizeof(double));
		vector<double> H;
		vector<double> I;
		// I = {0};

		pagerank_calculator::job::datasource_type datasource(H, I, (n*n)/reduce_tasks, n, 0, 0, n*n);


		// double* pr = (double*)malloc(n*sizeof(double));
		// m.pagerank(pr);

		// printMatrix(pr, n);
	}
	
	else{
		cerr << "Invalid Arguments!" << endl;
	}

	return 0;
}
