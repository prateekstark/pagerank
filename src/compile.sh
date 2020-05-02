g++-7 pagerank_mapreduce.cpp /usr/lib/x86_64-linux-gnu/libboost_system.a /usr/lib/x86_64-linux-gnu/libboost_iostreams.a /usr/lib/x86_64-linux-gnu/libboost_filesystem.a -pthread -o pagerank_mapreduce.o
mpic++ -std=c++11 pagerank_mapreduce_mpi_custom.cpp -o pagerank_mapreduce_mpi_custom.o
mpic++ -std=c++11 pagerank_mapreduce_mpi_custom_compact.cpp -o pagerank_mapreduce_mpi_custom_compact.o