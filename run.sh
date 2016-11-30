name=Max_Cluster_size_s0_161128.cpp
program=${name/.cpp/}
rm -vf read
g++ -std=c++11 $name -o $program
./$program
rm -vf $program
