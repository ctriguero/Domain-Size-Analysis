name=Max_Cluster_size_s0_161212.cpp
program=${name/.cpp/}
rm -vf read
g++ -O3 -std=c++11 $name -o $program
./$program
rm -vf $program
