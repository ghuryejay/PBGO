all:
	g++ -o min_cut -std=c++11 Louvain_lable.cpp -lNetworKit -fopenmp
	g++ create_links.cpp -o create_links
	g++ avg_neighrhood.cpp -o avg_neighrhood
	