#include <iostream>
#include <vector>
#include <ctime>
#include <fstream>
#include <string>
#include <numeric>
#include <map>
#include <algorithm>
#include <stdlib.h>
#include <cstring>
#include <set>

#include "cmdline.h"


#include <NetworKit/graph/Graph.h>
#include <NetworKit/io/GMLGraphReader.h>
#include <NetworKit/community/PLM.h>
#include <NetworKit/structures/Partition.h>	

using namespace std;
using namespace NetworKit;

char* getCharExpr(string s)
{
        char *a=new char[s.size()+1];
        a[s.size()]=0;
        memcpy(a,s.c_str(),s.size());
        return a;
}


int main(int argc, char* argv[])
{
	cmdline ::parser pr;
	pr.add<string>("dbg",'g',"de bruijn graph of reads",true,"");
	pr.add<string>("read_to_node",'r',"read to node mapping",true,"");
	pr.add<string>("label_file",'o',"outputfle",true,"");	

	pr.parse_check(argc,argv);	
	GMLGraphReader reader;
	Graph g = reader.read(getCharExpr(pr.get<string>("dbg")));

	map<string,unsigned long> read_to_node;
	map<unsigned long,vector<string> > node_to_reads;
	ifstream rfile;
	rfile.open(getCharExpr(pr.get<string>("read_to_node")));
	string line;
	/*
	while(getline(rfile,line))
	{
		istringstream iss(line);

		string a;
		unsigned long b;
		if(!(iss >> a >> b))
			break;
		read_to_node[a] = b;
		
		map<unsigned long, vector<string> > :: iterator it = node_to_reads.find(b);
		if(it == node_to_reads.end())
		{
			vector<string> ls;
			ls.push_back(a);
			node_to_reads[b] = ls;
		}
		else
		{
			node_to_reads[b].push_back(a);
		}
	}*/

	PLM c(g);
	c.run();
	cout<<c.toString()<<endl;
	Partition res = c.getPartition();
	set< set<unsigned long> > p = res.getSubsets();

	set< set<unsigned long> > :: iterator it;
	int cluster_no = 1;
	ofstream ofile;
	ofile.open(getCharExpr(pr.get<string>("label_file")));
	for(it = p.begin(); it != p.end(); ++it)
	{
		set<unsigned long> nums = *it;
		set<unsigned long> :: iterator it1;
		for(it1 = nums.begin(); it1 != nums.end(); ++it1)
		{
			//string kmer = kmers[i];
			ofile <<*it1<<"\t"<<cluster_no<<endl;
		}
		cluster_no += 1;
	}
	return 0;
}
