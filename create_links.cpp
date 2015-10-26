#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <map>
#include <algorithm>
#include <fstream>
#include <cctype>


#include "cmdline.h"

using namespace std;

char* getCharExpr(string s)
{
        char *a=new char[s.size()+1];
        a[s.size()]=0;
        memcpy(a,s.c_str(),s.size());
        return a;
}

string toLower(string s)
{
    string ret = "";
    for(int i = 0;i < s.length();i++)
    {
        ret += char(tolower(s[i]));
    }
    return ret;
}

int main(int argc, char* argv[])
{
	cmdline ::parser pr;
    pr.add<string>("reads",'r',"pacbio reads",true,"");
    pr.add<string>("edgefile",'e',"edgefile of kmers",true,"");
    pr.add<string>("nodefile",'n',"nodefile",true,"");
    pr.add<string>("kmercount",'c',"kmercount",true);
    pr.add<int>("kmersize",'k',"kmersize",true,15);
    pr.parse_check(argc,argv);

    int k = pr.get<int>("kmersize");
    //k -= 1;
    ifstream readfile(getCharExpr(pr.get<string>("reads")));
    if(!readfile.good()) 
    {
        cerr<<"Can not open reads file, Please give correct file, exit!!";
        return -1;
    }
    map<string,int> kmermap;
    ifstream kmerfile(getCharExpr(pr.get<string>("kmercount")));
    string fline;
    while(getline(kmerfile,fline))
    {
        istringstream iss(fline);
        string s;
        int u;
        if(!(iss >> s >> u))
        {
            //cout<<"here";
            break;
        }
        if(u >= 2)
        {
        	kmermap[s] = u;
        }
    }
    ofstream nofile(getCharExpr(pr.get<string>("nodefile")));
    ofstream eofile(getCharExpr(pr.get<string>("edgefile")));
    long readcount = 0;
    map<string,long> nodemap;
    //read fasta file
    std::string line, name, content;
    while( std::getline( readfile, line )){
        if( line.empty() || line[0] == '>' ){ // Identifier marker
            if( !name.empty() ){ // Print out what we read from the last entry
            	string prevkmer = "";
            	for(int i = 0;i < content.length() - k + 1;i++)
            	{
            		long curr_node_id;
            		string currkmer = content.substr(i,k);
            		map<string,long> :: iterator it = nodemap.find(currkmer);
            		if(it == nodemap.end())
            		{
            			nodemap[currkmer] = readcount;
            			curr_node_id = readcount;
            			readcount += 1;
            		}
            		else
            		{
            			curr_node_id = it->second;
            		}
            		if(prevkmer == "")
            		{
            			prevkmer = currkmer;
            			continue;
            		}
            		long prevkmerid = nodemap[prevkmer];
            		eofile << prevkmerid <<"\t"<<curr_node_id<<endl;

            	}
                name.clear();
            }
            if( !line.empty() ){
                name = line.substr(1);
            }
            content.clear();
        } else if( !name.empty() ){
            if( line.find(' ') != std::string::npos ){ // Invalid sequence--no spaces allowed
                name.clear();
                content.clear();
            } else {
                content += line;
            }
        }
    }
    if( !name.empty() ){ // Print out what we read from the last entry
        string prevkmer = "";
        for(int i = 0;i < content.length() - k + 1;i++)
    	{
    		long curr_node_id;
    		string currkmer = content.substr(i,k);
    		map<string,long> :: iterator it = nodemap.find(currkmer);
    		if(it == nodemap.end())
    		{
    			nodemap[currkmer] = readcount;
    			curr_node_id = readcount;
    			readcount += 1;
    		}
    		else
    		{
    			curr_node_id = it->second;
    		}
    		if(prevkmer == "")
    		{
    			prevkmer = currkmer;
    			continue;
    		}
    		long prevkmerid = nodemap[prevkmer];
    		eofile << prevkmerid <<"\t"<<curr_node_id<<endl;

    	}
    }

    map<string,long> ::iterator it;
    for(it = nodemap.begin();it != nodemap.end(); ++it)
    {
    	nofile << it->first << "\t" << it->second << endl;
    }
    return 0;
}