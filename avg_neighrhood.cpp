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

//data structures
map<string,string> reads;
map<string,long> nodes;
map<long,vector<long> > edges;

bool pairCompare(const std::pair<int, int>& firstElem, const std::pair<int, int>& secondElem) {
  return firstElem.first < secondElem.first;

}

string clean_read(string read)
{
    int x = read.find(" ");
    if(x == string::npos)
        return read;
    else
        return read.substr(0,x);
}

string reverse_complement(string s)
{
    map<char,char> m;
    m['a'] = 't';
    m['c'] = 'g';
    m['g'] = 'c';
    m['t'] = 'a';
    string ret = "";
    for(int i = 0;i < s.length();i++)
    {
        ret += m[s[i]];
    }
    return ret;
}

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

vector<pair <int, int> > getscore(vector<long> query, map<long,long> posmap)
{   
    map<long, long> :: iterator tmp;
        
    /*
    for(tmp = posmap.begin(); tmp != posmap.end(); ++tmp)
    {
        cout<<tmp->first<<"\t"<<tmp->second<<endl;
    }
    cout<<"--------------------------------"<<endl;
    */
    vector<pair <int, int> > sig;
    long score = 0;
    for(int i = 0;i < query.size();i++)
    {
        sig.push_back(make_pair(query[i],posmap[query[i]]));
        map<long,vector<long> > :: iterator hop_1 = edges.find(query[i]);
        if(hop_1 != edges.end())
        {
            vector<long> neighs = hop_1->second;
            for(int j = 0;j < neighs.size();j++)
            {
                long cur = neighs[i];
                if(posmap.find(cur) != posmap.end())
                {
                    sig.push_back(make_pair(cur,posmap[cur]));
                }
                else
                {
                    sig.push_back(make_pair(cur,posmap[query[i]]));
                }
                
                /*
                map<long,vector<long> > :: iterator hop_2 = edges.find(neighs[i]);
                if(hop_2 != edges.end())
                {
                    vector<long> nns = hop_2->second;
                    for(int k = 0;k < nns.size();k++)
                    {
                        //long cur = neighs[i];
                        if(posmap.find(nns[k]) != posmap.end())
                        {
                            sig.push_back(make_pair(cur,posmap[nns[k]]));
                        }
                        else
                        {
                            sig.push_back(make_pair(cur,posmap[query[i]]));
                        }
                    }
                } 
                */            
            }
            
        }
    }
    sort(sig.begin(),sig.end(),pairCompare);
    return sig;
}

pair<int,int> sort_merge(vector<pair<int, int> > first, vector<pair<int, int> > second)
{
    //cout<<"sortmerge"<<endl;
    int f = 0,s = 0;
    int match_score = 0;
    vector<int> diff;
    while(f < first.size() && s < second.size())
    {
        pair<int,int> first_pair = first[f];
        pair<int,int> second_pair = second[s];
        if(first_pair.first == second_pair.first)
        {
            match_score += 1;
            diff.push_back(first_pair.second - second_pair.second);
            f++;
            s++;
            continue;
        }
        if(first_pair.first < second_pair.first)
        {
            f++;
        }
        else
        {
            s++;
        }
    }
    int sz = diff.size();
    int median;
    if(sz == 0)
    {
        median = 0;
    }
    else
    {
        sort(diff.begin(),diff.end());
        if(sz %2 == 0)
        {
            median = (diff[sz/2-1]+diff[sz])/2;
        }     
        else
        {
            median = diff[sz/2-1];
        }
    }
    //cout<<median<<endl;
    pair<int,int> ret = make_pair(match_score,median);
    return ret;
}

pair<int,int> sort_merge_modified(vector<pair<int, int> > first, vector<pair<int, int> > second,map<long,long> firstpos, map<long,long> secondpos)
{
    //cout<<"sortmerge"<<endl;
    int f = 0,s = 0;
    int match_score = 0;
    vector<int> diff;
    while(f < first.size() && s < second.size())
    {
        pair<int,int> first_pair = first[f];
        pair<int,int> second_pair = second[s];
        if(first_pair.first == second_pair.first)
        {
            match_score += 1;
            if(firstpos.find(first_pair.first) != firstpos.end() && secondpos.find(second_pair.first)!= secondpos.end())
            {
                diff.push_back(firstpos[first_pair.first] - secondpos[second_pair.first]);
            }
            //diff.push_back(first_pair.second - second_pair.second);
            f++;
            s++;
            continue;
        }
        if(first_pair.first < second_pair.first)
        {
            f++;
        }
        else
        {
            s++;
        }
    }
    int sz = diff.size();
    int median;
    if(sz == 0)
    {
        median = 0;
    }
    else
    {
        sort(diff.begin(),diff.end());
        if(sz %2 == 0)
        {
            median = (diff[sz/2-1]+diff[sz])/2;
        }     
        else
        {
            median = diff[sz/2-1];
        }
    }
    //cout<<median<<endl;
    pair<int,int> ret = make_pair(match_score,median);
    return ret;
}


double find_jaccard(int read1_start, int read1_end, int read2_start,int read2_end, map<long,long> first, map<long,long> second)
{
    int common = 0;
    map<long,long> :: iterator it;
    for(it = first.begin(); it != first.end();++it)
    {
        long label1 = it->first;
        long pos1 = it->second;
        if(pos1 > read1_start && pos1 < read1_end)
        {
            if(second.find(label1) != second.end())
            {
                long pos2 = second[label1];
                if(pos2 > read2_start && pos2 < read2_end)
                {
                    common += 1;
                }
            }
        }
    }
    double jaccard = common*1.0/(first.size() + second.size() - common);
    return jaccard;
}


int main(int argc, char* argv[])
{
    cmdline ::parser pr;
    pr.add<string>("reads",'r',"pacbio reads",true,"");
    pr.add<string>("edgefile",'e',"edgefile of kmers",true,"");
    pr.add<string>("nodefile",'n',"nodemapping",true,"");
    pr.add<int>("kmersize",'k',"kmersize",true,15);
    pr.add<string>("output",'o',"output",true,"");
    pr.parse_check(argc,argv);

    int k = pr.get<int>("kmersize");
    //k -= 1;
    ifstream readfile(getCharExpr(pr.get<string>("reads")));
    if(!readfile.good()) 
    {
        cerr<<"Can not open reads file, Please give correct file, exit!!";
        return -1;
    }

    
    //read fasta file
    std::string line, name, content;
    long kmers = 0;
    while( std::getline( readfile, line )){
        if( line.empty() || line[0] == '>' ){ // Identifier marker
            if( !name.empty() ){ // Print out what we read from the last entry
                reads[name] = toLower(content);
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
        reads[name] = toLower(content);
    }
    
    map<string,string> :: iterator p;
    /*
    for(p = reads.begin();p != reads.end();++p)
    {
        cout<<p->first<<" : "<<p->second<<endl; 
    }
    */
    readfile.close();

    //load node file

    
    ifstream nodefile(getCharExpr(pr.get<string>("nodefile")));
    if(!nodefile.good()) 
    {
        cerr<<"Can not open edge file, Please give correct file, exit!!";
        return -1;
    }


    while(getline(nodefile,line))
    {
        istringstream iss(line);
        string s;
        long u;
        if(!(iss >> s >> u))
        {
            //cout<<"here";
            break;
        }
        nodes[s] = u;
        vector<long> adj;
        edges[u] = adj;
    }

    
    nodefile.close();

    //load edge file

    ifstream edgefile(getCharExpr(pr.get<string>("edgefile")));
    if(!edgefile.good()) 
    {
        cerr<<"Can not open edge file, Please give correct file, exit!!";
        return -1;
    }

    while(getline(edgefile,line))
    {
        istringstream iss(line);
        long u,v;
        //cout<<u<<" "<<v<<endl;
        if(!(iss >> u >> v))
            break;
        edges[u].push_back(v);
        edges[v].push_back(u);
    }
    /*
    for(map<long,vector<long> > ::iterator it = edges.begin();it != edges.end();it++)
    {
        cout<<it->first<<" : ";
        vector<long> x = it->second;
        for(int i = 0;i < x.size();i++)
        {
            cout<<x[i]<<" ";
        }
        cout<<endl;
    }*/
    edgefile.close();
    //iterate through each read, find its kmers and their labels, sort labels and for all remaining reads, do binary search from their adjs

    map<string,string> :: iterator it,it1;
    ofstream ofile(getCharExpr(pr.get<string>("output")));
    map<string,int> processed;
    long totscore = 0;
    map<string,vector< pair<int, int> > > readsig;
    map<string, vector< pair<int, int> > > revreadsig;
    map<string,map<long,long> > kmerposmaps;
    map<string,map<long,long> > revkmerposmaps;
    for(it = reads.begin(); it != reads.end(); ++it)
    {
        vector<long> kmernodes;
        map<long, long> kmerpos;
        map<long, long> revkmerpos;
        string seq = it->second;
        for(int i = 0;i < seq.length() - k + 1; i++)
        {
            string kmer = seq.substr(i,k);
            //cout<<kmer<<endl;
            map<string,long> :: iterator xt = nodes.find(kmer);
            if(xt != nodes.end())
            {
                long label = nodes[kmer];
                kmernodes.push_back(label);
                kmerpos[label] = i;
                //cout<<i<<endl;
            }
        }
        kmerposmaps[it->first] = kmerpos;
        /*
        map<long, long> :: iterator tmp;
        
        for(tmp = kmerpos.begin(); tmp != kmerpos.end(); ++tmp)
        {
            cout<<tmp->first<<"\t"<<tmp->second<<endl;
        }
        cout<<"--------------------------------"<<endl;*/
        //cout<<kmerpos.size()<<"\t"<<kmernodes.size()<<endl;
        
        string revseq = reverse_complement(seq);
        vector<long> revkmernodes;
        for(int i = 0;i < revseq.length() - k + 1; i++)
        {
            string kmer = revseq.substr(i,k);
            //cout<<kmer<<endl;
            map<string,long> :: iterator xt = nodes.find(kmer);
            if(xt != nodes.end())
            {
                long label = nodes[kmer];
                revkmernodes.push_back(label);
                revkmerpos[label] = i;
            }
        }
        revkmerposmaps[it->first] = revkmerpos;
        vector< pair<int, int> > res = getscore(kmernodes,kmerpos);
        vector< pair<int ,int> > revres = getscore(revkmernodes,revkmerpos);
        readsig[it->first] = res;
        revreadsig[it->first] = revres;
        
    }
    //cout<<"signature done"<<endl;
    map<string,vector< pair<int, int> > > :: iterator readiter;
    for(readiter = readsig.begin();readiter != readsig.end();++readiter)
    {
        string read = readiter->first;
        vector<pair<int,int> > currreadsig = readiter->second;
        map<string,vector<pair<int, int> > >:: iterator findrev = revreadsig.find(read);
        vector<pair<int,int> > currrevreadsig = findrev->second;
        map<string,vector< pair<int, int> > > :: iterator readiter1;
        for(readiter1 = readsig.begin();readiter1 != readsig.end();++readiter1)
        {
            if(readiter1->first != read)
            {
                string key;
                if(readiter1->first < read)
                {
                    key = readiter1->first + read;
                }
                else
                {
                    key = read + readiter1->first;
                }
                map<string,int> :: iterator findit = processed.find(key);
                if(findit != processed.end())
                {
                    continue;
                }
                processed[key] = 1;
                string read1 = readiter1->first;
                vector<pair<int,int> > readsig1 = readiter1->second;
                vector<pair<int,int> > revreadsig1 = revreadsig[read];

                string r1t = readiter->first;
                string r2t = readiter1->first;
                string r1_seq = reads[r1t];
                string r2_seq = reads[r2t];
                string r1 = clean_read(r1t);
                string r2 = clean_read(r2t);
                
                    
                pair<int,int> one = sort_merge(currreadsig,readsig1);
                pair<int,int> two = sort_merge(currreadsig,revreadsig1);
                pair<int,int> three = sort_merge(currrevreadsig,readsig1);
                
                /*
                pair<int,int> one = sort_merge_modified(currreadsig,readsig1,kmerposmaps[r1t],kmerposmaps[r2t]);
                pair<int,int> two = sort_merge_modified(currreadsig,revreadsig1,kmerposmaps[r1t],revkmerposmaps[r2t]);
                pair<int,int> three = sort_merge_modified(currrevreadsig,readsig1,revkmerposmaps[r1t],kmerposmaps[r2t]);
                */
                int max_score = max(one.first,max(two.first,three.first));
                int read1_start,read1_end,read2_start,read2_end;
                
                
                /*
                read1 read2 shared_labels read1_strand read1_start read1_end read2_strand read2_start read2_end read2_length
                */
                int MAX_OFFSET = 25000;
                if(max_score == one.first)
                {
                    //calculate offset and write to file
                    int offset = one.second;
                    if(offset > MAX_OFFSET)
                        continue;
                    if(offset < 0)
                    {
                        offset *= -1;
                        read1_start = 0;
                        read1_end = r2_seq.length()- offset;
                        if(read1_end > r1_seq.length())
                        {
                            read1_end = r1_seq.length();
                        }
                        read2_start = offset;
                        read2_end = r2_seq.length();
                    }
                    else
                    {
                        read2_start = 0;
                        read2_end = r1_seq.length()- offset;
                        if(read2_end > r2_seq.length())
                        {
                            read2_end = r2_seq.length();
                        }
                        read1_start = offset;
                        read1_end = r1_seq.length();

                    }
                    if(max_score != 0)
                    {
                        double jaccard = find_jaccard(read1_start,read1_end,read2_start,read2_end,kmerposmaps[r1t],kmerposmaps[r2t]);
                        ofile<<r1<<"\t"<<r2<<"\t"<<one.first<<"\t"<<jaccard<<"\t"<<0<<"\t"<<read1_start<<"\t"<<read1_end<<"\t"<<r1_seq.length()<<"\t"<<0<<"\t"<<read2_start<<"\t"<<read2_end<<"\t"<<r2_seq.length()<<endl;
                        continue;
                    }
                }
                if(max_score == two.first)
                {
                    int offset = two.second;
                    if(offset > MAX_OFFSET)
                        continue;
                    if(offset < 0)
                    {
                        offset *= -1;
                        read1_start = 0;
                        read1_end = r2_seq.length()- offset;
                        if(read1_end > r1_seq.length())
                        {
                            read1_end = r1_seq.length();
                        }
                        read2_start = offset;
                        read2_end = r2.length();
                    }
                    else
                    {
                        read2_start = 0;
                        read2_end = r1_seq.length()- offset;
                        if(read2_end > r2_seq.length())
                        {
                            read2_end = r2_seq.length();
                        }
                        read1_start = offset;
                        read1_end = r1_seq.length();
                    }
                    if(max_score != 0)
                    {
                        double jaccard = find_jaccard(read1_start,read1_end,read2_start,read2_end,kmerposmaps[r1t],revkmerposmaps[r2t]);
                        ofile<<r1<<"\t"<<r2<<"\t"<<one.first<<"\t"<<jaccard<<"\t"<<0<<"\t"<<read1_start<<"\t"<<read1_end<<"\t"<<r1_seq.length()<<"\t"<<1<<"\t"<<read2_start<<"\t"<<read2_end<<"\t"<<r2_seq.length()<<endl;
                        continue;
                    }
                }
                if(max_score == three.first)
                {
                    int offset = three.second;
                    if(offset > MAX_OFFSET)
                        continue;
                    if(offset < 0)
                    {
                        offset *= -1;
                        read1_start = 0;

                        read1_end = r2_seq.length()- offset;
                        read2_start = offset;
                        read2_end = r2_seq.length();
                        if(read1_end > r1_seq.length())
                        {
                            read1_end = r1_seq.length();
                        }
                    }
                    else
                    {
                        read2_start = 0;
                        read2_end = r1_seq.length()- offset;
                        read1_start = offset;
                        read1_end = r1_seq.length();
                        if(read2_end > r2_seq.length())
                        {
                            read2_end = r2_seq.length();
                        }
                    }
                    if(max_score != 0)
                    {
                        double jaccard = find_jaccard(read1_start,read1_end,read2_start,read2_end,revkmerposmaps[r1t],kmerposmaps[r2t]);
                        ofile<<r1<<"\t"<<r2<<"\t"<<one.first<<"\t"<<jaccard<<"\t"<<1<<"\t"<<read1_start<<"\t"<<read1_end<<"\t"<<r1_seq.length()<<"\t"<<1<<"\t"<<read2_start<<"\t"<<read2_end<<"\t"<<r2_seq.length()<<endl;
                        continue;
                    }
                }
            }
        }
    }
    return 0;
}
