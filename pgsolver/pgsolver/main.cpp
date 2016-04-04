//
//  CPPSolver.cpp
//  C++ Parity Games Solver
//

#include <array>
#include <chrono>
#include <fstream>
#include <future>
#include <iostream>
#include <sstream>
#include <map>
#include <vector>
#include <list>

class Node {
private:
    std::vector<int> adj;
    std::vector<int> inj;
    int priority;
    int player;
    int id;
public:
    std::mutex mutex;
    Node() {
        priority = -1;
        player = -1;
        id = -1;
    }
    
    Node(const Node& n)
    {
        priority = n.get_priority();
        player = n.get_player();
        id = n.get_id();

    }

    Node& operator=(Node& n)
    {
        Node m(n);
        return m;
    }
    void set_id(int j){id =j;}

    void set_priority(int pr) {
        priority = pr;
    }
    void set_player(int pl) {
        player = pl;
    }
    int get_id() const {return id;}
    int get_priority() const {
        return priority;
    }
    int get_player() const {
        return player;
    }
    std::vector<int> const &get_adj() {
        return adj;
    }
    std::vector<int> const &get_inj() {
        return inj;
    }
    void add_adj(int other) {
        mutex.lock();
        adj.push_back(other);
        mutex.unlock();
    }
    void add_inj(int other) {
        mutex.lock();
        inj.push_back(other);
        mutex.unlock();
    }
    
};

class Graph {
private:
    std::vector<Node> nodes;
    std::map<int, std::vector<int>> priorityMap;
public:
    Graph(int numNodes) {
        nodes = std::vector<Node>(numNodes);
    }
    Node &get(long n) {
        return nodes[n];
    }
    void addNode(int node, int priority, int player) {
        priorityMap[priority].push_back(node);
        nodes[node].set_priority(priority);
        nodes[node].set_player(player);
        nodes[node].set_id(node);
    }
    void addEdge(int origin, int destination) {
        nodes[origin].add_adj(destination);
        nodes[destination].add_inj(origin);
    }
    std::vector<Node> &getNodes(){
        return nodes;
    }
    long size() {
        return nodes.size();
    }
    std::map<int, std::vector<int>> &get_priority_map() {
        return priorityMap;
    }
};



void
display_graph(Graph* G) {
    std::cout << "Displaying Graph: " << std::endl;
    std::cout << "Size: " << G->size() << std::endl;
    int counter = 0;
    for (Node &n : G->getNodes()){
        std::cout << "Node:     " << counter <<std::endl;
        std::cout << " player:   " << n.get_player() << std::endl;
        std::cout << " priority: " << n.get_priority() << std::endl;
        counter++;
        for(auto & a : n.get_adj())
            std::cout << "  adj:      " << a << std::endl;
        for(auto & i : n.get_inj())
            std::cout << "  inj:      " << i << std::endl;
        
    }
}


std::vector<std::string>
&split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}


std::vector<std::string>
split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, elems);
    return elems;
}

void add_node_string(Graph *G, std::string line, std::mutex *gLock) {
    std::vector<std::string> x, edges;
    x= split( line , ' ');
    int node = std::atoi(x[0].c_str());
    gLock->lock();
    G->addNode(node, std::atoi(x[1].c_str()), std::atoi(x[2].c_str()));
    gLock->unlock();
    edges = split(x[3], ',');
    for (const auto& x : edges) {
        G->addEdge(node, atoi(x.c_str()));
    }
    
}

Graph
init_graph_from_file(std::string argf) {
    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();
    std::ifstream infile;
    infile.open(argf,std::ios::in);
    if (!infile.is_open()) {
        std::cerr << "Can't open " << &infile << "!\n";
    }
    std::string line;
    std::string first;
    std::getline(infile, first);
    int numNodes = 0;
    if (first.compare("parity") > -1) {
        std::vector<std::string> y;
        y = split(first, ' ');
        numNodes = atoi(y[1].substr(0, y[1].size()-1).c_str());
    } else {
        throw "Invalid file Passed as argument.";
    }
    Graph G(numNodes + 1);
    //Mutex pour eviter de se marcher dessus pendant le parsing
    std::mutex gLock;
    std::vector<std::future<void>> results;
    while (std::getline(infile, line)) {
        results.push_back(std::async(add_node_string, &G, std::string(line), &gLock));
    }
    infile.close();
    for (int i = 0; i < results.size(); i++) {
        results[i].wait();
    }
    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;
    std::cout << "Parsed in ........... " << elapsed_seconds.count() << "s" << std::endl;
    return G;
}

class Measure
{
private: 
    static bool initialized;
    std::map<int, int> mes;
    bool tau;
    static int size;
    static Measure* max_mes;

public:
    Measure():tau(false)
    {   
        for(int i = 0; i < size ; i ++)
        {
            mes[i] = 0;
        }
        
    }
    Measure(bool t)
    {
        tau = t;
    }
    Measure(const Measure& m)
    {
        mes = m.get_mes();
        tau = m.get_tau();
    }

    void set_tau(bool t){tau = t;}
    bool get_tau()const{return tau;}
    void addAt(int i, int j)
    {
        mes.at(i) += j;
    }

    bool eq(Measure& m)
    {
        if(m.get_tau() != this->tau)
        {
            return false;
        }
        if(m.get_tau() == this->tau && m.get_tau())
        {
            return true;
        }
        for(int i = 0; i < mes.size(); i++ )
        {
            if(mes[i] != m.get_mes()[i])
                return false;

        }
        
        return true;
    }
    static Measure* getMax()
    {
        if(initialized)
            return max_mes;
        throw std::exception();
    }

    static void init(std::map<int, std::vector<int> >& prio_map)
    {
        auto it = prio_map.end();
        it--;
        size = it->first + 1;
        
        std::list<int> comp;
        for(int i = 0 ; i < size; i++)
        {
            if(i%2==0)
            {
                comp.push_back(0);
                continue;   
            }
            if(prio_map.find(i)!= prio_map.end())
                comp.push_back(prio_map.at(i).size());
            else
                comp.push_back(0);
        }
        max_mes = new Measure(comp);
        initialized = true;
    }
   
    Measure(std::list<int> m)
    {
        tau = false;
        int i = 0;
        for(auto it = m.begin(); it != m.end(); it++)
        {
            mes.insert(std::pair<int,int>(i, *it));
            i++;
        }
    }
    bool less(Measure& m1, int prio)
    {
        if(this->tau == true && m1.get_tau() == false)
            return false;
        if(m1.get_tau() == true)
            return true;

        if(m1.getSize() != getSize())
        {
            std::cout << "Should not happen, measure of different sizes"<<std::endl;
            return false;
        }
            
        for(int i = 0; i <= prio; i++ )
        {
            if(getmAt(i) > m1.getmAt(i))
                return false;
        }
        return true;
    }
    int getmAt(int i)
    {
        return mes.at(i);
    }

    std::map<int,int> get_mes()const{return mes;}

    static int getSize(){return Measure::size;}

    std::string toString() const
    {
        if(tau)
            return " T ";
        std::string result = "(";
        for(auto n = mes.begin(); n != mes.end(); n++)
        {
            result = result + std::to_string(n->second) + ",";
        }
        result = result + ") ";
        return result;
    }

};

void display(std::map<int,Measure>& sig)
{
    for(auto it = sig.begin(); it != sig.end(); it++)
    {
        std::cout<<it->first<<" | "<<it->second.toString() <<std::endl;
    }
    std::cout<<std::endl;
}



Measure* prog(Node& v, Node& w, std::map<int,Measure>& sig)
{
  
    if(v.get_priority()%2 == 0)
    {

        return new Measure(sig[w.get_id()]);
    }
    else
    {
        if(sig[w.get_id()].less(*Measure::getMax(), v.get_priority()))
        {
            bool added = false;
            int j = v.get_priority();
            Measure* m = new Measure(sig[w.get_id()]);
            while(!added && j >= 0)
            {
                if(m->getmAt(j) < Measure::getMax()->getmAt(j))
                {
                    m->addAt(j,1);
                    added = true;
                }
                j--;

            }
            if(!added)
                return new Measure(true);
            return new Measure(*m);
        }
            
        
        return new Measure(true);
       
    }
    

}

void lift(Node& v, std::map<int,Measure>& sig, Graph& g)
{
    //display(sig);
    if(sig.at(v.get_id()).get_tau() )
        return;
    if(v.get_player() == 0)
    {
        Measure* min = Measure::getMax();
        bool assigned = false;
        for(auto n : v.get_adj())
        {
            Node w = g.get(n);
            Measure* p = prog(v,w, sig);
            if(p->get_tau())
                sig.at(w.get_id()).set_tau(true);

            if(p->less(*min, Measure::getSize()-1))
            {
                min = p;
                assigned = true;
            } 
        }
        if(!assigned)
            sig.at(v.get_id()).set_tau(true);
        else
        {
            sig.at(v.get_id()) = *min;
        }
            
    }
    else
    {
        Measure* max = new Measure(std::list<int>(Measure::getSize(), 0));
        bool assigned = false;
        for(auto n : v.get_adj())
        {
            Node w = g.get(n);
            if(sig.at(w.get_id()).eq(*(Measure::getMax())))
            {
                sig.at(v.get_id()).set_tau(true);
                sig.at(w.get_id()).set_tau(true);
                assigned = false;
                max = &sig.at(w.get_id());
                break;
            }
            
            Measure* p = prog(v,w, sig);
            if(p->get_tau()){
                sig.at(w.get_id()).set_tau(true);
                sig.at(v.get_id()).set_tau(true);
                assigned = false;
                break;
            }

            if(max->less(*p, Measure::getSize()-1))
            {
                max = p;

                assigned = true;
            }

        }
        if(assigned)
        {

            sig.at(v.get_id()) = *max;
        }
    }
}

bool comp(std::map<int,Measure>& sig1, std::map<int,Measure>& sig2)
{
    if(sig1.size() != sig2.size())
    {
        return false;
    }
    auto it1 = sig1.begin();
    auto it2 = sig2.begin();
    while(it1 != sig1.end() && it2 != sig2.end())
    {
        if(it1->first != it2->first)
        {
            return false;
        }
        
        if(!(it1->second).eq(it2->second))
        {
            return false;
        }
        
        it1++;
        it2++;
    }
    return true;
}


void spm(Graph& g)
{
    Measure::init(g.get_priority_map());
    std::map<int,Measure> sig;

    
    for(auto n : g.getNodes())
    {

        Measure i(std::list<int>(Measure::getSize(),0));
        sig[n.get_id()] = i;
    }
    std::map<int,Measure> prev_sig;
    while(!comp(prev_sig,sig))
    {

        prev_sig = sig;
        //@todo add other lifting techniques here
        for(Node& n : g.getNodes())
            lift(n, sig, g);
        //display(sig);
        //display(prev_sig);
    }
    bool sol = false;
    for(auto it = sig.begin(); it != sig.end(); it++)
    {
        if(!(it->second).get_tau())
        {
            std::cout<<"Solution for even from node "<<it->first<<std::endl;
            sol = true;
        }
    }
    if(!sol)
    {
        std::cout<<"No solutions for player even"<<std::endl;
    }

}


bool Measure::initialized = false;
int Measure::size = 0;
Measure* Measure::max_mes = nullptr;

int
main(int argc, const char * argv[]) {

    std::string file = argv[1];
    
    std::cout << "Parsing from ........ " << file << std::endl;
    auto G = init_graph_from_file(file);
    
    display_graph(&G);
    
    spm(G);
    
    return 0;
}