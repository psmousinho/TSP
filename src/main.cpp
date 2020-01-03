#include "data.h"
#include "hungarian.h"
#include <algorithm>
#include <iostream>
#include <numeric>
#include <vector>
#include <queue>
#include <chrono>

using namespace std;

Data* instance;
double** costMat;

void initSol(vector<int>&, double&);
void BB_best(vector<int>&, double&);
void BB_breadth(vector<int>&, double&);
void BB_depth(vector<int>&, double&);

void printTour(vector<int>&, bool);
void printMatrizAdj();


typedef struct node {
	vector<pair<int,int>> arcosProibidos;
	vector<vector<int>> subTours;
	double lb;
	int escolhido;
	bool podar;

	bool operator<(const node& rhs) const {
		return lb < rhs.lb;
	}

	void solve() {

		//Cria o problema
		hungarian_problem_t prob;
		int mode = HUNGARIAN_MODE_MINIMIZE_COST;
		int pSize = instance->getDimension();
		hungarian_init(&prob, costMat, pSize, pSize, mode);
		
		//Bloquia os pares no problema
		for(pair<int, int> p : arcosProibidos) {
			prob.cost[p.first][p.second] = INFINITE;
			//prob.cost[p.second][p.first] = INFINITE;
		}

		//Resolve o Problema
		lb = hungarian_solve(&prob);
		
		//Descobre os subtours retornados pelo alg. hungaro
		bool visited[pSize] = {false};
		int count = 0, current = 0;
		while(count < pSize) {
			for(int i = 0; i < pSize; i++) {
				if(!visited[i]) {
					current = i;
					break;
				}
			}

			std::vector<int> cicle = {current};

			while(!visited[current]) {
				for(int i = 0; i < pSize; i++) {
					//if(i == current) continue;
					if(prob.assignment[current][i] == 1) {
						visited[current] = true;
						cicle.push_back(i);
						current = i;
						break;
					}
				}
			}
			
			subTours.push_back(cicle);
	        count += cicle.size() -1;;
    	}	

		hungarian_free(&prob);

		//decide se vai podar ou nao
		podar = (subTours.size() == 1);
		
		//escolhe qual subcliclo bloquar nos filhos
		if(!podar) {
			int size = INFINITE;
			escolhido = 0;
			for(size_t i = 0; i < subTours.size(); i++) {
				if(subTours[i].size() < size) {
					size = subTours[i].size();
					escolhido = i;
				}
			}
		}

	}

	void printNode() {
		clog << endl;
		clog << "SubTours:\n";
		for(vector<int> s : subTours){
			printTour(s,false);
		}

		clog << "Proibidos: ";
		for(pair<int,int> p : arcosProibidos){
			clog << "(" << p.first << "," << p.second << ")," ;
		}
		clog << endl;
		
		clog << "lb: " << lb;
		clog << " escolhido: " << escolhido;
		clog << " podar: " << podar << endl;
	}
	
} Node;


int main(int argc, char** argv) {
	instance = new Data(argc, argv[1]);
	instance->readData();

	costMat = new double*[instance->getDimension()];
	for (int i = 0; i < instance->getDimension(); i++){
		costMat[i] = new double[instance->getDimension()];
		for (int j = 0; j < instance->getDimension(); j++){
			costMat[i][j] = instance->getDistance(i,j);
		}
	}
	
	vector<int> sol;
	double cost = 0;
	
	auto begin = chrono::steady_clock::now();
	BB_best(sol,cost);
	auto end = chrono::steady_clock::now();

	cout <<
	"INSTANCE: " << argv[1] << endl << 
	"COST: " << cost << endl <<
	"DURANTION(secs): " << chrono::duration_cast<chrono::seconds>(end-begin).count() << endl <<
	"ROUTE: ";
	printTour(sol,true);
	
	for (int i = 0; i < instance->getDimension(); i++) delete [] costMat[i];
	delete [] costMat;
	delete instance;
	
	return 0;	
}

void initSol(vector<int> &sol, double &cost) {
	int dimension = instance->getDimension();
	
	sol = {0,0};
	cost = 0;
	vector<int> candidates(dimension-1);
	iota(candidates.begin(),candidates.end(), 1);

	for(int i = 0; i < 3; i++) {
		int j = rand()% candidates.size();
		cost += costMat[sol[0]][candidates[j]] + costMat[candidates[j]][sol[1]]; 
		if(i != 0) cost -= costMat[sol[0]][sol[1]];
		sol.insert(sol.begin() + 1, candidates[j]);
		candidates.erase(candidates.begin() + j);
	}

	while(!candidates.empty()) {

		int candidatesSize = candidates.size();
		double delta = numeric_limits<double>::infinity();
		int pos1 = 0,pos2 = 0;
		for(size_t pos = 1; pos < sol.size(); pos++) {
			double aux = - costMat[sol[pos-1]][sol[pos]];
			for(int k = 0; k < candidatesSize; k++){
				double newDelta = costMat[sol[pos-1]][candidates[k]] + costMat[sol[pos]][candidates[k]] + aux;
				if(newDelta < delta) {
					pos1 = pos;
					pos2 = k;
					delta = newDelta;
				}
			}
		}

		sol.insert(sol.begin() + pos1, candidates[pos2]);
		cost += delta;
		candidates.erase(candidates.begin() + pos2);
  	}
}

void BB_best(vector<int> &bestSol, double &ub) {
	//Heuristic solution
	initSol(bestSol, ub);
	
	//Nodes queue sort by lowerBound
	priority_queue<Node> queue;
	
	//Init queue
	Node root;
	root.solve();
	queue.push(root);
	
	//Branch and Bound
	while(!queue.empty()) {
		
		Node node = queue.top();
		queue.pop();

		if(node.podar) {
			//new cheaper solution
			if(node.lb < ub) {
				bestSol = node.subTours[0];
				ub = node.lb;
				node.printNode();
			}

		} else {
			//branch
			for(size_t i = 0; i < node.subTours[node.escolhido].size() -1; i++) {
				Node newNode;
				
				newNode.arcosProibidos = node.arcosProibidos;
				pair<int,int> newArcoProibido(node.subTours[node.escolhido][i],node.subTours[node.escolhido][i+1]);	
				newNode.arcosProibidos.push_back(newArcoProibido);		

				newNode.solve();

				//bound
				if(newNode.lb < ub){
					queue.push(newNode);
				}
			}
		}

	}

}

void BB_breadth(vector<int> &bestSol, double &ub) {
	//Heuristic solution
	initSol(bestSol, ub);
	
	//Nodes queue sort by lowerBound
	vector<Node> queue;
	
	//Init queue
	Node root;
	root.solve();
	queue.push_back(root);
	
	//Branch and Bound
	while(!queue.empty()) {
		cout << queue.size() << endl;
		Node node = queue.front();
		queue.erase(queue.begin());

		if(node.podar) {
			//new cheaper solution
			if(node.lb < ub) {
				bestSol = node.subTours[0];
				ub = node.lb;
				node.printNode();
			}

		} else {
			//branch
			for(size_t i = 0; i < node.subTours[node.escolhido].size() -1; i++) {
				Node newNode;
				
				newNode.arcosProibidos = node.arcosProibidos;
				pair<int,int> newArcoProibido(node.subTours[node.escolhido][i],node.subTours[node.escolhido][i+1]);	
				newNode.arcosProibidos.push_back(newArcoProibido);		

				newNode.solve();

				//bound
				if(newNode.lb < ub){
					queue.push_back(newNode);
				}
			}
		}

	}

}

void BB_depth(vector<int> &bestSol, double &ub) {
	//Heuristic solution
	initSol(bestSol, ub);
	
	//Nodes queue sort by lowerBound
	vector<Node> queue;
	
	//Init queue
	Node root;
	root.solve();
	queue.push_back(root);
	
	//Branch and Bound
	while(!queue.empty()) {
		
		Node node = queue.back();
		queue.pop_back();

		if(node.podar) {
			//new cheaper solution
			if(node.lb < ub) {
				bestSol = node.subTours[0];
				ub = node.lb;
				node.printNode();
			}

		} else {
			//branch
			for(size_t i = 0; i < node.subTours[node.escolhido].size() -1; i++) {
				Node newNode;
				
				newNode.arcosProibidos = node.arcosProibidos;
				pair<int,int> newArcoProibido(node.subTours[node.escolhido][i],node.subTours[node.escolhido][i+1]);	
				newNode.arcosProibidos.push_back(newArcoProibido);		

				newNode.solve();

				//bound
				if(newNode.lb < ub){
					queue.push_back(newNode);
				}
			}
		}

	}

}

void printTour(vector<int> &tour, bool out) {
	if(out) {
		size_t i;
		cout << "{";
		for(i = 0; i < tour.size()-1;i++) {
			cout << tour[i] << ",";
		}
		cout << tour[i]<< "}\n\n"; 
	} else {
	 	size_t i;
		clog << "{";
		for(i = 0; i < tour.size()-1;i++) {
			clog << tour[i] << ",";
		}
		clog << tour[i]<< "}\n\n";
  	}
}

void printMatrizAdj() {
  for (size_t i = 0; i < instance->getDimension(); i++) {
    for (size_t j = 0; j < instance->getDimension(); j++) {
      cout << costMat[i][j] << " ";
    }
    cout << endl;
  }
}