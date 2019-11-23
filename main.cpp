#include "readData.h"
#include <fstream>
#include <iostream>
#include <algorithm>
#include <numeric>
#include <vector>
#include <math.h>
#include <tuple>

using namespace std;

double ** matrizAdj; // matriz de adjacencia
int dimension; // quantidade total de vertices
vector<int> bs; //melhor solucao encontrada globalmente
int bsCost; //custo da melhor solucao

void printMatrizAdj();
void printRoute(vector<int>);
double calCost(vector<int>);
void GILS_RVND(int, int, double, int);
vector<int> initSol(double&, int, double);
vector<int> RVND(double&, vector<int>, int);
void perturb(vector<int>&);


int main(int argc, char** argv) {
  int seed = time(nullptr);
  srand(seed);

  readData(argc, argv, &dimension, &matrizAdj);
  
  if(dimension >= 150)
    GILS_RVND(50,dimension/2,0.5,3);
  else
    GILS_RVND(1,dimension,0.5,3);

  //printRoute(bs);
  //cout << "custo: " << bsCost << "|seed: " << seed << endl;

  return 0;
}

void GILS_RVND(int Ig, int Iils, double alpha, int sizeInitSubTour) {
  for(int i = 0; i < Ig; i++) {
    
    double cost;
    vector<int> sol = initSol( cost, sizeInitSubTour, alpha);
    //printRoute(sol);
    
    for(int j = 0; j < Iils; j++) {
      double newCost;
      vector<int> newSol = RVND(newCost, sol, cost);
      if(newCost < cost) {
        sol = newSol;
        cost = newCost;
        j = 0;
      }
      perturb(sol);
    }

    if(bs.empty() || cost < bsCost) {
      bs = sol;
      bsCost = cost;
    }

  }
}

vector<int> initSol(double &cost, int sizeInitSubTour, double alpha) {
  vector<int> sol = {1,1};
  vector<int> candidates(dimension-1);
  iota(candidates.begin(),candidates.end(), 2);

  for(int i = 0; i < sizeInitSubTour; i++) {
    int j = rand()% candidates.size();
    sol.insert(sol.begin() + 1, candidates[j]);
    candidates.erase(candidates.begin() + j);
  }
  cost = calCost(sol);
  
  while(!candidates.empty()) {
    vector< tuple<double,int,int> > costIncert;
    
    for(int pos = 1; pos < sol.size(); pos++) {
      for(int k = 0; k < candidates.size(); k++){
        double delta = matrizAdj[sol[pos-1]][candidates[k]] + matrizAdj[sol[pos]][candidates[k]] - matrizAdj[sol[pos-1]][sol[pos]];
        costIncert.push_back(make_tuple(delta,pos,k));
      }
    }  
    
    sort(costIncert.begin(), costIncert.end());
    int index = rand()%  (int) floor(alpha*costIncert.size());
    sol.insert(sol.begin() + get<1>(costIncert[index]), candidates[get<2>(costIncert[index])]);
    cost += get<0>(costIncert[index]);
    candidates.erase(candidates.begin() + get<2>(costIncert[index]));
      
    #ifdef DEBUG_INIT
      cout << "\n";
      for(int i = 0; i < costIncert.size();i++) {
        cout <<  get<0>(costIncert[i]) << ";" << get<1>(costIncert[i]) << ";" << get<2>(costIncert[i]) << "||";
      }
      cout << "\n";
      cout <<  get<0>(costIncert[index]) << ";" << get<1>(costIncert[index]) << ";" << get<2>(costIncert[index]) << "||";
      cout << "\n";
      for(int i = 0; i < candidates.size();i++) {
        cout <<  candidates[i] << "|~|";
      }
      cout << "\n";
      printRoute();
    #endif
  }

  return sol;
}

vector<int> RVND(double &newCost, vector<int> sol, int cost) {
  vector<int> neighbors = {1,2,3};
  
  while(!neighbors.empty()) {
    int index = rand()% neighbors.size();
    
    vector<int> newSol;
    switch (neighbors[index])  {
      case 1:
        //newSol = 2_opt(newCost);
        break;
      case 2:
        //newSol = swap(newCost);
        break;
      case 3:
        //newSol = reinsertion(newCost);
        break;
    }

    if(newCost < cost) {
      sol = newSol;
      cost = newCost;
      neighbors = {1,2,3};
    } else {
      neighbors.erase(neighbors.begin()+index);
    }

    newCost = cost;
    return sol;
  }

}

void perturb(vector<int> &sol) {

}

double calCost(vector<int> route) {
  double cost = 0;
  for(size_t i = 0; i < route.size() -1; i++) {
    cost += matrizAdj[route[i]][route[i+1]];
  }
  return cost;
}

void printRoute(vector<int> route) {
  cout << "\n";
  for(int i = 0; i < route.size();i++) {
    cout << route[i] << "|";
  }
  cout << "\n";
}

void printMatrizAdj() {
  cout << "dimension: " << dimension << endl;
  for (size_t i = 1; i <= dimension; i++) {
    for (size_t j = 1; j <= dimension; j++) {
      cout << matrizAdj[i][j] << " ";
    }
    cout << endl;
  }
}
