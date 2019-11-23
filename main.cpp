#include "readData.h"
#include <algorithm>
#include <iostream>
#include <numeric>
#include <vector>
#include <math.h>
#include <tuple>

using namespace std;

double ** matrizAdj; // matriz de adjacencia
int dimension; // quantidade total de vertices
vector<int> bs; //melhor solucao encontrada globalmente
int bsCost; //custo da melhor solucao

void GILS_RVND(int, int, double, int);
vector<int> initSol(double&, int, double);
vector<int> RVND(double&, vector<int>, int);
vector<int> swap(double&, vector<int>);
vector<int> reinsertion(double&, vector<int>);
vector<int> two_opt(double&, vector<int>);
void perturb(vector<int>&);

void printMatrizAdj();
void printRoute(vector<int>);
double calCost(vector<int>);


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
    printRoute(sol);
    perturb(sol);
    
    /*for(int j = 0; j < Iils; j++) {
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
    }*/
  
  }
}

vector<int> initSol(double &cost, int sizeInitSubTour, double alpha) {
  cost = 0;
  vector<int> sol = {1,1};
  vector<int> candidates(dimension-1);
  iota(candidates.begin(),candidates.end(), 2);

  for(int i = 0; i < sizeInitSubTour; i++) {
    int j = rand()% candidates.size();
    cost += matrizAdj[sol[0]][candidates[j]] + matrizAdj[candidates[j]][sol[1]] - matrizAdj[sol[0]][sol[1]];
    sol.insert(sol.begin() + 1, candidates[j]);
    candidates.erase(candidates.begin() + j);
  }
  
  while(!candidates.empty()) {
    vector<tuple<double,int,int>> costIncert( (sol.size()-1) * candidates.size());
    
    for(int pos = 1, l = 0; pos < sol.size(); pos++) {
      for(int k = 0; k < candidates.size(); k++){
        double delta = matrizAdj[sol[pos-1]][candidates[k]] + matrizAdj[sol[pos]][candidates[k]] - matrizAdj[sol[pos-1]][sol[pos]];
        costIncert[l++] = make_tuple(delta,pos,k);
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
      cout << "\n\n";
      cout <<  get<0>(costIncert[index]) << ";" << get<1>(costIncert[index]) << ";" << get<2>(costIncert[index]) << "||";
      cout << "\n";
      for(int i = 0; i < candidates.size();i++) {
        cout <<  candidates[i] << "|~|";
      }
      cout << "\n";
      printRoute(sol);
    #endif
  }

  return sol;
}

vector<int> RVND(double &newCost, vector<int> sol, int cost) {
  vector<int> neighbors = {1,2,3};
  
  while(!neighbors.empty()) {
    int index = rand()% neighbors.size();
    
    vector<int> newSol;
    double delta;
    switch (neighbors[index])  {
      case 1:
        newSol = swap(delta, sol);
        break;
      case 2:
        newSol = reinsertion(delta, sol);
        break;
      case 3:
        newSol = two_opt(delta, sol);
        break;
    }

    if(cost + delta < cost) {
      sol = newSol;
      cost = cost+delta;
      neighbors = {1,2,3};
    } else {
      neighbors.erase(neighbors.begin()+index);
    }

  }

  newCost = cost;
  return sol;
}

vector<int> swap(double &delta, vector<int> sol) {
  int pos1, pos2;
  delta = numeric_limits<double>::infinity();
  for(int i = 0; i < sol.size(); i++) {
    for(int j = i+1; j < sol.size(); j++) {
      
      double newDelta;
      if(i != j+1){
        newDelta = matrizAdj[i][j-1] + matrizAdj[i][j+1] + matrizAdj[j][i-1] + matrizAdj[j][i+1]
                 - matrizAdj[i][i-1] - matrizAdj[i][i+1] - matrizAdj[j][j-1] - matrizAdj[j][j+1];
      } else {
        newDelta = matrizAdj[i-1][j] + matrizAdj[j+1][i] - matrizAdj[i-1][i] - matrizAdj[j+1][j];
      }

      if(newDelta < delta){
        pos1 = i;
        pos2 = j;
        delta = newDelta;
      }

    }
  }
  
  int aux = sol[pos1];
  sol[pos1] = sol[pos2];
  sol[pos2] = aux;

  return sol;
}

//todo
vector<int> reinsertion(double &delta, vector<int> sol) {
  
}

vector<int> two_opt(double &delta, vector<int> sol) {
  int pos1, pos2;
  delta = numeric_limits<double>::infinity();
  for(int i = 1; i < sol.size()-1; i++) {
    for(int j = i+1; j < sol.size()-1; j++) {
      double newDelta =  matrizAdj[i-1][j] + matrizAdj[j+1][i] - matrizAdj[i-1][i] - matrizAdj[j+1][j];
      if(newDelta < delta){
        pos1 = i;
        pos2 = j;
        delta = newDelta;
      }
    }
  }

  reverse(sol.begin() + pos1, sol.begin() + pos2);
  return sol;
}

//todo
void perturb(vector<int> &sol) {
  int intervalSize = sol.size()/4;
  vector<int> first(sol.begin()+1, sol.begin()+intervalSize);
  vector<int> second(sol.begin()+intervalSize, sol.begin()+intervalSize*2);
  vector<int> third(sol.begin()+intervalSize*2, sol.begin()+intervalSize*3);
  vector<int> fourth(sol.begin()+intervalSize*3, sol.end()-1);

  printRoute(first);
  printRoute(second);
  printRoute(third);
  printRoute(fourth);
  
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

void printRoute(vector<int> route) {
  cout << "\n";
  for(int i = 0; i < route.size();i++) {
    cout << route[i] << "|";
  }
  cout << "\n";
}

double calCost(vector<int> route) {
  double cost = 0;
  for(size_t i = 0; i < route.size() -1; i++) {
    cost += matrizAdj[route[i]][route[i+1]];
  }
  return cost;
}
