#include "readData.h"
#include <algorithm>
#include <iostream>
#include <numeric>
#include <vector>
#include <cmath>
#include <tuple>
#include<chrono>

using namespace std;

double ** matrizAdj; // matriz de adjacencia
int dimension; // quantidade total de vertices

vector<int> GILS_RVND(double&, int, int, double, int);
void initSol(vector<int>& ,double&, int, double);
void RVND(vector<int>&, double&);
void swap(vector<int>&, double&);
void two_opt(vector<int>&, double&);
void reinsertion(vector<int>&, double&, int);
void doubleBridge(vector<int>&, double&, const vector<int>&, double);

void printMatrizAdj();
void printRoute(vector<int>);

int main(int argc, char** argv) {
  if (argc < 2) {
   cout << "\nFaltando parametros\n";
   cout << " ./exec Instancia [Seed] "<< endl;
   exit(1);
  }
  if (argc > 3) {
    cout << "\nMuitos parametros\n";
    cout << " ./exec Instancia [Seed]" << endl;
    exit(1);
  }
  
  long int seed;
  if(argc == 3) {
    seed = atol(argv[2]);
  } else {
    seed = time(nullptr);
  }
  srand(seed);

  readData(argv[1], &dimension, &matrizAdj);

  double bestCost;
  vector<int> bestSol;
  if(dimension >= 150)
    bestSol = GILS_RVND(bestCost,25,dimension/2,0.5,3);
  else
    bestSol = GILS_RVND(bestCost,50,dimension,0.5,3);
  
  cout << 
  "INSTANCE: " << argv[1] << endl <<
  "COST: " << bestCost << endl <<
  "SEED: " << seed << endl <<
  "ROUTE: ";
  printRoute(bestSol);
  
  return 0;
}

vector<int> GILS_RVND(double &bestCost, int Ig, int Iils, double alpha, int sizeInitSubTour) {
  vector<int> bestSol;
  bestCost = numeric_limits<double>::infinity();

  for(int i = 0; i < Ig; i++) {
    double cost, newCost;
    vector<int> sol, newSol;
    cout << "reiniciando GRASP " << i << "/" << Ig << endl;
    initSol(sol, cost, sizeInitSubTour, alpha);
    newSol = sol;
    newCost = cost;

    for(int j = 0; j < Iils; j++) {
      RVND(newSol, newCost);
      if(newCost < cost) {
        sol = newSol;
        cost = newCost;
        j = 0;
      }
      doubleBridge(newSol, newCost, sol, cost);
    }

    if(cost < bestCost) {
      bestSol = sol;
      bestCost = cost;
      cout << "Atualizando melhor custo: " << bestCost << endl;
      printRoute(bestSol);
    }
  }

  return bestSol;
}

void initSol(vector<int> &sol, double &cost, int sizeInitSubTour, double alpha) {
  sol = {1,1};
  cost = 0;
  vector<int> candidates(dimension-1);
  iota(candidates.begin(),candidates.end(), 2);

  for(int i = 0; i < sizeInitSubTour; i++) {
    int j = rand()% candidates.size();
    cost += matrizAdj[sol[0]][candidates[j]] + matrizAdj[candidates[j]][sol[1]] - matrizAdj[sol[0]][sol[1]];
    sol.insert(sol.begin() + 1, candidates[j]);
    candidates.erase(candidates.begin() + j);
  }

  while(!candidates.empty()) {
    vector<tuple<double,int,int>> costInsert( (sol.size()-1) * candidates.size());

    int candidatesSize = candidates.size();
    for(int pos = 1, l = 0; pos < sol.size(); pos++) {
      for(int k = 0; k < candidatesSize; k++){
        double delta = matrizAdj[sol[pos-1]][candidates[k]] + matrizAdj[sol[pos]][candidates[k]] - matrizAdj[sol[pos-1]][sol[pos]];
        costInsert[l++] = make_tuple(delta,pos,k);
      }
    }

    sort(costInsert.begin(), costInsert.end());
    int index = rand()%  (int) floor(alpha*costInsert.size());
    sol.insert(sol.begin() + get<1>(costInsert[index]), candidates[get<2>(costInsert[index])]);
    cost += get<0>(costInsert[index]);
    candidates.erase(candidates.begin() + get<2>(costInsert[index]));
  }
}

void RVND(vector<int> &sol, double &cost) {
  vector<int> neighbors = {1,2,3,4,5};

  while(!neighbors.empty()) {
    int index = rand()% neighbors.size();
    double oldCost = cost;

    switch (neighbors[index])  {
      case 1:
        swap(sol, cost);
        break;
      case 2:
        two_opt(sol, cost);
        break;
      case 3:
        reinsertion(sol, cost,1);
        break;
      case 4:
        reinsertion(sol, cost,2);
        break;
      case 5:
        reinsertion(sol, cost,3);
        break;
    }

    if(cost < oldCost) {
      neighbors = {1,2,3,4,5};
    } else {
      neighbors.erase(neighbors.begin()+index);
    }
  }
}

void swap(vector<int> &sol, double &cost) {
  int pos1, pos2;
  double delta = numeric_limits<double>::infinity();

  int solSize = sol.size();
  for(int i = 1; i < solSize-1; i++) {
    double aux1 = - matrizAdj[sol[i]][sol[i-1]] - matrizAdj[sol[i]][sol[i+1]];
    for(int j = i+1; j < solSize-1; j++) {
      double newDelta;
      if(i != j-1){
        newDelta = aux1 + matrizAdj[sol[i]][sol[j-1]] + matrizAdj[sol[i]][sol[j+1]] + matrizAdj[sol[j]][sol[i-1]] + matrizAdj[sol[j]][sol[i+1]]
                  - matrizAdj[sol[j]][sol[j-1]] - matrizAdj[sol[j]][sol[j+1]] ;
      } else {
        newDelta = matrizAdj[sol[i-1]][sol[j]] + matrizAdj[sol[j+1]][sol[i]] - matrizAdj[sol[i-1]][sol[i]] - matrizAdj[sol[j+1]][sol[j]];
      }
      if(newDelta < delta){
        pos1 = i;
        pos2 = j;
        delta = newDelta;
      }
    }
  }

  if(delta < 0) {
    int aux = sol[pos1];
    sol[pos1] = sol[pos2];
    sol[pos2] = aux;
    cost += delta;
  }
}

void two_opt(vector<int> &sol, double &cost) {
  int pos1, pos2;
  double delta = numeric_limits<double>::infinity();

  int solSize = sol.size();
  for(int i = 1; i < solSize-1; i++) {
    double aux = - matrizAdj[sol[i-1]][sol[i]];
    for(int j = i+1; j < solSize-1; j++) {
      double newDelta =  aux + matrizAdj[sol[i-1]][sol[j]] + matrizAdj[sol[j+1]][sol[i]] - matrizAdj[sol[j+1]][sol[j]];
      if(newDelta < delta){
        pos1 = i;
        pos2 = j;
        delta = newDelta;
      }
    }
  }

  if(delta < 0) {
    reverse(sol.begin() + pos1, sol.begin() + pos2+1);
    cost += delta;
  }
}

void reinsertion(vector<int> &sol, double &cost, int tam) {
  int pos1, pos2;
  double delta = numeric_limits<double>::infinity();

  int solSize = sol.size(); 
  for(int i = 1; i < solSize-tam; i++) {
    double aux = + matrizAdj[sol[i-1]][sol[i+tam]] - matrizAdj[sol[i-1]][sol[i]];
    for(int j = 1; j < solSize; j++) {
      if(j >= i && j <= i+tam)
        continue;
      
      double newDelta = aux + matrizAdj[sol[j-1]][sol[i]] + matrizAdj[sol[i+tam-1]][sol[j]] - matrizAdj[sol[i+tam-1]][sol[i+tam]] - matrizAdj[sol[j-1]][sol[j]];
      if(newDelta < delta) {
        pos1 = i;
        pos2 = j;
        delta = newDelta;
      }
    }
  }

  if(delta < 0) {
    vector<int> vec;
    if(pos1 < pos2) {
      vec.insert(vec.end(),sol.begin()+pos1,sol.begin()+pos1+tam);
      sol.erase(sol.begin()+pos1,sol.begin()+pos1+tam);
      sol.insert(sol.begin()+pos2-tam,vec.begin(),vec.end());
    } else {
      vec.insert(vec.end(),sol.begin()+pos1,sol.begin()+pos1+tam);
      sol.erase(sol.begin()+pos1,sol.begin()+pos1+tam);
      sol.insert(sol.begin()+pos2,vec.begin(),vec.end());
    }

    cost += delta;
  }
}

void doubleBridge(vector<int> &newSol, double &newCost, const vector<int> &oldSol, double oldCost) {
  int oldSolSize = oldSol.size();
  int intervalSize = oldSolSize/4;
  int pos1 = 1 + rand()%intervalSize;
  int pos2 = pos1 + 1 + rand()%intervalSize;
  int pos3 = pos2 + 1 + rand()%intervalSize;

  newSol.clear();
  newSol.insert(newSol.end(),oldSol.begin(),oldSol.begin()+pos1);
  newSol.insert(newSol.end(),oldSol.begin()+pos3,oldSol.end()-1);
  newSol.insert(newSol.end(),oldSol.begin()+pos2,oldSol.begin()+pos3);
  newSol.insert(newSol.end(),oldSol.begin()+pos1,oldSol.begin()+pos2);
  newSol.insert(newSol.end(), oldSol[0]);

  newCost = matrizAdj[oldSol[pos1-1]][oldSol[pos3]] + matrizAdj[oldSol[oldSolSize-2]][oldSol[pos2]] + matrizAdj[oldSol[pos3-1]][oldSol[pos1]] + matrizAdj[oldSol[pos2-1]][oldSol[0]]
           -matrizAdj[oldSol[pos1-1]][oldSol[pos1]] - matrizAdj[oldSol[pos2-1]][oldSol[pos2]] - matrizAdj[oldSol[pos3-1]][oldSol[pos3]] - matrizAdj[oldSol[oldSolSize-2]][oldSol[0]];

  newCost += oldCost;
}

/*#################################################################################################*/

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
  int i;
  cout << "{";
  for(i = 0; i < route.size()-1;i++) {
    cout << route[i] << ",";
  }
  cout <<route[i]<< "}\n\n";
}
