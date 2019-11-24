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
double bsCost = numeric_limits<double>::infinity(); //custo da melhor solucao

void GILS_RVND(int, int, double, int);
vector<int> initSol(double&, int, double);
vector<int> RVND(double&, vector<int>, double);
vector<int> swap(double&, vector<int>, double);
vector<int> reinsertion(double&, vector<int>, double);
vector<int> two_opt(double&, vector<int>, double);
vector<int> doubleBridge(double&, const vector<int>&, double);

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

  printRoute(bs);
  cout << "custo: " << bsCost << "|seed: " << seed << endl;

  /*double cost;
  vector<int> sol = initSol(cost,3,0.5);
  double newCost= cost;
  vector<int> sol2 = sol;
  while(true) {
    sol2 = RVND(newCost,sol2,newCost);
    cout << cost << " " << newCost << " " << calCost(sol2) << endl; 
    printRoute(sol2);
  }*/
  
  
  return 0;
}

void GILS_RVND(int Ig, int Iils, double alpha, int sizeInitSubTour) {
  for(int i = 0; i < Ig; i++) {
    double cost, newCost;
    vector<int> sol, newSol;

    sol = initSol(cost, sizeInitSubTour, alpha);
    newSol = sol;
    newCost = cost;
    
    for(int j = 0; j < Iils; j++) {
      newSol = RVND(newCost, newSol, newCost);
      if(newCost < cost) {
        sol = newSol;
        cost = newCost;
        j = 0;
      }
      newSol = doubleBridge(newCost, sol, cost);
    }

    if(cost < bsCost) {
      bs = sol;
      bsCost = cost;
    }
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

vector<int> RVND(double &newCost, vector<int> sol, double cost) {
  vector<int> neighbors = {1,2};
  
  while(!neighbors.empty()) {
    int index = rand()% neighbors.size();
    double moveCost;
    vector<int> newSol;
    switch (neighbors[index])  {
      case 1:
        newSol = swap(moveCost, sol, cost);
        break;
      case 2:
        newSol = two_opt(moveCost, sol, cost);
        break;
      case 3:
        newSol = reinsertion(moveCost, sol, cost);
        break;
    }

    if(moveCost < cost) {
      sol = newSol;
      cost = moveCost;
      neighbors = {1,2};
    } else {
      neighbors.erase(neighbors.begin()+index);
    }
  }

  newCost = cost;
  return sol;
}

vector<int> swap(double &newCost, vector<int> sol, double oldCost) {
  int pos1, pos2;
  newCost = numeric_limits<double>::infinity();
  
  for(int i = 1; i < sol.size()-1; i++) {
    for(int j = i+1; j < sol.size()-1; j++) {
      double newDelta;
      if(i != j-1){
        newDelta = matrizAdj[sol[i]][sol[j-1]] + matrizAdj[sol[i]][sol[j+1]] + matrizAdj[sol[j]][sol[i-1]] + matrizAdj[sol[j]][sol[i+1]]
                 - matrizAdj[sol[i]][sol[i-1]] - matrizAdj[sol[i]][sol[i+1]] - matrizAdj[sol[j]][sol[j-1]] - matrizAdj[sol[j]][sol[j+1]];
      } else {
        newDelta = matrizAdj[sol[i-1]][sol[j]] + matrizAdj[sol[j+1]][sol[i]] - matrizAdj[sol[i-1]][sol[i]] - matrizAdj[sol[j+1]][sol[j]];
      }
      if(newDelta < newCost){
        pos1 = i;
        pos2 = j;
        newCost = newDelta;
      }
    }
  }
  
  int aux = sol[pos1];
  sol[pos1] = sol[pos2];
  sol[pos2] = aux;

  newCost += oldCost;
  return sol;
}

vector<int> two_opt(double &newCost, vector<int> sol, double oldCost) {
  int pos1, pos2;
  newCost = numeric_limits<double>::infinity();
  
  for(int i = 1; i < sol.size()-1; i++) {
    for(int j = i+1; j < sol.size()-1; j++) {
      double newDelta =  matrizAdj[sol[i-1]][sol[j]] + matrizAdj[sol[j+1]][sol[i]] - matrizAdj[sol[i-1]][sol[i]] - matrizAdj[sol[j+1]][sol[j]];
      if(newDelta < newCost){
        pos1 = i;
        pos2 = j;
        newCost = newDelta;
      }
    }
  }
  
  reverse(sol.begin() + pos1, sol.begin() + pos2+1);
  
  newCost += oldCost;
  return sol;
}

//todo
vector<int> reinsertion(double &newCost, vector<int> sol, double oldCost) {
  
}

vector<int> doubleBridge(double &newCost, const vector<int> &sol, double oldCost) {
  int intervalSize = sol.size()/4;
  int pos1 = 1 + rand()%intervalSize;
  int pos2 = pos1 + 1 + rand()%intervalSize;
  int pos3 = pos2 + 1 + rand()%intervalSize;
  
  vector<int> vec;
  vec.insert(vec.end(),sol.begin(),sol.begin()+pos1);
  vec.insert(vec.end(),sol.begin()+pos3,sol.end()-1);
  vec.insert(vec.end(),sol.begin()+pos2,sol.begin()+pos3);
  vec.insert(vec.end(),sol.begin()+pos1,sol.begin()+pos2);
  vec.insert(vec.end(), 1);

  newCost = matrizAdj[sol[pos1-1]][sol[pos3]] + matrizAdj[sol[sol.size()-2]][sol[pos2]] + matrizAdj[sol[pos3-1]][sol[pos1]] + matrizAdj[sol[pos2-1]][sol[0]]
           -matrizAdj[sol[pos1-1]][sol[pos1]] - matrizAdj[sol[pos2-1]][sol[pos2]] - matrizAdj[sol[pos3-1]][sol[pos3]] - matrizAdj[sol[sol.size()-2]][sol[0]];
  
  newCost += oldCost;
  return vec;
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
  cout << "\n";
  for(int i = 0; i < route.size();i++) {
    cout << route[i] << ",";
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

/*#################################################################################################*/

vector<int> cheapsterInsertion(double &cost, int sizeInitSubTour, double alpha) {
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
    int pos, cand;
    double delta = numeric_limits<double>::infinity();
    for(int i = 1, l = 0; i < sol.size(); i++) {
      for(int j = 0; j < candidates.size(); j++){
        double newDelta = matrizAdj[sol[i-1]][candidates[j]] + matrizAdj[sol[i]][candidates[j]] - matrizAdj[sol[i-1]][sol[i]];
        if(newDelta < delta) {
          delta = newDelta;
          pos = i;
          cand = j;
        }
      }
    }  
    
    sol.insert(sol.begin() + pos, candidates[cand]);
    cost += delta;
    candidates.erase(candidates.begin() + cand);
      
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

vector<int> nearestNeighbor(double &cost) {
  vector<int> sol = {1};
  vector<int> candidates(dimension-1);
  iota(candidates.begin(),candidates.end(), 2);
  
  while(!candidates.empty()) {
    int cand;
    double delta = numeric_limits<double>::infinity();
    for(int i = 0; i < candidates.size(); i++) {
      double newDelta = matrizAdj[sol.back()][i];
      if(newDelta < delta) {
        delta = newDelta;
        cand = i;
      }
    }

    sol.push_back(candidates[cand]);
    cost += delta;
    candidates.erase(candidates.begin() + cand);
  }

  cost += matrizAdj[sol.back()][1];
  sol.push_back(1);
  
  return sol;
}
