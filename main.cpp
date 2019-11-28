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
vector<int> initSol(double&, int, double);
vector<int> RVND(double&, vector<int>, double);
vector<int> swap(double&, vector<int>, double);
vector<int> reinsertion(double&, vector<int>, double);
vector<int> two_opt(double&, vector<int>, double);
vector<int> doubleBridge(double&, const vector<int>&, double);

void printMatrizAdj();
void printRoute(vector<int>);

int main(int argc, char** argv) {
  if (argc < 2) {
   cout << "\nFaltando parametros\n";
   cout << " ./exec [Instancia] <Seed> "<< endl;
   exit(1);
  }
  if (argc > 3) {
    cout << "\nMuitos parametros\n";
    cout << " ./exec [Instancia] <Seed>" << endl;
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
  
  auto start = chrono::steady_clock::now();

  double bestCost;
  vector<int> bestSol;
  if(dimension >= 150)
    bestSol = GILS_RVND(bestCost,50,dimension/2,0.5,3);
  else
    bestSol = GILS_RVND(bestCost,50,dimension,0.5,3);

  auto time = chrono::steady_clock::now() - start;
  
  cout << 
  "INSTANCE: " << argv[1] << endl <<
  "COST: " << bestCost << endl <<
  "SEED: " << seed << endl <<
  "EXECUTION TIME: " << time.count() << endl <<
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
    //cout << "reiniciando GRASP " << i << "/" << Ig << endl;
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

    if(cost < bestCost) {
      bestSol = sol;
      bestCost = cost;
      //cout << "Atualizando melhor custo: " << bestCost << endl;
    }
  }

  return bestSol;
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
    vector<tuple<double,int,int>> costInsert( (sol.size()-1) * candidates.size());

    for(int pos = 1, l = 0; pos < sol.size(); pos++) {
      for(int k = 0; k < candidates.size(); k++){
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

  return sol;
}

vector<int> RVND(double &newCost, vector<int> sol, double cost) {
  vector<int> neighbors = {1,2,3};

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
      neighbors = {1,2,3};
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

vector<int> reinsertion(double &newCost, vector<int> sol, double oldCost) {
  int pos1, pos2;
  newCost = numeric_limits<double>::infinity();

  int k = (random() %3) + 1;
  for(int i = 1; i < sol.size()-k-1; i++) {
    for(int j = 1; j < sol.size(); j++) {
      if(j >= i && j <= i+k)
        continue;
      
      double newDelta = matrizAdj[sol[j-1]][sol[i]] + matrizAdj[sol[i+k-1]][sol[j]] + matrizAdj[sol[i-1]][sol[i+k]]
                      - matrizAdj[sol[i-1]][sol[i]] - matrizAdj[sol[i+k-1]][sol[i+k]] - matrizAdj[sol[j-1]][sol[j]];
      if(newDelta < newCost) {
        pos1 = i;
        pos2 = j;
        newCost = newDelta;
      }
    }
  }

  vector<int> vec;
  if(pos1 < pos2) {
    vec.insert(vec.end(),sol.begin()+pos1,sol.begin()+pos1+k);
    sol.erase(sol.begin()+pos1,sol.begin()+pos1+k);
    sol.insert(sol.begin()+pos2-k,vec.begin(),vec.end());
  } else {
    vec.insert(vec.end(),sol.begin()+pos1,sol.begin()+pos1+k);
    sol.erase(sol.begin()+pos1,sol.begin()+pos1+k);
    sol.insert(sol.begin()+pos2,vec.begin(),vec.end());
  }

  newCost += oldCost;
  return sol;
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
  vec.insert(vec.end(), sol[0]);

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
  int i;
  cout << "{";
  for(i = 0; i < route.size()-1;i++) {
    cout << route[i] << ",";
  }
  cout <<route[i]<< "}\n\n";
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
      double newDelta = matrizAdj[sol.back()][candidates[i]];
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

vector<int> randomNeighborhood(double &newCost, vector<int> sol) {
  random_shuffle(sol.begin()+1,sol.end()-1);

  newCost = 0;
  for(size_t i = 0; i < sol.size() -1; i++) {
    newCost += matrizAdj[sol[i]][sol[i+1]];
  }

  return sol;
}
