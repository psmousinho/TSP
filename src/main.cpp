#include "readData.h"
#include <algorithm>
#include <iostream>
#include <numeric>
#include <vector>
#include <cmath>
#include <tuple>
#include <chrono>

using namespace std;

double *** subTourData; // dados auxiliares para calculo de custos
double ** matrizAdj; // matriz de adjacencia
int dimension; // quantidade total de vertices

vector<int> GILS_RVND(double&, int, int, double, int);
void initSol(vector<int>& ,double&, int, double);
void RVND(vector<int>&, double&);
void swap(vector<int>&, double&);
void two_opt(vector<int>&, double&);
void reinsertion(vector<int>&, double&, int);
void doubleBridge(vector<int>&, double&, const vector<int>&, double);

double calCost(vector<int> &);
inline void initSubTourData();
inline void updateSubTourData(vector<int>&);

void printMatrizAdj();
void printSubTourMat();
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
  
  //initializing dimension, matrizADJ
  readData(argv[1], &dimension, &matrizAdj);
  //initializing subTourData
  initSubTourData();
  
  auto begin = chrono::steady_clock::now();

  double bestCost;
  vector<int> bestSol;
  if(dimension >= 150)
    bestSol = GILS_RVND(bestCost,50,dimension/2,0.5,3);
  else
    bestSol = GILS_RVND(bestCost,50,dimension,0.5,3);
  
  auto end = chrono::steady_clock::now();

  cout << 
  "INSTANCE: " << argv[1] << endl <<
  "COST: " << bestCost << endl <<
  "SEED: " << seed << endl <<
  "DURATION(secs): " << chrono::duration_cast<chrono::seconds>(end - begin).count() << endl <<
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
    clog << "reiniciando GRASP " << i+1 << "/" << Ig << endl;
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
      clog << "Atualizando melhor custo: " << bestCost << endl;
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
    candidates.erase(candidates.begin() + get<2>(costInsert[index]));
  }
  cost = calCost(sol);
}

void RVND(vector<int> &sol, double &cost) {
  vector<int> neighbors = {1,2,3,4,5};

  updateSubTourData(sol);

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
      updateSubTourData(sol);
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
    for(int j = i+1; j < solSize-1; j++) {
      double C, T;

      if(i != j-1){
        C = subTourData[0][i-1][1] + subTourData[j][j][2] * (subTourData[0][i-1][1] + matrizAdj[sol[i-1]][sol[j]]) + subTourData[j][j][1];
        T = subTourData[0][i-1][0] + matrizAdj[sol[i-1]][sol[j]] + subTourData[j][j][0];

        C = C + subTourData[i+1][j-1][2] * (T + matrizAdj[sol[j]][sol[i+1]]) + subTourData[i+1][j-1][1];
        T = T + matrizAdj[sol[j]][sol[i+1]] + subTourData[i+1][j-1][0];

        C = C + subTourData[i][i][2] * (T + matrizAdj[sol[j-1]][sol[i]]) + subTourData[i][i][1];
        T = T + matrizAdj[sol[j-1]][sol[i]] + subTourData[i][i][0];

        C = C + subTourData[j+1][solSize-1][2] * (T + matrizAdj[sol[i]][sol[j+1]]) + subTourData[j+1][solSize-1][1];
        
      } else {
        C = subTourData[0][i-1][1] + subTourData[j][i][2] * (subTourData[0][i-1][0] + matrizAdj[sol[i-1]][sol[j]]) + subTourData[j][i][1];
        T = subTourData[0][i-1][0] + matrizAdj[sol[i-1]][sol[j]] + subTourData[j][i][0];

        C = C + subTourData[j+1][solSize-1][2] * (T + matrizAdj[sol[j+1]][sol[i]]) + subTourData[j+1][solSize-1][1];

      }
      if(C < delta){
        pos1 = i;
        pos2 = j;
        delta = C;
      }
    }
  }

  if(delta < cost) {
    int aux = sol[pos1];
    sol[pos1] = sol[pos2];
    sol[pos2] = aux;
    cost = delta;
  }
}

void two_opt(vector<int> &sol, double &cost) {
  int pos1, pos2;
  double delta = numeric_limits<double>::infinity();

  int solSize = sol.size();
  for(int i = 1; i < solSize-1; i++) {
    for(int j = i+1; j < solSize-1; j++) {

      double C = subTourData[0][i-1][1] + subTourData[j][i][2] * (subTourData[0][i-1][0] + matrizAdj[sol[i-1]][sol[j]]) + subTourData[j][i][1];
      double T = subTourData[0][i-1][0] + matrizAdj[sol[i-1]][sol[j]] + subTourData[j][i][0];
      
      C = C + subTourData[j+1][solSize-1][2] * (T + matrizAdj[sol[i]][sol[j+1]]) + subTourData[j+1][solSize-1][1];
      
      if(C < delta){
        pos1 = i;
        pos2 = j;
        delta = C;
      }
    }
  }

  if(delta < cost) {
    reverse(sol.begin() + pos1, sol.begin() + pos2+1);
    cost = delta;
  }
}

void reinsertion(vector<int> &sol, double &cost, int tam) {
  int pos1, pos2;
  double delta = numeric_limits<double>::infinity();

  int solSize = sol.size(); 
  for(int i = 1; i < solSize-tam; i++) {
    for(int j = 1; j < solSize; j++) {
      if(j >= i && j <= i+tam)
        continue;
      
      double C;
      double T;
      if(i < j) {
        C = subTourData[0][i-1][1] + subTourData[i+tam][j-1][2] * ( subTourData[0][i-1][0] + matrizAdj[sol[i-1]][sol[i+tam]] )+ subTourData[i+tam][j-1][1];
        T = subTourData[0][i-1][0] + matrizAdj[sol[i-1]][sol[i+tam]] + subTourData[i+tam][j-1][0];

        C = C + subTourData[i][i+tam-1][2] * ( T + matrizAdj[sol[j-1]][sol[i]] )+ subTourData[i][i+tam-1][1];
        T = T + matrizAdj[sol[j-1]][sol[i]] + subTourData[i][i+tam-1][0];

        C = C + subTourData[j][solSize-1][2] * ( T + matrizAdj[sol[i+tam-1]][sol[j]] )+ subTourData[j][solSize-1][1];
        
      } else {
        C = subTourData[0][j-1][1] + subTourData[i][i+tam-1][2] * ( subTourData[0][j-1][0] + matrizAdj[sol[j-1]][sol[i]] )+ subTourData[i][i+tam-1][1];
        T = subTourData[0][j-1][0] + matrizAdj[sol[j-1]][sol[i]] + subTourData[i][i+tam-1][0];

        C = C + subTourData[j][i-1][2] * ( T + matrizAdj[sol[i+tam-1]][sol[j]] )+ subTourData[j][i-1][1];
        T = T + matrizAdj[sol[i+tam-1]][sol[j]] + subTourData[j][i-1][0];

        C = C + subTourData[i+tam][solSize-1][2] * ( T + matrizAdj[sol[i-1]][sol[i+tam]] )+ subTourData[i+tam][solSize-1][1];
      }
      
      if(C < delta) {
        pos1 = i;
        pos2 = j;
        delta = C;
      }

    }
  }

  if(delta < cost) {
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

    cost = delta;
  }
}

void doubleBridge(vector<int> &newSol, double &newCost, vector<int> &oldSol, double oldCost) {
  updateSubTourData(oldSol);
  
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

  double C = subTourData[0][pos1-1][1] + subTourData[pos3][oldSolSize-2][2] * ( subTourData[0][pos1-1][0] + matrizAdj[oldSol[pos1-1]][oldSol[pos3]] )+ subTourData[pos3][oldSolSize-2][1];
  double T = subTourData[0][pos1-1][0] + matrizAdj[oldSol[pos1-1]][oldSol[pos3]] + subTourData[pos3][oldSolSize-2][0];
  
  C = C + subTourData[pos2][pos3-1][2] * ( T + matrizAdj[oldSol[oldSolSize-2]][oldSol[pos2]] )+ subTourData[pos2][pos3-1][1];
  T = T + matrizAdj[oldSol[oldSolSize-2]][oldSol[pos2]] + subTourData[pos2][pos3-1][0];
  
  C = C + subTourData[pos1][pos2-1][2] * ( T + matrizAdj[oldSol[pos3-1]][oldSol[pos1]] )+ subTourData[pos1][pos2-1][1];
  T = T + matrizAdj[oldSol[pos3-1]][oldSol[pos1]] + subTourData[pos1][pos2-1][0];

  C = C + subTourData[oldSolSize-1][oldSolSize-1][2] * ( T + matrizAdj[oldSol[pos2-1]][oldSol[0]] )+ subTourData[oldSolSize-1][oldSolSize-1][1];

  newCost = C;
}

double calCost(vector<int> &tour) {
  double cost = 0;
  double delay = 0;
  for(int i = 0; i < tour.size() -1 ; i++) {
    cost += delay + matrizAdj[tour[i]][tour[i+1]];
    delay += matrizAdj[tour[i]][tour[i+1]];
  }
  return cost;
}

inline void initSubTourData() { 
  subTourData = new double** [dimension+1];
  for(int i = 0; i < dimension+1; i++){
    subTourData[i] = new double* [dimension+1];
    for(int j = 0; j < dimension+1; j++) {
      subTourData[i][j] = new double[3];
    }
  }

  //origin
  subTourData[0][0][0] = 0; 
  subTourData[0][0][1] = 0;
  subTourData[0][0][2] = 0;
  //diagonal
  for(int i = 1; i <= dimension; i++) {
    subTourData[i][i][0] = 0;
    subTourData[i][i][1] = 0;
    subTourData[i][i][2] = 1;
  }
}

inline void updateSubTourData(vector<int> &tour) {
  int size = tour.size();

  for(int s = 2; s <= size; s++) {
    for(int i = 0; i < size - s +1; i++) {
      int j = i + s -1;
      //durantion
      subTourData[i][j][0] = subTourData[i][j-1][0] + matrizAdj[tour[j-1]][tour[j]] + subTourData[j][j][0];
      //cost
      subTourData[i][j][1] = subTourData[i][j-1][1]+ subTourData[j][j][2] * (subTourData[i][j-1][0] + matrizAdj[tour[j-1]][tour[j]]) + subTourData[j][j][1];
      //delay
      subTourData[i][j][2] = subTourData[i][j-1][2] + subTourData[j][j][2];

      //reverse durantion
      subTourData[j][i][0] = subTourData[j][j][0] + matrizAdj[tour[j-1]][tour[j]] + subTourData[j-1][i][0];
      //reverse cost
      subTourData[j][i][1] = subTourData[j][j][1]+ subTourData[j-1][i][2] * (subTourData[j][j][0] + matrizAdj[tour[j-1]][tour[j]]) + subTourData[j-1][i][1];
      //reverse delay
      subTourData[j][i][2] = subTourData[j][j][2] + subTourData[j-1][i][2];
    }
  }

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

void printSubTourMat() {
  for (size_t i = 0; i <= dimension; i++) {
    for (size_t j = 0; j <= dimension; j++) {
      cout <<"("<< subTourData[i][j][0] <<","<< subTourData[i][j][1] <<","<< subTourData[i][j][2] << ") ";
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