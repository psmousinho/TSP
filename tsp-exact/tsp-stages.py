#!/usr/bin/python
import sys
import cplex
from cplex.exceptions import CplexError
from readData import readInstance
from time import time

def createProblem(v, pesos):
    
    prob = cplex.Cplex()
    prob.objective.set_sense(prob.objective.sense.minimize)

    for i in range(v):
        for j in range(v):
            if(i == j):
                continue
            for t in range(v):
               prob.variables.add(obj=[pesos[i][j]], lb=[0], ub=[1], types="I", names=["y_" + str(i+1) + "_" + str(j+1) + "_" + str(t+1)]) 

    for i in range(v):
        var_list = []
        coeff_list = []
        for j in range(v):
            if(i == j):
                continue
            for t in range(v):
                var_list.append("y_" + str(i+1) + "_" + str(j+1) + "_"  +str(t+1))
                coeff_list.append(1)
        prob.linear_constraints.add(lin_expr=[[var_list, coeff_list]], senses = "E", rhs=[1], names = ["GRAU_OUT_v" + str(i+1)])
    
    for i in range(v):
        var_list = []
        coeff_list = []
        for j in range(v):
            if(i == j):
                continue
            for t in range(v):
                var_list.append("y_" + str(j+1) + "_" + str(i+1) + "_" + str(t+1))
                coeff_list.append(1)
        prob.linear_constraints.add(lin_expr=[[var_list, coeff_list]], senses = "E", rhs=[1], names = ["GRAU_IN_v" + str(i+1)])
    
    for t in range(v):
        var_list = []
        coeff_list = []
        for i in range(v):
            for j in range(v):
                if(i == j):
                    continue
                var_list.append("y_" + str(i+1) + "_" + str(j+1) + "_" + str(t+1))
                coeff_list.append(1)
        prob.linear_constraints.add(lin_expr=[[var_list, coeff_list]], senses = "E", rhs=[1], names = ["EXC_t" + str(t+1)])

    for i in range(1,v):
        var_list = []
        coeff_list = []
        for j in range(v):
            if(i==j):
                continue
            for t in range(v):
                 var_list.extend( [ "y_" + str(i+1) + "_" + str(j+1) + "_" + str(t+1), "y_" + str(j+1) + "_" + str(i+1) + "_" + str(t+1) ] ) 
                 coeff_list.extend([(t+1),-(t+1)])
        prob.linear_constraints.add(lin_expr=[[var_list, coeff_list]], senses = "E", rhs=[1], names = ["T_IN_T_OUT_v" + str(i+1)])

    return prob

def main():
    
    try:
        v,pesos = readInstance(sys.argv[1])
        prob = createProblem(v,pesos)
        prob.write("modelo.lp")
        st = time()
        prob.solve()
        dur = time() - st
    except CplexError as exc:
        print(exc)
        return

    sol = []
    for t in range(v):
        for i in range(v):
            for j in range(v):
                if(i==j):
                    continue
                if(prob.solution.get_values("y_" + str(i+1) + "_" + str(j+1) + "_" + str(t+1)) > 0.9):
                    sol.append(i+1)
    
    if(sys.argv.__contains__("-w")):
        f = open("results","a")
        f.write("Instance: " + sys.argv[1]+"\n")
        f.write("Solution status: " + str(prob.solution.get_status()) + " : " + prob.solution.status[prob.solution.get_status()] + "\n")
        f.write("Solution value: " + str(prob.solution.get_objective_value()) + "\n")
        f.write("Solution: " + str(sol) + "\n")
        f.write("Duration(secs): " + str(dur) + "\n\n")
        f.close()
    else:
        print("\nSolution status = ", prob.solution.get_status(), ":", prob.solution.status[prob.solution.get_status()])
        print("Solution value  = " , prob.solution.get_objective_value())
        print("Solution:", sol)
        print("Duration(secs):", dur)


if __name__ == "__main__":
   main()
