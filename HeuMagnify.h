#ifndef HEUMAG_H
#define HEUMAG_H

#ifdef LINUX
#include "cplex.h"
#else
#include <ilcplex/cplex.h>
#endif
#include <iostream>
#include <string>
#include <vector>
#include <algorithm> // sort

#include "QGAP.h"

using namespace std;

class HeuMagnify
{
public:
   QuadraticGAP* QGAP;

   HeuMagnify(QuadraticGAP*);
   ~HeuMagnify();
   void MagniGlass(CPXENVptr env, CPXLPptr lp, int maxnodes, int optimality_target, vector<int> sol);

private:
   double simpleContruct(double* x);
   void HeuMagnify::computeRegrets(double** c, int n, int m, vector<int> & regrets);
   double HeuMagnify::computeCost(double* x, int n, int m);
   double HeuMagnify::computeCost(int* sol, int n, int m);
   double HeuMagnify::fixVars(int * cnt, int m, int n, double *x, vector<int> &v_indices, vector<char> &v_lu, vector<double> &v_bd);
};


#endif // HEUMAG_H
