#ifndef HEUCHH_H
#define HEUCHH_H

#ifdef LINUX
#include "cplex.h"
#else
#include <ilcplex/cplex.h>
#endif
#include <iostream>
#include <string>
#include <vector>
#include <algorithm> // sort, shuffle
#include <random>    // std::default_random_engine

#include "QGAP.h"

using namespace std;

class HeuCHH
{
public:
   QuadraticGAP* QGAP;

   HeuCHH(QuadraticGAP*);
   ~HeuCHH();
   void goCHH(CPXENVptr env, CPXLPptr lp, int maxnodes, int optimality_target, vector<int> sol);

private:
};


#endif // HEUCHH_H
