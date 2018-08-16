#ifndef QGAP_H
#define QGAP_H

#ifdef LINUX
#include "/opt/ibm/ILOG/CPLEX_Studio128/cplex/include/ilcplex/cplex.h"
#else
#include <ilcplex/cplex.h>
#endif
#include <iostream>
#include <fstream>   // ifstream
#include <string>
#include <vector>

using namespace std;

class Config
{
public:
   string datapath;  // path to instances directory
   string datafile;  // the instance to solve
   string initsolfile;  // file with seed solutions
   int mode;         // usage mode: 1 optimize, 2 transcode
   int opt_target;   // optimality target
   int isverbose;    // console output detail
   int maxnodes;     // max num of nodes expanded in the tree search
   int maxiter;      // max num of iteration of the heuristic
   int selrule;      // rule to select variables to fix
   int fixperc;      // percentage of variables to fix at each MG iteration
};

class QuadraticGAP
{
public:
   QuadraticGAP();
   ~QuadraticGAP();

   string name;  // instance name
   int n;        // number of clients
   int m;        // number of servers
   //int v;        // unit traffic cost    <- factorize this into cqd
   double** cl;  // linear assignment costs
   double** cqf; // quadratic assignment costs, flows
   double** cqd; // quadratic assignment costs, distances
   int** req;    // client requests
   int*  cap;    // server capacities
   double EPS = 0.0001;

   double zub;                // best upper bound cost
   vector<int> solbest;       // best upper bound solution
   vector<vector<int>> zgl;   // linear cost estimations a' la gilmore-lawler

   Config* conf;

   int Qopt(void);
   int checkfeas(double* x, double solcost);
   int checkfeas(int* x, double solcost);
   int GLcosts();
   void saveZUB(int* sol, double solcost);
   void saveZUB(double* x, double solcost);

protected:
private:
   int setproblemdata(char **probname_p, int *numcols_p, int *numrows_p,
      int *objsen_p, double **obj_p, double **rhs_p,
      char **sense_p, int **matbeg_p, int **matcnt_p,
      int **matind_p, double **matval_p, double **lb_p,
      double **ub_p, char **ctype_p, int **qmatbeg_p, int **qmatcnt_p,
      int **qmatind_p, double **qmatval_p);
   double eigenValues(double *qmatval, int n);


   int optimality_target;  // convex function or not

};

// free
void free_and_null(char **ptr);

#endif // QGAP_H
