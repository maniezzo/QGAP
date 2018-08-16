#include <cfloat>
#include "HeuCHH.h"

HeuCHH::HeuCHH(QuadraticGAP* QGAPinst)
{  QGAP = QGAPinst;
}

HeuCHH::~HeuCHH()
{
}

// CHH heuristic
void HeuCHH::goCHH(CPXENVptr env, CPXLPptr lp, int maxnodes, int optimality_target, vector<int> sol)
{  int    i,j,n,m,status,solstat=0;
   double objval;
   double *x     = NULL;
   double *slack = NULL;
   vector<int> v_indices;
   vector<char> v_lu;
   vector<double> v_bd;

   cout << "CHH heuristic. Starting" << endl;
   n = QGAP->n;
   m = QGAP->m;

   if(n>0)  // DA IMPLEMENTARE!!!
   {  cout << "TODO, exiting" << endl;
      return;
   }

   int cur_numrows = CPXgetnumrows(env, lp);
   int cur_numcols = CPXgetnumcols(env, lp);
   slack = (double *)malloc(cur_numrows * sizeof(double));
   x     = (double *)malloc(cur_numcols * sizeof(double));   // initialize sol
   for(j=0;j<n;j++)  QGAP->solbest.push_back(-1);            // initialize ub sol

   if(sol.empty() || sol[0] == -1)
      objval = -1;
   else
   {
      //objval = computeCost(&sol[0],n,m);
      for(i=0;i<m;i++)
         for (j = 0; j<n; j++)
            x[i*n + j] = 0;
      for (j = 0; j<n; j++)
         x[sol[j]*n + j] = 1;
   }

   status = CPXsetintparam(env, CPXPARAM_MIP_Limits_Nodes, maxnodes);   // max num of nodes
   if (status)
   {  fprintf(stderr, "Failure to reset max tree search nodes, error %d.\n", status);
      goto TERMINATE;
   }
   cout << "Max num of search nodes: " << maxnodes << endl;

   status = CPXsetintparam(env, CPXPARAM_OptimalityTarget, optimality_target);
   if (status)
   {  fprintf(stderr, "Failure to reset optimality target, error %d.\n", status);
      goto TERMINATE;
   }
   cout << "Optimality target set to "<< optimality_target << endl;

   if(objval > DBL_MAX - QGAP->EPS)    // in case no solution was constructed
   {
      // Optimize the problem and obtain solution.
      if(optimality_target > 1)
         status = CPXqpopt(env, lp); // in case of non convex function (opt target > 1)
      if (status) 
      {  fprintf(stderr, "Failed to optimize QP.\n");
      //goto TERMINATE;
      }

      int nn = CPXgetsolnpoolnumsolns(env, lp);
      if(nn==0)
      {  fprintf(stderr, "Failed to find feasible solutions.\n");
         goto TERMINATE;
      }

      double mincost = DBL_MAX;
      int minind = -1;
      for(i=0;i<nn;i++)
      {  status = CPXgetsolnpoolobjval (env, lp, i, &objval);
         cout << "Solution " << i << " cost " << objval << endl;
         if(objval<mincost)
         {  minind = i;
            mincost = objval;
         }
      }

      status = CPXgetsolnpoolx (env, lp, minind, x, 0, CPXgetnumcols(env, lp)-1);

      status = CPXsolution(env, lp, &solstat, &objval, x, NULL, slack, NULL);
      if (status) 
      {  fprintf(stderr, "Failed to obtain solution.\n");
         goto TERMINATE;
      }
   }

   status = QGAP->checkfeas(x, objval);
   if (status)
   {  fprintf(stderr, "Solution infeasible.\n");
      cout << "Solution infeasible !!! status = " << status << endl;
      //goto TERMINATE;
   }
   else cout << "Solution checked, status " << status << " cost " << objval << endl;

   // Write the output to the screen.
   cout << "\nConstruction solution status = " << solstat << endl;
   cout << "Solution value  = " << objval << endl;
   cout << "Solution:" << endl;
   for(i=0;i<cur_numcols;i++)
      cout << x[i] << " ";
   cout << endl;

   // output to the screen
   status = CPXsetintparam(env, CPXPARAM_ScreenOutput,QGAP->conf->isverbose); // CPX_ON, CPX_OFF

   // lu[j] 	= 'L' 	bd[j] is a lower bound
   // lu[j] 	= 'U' 	bd[j] is an upper bound
   // lu[j] 	= 'B' 	bd[j] is the lower and upper bound
   int  iter,cnt,maxiter = QGAP->conf->maxiter;
   double p;

   iter = 0;
   do
   {  cout << "iter " << iter << " zub " << QGAP->zub << endl;

      //status = fixVars(&cnt, m, n, x, v_indices, v_lu, v_bd);

      status = CPXchgbds (env, lp, cnt, &v_indices[0], &v_lu[0], &v_bd[0]);

      v_indices.clear();
      v_lu.clear();
      v_bd.clear();

      status = CPXqpopt(env, lp); // in case of non convex function (opt target > 1)
      status = CPXsolution(env, lp, &solstat, &objval, x, NULL, slack, NULL);
      if(status)
         cout << "Failed to obtain solution.\n";
      else
      {  if(objval < QGAP->zub)
         {  cout << "New zub!! " << objval << endl;
            QGAP->saveZUB(x,objval);
         }
      }
      iter++;
   }
   while (iter < maxiter);

   cout << "Solution cost " << objval << " zub " << QGAP->zub << endl;

TERMINATE:
   free_and_null((char **)&x);
   free_and_null((char **)&slack);
   return;
}
