#include <cfloat>
#include "HeuMagnify.h"

HeuMagnify::HeuMagnify(QuadraticGAP* QGAPinst)
{  QGAP = QGAPinst;
}

HeuMagnify::~HeuMagnify()
{
}

// magnifying glass vlsn procedure
void HeuMagnify::MagniGlass(CPXENVptr env, CPXLPptr lp, int maxnodes, int optimality_target, vector<int> sol)
{  int    i,j,n,m,status,solstat=0;
   double objval;
   double *x     = NULL;
   double *slack = NULL;
   vector<int> v_indices;
   vector<char> v_lu;
   vector<double> v_bd;

   cout << "Magnifying glass heuristic. Starting" << endl;
   n = QGAP->n;
   m = QGAP->m;

   int cur_numrows = CPXgetnumrows(env, lp);
   int cur_numcols = CPXgetnumcols(env, lp);
   slack = (double *)malloc(cur_numrows * sizeof(double));
   x     = (double *)malloc(cur_numcols * sizeof(double));   // initialize sol
   for(j=0;j<n;j++)  QGAP->solbest.push_back(-1);            // initialize ub sol

   if(sol.empty())
      objval = simpleContruct(x);
   else
   {
      objval = computeCost(&sol[0],n,m);
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

      status = fixVars(&cnt, m, n, x, v_indices, v_lu, v_bd);

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
            QGAP->zub = objval;
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

// Chooses to variables to fix
double HeuMagnify::fixVars(int * cnt, int m, int n, double *x, vector<int> &v_indices, vector<char> &v_lu, vector<double> &v_bd)
{  int i,j;
   double p;

   *cnt = 0;

   for(i=0;i<m;i++)
      for(j=0;j<n;j++)
      {  p = rand()/(1.0*RAND_MAX) ;
         if(p < 0.3)
         {
            v_indices.push_back(i*n+j);
            v_lu.push_back('B');				//both bounds will be set
            v_bd.push_back(x[i*n+j]);				
            (*cnt)++;
         }
         else
         {
            v_indices.push_back(i*n+j);
            v_lu.push_back('U');				//upper bound
            v_bd.push_back(1);				
            (*cnt)++;

            v_indices.push_back(i*n+j);
            v_lu.push_back('L');				//lower bound
            v_bd.push_back(0);				
            (*cnt)++;
         }
      }

   return *cnt;
}

// constructive: each at the less requiring fecility
double HeuMagnify::simpleContruct(double* x)
{  int i,ii,j,jj,m,n;
   double zub;

   m = QGAP->m;
   n = QGAP->n;
   vector<int> cost(m),capleft(m),indReq(m);
   vector<int> regrets(n),indCost(n),sol(n);
   auto compCost = [&cost](int a, int b){ return cost[a] < cost[b]; };           // ASC order
   auto compRegr = [&regrets](int a, int b){ return regrets[a] > regrets[b]; };  // DESC order

   for(i=0;i<m;i++) capleft[i] = QGAP->cap[i];
   computeRegrets(QGAP->cl,n,m,regrets);        // just for ordering in some way

   for(j=0;j<n;j++) indCost[j] = j;
   std::sort(indCost.begin(), indCost.end(), compRegr);

   zub = 0;
   for(jj=0;jj<n;jj++)
   {  j = indCost[jj];  // client order by the function below
      for(i=0;i<m;i++)
      {  cost[i]= QGAP->req[i][j];
         indReq[i] = i;
      }

      std::sort(indReq.begin(), indReq.end(), compCost);

      ii=0;
      while(ii<m)
      {  i=indReq[ii];
         if(capleft[i]>=QGAP->req[i][j])
         {  sol[j]=i;
            capleft[i] -= QGAP->req[i][j];
            break;
         }
         ii++;
      }
      if(ii==m)
      {  cout << "[SimpleConstruct] Error. ii="+ii << endl;
         zub = DBL_MAX;
         break;
      }
   }

   for(i=0;i<m;i++)
      for(j=0;j<n;j++)
      {  x[i*n+j] = 0;
         if(sol[j]==i)
            x[i*n+j]=1;
      }
   zub = computeCost(x,n,m);
   cout << "Construction terminated. Zub = " << zub << endl;
   if(zub<DBL_MAX)
      for(int k=0;k<n;k++) QGAP->solbest[k] = sol[k];

   return zub;
}

// computes linear assignment regrets for each client
void HeuMagnify::computeRegrets(double** c, int n, int m, vector<int> & regrets)
{  int i,j;
   double first,second;

   for(j=0;j<n;j++)
   {  first = second = INT_MAX;
      for(i=0;i<m;i++)
      {  if(c[i][j] < first)
         {  second = first;
            first  = c[i][j];
         }
         else if (c[i][j] < second)
            second  = c[i][j];
      }
      regrets[j] = second - first;
   }
}

// computes the cost from a 0/1 array
double HeuMagnify::computeCost(double* x, int n, int m)
{
   double cost = 0;
   int i, j, h, k;

   for (i = 0; i<m; i++)
      for (j = 0; j<n; j++)
      {
         cost += x[i*n + j] * QGAP->cl[i][j]; // linear component
         for (h = 0; h<m; h++)
            for (k = 0; k<n; k++)
               cost += QGAP->cqd[i][h] * QGAP->cqf[j][k] * x[i*n + j] * x[h*n + k];  // quadratic component
      }

   return cost;
}

// computes the cost from an int valued array
double HeuMagnify::computeCost(int* sol, int n, int m)
{
   double cost = 0;
   int i, j, h, k;

   for (j = 0; j<n; j++)
   {  i = sol[j];
      cost += sol[j] * QGAP->cl[i][j]; // linear component
      for (h = 0; h<m; h++)
         for (k = 0; k<n; k++)
            cost += QGAP->cqd[i][h] * QGAP->cqf[j][k] * sol[j] * sol[k];  // quadratic component
   }

   return cost;
}
