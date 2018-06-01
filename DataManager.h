#ifndef DATA_H
#define DATA_H

#include <sstream>      // stringstream
#include <iomanip>      // std::setw
#include <functional>   // std::ptr_fun
#include <algorithm>    // std::find_if
#include <cctype>       // std::isspace

#include "json.h"
#include "QGAP.h"

class DataManager
{
public:
   DataManager();
   ~DataManager();

   Config* loadConfig();
   void readJSONdata(string);
   void transcode(string infile);
   int readAmplSol(string infile);
   vector<int> DataManager::readInitSol(string infile, string instance);

   QuadraticGAP* QGAP;

protected:
private:
   void json2ampl(string infile);
   void ampl2json(string infile);
   void leeMa2json(string infile);
};

#endif // DATA_H
