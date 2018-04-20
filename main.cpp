#include <iostream>
#include "QGAP.h"
#include "DataManager.h"

int main()
{  QuadraticGAP* QGAP = new QuadraticGAP();
   DataManager*  DM   = new DataManager();
   DM->QGAP = QGAP;
   DM->loadConfig();
   int res;

   switch(QGAP->conf->mode) 
   {
      case 1:     // optimizing instance
         DM->readJSONdata(QGAP->conf->datapath + QGAP->conf->datafile);
         QGAP->Qopt();
         break;
      case 2:     // transcoding instance format
         DM->transcode(QGAP->conf->datapath + QGAP->conf->datafile);
         break;
      case 3:     // checking ampl solution
      {  string str = QGAP->conf->datafile;
         DM->readJSONdata(QGAP->conf->datapath + str.replace(str.find("sol"), 3, "json"));
         res = DM->readAmplSol(QGAP->conf->datapath + QGAP->conf->datafile);
         break;
      }
      default:
         cout << "Usupported use mode" << endl;
         break;
   }

   cout << "\n<ENTER> to exit ..."; cin.get();

   delete QGAP;
   delete DM;
}
