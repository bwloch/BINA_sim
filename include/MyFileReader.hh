#include <list>
#include <fstream>
#include "G4LorentzVector.hh"
#include "G4ThreeVector.hh"

class MyFileReader
{
 public:
   MyFileReader(G4String fileName);
   ~MyFileReader();
   G4LorentzVector GetAnEvent();
 private:
   std::ifstream inputFile;
   std::list<G4LorentzVector> evList;
};
