#include "MyFileReader.hh"
MyFileReader::MyFileReader(G4String fileName)
{ inputFile.open(fileName.data()); }

MyFileReader::~MyFileReader()
{ inputFile.close(); }

G4LorentzVector MyFileReader::GetAnEvent()
{
//G4cout<<"\t MyfileReader: GetAnEvent()\n";
  if( evList.empty() )
  {
    for(int i=0;i<100;i++)
    {
      G4double ee, ex, ey, ez;
      inputFile >> ee >> ex >> ey >> ez;
      
      evList.push_back( G4LorentzVector(G4ThreeVector(ex,ey,ez),ee) );
          
    }
    //G4cout<<"\t MyfileReader: GetAnEvent() -> Reading from file\n";
  }
  G4LorentzVector ev = evList.front();
  evList.pop_front();
  //G4cout<<"\t MyfileReader: GetAnEvent() -> Reading="<<ev<<"\n";
  return ev;
}
