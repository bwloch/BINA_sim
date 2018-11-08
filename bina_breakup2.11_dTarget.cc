#include "Bina_DetectorConstruction.hh"
#include "Bina_PhysicsList.hh"
#include "Bina_PrimaryGeneratorAction.hh"
#include "Bina_EventAction.hh"
#include "Bina_SteppingAction.hh"
#include "Bina_RunAction.hh"

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"

#include "Randomize.hh"

#ifdef G4VIS_USE
#include "Bina_VisManager.hh"
#endif
#include <sys/stat.h>





int main(int argc,char** argv) 
{



 Bina_DetectorConstruction* Bina_detector = new Bina_DetectorConstruction;
 Bina_PhysicsList* Bina_physics=new Bina_PhysicsList;
 G4RunManager * runManager = new G4RunManager;

 runManager->SetUserInitialization(Bina_detector);
 runManager->SetUserInitialization(Bina_physics);

 if(argc<4)
 {
  G4String command = "/control/execute ";
  G4String fileInp;
  G4UImanager * UI = G4UImanager::GetUIpointer();
  
  (argc>2) ? fileInp = argv[2]:fileInp = "geo.mac";
  {
    struct stat fileData;
    if(stat(fileInp, &fileData)==-1)
      perror (fileInp), exit(1);
    G4cout<<"\n\t Dodano domyslny plik z geometria GEO.MAC\n";
    UI->ApplyCommand(command+fileInp);
  }

#ifdef G4VIS_USE
  G4VisManager* visManager = new Bina_VisManager;
  visManager->Initialize();
#endif
Bina_EventAction* EvnAct = new Bina_EventAction;

  // UserAction classes
  runManager->SetUserAction(new Bina_PrimaryGeneratorAction(Bina_detector));
  runManager->SetUserAction(EvnAct);
  runManager->SetUserAction(new Bina_SteppingAction(Bina_physics,EvnAct));
runManager->SetUserAction(new Bina_RunAction(EvnAct));

  //Initialize G4 kernel
  runManager->Initialize();
  G4UIsession* session = 0;
  session = new G4UIterminal(new G4UItcsh);
  

  if(argc<2)
  {
    fileInp = "run.mac";
    struct stat fileData;
    if(stat(fileInp, &fileData)!=-1)
    {G4cout<<"\n\t Dodano domyslny plik z instrukajami RUN.MAC\n";
     UI->ApplyCommand(command+fileInp);
    }
  }

  
  if(argc==3||argc==2) 
  {
    fileInp = argv[1];
    UI->ApplyCommand(command+fileInp);
  }
    session->SessionStart();
    delete session;

#ifdef G4VIS_USE
    delete visManager;
#endif
  delete runManager;


  return 0;
 }
 else
 {
  G4cout<<"Bad argument list"<<G4endl;
  exit(1);
  return 1;
 }
}
