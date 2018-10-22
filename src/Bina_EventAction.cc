#include "Bina_EventAction.hh"

//#include "Bina_CalorHit.hh"
//#include "Bina_EventActionMessenger.hh"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4VHitsCollection.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4VVisManager.hh"
#include "G4SDManager.hh"
#include "G4UImanager.hh"
#include "G4ios.hh"
#include "G4UnitsTable.hh"


Bina_EventAction::Bina_EventAction()
{;}

Bina_EventAction::~Bina_EventAction()
{;}

void Bina_EventAction::BeginOfEventAction(const G4Event* evt)
{
 G4int evtNb = evt->GetEventID();
 getNb(evtNb);
 if (!(evtNb%100)) G4cout << "\n--> Begin of event: " << evtNb <<G4endl;
 fEnergy=0;
}

void Bina_EventAction::EndOfEventAction(const G4Event* evt)
{
G4int evt_Num = evt->GetEventID();
G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
analysisManager->FillNtupleIColumn(1,0,evt_Num);
analysisManager->FillNtupleDColumn(1,1,fEnergy);
analysisManager->AddNtupleRow(1);


}
void Bina_EventAction::AddEnergy(G4double En){
	fEnergy+=En;
}


