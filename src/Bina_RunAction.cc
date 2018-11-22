#include "Bina_RunAction.hh"
#include "Bina_EventAction.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

Bina_RunAction::Bina_RunAction(Bina_EventAction* eventAction)
: fEventAction(eventAction)
{ 
G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
G4cout << "\n\t MyLog: Using " << analysisManager->GetType() << G4endl;
//Creating Ntuples
analysisManager->SetFileName("test.root");
if( fEventAction){
analysisManager->SetFirstNtupleId(1);
analysisManager->CreateNtuple("T","Title");
analysisManager->CreateNtupleIColumn(1,"Event_evNum");
analysisManager->CreateNtupleDColumn(1,"Event_fX1");
analysisManager->CreateNtupleDColumn(1,"Event_fX2");
analysisManager->CreateNtupleDColumn(1,"Event_fX3");
analysisManager->CreateNtupleDColumn(1,"Event_fY1");
analysisManager->CreateNtupleDColumn(1,"Event_fY2");
analysisManager->CreateNtupleDColumn(1,"Event_fY3");
analysisManager->CreateNtupleIColumn(1,"Event_fP1Type");
analysisManager->CreateNtupleIColumn(1,"Event_fP2Type");
analysisManager->CreateNtupleIColumn(1,"Event_fP3Type");
analysisManager->CreateNtupleIColumn(1,"Event_fE1");
analysisManager->CreateNtupleIColumn(1,"Event_fE2");
analysisManager->CreateNtupleIColumn(1,"Event_fE3");
analysisManager->CreateNtupleIColumn(1,"Event_fdE1");
analysisManager->CreateNtupleIColumn(1,"Event_fdE2");
analysisManager->CreateNtupleIColumn(1,"Event_fdE3");
analysisManager->CreateNtupleIColumn(1,"Event_fN");
analysisManager->CreateNtupleDColumn(1,"Event_fEn1");
analysisManager->CreateNtupleDColumn(1,"Event_fEn2");
analysisManager->CreateNtupleDColumn(1,"Event_fEn3");
analysisManager->CreateNtupleDColumn(1,"Event_fEd1");
analysisManager->CreateNtupleDColumn(1,"Event_fEd2");
analysisManager->CreateNtupleDColumn(1,"Event_fEd3");
analysisManager->CreateNtupleDColumn(1,"Event_fTh1");
analysisManager->CreateNtupleDColumn(1,"Event_fTh2");
analysisManager->CreateNtupleDColumn(1,"Event_fTh3");
analysisManager->CreateNtupleDColumn(1,"Event_fPhi1");
analysisManager->CreateNtupleDColumn(1,"Event_fPhi2");
analysisManager->CreateNtupleDColumn(1,"Event_fPhi3");
analysisManager->CreateNtupleDColumn(1,"Event_fXz");
analysisManager->CreateNtupleDColumn(1,"Event_fYz");
analysisManager->CreateNtupleDColumn(1,"Event_fZv");
analysisManager->CreateNtupleDColumn(1,"Event_fVexPhi1");
analysisManager->CreateNtupleDColumn(1,"Event_fVexfPhi2");
analysisManager->CreateNtupleDColumn(1,"Event_fVexfPhi3");
analysisManager->CreateNtupleDColumn(1,"Event_fVexTh1");
analysisManager->CreateNtupleDColumn(1,"Event_fVexTh2");
analysisManager->CreateNtupleDColumn(1,"Event_fVexTh3");
analysisManager->CreateNtupleDColumn(1,"Event_fVexEn1");
analysisManager->CreateNtupleDColumn(1,"Event_fVexEn2");
analysisManager->CreateNtupleDColumn(1,"Event_fVexEn3");
analysisManager->FinishNtuple();
}
}
Bina_RunAction::~Bina_RunAction(){
 delete G4AnalysisManager::Instance(); 
}
void Bina_RunAction::BeginOfRunAction(const G4Run* ){
fbegin = clock(); //Start clock
G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
analysisManager->OpenFile();
}

void Bina_RunAction::EndOfRunAction(const G4Run* ){
//Writing and closing root file
G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
analysisManager->Write();
analysisManager->CloseFile();
//calculating time
fend = clock();
double time_spent;
time_spent = (double)(fend - fbegin) / CLOCKS_PER_SEC;
G4cout<<"\n\t MyLog: czas wykonywania programu ="<<time_spent<<G4endl;
}
