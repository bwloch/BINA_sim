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
analysisManager->CreateNtupleIColumn(1,"evNum");
analysisManager->CreateNtupleDColumn(1,"fX1");
analysisManager->CreateNtupleDColumn(1,"fX2");
analysisManager->CreateNtupleDColumn(1,"fX3");
analysisManager->CreateNtupleDColumn(1,"fY1");
analysisManager->CreateNtupleDColumn(1,"fY2");
analysisManager->CreateNtupleDColumn(1,"fY3");
analysisManager->CreateNtupleIColumn(1,"fP1Type");
analysisManager->CreateNtupleIColumn(1,"fP2Type");
analysisManager->CreateNtupleIColumn(1,"fP3Type");
analysisManager->CreateNtupleIColumn(1,"fE1");
analysisManager->CreateNtupleIColumn(1,"fE2");
analysisManager->CreateNtupleIColumn(1,"fE3");
analysisManager->CreateNtupleIColumn(1,"fdE1");
analysisManager->CreateNtupleIColumn(1,"fdE2");
analysisManager->CreateNtupleIColumn(1,"fdE3");
analysisManager->CreateNtupleIColumn(1,"fN");
analysisManager->CreateNtupleDColumn(1,"fEn1");
analysisManager->CreateNtupleDColumn(1,"fEn2");
analysisManager->CreateNtupleDColumn(1,"fEn3");
analysisManager->CreateNtupleDColumn(1,"fEd1");
analysisManager->CreateNtupleDColumn(1,"fEd2");
analysisManager->CreateNtupleDColumn(1,"fEd3");
analysisManager->CreateNtupleDColumn(1,"fTh1");
analysisManager->CreateNtupleDColumn(1,"fTh2");
analysisManager->CreateNtupleDColumn(1,"fTh3");
analysisManager->CreateNtupleDColumn(1,"fPhi1");
analysisManager->CreateNtupleDColumn(1,"fPhi2");
analysisManager->CreateNtupleDColumn(1,"fPhi3");
analysisManager->CreateNtupleDColumn(1,"fXz");
analysisManager->CreateNtupleDColumn(1,"fYz");
analysisManager->CreateNtupleDColumn(1,"fZv");
analysisManager->CreateNtupleDColumn(1,"fVexPhi1");
analysisManager->CreateNtupleDColumn(1,"fVexfPhi2");
analysisManager->CreateNtupleDColumn(1,"fVexfPhi3");
analysisManager->CreateNtupleDColumn(1,"fVexTh1");
analysisManager->CreateNtupleDColumn(1,"fVexTh2");
analysisManager->CreateNtupleDColumn(1,"fVexTh3");
analysisManager->CreateNtupleDColumn(1,"fVexEn1");
analysisManager->CreateNtupleDColumn(1,"fVexEn2");
analysisManager->CreateNtupleDColumn(1,"fVexEn3");
analysisManager->CreateNtupleIColumn(1,"fFlagMWPC1");
analysisManager->CreateNtupleIColumn(1,"fFlagMWPC2");
analysisManager->CreateNtupleIColumn(1,"fFlagMWPC3");
analysisManager->CreateNtupleIColumn(1,"fFlagE1");
analysisManager->CreateNtupleIColumn(1,"fFlagE2");
analysisManager->CreateNtupleIColumn(1,"fFlagE3");
analysisManager->CreateNtupleIColumn(1,"fFlagdE1");
analysisManager->CreateNtupleIColumn(1,"fFlagdE2");
analysisManager->CreateNtupleIColumn(1,"fFlagdE3");
analysisManager->CreateNtupleDColumn(1,"fX4");
analysisManager->CreateNtupleDColumn(1,"fY4");
analysisManager->CreateNtupleIColumn(1,"fP4Type");
analysisManager->CreateNtupleIColumn(1,"fE4");
analysisManager->CreateNtupleIColumn(1,"fdE4");
analysisManager->CreateNtupleDColumn(1,"fEn4");
analysisManager->CreateNtupleDColumn(1,"fEd4");
analysisManager->CreateNtupleDColumn(1,"fTh4");
analysisManager->CreateNtupleDColumn(1,"fPhi4");
analysisManager->CreateNtupleIColumn(1,"fFlagMWPC4");
analysisManager->CreateNtupleIColumn(1,"fFlagE4");
analysisManager->CreateNtupleIColumn(1,"fFlagdE4");
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
