#include "Bina_RunAction.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

Bina_RunAction::Bina_RunAction()
{ 
}
Bina_RunAction::~Bina_RunAction(){

}
void Bina_RunAction::BeginOfRunAction(const G4Run* ){
fbegin = clock(); //Start clock
G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
G4cout << "\n\t MyLog: Using " << analysisManager->GetType() << G4endl;
analysisManager->OpenFile("test.root");
//Creating Ntuples
analysisManager->SetFirstNtupleId(1);
analysisManager->CreateNtuple("T","Title");
analysisManager->CreateNtupleIColumn("Event_evNum");
analysisManager->CreateNtupleDColumn("Event_fX1");
analysisManager->CreateNtupleDColumn("Event_fX2");
analysisManager->CreateNtupleDColumn("Event_fX3");
analysisManager->CreateNtupleDColumn("Event_fY1");
analysisManager->CreateNtupleDColumn("Event_fY2");
analysisManager->CreateNtupleDColumn("Event_fY3");
analysisManager->CreateNtupleIColumn("Event_fP1Type");
analysisManager->CreateNtupleIColumn("Event_fP2Type");
analysisManager->CreateNtupleIColumn("Event_fP3Type");
analysisManager->CreateNtupleIColumn("Event_fE1");
analysisManager->CreateNtupleIColumn("Event_fE2");
analysisManager->CreateNtupleIColumn("Event_fE3");
analysisManager->CreateNtupleIColumn("Event_fdE1");
analysisManager->CreateNtupleIColumn("Event_fdE2");
analysisManager->CreateNtupleIColumn("Event_fdE3");
analysisManager->CreateNtupleIColumn("Event_fN");
analysisManager->CreateNtupleDColumn("Event_fEn1");
analysisManager->CreateNtupleDColumn("Event_fEn2");
analysisManager->CreateNtupleDColumn("Event_fEn3");
analysisManager->CreateNtupleDColumn("Event_fEd1");
analysisManager->CreateNtupleDColumn("Event_fEd2");
analysisManager->CreateNtupleDColumn("Event_fEd3");
analysisManager->CreateNtupleDColumn("Event_fTh1");
analysisManager->CreateNtupleDColumn("Event_fTh2");
analysisManager->CreateNtupleDColumn("Event_fTh3");
analysisManager->CreateNtupleDColumn("Event_fPhi1");
analysisManager->CreateNtupleDColumn("Event_fPhi2");
analysisManager->CreateNtupleDColumn("Event_fPhi3");
analysisManager->CreateNtupleDColumn("Event_fXz");
analysisManager->CreateNtupleDColumn("Event_fYz");
analysisManager->CreateNtupleDColumn("Event_fZv");
analysisManager->CreateNtupleDColumn("Event_fVexPhi1");
analysisManager->CreateNtupleDColumn("Event_fVexfPhi2");
analysisManager->CreateNtupleDColumn("Event_fVexfPhi3");
analysisManager->CreateNtupleDColumn("Event_fVexTh1");
analysisManager->CreateNtupleDColumn("Event_fVexTh2");
analysisManager->CreateNtupleDColumn("Event_fVexTh3");
analysisManager->CreateNtupleDColumn("Event_fVexEn1");
analysisManager->CreateNtupleDColumn("Event_fVexEn2");
analysisManager->CreateNtupleDColumn("Event_fVexEn3");

  
analysisManager->FinishNtuple();
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
