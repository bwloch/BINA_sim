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
analysisManager->CreateNtuple("Bina_sim","Title");
analysisManager->CreateNtupleIColumn("evNum");
analysisManager->CreateNtupleDColumn("fX1");
analysisManager->CreateNtupleDColumn("fX2");
analysisManager->CreateNtupleDColumn("fX3");
analysisManager->CreateNtupleDColumn("fY1");
analysisManager->CreateNtupleDColumn("fY2");
analysisManager->CreateNtupleDColumn("fY3");
analysisManager->CreateNtupleIColumn("fP1Type");
analysisManager->CreateNtupleIColumn("fP2Type");
analysisManager->CreateNtupleIColumn("fP3Type");
analysisManager->CreateNtupleIColumn("fE1");
analysisManager->CreateNtupleIColumn("fE2");
analysisManager->CreateNtupleIColumn("fE3");
analysisManager->CreateNtupleIColumn("fdE1");
analysisManager->CreateNtupleIColumn("fdE2");
analysisManager->CreateNtupleIColumn("fdE3");
analysisManager->CreateNtupleIColumn("fN");
analysisManager->CreateNtupleDColumn("fEn1");
analysisManager->CreateNtupleDColumn("fEn2");
analysisManager->CreateNtupleDColumn("fEn3");
analysisManager->CreateNtupleDColumn("fEd1");
analysisManager->CreateNtupleDColumn("fEd2");
analysisManager->CreateNtupleDColumn("fEd3");
analysisManager->CreateNtupleDColumn("fTh1");
analysisManager->CreateNtupleDColumn("fTh2");
analysisManager->CreateNtupleDColumn("fTh3");
analysisManager->CreateNtupleDColumn("fPhi1");
analysisManager->CreateNtupleDColumn("fPhi2");
analysisManager->CreateNtupleDColumn("fPhi3");
analysisManager->CreateNtupleDColumn("fXz");
analysisManager->CreateNtupleDColumn("fYz");
analysisManager->CreateNtupleDColumn("fZv");
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
