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
fbegin = clock();
G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
G4cout << "\n\t MyLog: Using " << analysisManager->GetType() << G4endl;
analysisManager->OpenFile("test.root");
analysisManager->SetFirstNtupleId(1);
analysisManager->CreateNtuple("name","Title");
analysisManager->CreateNtupleIColumn("evNum");
analysisManager->CreateNtupleDColumn("Energy");
analysisManager->FinishNtuple();
}
void Bina_RunAction::EndOfRunAction(const G4Run* ){
G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
analysisManager->Write();
analysisManager->CloseFile();
fend = clock();
double time_spent;
time_spent = (double)(fend - fbegin) / CLOCKS_PER_SEC;
G4cout<<"\n\t MyLog: czas wykonywania programu ="<<time_spent<<G4endl;
}
