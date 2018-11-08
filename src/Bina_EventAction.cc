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
//Initialization
Part_num=0;
fN=0;
fX1=-999;
fY1=-999;
fTh1=-999;
fPhi1=-999;
fEn1=-999;
fEd1=-999;
fE1=-999;
fdE1=-999;
fP1Type=-999;
fX2=-999;
fY2=-999;
fTh2=-999;
fPhi2=-999;
fEn2=-999;
fEd2=-999;
fE2=-999;
fdE2=-999;
fP2Type=-999;
fX3=-999;
fY3=-999;
fTh3=-999;
fPhi3=-999;
fEn3=-999;
fEd3=-999;
fE3=-999;
fdE3=-999;
fP3Type=-999;
fXv=-999;
fYv=-999;
fZv=-999;
fX1vec.clear();
//G4cout<<"\n\t MyLog: BeginOfEventAction";

 G4int evtNb = evt->GetEventID();
 getNb(evtNb);
 if (!(evtNb%100)) G4cout << "\n--> Begin of event: " << evtNb <<G4endl;
 //G4cout<<"\n\t MyLog: początek Eventu"<<G4endl;
}

void Bina_EventAction::EndOfEventAction(const G4Event* evt)
{
//G4cout<<"\n\t MyLog: koniec Eventu"<<G4endl;
G4int evt_Num = evt->GetEventID();
G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

//Filling Ntuples
analysisManager->FillNtupleIColumn(1,0,evt_Num);
//analysisManager->FillNtupleDColumn(1,1,fX1vec[0]);
analysisManager->FillNtupleDColumn(1,2,fX2);
analysisManager->FillNtupleDColumn(1,3,fX3);
analysisManager->FillNtupleDColumn(1,4,fY1);
analysisManager->FillNtupleDColumn(1,5,fY2);
analysisManager->FillNtupleDColumn(1,6,fY3);
analysisManager->FillNtupleIColumn(1,7,fP1Type);
analysisManager->FillNtupleIColumn(1,8,fP2Type);
analysisManager->FillNtupleIColumn(1,9,fP3Type);
analysisManager->FillNtupleIColumn(1,10,fE1);
analysisManager->FillNtupleIColumn(1,11,fE2);
analysisManager->FillNtupleIColumn(1,12,fE3);
analysisManager->FillNtupleIColumn(1,13,fdE1);
analysisManager->FillNtupleIColumn(1,14,fdE2);
analysisManager->FillNtupleIColumn(1,15,fdE3);
analysisManager->FillNtupleIColumn(1,16,fN);
analysisManager->FillNtupleDColumn(1,17,fEn1);
analysisManager->FillNtupleDColumn(1,18,fEn2);
analysisManager->FillNtupleDColumn(1,19,fEn3);
analysisManager->FillNtupleDColumn(1,20,fEd1);
analysisManager->FillNtupleDColumn(1,21,fEd2);
analysisManager->FillNtupleDColumn(1,22,fEd3);
analysisManager->FillNtupleDColumn(1,23,fTh1);
analysisManager->FillNtupleDColumn(1,24,fTh2);
analysisManager->FillNtupleDColumn(1,25,fTh3);
analysisManager->FillNtupleDColumn(1,26,fPhi1);
analysisManager->FillNtupleDColumn(1,27,fPhi2);
analysisManager->FillNtupleDColumn(1,28,fPhi3);
analysisManager->FillNtupleDColumn(1,29,fXv);
analysisManager->FillNtupleDColumn(1,30,fYv);
analysisManager->FillNtupleDColumn(1,31,fZv);
analysisManager->FillNtupleDColumn(1,32,fPhi1);
analysisManager->FillNtupleDColumn(1,33,fPhi2);
analysisManager->FillNtupleDColumn(1,34,fPhi3);
analysisManager->FillNtupleDColumn(1,35,fTh1);
analysisManager->FillNtupleDColumn(1,36,fTh2);
analysisManager->FillNtupleDColumn(1,37,fTh3);
analysisManager->FillNtupleDColumn(1,38,fEn1);
analysisManager->FillNtupleDColumn(1,39,fEn2);
analysisManager->FillNtupleDColumn(1,40,fEn3);
//analysisManager->FillNtupleDColumn(1,1,fEnergy);
analysisManager->AddNtupleRow(1);
//G4cout<<" \n\t Myog: EventAction En1="<<fEn1<<" \t En2="<<fEn2<<"\t En1+En2="<<fEn1+fEn2<<G4endl;

}

void Bina_EventAction::AddHits(G4int Ptype, G4double X, G4double Y, G4double Th, G4double phi, G4double En, G4double Ed, G4double E, G4double dE, G4double Xv, G4double Yv, G4double Zv){
G4int pos_num=0;

//Checking particle type and its position in FBEvent structure
if(Ptype==2){//Proton
	if(Part_num==0) {
		pos_num=1;
		Part_num++;
		}
	else pos_num=2;
	}
if(Ptype==1) pos_num=3; //Neutron
if(Ptype==3) pos_num=2; //Deuteron

//Filling variables
if(pos_num==1){
	fX1=X;
	fY1=Y;
	fTh1=Th;
	fPhi1=phi;
	fEn1=En;
	fEd1=Ed;
	fE1=E;
	fdE1=dE;
	fP1Type=1;
	fX1vec.push_back(X);
}
if(pos_num==2){
	fX2=X;
	fY2=Y;
	fTh2=Th;
	fPhi2=phi;
	fEn2=En;
	fEd2=Ed;
	fE2=E;
	fdE2=dE;
	fP2Type=2;
}
if(pos_num==3){
	fX3=X;
	fY3=Y;
	fTh3=Th;
	fPhi3=phi;
	fEn3=En;
	fEd3=Ed;
	fE3=E;
	fdE3=dE;
	fP3Type=3;
}
fXv=Xv;
fYv=Yv;
fZv=Zv;
fN++;
//G4cout<<"\n\t Mylog: AddHits En1="<<En<<" \t Ptype="<<Ptype<<"";
}


