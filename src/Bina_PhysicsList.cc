//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
/// \file hadronic/Hadr03/src/PhysicsList.cc
/// \brief Implementation of the PhysicsList class
//
// $Id: PhysicsList.cc 70268 2013-05-28 14:17:50Z maire $

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "Bina_PhysicsList.hh"

#include "G4SystemOfUnits.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"

#include "G4HadronElasticPhysics.hh"
#include "G4HadronElasticPhysicsHP.hh"
#include "G4HadronPhysicsQGSP_BERT_HP.hh"
#include "G4HadronPhysicsQGSP_BERT.hh"
//#include "G4HadronPhysicsQGSP_BIC_HP.hh"
#include "G4HadronInelasticQBBC.hh"
#include "G4HadronPhysicsINCLXX.hh"
#include "G4IonPhysics.hh"
#include "G4IonINCLXXPhysics.hh"
#include "G4EmExtraPhysics.hh"

#include "G4EmStandardPhysics_option4.hh"
#include "G4DecayPhysics.hh"
#include "G4RadioactiveDecayPhysics.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
Bina_PhysicsList::Bina_PhysicsList()
:G4VModularPhysicsList()
{
  physicsMessenger = new Bina_PhysicsMessenger(this);
  //add new units
  //
  new G4UnitDefinition( "millielectronVolt", "meV", "Energy", 1.e-3*eV);   
  new G4UnitDefinition( "mm2/g",  "mm2/g", "Surface/Mass", mm2/g);
  new G4UnitDefinition( "um2/mg", "um2/mg","Surface/Mass", um*um/mg);
    
    G4cout<<"Bina_PhysicsList::file_outputs="<<file_outputs<<"===========================================-==-==-=-=-==-\n";


/*
  G4int verb = 1;
  SetVerboseLevel(verb);

  // Hadron Elastic scattering

  
  // Hadron Inelastic Physics
  ////RegisterPhysics( new G4HadronPhysicsFTFP_BERT_HP(verb));
//  RegisterPhysics( new G4HadronPhysicsQGSP_BIC_HP(verb));

  if (GetNeutronModel()==0) 
    {
    RegisterPhysics( new G4HadronElasticPhysicsHP(verb) );
    RegisterPhysics( new G4HadronPhysicsQGSP_BERT_HP(verb));  
    }
  else if (GetNeutronModel()==1) 
    {
    RegisterPhysics( new G4HadronElasticPhysics(verb) );
    RegisterPhysics( new G4HadronPhysicsQGSP_BERT(verb));  
    }*/
  ////RegisterPhysics( new G4HadronInelasticQBBC(verb));        
  ////RegisterPhysics( new G4HadronPhysicsINCLXX(verb));
  
  // Ion Physics
//  RegisterPhysics( new G4IonPhysics(verb));
  ////RegisterPhysics( new G4IonINCLXXPhysics(verb));
    
  // Gamma-Nuclear Physics
  RegisterPhysics( new G4EmExtraPhysics());
  
  // EM physics
  RegisterPhysics(new G4EmStandardPhysics_option4());
  
  // Decay
  RegisterPhysics(new G4DecayPhysics());

  // Radioactive decay
//  RegisterPhysics(new G4RadioactiveDecayPhysics());

}


void Bina_PhysicsList::RegisterHadrons(G4String option) {
  G4int verb = 2;
  SetVerboseLevel(verb);
  // Hadron Elastic scattering

  if (option=='0')   {
  // Hadron Inelastic Physics
  ////RegisterPhysics( new G4HadronPhysicsFTFP_BERT_HP(verb));
//  RegisterPhysics( new G4HadronPhysicsQGSP_BIC_HP(verb));
    RegisterPhysics( new G4HadronElasticPhysicsHP(verb) );
    RegisterPhysics( new G4HadronPhysicsQGSP_BERT_HP(verb));  
    }

  else if (option=='1') 
    {
    RegisterPhysics( new G4HadronElasticPhysics(verb) );
    RegisterPhysics( new G4HadronPhysicsQGSP_BERT(verb));  
    }
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Bina_PhysicsList::~Bina_PhysicsList()
{
delete physicsMessenger;
 }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Bina_PhysicsList::SetCuts()
{
//SetCutsWithDefault();
  SetCutValue(0.01*mm, "proton");
  SetCutValue(0.01*mm, "e-");
  SetCutValue(0.01*mm, "e+");
  SetCutValue(0.01*mm, "gamma");      

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/*
void Bina_PhysicsList::SetParamUpdate() 
  {G4RunManager::GetRunManager()->SetUserAction(this);
  }
*/
