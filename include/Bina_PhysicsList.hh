#ifndef Bina_PhysicsList_h
#define Bina_PhysicsList_h 1

#include "G4VUserPhysicsList.hh"
#include "globals.hh"
#include "G4VModularPhysicsList.hh"
#include "Bina_PhysicsMessenger.hh"
class G4VPhysicsConstructor;
class Bina_PhysicsMessenger;

class Bina_PhysicsList: public G4VModularPhysicsList, public G4UImessenger
{
  public:
    Bina_PhysicsList();
   ~Bina_PhysicsList();
    void SetNeutronModel(G4int neutrony)       {neutrons_model=neutrony;};
//    void SetNeutronType(G4int type)              {neutrons_type=type;};
    void SetBroadening(G4int broad)              {broadening=broad;};
    void SetFileOutputs(G4int fileOut)           {file_outputs=fileOut;};
    int  GetNeutronModel() {return neutrons_model;};
//    int  GetNeutronType()    {return neutrons_type;};
    int  GetBroadening()     {return broadening;};
    double GetKinematicsMin(){return kinematics_min;};
    double GetKinematicsMax(){return kinematics_max;};
    int  GetFileOutputs()    {return file_outputs;};

    void SetParamUpdate();
    void RegisterHadrons (G4String option);    
  private:
    G4int neutrons_model;
//    G4int neutrons_type;
    G4int broadening;
    G4double kinematics_min;
    G4double kinematics_max;
    G4int file_outputs;
    G4String particleName;
    Bina_PhysicsMessenger* physicsMessenger;
    G4VPhysicsConstructor*  emPhysicsList;
 Bina_PhysicsList* phyList;
    G4VPhysicsConstructor*  particleList;

  protected:
    // Construct particle and physics
    void SetCuts();


  protected:
    // these methods Construct particles 
    void ConstructBosons();
    void ConstructLeptons();
    void ConstructMesons();
    void ConstructBaryons();

  protected:
  // these methods Construct physics processes and register them
    void ConstructGeneral();
    void ConstructEM();
    void ConstructOp();
};

#endif
