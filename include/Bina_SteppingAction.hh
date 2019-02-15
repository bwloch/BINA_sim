#ifndef Bina_SteppingAction_H
#define Bina_SteppingAction_H 1

#include "globals.hh"
#include "G4UserSteppingAction.hh"
class Bina_PhysicsMessenger;
class G4Step;
class G4Track;
class Bina_EventAction;
class Bina_PrimaryGeneratorAction;

class Bina_SteppingAction : public G4UserSteppingAction
{
  public:
    Bina_SteppingAction(Bina_EventAction*, Bina_PrimaryGeneratorAction* );
    virtual ~Bina_SteppingAction();
    virtual void UserSteppingAction(const G4Step*);
    G4Track* fTrack;
    G4Step* fStep;
    int energy_broadening;
    int file_types;
  private:
  double *energy, *theta, *phi ,*position;
  Bina_EventAction* fEventAction;
  Bina_PrimaryGeneratorAction* fPrimaryGeneratorAction;
};

#endif

