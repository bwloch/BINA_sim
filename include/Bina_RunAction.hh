#ifndef Bina_RunAction_h
#define Bina_RunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"
#include "time.h"
#include "g4root.hh"
class Bina_EventAction;
class G4Run;
class Bina_RunAction : public G4UserRunAction
{
  public:
    Bina_RunAction(Bina_EventAction* eventAction);
    virtual ~Bina_RunAction();

    virtual void BeginOfRunAction(const G4Run* );
    virtual void   EndOfRunAction(const G4Run* );

    private:
    clock_t fbegin;
    clock_t fend;
    Bina_EventAction* fEventAction;



};



#endif

