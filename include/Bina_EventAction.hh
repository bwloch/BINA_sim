
#ifndef Bina_EventAction_h
#define Bina_EventAction_h 1

#include "G4UserEventAction.hh"
#include "g4root.hh"
#include "globals.hh"


class Bina_EventAction : public G4UserEventAction
{
  public:
    Bina_EventAction();
    virtual ~Bina_EventAction();

  public:
    virtual void   BeginOfEventAction(const G4Event*);
    virtual void   EndOfEventAction(const G4Event*);
    void AddEnergy(G4double En);
    inline static int getNb(int num = -1 )
    {
      static int temp;
      if (num != -1) temp = num;
      return temp;
    };

    G4double fEnergy;
};


#endif


