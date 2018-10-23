
#ifndef Bina_EventAction_h
#define Bina_EventAction_h 1

#include "G4UserEventAction.hh"
#include "g4root.hh"
#include "globals.hh"


class Bina_EventAction : public G4UserEventAction
{
 private:
  G4double  fX1, fX2, fX3;	         // Coord. on MWPC
  G4double  fY1, fY2, fY3;             // -||-
  G4int    fP1Type,fP2Type,fP3Type;      // Kind of particle
  G4int    fE1, fE2, fE3;             // Which E det
  G4int    fdE1, fdE2, fdE3;           // Which dE det
  G4int    fN;                   // number of particles in event
  G4double fEn1, fEn2, fEn3;            // Energies from vert
  G4double fEd1, fEd2, fEd3;            // Energies dep. in E
  G4double fTh1, fTh2, fTh3;            // Theta angles of particles (vert)
  G4double fPhi1, fPhi2, fPhi3;          // Phi angles of particles (vert)
  G4double  fXv, fYv, fZv;			// Vertex position
  G4int Part_num;
  
  
  public:
    Bina_EventAction();
    virtual ~Bina_EventAction();

  public:
    virtual void   BeginOfEventAction(const G4Event*);
    virtual void   EndOfEventAction(const G4Event*);
    void AddHits(G4int Ptype, G4double X, G4double Y, G4double Th, G4double phi, G4double En, G4double Ed, G4double E, G4double dE, G4double Xv, G4double Yv, G4double Zv);
    inline static int getNb(int num = -1 )
    {
      static int temp;
      if (num != -1) temp = num;
      return temp;
    };

	
};


#endif


