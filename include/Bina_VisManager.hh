#ifndef Bina_VisManager_h
#define Bina_VisManager_h 1

#ifdef G4VIS_USE

#include "G4VisManager.hh"


class Bina_VisManager: public G4VisManager {

public:

  Bina_VisManager ();


private:

  void RegisterGraphicsSystems ();

};

#endif

#endif
