#ifndef Bina_ActionInitialization_h
#define Bina_ActionInitialization_h 1

#include "G4VUserActionInitialization.hh"

class Bina_ActionInitialization : public G4VUserActionInitialization 
{
	public:
		Bina_ActionInitialization();
		virtual ~Bina_ActionInitialization();
		virtual void Build() const;
		virtual void BuildForMaster() const;
	
	private:
	
	
};
#endif
