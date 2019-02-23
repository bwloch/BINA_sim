#include "Bina_ActionInitialization.hh"
#include "Bina_PrimaryGeneratorAction.hh"
#include "Bina_RunAction.hh"
#include "Bina_SteppingAction.hh"
#include "Bina_EventAction.hh"


Bina_ActionInitialization::Bina_ActionInitialization()
 : G4VUserActionInitialization()
{

}

Bina_ActionInitialization::~Bina_ActionInitialization()
{}

void Bina_ActionInitialization::Build() const
{
	Bina_EventAction* EvnAct = new Bina_EventAction();
	Bina_RunAction* runAction= new Bina_RunAction();
	Bina_PrimaryGeneratorAction* GenAct = new Bina_PrimaryGeneratorAction();
	SetUserAction(GenAct);
	SetUserAction(EvnAct);
	SetUserAction(new Bina_SteppingAction(EvnAct,GenAct));
	SetUserAction( runAction);



}

void Bina_ActionInitialization::BuildForMaster() const
{
	Bina_RunAction* runAction= new Bina_RunAction;
	SetUserAction(runAction);
}
