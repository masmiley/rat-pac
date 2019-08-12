#include <string>
#include <G4UIcmdWithAString.hh>
#include <G4UIcmdWithABool.hh>
#include <RAT/PhysicsListMessenger.hh>
#include <RAT/Gsim.hh>
#include <RAT/PhysicsList.hh>
#include <RAT/Log.hh>

namespace RAT {

PhysicsListMessenger::PhysicsListMessenger(PhysicsList* physicsList)
    : fPhysicsList(physicsList) {
  fSetOpWLSCmd = new G4UIcmdWithAString("/PhysicsList/setOpWLS", this);
  fSetOpWLSCmd->SetParameterName("model", false);
  fSetOpWLSCmd->SetGuidance("Select a WLS model (g4|bnl)");
  fSetOpWLSCmd->SetDefaultValue("g4");
  fDisableCerenkovCmd = new G4UIcmdWithABool("/PhysicsList/disableCerenkov", this);
  fDisableCerenkovCmd->SetParameterName("disableCerenkov", true, false);
  fDisableCerenkovCmd->SetGuidance("Disable cherenkov photons (true|false");
  fDisableCerenkovCmd->SetDefaultValue(false);
}

PhysicsListMessenger::~PhysicsListMessenger() {
  delete fSetOpWLSCmd;
  delete fDisableCerenkovCmd;
}

G4String PhysicsListMessenger::GetCurrentValue(G4UIcommand* command) {
  if (command == fSetOpWLSCmd) {
    return G4String(fPhysicsList->GetOpWLSModelName().c_str());
  }
  else if (command == fDisableCerenkovCmd) {
    return fPhysicsList->GetDisableCerenkov();
  }
  else {
    Log::Die(
      dformat("PhysicsListMessenger::GetCurrentValue: Unknown command %s",
      command->GetCommandPath().data()));
  }

  return "";
}

void PhysicsListMessenger::SetNewValue(G4UIcommand* command,
                                       G4String newValue) {
  info << "PhysicsListMessenger: Setting WLS model to " << newValue << newline;
  if (command == fSetOpWLSCmd) {
    fPhysicsList->SetOpWLSModel(std::string(newValue.data()));
  }
  else if (command == fDisableCerenkovCmd){
    info << "DISABLING CERENKOV" << newline;
    fPhysicsList->SetDisableCerenkov(fDisableCerenkovCmd->GetNewBoolValue(newValue));
  }
  else {
    Log::Die(
      dformat("PhysicsListMessenger::SetCurrentValue: Unknown command %s",
              command->GetCommandPath().data()));
  }
}

}  // namespace RAT

