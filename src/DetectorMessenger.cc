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
/// \file electromagnetic/TestEm18/src/DetectorMessenger.cc
/// \brief Implementation of the DetectorMessenger class
//
// $Id: DetectorMessenger.cc 67268 2013-02-13 11:38:40Z ihrivnac $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorMessenger.hh"

#include "DetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::DetectorMessenger(DetectorConstruction * Det)
:G4UImessenger(),fDetector(Det),
 fTestemDir(0),
 fDetDir(0),    
 fMaterWCmd(0),
 fMaterDCmd(0),
 fSizeZWCmd(0),
 fSizeZ1Cmd(0),
 fSizeZ2Cmd(0),
 fUpdateCmd(0)
{ 
  fTestemDir = new G4UIdirectory("/testem/");
  fTestemDir->SetGuidance("commands specific to this example");
  
  fDetDir = new G4UIdirectory("/testem/det/");
  fDetDir->SetGuidance("detector construction");
        
  fMaterWCmd = new G4UIcmdWithAString("/testem/det/setMatW",this);
  fMaterWCmd->SetGuidance("Select material of the world.");
  fMaterWCmd->SetParameterName("choice",false);
  fMaterWCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fMaterDCmd = new G4UIcmdWithAString("/testem/det/setMatD",this);
  fMaterDCmd->SetGuidance("Select material of the detectors.");
  fMaterDCmd->SetParameterName("choice",false);
  fMaterDCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  fSizeZWCmd = new G4UIcmdWithADoubleAndUnit("/testem/det/setSizeZW",this);
  fSizeZWCmd->SetGuidance("Set size-Z of the world");
  fSizeZWCmd->SetParameterName("SizeZW",false);
  fSizeZWCmd->SetRange("SizeZW>0.");
  fSizeZWCmd->SetUnitCategory("Length");
  fSizeZWCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fSizeZ1Cmd = new G4UIcmdWithADoubleAndUnit("/testem/det/setSizeZ1",this);
  fSizeZ1Cmd->SetGuidance("Set size-Z of the detector 1");
  fSizeZ1Cmd->SetParameterName("SizeZ1",false);
  fSizeZ1Cmd->SetRange("SizeZ1>0.");
  fSizeZ1Cmd->SetUnitCategory("Length");
  fSizeZ1Cmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fSizeZ2Cmd = new G4UIcmdWithADoubleAndUnit("/testem/det/setSizeZ2",this);
  fSizeZ2Cmd->SetGuidance("Set size-Z of the detector 2");
  fSizeZ2Cmd->SetParameterName("SizeZ2",false);
  fSizeZ2Cmd->SetRange("SizeZ2>0.");
  fSizeZ2Cmd->SetUnitCategory("Length");
  fSizeZ2Cmd->AvailableForStates(G4State_PreInit,G4State_Idle);
    
  fUpdateCmd = new G4UIcmdWithoutParameter("/testem/det/update",this);
  fUpdateCmd->SetGuidance("Update calorimeter geometry.");
  fUpdateCmd->SetGuidance("This command MUST be applied before \"beamOn\" ");
  fUpdateCmd->SetGuidance("if you changed geometrical value(s).");
  fUpdateCmd->AvailableForStates(G4State_Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::~DetectorMessenger()
{
  delete fMaterWCmd;
  delete fMaterDCmd;
  delete fSizeZWCmd;
  delete fSizeZ1Cmd;
  delete fSizeZ2Cmd; 
  delete fUpdateCmd;
  delete fDetDir;  
  delete fTestemDir;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{ 
  if( command == fMaterWCmd )
   { fDetector->SetMaterialW(newValue);}

  if( command == fMaterDCmd )
   { fDetector->SetMaterialD(newValue);}
   
  if( command == fSizeZWCmd )
   { fDetector->SetSizeW(fSizeZWCmd->GetNewDoubleValue(newValue));}

  if( command == fSizeZ1Cmd )
   { fDetector->SetSize1(fSizeZ1Cmd->GetNewDoubleValue(newValue));}

  if( command == fSizeZ2Cmd )
   { fDetector->SetSize2(fSizeZ2Cmd->GetNewDoubleValue(newValue));}
     
  if( command == fUpdateCmd )
   { fDetector->UpdateGeometry(); }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
