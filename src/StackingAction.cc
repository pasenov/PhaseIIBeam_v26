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
/// \file electromagnetic/TestEm18/src/StackingAction.cc
/// \brief Implementation of the StackingAction class
//
// $Id: StackingAction.cc 67268 2013-02-13 11:38:40Z ihrivnac $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "StackingAction.hh"

#include "DetectorConstruction.hh"
#include "EventAction.hh"
#include "HistoManager.hh"
#include "RunAction.hh"

#include "G4Step.hh"
#include "G4Track.hh"
#include "G4VTrajectory.hh"
#include "G4TransportationManager.hh"
#include "G4Navigator.hh"
#include "G4PVPlacement.hh"
#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4Positron.hh"
#include "G4Electron.hh"
#include "G4PionPlus.hh"
#include "G4PionMinus.hh"
#include "G4MuonPlus.hh"
#include "G4MuonMinus.hh"
#include "G4Gamma.hh"

#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4VisAttributes.hh"

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

StackingAction::StackingAction(DetectorConstruction* DA, RunAction* RA, EventAction* EA)
:G4UserStackingAction(), fDetectorconstruction(DA), fRunaction(RA), fEventaction(EA)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

StackingAction::~StackingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ClassificationOfNewTrack
StackingAction::ClassifyNewTrack(const G4Track* track)
{
  //define speed unit
  new G4UnitDefinition ("meter/second", "m/s", "Speed", m/s);
  
  //energy spectrum, polarization, initial momentum direction, track weight, velocity and local time of secondaries and tertiaries
  //
  energy = track->GetKineticEnergy();
  xpolarization = track->GetPolarization().x();
  ypolarization = track->GetPolarization().y();
  zpolarization = track->GetPolarization().z();
  charged = (track->GetDefinition()->GetPDGCharge() != 0.);
  G4VPhysicalVolume* volume = track->GetVolume();
  G4double trackWeight = track->GetWeight();
  G4double initVelocity = track->GetVelocity();
  G4double kinEnergy = track->GetKineticEnergy();
  G4double time = track->GetLocalTime();

  //geometry parameters
  G4double Det1SizeZ, Det2SizeZ, Dist, Strip1Depth, Strip1Length, Strip2Depth, Strip2Length, StripDist, StripWidth, xAbs, XangleDUT, XangleBPIX, YangleBPIX, theta;
  G4int div;
  G4int totalNbStrips, stripNo;
  Det1SizeZ=fDetectorconstruction->GetSize1();
  Det2SizeZ=fDetectorconstruction->GetSize2();
  Dist=fDetectorconstruction->GetDist();
  Strip1Depth=fDetectorconstruction->GetStrip1Depth();
  Strip1Length=fDetectorconstruction->GetStrip1Length();
  Strip2Depth=fDetectorconstruction->GetStrip2Depth();
  Strip2Length=fDetectorconstruction->GetStrip2Length();
  StripDist=fDetectorconstruction->GetStripDist();
  StripWidth=fDetectorconstruction->GetStripWidth();
  XangleDUT=fDetectorconstruction->GetDUTangleX();
  XangleBPIX=fDetectorconstruction->GetBPIXangleX();
  YangleBPIX=fDetectorconstruction->GetBPIXangleY();
  totalNbStrips = 1016;

  //G4ThreeVector primVert = track->GetVertexPosition();
  G4double x_prim, y_prim, z_prim, xPix1, yPix1, zPix1, xPix, yPix, zPix, x_after, y_after, z_after;
  x_after = track->GetVertexPosition().x();
  y_after = track->GetVertexPosition().y();
  z_after = track->GetVertexPosition().z();
  //position if the DUT wasn't rotated
  x_prim = x_after;
  y_prim = y_after*cos(XangleDUT) - z_after*sin(XangleDUT);
  z_prim = y_after*sin(XangleDUT) + z_after*cos(XangleDUT);
  //position if the BPIX modules weren't rotated
  xPix1 = x_after;
  yPix1 = y_after*cos(XangleBPIX) - z_after*sin(XangleBPIX);
  zPix1 = y_after*sin(XangleBPIX) + z_after*cos(XangleBPIX);
  xPix = xPix1*cos(YangleBPIX) + zPix1*sin(YangleBPIX);
  yPix = yPix1;
  zPix = -xPix1*sin(YangleBPIX) + zPix1*cos(YangleBPIX);


  //secret number
  G4int iSecret, jSecret;

  //G4cout 
      //<< "\n Local time: " 
      //<< G4BestUnit(time, "Time")
      //<< "\n Initial velocity: " 
      //<< G4BestUnit(initVelocity, "Speed")
      //<< "\n Kinetic energy: "
      //<< G4BestUnit(kinEnergy, "Energy")
      //<< "\n Track weight: " 
      //<< trackWeight
      //<< G4endl;

  //particle definition
  const G4ParticleDefinition* particleDefinition = track->GetParticleDefinition();

  G4ThreeVector initMom = track->GetMomentum();

  //polarization of secondaries and tertiaries calculation
  if (track->GetParentID() == 1)  {
    fEventaction->AddSecondaryxPolarization(xpolarization);
    fEventaction->AddSecondaryyPolarization(ypolarization ); 
    fEventaction->AddSecondaryzPolarization(zpolarization); 
  }

  if (track->GetParentID() > 1)  {
    fEventaction->AddTertiaryxPolarization(xpolarization); 
    fEventaction->AddTertiaryyPolarization(ypolarization); 
    fEventaction->AddTertiaryzPolarization(zpolarization); 
  }

  //kinetic energy of secondaries calculation 
  if (track->GetParentID() == 1)  {

    if (volume == fDetectorconstruction->GetWorld()) {

      if (z_prim < -Dist/2 - Det1SizeZ)  {

      	if(particleDefinition == G4Positron::Definition())   {
          fRunaction->AddNbPositronsB1();
      	}

      	if(particleDefinition == G4Electron::Definition())   {
          fRunaction->AddNbElectronsB1();
      	}

      }
      if ((z_prim > -Dist/2) && (z_prim < Dist/2))  {

      	if(particleDefinition == G4Positron::Definition())   {
          fRunaction->AddNbPositrons12();
      	}

      	if(particleDefinition == G4Electron::Definition())   {
          fRunaction->AddNbElectrons12();
      	}

      }
      if (z_prim > Dist/2 + Det2SizeZ)  {

      	if(particleDefinition == G4Positron::Definition())   {
          fRunaction->AddNbPositronsA2();
      	}

      	if(particleDefinition == G4Electron::Definition())   {
          fRunaction->AddNbElectronsA2();
      	}

      }

      //G4cout 
          //<< "\n Secondary produced with initial momentum direction:" 
          //<< initMom
          //<< G4endl;
    }

    if (volume == fDetectorconstruction->GetDet1()) {
      G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  
      if (charged) {
       fRunaction->AddChargedSecondary1(energy);
       analysisManager->FillH1(9,energy);
      } else {
       fRunaction->AddNeutralSecondary1(energy);
       analysisManager->FillH1(11,energy);
      }

      fEventaction->AddSecondary1(energy); 

      if(particleDefinition == G4Positron::Definition())   {
        fRunaction->AddNbPositrons1();
      }

      if(particleDefinition == G4Electron::Definition())   {
        fRunaction->AddNbElectrons1();
	if ((z_prim >= -Dist/2 - Det1SizeZ) && (z_prim <= -Dist/2 - Det1SizeZ + Strip1Depth) && (x_prim >= (-Strip1Length/2)) && (x_prim <= Strip1Length/2)) {
            div = y_prim/(StripWidth+StripDist);
	    if (y_prim == 0)   {
	       iSecret = rand() % 99;
	       if (iSecret < 50)   {
		  stripNo = totalNbStrips/2 + div;
	       }
	       if (iSecret >= 50)   {
	       	  stripNo = totalNbStrips/2 + 1 + div;
	       }
	    }	       
            if (y_prim > 0)    {
	       stripNo = totalNbStrips/2 + 1 + div;
            }
            if (y_prim < 0)    {
               stripNo = totalNbStrips/2 + div;
            }
            if ((stripNo <= totalNbStrips) && (stripNo > 0))  {
     	    //G4cout
	       //<< "\n Secondary electron created below strip No "
               //<< stripNo
               //<< "\n y = "
	       //<< y_prim
               //<< "\n div = "
	       //<< div
               //<< G4endl;
               if (x_prim > 0)  {
	          fRunaction->AddNbElectronsStrip1a(stripNo);
               }
               if (x_prim < 0)  {
	          fRunaction->AddNbElectronsStrip1b(stripNo);
               }
               if (x_prim == 0)  {
	       	jSecret = rand() % 99;
	       	if (jSecret < 50)   {
	            fRunaction->AddNbElectronsStrip1a(stripNo);
	       	}
	       	if (jSecret >= 50)   {
	            fRunaction->AddNbElectronsStrip1b(stripNo);
	       	}
               }
            }
        }
      }

      if(particleDefinition == G4PionPlus::Definition())   {
        fRunaction->AddNbPionsPlus1();
      }

      if(particleDefinition == G4PionMinus::Definition())   {
        fRunaction->AddNbPionsMinus1();
      }

      if(particleDefinition == G4MuonPlus::Definition())   {
        fRunaction->AddNbMuonsPlus1();
      }

      if(particleDefinition == G4MuonMinus::Definition())   {
        fRunaction->AddNbMuonsMinus1();
      }

      //G4cout 
          //<< "\n Secondary produced with initial momentum direction:" 
          //<< initMom
          //<< G4endl;
    }

    if (volume == fDetectorconstruction->GetDet2()) {
      G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  
      if (charged) {
       fRunaction->AddChargedSecondary2(energy);
       analysisManager->FillH1(10,energy);
      } else {
       fRunaction->AddNeutralSecondary2(energy);
       analysisManager->FillH1(12,energy);
      }

      fEventaction->AddSecondary2(energy); 

      if(particleDefinition == G4Positron::Definition())   {
        fRunaction->AddNbPositrons2();
      }

      if(particleDefinition == G4Electron::Definition())   {
        fRunaction->AddNbElectrons2();
	if ((z_prim >= Dist/2) && (z_prim <= (Dist/2 + Strip2Depth))&& (x_prim >= (-Strip2Length/2)) && (x_prim <= Strip2Length/2)) {
            div = y_prim/(StripWidth+StripDist);
	    if (y_prim == 0)   {
	       iSecret = rand() % 99;
	       if (iSecret < 50)   {
		  stripNo = totalNbStrips/2 + div;
	       }
	       if (iSecret >= 50)   {
	       	  stripNo = totalNbStrips/2 + 1 + div;
	       }
	    }
            if (y_prim > 0)    {
	       stripNo = totalNbStrips/2 + 1 + div;
            }
            if (y_prim < 0)    {
               stripNo = totalNbStrips/2 + div;
            }
            if ((stripNo <= totalNbStrips) && (stripNo > 0))  {
     	    //G4cout
	       //<< "\n Secondary electron created below strip No "
               //<< stripNo
               //<< "\n y = "
	       //<< y_prim
               //<< "\n div = "
	       //<< div
               //<< G4endl;
               if (x_prim > 0)  {
	          fRunaction->AddNbElectronsStrip2a(stripNo);
               }
               if (x_prim < 0)  {
	          fRunaction->AddNbElectronsStrip2b(stripNo);
               }
               if (x_prim == 0)  {
	       	jSecret = rand() % 99;
	       	if (jSecret < 50)   {
	            fRunaction->AddNbElectronsStrip2a(stripNo);
	       	}
	       	if (jSecret >= 50)   {
	            fRunaction->AddNbElectronsStrip2b(stripNo);
	       	}
               }
            }
        }
      }

      if(particleDefinition == G4PionPlus::Definition())   {
        fRunaction->AddNbPionsPlus2();
      }

      if(particleDefinition == G4PionMinus::Definition())   {
        fRunaction->AddNbPionsMinus2();
      }

      if(particleDefinition == G4MuonPlus::Definition())   {
        fRunaction->AddNbMuonsPlus2();
      }

      if(particleDefinition == G4MuonMinus::Definition())   {
        fRunaction->AddNbMuonsMinus2();
      }

      //G4cout 
          //<< "\n Secondary produced with initial momentum direction:" 
          //<< initMom
          //<< G4endl;
    }

  }

  //kinetic energy of tertiaries calculation 
  if (track->GetParentID() > 1)  {
    if (volume == fDetectorconstruction->GetDet1()) {
      G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  
      if (charged) {
       fRunaction->AddChargedTertiary1(energy);
       analysisManager->FillH1(13,energy);
      } else {
       fRunaction->AddNeutralTertiary1(energy);
       analysisManager->FillH1(15,energy);
      }

      fEventaction->AddTertiary1(energy);

      if(particleDefinition == G4Positron::Definition())   {
        fRunaction->AddNbPositrons1();
      }

      if(particleDefinition == G4Electron::Definition())   {
        fRunaction->AddNbElectrons1();
      }

      if(particleDefinition == G4PionPlus::Definition())   {
        fRunaction->AddNbPionsPlus1();
      }

      if(particleDefinition == G4PionMinus::Definition())   {
        fRunaction->AddNbPionsMinus1();
      }

      if(particleDefinition == G4MuonPlus::Definition())   {
        fRunaction->AddNbMuonsPlus1();
      }

      if(particleDefinition == G4MuonMinus::Definition())   {
        fRunaction->AddNbMuonsMinus1();
      } 
    }

    if (volume == fDetectorconstruction->GetDet2()) {
      G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  
      if (charged) {
       fRunaction->AddChargedTertiary2(energy);
       analysisManager->FillH1(14,energy);
      } else {
       fRunaction->AddNeutralTertiary2(energy);
       analysisManager->FillH1(16,energy);
      }

      fEventaction->AddTertiary2(energy); 

      if(particleDefinition == G4Positron::Definition())   {
        fRunaction->AddNbPositrons2();
      }

      if(particleDefinition == G4Electron::Definition())   {
        fRunaction->AddNbElectrons2();
      }

      if(particleDefinition == G4PionPlus::Definition())   {
        fRunaction->AddNbPionsPlus2();
      }

      if(particleDefinition == G4PionMinus::Definition())   {
        fRunaction->AddNbPionsMinus2();
      }

      if(particleDefinition == G4MuonPlus::Definition())   {
        fRunaction->AddNbMuonsPlus2();
      }

      if(particleDefinition == G4MuonMinus::Definition())   {
        fRunaction->AddNbMuonsMinus2();
      }
    }
  }

  return fUrgent;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
