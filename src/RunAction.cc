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
/// \file electromagnetic/TestEm18/src/RunAction.cc
/// \brief Implementation of the RunAction class
//
// $Id: RunAction.cc 82401 2014-06-18 14:43:54Z gcosmo $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "RunAction.hh"
#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"
#include "HistoManager.hh"
#include "EventAction.hh"

#include "G4Event.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4VisAttributes.hh"
#include "G4UImanager.hh"
#include "G4EmCalculator.hh"
#include "G4Step.hh"

#include "Randomize.hh"
#include <iomanip>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction(DetectorConstruction* det, PrimaryGeneratorAction* kin)
:G4UserRunAction(),fDetector(det), fPrimary(kin), fHistoManager(0)
{ 
  fHistoManager = new HistoManager(); 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{ 
  delete fHistoManager; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run*)
{
  //initialisation
  //
  fEnergyDeposit1 = 0.;
  fEnergyDeposit2 = 0.;
  
  fNbCharged1 = fNbCharged2 = fNbNeutral1 = fNbNeutral2 = 0;
  fNbPositronsB1 = fNbPositrons1 = fNbPositrons12 = fNbPositrons2 = fNbPositronsA2 = fNbElectronsB1 = fNbElectrons1 = fNbElectrons12 = fNbElectrons2 = fNbElectronsA2 = fNbPionsPlus1 = fNbPionsPlus2 = fNbPionsMinus1 = fNbPionsMinus2 = fNbMuonsPlus1 = fNbMuonsPlus2 = fNbMuonsMinus1 = fNbMuonsMinus2 = 0;
  fNbDelta12 = fNbDelta21 = 0;
  fEnergyCharged1 = fEnergyCharged2 = fEnergyNeutral1 = fEnergyNeutral2 = 0.;  
  fEmin1[0] = fEmin1[1] = DBL_MAX;
  fEmax1[0] = fEmax1[1] = 0.;    
  fEmin2[0] = fEmin2[1] = DBL_MAX;
  fEmax2[0] = fEmax2[1] = 0.;   

  for (G4int k = 0; k < 1017; k = k + 1) {
    fNbElectronsStrip1a[k] = 0;
    fNbElectronsStrip1b[k] = 0;
    fNbElectronsStrip2a[k] = 0;
    fNbElectronsStrip2b[k] = 0;

    fNbHitsStrip1a[k] = 0;
    fNbHitsStrip1b[k] = 0;
    fNbHitsStrip2a[k] = 0;
    fNbHitsStrip2b[k] = 0;
  }

  for (G4int ke1 = 0; ke1 < 3; ke1++) {
     for (G4int ke2 = 0; ke2 < 17; ke2++) {
  	 for (G4int ke3 = 0; ke3 < 81; ke3++) {
  	    for (G4int ke4 = 0; ke4 < 53; ke4++) {
		fNbHitsPixel[ke1][ke2][ke3][ke4] = 0;
 	    }
 	 }
     }
  }
    
  fNbSteps = 0;
  fTrackLength = fSecTrackLength = fTertTrackLength = fSecDetTrackLength = 0.;
   
  //histograms
  //
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  if ( analysisManager->IsActive() ) {
    analysisManager->OpenFile();
  }       

  // do not save Rndm status
  G4RunManager::GetRunManager()->SetRandomNumberStore(false);
  CLHEP::HepRandom::showEngineStatus();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run* aRun)
{
  G4int nbEvents = aRun->GetNumberOfEvent();
  if (nbEvents == 0) return;
  
  G4Material* material = fDetector->GetMaterialD();
  G4double length  = fDetector->GetSizeW();
   
  G4ParticleDefinition* particle = fPrimary->GetParticleGun()
                                          ->GetParticleDefinition();
  G4String partName = particle->GetParticleName();
  G4double ePrimary = fPrimary->GetParticleGun()->GetParticleEnergy();
  
  G4int prec = G4cout.precision(3);
  G4cout << "\n ======================== run summary ======================\n";
  G4cout << "\n The run was " << nbEvents << " " << partName << " of "
         << G4BestUnit(ePrimary,"Energy") << " through " 
         << G4BestUnit(length,"Length") << ".";
  G4cout << "\n ===========================================================\n";
  G4cout << G4endl;
  
  //save histograms      
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->FillH1(32, fNbPositrons1);
  analysisManager->FillH1(33, fNbPositrons2);
  analysisManager->FillH1(34, fNbElectrons1);
  analysisManager->FillH1(35, fNbElectrons2);
  analysisManager->FillH1(36, fNbPionsPlus1);
  analysisManager->FillH1(37, fNbPionsPlus2);
  analysisManager->FillH1(38, fNbPionsMinus1);
  analysisManager->FillH1(39, fNbPionsMinus2);
  analysisManager->FillH1(40, fNbMuonsPlus1);
  analysisManager->FillH1(41, fNbMuonsPlus2);
  analysisManager->FillH1(42, fNbMuonsMinus1);
  analysisManager->FillH1(43, fNbMuonsMinus2);  
  analysisManager->FillH1(44, fNbDelta12);
  analysisManager->FillH1(45, fNbDelta21); 
  analysisManager->FillH1(46, fNbPositronsB1);
  analysisManager->FillH1(47, fNbElectronsB1);
  analysisManager->FillH1(48, fNbPositrons12);
  analysisManager->FillH1(49, fNbElectrons12);  
  analysisManager->FillH1(50, fNbPositronsA2);
  analysisManager->FillH1(51, fNbElectronsA2); 
  for (G4int i = 84; i < 92; i = i + 1) {
    analysisManager->FillH1(i, fNbElectronsStrip1a[i+421]);
  }
  for (G4int j = 92; j < 100; j = j + 1) {
    analysisManager->FillH1(j, fNbElectronsStrip1b[j+413]);
  }
  for (G4int k = 100; k < 108; k = k + 1) {
    analysisManager->FillH1(k, fNbElectronsStrip2a[k+405]);
  }
  for (G4int l = 108; l < 116; l = l + 1) {
    analysisManager->FillH1(l, fNbElectronsStrip2b[l+397]);
  }

  for (G4int i1 = 123; i1 < 131; i1 = i1 + 1) {
    analysisManager->FillH1(i1, fNbHitsStrip1a[i1+382]);
  }
  for (G4int p1 = 131; p1 < 139; p1 = p1 + 1) {
    analysisManager->FillH1(p1, fNbHitsStrip1b[p1+374]);
  }
  for (G4int l1 = 139; l1 < 147; l1 = l1 + 1) {
    analysisManager->FillH1(l1, fNbHitsStrip2a[l1+366]);
  }
  for (G4int n1 = 147; n1 < 155; n1 = n1 + 1) {
    analysisManager->FillH1(n1, fNbHitsStrip2b[n1+358]);
  }

  // fill number of hits per strip ntuple
  for (G4int k2 = 0; k2 < 1017; k2++) {
    analysisManager->FillNtupleDColumn(k2+4*1017, fNbHitsStrip1a[k2]);
    analysisManager->FillNtupleDColumn(k2+5*1017, fNbHitsStrip1b[k2]);
    analysisManager->FillNtupleDColumn(k2+6*1017, fNbHitsStrip2a[k2]);
    analysisManager->FillNtupleDColumn(k2+7*1017, fNbHitsStrip2b[k2]);
  }
  analysisManager->AddNtupleRow();


  if ( analysisManager->IsActive() ) {
    analysisManager->Write();
    analysisManager->CloseFile();
  }      

  for (G4int k1 = 505; k1 < 513; k1 = k1 + 1) {
    G4cout 
      << "\n Number of hits in strip " << k1 
      << "a of Sensor 1 = " <<fNbHitsStrip1a[k1]
      << "\n Number of hits in strip " << k1 
      << "b of Sensor 1 = " <<fNbHitsStrip1b[k1]
      << "\n Number of hits in strip " << k1 
      << "a of Sensor 2 = " <<fNbHitsStrip2a[k1]
      << "\n Number of hits in strip " << k1 
      << "b of Sensor 2 = " <<fNbHitsStrip2b[k1]
      << G4endl; 
  }  
    
  if (particle->GetPDGCharge() == 0.) return;
   
  G4cout.precision(5);
  
  //track length
  //
  G4double trackLPerEvent = fTrackLength/nbEvents;
  G4double trackLS = fSecTrackLength/nbEvents;
  G4double trackLT = fTertTrackLength/nbEvents;
  G4double nbStepPerEvent = double(fNbSteps)/nbEvents;
  G4double stepSize = fTrackLength/fNbSteps;
  G4double trackSecDet = fSecDetTrackLength/nbEvents;
  
  G4cout 
    << "\n TrackLength= " 
    << G4BestUnit(trackLPerEvent, "Length")
    << "\n Secondary TrackLength= " 
    << G4BestUnit(trackLS, "Length")
    << "\n Tertiary TrackLength= " 
    << G4BestUnit(trackLT, "Length")
    << "\t nb of steps = " << nbStepPerEvent
    << "  stepSize= " << G4BestUnit(stepSize, "Length")
    << "  stepSize for secondaries created in detectors = " << G4BestUnit(trackSecDet, "Length")
    << G4endl;
      
  //charged secondaries (ionization, direct pair production)
  //
  G4double energy1PerEvent = fEnergyCharged1/nbEvents;
  G4double energy2PerEvent = fEnergyCharged2/nbEvents;
  G4double nb1PerEvent = double(fNbCharged1)/nbEvents;
  G4double nb2PerEvent = double(fNbCharged2)/nbEvents;
  G4double meanEkin1 = 0.;
  if (fNbCharged1) meanEkin1 = fEnergyCharged1/fNbCharged1;
  G4double meanEkin2 = 0.;
  if (fNbCharged2) meanEkin2 = fEnergyCharged2/fNbCharged2;
  
  G4cout 
    << "\n d-rays  : eLoss1/primary= " 
    << G4BestUnit(energy1PerEvent, "Energy")
    << "\n d-rays  : eLoss2/primary= " 
    << G4BestUnit(energy2PerEvent, "Energy")
    << "\t  nb1 of d-rays= " << nb1PerEvent
    << "\t  nb2 of d-rays= " << nb2PerEvent
    << G4endl;
         
  //neutral secondaries (bremsstrahlung, pixe)
  //
  energy1PerEvent = fEnergyNeutral1/nbEvents;
  energy2PerEvent = fEnergyNeutral2/nbEvents;
  nb1PerEvent = double(fNbNeutral1)/nbEvents;
  nb2PerEvent = double(fNbNeutral2)/nbEvents;
  meanEkin1 = 0.;
  if (fNbNeutral1) meanEkin1 = fEnergyNeutral1/fNbNeutral1;
  meanEkin2 = 0.;
  if (fNbNeutral2) meanEkin2 = fEnergyNeutral2/fNbNeutral2;
  
  G4cout 
    << "\n gamma   : eLoss1/primary= " 
    << G4BestUnit(energy1PerEvent, "Energy")
    << "\n gamma   : eLoss2/primary= " 
    << G4BestUnit(energy2PerEvent, "Energy")
    << "\t  nb1 of gammas= " << nb1PerEvent
    << "\t  nb2 of gammas= " << nb2PerEvent
    << "  <Tkin1>= " << G4BestUnit(meanEkin1, "Energy")
    << "  Tmin1= "   << G4BestUnit(fEmin1[1],  "Energy")
    << "  Tmax1= "   << G4BestUnit(fEmax1[1],  "Energy")
    << "  <Tkin2>= " << G4BestUnit(meanEkin2, "Energy")
    << "  Tmin2= "   << G4BestUnit(fEmin2[1],  "Energy")
    << "  Tmax2= "   << G4BestUnit(fEmax2[1],  "Energy")
    << G4endl;
    

  G4EmCalculator emCal;
  
  //local energy deposit for Detector 1
  //
  energy1PerEvent = fEnergyDeposit1/nbEvents;
  //
  G4double r01  = emCal.GetRangeFromRestricteDEDX(ePrimary,particle,material);  
  G4double r11 = r01 - trackLPerEvent;
  G4double etry1 = ePrimary - energy1PerEvent;  
  G4double efinal1 = 0.;
  if (r11 > 0.) efinal1 = GetEnergyFromRestrictedRange(r11,particle,material,etry1);
  G4double dEtable1 = ePrimary - efinal1;
  G4double ratio1 = 0.;
  if (dEtable1 > 0.) ratio1 = energy1PerEvent/dEtable1;
    
  G4cout 
    << "\n deposit1 : eLoss1/primary= " 
    << G4BestUnit(energy1PerEvent, "Energy")
    << "\t <dEcut1 > table= " 
    << G4BestUnit(dEtable1, "Energy")
    << "   ---> simul1/reference= " << ratio1          
    << G4endl;
    
  //total energy transferred for Detector 1
  //
  G4double energyTotal1 = fEnergyDeposit1 + fEnergyCharged1 + fEnergyNeutral1;
  energy1PerEvent = energyTotal1/nbEvents;
  //
  r01  = emCal.GetCSDARange(ePrimary,particle,material);  
  r11 = r01 - trackLPerEvent;
  etry1 = ePrimary - energy1PerEvent;
  efinal1 = 0.;
  if (r11 > 0.) efinal1 = GetEnergyFromCSDARange(r11,particle,material,etry1);
  dEtable1 = ePrimary - efinal1;
  ratio1 = 0.;
  if (dEtable1 > 0.) ratio1 = energy1PerEvent/dEtable1;
    
  G4cout 
    << "\n total   : eLoss1/primary= " 
    << G4BestUnit(energy1PerEvent, "Energy")
    << "\t <dEfull1> table= " 
    << G4BestUnit(dEtable1, "Energy")
    << "   ---> simul1/reference= " << ratio1           
    << G4endl; 

  //local energy deposit for Detector 2
  //
  energy2PerEvent = fEnergyDeposit2/nbEvents;
  //
  G4double r02  = emCal.GetRangeFromRestricteDEDX(ePrimary,particle,material);  
  G4double r12 = r02 - trackLPerEvent;
  G4double etry2 = ePrimary - energy2PerEvent;  
  G4double efinal2 = 0.;
  if (r12 > 0.) efinal2 = GetEnergyFromRestrictedRange(r12,particle,material,etry2);
  G4double dEtable2 = ePrimary - efinal2;
  G4double ratio2 = 0.;
  if (dEtable2 > 0.) ratio2 = energy2PerEvent/dEtable2;
    
  G4cout 
    << "\n deposit2 : eLoss2/primary= " 
    << G4BestUnit(energy2PerEvent, "Energy")
    << "\t <dEcut2 > table= " 
    << G4BestUnit(dEtable2, "Energy")
    << "   ---> simul2/reference= " << ratio2          
    << G4endl;
    
  //total energy transferred for Detector 2
  //
  G4double energyTotal2 = fEnergyDeposit2 + fEnergyCharged2 + fEnergyNeutral2;
  energy2PerEvent = energyTotal2/nbEvents;
  //
  r02  = emCal.GetCSDARange(ePrimary,particle,material);  
  r12 = r02 - trackLPerEvent;
  etry2 = ePrimary - energy2PerEvent;
  efinal2 = 0.;
  if (r12 > 0.) efinal2 = GetEnergyFromCSDARange(r12,particle,material,etry2);
  dEtable2 = ePrimary - efinal2;
  ratio2 = 0.;
  if (dEtable2 > 0.) ratio2 = energy2PerEvent/dEtable2;
    
  G4cout 
    << "\n total   : eLoss2/primary= " 
    << G4BestUnit(energy2PerEvent, "Energy")
    << "\t <dEfull2> table= " 
    << G4BestUnit(dEtable2, "Energy")
    << "   ---> simul2/reference= " << ratio2           
    << G4endl; 

  G4cout.precision(prec);

  // show Rndm status
  CLHEP::HepRandom::showEngineStatus();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double RunAction::GetEnergyFromRestrictedRange(G4double range,
            G4ParticleDefinition* particle, G4Material* material, G4double Etry)
{
  G4EmCalculator emCal;
    
  G4double Energy = Etry, dE = 0., dEdx;
  G4double r, dr;
  G4double err  = 1., errmax = 0.00001;
  G4int    iter = 0 , itermax = 10;  
  while (err > errmax && iter < itermax) {
    iter++;
    Energy -= dE;
    r = emCal.GetRangeFromRestricteDEDX(Energy,particle,material);
    dr = r - range;          
    dEdx = emCal.GetDEDX(Energy,particle,material);
    dE = dEdx*dr;
    err = std::abs(dE)/Energy;    
  }
  if (iter == itermax) {
    G4cout 
    << "\n  ---> warning: RunAction::GetEnergyFromRestRange() did not converge"
    << "   Etry = " << G4BestUnit(Etry,"Energy")
    << "   Energy = " << G4BestUnit(Energy,"Energy")
    << "   err = " << err
    << "   iter = " << iter << G4endl;
  }         
         
  return Energy;         
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double RunAction::GetEnergyFromCSDARange(G4double range,
            G4ParticleDefinition* particle, G4Material* material, G4double Etry)
{
  G4EmCalculator emCal;
    
  G4double Energy = Etry, dE = 0., dEdx;
  G4double r, dr;
  G4double err  = 1., errmax = 0.00001;
  G4int    iter = 0 , itermax = 10;  
  while (err > errmax && iter < itermax) {
    iter++;
    Energy -= dE;
    r = emCal.GetCSDARange(Energy,particle,material);
    dr = r - range;          
    dEdx = emCal.ComputeTotalDEDX(Energy,particle,material);
    dE = dEdx*dr;
    err = std::abs(dE)/Energy;
  }
  if (iter == itermax) {
    G4cout 
    << "\n  ---> warning: RunAction::GetEnergyFromCSDARange() did not converge"
    << "   Etry = " << G4BestUnit(Etry,"Energy")
    << "   Energy = " << G4BestUnit(Energy,"Energy")
    << "   err = " << err
    << "   iter = " << iter << G4endl;
  }         
         
  return Energy;         
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
