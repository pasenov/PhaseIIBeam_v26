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
/// \file electromagnetic/TestEm18/include/EventAction.hh
/// \brief Definition of the EventAction class
//
// $Id: EventAction.hh 82401 2014-06-18 14:43:54Z gcosmo $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef EventAction_h
#define EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

#include "RunAction.hh"
#include "HistoManager.hh"

#include "G4Event.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4EmCalculator.hh"
#include "G4Step.hh"

#include "Randomize.hh"
#include <iomanip>

class RunAction;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class EventAction : public G4UserEventAction
{
  public:
    EventAction(RunAction*);
   ~EventAction();

  public:
    virtual void BeginOfEventAction(const G4Event* evt);
    virtual void   EndOfEventAction(const G4Event*);

    
    void AddEnergyDeposit1(G4double edep)                     {fEnergyDeposit1  += edep;};
    void AddEnergyDeposit2(G4double edep)                     {fEnergyDeposit2  += edep;};
    void AddPrimaryTrackLength(G4double track)                {fPrimaryTrackLength  += track;};
    void AddSecondary1(G4double ekin)                         {fEnergySecondary1  += ekin;};
    void AddSecondary2(G4double ekin)                         {fEnergySecondary2  += ekin;};
    void AddSecondaryxPolarization(G4double spolarization)    {fSecondaryxPolarization  += spolarization;};
    void AddSecondaryyPolarization(G4double spolarization)    {fSecondaryyPolarization  += spolarization;};
    void AddSecondaryzPolarization(G4double spolarization)    {fSecondaryzPolarization  += spolarization;};
    void AddSecondaryTrackLength(G4double track)              {fSecondaryTrackLength  += track;};
    void AddSecondaryDetTrackLength(G4double track)           {fSecondaryDetTrackLength  += track;};
    void AddTertiary1(G4double ekin)                          {fEnergyTertiary1  += ekin;};
    void AddTertiary2(G4double ekin)                          {fEnergyTertiary2  += ekin;};
    void AddTertiaryxPolarization(G4double spolarization)     {fTertiaryxPolarization  += spolarization;};
    void AddTertiaryyPolarization(G4double spolarization)     {fTertiaryyPolarization  += spolarization;};
    void AddTertiaryzPolarization(G4double spolarization)     {fTertiaryzPolarization  += spolarization;};
    void AddTertiaryTrackLength(G4double track)               {fTertiaryTrackLength  += track;};
    void TrackCheck(G4double trackcheck)                      {fTrack1 = fTrack2; fTrack2 = trackcheck; 
                                                               if (fTrack1 > fTrack2) {
                                                                   G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
                                                                   analysisManager->FillH1(31, fTrack1);
                                                                   //G4cout 
     									//<< "\n Last track length of secondary created in detector calculated = " 
     									//<< G4BestUnit(fTrack1, "Length")
     									//<< G4endl;
                                                               }};
    void AddEnergyStrip1a(G4double edep, G4int i)  {fEnergyStrip1a[i] += edep;};
    void AddEnergyStrip1b(G4double edep, G4int i)  {fEnergyStrip1b[i] += edep;};
    void AddEnergyStrip2a(G4double edep, G4int i)  {fEnergyStrip2a[i] += edep;};
    void AddEnergyStrip2b(G4double edep, G4int i)  {fEnergyStrip2b[i] += edep;};

    void AddEnergyPixel(G4double edep, G4int i, G4int j, G4int k, G4int l)  {fEnergyPixel[i][j][k][l] += edep;};

    void AddMomentumDirection1(G4double dirx, G4double diry, G4double dirz)  {fMomDir1x += dirx; fMomDir1y += diry; fMomDir1z += dirz;};
    void AddMomentumDirection2(G4double dirx, G4double diry, G4double dirz)  {fMomDir2x += dirx; fMomDir2y += diry; fMomDir2z += dirz;};

    void AddPointPix1ent(G4ThreeVector vec)  {fPointPix1ent += vec;};
    void AddPointPix1ex(G4ThreeVector vec)  {fPointPix1ex += vec;};
    void AddPointDet1ent(G4ThreeVector vec)  {fPointDet1ent += vec;};
    void AddPointDet1ex(G4ThreeVector vec)  {fPointDet1ex += vec;};
    void AddPointDet2ent(G4ThreeVector vec)  {fPointDet2ent += vec;};
    void AddPointDet2ex(G4ThreeVector vec)  {fPointDet2ex += vec;};
    void AddPointPix2ent(G4ThreeVector vec)  {fPointPix2ent += vec;};
    void AddPointPix2ex(G4ThreeVector vec)  {fPointPix2ex += vec;};
        
  private:
    RunAction*    fRunAction;
    
    G4double      fEnergyDeposit1;
    G4double      fEnergyDeposit2;
    G4double      fPrimaryTrackLength;
    G4double      fEnergySecondary1;
    G4double      fEnergySecondary2;       
    G4double      fSecondaryxPolarization;  
    G4double      fSecondaryyPolarization;  
    G4double      fSecondaryzPolarization; 
    G4double      fSecondaryTrackLength; 
    G4double      fSecondaryDetTrackLength;
    G4double      fEnergyTertiary1;
    G4double      fEnergyTertiary2;   
    G4double      fTertiaryxPolarization;   
    G4double      fTertiaryyPolarization; 
    G4double      fTertiaryzPolarization; 
    G4double      fTertiaryTrackLength;
    G4double      fEnergyStrip1a[1017];
    G4double      fEnergyStrip1b[1017];
    G4double      fEnergyStrip2a[1017];
    G4double      fEnergyStrip2b[1017];
    G4double      fMomDir1x, fMomDir1y, fMomDir1z;
    G4double      fMomDir2x, fMomDir2y, fMomDir2z;
    G4double      fEnergyPixel[3][17][81][53];

    G4int         fCharge1, fCharge2, fCS1a, fCS1b, fCS2a, fCS2b, fChargePix;

    G4int         fNbHitsStrip1a[1017];
    G4int         fNbHitsStrip1b[1017];
    G4int         fNbHitsStrip2a[1017];
    G4int         fNbHitsStrip2b[1017];
    G4int         fChargeStrip1a[1017];
    G4int         fChargeStrip1b[1017];
    G4int         fChargeStrip2a[1017];
    G4int         fChargeStrip2b[1017];
    G4int         fChargePixel[3][17][81][53];
    G4int         fNbHitsPixel[3][17][81][53];

    G4int         fHitsMultiplicity1b;
    G4int         fHitsMultiplicity2b;
    G4int         fHitsMultiplicityPix[3];

    G4ThreeVector fPointPix1ent, fPointPix1mid, fPointPix1ex, fPointDet1ent, fPointDet1mid, fPointDet1ex, fPointDet2ent, fPointDet2mid, fPointDet2ex, fPointPix2ent, fPointPix2mid, fPointPix2ex;

    //Let A = fPointPix1mid, B = fPointDet1mid, C = fPointDet2mid, D = fPointPix2mid. Let dB = distance between B and AD and dC = distance between C and AD.
    //Let B1 be the AD point with same z as B. Let C1 be the AD point with the same z as C.
    G4double dB, dC;
    G4double dBAD, dCAD, dAD;
    G4ThreeVector xBA, xBD, xCA, xCD, xAD;
    G4ThreeVector xBAD, xCAD;
    G4ThreeVector B1, C1;
    G4double B1Bx, B1By, C1Cx, C1Cy;

    //primary track's deflection angle in radians
    G4double ftheta;

    G4double      fTrack1;
    G4double      fTrack2;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

    
