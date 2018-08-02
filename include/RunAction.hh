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
/// \file electromagnetic/TestEm18/include/RunAction.hh
/// \brief Definition of the RunAction class
//
// $Id: RunAction.hh 66241 2012-12-13 18:34:42Z gunter $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef RunAction_h
#define RunAction_h 1

#include "G4UserRunAction.hh"

#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4Run;
class G4ParticleDefinition;
class G4Material;

class DetectorConstruction;
class PrimaryGeneratorAction;
class HistoManager;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class RunAction : public G4UserRunAction
{
  public:
    RunAction(DetectorConstruction*, PrimaryGeneratorAction*);
   ~RunAction();

  public:
    virtual void BeginOfRunAction(const G4Run*);
    virtual void   EndOfRunAction(const G4Run*);

    void AddEnergyDeposit1 (G4double edep)
                   {fEnergyDeposit1 += edep;};

    void AddEnergyDeposit2 (G4double edep)
                   {fEnergyDeposit2 += edep;};

    void AddTrackLength (G4double step)
                 {fTrackLength += step; fNbSteps++;};

    void AddSecTrackLength (G4double step)
                 {fSecTrackLength += step; fNbSteps++;};

    void AddSecDetTrackLength (G4double step)
                 {fSecDetTrackLength += step; fNbSteps++;};

    void AddTertTrackLength (G4double step)
                 {fTertTrackLength += step; fNbSteps++;};
                 
    void AddChargedSecondary1 (G4double ekin)
                 {fEnergyCharged1 += ekin; fNbCharged1++;
                  if (ekin<fEmin1[0]) fEmin1[0] = ekin;
                  if (ekin>fEmax1[0]) fEmax1[0] = ekin;
                 };

    void AddChargedSecondary2 (G4double ekin)
                 {fEnergyCharged2 += ekin; fNbCharged2++;
                  if (ekin<fEmin2[0]) fEmin2[0] = ekin;
                  if (ekin>fEmax2[0]) fEmax2[0] = ekin;
                 };
                 
    void AddNeutralSecondary1 (G4double ekin)
                 {fEnergyNeutral1 += ekin; fNbNeutral1++;
                  if (ekin<fEmin1[1]) fEmin1[1] = ekin;
                  if (ekin>fEmax1[1]) fEmax1[1] = ekin;
                 };

    void AddNeutralSecondary2 (G4double ekin)
                 {fEnergyNeutral2 += ekin; fNbNeutral2++;
                  if (ekin<fEmin2[1]) fEmin2[1] = ekin;
                  if (ekin>fEmax2[1]) fEmax2[1] = ekin;
                 };

    void AddChargedTertiary1 (G4double ekin)
                 {fEnergyCharged1 += ekin; fNbCharged1++;
                  if (ekin<fEmin1[0]) fEmin1[0] = ekin;
                  if (ekin>fEmax1[0]) fEmax1[0] = ekin;
                 };

    void AddChargedTertiary2 (G4double ekin)
                 {fEnergyCharged2 += ekin; fNbCharged2++;
                  if (ekin<fEmin2[0]) fEmin2[0] = ekin;
                  if (ekin>fEmax2[0]) fEmax2[0] = ekin;
                 };
                 
    void AddNeutralTertiary1 (G4double ekin)
                 {fEnergyNeutral1 += ekin; fNbNeutral1++;
                  if (ekin<fEmin1[1]) fEmin1[1] = ekin;
                  if (ekin>fEmax1[1]) fEmax1[1] = ekin;
                 };

    void AddNeutralTertiary2 (G4double ekin)
                 {fEnergyNeutral2 += ekin; fNbNeutral2++;
                  if (ekin<fEmin2[1]) fEmin2[1] = ekin;
                  if (ekin>fEmax2[1]) fEmax2[1] = ekin;
                 };

    void AddNbPositronsB1 ()
                 {fNbPositronsB1++;};

    void AddNbPositrons1 ()
                 {fNbPositrons1++;};

    void AddNbPositrons12 ()
                 {fNbPositrons12++;};

    void AddNbPositrons2 ()
                 {fNbPositrons2++;};

    void AddNbPositronsA2 ()
                 {fNbPositronsA2++;};

    void AddNbElectronsB1 ()
                 {fNbElectronsB1++;};

    void AddNbElectrons1 ()
                 {fNbElectrons1++;};

    void AddNbElectrons12 ()
                 {fNbElectrons12++;};

    void AddNbElectrons2 ()
                 {fNbElectrons2++;};

    void AddNbElectronsA2 ()
                 {fNbElectronsA2++;};

    void AddNbPionsPlus1 ()
                 {fNbPionsPlus1++;};

    void AddNbPionsPlus2 ()
                 {fNbPionsPlus2++;};

    void AddNbPionsMinus1 ()
                 {fNbPionsMinus1++;};

    void AddNbPionsMinus2 ()
                 {fNbPionsMinus2++;};

    void AddNbMuonsPlus1 ()
                 {fNbMuonsPlus1++;};

    void AddNbMuonsPlus2 ()
                 {fNbMuonsPlus2++;};

    void AddNbMuonsMinus1 ()
                 {fNbMuonsMinus1++;};

    void AddNbMuonsMinus2 ()
                 {fNbMuonsMinus2++;};

    void AddNbDelta12 ()
                 {fNbDelta12++;};

    void AddNbDelta21 ()
                 {fNbDelta21++;};

    void AddNbElectronsStrip1a(G4int i)  {fNbElectronsStrip1a[i] ++;};
    void AddNbElectronsStrip1b(G4int j)  {fNbElectronsStrip1b[j] ++;};
    void AddNbElectronsStrip2a(G4int k)  {fNbElectronsStrip2a[k] ++;};
    void AddNbElectronsStrip2b(G4int l)  {fNbElectronsStrip2b[l] ++;};

    void AddNbHitsStrip1a(G4int i1)  {fNbHitsStrip1a[i1] ++;};
    void AddNbHitsStrip1b(G4int j1)  {fNbHitsStrip1b[j1] ++;};
    void AddNbHitsStrip2a(G4int k1)  {fNbHitsStrip2a[k1] ++;};
    void AddNbHitsStrip2b(G4int l1)  {fNbHitsStrip2b[l1] ++;};

    void AddNbHitsPixel(G4int le1, G4int le2, G4int le3, G4int le4)  {fNbHitsPixel[le1][le2][le3][le4] ++;};

               
  public:
    G4double GetEnergyFromRestrictedRange
             (G4double,G4ParticleDefinition*,G4Material*,G4double);
                       
    G4double GetEnergyFromCSDARange
             (G4double,G4ParticleDefinition*,G4Material*,G4double);                 
                 
  private:
    G4double fEnergyDeposit1, fEnergyDeposit2;
    G4double fTrackLength, fSecTrackLength, fTertTrackLength, fSecDetTrackLength;
    G4double fEnergyCharged1, fEnergyCharged2, fEnergyNeutral1, fEnergyNeutral2;
    G4double fEmin1[2], fEmax1[2], fEmin2[2], fEmax2[2];
    G4int fNbElectronsStrip1a[1017];
    G4int fNbElectronsStrip1b[1017];
    G4int fNbElectronsStrip2a[1017];
    G4int fNbElectronsStrip2b[1017];
    G4int fNbHitsStrip1a[1017];
    G4int fNbHitsStrip1b[1017];
    G4int fNbHitsStrip2a[1017];
    G4int fNbHitsStrip2b[1017];
    G4int fNbHitsPixel[3][17][81][53];
    
    G4long   fNbSteps;
    G4int    fNbCharged1, fNbCharged2, fNbNeutral1, fNbNeutral2;
    G4int    fNbPositronsB1, fNbPositrons1, fNbPositrons12, fNbPositrons2, fNbPositronsA2, fNbElectronsB1, fNbElectrons1, fNbElectrons12, fNbElectrons2, fNbElectronsA2, fNbPionsPlus1, fNbPionsPlus2, fNbPionsMinus1, fNbPionsMinus2, fNbMuonsPlus1, fNbMuonsPlus2, fNbMuonsMinus1, fNbMuonsMinus2;
    G4int    fNbDelta12, fNbDelta21;

    DetectorConstruction*   fDetector;
    PrimaryGeneratorAction* fPrimary;
    HistoManager*           fHistoManager;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

