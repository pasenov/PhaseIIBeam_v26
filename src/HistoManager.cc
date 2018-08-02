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
//
// $Id: HistoManager.cc 72242 2013-07-12 08:44:19Z gcosmo $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include <CLHEP/Units/SystemOfUnits.h>
#include <sstream>

#include "HistoManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HistoManager::HistoManager()
  :fFileName("phaseiiv1")
{
  Book();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HistoManager::~HistoManager()
{
  delete G4AnalysisManager::Instance();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::Book()
{
  // Create or get analysis manager
  // The choice of analysis technology is done via selection of a namespace
  // in HistoManager.hh
  // Creating a tree container to handle histograms and ntuples.
  // This tree is associated to an output file.
  //
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->SetFileName(fFileName);
  analysisManager->SetVerboseLevel(1);
  analysisManager->SetActivation(true);   //enable inactivation of histograms

  
  // Define histograms start values
  // const G4int kNbofStrips = 1016;
  const G4int kMaxHisto = 199;

  const G4String id[] = { "0", "1", "2", "3" , "4", "5", "6", "7", "8", "9" , "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29", "30", "31", "32", "33", "34", "35", "36", "37", "38", "39", "40", "41", "42", "43", "44", "45", "46", "47", "48", "49", "50", "51", "52", "53", "54", "55", "56", "57", "58", "59", "60", "61", "62", "63", "64", "65", "66", "67", "68", "69", "70", "71", "72", "73", "74", "75", "76", "77", "78", "79", "80", "81", "82", "83", "84", "85", "86", "87", "88", "89", "90", "91", "92", "93", "94", "95", "96", "97", "98", "99", "100", "101", "102", "103", "104", "105", "106", "107", "108", "109", "110", "111", "112", "113", "114", "115", "116", "117", "118", "119", "120", "121", "122", "123", "124", "125", "126", "127", "128", "129", "130", "131", "132", "133", "134", "135", "136", "137", "138", "139", "140", "141", "142", "143", "144", "145", "146", "147", "148", "149", "150", "151", "152", "153", "154", "155", "156", "157", "158", "159", "160", "161", "162", "163", "164", "165", "166", "167", "168", "169", "170", "171", "172", "173", "174", "175", "176", "177", "178", "179", "180", "181", "182", "183", "184", "185", "186", "187", "188", "189", "190", "191", "192", "193", "194", "195", "196", "197", "198" };
  const G4String title[] =
                { "dummy",                                                       //0
                  "continuous energy loss along primary track for Sensor 1",   //1
                  "continuous energy loss along primary track for Sensor 2",   //2
                  "energy from secondaries for Sensor 1",                      //3
                  "energy from secondaries for Sensor 2",                      //4
                  "energy from tertiaries for Sensor 1",                       //5
                  "energy from tertiaries for Sensor 2",                       //6
                  "total energy lost by primary track for Sensor 1",           //7
                  "total energy lost by primary track for Sensor 2",           //8
                  "energy spectrum of secondary e-+ for Sensor 1",             //9
                  "energy spectrum of secondary e-+ for Sensor 2",             //10
                  "energy spectrum of secondary gamma for Sensor 1",           //11
                  "energy spectrum of secondary gamma for Sensor 2",           //12
                  "energy spectrum of tertiary e-+ for Sensor 1",              //13
                  "energy spectrum of tertiary e-+ for Sensor 2",              //14
                  "energy spectrum of tertiary gamma for Sensor 1",            //15
                  "energy spectrum of tertiary gamma for Sensor 2",            //16
                  "step size for primary",                                       //17
                  "step size for secondaries",                                   //18
                  "step size for tertiaries",                                    //19
                  "step size for secondaries created in sensors",              //20
                  "x-polarization of secondaries",                               //21
                  "y-polarization of secondaries",                               //22
                  "z-polarization of secondaries",                               //23
                  "x-polarization of tertiaries",                                //24
                  "y-polarization of tertiaries",                                //25
                  "z-polarization of tertiaries",                                //26
                  "track length for primary",                                    //27
                  "total track length of secondaries",                           //28
                  "total track length of tertiaries",                            //29
                  "total track length of secondaries created in sensors",      //30
                  "track length of secondaries created in sensors",            //31
                  "number of positrons created in Sensor 1",       		 //32
                  "number of positrons created in Sensor 2",       		 //33
                  "number of electrons created in Sensor 1",       		 //34
                  "number of electrons created in Sensor 2",       		 //35
                  "number of pi+ created in Sensor 1",             		 //36
                  "number of pi+ created in Sensor 2",             		 //37
                  "number of pi- created in Sensor 1",              		 //38
                  "number of pi- created in Sensor 2",             		 //39
                  "number of mu+ created in Sensor 1",             		 //40
                  "number of mu+ created in Sensor 2",             		 //41
                  "number of mu- created in Sensor 1",             		 //42
                  "number of mu- created in Sensor 2",              		 //43
		  "number of delta e- created in Sensor 1 that have reached Sensor 2", //44
		  "number of delta e- created in Sensor 2 that have reached Sensor 1", //45
                  "number of e+ created before Sensor 1",               	 //46
                  "number of e- created before Sensor 1",               	 //47
		  "number of e+ created between Sensor 1 and Sensor 2", 	 //48
		  "number of e- created between Sensor 1 and Sensor 2", 	 //49
                  "number of e+ created after Sensor 2",               	 //50
                  "number of e- created after Sensor 2",                	 //51
                  "continuous energy loss along primary track for strip 505 a Sensor 1",   //52
                  "continuous energy loss along primary track for strip 506 a Sensor 1",   //53
                  "continuous energy loss along primary track for strip 507 a Sensor 1",   //54
                  "continuous energy loss along primary track for strip 508 a Sensor 1",   //55
                  "continuous energy loss along primary track for strip 509 a Sensor 1",   //56
                  "continuous energy loss along primary track for strip 510 a Sensor 1",   //57
                  "continuous energy loss along primary track for strip 511 a Sensor 1",   //58
                  "continuous energy loss along primary track for strip 512 a Sensor 1",   //59
                  "continuous energy loss along primary track for strip 505 b Sensor 1",   //60
                  "continuous energy loss along primary track for strip 506 b Sensor 1",   //61
                  "continuous energy loss along primary track for strip 507 b Sensor 1",   //62
                  "continuous energy loss along primary track for strip 508 b Sensor 1",   //63
                  "continuous energy loss along primary track for strip 509 b Sensor 1",   //64
                  "continuous energy loss along primary track for strip 510 b Sensor 1",   //65
                  "continuous energy loss along primary track for strip 511 b Sensor 1",   //66
                  "continuous energy loss along primary track for strip 512 b Sensor 1",   //67
                  "continuous energy loss along primary track for strip 505 a Sensor 2",   //68
                  "continuous energy loss along primary track for strip 506 a Sensor 2",   //69
                  "continuous energy loss along primary track for strip 507 a Sensor 2",   //70
                  "continuous energy loss along primary track for strip 508 a Sensor 2",   //71
                  "continuous energy loss along primary track for strip 509 a Sensor 2",   //72
                  "continuous energy loss along primary track for strip 510 a Sensor 2",   //73
                  "continuous energy loss along primary track for strip 511 a Sensor 2",   //74
                  "continuous energy loss along primary track for strip 512 a Sensor 2",   //75
                  "continuous energy loss along primary track for strip 505 b Sensor 2",   //76
                  "continuous energy loss along primary track for strip 506 b Sensor 2",   //77
                  "continuous energy loss along primary track for strip 507 b Sensor 2",   //78
                  "continuous energy loss along primary track for strip 508 b Sensor 2",   //79
                  "continuous energy loss along primary track for strip 509 b Sensor 2",   //80
                  "continuous energy loss along primary track for strip 510 b Sensor 2",   //81
                  "continuous energy loss along primary track for strip 511 b Sensor 2",   //82
                  "continuous energy loss along primary track for strip 512 b Sensor 2",   //83
                  "number of electrons created in strip 505 a Sensor 1",   //84
                  "number of electrons created in strip 506 a Sensor 1",   //85
                  "number of electrons created in strip 507 a Sensor 1",   //86
                  "number of electrons created in strip 508 a Sensor 1",   //87
                  "number of electrons created in strip 509 a Sensor 1",   //88
                  "number of electrons created in strip 510 a Sensor 1",   //89
                  "number of electrons created in strip 511 a Sensor 1",   //90
                  "number of electrons created in strip 512 a Sensor 1",   //91
                  "number of electrons created in strip 505 b Sensor 1",   //92
                  "number of electrons created in strip 506 b Sensor 1",   //93
                  "number of electrons created in strip 507 b Sensor 1",   //94
                  "number of electrons created in strip 508 b Sensor 1",   //95
                  "number of electrons created in strip 509 b Sensor 1",   //96
                  "number of electrons created in strip 510 b Sensor 1",   //97
                  "number of electrons created in strip 511 b Sensor 1",   //98
                  "number of electrons created in strip 512 b Sensor 1",   //99
                  "number of electrons created in strip 505 a Sensor 2",   //100
                  "number of electrons created in strip 506 a Sensor 2",   //101
                  "number of electrons created in strip 507 a Sensor 2",   //102
                  "number of electrons created in strip 508 a Sensor 2",   //103
                  "number of electrons created in strip 509 a Sensor 2",   //104
                  "number of electrons created in strip 510 a Sensor 2",   //105
                  "number of electrons created in strip 511 a Sensor 2",   //106
                  "number of electrons created in strip 512 a Sensor 2",   //107
                  "number of electrons created in strip 505 b Sensor 2",   //108
                  "number of electrons created in strip 506 b Sensor 2",   //109
                  "number of electrons created in strip 507 b Sensor 2",   //110
                  "number of electrons created in strip 508 b Sensor 2",   //111
                  "number of electrons created in strip 509 b Sensor 2",   //112
                  "number of electrons created in strip 510 b Sensor 2",   //113
                  "number of electrons created in strip 511 b Sensor 2",   //114
                  "number of electrons created in strip 512 b Sensor 2",   //115
                  "deflection angle of primary particle (rad)",              //116
                  "distance between midpoint of primary inside Sensor 1 and AD",   //117
                  "distance between midpoint of primary inside Sensor 2 and AD",   //118
                  "2S first-plane residuals in x-direction",   //119
                  "2S first-plane residuals in y-direction",   //120
                  "2S second-plane residuals in x-direction",   //121
                  "2S second-plane residuals in y-direction",   //122
                  "number of hits for strip 505 a Sensor 1",   //123
                  "number of hits for strip 506 a Sensor 1",   //124
                  "number of hits for strip 507 a Sensor 1",   //125
                  "number of hits for strip 508 a Sensor 1",   //126
                  "number of hits for strip 509 a Sensor 1",   //127
                  "number of hits for strip 510 a Sensor 1",   //128
                  "number of hits for strip 511 a Sensor 1",   //129
                  "number of hits for strip 512 a Sensor 1",   //130
                  "number of hits for strip 505 b Sensor 1",   //131
                  "number of hits for strip 506 b Sensor 1",   //132
                  "number of hits for strip 507 b Sensor 1",   //133
                  "number of hits for strip 508 b Sensor 1",   //134
                  "number of hits for strip 509 b Sensor 1",   //135
                  "number of hits for strip 510 b Sensor 1",   //136
                  "number of hits for strip 511 b Sensor 1",   //137
                  "number of hits for strip 512 b Sensor 1",   //138
                  "number of hits for strip 505 a Sensor 2",   //139
                  "number of hits for strip 506 a Sensor 2",   //140
                  "number of hits for strip 507 a Sensor 2",   //141
                  "number of hits for strip 508 a Sensor 2",   //142
                  "number of hits for strip 509 a Sensor 2",   //143
                  "number of hits for strip 510 a Sensor 2",   //144
                  "number of hits for strip 511 a Sensor 2",   //145
                  "number of hits for strip 512 a Sensor 2",   //146
                  "number of hits for strip 505 b Sensor 2",   //147
                  "number of hits for strip 506 b Sensor 2",   //148
                  "number of hits for strip 507 b Sensor 2",   //149
                  "number of hits for strip 508 b Sensor 2",   //150
                  "number of hits for strip 509 b Sensor 2",   //151
                  "number of hits for strip 510 b Sensor 2",   //152
                  "number of hits for strip 511 b Sensor 2",   //153
                  "number of hits for strip 512 b Sensor 2",   //154
                  "charge in Sensor 1 corresponding to the deposited energy in Sensor 1",   //155
                  "charge in Sensor 2 corresponding to the deposited energy in Sensor 2",   //156
                  "charge in strip 505 a Sensor 1",   //157
                  "charge in strip 506 a Sensor 1",   //158
                  "charge in strip 507 a Sensor 1",   //159
                  "charge in strip 508 a Sensor 1",   //160
                  "charge in strip 509 a Sensor 1",   //161
                  "charge in strip 510 a Sensor 1",   //162
                  "charge in strip 511 a Sensor 1",   //163
                  "charge in strip 512 a Sensor 1",   //164
                  "charge in strip 505 b Sensor 1",   //165
                  "charge in strip 506 b Sensor 1",   //166
                  "charge in strip 507 b Sensor 1",   //167
                  "charge in strip 508 b Sensor 1",   //168
                  "charge in strip 509 b Sensor 1",   //169
                  "charge in strip 510 b Sensor 1",   //170
                  "charge in strip 511 b Sensor 1",   //171
                  "charge in strip 512 b Sensor 1",   //172
                  "charge in strip 505 a Sensor 2",   //173
                  "charge in strip 506 a Sensor 2",   //174
                  "charge in strip 507 a Sensor 2",   //175
                  "charge in strip 508 a Sensor 2",   //176
                  "charge in strip 509 a Sensor 2",   //177
                  "charge in strip 510 a Sensor 2",   //178
                  "charge in strip 511 a Sensor 2",   //179
                  "charge in strip 512 a Sensor 2",   //180
                  "charge in strip 505 b Sensor 2",   //181
                  "charge in strip 506 b Sensor 2",   //182
                  "charge in strip 507 b Sensor 2",   //183
                  "charge in strip 508 b Sensor 2",   //184
                  "charge in strip 509 b Sensor 2",   //185
                  "charge in strip 510 b Sensor 2",   //186
                  "charge in strip 511 b Sensor 2",   //187
                  "charge in strip 512 b Sensor 2",   //188
                  "hit multiplicity per event: Sensor 1",   //189
                  "hit multiplicity per event: Sensor 2",   //190
                  "average x-position of primary in sensor 1",    //191
                  "average y-position of primary in sensor 1",    //192
                  "average z-position of primary in sensor 1",    //193
                  "average x-position of primary in sensor 1",    //194
                  "average y-position of primary in sensor 2",    //195
                  "average z-position of primary in sensor 2",    //196
                  "hit multiplicity per event: BPIX module 1",    //197
                  "hit multiplicity per event: BPIX module 2"     //198
                  //"midpoint of primary particle inside Pixel Detector 1",
                  //"midpoint of primary particle inside Strip Detector 1",
                  //"midpoint of primary particle inside Strip Detector 2",
                  //"midpoint of primary particle inside Pixel Detector 2"
                 };
            
  // Default values (to be reset via /analysis/h1/set command)               
  G4int nbins = 100;
  G4double vmin = -1.5;
  G4double vmax = 100.;

  // Create all histograms as inactivated 
  // as we have not yet set nbins, vmin, vmax
  for (G4int k=0; k<kMaxHisto; k++) {
    G4int ih = analysisManager->CreateH1(id[k], title[k], nbins, vmin, vmax);
    analysisManager->SetH1Activation(ih, false);
  }

  //Create H3 histograms
  //G4int ih1 = analysisManager->CreateH3("85", "midpoint of primary particle inside Pixel Detector 1", 10000, -2*cm, 2*cm, 10000, -2*cm, 2*cm, 10000, -2*cm, 2*cm, "cm", "cm", "cm");
  //analysisManager->SetH3Activation(ih1, false);
  //G4int ih2 = analysisManager->CreateH3("86", "midpoint of primary particle inside Strip Detector 1", 1000, -2*mm, 2*mm, 1000, -2*mm, 2*mm, 1000, -2*mm, 2*mm, "mm", "mm", "mm");
  //analysisManager->SetH3Activation(ih2, false);
  //G4int ih3 = analysisManager->CreateH3("87", "midpoint of primary particle inside Strip Detector 2", 1000, -2*mm, 2*mm, 1000, -2*mm, 2*mm, 1000, -2*mm, 2*mm, "mm", "mm", "mm");
  //analysisManager->SetH3Activation(ih3, false);
  //G4int ih4 = analysisManager->CreateH3("88", "midpoint of primary particle inside Pixel Detector 2", 10000, -2*cm, 2*cm, 10000, -2*cm, 2*cm, 10000, -2*cm, 2*cm, "cm", "cm", "cm");
  //analysisManager->SetH3Activation(ih4, false);

  // nTuples, column l: 1 to 1016 correspond to Detector 1, 2 to 2032 correspond to Detector 2
  //
  analysisManager->SetNtupleDirectoryName("ntuple");
  analysisManager->SetFirstNtupleId(1);       
  analysisManager->CreateNtuple("1", "Primary Particle Tuple 1");
  for (G4int l=0; l<1017; l++) {
     std::ostringstream os1;
     os1 <<"EdepStripSensor1a_" << l;
     analysisManager->CreateNtupleDColumn(os1.str());          //
  }
  for (G4int m1=0; m1<1017; m1++) {
     std::ostringstream os2;
     os2 <<"EdepStripSensor1b_" << m1;
     analysisManager->CreateNtupleDColumn(os2.str());          //
  }
  for (G4int ma=0; ma<1017; ma++) {
     std::ostringstream os3;
     os3 <<"EdepStripSensor2a_" << ma;
     analysisManager->CreateNtupleDColumn(os3.str());          //
  }
  for (G4int mb=0; mb<1017; mb++) {
     std::ostringstream os4;
     os4 <<"EdepStripSensor2b_" << mb;
     analysisManager->CreateNtupleDColumn(os4.str());          //
  }

  for (G4int l1=0; l1<1017; l1++) {
     std::ostringstream os5;
     os5 <<"NbHitsSensor1a_" << l1;
     analysisManager->CreateNtupleDColumn(os5.str());          //
  }
  for (G4int l2=0; l2<1017; l2++) {
     std::ostringstream os6;
     os6 <<"NbHitsSensor1b_" << l2;
     analysisManager->CreateNtupleDColumn(os6.str());          //
  }
  for (G4int l3=0; l3<1017; l3++) {
     std::ostringstream os7;
     os7 <<"NbHitsSensor2a_" << l3;
     analysisManager->CreateNtupleDColumn(os7.str());          //
  }
  for (G4int l4=0; l4<1017; l4++) {
     std::ostringstream os8;
     os8 <<"NbHitsSensor2b_" << l4;
     analysisManager->CreateNtupleDColumn(os8.str());          //
  }

  G4int NbBPIX = 2;
  G4int NbROC = 16;
  G4int NbRows = 80;
  G4int NbCols = 52;

  for (G4int ne1=0; ne1<=NbBPIX; ne1++) {
     for (G4int ne2=0; ne2<=NbROC; ne2++) {
        for (G4int ne3=0; ne3<=NbRows; ne3++) {
           for (G4int ne4=0; ne4<=NbCols; ne4++) {
              std::ostringstream os9;
     	      os9 <<"EdepPixel_" << ne1 <<"_" << ne2 <<"_" << ne3 << "_" << ne4;
     	      //analysisManager->CreateNtupleDColumn(os9.str());  //Name convention: EdepPixel_NbBPIX_NbROC_NbRow_NbCol, zeros are dummy
           }
        }
     }
  }

  for (G4int n1=0; n1<=NbBPIX; n1++) {
     for (G4int n2=0; n2<=NbROC; n2++) {
        for (G4int n3=0; n3<=NbRows; n3++) {
           for (G4int n4=0; n4<=NbCols; n4++) {
              std::ostringstream os10;
     	      os10 <<"NbHitsPixel_" << n1 <<"_" << n2 <<"_" << n3 << "_" << n4;
     	      //analysisManager->CreateNtupleDColumn(os10.str());  //Name convention: NbHitsPixel_NbBPIX_NbROC_NbRow_NbCol, zeros are dummy
           }
        }
     }
  }

  //for (G4int n=0; n<1017; n++) {
     //std::ostringstream os3;
     //os3 <<"NbElectronsStripSensor1_" << n;
     //analysisManager->CreateNtupleDColumn(os3.str());          //
  //}
  //for (G4int o=0; o<1017; o++) {
     //std::ostringstream os4;
     //os4 <<"NbElectronsStripSensor2_" << o;
     //analysisManager->CreateNtupleDColumn(os4.str());          //
  //}
  analysisManager->FinishNtuple();
  
  //analysisManager->SetNtupleActivation(false);  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
