#include "./SetupTree.C"

#include <TStyle.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1D.h>
#include <TGraphAsymmErrors.h>
#include <TMultiGraph.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TROOT.h>
#include <THStack.h>
#include <TMath.h>
#include <TText.h>
#include <TPad.h>
#include <TPaveText.h>
#include <TLorentzVector.h>

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>

using namespace std;
TString _dataset = "SingleMu";
TString _trig = "HLT_IsoMu24_eta2p1_v*"; //single mu trigger for mumu, 4mu, 3mu1e
double intLumi = 14817;
TH1D* fillHist( TString, TString, TString, TString, double, int, double, double, bool );
bool getDYSpectrum( TString, double );
void PrintAllDetails( void );

bool reorder(const TLorentzVector &a, const TLorentzVector &b) {
  return a.Pt() > b.Pt();
}

double deltaPhi(double phi1, double phi2) {
  double result = phi1 - phi2;
  while (result > M_PI) result -= 2*M_PI;
  while (result <= -M_PI) result += 2*M_PI;
  return result;
}

double deltaR(double eta1, double phi1, double eta2, double phi2) {
  double deta = eta1 - eta2;
  double dphi = deltaPhi(phi1, phi2);
  return sqrt(deta*deta + dphi*dphi);
}

ofstream fout;
bool isDataPrint = true;
TString optPlot = "mass";

void massSpec( double _ptCut = 30.0, TString _channel = "mumu" ) 
{
    TString _ptCutname = "";
    ostringstream name0;
    name0 << "_pt" << (int)_ptCut;
    _ptCutname = (TString) name0.str();
    cout << "_ptCutname = " << _ptCutname << endl;

    TFile* f = new TFile("rootfiles/massSpec_"+_channel+"_"+optPlot+_ptCutname+".root", "recreate");
    f->cd();
    gStyle->SetOptStat(0);

    const int nbin = 120;
    double xbin1 = 60;
    double xbin2 = 120;

    // data
    TH1D* hdata1 = fillHist(_dataset+"_Run2011B", _trig, optPlot, _channel, _ptCut, nbin, xbin1, xbin2, false);
    // DY MC 
    TH1D* hDY = fillHist("DYJets", _trig, optPlot, _channel, _ptCut, nbin, xbin1, xbin2, true);

    f->Write();
}

TH1D* fillHist( TString _fname, TString trig, TString _optPlot, TString _channel, double _ptCut, int nbin, double xbin1, double xbin2, bool isMC )
{
    TChain* sample = new TChain("recoTree/DYTree");
    SetupTree(_fname, sample);

    TH1D* _h1 = new TH1D("hist_"+_fname, "hist_"+_fname, nbin, xbin1, xbin2);

    double _etaCut = 2.4;

    for( int i = 0; i < sample->GetEntries(); i++ ) {
        sample->GetEntry(i);

        if( i % 100000 == 0 ) cout << _fname << " = " << i << endl;
	//if( i > 1000 ) break;

	/*
	if( isMC && isDY ) {
          // gen. mass
	  vector<TLorentzVector> genLeptons;
	  TLorentzVector tmpGen;
          for( int j = 0; j < gnpair; j++ ) {
	    if( fabs(genLepton_mother[j]) == 23 && genLepton_status[j] == 3 ) {
	      double parMass = -1;
	      if( fabs(genLepton_ID[j]) == 11 ) {
	        parMass = 0.000511;
	      }
	      if( fabs(genLepton_ID[j]) == 13 ) {
	        parMass = 0.105658;
	      }
	      double rpx = genLepton_px[j];
	      double rpy = genLepton_py[j];
	      double rpz = genLepton_pz[j];
	      double rp = rpx*rpx + rpy*rpy + rpz*rpz;
	      double re = sqrt(parMass*parMass + rp);

	      tmpGen.SetPxPyPzE(rpx, rpy, rpz, re);
	      genLeptons.push_back(tmpGen);
	      if( genLeptons.size() == 2 ) break;
	    }
	  }
          TLorentzVector genCand = genLeptons[0] + genLeptons[1];
	  double genMass = genCand.M();
	  if( getDYSpectrum(_fname, genMass) ) continue;
        }
	*/

	bool isTriggered = false;
	for( int k = 0; k < hlt_ntrig; k++ ) {
	  if( (hlt_trigName->at((unsigned int)k)) == trig ) {
	    if( hlt_trigFired[k] == 1 ) {
	      isTriggered = true;
	      break;
	    }
	  }
	}
        if( !isTriggered ) continue;
	//cout << "evt = " << evtnum << endl;

      if( _channel == "mumu" ) {
	vector<TLorentzVector> muPlus, muMinus;
        for( int j = 0; j < nMuon; j++ ) {
	  // for check muon charge
	  TLorentzVector tmpMu;
	  double mu_mass = 0.105658;
	  double rpx = muon_px[j];
	  double rpy = muon_py[j];
	  double rpz = muon_pz[j];
	  double rp = rpx*rpx + rpy*rpy + rpz*rpz;
	  double re = sqrt(mu_mass*mu_mass + rp);
	  tmpMu.SetPxPyPzE(rpx, rpy, rpz, re);

	  // GLB muon only
	  if( muon_type[j] != 0 && muon_type[j] != 1 ) continue;
	  // relative track isolation 
	  double _muIso = muon_trkiso[j];
	  double _muRelIso = (_muIso)/muon_cktpt[j];

   	  // good muon selection
	  if( muon_cktpt[j] > _ptCut 
              && fabs(muon_eta[j]) < _etaCut 
              && muon_trackerLayers[j] > 5
              && muon_pixelHits[j] > 0 
	      && _muRelIso < 0.10
              && fabs(muon_dxyVTX[j]) < 0.2 && fabs(muon_dzVTX[j]) < 0.5 
              && muon_muonHits[j] > 0 
	      && muon_nMatches[j] > 1 
	      && muon_cktptError[j] / muon_cktpt[j] < 0.3 ) {
		if( muon_charge[j] > 0 ) muPlus.push_back(tmpMu);
		if( muon_charge[j] < 0 ) muMinus.push_back(tmpMu);
	  }
	}

	// dilepton candidate
	TLorentzVector recoCand;
	double recoMass = 0;
	int iMuPlus = muPlus.size();
	int iMuMinus = muMinus.size();
	if( iMuPlus == 1 && iMuMinus == 1 ) {
	  recoCand = muPlus[0] + muMinus[0];
	}
	else if( iMuPlus+iMuMinus > 2 && iMuPlus >= 1 && iMuMinus >= 1 ) {
	  if( isDataPrint ) PrintAllDetails();
	}
	recoMass =  recoCand.M();
	if( recoMass > 0 ) _h1->Fill(recoMass);

	if( isDataPrint && recoMass > 0 ) {
	  //fout << "###########################" << endl;
	  if( !isMC ) {
	    //fout << _fname << " " << runnum << " " << lumisection << " " << evtnum << " " << endl;
	    //PrintAllDetails();
	  }
	}
      }
    }
    cout << "result = " << _h1->Integral() << " " << endl;
    return _h1;
}

bool getDYSpectrum( TString fname, double genMass )
{
  bool isPass = false;
  // muon
  if( fname == "DYMuMu_M20" ) 
    if( genMass > 200 ) isPass = true;
  if( fname == "DYMuMu_M200" ) 
    if( genMass > 400 ) isPass = true;
  if( fname == "DYMuMu_M400" ) 
    if( genMass > 500 ) isPass = true;
  if( fname == "DYMuMu_M500" ) 
    if( genMass > 700 ) isPass = true;
  if( fname == "DYMuMu_M700" ) 
    if( genMass > 800 ) isPass = true;
  if( fname == "DYMuMu_M800" ) 
    if( genMass > 1000 ) isPass = true;
  if( fname == "DYMuMu_M1000" ) 
    if( genMass > 1500 ) isPass = true;
  if( fname == "DYMuMu_M1500" ) 
    if( genMass > 2000 ) isPass = true;

  // electron
  if( fname == "DYEE_M20" ) 
    if( genMass > 200 ) isPass = true;
  if( fname == "DYEE_M200" ) 
    if( genMass > 400 ) isPass = true;
  if( fname == "DYEE_M400" ) 
    if( genMass > 500 ) isPass = true;
  if( fname == "DYEE_M500" ) 
    if( genMass > 700 ) isPass = true;
  if( fname == "DYEE_M700" ) 
    if( genMass > 800 ) isPass = true;
  if( fname == "DYEE_M800" ) 
    if( genMass > 1000 ) isPass = true;
  if( fname == "DYEE_M1000" ) 
    if( genMass > 1500 ) isPass = true;
  if( fname == "DYEE_M1500" ) 
    if( genMass > 2000 ) isPass = true;

  return isPass;
}


void PrintAllDetails( void )
{
  fout << "muons: " << endl;
  for( int j = 0; j < nMuon; j++ ) {
    fout << j << " = " 
    << muon_cktpt[j] << "(pt), " 
    << muon_eta[j] << "(eta), " 
    << muon_type[j] << "(type), "
    << muon_trackerLayers[j] << "(tkLayers), "
    << muon_pixelHits[j] << "(pHits), " 
    << muon_muonHits[j] << "(mHits), "
    << muon_nMatches[j] << "(nMatches), "
    << muon_cktptError[j] / muon_cktpt[j] << "(sigPt), "
    << (muon_trkiso[j])/muon_cktpt[j] << "(iso), " 
    << muon_dxyVTX[j] << "(dxy), "
    << muon_dzVTX[j] << "(dz), "
    << endl;
  }

  fout << "electrons: " << endl;
  //int _nMerged = 0;
  for( int j = 0; j < nelec; j++ ) {
    // check additional gsfTrack in a cone (dR < 0.15)
    /*
    bool isMoreGsfTrack = false;
    for( int l = 0; l < nTT; l++ ) {
      if( track_pT[l] < 30.0 ) continue;
      double dR = deltaR(elec_etaSC[j], elec_phiSC[j], track_eta[l], track_phi[l]);
      if( dR < 0.30 && elec_eoverp[j] > 0.7 ) isMoreGsfTrack = true;
    }
    if( isMoreGsfTrack ) _nMerged++;
    */
    fout << j << " = " 
    << elec_et[j] << "(et), "
    << elec_etaSC[j] << "(eta), "
    << elec_ecalDriven[j] << "(eDriven), "
    << elec_isoPtTrks[j] << "(isoPtTrk), "
    << elec_dPhiIn[j] << "(dPhiIn), "
    << elec_HoverE[j] << "(HoverE), "
    << elec_mHits[j] << "(mHits), "
    << elec_dEtaIn[j] << "(dEtaIn), "
    << elec_25over55[j] << "(25over55), "
    << elec_15over55[j] << "(15over55), "
    << elec_dxyVTX[j] << "(dxy), "
    << elec_sigmaIEtaIEta[j] << "(sigmaIetaIeta), "
    << elec_isoEMHADDepth1[j] << "(isoEMHADDepth1), "
    << endl;
  }
  fout << "###########################" << endl;
}

