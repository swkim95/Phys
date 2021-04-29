#include "chain_SingleMu_Run2011B.C"
#include "chain_DYJets.C"

#include "SetupTree.h"
#include <TLorentzVector.h>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>

void LoadBranches( TChain* chains )
{
    chains->SetBranchStatus("*", 0);
    chains->SetBranchStatus("nTotal", 1);
    chains->SetBranchStatus("runNum", 1);
    chains->SetBranchStatus("lumiBlock", 1);
    chains->SetBranchStatus("evtNum", 1);
    chains->SetBranchStatus("nVertices", 1);
    chains->SetBranchStatus("pileUpReweight", 1);
    chains->SetBranchStatus("HLT_ntrig", 1);
    chains->SetBranchStatus("HLT_trigType", 1);
    chains->SetBranchStatus("HLT_trigFired", 1);
    chains->SetBranchStatus("HLT_trigName", 1);
    chains->SetBranchStatus("HLT_trigPS", 1);
    chains->SetBranchStatus("HLT_trigPt", 1);
    chains->SetBranchStatus("HLT_trigEta", 1);
    chains->SetBranchStatus("HLT_trigPhi", 1);

    chains->SetBranchStatus("nMuon", 1);
    chains->SetBranchStatus("Muon_pT", 1);
    chains->SetBranchStatus("Muon_eta", 1);
    chains->SetBranchStatus("Muon_charge", 1);
    chains->SetBranchStatus("Muon_muonType", 1);
    chains->SetBranchStatus("Muon_phi", 1);
    chains->SetBranchStatus("Muon_chi2dof", 1);
    chains->SetBranchStatus("Muon_nhits", 1);
    chains->SetBranchStatus("Muon_trackerHits", 1);
    chains->SetBranchStatus("Muon_trackerLayers", 1);
    chains->SetBranchStatus("Muon_pixelHits", 1);
    chains->SetBranchStatus("Muon_muonHits", 1);
    chains->SetBranchStatus("Muon_nMatches", 1);
    chains->SetBranchStatus("Muon_dxyVTX", 1);
    chains->SetBranchStatus("Muon_dzVTX", 1);
    chains->SetBranchStatus("Muon_dxyBS", 1);
    chains->SetBranchStatus("Muon_dzBS", 1);
    chains->SetBranchStatus("Muon_trkiso", 1);
    chains->SetBranchStatus("Muon_PfChargedHadronIsoR04", 1);
    chains->SetBranchStatus("Muon_PfNeutralHadronIsoR04", 1);
    chains->SetBranchStatus("Muon_PfGammaIsoR04", 1);
    chains->SetBranchStatus("Muon_Px", 1);
    chains->SetBranchStatus("Muon_Py", 1);
    chains->SetBranchStatus("Muon_Pz", 1);
    chains->SetBranchStatus("Muon_cktpT", 1);
    chains->SetBranchStatus("Muon_cktPx", 1);
    chains->SetBranchStatus("Muon_cktPy", 1);
    chains->SetBranchStatus("Muon_cktPz", 1);
    chains->SetBranchStatus("Muon_cktpTError", 1);
    chains->SetBranchStatus("Muon_dxycktVTX", 1);
    chains->SetBranchStatus("Muon_dzcktVTX", 1);
    chains->SetBranchStatus("CosAngle", 1);

    //chains->SetBranchStatus("Njets", 1);
    //chains->SetBranchStatus("JETpt", 1);
    //chains->SetBranchStatus("JETeta", 1);
    //chains->SetBranchStatus("JETphi", 1);

    chains->SetBranchAddress("nTotal", &nTotal);
    chains->SetBranchAddress("runNum", &runnum);
    chains->SetBranchAddress("evtNum", &evtnum);
    chains->SetBranchAddress("lumiBlock", &lumisection);
    chains->SetBranchAddress("nMuon", &nMuon);
    chains->SetBranchAddress("pileUpReweight", &PUweight);
    chains->SetBranchAddress("nVertices", &nvertices);

    chains->SetBranchAddress("HLT_ntrig", &hlt_ntrig);
    chains->SetBranchAddress("HLT_trigType", &hlt_trigType);
    chains->SetBranchAddress("HLT_trigName", &hlt_trigName);
    chains->SetBranchAddress("HLT_trigPS", &hlt_trigPS);
    chains->SetBranchAddress("HLT_trigFired", &hlt_trigFired);
    chains->SetBranchAddress("HLT_trigPt", &hlt_trigPt);
    chains->SetBranchAddress("HLT_trigEta", &hlt_trigEta);
    chains->SetBranchAddress("HLT_trigPhi", &hlt_trigPhi);

    chains->SetBranchAddress("Muon_pT", &muon_pt);
    chains->SetBranchAddress("Muon_eta", &muon_eta);
    chains->SetBranchAddress("Muon_charge", &muon_charge);
    chains->SetBranchAddress("Muon_muonType", &muon_type);
    chains->SetBranchAddress("Muon_phi", &muon_phi);
    chains->SetBranchAddress("Muon_chi2dof", &muon_chi2dof);
    chains->SetBranchAddress("Muon_nhits", &muon_nhits);
    chains->SetBranchAddress("Muon_trackerHits", &muon_trackerHits);
    chains->SetBranchAddress("Muon_trackerLayers", &muon_trackerLayers);
    chains->SetBranchAddress("Muon_pixelHits", &muon_pixelHits);
    chains->SetBranchAddress("Muon_muonHits", &muon_muonHits);
    chains->SetBranchAddress("Muon_nMatches", &muon_nMatches);
    chains->SetBranchAddress("Muon_dxyVTX", &muon_dxyVTX);
    chains->SetBranchAddress("Muon_dzVTX", &muon_dzVTX);
    chains->SetBranchAddress("Muon_dxyBS", &muon_dxyBS);
    chains->SetBranchAddress("Muon_dzBS", &muon_dzBS);
    chains->SetBranchAddress("Muon_trkiso", &muon_trkiso);
    chains->SetBranchAddress("Muon_PfChargedHadronIsoR04", &muon_PfChargedHadronIso);
    chains->SetBranchAddress("Muon_PfNeutralHadronIsoR04", &muon_PfNeutralHadronIso);
    chains->SetBranchAddress("Muon_PfGammaIsoR04", &muon_PfGammaIso);
    chains->SetBranchAddress("Muon_Px", &muon_px);
    chains->SetBranchAddress("Muon_Py", &muon_py);
    chains->SetBranchAddress("Muon_Pz", &muon_pz);

    chains->SetBranchAddress("Muon_cktpT", &muon_cktpt);
    chains->SetBranchAddress("Muon_cktPx", &muon_cktpx);
    chains->SetBranchAddress("Muon_cktPy", &muon_cktpy);
    chains->SetBranchAddress("Muon_cktPz", &muon_cktpz);
    chains->SetBranchAddress("Muon_cktpTError", &muon_cktptError);
    chains->SetBranchAddress("Muon_dxycktVTX", &muon_dxycktVTX);
    chains->SetBranchAddress("Muon_dzcktVTX", &muon_dzcktVTX);

    //chains->SetBranchAddress("Njets", &njets);
    //chains->SetBranchAddress("JETpt", &jet_pt);
    //chains->SetBranchAddress("JETeta", &jet_eta);
    //chains->SetBranchAddress("JETphi", &jet_phi);
}

void LoadBranchesMC( TChain* chains )
{
    LoadBranches(chains);

    chains->SetBranchStatus("GENnPair", 1);
    chains->SetBranchStatus("GENLepton_Px", 1);
    chains->SetBranchStatus("GENLepton_Py", 1);
    chains->SetBranchStatus("GENLepton_Pz", 1);
    chains->SetBranchStatus("GENLepton_pT", 1);
    chains->SetBranchStatus("GENLepton_eta", 1);
    chains->SetBranchStatus("GENLepton_phi", 1);
    chains->SetBranchStatus("GENLepton_mother", 1);
    chains->SetBranchStatus("GENLepton_charge", 1);
    chains->SetBranchStatus("GENLepton_status", 1);
    chains->SetBranchStatus("GENLepton_ID", 1);

    chains->SetBranchAddress("GENnPair", &gnpair);
    chains->SetBranchAddress("GENLepton_Px", &genLepton_px);
    chains->SetBranchAddress("GENLepton_Py", &genLepton_py);
    chains->SetBranchAddress("GENLepton_Pz", &genLepton_pz);
    chains->SetBranchAddress("GENLepton_mother", &genLepton_mother);
    chains->SetBranchAddress("GENLepton_pT", &genLepton_pt);
    chains->SetBranchAddress("GENLepton_eta", &genLepton_eta);
    chains->SetBranchAddress("GENLepton_phi", &genLepton_phi);
    chains->SetBranchAddress("GENLepton_charge", &genLepton_charge);
    chains->SetBranchAddress("GENLepton_status", &genLepton_status);
    chains->SetBranchAddress("GENLepton_ID", &genLepton_ID);
}

void SetupTree( TString sample, TChain* chains )
{
    if( sample == "SingleMu_Run2011B" ) {
        chain_SingleMu_Run2011B(chains);
        LoadBranches(chains);
    }
    if( sample == "DYJets" ) {
        chain_DYJets(chains);
        LoadBranchesMC(chains);
    }
}
