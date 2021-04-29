#ifndef DYntupleMaker_H
#define DYntupleMaker_H

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"
//#include "SimDataFormats/HepMCProduct/interface/HepMCProduct.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/GeometrySurface/interface/GloballyPositioned.h"

#include "Geometry/Records/interface/GlobalTrackingGeometryRecord.h"
#include "DataFormats/Common/interface/RefToBaseVector.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/Candidate/interface/NamedCompositeCandidateFwd.h"
#include "DataFormats/Candidate/interface/NamedCompositeCandidate.h"
#include "DataFormats/Candidate/interface/CompositeCandidateFwd.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "DataFormats/Common/interface/View.h"

#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidateFwd.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "PhysicsTools/RecoUtils/interface/CandCommonVertexFitter.h"

#include "DataFormats/MuonReco/interface/MuonCosmicCompatibility.h"

#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

// PU reweight
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"

#include <iostream>

#include <TFile.h>
#include <TTree.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TFrame.h>
#include <TMath.h>


//
// class decleration
//

using namespace std;
namespace reco { class CandCommonVertexFitterBase; class VertexCompositeCandidate; class CandCommonVertexFitter; }

class DYntupleMaker : public edm::EDAnalyzer {
   public:
      explicit DYntupleMaker(const edm::ParameterSet&);
      ~DYntupleMaker();

   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      virtual void beginRun(const edm::Run &, const edm::EventSetup & );

      virtual void fillPrimaryVertex(const edm::Event &iEvent);  // fill primary vertex information
      virtual void hltReport(const edm::Event &iEvent);          // fill list of triggers fired in an event
      //virtual void fillGENInfo(const edm::Event &iEvent);            // fill MET information
      virtual void fillTT(const edm::Event&);

      bool reorder(double &a, double &b) {
	return a > b;
      }
      
      int theDebugLevel;                   // 0 no prints, 1 some, 2 lots
      std::string processName;
      std::string theElectronID;

      HLTConfigProvider hltConfig_;

      edm::InputTag theElectronLabel;
      edm::InputTag theElectronCalibLabel;
      edm::InputTag theMuonLabel;
      edm::InputTag theMETLabel;
      edm::InputTag pfMetCollection_;
      edm::InputTag theJetLabel;
      edm::InputTag theTimeTags;
      edm::InputTag theBeamSpot;
      edm::InputTag theTrackLabel;

      edm::InputTag theConversionsInputTag;

      double theCrossSection;
      double theFilterEfficiency;
      double theTotalNevents;
      double theIntegratedLumi;
      int theNHighQualLeptons;
      
      bool theStorePriVtxFlag;                // Yes or No to store primary vertex
      bool theStoreHLTReportFlag;             // Yes or No to store HLT reuslts (list of triggers fired)
      bool theStoreGENFlag;
      bool theStoreTTFlag;
      bool isMC;  //turn gen on and off

      edm::ESHandle<TransientTrackBuilder> theTTBuilder;

      std::vector<std::string > MuonHLT;
      std::vector<int > MuonHLTPS;
      std::vector<std::string > trigModuleNames;
      std::vector<std::string > trigModuleNames_preFil;

      // Pile-up Reweight
      edm::LumiReWeighting LumiWeights_;
      reweight::PoissonMeanShifter PShiftUp_;
      reweight::PoissonMeanShifter PShiftDown_;
      edm::LumiReWeighting LumiWeightsMuonPhys_;
      reweight::PoissonMeanShifter PShiftUpMuonPhys_;
      reweight::PoissonMeanShifter PShiftDownMuonPhys_;

      std::vector<double> PileUpRD_;
      std::vector<double> PileUpRDMuonPhys_;
      std::vector<double> PileUpMC_;

      unsigned int nPileUp;
      double pileUpReweightIn;
      double pileUpReweight;
      double pileUpReweightPlus;
      double pileUpReweightMinus;

      double pileUpReweightInMuonPhys;
      double pileUpReweightMuonPhys;
      double pileUpReweightPlusMuonPhys;
      double pileUpReweightMinusMuonPhys;

      TTree *DYTree;

      int nEvt;

      //FILL a tree
      static const int MPSIZE = 2000;

      // Invariant Mass distribution of SS(OS) di-muon events GG (2 global)
      // GlbTree
      int runNum;
      int evtNum;
      int lumiBlock;
      int nMuon;
      double PUweight;
      int Nmuons;

      // PV
      int nVertices;
      int PVtrackSize;
      double PVchi2;
      double PVndof;
      double PVnormalizedChi2;
      double PVx;
      double PVy;
      double PVz;
      double PVprob;

      // trigger object
      int _HLT_ntrig;
      int _HLT_trigType[MPSIZE];
      int _HLT_trigFired[MPSIZE];
      std::vector<std::string> _HLT_trigName;
      std::vector<int> _HLT_trigPS;
      double _HLT_trigPt[MPSIZE];
      double _HLT_trigEta[MPSIZE];
      double _HLT_trigPhi[MPSIZE];

      // Pat Muon
      //pf isolations
      double Muon_PfChargedHadronIsoR05[MPSIZE];
      double Muon_PfNeutralHadronIsoR05[MPSIZE];
      double Muon_PfGammaIsoR05[MPSIZE];
      double Muon_PfChargedHadronIsoR04[MPSIZE];
      double Muon_PfNeutralHadronIsoR04[MPSIZE];
      double Muon_PfGammaIsoR04[MPSIZE];
      double Muon_PfChargedHadronIsoR03[MPSIZE];
      double Muon_PfNeutralHadronIsoR03[MPSIZE];
      double Muon_PfGammaIsoR03[MPSIZE];

      int isPFmuon[MPSIZE];
      int isGLBmuon[MPSIZE];
      int isTRKmuon[MPSIZE];
      int isSTAmuon[MPSIZE];
      int Muon_muonType[MPSIZE];
      int Muon_nTrig[MPSIZE];
      int Muon_triggerObjectType[MPSIZE];
      int Muon_filterName[MPSIZE];
      double Muon_phi[MPSIZE];
      double Muon_eta[MPSIZE];
      double Muon_pT[MPSIZE];
      double Muon_cktpT[MPSIZE];
      double Muon_cktPx[MPSIZE];
      double Muon_cktPy[MPSIZE];
      double Muon_cktPz[MPSIZE];
      double Muon_cktpTError[MPSIZE];
      double Muon_Px[MPSIZE];
      double Muon_Py[MPSIZE];
      double Muon_Pz[MPSIZE];
      double Muon_sumtrkpt[MPSIZE];
      double Muon_trkiso[MPSIZE];
      double Muon_hcaliso[MPSIZE];
      double Muon_ecaliso[MPSIZE];
      double Muon_trkisoR05[MPSIZE];
      double Muon_hcalisoR05[MPSIZE];
      double Muon_ecalisoR05[MPSIZE];
      int Muon_charge[MPSIZE];
      int Muon_nChambers[MPSIZE];
      int Muon_nMatches[MPSIZE];
      int Muon_stationMask[MPSIZE];
      int Muon_nSegments[MPSIZE];
      double Muon_chi2dof[MPSIZE];
      int Muon_nhits[MPSIZE];
      int Muon_trackerHits[MPSIZE];
      int Muon_pixelHits[MPSIZE];
      int Muon_muonHits[MPSIZE];
      int Muon_trackerLayers[MPSIZE];
      double Muon_qoverp[MPSIZE];
      double Muon_theta[MPSIZE];
      double Muon_lambda[MPSIZE];
      double Muon_dxy[MPSIZE];
      double Muon_d0[MPSIZE];
      double Muon_dsz[MPSIZE];
      double Muon_dz[MPSIZE];
      double Muon_dxyBS[MPSIZE];
      double Muon_dzBS[MPSIZE];
      double Muon_dszBS[MPSIZE];
      double Muon_dxyVTX[MPSIZE];
      double Muon_dzVTX[MPSIZE];
      double Muon_dszVTX[MPSIZE];
      double Muon_dxycktVTX[MPSIZE];
      double Muon_dzcktVTX[MPSIZE];
      double Muon_dszcktVTX[MPSIZE];
      double Muon_vx[MPSIZE];
      double Muon_vy[MPSIZE];
      double Muon_vz[MPSIZE];
      std::vector<double> CosAngle;
      std::vector<double> vtxTrkChi2;
      std::vector<double> vtxTrkProb;
      std::vector<double> vtxTrkNdof;
      std::vector<double> vtxTrkCkt1Pt;
      std::vector<double> vtxTrkCkt2Pt;
      std::vector<double> vtxTrkDiEChi2;
      std::vector<double> vtxTrkDiEProb;
      std::vector<double> vtxTrkDiENdof;
      std::vector<double> vtxTrkDiE1Pt;
      std::vector<double> vtxTrkDiE2Pt;
      std::vector<double> vtxTrkEMuChi2;
      std::vector<double> vtxTrkEMuProb;
      std::vector<double> vtxTrkEMuNdof;
      std::vector<double> vtxTrkEMu1Pt;
      std::vector<double> vtxTrkEMu2Pt;

      //Various track informations
      //MuonBestTrack
      double Muon_Best_pT[MPSIZE];
      double Muon_Best_pTError[MPSIZE];
      double Muon_Best_Px[MPSIZE];
      double Muon_Best_Py[MPSIZE];
      double Muon_Best_Pz[MPSIZE];
      double Muon_Best_eta[MPSIZE];
      double Muon_Best_phi[MPSIZE];
      //Inner Track
      double Muon_Inner_pT[MPSIZE];
      double Muon_Inner_pTError[MPSIZE];
      double Muon_Inner_Px[MPSIZE];
      double Muon_Inner_Py[MPSIZE];
      double Muon_Inner_Pz[MPSIZE];
      double Muon_Inner_eta[MPSIZE];
      double Muon_Inner_phi[MPSIZE];
      //Outer Track
      double Muon_Outer_pT[MPSIZE];
      double Muon_Outer_pTError[MPSIZE];
      double Muon_Outer_Px[MPSIZE];
      double Muon_Outer_Py[MPSIZE];
      double Muon_Outer_Pz[MPSIZE];
      double Muon_Outer_eta[MPSIZE];
      double Muon_Outer_phi[MPSIZE];
      //Global Track
      double Muon_GLB_pT[MPSIZE];
      double Muon_GLB_pTError[MPSIZE];
      double Muon_GLB_Px[MPSIZE];
      double Muon_GLB_Py[MPSIZE];
      double Muon_GLB_Pz[MPSIZE];
      double Muon_GLB_eta[MPSIZE];
      double Muon_GLB_phi[MPSIZE];

      //tuneP MuonBestTrack
      double Muon_TuneP_pT[MPSIZE];
      double Muon_TuneP_pTError[MPSIZE];
      double Muon_TuneP_Px[MPSIZE];
      double Muon_TuneP_Py[MPSIZE];
      double Muon_TuneP_Pz[MPSIZE];
      double Muon_TuneP_eta[MPSIZE];
      double Muon_TuneP_phi[MPSIZE];

      // GEN
      int GENnPair;
      double GENLepton_phi[MPSIZE];
      double GENLepton_eta[MPSIZE];
      double GENLepton_pT[MPSIZE];
      double GENLepton_Px[MPSIZE];
      double GENLepton_Py[MPSIZE];
      double GENLepton_Pz[MPSIZE];
      double GENLepton_mother[MPSIZE];
      int GENLepton_charge[MPSIZE];
      int GENLepton_status[MPSIZE];
      int GENLepton_ID[MPSIZE];
      int GENLepton_isPrompt[MPSIZE];
      int GENLepton_isPromptFinalState[MPSIZE];
      int GENLepton_isTauDecayProduct[MPSIZE];
      int GENLepton_isPromptTauDecayProduct[MPSIZE];
      int GENLepton_isDirectPromptTauDecayProductFinalState[MPSIZE];
      int GENLepton_isHardProcess[MPSIZE];
      int GENLepton_isLastCopy[MPSIZE];
      int GENLepton_isLastCopyBeforeFSR[MPSIZE];
      int GENLepton_isPromptDecayed[MPSIZE];
      int GENLepton_isDecayedLeptonHadron[MPSIZE];
      int GENLepton_fromHardProcessBeforeFSR[MPSIZE];
      int GENLepton_fromHardProcessDecayed[MPSIZE];
      int GENLepton_fromHardProcessFinalState[MPSIZE];
      double GENEvt_weight;

      int NTT;
      double TTrack_dxy[MPSIZE];
      double TTrack_dxyErr[MPSIZE];
      double TTrack_d0[MPSIZE];
      double TTrack_d0Err[MPSIZE];
      double TTrack_dsz[MPSIZE];
      double TTrack_dszErr[MPSIZE];
      double TTrack_dz[MPSIZE];
      double TTrack_dzErr[MPSIZE];
      double TTrack_dxyBS[MPSIZE];
      double TTrack_dszBS[MPSIZE];
      double TTrack_dzBS[MPSIZE];
      double TTrack_pT[MPSIZE];
      double TTrack_Px[MPSIZE];
      double TTrack_Py[MPSIZE];
      double TTrack_Pz[MPSIZE];
      double TTrack_eta[MPSIZE];
      double TTrack_phi[MPSIZE];
      double TTrack_charge[MPSIZE];
};
#endif
