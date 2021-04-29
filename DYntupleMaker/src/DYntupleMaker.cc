//-------------------------------------------------
//
//   Class: DYntupleMaker
//
//   Description: Ntuple maker for DY analysis
//
//
//   Author:
//   H.D. Yoo           Purdue University
//
//--------------------------------------------------


#include "Phys/DYntupleMaker/interface/DYntupleMaker.h"
// system include files
#include <memory>

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/METReco/interface/PFMETFwd.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Vertexing.h"
//#include "DataFormats/PatCandidates/interface/TriggerPrimitive.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

// more info
#include "DataFormats/L1GlobalMuonTrigger/interface/L1MuGMTCand.h"
#include "DataFormats/GeometrySurface/interface/GloballyPositioned.h"
#include "DataFormats/Candidate/interface/CandMatchMap.h"
#include "Geometry/Records/interface/GlobalTrackingGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/GlobalTrackingGeometry.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"

// trigger
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"

#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/ParameterSet/interface/ParameterSetfwd.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"

#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidateFwd.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/VertexTools/interface/InvariantMassFromVertex.h"

#include "PhysicsTools/RecoUtils/interface/CandCommonVertexFitter.h"

#include "DataFormats/MuonReco/interface/MuonTimeExtra.h"
#include "DataFormats/MuonReco/interface/MuonTimeExtraMap.h"

#include "DataFormats/MuonReco/interface/MuonCosmicCompatibility.h"

#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"
#include "DataFormats/RecoCandidate/interface/IsoDepositFwd.h"
#include "DataFormats/RecoCandidate/interface/IsoDepositVetos.h"
#include "PhysicsTools/IsolationAlgos/interface/IsoDepositVetoFactory.h"

#include "DataFormats/RecoCandidate/interface/IsoDepositDirection.h"
#include "DataFormats/PatCandidates/interface/Isolation.h"
#include "DataFormats/Math/interface/deltaPhi.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "FWCore/Utilities/interface/RegexMatch.h"

// high pT lepton selection
#include "DataFormats/MuonReco/interface/MuonCocktails.h"

//for TransientTrackBuilder
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

#include "MagneticField/Engine/interface/MagneticField.h"

#include <iostream>

#include <TTree.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TFile.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TFrame.h>
#include <TMath.h>
#include <TF1.h>
#include <TLorentzVector.h>
#include <TVector3.h>
#include <boost/foreach.hpp>

//
// class decleration
//

using namespace std;
using namespace reco;
using namespace edm;
using namespace pat;
using namespace isodeposit;


// Constructor
DYntupleMaker::DYntupleMaker(const edm::ParameterSet& iConfig)
{

  nEvt = 0;

  processName 			    = iConfig.getUntrackedParameter<string>("processName", "HLT");
  theDebugLevel                     = iConfig.getUntrackedParameter<int>("DebugLevel", 0);

  // RECO
  theMuonLabel                      = iConfig.getUntrackedParameter<edm::InputTag>("Muon", edm::InputTag("selectedPatMuonsPFlow"));

  theCrossSection                   = iConfig.getUntrackedParameter<double>("CrossSection", 1.0);
  theFilterEfficiency               = iConfig.getUntrackedParameter<double>("FilterEfficiency", 1.0);
  theTotalNevents                   = iConfig.getUntrackedParameter<double>("TotalNevents", 1000);
  theIntegratedLumi                 = iConfig.getUntrackedParameter<double>("IntLuminosity", 100);
  isMC                              =  iConfig.getUntrackedParameter<bool>("isMC");
  theStoreHLTReportFlag             =  iConfig.getUntrackedParameter<bool>("StoreHLTReportFlag", true);

/*
  theStorePriVtxFlag                =  iConfig.getUntrackedParameter<bool>("StorePriVtxFlag", true);
  theStoreGENFlag                   =  iConfig.getUntrackedParameter<bool>("StoreGENFlag", true);
  theStoreTTFlag                    =  iConfig.getUntrackedParameter<bool>("StoreTTFlag", false);

  if( isMC ) {
    PileUpRD_ = iConfig.getParameter< std::vector<double> >("PileUpRD");
    PileUpRDMuonPhys_ = iConfig.getParameter< std::vector<double> >("PileUpRDMuonPhys");
    PileUpMC_ = iConfig.getParameter< std::vector<double> >("PileUpMC");
  }
*/
}


DYntupleMaker::~DYntupleMaker()
{
}


//
// member functions
//

// ------------ method called to for each event  ------------
void DYntupleMaker::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theTTBuilder);
    
   // initialize for ntuple variables
   nMuon = -1; 
   Nmuons = -1;

   nVertices = -1;
   PVtrackSize = -1;
   PVchi2 = -1;
   PVndof = -1;
   PVnormalizedChi2 = -1;
   PVx = -1000;
   PVy = -1000;
   PVz = -1000;
   PVprob = -1;

   GENnPair = -1;

   //PU reweighting
   PUweight = -1;

   // trigger object
   _HLT_ntrig = -1;

   // PU reweight
   pileUpReweightIn = pileUpReweight = 1.0;
   pileUpReweightPlus = pileUpReweightMinus = 0.0;
   pileUpReweightInMuonPhys = pileUpReweightMuonPhys = 1.0;
   pileUpReweightPlusMuonPhys = pileUpReweightMinusMuonPhys = 0.0;

   _HLT_trigName.clear();
   _HLT_trigPS.clear();
   CosAngle.clear();
   vtxTrkCkt1Pt.clear();
   vtxTrkCkt2Pt.clear();
   vtxTrkProb.clear();
   vtxTrkNdof.clear();
   vtxTrkChi2.clear();
   vtxTrkDiE1Pt.clear();
   vtxTrkDiE2Pt.clear();
   vtxTrkDiEProb.clear();
   vtxTrkDiENdof.clear();
   vtxTrkDiEChi2.clear();
   vtxTrkEMu1Pt.clear();
   vtxTrkEMu2Pt.clear();
   vtxTrkEMuProb.clear();
   vtxTrkEMuNdof.clear();
   vtxTrkEMuChi2.clear();
    
   for( int i = 0; i < MPSIZE; i++ ) {
       // Trigger
       _HLT_trigType[i] = -1;
       _HLT_trigFired[i] = -1;
       _HLT_trigPt[i] = _HLT_trigEta[i] = _HLT_trigPhi[i] = -100;

	//pf isolations
       Muon_PfChargedHadronIsoR05[i] = -1;
       Muon_PfNeutralHadronIsoR05[i] = -1;
       Muon_PfGammaIsoR05[i] = -1;
       Muon_PfChargedHadronIsoR04[i] = -1;
       Muon_PfNeutralHadronIsoR04[i] = -1;
       Muon_PfGammaIsoR04[i] = -1;
       Muon_PfChargedHadronIsoR03[i] = -1;
       Muon_PfNeutralHadronIsoR03[i] = -1;
       Muon_PfGammaIsoR03[i] = -1;

       // muon type
       Muon_muonType[i] = -1;
       isPFmuon[i] = 0;
       isGLBmuon[i] = 0;
       isTRKmuon[i] = 0;
       isSTAmuon[i] = 0;

       // trigger
       Muon_nTrig[i] = Muon_triggerObjectType[i] = Muon_filterName[i] = -1;

       // Muon kinematics
       Muon_phi[i] = Muon_eta[i] = Muon_cktpT[i] = Muon_pT[i] = -100;
       Muon_cktPx[i] = Muon_cktPy[i] = Muon_cktPz[i] = -999;
       Muon_cktpTError[i] = -999;
       Muon_Px[i] = Muon_Py[i] = Muon_Pz[i] = -100;
       Muon_trkiso[i] = Muon_hcaliso[i] = Muon_ecaliso[i] = Muon_chi2dof[i] = -100;
       Muon_trkisoR05[i] = Muon_hcalisoR05[i] = Muon_ecalisoR05[i] = -100;
       Muon_nChambers[i] = Muon_nMatches[i] = Muon_stationMask[i] = Muon_nSegments[i] =  -1;
       Muon_charge[i] = Muon_nhits[i] = -100;
       Muon_trackerHits[i] = Muon_pixelHits[i] = Muon_muonHits[i] = -1;
       Muon_trackerLayers[i] = -1;

       Muon_qoverp[i] = Muon_theta[i] = Muon_lambda[i] = -100;
       Muon_dxy[i] = Muon_d0[i] = Muon_dsz[i] = Muon_dz[i] = -100;
       Muon_vx[i] = Muon_vy[i] = Muon_vz[i] = -100;
       Muon_dxyBS[i] = Muon_dszBS[i] = Muon_dzBS[i] = -100;
       Muon_dxyVTX[i] = Muon_dszVTX[i] = Muon_dzVTX[i] = -100;
       Muon_dxycktVTX[i] = Muon_dszcktVTX[i] = Muon_dzcktVTX[i] = -100;

       //Various track informations
       //MuonBestTrack
       Muon_Best_pT[i] = -9999;
       Muon_Best_pTError[i] = -9999;
       Muon_Best_Px[i] = -9999;
       Muon_Best_Py[i] = -9999;
       Muon_Best_Pz[i] = -9999;
       Muon_Best_eta[i] = -9999;
       Muon_Best_phi[i] = -9999;
       //Inner Track
       Muon_Inner_pT[i] = -9999;
       Muon_Inner_pTError[i] = -9999;
       Muon_Inner_Px[i] = -9999;
       Muon_Inner_Py[i] = -9999;
       Muon_Inner_Pz[i] = -9999;
       Muon_Inner_eta[i] = -9999;
       Muon_Inner_phi[i] = -9999;
       //Outer Track
       Muon_Outer_pT[i] = -9999;
       Muon_Outer_pTError[i] = -9999;
       Muon_Outer_Px[i] = -9999;
       Muon_Outer_Py[i] = -9999;
       Muon_Outer_Pz[i] = -9999;
       Muon_Outer_eta[i] = -9999;
       Muon_Outer_phi[i] = -9999;
       //Global Track
       Muon_GLB_pT[i] = -9999;
       Muon_GLB_pTError[i] = -9999;
       Muon_GLB_Px[i] = -9999;
       Muon_GLB_Py[i] = -9999;
       Muon_GLB_Pz[i] = -9999;
       Muon_GLB_eta[i] = -9999;
       Muon_GLB_phi[i] = -9999;

       //tuneP MuonBestTrack
       Muon_TuneP_pT[i] = -9999;
       Muon_TuneP_pTError[i] = -9999;
       Muon_TuneP_Px[i] = -9999;
       Muon_TuneP_Py[i] = -9999;
       Muon_TuneP_Pz[i] = -9999;
       Muon_TuneP_eta[i] = -9999;
       Muon_TuneP_phi[i] = -9999;
   
       // GEN
       GENLepton_phi[i] = GENLepton_eta[i] = GENLepton_pT[i] = GENLepton_mother[i] = -100;
       GENLepton_Px[i] = GENLepton_Py[i] = GENLepton_Pz[i] = -100;
       GENLepton_charge[i] = GENLepton_status[i] = GENLepton_ID[i] = -100;
       GENLepton_isPrompt[i] = 0;
       GENLepton_isPromptFinalState[i] = 0;
       GENLepton_isTauDecayProduct[i] = 0;
       GENLepton_isPromptTauDecayProduct[i] = 0;
       GENLepton_isDirectPromptTauDecayProductFinalState[i] = 0;
       GENLepton_isHardProcess[i] = 0;
       GENLepton_isLastCopy[i] = 0;
       GENLepton_isLastCopyBeforeFSR[i] = 0;
       GENLepton_isPromptDecayed[i] = 0;
       GENLepton_isDecayedLeptonHadron[i] = 0;
       GENLepton_fromHardProcessBeforeFSR[i] = 0;
       GENLepton_fromHardProcessDecayed[i] = 0;
       GENLepton_fromHardProcessFinalState[i] = 0;
       GENEvt_weight = 0; //Weights for NLO generated events

   }

   nEvt++;
   // run number & event number
   runNum = iEvent.id().run();
   evtNum = iEvent.id().event();
   lumiBlock = iEvent.id().luminosityBlock();
   //const bool isRD = iEvent.isRealData();

   edm::Handle<double> weight_;
   iEvent.getByLabel("PUweight", weight_);

   if(weight_.isValid()){
     PUweight = *weight_;
   }else{
     PUweight = 1.0;
   }

   //get the geometry
   edm::ESHandle<GlobalTrackingGeometry> glbTrackingGeometry;
   iSetup.get<GlobalTrackingGeometryRecord>().get(glbTrackingGeometry);

   // // Pile-up reweight
   // if( isMC ) {
   //   edm::Handle<std::vector< PileupSummaryInfo > >  PupInfo;
   //   iEvent.getByLabel(edm::InputTag("addPileupInfo"), PupInfo);
   //   std::vector<PileupSummaryInfo>::const_iterator PVI;

   //   int npv = -1;
   //   int npvin = -1;
   //   for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {

   //     int BX = PVI->getBunchCrossing();

   //     if(BX == 0) {
   //       npvin = PVI->getPU_NumInteractions(); // in time only
   //       npv = PVI->getTrueNumInteractions(); // in and out of time
   //       continue;
   //     }
   //   }

   //   nPileUp = npv;
   //   pileUpReweightIn = LumiWeights_.weight( npvin );
   //   pileUpReweight = LumiWeights_.weight( npv );
   
   //   pileUpReweightPlus  = PShiftUp_.ShiftWeight( npv );
   //   pileUpReweightMinus = PShiftDown_.ShiftWeight( npv );

   //   pileUpReweightInMuonPhys = LumiWeightsMuonPhys_.weight( npvin );
   //   pileUpReweightMuonPhys = LumiWeightsMuonPhys_.weight( npv );
   
   //   pileUpReweightPlusMuonPhys  = PShiftUpMuonPhys_.ShiftWeight( npv );
   //   pileUpReweightMinusMuonPhys = PShiftDownMuonPhys_.ShiftWeight( npv );
   // }

   // call pat objects
   edm::Handle< pat::MuonCollection > muonHandle;
   iEvent.getByLabel(theMuonLabel, muonHandle);
   using reco::MuonCollection;
   MuonCollection::const_iterator imuon;

   // //vertex
   math::XYZPoint RefVtx;
   RefVtx.SetXYZ(0, 0, 0);

   edm::Handle<reco::VertexCollection> pvHandle;
   iEvent.getByLabel("offlinePrimaryVerticesWithBS", pvHandle);
   const reco::VertexCollection &vertices = *pvHandle.product();
   nVertices = pvHandle->size();
   for(reco::VertexCollection::const_iterator it=vertices.begin() ; it!=vertices.end() ; ++it) {
     RefVtx = it->position();
     break;
   }
   const reco::Vertex &vtx = pvHandle->front();


   //
   // DYTree 
   //

   // muons
   ESHandle<MagneticField> B;
   iSetup.get<IdealMagneticFieldRecord>().get(B);
          
   int _nMuon = 0;
   Nmuons = muonHandle->size();
   for( unsigned i = 0; i != muonHandle->size(); i++ ) {
        const pat::Muon imuon = muonHandle->at(i);
       
    if( imuon.isStandAloneMuon() )
       isSTAmuon[_nMuon] = 1;
       
    if( imuon.isGlobalMuon() )
       isGLBmuon[_nMuon] = 1;
       
    if( imuon.isTrackerMuon() )
       isTRKmuon[_nMuon] = 1;
    
    if( imuon.isPFMuon() )
       isPFmuon[_nMuon] = 1;
    
	if( imuon.isStandAloneMuon() ) {
	  if( imuon.isGlobalMuon() ) {
	    if( imuon.isTrackerMuon() ) Muon_muonType[_nMuon] = 0; // STA+GLB+TM
	    else Muon_muonType[_nMuon] = 1; // STA+GLB
	  }
	  else {
	    if( imuon.isTrackerMuon() ) Muon_muonType[_nMuon] = 2; // STA+TM
	    else Muon_muonType[_nMuon] = 3; // STA
	  }
	}
	else {
	  if( imuon.isTrackerMuon() ) Muon_muonType[_nMuon] = 4; // TM
	}
	if( Muon_muonType[_nMuon] == 3 ) continue;

        // bits 0-1-2-3 = DT stations 1-2-3-4
        // bits 4-5-6-7 = CSC stations 1-2-3-4
        int _segments = 0;
        for( int idet = 1; idet < 4; idet++ ) {
          // DT (1), CSC (2), RPC (3)
          for( int istation = 1; istation < 5; istation++ ) {
            // station 1, 2, 3, 4
            _segments += imuon.numberOfSegments(istation, idet);
          }
        }
        Muon_nSegments[_nMuon] = _segments;

        // reco track information
        reco::TrackRef trackerTrack = imuon.innerTrack();
        reco::TrackRef muonTrack    = imuon.outerTrack();
        reco::TrackRef glbTrack     = imuon.globalTrack();

	reco::TrackRef cktTrack = (muon::tevOptimized(imuon, 200, 17., 40., 0.25)).first;

	if( glbTrack.isNonnull() ) {
          Muon_chi2dof[_nMuon] = glbTrack->normalizedChi2();
          Muon_nhits[_nMuon] = glbTrack->numberOfValidHits();
          
	  const reco::HitPattern & glbhit = glbTrack->hitPattern();
          Muon_muonHits[_nMuon] = glbhit.numberOfValidMuonHits();

          Muon_qoverp[_nMuon] = glbTrack->qoverp();
          Muon_theta[_nMuon] = glbTrack->theta();
          Muon_lambda[_nMuon] = glbTrack->lambda();
          Muon_dxy[_nMuon] = glbTrack->dxy();
          Muon_d0[_nMuon] = glbTrack->d0();
          Muon_dsz[_nMuon] = glbTrack->dsz();
          Muon_dz[_nMuon] = glbTrack->dz();
          
	  Muon_vx[_nMuon] = glbTrack->vx();
	  Muon_vy[_nMuon] = glbTrack->vy();
	  Muon_vz[_nMuon] = glbTrack->vz();
	}
	else {
	  if( trackerTrack.isNonnull() ) {
            Muon_chi2dof[_nMuon] = trackerTrack->normalizedChi2();
            Muon_nhits[_nMuon] = trackerTrack->numberOfValidHits();
            
	    if( muonTrack.isNonnull() ) {
	      const reco::HitPattern & muonhit = muonTrack->hitPattern();
              Muon_muonHits[_nMuon] = muonhit.numberOfValidMuonHits();
	    }
	    else {
              Muon_muonHits[_nMuon] = 0;
	    }
  
            Muon_qoverp[_nMuon] = trackerTrack->qoverp();
            Muon_theta[_nMuon] = trackerTrack->theta();
            Muon_lambda[_nMuon] = trackerTrack->lambda();
            Muon_dxy[_nMuon] = trackerTrack->dxy();
            Muon_d0[_nMuon] = trackerTrack->d0();
            Muon_dsz[_nMuon] = trackerTrack->dsz();
            Muon_dz[_nMuon] = trackerTrack->dz();
          
	    Muon_vx[_nMuon] = trackerTrack->vx();
	    Muon_vy[_nMuon] = trackerTrack->vy();
	    Muon_vz[_nMuon] = trackerTrack->vz();
	  }
	}

	if( trackerTrack.isNonnull() ) {
          const reco::HitPattern & inhit = trackerTrack->hitPattern();

          Muon_trackerHits[_nMuon] = inhit.numberOfValidTrackerHits();
          Muon_pixelHits[_nMuon] = inhit.numberOfValidPixelHits();
          Muon_trackerLayers[_nMuon] = inhit.trackerLayersWithMeasurement();
	}

        if( !pvHandle->empty() && !pvHandle->front().isFake() ) {
          Muon_dxyVTX[_nMuon] = imuon.muonBestTrack()->dxy(vtx.position());
          Muon_dszVTX[_nMuon] = imuon.muonBestTrack()->dsz(vtx.position());
          Muon_dzVTX[_nMuon] = imuon.muonBestTrack()->dz(vtx.position());

          Muon_dxycktVTX[_nMuon] = cktTrack->dxy(vtx.position());
          Muon_dszcktVTX[_nMuon] = cktTrack->dsz(vtx.position());
          Muon_dzcktVTX[_nMuon] = cktTrack->dz(vtx.position());
        }

        // muon1 kinematics
        Muon_cktpT[_nMuon] = cktTrack->pt();
        Muon_cktPx[_nMuon] = cktTrack->px();
        Muon_cktPy[_nMuon] = cktTrack->py();
        Muon_cktPz[_nMuon] = cktTrack->pz();
        Muon_pT[_nMuon] = imuon.pt();
        Muon_Px[_nMuon] = imuon.px();
        Muon_Py[_nMuon] = imuon.py();
        Muon_Pz[_nMuon] = imuon.pz();
        Muon_eta[_nMuon] = imuon.eta();
        Muon_phi[_nMuon] = imuon.phi();

        //Various track informations
        //MuonBestTrack
        if( imuon.muonBestTrack().isNonnull() )
        {
          Muon_Best_pT[_nMuon] = imuon.muonBestTrack()->pt();
          Muon_Best_pTError[_nMuon] = imuon.muonBestTrack()->ptError();
          Muon_Best_Px[_nMuon] = imuon.muonBestTrack()->px();
          Muon_Best_Py[_nMuon] = imuon.muonBestTrack()->py();
          Muon_Best_Pz[_nMuon] = imuon.muonBestTrack()->pz();
          Muon_Best_eta[_nMuon] = imuon.muonBestTrack()->eta();
          Muon_Best_phi[_nMuon] = imuon.muonBestTrack()->phi();
        }

        //Inner Track
        if( imuon.innerTrack().isNonnull() )
        {
          Muon_Inner_pT[_nMuon] = imuon.innerTrack()->pt();
          Muon_Inner_pTError[_nMuon] = imuon.innerTrack()->ptError();
          Muon_Inner_Px[_nMuon] = imuon.innerTrack()->px();
          Muon_Inner_Py[_nMuon] = imuon.innerTrack()->py();
          Muon_Inner_Pz[_nMuon] = imuon.innerTrack()->pz();
          Muon_Inner_eta[_nMuon] = imuon.innerTrack()->eta();
          Muon_Inner_phi[_nMuon] = imuon.innerTrack()->phi();
        }

        //Outer Track
        if( imuon.outerTrack().isNonnull() )
        {
          Muon_Outer_pT[_nMuon] = imuon.outerTrack()->pt();
          Muon_Outer_pTError[_nMuon] = imuon.outerTrack()->ptError();
          Muon_Outer_Px[_nMuon] = imuon.outerTrack()->px();
          Muon_Outer_Py[_nMuon] = imuon.outerTrack()->py();
          Muon_Outer_Pz[_nMuon] = imuon.outerTrack()->pz();
          Muon_Outer_eta[_nMuon] = imuon.outerTrack()->eta();
          Muon_Outer_phi[_nMuon] = imuon.outerTrack()->phi();
        }

        //Global Track
        if( imuon.globalTrack().isNonnull() )
        {
          Muon_GLB_pT[_nMuon] = imuon.globalTrack()->pt();
          Muon_GLB_pTError[_nMuon] = imuon.globalTrack()->ptError();
          Muon_GLB_Px[_nMuon] = imuon.globalTrack()->px();
          Muon_GLB_Py[_nMuon] = imuon.globalTrack()->py();
          Muon_GLB_Pz[_nMuon] = imuon.globalTrack()->pz();
          Muon_GLB_eta[_nMuon] = imuon.globalTrack()->eta();
          Muon_GLB_phi[_nMuon] = imuon.globalTrack()->phi();
        }

        //tuneP MuonBestTrack
/*
        if( imuon.tunePMuonBestTrack().isNonnull() )
        {
          Muon_TuneP_pT[_nMuon] = imuon.tunePMuonBestTrack()->pt();
          Muon_TuneP_pTError[_nMuon] = imuon.tunePMuonBestTrack()->ptError();
          Muon_TuneP_Px[_nMuon] = imuon.tunePMuonBestTrack()->px();
          Muon_TuneP_Py[_nMuon] = imuon.tunePMuonBestTrack()->py();
          Muon_TuneP_Pz[_nMuon] = imuon.tunePMuonBestTrack()->pz();
          Muon_TuneP_eta[_nMuon] = imuon.tunePMuonBestTrack()->eta();
          Muon_TuneP_phi[_nMuon] = imuon.tunePMuonBestTrack()->phi();
        }
*/



        //ISOLATIONS GO HERE
        //detector based
        Muon_trkiso[_nMuon] = imuon.isolationR03().sumPt;
        Muon_hcaliso[_nMuon] = imuon.isolationR03().hadEt;
        Muon_ecaliso[_nMuon] = imuon.isolationR03().emEt;
        Muon_trkisoR05[_nMuon] = imuon.isolationR05().sumPt;
        Muon_hcalisoR05[_nMuon] = imuon.isolationR05().hadEt;
        Muon_ecalisoR05[_nMuon] = imuon.isolationR05().emEt; 
        //pf isolation 
        Muon_PfChargedHadronIsoR04[_nMuon] = imuon.pfIsolationR04().sumChargedHadronPt;
        Muon_PfNeutralHadronIsoR04[_nMuon] = imuon.pfIsolationR04().sumNeutralHadronEt;
        Muon_PfGammaIsoR04[_nMuon] = imuon.pfIsolationR04().sumPhotonEt;
        Muon_PfChargedHadronIsoR03[_nMuon] = imuon.pfIsolationR03().sumChargedHadronPt;
        Muon_PfNeutralHadronIsoR03[_nMuon] = imuon.pfIsolationR03().sumNeutralHadronEt;
        Muon_PfGammaIsoR03[_nMuon] = imuon.pfIsolationR03().sumPhotonEt;

        Muon_charge[_nMuon] = imuon.charge();
        Muon_nChambers[_nMuon] = imuon.numberOfChambers(); // # of chambers
        Muon_nMatches[_nMuon] = imuon.numberOfMatchedStations(); // # of chambers with matched segments
        Muon_stationMask[_nMuon] = imuon.stationMask(); // bit map of stations with matched segments

	Muon_cktpTError[_nMuon] = cktTrack->ptError();

	// filter for high pt Tight muon
	// https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonId#New_Version_recommended
	_nMuon++;

	// dimuon variables
        for( unsigned j = 0; j != muonHandle->size(); j++ ) {
          const pat::Muon imuon2 = muonHandle->at(j);
	  int index_type = -1;
	  if( imuon2.isStandAloneMuon() ) {
	    if( imuon2.isGlobalMuon() ) {
	      if( imuon2.isTrackerMuon() ) index_type = 0; // STA+GLB+TM
	      else index_type = 1; // STA+GLB
	    }
	    else {
	      if( imuon2.isTrackerMuon() ) index_type = 2; // STA+TM
	      else index_type = 3; // STA
	    }
	  }
	  else {
	    if( imuon2.isTrackerMuon() ) index_type = 4; // TM
  	  }
  	  if( index_type == 3 ) continue;

	  reco::TrackRef cktTrack2 = (muon::tevOptimized(imuon2, 200, 17., 40., 0.25)).first;
     
	  if( i <= j ) continue;
	 
	  if( cktTrack.isNonnull() && cktTrack2.isNonnull() ) {
            reco::TransientTrack muTransient1(cktTrack, B.product());
            reco::TransientTrack muTransient2(cktTrack2, B.product());

            vector<reco::TransientTrack> dimuonTracksTrk;
            dimuonTracksTrk.push_back(muTransient1);
            dimuonTracksTrk.push_back(muTransient2);
            KalmanVertexFitter KalmanFitterTrk(true);
            CachingVertex<5> vertexTrk;
            TransientVertex vtxtmpTrk;
            bool isVertexTrk = true;
            try {
              vertexTrk = KalmanFitterTrk.vertex(dimuonTracksTrk);
              vtxtmpTrk = KalmanFitterTrk.vertex(dimuonTracksTrk);
            }
            catch( exception & err ) {
              isVertexTrk = false;
            }
            if( isVertexTrk && vertexTrk.isValid() ) {
	      // inv. mass refit using the dimuon vtx
	      InvariantMassFromVertex imfvTrk;
	      static const double muon_mass = 0.1056583;
	      const CachingVertex<5>& vtxTrk = vertexTrk;
	      Measurement1D new_massTrk = imfvTrk.invariantMass(vtxTrk, muon_mass);
 
              vtxTrkCkt1Pt.push_back(cktTrack->pt());
              vtxTrkCkt2Pt.push_back(cktTrack2->pt());
 	      vtxTrkChi2.push_back(vtxTrk.totalChiSquared());
 	      vtxTrkNdof.push_back(vtxTrk.degreesOfFreedom());
 	      vtxTrkProb.push_back(TMath::Prob(vtxTrk.totalChiSquared(),(int)vtxTrk.degreesOfFreedom()));
	    }

            // cosmic variable
            double cosine = acos(-cktTrack->momentum().Dot(cktTrack2->momentum()/
                                   cktTrack->p()/cktTrack2->p()));
            CosAngle.push_back(cosine);
	  }
	}
   }
   nMuon = _nMuon;

   // fills
   if( theStoreHLTReportFlag ) hltReport(iEvent);
   //if( theStorePriVtxFlag ) fillPrimaryVertex(iEvent);
   //if( !isRD && theStoreGENFlag ) fillGENInfo(iEvent);
   //if( theStoreTTFlag ) fillTT(iEvent);
   DYTree->Fill();
}


// ------------ method called once each job just before starting event loop  ------------
void 
DYntupleMaker::beginJob()
{

  // if( isMC ) {
  //   // Pileup Reweight: 2012, Summer12_S10
  //   std::vector< float > _PUreweightRun2012 ;
  //   std::vector< float > _PUreweightRun2012MuonPhys ;
  //   std::vector< float > _MC2012;

  //   for( int i = 0; i < 100; ++i) {  
  //     _PUreweightRun2012.push_back((float)PileUpRD_[i]);
  //     _PUreweightRun2012MuonPhys.push_back((float)PileUpRDMuonPhys_[i]);
  //     _MC2012.push_back((float)PileUpMC_[i]);
  //   }

  //   LumiWeights_ = edm::LumiReWeighting(_MC2012, _PUreweightRun2012);
  //   PShiftDown_ = reweight::PoissonMeanShifter(-0.5);
  //   PShiftUp_ = reweight::PoissonMeanShifter(0.5);

  //   LumiWeightsMuonPhys_ = edm::LumiReWeighting(_MC2012, _PUreweightRun2012MuonPhys);
  //   PShiftDownMuonPhys_ = reweight::PoissonMeanShifter(-0.5);
  //   PShiftUpMuonPhys_ = reweight::PoissonMeanShifter(0.5);
  // }

  edm::Service<TFileService> fs;
  DYTree = fs->make<TTree>("DYTree","DYTree");

  // global event variables
  DYTree->Branch("nTotal",&nEvt,"nTotal/I");
  DYTree->Branch("runNum",&runNum,"runNum/I");
  DYTree->Branch("evtNum",&evtNum,"evtNum/I");
  DYTree->Branch("lumiBlock",&lumiBlock,"lumiBlock/I");
  DYTree->Branch("nMuon",&nMuon,"nMuon/I");
  DYTree->Branch("Nmuons",&Nmuons,"Nmuons/I");
  DYTree->Branch("PUweight",&PUweight,"PUweight/D");
  DYTree->Branch("nVertices",&nVertices,"nVertices/I");
 
  if(theStorePriVtxFlag){
    DYTree->Branch("PVtrackSize", &PVtrackSize,"PVtrackSize/I");
    DYTree->Branch("PVchi2", &PVchi2,"PVchi2/D");
    DYTree->Branch("PVndof", &PVndof,"PVndof/D");
    DYTree->Branch("PVnormalizedChi2", &PVnormalizedChi2,"PVnormalizedChi2/D");
    DYTree->Branch("PVx", &PVx,"PVx/D");
    DYTree->Branch("PVy", &PVy,"PVy/D");
    DYTree->Branch("PVz", &PVz,"PVz/D");
  }

  if (theStoreHLTReportFlag){
    DYTree->Branch("HLT_ntrig", &_HLT_ntrig,"HLT_ntrig/I");
    DYTree->Branch("HLT_trigType", &_HLT_trigType,"HLT_trigType[HLT_ntrig]/I");
    DYTree->Branch("HLT_trigFired", &_HLT_trigFired,"HLT_trigFired[HLT_ntrig]/I");
    DYTree->Branch("HLT_trigName", &_HLT_trigName);
    DYTree->Branch("HLT_trigPS", &_HLT_trigPS);
    DYTree->Branch("HLT_trigPt", &_HLT_trigPt,"HLT_trigPt[HLT_ntrig]/D");
    DYTree->Branch("HLT_trigEta", &_HLT_trigEta,"HLT_trigEta[HLT_ntrig]/D");
    DYTree->Branch("HLT_trigPhi", &_HLT_trigPhi,"HLT_trigPhi[HLT_ntrig]/D");
  }


  // object variables
  DYTree->Branch("Muon_muonType", &Muon_muonType,"Muon_muonType[nMuon]/I");
  DYTree->Branch("isPFmuon", &isPFmuon, "isPFmuon[nMuon]/I");
  DYTree->Branch("isGLBmuon", &isGLBmuon, "isGLBmuon[nMuon]/I");
  DYTree->Branch("isTRKmuon", &isTRKmuon, "isTRKmuon[nMuon]/I");
  DYTree->Branch("isSTAmuon", &isSTAmuon, "isSTAmuon[nMuon]/I");
  DYTree->Branch("Muon_nTrig", &Muon_nTrig,"Muon_nTrig[nMuon]/I");
  DYTree->Branch("Muon_triggerObjectType", &Muon_triggerObjectType,"Muon_triggerObjectType[nMuon]/I");
  DYTree->Branch("Muon_filterName", &Muon_filterName,"Muon_filterName[nMuon]/I");
  DYTree->Branch("Muon_phi", &Muon_phi,"Muon_phi[nMuon]/D");
  DYTree->Branch("Muon_eta", &Muon_eta,"Muon_eta[nMuon]/D");
  DYTree->Branch("Muon_pT", &Muon_pT,"Muon_pT[nMuon]/D");
  DYTree->Branch("Muon_cktpT", &Muon_cktpT,"Muon_cktpT[nMuon]/D");
  DYTree->Branch("Muon_cktPx", &Muon_cktPx,"Muon_cktPx[nMuon]/D");
  DYTree->Branch("Muon_cktPy", &Muon_cktPy,"Muon_cktPy[nMuon]/D");
  DYTree->Branch("Muon_cktPz", &Muon_cktPz,"Muon_cktPz[nMuon]/D");
  DYTree->Branch("Muon_cktpTError", &Muon_cktpTError,"Muon_cktpTError[nMuon]/D");
  DYTree->Branch("Muon_Px", &Muon_Px,"Muon_Px[nMuon]/D");
  DYTree->Branch("Muon_Py", &Muon_Py,"Muon_Py[nMuon]/D");
  DYTree->Branch("Muon_Pz", &Muon_Pz,"Muon_Pz[nMuon]/D");
  DYTree->Branch("Muon_trkiso", &Muon_trkiso,"Muon_trkiso[nMuon]/D");
  DYTree->Branch("Muon_hcaliso", &Muon_hcaliso,"Muon_hcaliso[nMuon]/D");
  DYTree->Branch("Muon_ecaliso", &Muon_ecaliso,"Muon_ecaliso[nMuon]/D");
  DYTree->Branch("Muon_trkisoR05", &Muon_trkisoR05,"Muon_trkisoR05[nMuon]/D");
  DYTree->Branch("Muon_hcalisoR05", &Muon_hcalisoR05,"Muon_hcalisoR05[nMuon]/D");
  DYTree->Branch("Muon_ecalisoR05", &Muon_ecalisoR05,"Muon_ecalisoR05[nMuon]/D");

  //Various track informations
  DYTree->Branch("Muon_Best_pT", &Muon_Best_pT, "Muon_Best_pT[nMuon]/D");
  DYTree->Branch("Muon_Best_pTError", &Muon_Best_pTError, "Muon_Best_pTError[nMuon]/D");
  DYTree->Branch("Muon_Best_Px", &Muon_Best_Px, "Muon_Best_Px[nMuon]/D");
  DYTree->Branch("Muon_Best_Py", &Muon_Best_Py, "Muon_Best_Py[nMuon]/D");
  DYTree->Branch("Muon_Best_Pz", &Muon_Best_Pz, "Muon_Best_Pz[nMuon]/D");
  DYTree->Branch("Muon_Best_eta", &Muon_Best_eta, "Muon_Best_eta[nMuon]/D");
  DYTree->Branch("Muon_Best_phi", &Muon_Best_phi, "Muon_Best_phi[nMuon]/D");

  DYTree->Branch("Muon_Inner_pT", &Muon_Inner_pT, "Muon_Inner_pT[nMuon]/D");
  DYTree->Branch("Muon_Inner_pTError", &Muon_Inner_pTError, "Muon_Inner_pTError[nMuon]/D");
  DYTree->Branch("Muon_Inner_Px", &Muon_Inner_Px, "Muon_Inner_Px[nMuon]/D");
  DYTree->Branch("Muon_Inner_Py", &Muon_Inner_Py, "Muon_Inner_Py[nMuon]/D");
  DYTree->Branch("Muon_Inner_Pz", &Muon_Inner_Pz, "Muon_Inner_Pz[nMuon]/D");
  DYTree->Branch("Muon_Inner_eta", &Muon_Inner_eta, "Muon_Inner_eta[nMuon]/D");
  DYTree->Branch("Muon_Inner_phi", &Muon_Inner_phi, "Muon_Inner_phi[nMuon]/D");

  DYTree->Branch("Muon_Outer_pT", &Muon_Outer_pT, "Muon_Outer_pT[nMuon]/D");
  DYTree->Branch("Muon_Outer_pTError", &Muon_Outer_pTError, "Muon_Outer_pTError[nMuon]/D");
  DYTree->Branch("Muon_Outer_Px", &Muon_Outer_Px, "Muon_Outer_Px[nMuon]/D");
  DYTree->Branch("Muon_Outer_Py", &Muon_Outer_Py, "Muon_Outer_Py[nMuon]/D");
  DYTree->Branch("Muon_Outer_Pz", &Muon_Outer_Pz, "Muon_Outer_Pz[nMuon]/D");
  DYTree->Branch("Muon_Outer_eta", &Muon_Outer_eta, "Muon_Outer_eta[nMuon]/D");
  DYTree->Branch("Muon_Outer_phi", &Muon_Outer_phi, "Muon_Outer_phi[nMuon]/D");

  DYTree->Branch("Muon_GLB_pT", &Muon_GLB_pT, "Muon_GLB_pT[nMuon]/D");
  DYTree->Branch("Muon_GLB_pTError", &Muon_GLB_pTError, "Muon_GLB_pTError[nMuon]/D");
  DYTree->Branch("Muon_GLB_Px", &Muon_GLB_Px, "Muon_GLB_Px[nMuon]/D");
  DYTree->Branch("Muon_GLB_Py", &Muon_GLB_Py, "Muon_GLB_Py[nMuon]/D");
  DYTree->Branch("Muon_GLB_Pz", &Muon_GLB_Pz, "Muon_GLB_Pz[nMuon]/D");
  DYTree->Branch("Muon_GLB_eta", &Muon_GLB_eta, "Muon_GLB_eta[nMuon]/D");
  DYTree->Branch("Muon_GLB_phi", &Muon_GLB_phi, "Muon_GLB_phi[nMuon]/D");

  DYTree->Branch("Muon_TuneP_pT", &Muon_TuneP_pT, "Muon_TuneP_pT[nMuon]/D");
  DYTree->Branch("Muon_TuneP_pTError", &Muon_TuneP_pTError, "Muon_TuneP_pTError[nMuon]/D");
  DYTree->Branch("Muon_TuneP_Px", &Muon_TuneP_Px, "Muon_TuneP_Px[nMuon]/D");
  DYTree->Branch("Muon_TuneP_Py", &Muon_TuneP_Py, "Muon_TuneP_Py[nMuon]/D");
  DYTree->Branch("Muon_TuneP_Pz", &Muon_TuneP_Pz, "Muon_TuneP_Pz[nMuon]/D");
  DYTree->Branch("Muon_TuneP_eta", &Muon_TuneP_eta, "Muon_TuneP_eta[nMuon]/D");
  DYTree->Branch("Muon_TuneP_phi", &Muon_TuneP_phi, "Muon_TuneP_phi[nMuon]/D");
         
  //pf iso
  DYTree->Branch("Muon_PfChargedHadronIsoR05", &Muon_PfChargedHadronIsoR05,"Muon_PfChargedHadronIsoR05[nMuon]/D");
  DYTree->Branch("Muon_PfNeutralHadronIsoR05", &Muon_PfNeutralHadronIsoR05,"Muon_PfNeutralHadronIsoR05[nMuon]/D");
  DYTree->Branch("Muon_PfGammaIsoR05", &Muon_PfGammaIsoR05,"Muon_PfGammaIsoR05[nMuon]/D");
  DYTree->Branch("Muon_PfChargedHadronIsoR04", &Muon_PfChargedHadronIsoR04,"Muon_PfChargedHadronIsoR04[nMuon]/D");
  DYTree->Branch("Muon_PfNeutralHadronIsoR04", &Muon_PfNeutralHadronIsoR04,"Muon_PfNeutralHadronIsoR04[nMuon]/D");
  DYTree->Branch("Muon_PfGammaIsoR04", &Muon_PfGammaIsoR04,"Muon_PfGammaIsoR04[nMuon]/D");
  DYTree->Branch("Muon_PfChargedHadronIsoR03", &Muon_PfChargedHadronIsoR03,"Muon_PfChargedHadronIsoR03[nMuon]/D");
  DYTree->Branch("Muon_PfNeutralHadronIsoR03", &Muon_PfNeutralHadronIsoR03,"Muon_PfNeutralHadronIsoR03[nMuon]/D");
  DYTree->Branch("Muon_PfGammaIsoR03", &Muon_PfGammaIsoR03,"Muon_PfGammaIsoR03[nMuon]/D");
  DYTree->Branch("Muon_charge", &Muon_charge,"Muon_charge[nMuon]/I");
  DYTree->Branch("Muon_nChambers", &Muon_nChambers,"Muon_nChambers[nMuon]/I");
  DYTree->Branch("Muon_nMatches", &Muon_nMatches,"Muon_nMatches[nMuon]/I");
  DYTree->Branch("Muon_stationMask", &Muon_stationMask,"Muon_stationMask[nMuon]/I");
  DYTree->Branch("Muon_nSegments", &Muon_nSegments,"Muon_nSegments[nMuon]/I");
  DYTree->Branch("Muon_chi2dof", &Muon_chi2dof,"Muon_chi2dof[nMuon]/D");
  DYTree->Branch("Muon_nhits", &Muon_nhits,"Muon_nhits[nMuon]/I");
  DYTree->Branch("Muon_trackerHits", &Muon_trackerHits,"Muon_trackerHits[nMuon]/I");
  DYTree->Branch("Muon_trackerLayers", &Muon_trackerLayers,"Muon_trackerLayers[nMuon]/I");
  DYTree->Branch("Muon_pixelHits", &Muon_pixelHits,"Muon_pixelHits[nMuon]/I");
  DYTree->Branch("Muon_muonHits", &Muon_muonHits,"Muon_muonHits[nMuon]/I");
  DYTree->Branch("Muon_qoverp", &Muon_qoverp,"Muon_qoverp[nMuon]/D");
  DYTree->Branch("Muon_theta", &Muon_theta,"Muon_theta[nMuon]/D");
  DYTree->Branch("Muon_lambda", &Muon_lambda,"Muon_lambda[nMuon]/D");
  DYTree->Branch("Muon_dxy", &Muon_dxy,"Muon_dxy[nMuon]/D");
  DYTree->Branch("Muon_d0", &Muon_d0,"Muon_d0[nMuon]/D");
  DYTree->Branch("Muon_dsz", &Muon_dsz,"Muon_dsz[nMuon]/D");
  DYTree->Branch("Muon_dz", &Muon_dz,"Muon_dz[nMuon]/D");
  DYTree->Branch("Muon_dxyBS", &Muon_dxyBS,"Muon_dxyBS[nMuon]/D");
  DYTree->Branch("Muon_dszBS", &Muon_dszBS,"Muon_dszBS[nMuon]/D");
  DYTree->Branch("Muon_dzBS", &Muon_dzBS,"Muon_dzBS[nMuon]/D");
  DYTree->Branch("Muon_dxyVTX", &Muon_dxyVTX,"Muon_dxyVTX[nMuon]/D");
  DYTree->Branch("Muon_dszVTX", &Muon_dszVTX,"Muon_dszVTX[nMuon]/D");
  DYTree->Branch("Muon_dzVTX", &Muon_dzVTX,"Muon_dzVTX[nMuon]/D");
  DYTree->Branch("Muon_dxycktVTX", &Muon_dxycktVTX,"Muon_dxycktVTX[nMuon]/D");
  DYTree->Branch("Muon_dszcktVTX", &Muon_dszcktVTX,"Muon_dszcktVTX[nMuon]/D");
  DYTree->Branch("Muon_dzcktVTX", &Muon_dzcktVTX,"Muon_dzcktVTX[nMuon]/D");
  DYTree->Branch("Muon_vx", &Muon_vx,"Muon_vx[nMuon]/D");
  DYTree->Branch("Muon_vy", &Muon_vy,"Muon_vy[nMuon]/D");
  DYTree->Branch("Muon_vz", &Muon_vz,"Muon_vz[nMuon]/D");
  DYTree->Branch("CosAngle", &CosAngle);
  DYTree->Branch("vtxTrkCkt1Pt", &vtxTrkCkt1Pt);
  DYTree->Branch("vtxTrkCkt2Pt", &vtxTrkCkt2Pt);
  DYTree->Branch("vtxTrkChi2", &vtxTrkChi2);
  DYTree->Branch("vtxTrkProb", &vtxTrkProb);
  DYTree->Branch("vtxTrkNdof", &vtxTrkNdof);
  DYTree->Branch("vtxTrkDiE1Pt", &vtxTrkDiE1Pt);
  DYTree->Branch("vtxTrkDiE2Pt", &vtxTrkDiE2Pt);
  DYTree->Branch("vtxTrkDiEChi2", &vtxTrkDiEChi2);
  DYTree->Branch("vtxTrkDiEProb", &vtxTrkDiEProb);
  DYTree->Branch("vtxTrkDiENdof", &vtxTrkDiENdof);
  DYTree->Branch("vtxTrkEMu1Pt", &vtxTrkEMu1Pt);
  DYTree->Branch("vtxTrkEMu2Pt", &vtxTrkEMu2Pt);
  DYTree->Branch("vtxTrkEMuChi2", &vtxTrkEMuChi2);
  DYTree->Branch("vtxTrkEMuProb", &vtxTrkEMuProb);
  DYTree->Branch("vtxTrkEMuNdof", &vtxTrkEMuNdof);

  // GEN info
  DYTree->Branch("GENnPair",&GENnPair,"GENnPair/I");
  DYTree->Branch("GENLepton_phi", &GENLepton_phi,"GENLepton_phi[GENnPair]/D");
  DYTree->Branch("GENLepton_eta", &GENLepton_eta,"GENLepton_eta[GENnPair]/D");
  DYTree->Branch("GENLepton_pT", &GENLepton_pT,"GENLepton_pT[GENnPair]/D");
  DYTree->Branch("GENLepton_Px", &GENLepton_Px,"GENLepton_Px[GENnPair]/D");
  DYTree->Branch("GENLepton_Py", &GENLepton_Py,"GENLepton_Py[GENnPair]/D");
  DYTree->Branch("GENLepton_Pz", &GENLepton_Pz,"GENLepton_Pz[GENnPair]/D");
  DYTree->Branch("GENLepton_mother", &GENLepton_mother,"GENLepton_mother[GENnPair]/D");
  DYTree->Branch("GENLepton_charge", &GENLepton_charge,"GENLepton_charge[GENnPair]/I");
  DYTree->Branch("GENLepton_status", &GENLepton_status,"GENLepton_status[GENnPair]/I");
  DYTree->Branch("GENLepton_ID", &GENLepton_ID,"GENLepton_ID[GENnPair]/I");
  DYTree->Branch("GENLepton_isPrompt", &GENLepton_isPrompt,"GENLepton_isPrompt[GENnPair]/I");
  DYTree->Branch("GENLepton_isPromptFinalState", &GENLepton_isPromptFinalState,"GENLepton_isPromptFinalState[GENnPair]/I");
  DYTree->Branch("GENLepton_isTauDecayProduct", &GENLepton_isTauDecayProduct,"GENLepton_isTauDecayProduct[GENnPair]/I");
  DYTree->Branch("GENLepton_isPromptTauDecayProduct", &GENLepton_isPromptTauDecayProduct,"GENLepton_isPromptTauDecayProduct[GENnPair]/I");
  DYTree->Branch("GENLepton_isDirectPromptTauDecayProductFinalState", &GENLepton_isDirectPromptTauDecayProductFinalState,"GENLepton_isDirectPromptTauDecayProductFinalState[GENnPair]/I");
  DYTree->Branch("GENLepton_isHardProcess",&GENLepton_isHardProcess,"GENLepton_isHardProcess[GENnPair]/I");
  DYTree->Branch("GENLepton_isLastCopy",&GENLepton_isLastCopy,"GENLepton_isLastCopy[GENnPair]/I");
  DYTree->Branch("GENLepton_isLastCopyBeforeFSR",&GENLepton_isLastCopyBeforeFSR,"GENLepton_isLastCopyBeforeFSR[GENnPair]/I");
  DYTree->Branch("GENLepton_isPromptDecayed",&GENLepton_isPromptDecayed,"GENLepton_isPromptDecayed[GENnPair]/I");
  DYTree->Branch("GENLepton_isDecayedLeptonHadron",&GENLepton_isDecayedLeptonHadron,"GENLepton_isDecayedLeptonHadron[GENnPair]/I");
  DYTree->Branch("GENLepton_fromHardProcessBeforeFSR",&GENLepton_fromHardProcessBeforeFSR,"GENLepton_fromHardProcessBeforeFSR[GENnPair]/I");
  DYTree->Branch("GENLepton_fromHardProcessDecayed",&GENLepton_fromHardProcessDecayed,"GENLepton_fromHardProcessDecayed[GENnPair]/I");
  DYTree->Branch("GENLepton_fromHardProcessFinalState",&GENLepton_fromHardProcessFinalState,"GENLepton_fromHardProcessFinalState[GENnPair]/I");
  DYTree->Branch("GENEvt_weight",&GENEvt_weight,"GENEvt_weight/D");



  // Pile-up Reweight
  DYTree->Branch("nPileUp",&nPileUp,"nPileUp/I");
  DYTree->Branch("pileUpReweightIn",&pileUpReweightIn,"pileUpReweightIn/D");
  DYTree->Branch("pileUpReweight",&pileUpReweight,"pileUpReweight/D");
  DYTree->Branch("pileUpReweightPlus",&pileUpReweightPlus,"pileUpReweightPlus/D");
  DYTree->Branch("pileUpReweightMinus",&pileUpReweightMinus,"pileUpReweightMinus/D");
  DYTree->Branch("pileUpReweightInMuonPhys",&pileUpReweightInMuonPhys,"pileUpReweightInMuonPhys/D");
  DYTree->Branch("pileUpReweightMuonPhys",&pileUpReweightMuonPhys,"pileUpReweightMuonPhys/D");
  DYTree->Branch("pileUpReweightPlusMuonPhys",&pileUpReweightPlusMuonPhys,"pileUpReweightPlusMuonPhys/D");
  DYTree->Branch("pileUpReweightMinusMuonPhys",&pileUpReweightMinusMuonPhys,"pileUpReweightMinusMuonPhys/D");

  if( theStoreTTFlag ) {
    DYTree->Branch("NTT", &NTT,"NTT/I");
    DYTree->Branch("TTrack_dxy", &TTrack_dxy,"TTrack_dxy[NTT]/D");
    DYTree->Branch("TTrack_dxyErr", &TTrack_dxyErr,"TTrack_dxyErr[NTT]/D");
    DYTree->Branch("TTrack_d0", &TTrack_d0,"TTrack_d0[NTT]/D");
    DYTree->Branch("TTrack_d0Err", &TTrack_d0Err,"TTrack_d0Err[NTT]/D");
    DYTree->Branch("TTrack_dsz", &TTrack_dsz,"TTrack_dsz[NTT]/D");
    DYTree->Branch("TTrack_dszErr", &TTrack_dszErr,"TTrack_dszErr[NTT]/D");
    DYTree->Branch("TTrack_dz", &TTrack_dz,"TTrack_dz[NTT]/D");
    DYTree->Branch("TTrack_dzErr", &TTrack_dzErr,"TTrack_dzErr[NTT]/D");
    DYTree->Branch("TTrack_dxyBS", &TTrack_dxyBS,"TTrack_dxyBS[NTT]/D");
    DYTree->Branch("TTrack_dszBS", &TTrack_dszBS,"TTrack_dszBS[NTT]/D");
    DYTree->Branch("TTrack_dzBS", &TTrack_dzBS,"TTrack_dzBS[NTT]/D");
    DYTree->Branch("TTrack_pT", &TTrack_pT,"TTrack_pT[NTT]/D");
    DYTree->Branch("TTrack_Px", &TTrack_Px,"TTrack_Px[NTT]/D");
    DYTree->Branch("TTrack_Py", &TTrack_Py,"TTrack_Py[NTT]/D");
    DYTree->Branch("TTrack_Pz", &TTrack_Pz,"TTrack_Pz[NTT]/D");
    DYTree->Branch("TTrack_eta", &TTrack_eta,"TTrack_eta[NTT]/D");
    DYTree->Branch("TTrack_phi", &TTrack_phi,"TTrack_phi[NTT]/D");
    DYTree->Branch("TTrack_charge", &TTrack_charge,"TTrack_charge[NTT]/D");
  }

}

void
DYntupleMaker::beginRun(const Run & iRun, const EventSetup & iSetup)
{
   const int nTrigName = 16;
   string trigs[nTrigName] = {
     "HLT_IsoMu24_eta2p1_v*",
     "HLT_Mu50_v*",
     "HLT_IsoMu24_eta2p1_v*",
     "HLT_Mu17_Mu8_v*", 
     "HLT_Mu17_TkMu8_v*",
     "HLT_Mu22_TkMu8_v*",
     "HLT_Mu24_v*",
     "HLT_DoubleEle33_CaloIdL_v*",
     "HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v*",
     "HLT_Ele27_WP80_v*",
     "HLT_IsoMu27_v*",
     "HLT_Mu45_eta2p1_v*",
     "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*",
     "HLT_IsoMu20_v*",
     "HLT_IsoMu24_IterTrk02_v*",
     "HLT_IsoMu20_eta2p1_IterTrk02_v*"
   };
   MuonHLT.clear();
   MuonHLTPS.clear();

   for( int i = 0; i < nTrigName; i++ ) {
     MuonHLT.push_back(trigs[i]);
   }

  int listRemoval[nTrigName] = {-1};
  int ntrigName = MuonHLT.size();
  bool changedConfig;
  if (!hltConfig_.init(iRun, iSetup, processName, changedConfig)) {
    LogError("HLTMuonVal") << "Initialization of HLTConfigProvider failed!!";
    return;
  }
  else {
    std::vector<std::string> triggerNames = hltConfig_.triggerNames();

    for( int itrigName = 0; itrigName < ntrigName; itrigName++ ) {
      listRemoval[itrigName] = 0;

      // check list of triggers
      //cout << "trigger = " << itrigName << " " << MuonHLT[itrigName] << endl;
      bool isMatched = false;
      for( size_t i = 0; i < triggerNames.size(); i++) {
        std::vector<std::string> moduleNames = hltConfig_.moduleLabels( triggerNames[i] );
        //cout << "Trigger Path: " << triggerNames[i] << endl;
	// find using wild card
        std::vector<std::vector<std::string>::const_iterator> matches = edm::regexMatch(triggerNames, MuonHLT[itrigName]);
        int _preScaleValue = hltConfig_.prescaleValue(0, triggerNames[i]);
	//cout << "prescale = " << _preScaleValue << endl;
        if( !matches.empty() ) {
	  BOOST_FOREACH(std::vector<std::string>::const_iterator match, matches) {
	    //cout << "trigger match = " << *match << endl;
	    if( *match == triggerNames[i] ) {
	      isMatched = true;
              MuonHLTPS.push_back(_preScaleValue);
	      /*
              cout << "Filter name: " << trigModuleNames[moduleNames.size()-2] << endl;
              for( size_t j = 0; j < moduleNames.size(); j++) {
                TString name = moduleNames[j];
                cout << "\t  Fliter Name: "<<moduleNames[j] << endl;
              }
	      */
	      int nsize = moduleNames.size();
	      if( nsize-2 >= 0 ) {
		//cout << "module names = " << moduleNames[nsize-2] << " " << moduleNames[nsize-3] << endl;
                trigModuleNames.push_back(moduleNames[nsize-2]);
                //cout << "Filter name: " << trigModuleNames[trigModuleNames.size()-1] << endl;
	        if( nsize-3 >= 0 ) {
                  trigModuleNames_preFil.push_back(moduleNames[nsize-3]);
	        }
		else {
                  trigModuleNames_preFil.push_back("");
		}
	      }
	      break;
	    }
	  }
	}
      }
      if( !isMatched ) listRemoval[itrigName] = 1;
    }
  }

   for( int i = 0; i < nTrigName; i++ ) {
     //cout << "prescale = " << MuonHLT[i] << " " << MuonHLTPS[i] << endl;
   }

  // remove unavailable triggers
  int itmp = 0;
  for( vector<string>::iterator iter = MuonHLT.begin(); iter != MuonHLT.end(); ) {
    if( listRemoval[itmp] > 0 ) 
      iter = MuonHLT.erase(iter);
    else 
      ++iter;
    itmp++;
  }
  ntrigName = MuonHLT.size();

  // trigger filters
  for( int itrig = 0; itrig < ntrigName; itrig++ ) {
    //cout << "Filter name: " << itrig << " " << MuonHLT[itrig] << " " << trigModuleNames[itrig] << " " << trigModuleNames_preFil[itrig] << endl;
  }
//cout << "beginRun: 2" << endl;

}

// ------------ method called once each job just after ending the event loop  ------------
void 
DYntupleMaker::endJob() {
  std::cout <<"++++++++++++++++++++++++++++++++++++++" << std::endl;
  std::cout <<"analyzed " << nEvt << " events: " << std::endl;
  std::cout <<"++++++++++++++++++++++++++++++++++++++" << std::endl;

}

///////////////////////////////////////////////////////////////
// makes hlt report and fills it to the ntuple
///////////////////////////////////////////////////////////////
void DYntupleMaker::hltReport(const edm::Event &iEvent) {
   //
   // trigger
   //

   int ntrigName = MuonHLT.size();
   for( int itrig = 0; itrig < ntrigName; itrig++ ) {
   }

   // read the whole HLT trigger lists fired in an event
   bool *trigFired = new bool[ntrigName];
   for( int i = 0; i < ntrigName; i++ ) trigFired[i] = false;
   Handle<TriggerResults> trigResult;
   iEvent.getByLabel(edm::InputTag("TriggerResults","",processName), trigResult);
   if( !trigResult.failedToGet() ) {
       int ntrigs = trigResult->size();
       const edm::TriggerNames trigName = iEvent.triggerNames(*trigResult);
       for( int itrig = 0; itrig != ntrigs; ++itrig) {
	   //cout << "trigName = " << trigName.triggerName(itrig) << " " << itrig << endl;
	   for( int itrigName = 0; itrigName < ntrigName; itrigName++ ) {
             std::vector<std::vector<std::string>::const_iterator> matches = edm::regexMatch(trigName.triggerNames(), MuonHLT[itrigName]);
             if( !matches.empty() ) {
	       BOOST_FOREACH(std::vector<std::string>::const_iterator match, matches) {
		 //cout << "trigger match = " << *match << endl;
	         if( trigName.triggerIndex(*match) >= (unsigned int)ntrigs ) continue;
	         if( trigResult->accept(trigName.triggerIndex(*match)) ) trigFired[itrigName] = true;
  	         else break;
                 //cout << "trigger fire = " << trigFired[itrigName] << endl;
	       }
	     }
	   }
       }
   }

   edm::Handle< trigger::TriggerEvent > triggerObject;
   iEvent.getByLabel(edm::InputTag("hltTriggerSummaryAOD","",processName), triggerObject);
   const trigger::TriggerObjectCollection & toc(triggerObject->getObjects());
   int ntrigTot = 0;
   //cout << "size filter = " << triggerObject->sizeFilters() << endl;
   for( size_t k = 0; k < triggerObject->sizeFilters(); ++ k ) {
       std::string fullname = triggerObject->filterTag(k).encode();
       std::string filterName;
       //cout << "filterName = " << filterName << endl;
       size_t m = fullname.find_first_of(':');
       if( m != std::string::npos ) {
	 filterName = fullname.substr(0, m);
       }
       else {
	 filterName = fullname;
       }

       if( &toc != 0 ) {
	 const trigger::Keys & it = triggerObject->filterKeys(k);
	 for( trigger::Keys::const_iterator ky = it.begin(); ky !=it.end(); ++ky ) {
	   double hlt_pt = toc[*ky].pt();
	   //double hlt_eta = toc[*ky].eta();
	   //cout << "hlt kine = " << hlt_pt << " " << hlt_eta << endl;

	   for( int itf = 0; itf < ntrigName; itf++ ) {
	     string names = "";
             //cout << "filterName = " << k << " " << filterName << " " << trigModuleNames[itf] << " " << trigModuleNames_preFil[itf] << " " << _HLT_trigName.size() << endl;
	       if( filterName == trigModuleNames[itf] ) { 
		 //((itf == index_hlt_mu17_mu8 || itf == index_hlt_mu17_tkmu8 || itf == index_hlt_mu22_tkmu8 || itf == index_hlt_mu22_tkmu22) && filterName == trigModuleNames_preFil[itf]) ) {
	         //cout << "filter name = " << filterName << " " << ntrigTot << " " << itf << " " << hlt_pt << " " << toc[*ky].eta() << " " << toc[*ky].phi() << endl;
		 names = MuonHLT[itf];
		 int _ps = MuonHLTPS[itf];
	         _HLT_trigType[ntrigTot] = itf;
	         _HLT_trigFired[ntrigTot] = trigFired[itf];
	         _HLT_trigPt[ntrigTot] = hlt_pt;
	         _HLT_trigEta[ntrigTot] = toc[*ky].eta();
	         _HLT_trigPhi[ntrigTot] = toc[*ky].phi();
		 _HLT_trigName.push_back(names);
		 _HLT_trigPS.push_back(_ps);
	         ntrigTot++;
	       }
	   }
	 }
       }
   }
   _HLT_ntrig = ntrigTot;
   
   for( int i = 0; i < _HLT_ntrig; i++ ) {
     //cout << "trig = " << i << " " << _HLT_trigType[i] << " " << _HLT_trigPt[i] << " " << _HLT_trigEta[i] << " " << _HLT_trigPhi[i] << " " << _HLT_trigName[i] << " " << _HLT_trigPS[i] <<  endl;
   }
}

//////////////////////////////////////////////////////////////
// Get Primary vertex info
///////////////////////////////////////////////////////////////
void DYntupleMaker::fillPrimaryVertex(const edm::Event &iEvent) {

   edm::Handle<reco::VertexCollection> pvHandle;
   iEvent.getByLabel("offlinePrimaryVertices", pvHandle);
   const reco::VertexCollection vtx = *(pvHandle.product());

   if( vtx.size() > 2 && theDebugLevel > 0) cout << "Reconstructed "<< vtx.size() << " vertices" << endl;
   if (vtx.size() > 0 ){
       PVtrackSize = vtx.front().tracksSize();
       PVchi2 = vtx.front().chi2();
       PVndof = vtx.front().ndof();
       PVnormalizedChi2 = vtx.front().normalizedChi2();
       PVx = vtx.front().x();
       PVy = vtx.front().y();
       PVz = vtx.front().z();
       PVprob = TMath::Prob(PVchi2,(int)PVndof) ;
   }
}

//////////////////////////////////////////////////////////////
// Get GEN info
///////////////////////////////////////////////////////////////
/*
void DYntupleMaker::fillGENInfo(const edm::Event &iEvent) {
  edm::Handle <reco::GenParticleCollection> particles;
  iEvent.getByLabel("genParticles", particles);

  int _GennPair = 0;
  for( size_t ipar = 0; ipar < particles->size(); ipar++ ) {
    const reco::GenParticle &parCand = (*particles)[ipar];
    if( abs(parCand.pdgId()) == 13 || abs(parCand.pdgId()) == 11 || abs(parCand.pdgId()) == 15 ) {
      GENLepton_ID[_GennPair] = parCand.pdgId(); 
      GENLepton_pT[_GennPair] = parCand.pt(); 
      GENLepton_Px[_GennPair] = parCand.px();
      GENLepton_Py[_GennPair] = parCand.py();
      GENLepton_Pz[_GennPair] = parCand.pz();
      GENLepton_eta[_GennPair] = parCand.eta();
      GENLepton_phi[_GennPair] = parCand.phi();
      GENLepton_charge[_GennPair] = parCand.charge();
      GENLepton_status[_GennPair] = parCand.status();
      GENLepton_mother[_GennPair] = parCand.mother(0)->pdgId();

      //Flags (Ref: https://indico.cern.ch/event/402279/contribution/5/attachments/805964/1104514/mcaod-Jun17-2015.pdf)
      GENLepton_isPrompt[_GennPair] = parCand.statusFlags().isPrompt(); //not from hadron, muon or tau decay
      GENLepton_isPromptFinalState[_GennPair] = parCand.isPromptFinalState(); //isPrompt && final state (status==1)
      GENLepton_isTauDecayProduct[_GennPair] = parCand.statusFlags().isTauDecayProduct(); //is directly or indirectly from a tau decay
      GENLepton_isPromptTauDecayProduct[_GennPair] = parCand.statusFlags().isPromptTauDecayProduct(); //is directly or indirectly from a tau decay, where the tau did not come from a hadron decay
      GENLepton_isDirectPromptTauDecayProductFinalState[_GennPair] = parCand.isDirectPromptTauDecayProductFinalState(); // is the direct decay product from a tau decay (ie no intermediate hadron), where the tau did not come from a hadron decay && final state
      GENLepton_isHardProcess[_GennPair] = parCand.isHardProcess();
      GENLepton_isLastCopy[_GennPair] = parCand.isLastCopy();
      GENLepton_isLastCopyBeforeFSR[_GennPair] = parCand.isLastCopyBeforeFSR();
      GENLepton_isPromptDecayed[_GennPair] = parCand.isPromptDecayed();
      GENLepton_isDecayedLeptonHadron[_GennPair] = parCand.statusFlags().isDecayedLeptonHadron();
      GENLepton_fromHardProcessBeforeFSR[_GennPair] = parCand.fromHardProcessBeforeFSR();
      GENLepton_fromHardProcessDecayed[_GennPair] = parCand.fromHardProcessDecayed();
      GENLepton_fromHardProcessFinalState[_GennPair] = parCand.fromHardProcessFinalState();

      _GennPair++;
    }
  }
  GENnPair = _GennPair;

  edm::Handle<GenEventInfoProduct> genEvtInfo;
  iEvent.getByLabel("generator", genEvtInfo);
  GENEvt_weight = genEvtInfo->weight();
  // std::cout << "Evt Weight: " << GENEvt_weight<< std::endl;
  }
}
*/



//////////////////////////////////////////////////////////////
// Get Tracker track info (all single tracks)
///////////////////////////////////////////////////////////////

void DYntupleMaker::fillTT(const edm::Event &iEvent) {

/*
   edm::Handle<edm::View<reco::Track> > trackHandle;
   iEvent.getByLabel(theTrackLabel, trackHandle);

   edm::Handle<reco::BeamSpot> beamSpotHandle;
   iEvent.getByLabel(theBeamSpot, beamSpotHandle);
   reco::BeamSpot beamSpot = (*beamSpotHandle);

   edm::Handle<reco::VertexCollection> _pvHandle;
   iEvent.getByLabel("offlinePrimaryVerticesWithBS", _pvHandle);
   const reco::Vertex &vtx = _pvHandle->front();

   // to store gsftracks close to electrons
   edm::Handle<reco::GsfTrackCollection> gsfTracks; 
   iEvent.getByLabel("electronGsfTracks", gsfTracks);

   int _nTT = 0;
   //for(edm::View<reco::Track>::const_iterator iTT = trackHandle->begin(); iTT != trackHandle->end(); ++iTT) {
   for(unsigned igsf = 0; igsf < gsfTracks->size(); igsf++ ) {
     GsfTrackRef iTT(gsfTracks, igsf);
     //if( iTT->pt() < 1.0 || iTT->pt() > 100000 ) continue;
     bool _isMatch = false;
     for( int i = 0; i < Nelectrons; i++ ) {
       double dpT = fabs(iTT->pt() - Electron_gsfpT[i]);
       double dR = deltaR(iTT->eta(), iTT->phi(), Electron_gsfEta[i], Electron_gsfPhi[i]);
       //cout << "elec = " << i << " " << Electron_gsfpT[i] << " " << Electron_gsfEta[i] << " " << Electron_gsfPhi[i] << " " << dR << endl;
       if( dR < 0.001 && dpT < 1.0 ) _isMatch = true;
     }
     if( _isMatch ) continue;

     TTrack_dxy[_nTT] = iTT->dxy(vtx.position());
     TTrack_dxyErr[_nTT] = iTT->dxyError();
     TTrack_d0[_nTT] = iTT->d0();
     TTrack_d0Err[_nTT] = iTT->d0Error(); 
     TTrack_dsz[_nTT] = iTT->dsz(vtx.position());
     TTrack_dszErr[_nTT] = iTT->dszError();
     TTrack_dz[_nTT] = iTT->dz(vtx.position());
     TTrack_dzErr[_nTT] = iTT->dzError();
     TTrack_dxyBS[_nTT] = iTT->dxy(beamSpot.position());
     TTrack_dszBS[_nTT] = iTT->dsz(beamSpot.position());
     TTrack_dzBS[_nTT] = iTT->dz(beamSpot.position());
     TTrack_pT[_nTT] = iTT->pt();
     TTrack_Px[_nTT] = iTT->px();
     TTrack_Py[_nTT] = iTT->py();
     TTrack_Pz[_nTT] = iTT->pz();
     TTrack_eta[_nTT] = iTT->eta();
     TTrack_phi[_nTT] = iTT->phi();
     TTrack_charge[_nTT] = iTT->charge();
     _nTT++;
   }
   NTT = _nTT;
*/

}

//define this as a plug-in
DEFINE_FWK_MODULE(DYntupleMaker);


