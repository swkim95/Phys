import FWCore.ParameterSet.Config as cms

#from SHarper.HEEPAnalyzer.HEEPSelectionCuts_cfi import *

DYntupleMaker = cms.EDAnalyzer("DYntupleMaker",
        processName = cms.untracked.string("HLT"),
        isMC = cms.untracked.bool(True),
        ApplyFilter = cms.untracked.bool(True),
        FilterType = cms.untracked.int32(0),
        Muon = cms.untracked.InputTag("selectedPatMuons"),
        IntLuminosity = cms.untracked.double(0.25),
        nHighQualLeptons = cms.untracked.int32(2),
        CrossSection = cms.untracked.double(1.0),
        FilterEfficiency = cms.untracked.double(1.0),
        TotalNevents = cms.untracked.double(1),
        #StoreGENFlag = cms.untracked.bool(False),
        #StoreJetFlag = cms.untracked.bool(False),
        #StoreTTFlag = cms.untracked.bool(True),
        #StoreElectronFlag = cms.untracked.bool(True),
)
