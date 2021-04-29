# "Phys" package is to make ntuples and dimuon analysis to select events for Z peak region (60 < M_ll < 120 GeV)


# setup envirnment 
###this runs with the following recipe (with csh or tcsh based)


setenv SCRAM_ARCH slc6_amd64_gcc472


setenv VO_CMS_SW_DIR /cvmfs/cms.cern.ch

source $VO_CMS_SW_DIR/cmsset_default.csh

# setup working directory with cmssw release
mkdir test
cd test
cmsrel CMSSW_5_3_32
cd CMSSW_5_3_32/src/
cmsenv

# compile
git clone git@github.com:hdyoo/Phys.git ### replaced by downloading (SW)
cd Phys/DYntupleMaker/
scramv1 b

# produce ntuples
cd ntuples
replace input root files by yours in data_cfg.py and mc_cfg.py
cmsRun data_cfg.py >&log_data& ### if the muon server is busy, it takes 1-2 days for both data and MC
cmsRun mc_cfg.py >&log_mc&
top ## make sure your codes are running and finishes

# produce the dimuon invariant mass distribution of selected events: Z boson candidates
cd ../Mll/
vi chain_DYJets.C;; vi chain_SingleMu_Run2011B.C ### modify the root file by your ntuple produced by above recipe
root -b massSpec.C+ ### compile your macro in ROOT
root getControlPlots.C ### you should see outputs
