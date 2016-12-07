#!/usr/bin/env python

# standard imports
import ROOT
import sys
import os
import copy
import random
import subprocess
import datetime
import shutil

from array import array
from operator import mul
from math import sqrt, atan2, sin, cos

# RootTools
from RootTools.core.standard import *

# User specific
import StopsDilepton.tools.user as user

# Tools for systematics
from StopsDilepton.tools.mt2Calculator import mt2Calculator
mt2Calc = mt2Calculator()  #smth smarter possible?
from StopsDilepton.tools.helpers import closestOSDLMassToMZ, checkRootFile, writeObjToFile, m3, deltaR, bestDRMatchInCollection
from StopsDilepton.tools.addJERScaling import addJERScaling
from StopsDilepton.tools.objectSelection import getMuons, getElectrons, muonSelector, eleSelector, getGoodLeptons, getGoodAndOtherLeptons,  getGoodBJets, getGoodJets, isBJet, jetId, isBJet, getGoodPhotons, getGenPartsAll, multiIsoWPInt
from StopsDilepton.tools.overlapRemovalTTG import getTTGJetsEventType
from StopsDilepton.tools.getGenBoson import getGenZ, getGenPhoton

from StopsDilepton.tools.triggerEfficiency import triggerEfficiency
triggerEff_withBackup = triggerEfficiency(with_backup_triggers = True)
triggerEff            = triggerEfficiency(with_backup_triggers = False)

from StopsDilepton.tools.leptonHIPEfficiency import leptonHIPEfficiency
leptonHIPSF = leptonHIPEfficiency()

#MC tools
from StopsDilepton.tools.mcTools import pdgToName, GenSearch, B_mesons, D_mesons, B_mesons_abs, D_mesons_abs
genSearch = GenSearch()

# central configuration
targetLumi = 1000 #pb-1 Which lumi to normalize to

def get_parser():
    ''' Argument parser for post-processing module.
    '''
    import argparse
    argParser = argparse.ArgumentParser(description = "Argument parser for cmgPostProcessing")

    argParser.add_argument('--logLevel',
        action='store',
        nargs='?',
        choices=['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'TRACE', 'NOTSET'],
        default='INFO',
        help="Log level for logging"
        )

    argParser.add_argument('--overwrite',
        action='store_true',
        help="Overwrite existing output files, bool flag set to True  if used")

    argParser.add_argument('--samples',
        action='store',
        nargs='*',
        type=str,
#        default=['MuonEG_Run2015D_16Dec'],
        default=['TTJets'],
        help="List of samples to be post-processed, given as CMG component name"
        )

    argParser.add_argument('--triggerSelection',
        action='store',
        nargs='?',
        type=str,
        default=None,
        choices=['mumu', 'ee', 'mue', 'mu_for_mumu', 'e_for_ee', 'mu_for_mue', 'e_for_mue'],
        help="Trigger selection?"
        )

    argParser.add_argument('--eventsPerJob',
        action='store',
        nargs='?',
        type=int,
        default=300000,
        help="Maximum number of events per job (Approximate!)."
        )

    argParser.add_argument('--nJobs',
        action='store',
        nargs='?',
        type=int,
        default=1,
        help="Maximum number of simultaneous jobs."
        )
    argParser.add_argument('--job',
        action='store',
        nargs='*',
        type=int,
        default=[],
        help="Run only jobs i"
        )

    argParser.add_argument('--minNJobs',
        action='store',
        nargs='?',
        type=int,
        default=1,
        help="Minimum number of simultaneous jobs."
        )

    argParser.add_argument('--dataDir',
        action='store',
        nargs='?',
        type=str,
        default=user.cmg_directory,
        help="Name of the directory where the input data is stored (for samples read from Heppy)."
        )

    argParser.add_argument('--targetDir',
        action='store',
        nargs='?',
        type=str,
        default=user.data_output_directory,
        help="Name of the directory the post-processed files will be saved"
        )

    argParser.add_argument('--processingEra',
        action='store',
        nargs='?',
        type=str,
        default='postProcessed_80X_v16',
        help="Name of the processing era"
        )

    argParser.add_argument('--skim',
        action='store',
        nargs='?',
        type=str,
        default='dilepTiny',
        help="Skim conditions to be applied for post-processing"
        )

    argParser.add_argument('--LHEHTCut',
        action='store',
        nargs='?',
        type=int,
        default=-1,
        help="LHE cut."
        )

    argParser.add_argument('--keepForwardJets',
        action='store_true',
        help="Is T2tt signal?"
        )

    argParser.add_argument('--small',
        action='store_true',
        help="Run the file on a small sample (for test purpose), bool flag set to True if used",
        #default = True
        )

    argParser.add_argument('--T2tt',
        action='store_true',
        help="Is T2tt signal?"
        )

    argParser.add_argument('--TTDM',
        action='store_true',
        help="Is TTDM signal?"
        )

    argParser.add_argument('--fastSim',
        action='store_true',
        help="FastSim?"
        )

    argParser.add_argument('--skipGenLepMatching',
        action='store_true',
        help="skip matched genleps??"
        )

    argParser.add_argument('--keepLHEWeights',
        action='store_true',
        help="Keep LHEWeights?"
        )

    argParser.add_argument('--checkTTGJetsOverlap',
        action='store_true',
        default=True,
        help="Keep TTGJetsEventType which can be used to clean TTG events from TTJets samples"
        )

    argParser.add_argument('--skipSystematicVariations',
        action='store_true',
        help="Don't calulcate BTag, JES and JER variations."
        )

    argParser.add_argument('--noTopPtReweighting',
        action='store_true',
        help="Skip top pt reweighting.")

    return argParser

options = get_parser().parse_args()

# Logging
import StopsDilepton.tools.logger as logger
logger = logger.get_logger(options.logLevel, logFile ='/tmp/%s_%s.txt'%(options.skim, '_'.join(options.samples) ) )
import RootTools.core.logger as logger_rt
logger_rt = logger_rt.get_logger(options.logLevel, logFile = None )

# flags (I think string searching is slow, so let's not do it in the filler function)
isDiLep     =   options.skim.lower().startswith('dilep')
isTriLep     =   options.skim.lower().startswith('trilep')
isSingleLep =   options.skim.lower().startswith('singlelep')
isTiny      =   options.skim.lower().count('tiny') 
isSmall      =   options.skim.lower().count('small')
isInclusive  = options.skim.lower().count('inclusive') 
isVeryLoose =  'veryloose' in options.skim.lower()
isVeryLoosePt10 =  'veryloosept10' in options.skim.lower()
isLoose     =  'loose' in options.skim.lower() and not isVeryLoose
isJet250    = 'jet250' in options.skim.lower()

# Skim condition
skimConds = []
if isDiLep:
    skimConds.append( "Sum$(LepGood_pt>20&&abs(LepGood_eta)<2.5) + Sum$(LepOther_pt>20&&abs(LepOther_eta)<2.5)>=2" )
if isTriLep:
    raise NotImplementedError
    # skimConds.append( "Sum$(LepGood_pt>20&&abs(LepGood_eta)&&LepGood_miniRelIso<0.4) + Sum$(LepOther_pt>20&&abs(LepOther_eta)<2.5&&LepGood_miniRelIso<0.4)>=2 && Sum$(LepOther_pt>10&&abs(LepOther_eta)<2.5)+Sum$(LepGood_pt>10&&abs(LepGood_eta)<2.5)>=3" )
elif isSingleLep:
    skimConds.append( "Sum$(LepGood_pt>20&&abs(LepGood_eta)<2.5) + Sum$(LepOther_pt>20&&abs(LepOther_eta)<2.5)>=1" )
elif isJet250:
    skimConds.append( "Sum$(Jet_pt>250) +  Sum$(DiscJet_pt>250) + Sum$(JetFailId_pt>250) + Sum$(gamma_pt>250) > 0" )

if isInclusive:
    skimConds = []

#Samples: Load samples
maxN = 2 if options.small else None
from StopsDilepton.samples.helpers import fromHeppySample
if options.T2tt:
    samples = [ fromHeppySample(s, data_path = options.dataDir, maxN = maxN) for s in options.samples ]
    from StopsDilepton.samples.helpers import getT2ttSignalWeight
    logger.info( "T2tt signal samples to be processed: %s", ",".join(s.name for s in samples) )
    # FIXME I'm forcing ==1 signal sample because I don't have a good idea how to construct a sample name from the complicated T2tt_x_y_z_... names
    assert len(samples)==1, "Can only process one T2tt sample at a time."
    samples[0].files = samples[0].files[:maxN]
    logger.debug( "Fetching signal weights..." )
    signalWeight = getT2ttSignalWeight( samples[0], lumi = targetLumi )
    logger.debug("Done fetching signal weights.")
#elif options.TTDM:
#    samples = [ fromHeppySample(s, data_path = "/scratch/rschoefbeck/cmgTuples/80X_0l_TTDM/", \
#                    module = "CMGTools.StopsDilepton.TTbarDMJets_signals_RunIISpring16MiniAODv2",  
#                    maxN = maxN)\
#                for s in options.samples ]
else:
    samples = [ fromHeppySample(s, data_path = options.dataDir, maxN = maxN) for s in options.samples ]

if len(samples)==0:
    logger.info( "No samples found. Was looking for %s. Exiting" % options.samples )
    sys.exit(-1)

isData = False not in [s.isData for s in samples]
isMC   =  True not in [s.isData for s in samples]

if options.T2tt:
    xSection = None
    # special filet for bad jets in FastSim: https://twiki.cern.ch/twiki/bin/viewauth/CMS/SUSRecommendationsICHEP16#Cleaning_up_of_fastsim_jets_from
    skimConds.append( "Sum$(JetFailId_pt>30&&abs(JetFailId_eta)<2.5&&JetFailId_mcPt==0&&JetFailId_chHEF<0.1)+Sum$(Jet_pt>30&&abs(Jet_eta)<2.5&&Jet_mcPt==0&&Jet_chHEF<0.1)==0" )
else:
    # Check that all samples which are concatenated have the same x-section.
    assert isData or len(set([s.heppy.xSection for s in samples]))==1, "Not all samples have the same xSection: %s !"%(",".join([s.name for s in samples]))
    assert isMC or len(samples)==1, "Don't concatenate data samples"

    xSection = samples[0].heppy.xSection if isMC else None


if isData and options.triggerSelection is not None:
    if options.triggerSelection == 'mumu':
        skimConds.append( "(HLT_mumuIso||HLT_mumuNoiso)" )
    elif options.triggerSelection == 'ee':
        skimConds.append( "(HLT_ee_DZ||HLT_ee_33||HLT_ee_33_MW)" )
    elif options.triggerSelection == 'mue':
        skimConds.append( "(HLT_mue||HLT_mu30e30)" )
    elif options.triggerSelection == 'mu_for_mumu':
        skimConds.append( "HLT_SingleMu_noniso&&(!(HLT_mumuIso||HLT_mumuNoiso))" )
    elif options.triggerSelection == 'e_for_ee':
        skimConds.append( "HLT_SingleEle_noniso&&(!(HLT_ee_DZ||HLT_ee_33||HLT_ee_33_MW))" )
    elif options.triggerSelection == 'e_for_mue':
        skimConds.append( "HLT_SingleEle_noniso && (!(HLT_mue||HLT_mu30e30))" )
    elif options.triggerSelection == 'mu_for_mue':
        skimConds.append( "HLT_SingleMu_noniso && (!(HLT_mue||HLT_mu30e30)) && (!HLT_SingleEle_noniso)" )
    else:
        raise ValueError( "Don't know about triggerSelection %s"%options.triggerSelection )
    sample_name_postFix = "_Trig_"+options.triggerSelection
    logger.info( "Added trigger selection %s and postFix %s", options.triggerSelection, sample_name_postFix )
else:
    sample_name_postFix = ""

#Samples: combine if more than one
if len(samples)>1:
    sample_name =  samples[0].name+"_comb"
    logger.info( "Combining samples %s to %s.", ",".join(s.name for s in samples), sample_name )
    sample = Sample.combine(sample_name, samples, maxN = maxN)
    # Clean up
    for s in samples:
        sample.clear()
elif len(samples)==1:
    sample = samples[0]
    sample.name+=sample_name_postFix
else:
    raise ValueError( "Need at least one sample. Got %r",samples )

if isMC:
    from StopsDilepton.tools.puReweighting import getReweightingFunction
    # nTrueIntReweighting
    nTrueInt27fb_puRW        = getReweightingFunction(data="PU_2016_27000_XSecCentral", mc="Spring16")
    nTrueInt27fb_puRWDown    = getReweightingFunction(data="PU_2016_27000_XSecDown", mc="Spring16")
    nTrueInt27fb_puRWUp      = getReweightingFunction(data="PU_2016_27000_XSecUp", mc="Spring16")
    nTrueInt12fb_puRW        = getReweightingFunction(data="PU_2016_12000_XSecCentral", mc="Spring16")
    nTrueInt12fb_puRWDown    = getReweightingFunction(data="PU_2016_12000_XSecDown", mc="Spring16")
    nTrueInt12fb_puRWUp      = getReweightingFunction(data="PU_2016_12000_XSecUp", mc="Spring16")
    nTrueInt_puRW        = getReweightingFunction(data="PU_2016_5300_XSecCentral", mc="Spring16")
    nTrueInt_puRWDown    = getReweightingFunction(data="PU_2016_5300_XSecDown", mc="Spring16")
    nTrueInt_puRWUp      = getReweightingFunction(data="PU_2016_5300_XSecUp", mc="Spring16")
    ## 2016 NVTX reweighting
    #from StopsDilepton.tools.puReweighting import getNVTXReweightingFunction
    #nvtx_reweight_central = getNVTXReweightingFunction(key = 'rw', filename = "dilepton_allZ_isOS_5300pb.pkl")
    #nvtx_reweight_up      = getNVTXReweightingFunction(key = 'up', filename = "dilepton_allZ_isOS_5300pb.pkl")
    #nvtx_reweight_down    = getNVTXReweightingFunction(key = 'down', filename = "dilepton_allZ_isOS_5300pb.pkl")
    #nvtx_reweight_vup      = getNVTXReweightingFunction(key = 'vup', filename = "dilepton_allZ_isOS_5300pb.pkl")
    #nvtx_reweight_vdown    = getNVTXReweightingFunction(key = 'vdown', filename = "dilepton_allZ_isOS_5300pb.pkl")
        
# top pt reweighting
from StopsDilepton.tools.topPtReweighting import getUnscaledTopPairPtReweightungFunction, getTopPtDrawString, getTopPtsForReweighting
# Decision based on sample name -> whether TTJets or TTLep is in the sample name
isTT = sample.name.startswith("TTJets") or sample.name.startswith("TTLep") or sample.name.startswith("TT_pow")
doTopPtReweighting = isTT and not options.noTopPtReweighting
if doTopPtReweighting:
    logger.info( "Sample will have top pt reweighting." )
    topPtReweightingFunc = getUnscaledTopPairPtReweightungFunction(selection = "dilep")
    # Compute x-sec scale factor on unweighted events
    selectionString = "&&".join(skimConds)
    topScaleF = sample.getYieldFromDraw( selectionString = selectionString, weightString = getTopPtDrawString(selection = "dilep"))
    topScaleF = topScaleF['val']/float(sample.chain.GetEntries(selectionString))
    logger.info( "Found topScaleF %f", topScaleF )
else:
    topScaleF = 1
    logger.info( "Sample will NOT have top pt reweighting. topScaleF=%f",topScaleF )

if options.fastSim:
   from StopsDilepton.tools.leptonFastSimSF import leptonFastSimSF as leptonFastSimSF_
   leptonFastSimSF = leptonFastSimSF_()

from StopsDilepton.tools.leptonSF import leptonSF as leptonSF_
leptonSF = leptonSF_()


# systematic variations
addSystematicVariations = (not isData) and (not options.skipSystematicVariations)
if addSystematicVariations:
    # B tagging SF
    from StopsDilepton.tools.btagEfficiency import btagEfficiency
    btagEff = btagEfficiency( fastSim = options.fastSim )

# LHE cut (DY samples)
if options.LHEHTCut>0:
    sample.name+="_lheHT"+str(options.LHEHTCut)
    logger.info( "Adding upper LHE cut at %f", options.LHEHTCut )
    skimConds.append( "lheHTIncoming<%f"%options.LHEHTCut )

# output directory
outDir = os.path.join(options.targetDir, options.processingEra, options.skim, sample.name)

# Directory for individual signal files
if options.T2tt:
    signalDir = os.path.join(options.targetDir, options.processingEra, options.skim, "T2tt")
    if not os.path.exists(signalDir): os.makedirs(signalDir)

if os.path.exists(outDir) and options.overwrite:
    if options.nJobs > 1:
        logger.warning( "NOT removing directory %s because nJobs = %i", outDir, options.nJobs )
    else:
        logger.info( "Output directory %s exists. Deleting.", outDir )
        shutil.rmtree(outDir)

try:    #Avoid trouble with race conditions in multithreading
    os.makedirs(outDir)
    logger.info( "Created output directory %s.", outDir )
except:
    pass

if isTiny:
    #branches to be kept for data and MC
    branchKeepStrings_DATAMC = \
       ["run", "lumi", "evt", "isData", "nVert",
        "met_pt", "met_phi",
#        "puppiMet_pt","puppiMet_phi",
        "Flag_*", "HLT_*",
        #"HLT_mumuIso", "HLT_ee_DZ", "HLT_mue",
        #"HLT_3mu", "HLT_3e", "HLT_2e1mu", "HLT_2mu1e",
        "LepGood_eta", "LepGood_etaSc", "LepGood_pt","LepGood_phi", "LepGood_dxy", "LepGood_dz","LepGood_tightId", "LepGood_pdgId",
        "LepGood_mediumMuonId", "LepGood_miniRelIso", "LepGood_relIso03", "LepGood_sip3d", "LepGood_mvaIdSpring15", "LepGood_convVeto", "LepGood_lostHits","LepGood_jetPtRelv2", "LepGood_jetPtRatiov2", "LepGood_eleCutId_Spring2016_25ns_v1_ConvVetoDxyDz"
        ]

    #branches to be kept for MC samples only
    branchKeepStrings_MC = [ \
        "nTrueInt", "genWeight", "xsec", "met_genPt", "met_genPhi", "lheHTIncoming"
    ]

    #branches to be kept for data only
    branchKeepStrings_DATA = [ ]
elif isSmall:
    #branches to be kept for data and MC
    branchKeepStrings_DATAMC = [\
        "run", "lumi", "evt", "isData", "rho", "nVert",
        "met_pt", "met_phi","met_Jet*", "met_Unclustered*", "met_sumEt", "met_rawPt","met_rawPhi", "met_rawSumEt",
#        "metNoHF_pt", "metNoHF_phi",
#        "puppiMet_pt","puppiMet_phi","puppiMet_sumEt","puppiMet_rawPt","puppiMet_rawPhi","puppiMet_rawSumEt",
        "Flag_*","HLT_*",
#        "nDiscJet", "DiscJet_*",
#        "nJetFailId", "JetFailId_*",
        "nLepGood", "LepGood_*",
        "nLepOther", "LepOther_*",
#        "nTauGood", "TauGood_*",
    ]
    #branches to be kept for MC samples only
    branchKeepStrings_MC = [\
        "nTrueInt", "genWeight", "xsec", "met_gen*", "lheHTIncoming",
#        "ngenPartAll","genPartAll_*","ngenLep","genLep_*"
    ]

    #branches to be kept for data only
    branchKeepStrings_DATA = [ ]
else:
    #branches to be kept for data and MC
    branchKeepStrings_DATAMC = [\
        "run", "lumi", "evt", "isData", "rho", "nVert",
        "met_pt", "met_phi","met_Jet*", "met_Unclustered*", "met_sumEt", "met_rawPt","met_rawPhi", "met_rawSumEt",
#        "metNoHF_pt", "metNoHF_phi",
#        "puppiMet_pt","puppiMet_phi","puppiMet_sumEt","puppiMet_rawPt","puppiMet_rawPhi","puppiMet_rawSumEt",
        "Flag_*","HLT_*",
        "nDiscJet", "DiscJet_*",
        "nJetFailId", "JetFailId_*",
        "nLepGood", "LepGood_*",
        "nLepOther", "LepOther_*",
        "nTauGood", "TauGood_*",
    ]
    #branches to be kept for MC samples only
    branchKeepStrings_MC = [\
        "nTrueInt", "genWeight", "xsec", "met_gen*", "lheHTIncoming",
        "ngenPartAll","genPartAll_*","ngenLep","genLep_*"
    ]

    #branches to be kept for data only
    branchKeepStrings_DATA = [ ]

if options.T2tt or options.TTDM:
    branchKeepStrings_MC += ["nIsr"]
if options.keepLHEWeights or options.T2tt or options.TTDM:
        branchKeepStrings_MC+=["nLHEweight", "LHEweight_id", "LHEweight_wgt", "LHEweight_original"]

if isSingleLep:
    branchKeepStrings_DATAMC += ['HLT_*']

if options.T2tt: branchKeepStrings_MC += ['GenSusyMStop', 'GenSusyMNeutralino']

# Jet variables to be read from chain
jetCorrInfo = ['corr/F', 'corr_JECUp/F', 'corr_JECDown/F'] if addSystematicVariations else []
if isMC:
    if isTiny or isSmall:
        jetMCInfo = ['mcPt/F', 'hadronFlavour/I']
    else:
        jetMCInfo = ['mcMatchFlav/I', 'partonId/I', 'partonMotherId/I', 'mcPt/F', 'mcFlavour/I', 'hadronFlavour/I', 'mcMatchId/I']
        if not options.T2tt:
            jetMCInfo.append('partonFlavour/I')
else:
    jetMCInfo = []

if not (isTiny or isSmall):
    branchKeepStrings_DATAMC+=[
        "ngamma", "gamma_idCutBased", "gamma_hOverE", "gamma_r9", "gamma_sigmaIetaIeta", "gamma_chHadIso04", "gamma_chHadIso", "gamma_phIso",
        "gamma_neuHadIso", "gamma_relIso", "gamma_pdgId", "gamma_pt", "gamma_eta", "gamma_phi", "gamma_mass",
        "gamma_chHadIsoRC04", "gamma_chHadIsoRC"]
    if isMC: branchKeepStrings_DATAMC+=[ "gamma_mcMatchId", "gamma_mcPt", "gamma_genIso04", "gamma_genIso03", "gamma_drMinParton"]

if sample.isData:
    lumiScaleFactor=None
    branchKeepStrings = branchKeepStrings_DATAMC + branchKeepStrings_DATA
    from FWCore.PythonUtilities.LumiList import LumiList
    # Apply golden JSON
    sample.heppy.json = '$CMSSW_BASE/src/CMGTools/TTHAnalysis/data/json/Cert_271036-284044_13TeV_PromptReco_Collisions16_JSON_NoL1T.txt'
    lumiList = LumiList(os.path.expandvars(sample.heppy.json))
    logger.info( "Loaded json %s", sample.heppy.json )
else:
    lumiScaleFactor = xSection*targetLumi/float(sample.normalization) if xSection is not None else None
    branchKeepStrings = branchKeepStrings_DATAMC + branchKeepStrings_MC

jetVars = ['pt/F', 'rawPt/F', 'eta/F', 'phi/F', 'id/I', 'btagCSV/F'] + jetCorrInfo + jetMCInfo
jetVarNames = [x.split('/')[0] for x in jetVars]
genLepVars      = ['pt/F', 'phi/F', 'eta/F', 'pdgId/I', 'index/I', 'lepGoodMatchIndex/I', 'matchesPromptGoodLepton/I', 'n_t/I','n_W/I', 'n_B/I', 'n_D/I', 'n_tau/I']
genLepVarNames  = [x.split('/')[0] for x in genLepVars]

read_variables = map(TreeVariable.fromString, ['met_pt/F', 'met_phi/F', 'run/I', 'lumi/I', 'evt/l', 'nVert/I'] )
read_variables += [ TreeVariable.fromString('ngamma/I'),
		    VectorTreeVariable.fromString('gamma[pt/F,eta/F,phi/F,mass/F,idCutBased/I,pdgId/I]') ]

new_variables = [ 'weight/F']
if isMC:
    read_variables+= [TreeVariable.fromString('nTrueInt/F')]
    # reading gen particles for top pt reweighting
    read_variables.append( TreeVariable.fromString('ngenPartAll/I') )
    read_variables.append( VectorTreeVariable.fromString('genPartAll[pt/F,eta/F,phi/F,pdgId/I,status/I,charge/I,motherId/I,grandmotherId/I,nMothers/I,motherIndex1/I,motherIndex2/I,nDaughters/I,daughterIndex1/I,daughterIndex2/I,isPromptHard/I]', nMax=200 )) # default nMax is 100, which would lead to corrupt values in this case
    read_variables.append( TreeVariable.fromString('genWeight/F') )
    read_variables.append( VectorTreeVariable.fromString('gamma[mcPt/F]') )

    new_variables.extend([ 'reweightTopPt/F', 'reweightPU/F','reweightPUUp/F','reweightPUDown/F', 'reweightPU12fb/F','reweightPU12fbUp/F','reweightPU12fbDown/F','reweightPU27fb/F','reweightPU27fbUp/F','reweightPU27fbDown/F'])
    if not options.skipGenLepMatching:
        TreeVariable.fromString( 'nGenLep/I' ),
        new_variables.append( 'GenLep[%s]'% ( ','.join(genLepVars) ) )

read_variables += [\
    TreeVariable.fromString('nLepGood/I'),
    VectorTreeVariable.fromString('LepGood[pt/F,eta/F,etaSc/F,phi/F,pdgId/I,tightId/I,miniRelIso/F,relIso03/F,sip3d/F,mediumMuonId/I,mvaIdSpring15/F,lostHits/I,convVeto/I,dxy/F,dz/F,jetPtRelv2/F,jetPtRatiov2/F,eleCutId_Spring2016_25ns_v1_ConvVetoDxyDz/I]'),
    TreeVariable.fromString('nJet/I'),
    VectorTreeVariable.fromString('Jet[%s]'% ( ','.join(jetVars) ) )
]
if isVeryLoose:
    read_variables += [\
        TreeVariable.fromString('nLepOther/I'),
        VectorTreeVariable.fromString('LepOther[pt/F,eta/F,etaSc/F,phi/F,pdgId/I,tightId/I,miniRelIso/F,relIso03/F,sip3d/F,mediumMuonId/I,mvaIdSpring15/F,lostHits/I,convVeto/I,dxy/F,dz/F,jetPtRelv2/F,jetPtRatiov2/F,eleCutId_Spring2016_25ns_v1_ConvVetoDxyDz/I]'),
    ]
new_variables += [\
    'JetGood[%s]'% ( ','.join(jetVars) )
]

if isData: new_variables.extend( ['jsonPassed/I'] )
new_variables.extend( ['nJetGood/I','nBTag/I', 'ht/F', 'metSig/F'] )

if isSingleLep:
    new_variables.extend( ['m3/F', 'm3_ind1/I', 'm3_ind2/I', 'm3_ind3/I'] )
if isTriLep or isDiLep or isSingleLep:
    new_variables.extend( ['nGoodMuons/I', 'nGoodElectrons/I' ] )
    new_variables.extend( ['l1_pt/F', 'l1_eta/F', 'l1_phi/F', 'l1_pdgId/I', 'l1_index/I', 'l1_jetPtRelv2/F', 'l1_jetPtRatiov2/F', 'l1_miniRelIso/F', 'l1_relIso03/F', 'l1_dxy/F', 'l1_dz/F', 'l1_mIsoWP/I' ] )
    # new_variables.extend( ['mt/F', 'mlmZ_mass/F'] )
    new_variables.extend( ['mlmZ_mass/F'] )
    new_variables.extend( ['mt_photonEstimated/F'] )
if isTriLep or isDiLep:
    new_variables.extend( ['l2_pt/F', 'l2_eta/F', 'l2_phi/F', 'l2_pdgId/I', 'l2_index/I', 'l2_jetPtRelv2/F', 'l2_jetPtRatiov2/F', 'l2_miniRelIso/F', 'l2_relIso03/F', 'l2_dxy/F', 'l2_dz/F', 'l2_mIsoWP/I' ] )
    new_variables.extend( ['isEE/I', 'isMuMu/I', 'isEMu/I', 'isOS/I' ] )
    new_variables.extend( ['dl_pt/F', 'dl_eta/F', 'dl_phi/F', 'dl_mass/F'] )
    new_variables.extend( ['dl_mt2ll/F', 'dl_mt2bb/F', 'dl_mt2blbl/F' ] )
    if isMC: new_variables.extend( \
        [   'zBoson_genPt/F', 'zBoson_genEta/F', 
            'reweightDilepTrigger/F', 'reweightDilepTriggerUp/F', 'reweightDilepTriggerDown/F', 'reweightDilepTriggerBackup/F', 'reweightDilepTriggerBackupUp/F', 'reweightDilepTriggerBackupDown/F',
            'reweightLeptonSF/F', 'reweightLeptonSFUp/F', 'reweightLeptonSFDown/F',
            'reweightLeptonHIPSF/F',
         ] )
    if options.T2tt or options.TTDM:
        new_variables.extend( ['dl_mt2ll_gen/F', 'dl_mt2bb_gen/F', 'dl_mt2blbl_gen/F' ] )
new_variables.extend( ['nPhotonGood/I','photon_pt/F','photon_eta/F','photon_phi/F','photon_idCutBased/I'] )
if isMC: new_variables.extend( ['photon_genPt/F', 'photon_genEta/F'] )
new_variables.extend( ['met_pt_photonEstimated/F','met_phi_photonEstimated/F','metSig_photonEstimated/F'] )
new_variables.extend( ['photonJetdR/F','photonLepdR/F'] )
if isTriLep or isDiLep:
  new_variables.extend( ['dlg_mass/F','dl_mt2ll_photonEstimated/F', 'dl_mt2bb_photonEstimated/F', 'dl_mt2blbl_photonEstimated/F' ] )

if options.checkTTGJetsOverlap:
    new_variables.extend( ['TTGJetsEventType/I'] )

if addSystematicVariations:
    read_variables += map(TreeVariable.fromString, [\
    "met_JetEnUp_Pt/F", "met_JetEnUp_Phi/F", "met_JetEnDown_Pt/F", "met_JetEnDown_Phi/F", "met_JetResUp_Pt/F", "met_JetResUp_Phi/F", "met_JetResDown_Pt/F", "met_JetResDown_Phi/F", 
    "met_UnclusteredEnUp_Pt/F", "met_UnclusteredEnUp_Phi/F", "met_UnclusteredEnDown_Pt/F", "met_UnclusteredEnDown_Phi/F", 
    ] )

    for var in ['JECUp', 'JECDown', 'JERUp', 'JERDown', 'UnclusteredEnUp', 'UnclusteredEnDown']:
        if 'Unclustered' not in var: new_variables.extend( ['nJetGood_'+var+'/I', 'nBTag_'+var+'/I','ht_'+var+'/F'] )
        new_variables.extend( ['met_pt_'+var+'/F', 'met_phi_'+var+'/F', 'metSig_'+var+'/F'] )
        if isTriLep or isDiLep:
            new_variables.extend( ['dl_mt2ll_'+var+'/F', 'dl_mt2bb_'+var+'/F', 'dl_mt2blbl_'+var+'/F'] )
	new_variables.extend( ['met_pt_photonEstimated_'+var+'/F', 'met_phi_photonEstimated_'+var+'/F', 'metSig_photonEstimated_'+var+'/F'] )
	if isTriLep or isDiLep:
	    new_variables.extend( ['dl_mt2ll_photonEstimated_'+var+'/F', 'dl_mt2bb_photonEstimated_'+var+'/F', 'dl_mt2blbl_photonEstimated_'+var+'/F'] )
    # Btag weights Method 1a
    for var in btagEff.btagWeightNames:
        if var!='MC':
            new_variables.append('reweightBTag_'+var+'/F')

if options.T2tt or options.TTDM:
    read_variables += map(TreeVariable.fromString, ['met_genPt/F', 'met_genPhi/F'] )
if options.T2tt:
    read_variables += map(TreeVariable.fromString, ['GenSusyMStop/I', 'GenSusyMNeutralino/I'] )
    new_variables  += ['reweightXSecUp/F', 'reweightXSecDown/F', 'mStop/I', 'mNeu/I']

if options.fastSim and (isTriLep or isDiLep):
    new_variables  += ['reweightLeptonFastSimSF/F', 'reweightLeptonFastSimSFUp/F', 'reweightLeptonFastSimSFDown/F']


# Define a reader
reader = sample.treeReader( \
    variables = read_variables ,
    selectionString = "&&".join(skimConds)
    )

# Calculate photonEstimated met
def getMetPhotonEstimated(met_pt, met_phi, photon):
  met = ROOT.TLorentzVector()
  met.SetPtEtaPhiM(met_pt, 0, met_phi, 0 )
  gamma = ROOT.TLorentzVector()
  gamma.SetPtEtaPhiM(photon['pt'], photon['eta'], photon['phi'], photon['mass'] )
  metGamma = met + gamma
  return (metGamma.Pt(), metGamma.Phi())


## Calculate corrected met pt/phi using systematics for jets
def getMetJetCorrected(met_pt, met_phi, jets, var):
  met_corr_px  = met_pt*cos(met_phi) + sum([(j['pt']-j['pt_'+var])*cos(j['phi']) for j in jets])
  met_corr_py  = met_pt*sin(met_phi) + sum([(j['pt']-j['pt_'+var])*sin(j['phi']) for j in jets])
  met_corr_pt  = sqrt(met_corr_px**2 + met_corr_py**2)
  met_corr_phi = atan2(met_corr_py, met_corr_px)
  return (met_corr_pt, met_corr_phi)

def getMetCorrected(r, var, addPhoton = None):
    if var == "":
      if addPhoton is not None: return getMetPhotonEstimated(r.met_pt, r.met_phi, addPhoton)
      else:                     return (r.met_pt, r.met_phi)

    elif var in [ "JECUp", "JECDown", "UnclusteredEnUp", "UnclusteredEnDown" ]:
      var_ = var
      var_ = var_.replace("JEC", "JetEn")
      var_ = var_.replace("JER", "JetRes")
      met_pt  = getattr(r, "met_" + var_ + "_Pt")
      met_phi = getattr(r, "met_" + var_ + "_Phi")
      if addPhoton is not None: return getMetPhotonEstimated(met_pt, met_phi, addPhoton)
      else:                     return (met_pt, met_phi)

    else:
        raise ValueError

mothers = {"D":0, "B":0}
grannies_D = {}
grannies_B = {}

def filler( event ):
    # shortcut
    r = reader.event
    if isMC: gPart = getGenPartsAll(r)

    # weight
    if options.T2tt:
        event.weight=signalWeight[(r.GenSusyMStop, r.GenSusyMNeutralino)]['weight']
        event.mStop = r.GenSusyMStop
        event.mNeu  = r.GenSusyMNeutralino
        event.reweightXSecUp    = signalWeight[(r.GenSusyMStop, r.GenSusyMNeutralino)]['xSecFacUp']
        event.reweightXSecDown  = signalWeight[(r.GenSusyMStop, r.GenSusyMNeutralino)]['xSecFacDown']
    elif isMC:
        event.weight = lumiScaleFactor*r.genWeight if lumiScaleFactor is not None else 1
    elif isData:
        event.weight = 1
    else:
        raise NotImplementedError( "isMC %r isData %r T2tt? %r TTDM?" % (isMC, isData, options.T2tt, options.TTDM) )

    # lumi lists and vetos
    if isData:
        #event.vetoPassed  = vetoList.passesVeto(r.run, r.lumi, r.evt)
        event.jsonPassed  = lumiList.contains(r.run, r.lumi)
        # store decision to use after filler has been executed
        event.jsonPassed_ = event.jsonPassed

    if isMC:
        event.reweightPU     = nTrueInt_puRW       ( r.nTrueInt ) 
        event.reweightPUDown = nTrueInt_puRWDown   ( r.nTrueInt ) 
        event.reweightPUUp   = nTrueInt_puRWUp     ( r.nTrueInt ) 
        event.reweightPU12fb     = nTrueInt12fb_puRW       ( r.nTrueInt ) 
        event.reweightPU12fbDown = nTrueInt12fb_puRWDown   ( r.nTrueInt ) 
        event.reweightPU12fbUp   = nTrueInt12fb_puRWUp     ( r.nTrueInt ) 
        event.reweightPU27fb     = nTrueInt27fb_puRW       ( r.nTrueInt )
        event.reweightPU27fbDown = nTrueInt27fb_puRWDown   ( r.nTrueInt )
        event.reweightPU27fbUp   = nTrueInt27fb_puRWUp     ( r.nTrueInt )

    # top pt reweighting
    if isMC: event.reweightTopPt = topPtReweightingFunc(getTopPtsForReweighting(r))/topScaleF if doTopPtReweighting else 1.

    # jet/met related quantities, also load the leptons already
    if options.keepForwardJets:
        jetAbsEtaCut = 99.
    else:
        jetAbsEtaCut = 2.4
    allJets      = getGoodJets(r, ptCut=0, jetVars = jetVarNames, absEtaCut=jetAbsEtaCut)
    jets         = filter(lambda j:jetId(j, ptCut=30, absEtaCut=jetAbsEtaCut), allJets)
    bJets        = filter(lambda j:isBJet(j), jets)
    nonBJets     = filter(lambda j:not isBJet(j), jets)
    if isVeryLoose:
        raise NotImplementedError
        ## all leptons up to relIso 0.5
        #mu_selector = muonSelector( iso = 999., dxy = 1., dz = 0.1 )
        #ele_selector = eleSelector( iso = 999., dxy = 1., dz = 0.1 )
        #leptons_pt10 = getGoodAndOtherLeptons(r, ptCut=10, mu_selector = mu_selector, ele_selector = ele_selector)
        #ptCut = 20 if not isVeryLoosePt10 else 10 
        #leptons      = filter(lambda l:l['pt']>ptCut, leptons_pt10)
    elif isLoose:
        raise NotImplementedError
        ## loose relIso lepton selection
        #mu_selector = muonSelector( iso = 0.4 )
        #ele_selector = eleSelector( iso = 0.4 )

        #leptons_pt10 = getGoodLeptons(r, ptCut=10, mu_selector = mu_selector, ele_selector = ele_selector)
        #leptons      = filter(lambda l:l['pt']>20, leptons_pt10)
    else:
        # using miniRelIso 0.2 as baseline 
        mu_selector = muonSelector( relIso03 = 0.2 )
        ele_selector = eleSelector( relIso03 = 0.2 )
        leptons_pt10 = getGoodLeptons(r, ptCut=10, mu_selector = mu_selector, ele_selector = ele_selector)
        leptons      = filter(lambda l:l['pt']>20, leptons_pt10)

    leptons.sort(key = lambda p:-p['pt'])
    
    event.met_pt  = r.met_pt
    event.met_phi = r.met_phi

    # Filling jets
    event.nJetGood   = len(jets)
    for iJet, jet in enumerate(jets):
        for b in jetVarNames:
            getattr(event, "JetGood_"+b)[iJet] = jet[b]
    if isSingleLep:
        # Compute M3 and the three indiced of the jets entering m3
        event.m3, event.m3_ind1, event.m3_ind2, event.m3_ind3 = m3( jets )

    event.ht         = sum([j['pt'] for j in jets])
    event.metSig     = event.met_pt/sqrt(event.ht) if event.ht>0 else float('nan')
    event.nBTag      = len(bJets)

    jets_sys      = {}
    bjets_sys     = {}
    nonBjets_sys  = {}

    metVariants = [''] # default

    # Keep photons and estimate met including (leading pt) photon
    photons = getGoodPhotons(r, ptCut=20, idLevel="loose", isData=isData)
    event.nPhotonGood = len(photons)
    if event.nPhotonGood > 0:
      metVariants += ['_photonEstimated']  # do all met calculations also for the photonEstimated variant
      event.photon_pt         = photons[0]['pt']
      event.photon_eta        = photons[0]['eta']
      event.photon_phi        = photons[0]['phi']
      event.photon_idCutBased = photons[0]['idCutBased']
      if isMC:
        genPhoton       = getGenPhoton(gPart)
        event.photon_genPt  = genPhoton['pt']  if genPhoton is not None else float('nan')
        event.photon_genEta = genPhoton['eta'] if genPhoton is not None else float('nan')

      event.met_pt_photonEstimated, event.met_phi_photonEstimated = getMetPhotonEstimated(r.met_pt, r.met_phi, photons[0])
      event.metSig_photonEstimated = event.met_pt_photonEstimated/sqrt(event.ht) if event.ht>0 else float('nan')

      event.photonJetdR = min(deltaR(photons[0], j) for j in jets) if len(jets) > 0 else 999
      event.photonLepdR = min(deltaR(photons[0], l) for l in leptons_pt10) if len(leptons_pt10) > 0 else 999

    if options.checkTTGJetsOverlap and isMC:
       event.TTGJetsEventType = getTTGJetsEventType(r)

    if addSystematicVariations:
        for j in allJets:
            j['pt_JECUp']   =j['pt']/j['corr']*j['corr_JECUp']
            j['pt_JECDown'] =j['pt']/j['corr']*j['corr_JECDown']
            # JERUp, JERDown, JER
            addJERScaling(j)
        for var in ['JECUp', 'JECDown', 'JERUp', 'JERDown']:
            jets_sys[var]       = filter(lambda j:jetId(j, ptCut=30, absEtaCut=jetAbsEtaCut, ptVar='pt_'+var), allJets)
            bjets_sys[var]      = filter(isBJet, jets_sys[var])
            nonBjets_sys[var]   = filter(lambda j: not isBJet(j), jets_sys[var])

            setattr(event, "nJetGood_"+var, len(jets_sys[var]))
            setattr(event, "ht_"+var,       sum([j['pt_'+var] for j in jets_sys[var]]))
            setattr(event, "nBTag_"+var,    len(bjets_sys[var]))

        for var in ['JECUp', 'JECDown', 'JERUp', 'JERDown', 'UnclusteredEnUp', 'UnclusteredEnDown']:
            for i in metVariants:
                # use cmg MET correction values ecept for JER where it is zero. There, propagate jet variations.
                if 'JER' in var or 'JECV' in var:
                  (met_corr_pt, met_corr_phi) = getMetJetCorrected(getattr(event, "met_pt" + i), getattr(event,"met_phi" + i), jets_sys[var], var)
                else:
                  (met_corr_pt, met_corr_phi) = getMetCorrected(r, var, photons[0] if i.count("photonEstimated") else None)

                setattr(event, "met_pt" +i+"_"+var, met_corr_pt)
                setattr(event, "met_phi"+i+"_"+var, met_corr_phi)
                ht = getattr(event, "ht_"+var) if 'Unclustered' not in var else event.ht 
                setattr(event, "metSig" +i+"_"+var, getattr(event, "met_pt"+i+"_"+var)/sqrt( ht ) if ht>0 else float('nan') )

    if isSingleLep or isTriLep or isDiLep:
        event.nGoodMuons      = len(filter( lambda l:abs(l['pdgId'])==13, leptons))
        event.nGoodElectrons  = len(filter( lambda l:abs(l['pdgId'])==11, leptons))

        if len(leptons)>=1:
            event.l1_pt     = leptons[0]['pt']
            event.l1_eta    = leptons[0]['eta']
            event.l1_phi    = leptons[0]['phi']
            event.l1_pdgId  = leptons[0]['pdgId']
            event.l1_index  = leptons[0]['index']
            event.l1_jetPtRatiov2  = leptons[0]['jetPtRatiov2']
            event.l1_jetPtRelv2    = leptons[0]['jetPtRelv2']
            event.l1_miniRelIso = leptons[0]['miniRelIso']
            event.l1_relIso03 = leptons[0]['relIso03']
            event.l1_dxy = leptons[0]['dxy']
            event.l1_dz = leptons[0]['dz']
            event.l1_mIsoWP = multiIsoWPInt(leptons[0])

        # For TTZ studies: find Z boson candidate, and use third lepton to calculate mt
        (event.mlmZ_mass, zl1, zl2) = closestOSDLMassToMZ(leptons_pt10)
#        if len(leptons_pt10) >= 3:
#            thirdLepton = leptons_pt10[[x for x in range(len(leptons_pt10)) if x != zl1 and x != zl2][0]]
#            for i in metVariants:
#              setattr(event, "mt"+i, sqrt(2*thirdLepton['pt']*getattr(event, "met_pt"+i)*(1-cos(thirdLepton['phi']-getattr(event, "met_phi"+i)))))

        if options.fastSim:
            event.reweightLeptonFastSimSF     = reduce(mul, [leptonFastSimSF.get2DSF(pdgId=l['pdgId'], pt=l['pt'], eta=l['eta'] , nvtx = r.nVert) for l in leptons], 1)
            event.reweightLeptonFastSimSFUp   = reduce(mul, [leptonFastSimSF.get2DSF(pdgId=l['pdgId'], pt=l['pt'], eta=l['eta'] , nvtx = r.nVert, sigma = +1) for l in leptons], 1)
            event.reweightLeptonFastSimSFDown = reduce(mul, [leptonFastSimSF.get2DSF(pdgId=l['pdgId'], pt=l['pt'], eta=l['eta'] , nvtx = r.nVert, sigma = -1) for l in leptons], 1)

    if isTriLep or isDiLep:
        if isMC:
            event.reweightDilepTrigger       = 0 
            event.reweightDilepTriggerUp     = 0 
            event.reweightDilepTriggerDown   = 0 
            event.reweightDilepTriggerBackup       = 0 
            event.reweightDilepTriggerBackupUp     = 0 
            event.reweightDilepTriggerBackupDown   = 0 

            leptonsForSF = (leptons[:2] if isDiLep else (leptons[:3] if isTriLep else leptons))
            event.reweightLeptonSF           = reduce(mul, [leptonSF.getSF(pdgId=l['pdgId'], pt=l['pt'], eta=l['eta']) for l in leptonsForSF], 1)
            event.reweightLeptonSFUp         = reduce(mul, [leptonSF.getSF(pdgId=l['pdgId'], pt=l['pt'], eta=l['eta'] , sigma = +1) for l in leptonsForSF], 1)
            event.reweightLeptonSFDown       = reduce(mul, [leptonSF.getSF(pdgId=l['pdgId'], pt=l['pt'], eta=l['eta'] , sigma = -1) for l in leptonsForSF], 1)

            event.reweightLeptonHIPSF      = reduce(mul, [leptonHIPSF.getSF( \
                pdgId = l['pdgId'],
                pt  =   l['pt'], 
                eta =   (l['etaSc'] if abs(l['pdgId'])==11 else l['eta'])
                )[0]  for l in leptonsForSF], 1)

        if len(leptons)>=2:
            mt2Calc.reset()
            event.l2_pt     = leptons[1]['pt']
            event.l2_eta    = leptons[1]['eta']
            event.l2_phi    = leptons[1]['phi']
            event.l2_pdgId  = leptons[1]['pdgId']
            event.l2_index  = leptons[1]['index']
            event.l2_jetPtRatiov2  = leptons[1]['jetPtRatiov2']
            event.l2_jetPtRelv2    = leptons[1]['jetPtRelv2']
            event.l2_miniRelIso = leptons[1]['miniRelIso']
            event.l2_relIso03 = leptons[1]['relIso03']
            event.l2_dxy = leptons[1]['dxy']
            event.l2_dz = leptons[1]['dz']
            event.l2_mIsoWP = multiIsoWPInt(leptons[1])
            

            l_pdgs = [abs(leptons[0]['pdgId']), abs(leptons[1]['pdgId'])]
            l_pdgs.sort()
            event.isMuMu = l_pdgs==[13,13]
            event.isEE   = l_pdgs==[11,11]
            event.isEMu  = l_pdgs==[11,13]
            event.isOS   = event.l1_pdgId*event.l2_pdgId<0

            l1 = ROOT.TLorentzVector()
            l1.SetPtEtaPhiM(leptons[0]['pt'], leptons[0]['eta'], leptons[0]['phi'], 0 )
            l2 = ROOT.TLorentzVector()
            l2.SetPtEtaPhiM(leptons[1]['pt'], leptons[1]['eta'], leptons[1]['phi'], 0 )
            dl = l1+l2
            event.dl_pt   = dl.Pt()
            event.dl_eta  = dl.Eta()
            event.dl_phi  = dl.Phi()
            event.dl_mass = dl.M()
            mt2Calc.setLeptons(event.l1_pt, event.l1_eta, event.l1_phi, event.l2_pt, event.l2_eta, event.l2_phi)

            # To check MC truth when looking at the TTZToLLNuNu sample
            if isMC:
              trig_eff, trig_eff_err =  triggerEff.getSF(event.l1_pt, event.l1_eta, event.l1_pdgId, event.l2_pt, event.l2_eta, event.l2_pdgId)
              event.reweightDilepTrigger       = trig_eff 
              event.reweightDilepTriggerUp     = trig_eff + trig_eff_err
              event.reweightDilepTriggerDown   = trig_eff - trig_eff_err
              trig_eff, trig_eff_err =  triggerEff_withBackup.getSF(event.l1_pt, event.l1_eta, event.l1_pdgId, event.l2_pt, event.l2_eta, event.l2_pdgId)
              event.reweightDilepTriggerBackup       = trig_eff 
              event.reweightDilepTriggerBackupUp     = trig_eff + trig_eff_err
              event.reweightDilepTriggerBackupDown   = trig_eff - trig_eff_err

              zBoson          = getGenZ(gPart)
              event.zBoson_genPt  = zBoson['pt']  if zBoson is not None else float('nan')
              event.zBoson_genEta = zBoson['eta'] if zBoson is not None else float('nan')

            if event.nPhotonGood > 0:
              gamma = ROOT.TLorentzVector()
              gamma.SetPtEtaPhiM(photons[0]['pt'], photons[0]['eta'], photons[0]['phi'], photons[0]['mass'] )
              dlg = dl + gamma
              event.dlg_mass = dlg.M()

            if options.T2tt or options.TTDM:
                mt2Calc.setMet(getattr(r, 'met_genPt'), getattr(r, 'met_genPhi'))
                setattr(event, "dl_mt2ll_gen", mt2Calc.mt2ll())
                if len(jets)>=2:
                    bj0, bj1 = (bJets+nonBJets)[:2]
                    mt2Calc.setBJets(bj0['pt'], bj0['eta'], bj0['phi'], bj1['pt'], bj1['eta'], bj1['phi'])
                    setattr(event, "dl_mt2bb_gen",   mt2Calc.mt2bb())
                    setattr(event, "dl_mt2blbl_gen", mt2Calc.mt2blbl())
                
            for i in metVariants:
                mt2Calc.setMet(getattr(event, 'met_pt'+i), getattr(event, 'met_phi'+i))
                setattr(event, "dl_mt2ll"+i, mt2Calc.mt2ll())

                bj0, bj1 = None, None
                if len(jets)>=2:
                    bj0, bj1 = (bJets+nonBJets)[:2]
                    mt2Calc.setBJets(bj0['pt'], bj0['eta'], bj0['phi'], bj1['pt'], bj1['eta'], bj1['phi'])
                    setattr(event, "dl_mt2bb"+i,   mt2Calc.mt2bb())
                    setattr(event, "dl_mt2blbl"+i, mt2Calc.mt2blbl())

                if addSystematicVariations:
                    for var in ['JECUp', 'JECDown', 'JERUp', 'JERDown', 'UnclusteredEnUp', 'UnclusteredEnDown']:
                        mt2Calc.setMet( getattr(event, "met_pt"+i+"_"+var), getattr(event, "met_phi"+i+"_"+var) )
                        setattr(event, "dl_mt2ll"+i+"_"+var,  mt2Calc.mt2ll())
                        if not 'Unclustered' in var:
                            if len(jets_sys[var])>=2:
                                bj0_, bj1_ = (bjets_sys[var]+nonBjets_sys[var])[:2]
                            else: 
                                bj0_, bj1_ = None, None
                        else:
                            bj0_, bj1_ = bj0, bj1
                        if bj0_ and bj1_:
                            mt2Calc.setBJets(bj0_['pt'], bj0_['eta'], bj0_['phi'], bj1_['pt'], bj1_['eta'], bj1_['phi'])
                            setattr(event, 'dl_mt2bb'  +i+'_'+var, mt2Calc.mt2bb())
                            setattr(event, 'dl_mt2blbl'+i+'_'+var, mt2Calc.mt2blbl())

    if addSystematicVariations:
        # B tagging weights method 1a
        for j in jets:
            btagEff.addBTagEffToJet(j)
        for var in btagEff.btagWeightNames:
            if var!='MC':
                setattr(event, 'reweightBTag_'+var, btagEff.getBTagSF_1a( var, bJets, nonBJets ) )

    # gen information on extra leptons
    if isMC and not options.skipGenLepMatching:
        genSearch.init( gPart )
        # Start with status 1 gen leptons in acceptance
        gLep = filter( lambda p:abs(p['pdgId']) in [11, 13] and p['status']==1 and p['pt']>20 and abs(p['eta'])<2.5, gPart )
        for l in gLep:
            ancestry = [ gPart[x]['pdgId'] for x in genSearch.ancestry( l ) ]
            l["n_D"]   =  sum([ancestry.count(p) for p in D_mesons])
            l["n_B"]   =  sum([ancestry.count(p) for p in B_mesons])
            l["n_W"]   =  sum([ancestry.count(p) for p in [24, -24]])
            l["n_t"]   =  sum([ancestry.count(p) for p in [6, -6]])
            l["n_tau"] =  sum([ancestry.count(p) for p in [15, -15]])
            matched_lep = bestDRMatchInCollection(l, leptons_pt10)
            if matched_lep:
                l["lepGoodMatchIndex"] = matched_lep['index']
                if isSingleLep:
                    l["matchesPromptGoodLepton"] = l["lepGoodMatchIndex"] in [event.l1_index]
                elif isTriLep or isDiLep:
                    l["matchesPromptGoodLepton"] = l["lepGoodMatchIndex"] in [event.l1_index, event.l2_index]
                else:
                    l["matchesPromptGoodLepton"] = 0
            else:
                l["lepGoodMatchIndex"] = -1
                l["matchesPromptGoodLepton"] = 0
#            if      l["n_t"]>0 and l["n_W"]>0 and l["n_B"]==0 and l["n_D"]==0 and l["n_tau"]==0:
#                print "t->W->l"
#            elif    l["n_t"]>0 and l["n_W"]==0 and l["n_B"]>0 and l["n_D"]==0 and l["n_tau"]==0:
#                print "t->b->B->l"
#            elif    l["n_t"]>0 and l["n_W"]==0 and l["n_B"]>0 and l["n_D"]>0 and l["n_tau"]==0:
#                print "t->b->B->D->l"
#            elif    l["n_t"]>0 and l["n_W"]>0 and l["n_B"]==0 and l["n_D"]==0 and l["n_tau"]>0 :
#                print "t->W->tau->l"
#            elif    l["n_t"]>0 and l["n_W"]>0 and l["n_B"]==0 and l["n_D"]>0 and l["n_tau"]==0:
#                print "t->W->c->D->l"
#            elif    l["n_t"]==0 and l["n_W"]==0 and l["n_B"]>0 and l["n_D"]>=0 and l["n_tau"]==0:
#                print l['pdgId'], l['pt'], l['phi'], l['eta'], ",".join(pdgToName( gPart[x]['pdgId']) for x in genSearch.ancestry(l) )
#                for p in genSearch.ancestry(l):
#                    print p, gPart[p]
#            else:
#                pass
                # print l['pdgId'], l['pt'], l['phi'], l['eta'], ",".join(pdgToName(gPart[x]['pdgId']) for x in genSearch.ancestry(l))
        event.nGenLep   = len(gLep)
        for iLep, lep in enumerate(gLep):
            for b in genLepVarNames:
                getattr(event, "GenLep_"+b)[iLep] = lep[b]

# Create a maker. Maker class will be compiled. This instance will be used as a parent in the loop
treeMaker_parent = TreeMaker(
    sequence  = [ filler ],
    variables = [ TreeVariable.fromString(x) for x in new_variables ],
    treeName = "Events"
    )

# Split input in ranges
if options.nJobs>1:
    eventRanges = reader.getEventRanges( nJobs = options.nJobs )
else:
    eventRanges = reader.getEventRanges( maxNEvents = options.eventsPerJob, minJobs = options.minNJobs )

logger.info( "Splitting into %i ranges of %i events on average.",  len(eventRanges), (eventRanges[-1][1] - eventRanges[0][0])/len(eventRanges) )

#Define all jobs
jobs = [(i, eventRanges[i]) for i in range(len(eventRanges))]

filename, ext = os.path.splitext( os.path.join(outDir, sample.name + '.root') )

clonedEvents = 0
convertedEvents = 0
outputLumiList = {}
for ievtRange, eventRange in enumerate( eventRanges ):

    if len(options.job)>0 and not ievtRange in options.job: continue

    logger.info( "Processing range %i/%i from %i to %i which are %i events.",  ievtRange, len(eventRanges), eventRange[0], eventRange[1], eventRange[1]-eventRange[0] )

    # Check whether file exists
    outfilename = filename+'_'+str(ievtRange)+ext
    if os.path.isfile(outfilename):
        logger.info( "Output file %s found.", outfilename)
        if not checkRootFile(outfilename, checkForObjects=["Events"]):
            logger.info( "File %s is broken. Overwriting.", outfilename)
        elif not options.overwrite:
            logger.info( "Skipping.")
            continue
        else:
            logger.info( "Overwriting.")

    tmp_directory = ROOT.gDirectory
    outputfile = ROOT.TFile.Open(outfilename, 'recreate')
    tmp_directory.cd()

    # Set the reader to the event range
    reader.setEventRange( eventRange )
    clonedTree = reader.cloneTree( branchKeepStrings, newTreename = "Events", rootfile = outputfile )
    clonedEvents += clonedTree.GetEntries()

    # Clone the empty maker in order to avoid recompilation at every loop iteration
    maker = treeMaker_parent.cloneWithoutCompile( externalTree = clonedTree )

    maker.start()
    # Do the thing
    reader.start()

    while reader.run():
        maker.run()
        if isData:
            if maker.event.jsonPassed_:
                if reader.event.run not in outputLumiList.keys():
                    outputLumiList[reader.event.run] = {reader.event.lumi}
                else:
                    if reader.event.lumi not in outputLumiList[reader.event.run]:
                        outputLumiList[reader.event.run].add(reader.event.lumi)

    convertedEvents += maker.tree.GetEntries()
    maker.tree.Write()
    outputfile.Close()
    logger.info( "Written %s", outfilename)

  # Destroy the TTree
    maker.clear()


logger.info( "Converted %i events of %i, cloned %i",  convertedEvents, reader.nEvents , clonedEvents )

# Storing JSON file of processed events
if isData:
    jsonFile = filename+'.json'
    LumiList( runsAndLumis = outputLumiList ).writeJSON(jsonFile)
    logger.info( "Written JSON file %s",  jsonFile )

# Write one file per mass point for T2tt
if options.T2tt:
    output = Sample.fromDirectory("T2tt_output", outDir)
    for s in signalWeight.keys():
        cut = "GenSusyMStop=="+str(s[0])+"&&GenSusyMNeutralino=="+str(s[1])
        signalFile = os.path.join(signalDir, 'T2tt_'+str(s[0])+'_'+str(s[1])+'.root' )
        if not os.path.exists(signalFile) or options.overwrite:
            t = output.chain.CopyTree(cut)
            writeObjToFile(signalFile, t)
            logger.info( "Written signal file for masses mStop %i mNeu %i to %s", s[0], s[1], signalFile)
        else:
            logger.info( "Found file %s -> Skipping"%(signalFile) )

    output.clear()

