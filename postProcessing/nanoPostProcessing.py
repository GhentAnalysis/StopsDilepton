#!/usr/bin/env python

# standard imports
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
import sys
import os
import copy
import random
import subprocess
import datetime
import shutil
import uuid

theano_compile_dir = '/tmp/%s'%str(uuid.uuid4())
if not os.path.exists( theano_compile_dir ):
    os.makedirs( theano_compile_dir )
os.environ['THEANO_FLAGS'] = 'base_compiledir=%s'%theano_compile_dir 

from array import array
from operator import mul
from math import sqrt, atan2, sin, cos

# RootTools
from RootTools.core.standard import *

# User specific
import StopsDilepton.tools.user as user

# Tools for systematics
from StopsDilepton.tools.mt2Calculator       import mt2Calculator
from StopsDilepton.tools.helpers             import closestOSDLMassToMZ, writeObjToFile, m3, deltaR, bestDRMatchInCollection, deltaPhi, nonEmptyFile, getSortedZCandidates, getMinDLMass
from StopsDilepton.tools.addJERScaling       import addJERScaling
from StopsDilepton.tools.objectSelection     import muonSelector, eleSelector, getGoodMuons, getGoodElectrons,  getGoodJets, isBJet, jetId, isBJet, getGoodPhotons, getGenPartsAll, getJets, getPhotons, getAllJets, filterGenPhotons, genPhotonSelector, mergeCollections, genLepFromZ
from StopsDilepton.tools.getGenBoson         import getGenZ, getGenPhoton
from StopsDilepton.tools.polReweighting      import getPolWeights
from StopsDilepton.tools.triggerEfficiency   import triggerEfficiency
from StopsDilepton.tools.leptonSF            import leptonSF as leptonSF_
from StopsDilepton.tools.leptonHit0SF        import leptonHit0SF as leptonHit0SF_
from StopsDilepton.tools.leptonSip3dSF       import leptonSip3dSF as leptonSip3dSF_
from StopsDilepton.tools.leptonFastSimSF     import leptonFastSimSF as leptonFastSimSF_
from StopsDilepton.samples.helpers           import getTTDMSignalWeightForEvent, getTTDMSignalWeight, getTTDMBranchNames
from Analysis.Tools.overlapRemovalTTG        import photonFromTopDecay, hasMesonMother, getParentIds, isIsolatedPhoton, getPhotonCategory
from Analysis.Tools.puProfileCache           import *
from Analysis.Tools.L1PrefireWeight          import L1PrefireWeight
from Analysis.Tools.LeptonTrackingEfficiency import LeptonTrackingEfficiency
from Analysis.Tools.isrWeight                import ISRweight
from Analysis.Tools.helpers                  import checkRootFile, deepCheckRootFile, deepCheckWeight

#MC tools
from StopsDilepton.tools.mcTools import pdgToName, GenSearch, B_mesons, D_mesons, B_mesons_abs, D_mesons_abs
genSearch = GenSearch()

xSecModifications = {
    # low mass DY k-factor. analysis specific modification, won't change the central sample file.
    'DYJetsToLL_M4to':     {'kFactor': 1.23, 'xsec':None},
    'DYJetsToLL_M5to':     {'kFactor': 1.23, 'xsec':None},
    'DYJetsToLL_M10to':    {'kFactor': 1.23, 'xsec':None},
    }
    
# central configuration
targetLumi = 1000 #pb-1 Which lumi to normalize to

def extractEra(sampleName):
    return sampleName[sampleName.find("Run"):sampleName.find("Run")+len('Run2000A')]

def get_parser():
    ''' Argument parser for post-processing module.
    '''
    import argparse
    argParser = argparse.ArgumentParser(description = "Argument parser for cmgPostProcessing")

    argParser.add_argument('--logLevel',    action='store',         nargs='?',  choices=['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'TRACE', 'NOTSET'],   default='INFO', help="Log level for logging" )
    argParser.add_argument('--samples',     action='store',         nargs='*',  type=str, default=['TTZToLLNuNu_ext'],                  help="List of samples to be post-processed, given as CMG component name" )
    argParser.add_argument('--eventsPerJob',action='store',         nargs='?',  type=int, default=30000000,                             help="Maximum number of events per job (Approximate!)." )
    argParser.add_argument('--nJobs',       action='store',         nargs='?',  type=int, default=1,                                    help="Maximum number of simultaneous jobs." )
    argParser.add_argument('--job',         action='store',                     type=int, default=0,                                    help="Run only jobs i" )
    argParser.add_argument('--minNJobs',    action='store',         nargs='?',  type=int, default=1,                                    help="Minimum number of simultaneous jobs." )
    argParser.add_argument('--targetDir',   action='store',         nargs='?',  type=str, default=user.postprocessing_output_directory, help="Name of the directory the post-processed files will be saved" )
    argParser.add_argument('--processingEra', action='store',       nargs='?',  type=str, default='postProcessed_80X_v22',              help="Name of the processing era" )
    argParser.add_argument('--writeToDPM', action='store_true',                                                                         help="Write output to DPM?")
    argParser.add_argument('--skim',        action='store',         nargs='?',  type=str, default='dilepTiny',                          help="Skim conditions to be applied for post-processing" )
    argParser.add_argument('--massSkim',    action='store',         nargs='?',  type=str, default='',                                   help="Mass skim conditions to be applied for post-processing" )
    argParser.add_argument('--LHEHTCut',    action='store',         nargs='?',  type=int, default=-1,                                   help="LHE cut." )
    argParser.add_argument('--year',        action='store',                     type=int,                                               help="Which year?" )
    argParser.add_argument('--overwriteJEC',action='store',                               default=None,                                 help="Overwrite JEC?" )
    argParser.add_argument('--overwrite',   action='store_true',                                                                        help="Overwrite existing output files, bool flag set to True  if used" )
    argParser.add_argument('--runOnLxPlus', action='store_true',                                                                        help="Change the global redirector of samples to run on lxplus")
    argParser.add_argument('--keepAllJets', action='store_true',                                                                        help="Keep also forward jets?" )
    argParser.add_argument('--small',       action='store_true',                                                                        help="Run the file on a small sample (for test purpose), bool flag set to True if used" )
    argParser.add_argument('--susySignal',  action='store_true',                                                                        help="Is SUSY signal?" )
    argParser.add_argument('--flagTTGamma', action='store_true',                                                                        help="Is ttgamma?" )
    argParser.add_argument('--flagTTBar',   action='store_true',                                                                        help="Is ttbar?" )
    argParser.add_argument('--TTDM',        action='store_true',                                                                        help="Is TTDM signal?" )
    argParser.add_argument('--fastSim',     action='store_true',                                                                        help="FastSim?" )
    argParser.add_argument('--triggerSelection',            action='store_true',                                                        help="Trigger selection?" ) 
    argParser.add_argument('--skipGenLepMatching',          action='store_true',                                                        help="skip matched genleps??" )
    argParser.add_argument('--keepLHEWeights',              action='store_true',                                                        help="Keep LHEWeights?" )
    argParser.add_argument('--checkTTGJetsOverlap',         action='store_true',                                                        help="Keep TTGJetsEventType which can be used to clean TTG events from TTJets samples" )
    argParser.add_argument('--skipSystematicVariations',    action='store_true',                                                        help="Don't calulcate BTag, JES and JER variations.")
    argParser.add_argument('--noTopPtReweighting',          action='store_true',                                                        help="Skip top pt reweighting.")
    argParser.add_argument('--forceProxy',                  action='store_true',                                                        help="Don't check certificate")
    argParser.add_argument('--skipNanoTools',               action='store_true',                                                        help="Skipt the nanoAOD tools step for computing JEC/JER/MET etc uncertainties")
    argParser.add_argument('--keepNanoAOD',                 action='store_true',                                                        help="Keep nanoAOD output?")
    argParser.add_argument('--reuseNanoAOD',                action='store_true',                                                        help="Keep nanoAOD output?")
    argParser.add_argument('--reapplyJECS',                 action='store_true',                                                        help="Reapply JECs to data?")
    argParser.add_argument('--reduceSizeBy',                action='store',     type=int,                                               help="Reduce the size of the sample by a factor of...")
    argParser.add_argument('--event',                       action='store',     type=int, default=-1,                                   help="Just process event no")

    return argParser

options = get_parser().parse_args()

# Logging
import StopsDilepton.tools.logger as _logger
logFile = '/tmp/%s_%s_%s_njob%s.txt'%(options.skim, '_'.join(options.samples), os.environ['USER'], str(0 if options.nJobs==1 else options.job))
logger  = _logger.get_logger(options.logLevel, logFile = logFile)

import Analysis.Tools.logger as _logger_an
logger_an = _logger_an.get_logger('DEBUG', logFile = logFile )

import RootTools.core.logger as _logger_rt
logger_rt = _logger_rt.get_logger(options.logLevel, logFile = logFile )

def fill_vector_collection( event, collection_name, collection_varnames, objects, maxN = 100):
    setattr( event, "n"+collection_name, len(objects) )
    for i_obj, obj in enumerate(objects[:maxN]):
        for var in collection_varnames:
            if var in obj.keys():
                if type(obj[var]) == type("string"):
                    obj[var] = int(ord(obj[var]))
                if type(obj[var]) == type(True):
                    obj[var] = int(obj[var])
                getattr(event, collection_name+"_"+var)[i_obj] = obj[var]

#_logger.   add_fileHandler( user.data_output_directory + '/logs/%s_%s_debug.txt'%(options.samples[0], options.job), options.logLevel )

# Flags 
isDiLep         = options.skim.lower().startswith('dilep')
isTriLep        = options.skim.lower().startswith('trilep')
isSingleLep     = options.skim.lower().startswith('singlelep')
isTiny          = options.skim.lower().count('tiny') 
isSmall         = options.skim.lower().count('small')
isInclusive     = options.skim.lower().count('inclusive') 

# Skim condition
skimConds = []

if options.event > 0:
    skimConds.append( "event == %s"%options.event )

if isDiLep:
    skimConds.append( "Sum$(Electron_pt>20&&abs(Electron_eta)<2.4) + Sum$(Muon_pt>20&&abs(Muon_eta)<2.4)>=2" )
if isTriLep:
    skimConds.append( "Sum$(Electron_pt>20&&abs(Electron_eta)&&Electron_pfRelIso03_all<0.4) + Sum$(Muon_pt>20&&abs(Muon_eta)<2.5&&Muon_pfRelIso03_all<0.4)>=2 && Sum$(Electron_pt>10&&abs(Electron_eta)<2.5)+Sum$(Muon_pt>10&&abs(Muon_eta)<2.5)>=3" )
elif isSingleLep:
    skimConds.append( "Sum$(Electron_pt>20&&abs(Electron_eta)<2.5) + Sum$(Muon_pt>20&&abs(Muon_eta)<2.5)>=1" )

if isInclusive:
    skimConds.append('(1)')
    isSingleLep = True #otherwise no lepton variables?!
    isDiLep     = True
    isTriLep    = True

if options.massSkim:
    stop, lsp = options.massSkim.split(',')
    logger.info("Selecting mass point: %s, %s", stop, lsp)
    skimConds.append( "Max$(GenPart_mass*(abs(GenPart_pdgId)==1000006))=="+stop+"&&Max$(GenPart_mass*(abs(GenPart_pdgId)==1000022))=="+lsp)

print skimConds
#Samples: Load samples
maxN = 1 if options.small else None
if options.small:
    options.job = 0
    options.nJobs = 10000 # set high to just run over 1 input file

if options.runOnLxPlus:
    # Set the redirector in the samples repository to the global redirector
    from Samples.Tools.config import redirector_global as redirector

if options.year == 2016:
    from Samples.nanoAOD.Summer16_private_legacy_v1 import allSamples as mcSamples
    from Samples.nanoAOD.Run2016_17Jul2018_private  import allSamples as dataSamples
    from StopsDilepton.samples.nanoAOD_TTDM_2016    import allSamples as TTDMSamples
    allSamples = mcSamples + dataSamples + TTDMSamples
elif options.year == 2017:
    from Samples.nanoAOD.Fall17_private_legacy_v1   import allSamples as mcSamples
    from Samples.nanoAOD.Run2017_31Mar2018_private  import allSamples as dataSamples
    from StopsDilepton.samples.nanoAOD_TTDM_2017    import allSamples as TTDMSamples
    allSamples = mcSamples + dataSamples + TTDMSamples
elif options.year == 2018:
    from Samples.nanoAOD.Spring18_private           import allSamples as HEMSamples
    from Samples.nanoAOD.Run2018_26Sep2018_private  import allSamples as HEMDataSamples
    from Samples.nanoAOD.Autumn18_private_legacy_v1 import allSamples as mcSamples
    from Samples.nanoAOD.Run2018_17Sep2018_private  import allSamples as dataSamples
    from StopsDilepton.samples.nanoAOD_TTDM_2018    import allSamples as TTDMSamples
    allSamples = HEMSamples + HEMDataSamples + mcSamples + dataSamples + TTDMSamples
else:
    raise NotImplementedError

samples = []
for selectedSamples in options.samples:
    for sample in allSamples:
        if selectedSamples == sample.name:
            samples.append(sample)

if len(samples)==0:
    logger.info( "No samples found. Was looking for %s. Exiting" % options.samples )
    sys.exit(-1)

isData = False not in [s.isData for s in samples]
isMC   =  True not in [s.isData for s in samples]


if options.susySignal:
    xSection = None
    ## special filet for bad jets in FastSim: https://twiki.cern.ch/twiki/bin/viewauth/CMS/SUSRecommendationsICHEP16#Cleaning_up_of_fastsim_jets_from
    #skimConds.append( "Sum$(JetFailId_pt>30&&abs(JetFailId_eta)<2.5&&JetFailId_mcPt==0&&JetFailId_chHEF<0.1)+Sum$(Jet_pt>30&&abs(Jet_eta)<2.5&&Jet_mcPt==0&&Jet_chHEF<0.1)==0" )
else:
    # Check that all samples which are concatenated have the same x-section.
    assert isData or len(set([s.xSection for s in samples]))==1, "Not all samples have the same xSection: %s !"%(",".join([s.name for s in samples]))
    assert isMC or len(samples)==1, "Don't concatenate data samples"

    xSection = samples[0].xSection if isMC else None

## Trigger selection

L1PW = L1PrefireWeight(options.year)

# set default era
era = None

if isData and options.triggerSelection:
    from StopsDilepton.tools.triggerSelector import triggerSelector
    era = extractEra(samples[0].name)[-1]
    logger.info( "######### Era %s ########", era )
    ts = triggerSelector(options.year, era=era)
    triggerCond  = ts.getSelection(options.samples[0] if isData else "MC")
    treeFormulas = {"triggerDecision": {'string':triggerCond} }

    logger.info("Sample will have the following trigger skim: %s"%triggerCond)
    skimConds.append( triggerCond )
elif isData and not options.triggerSelection:
    raise Exception( "Data should have a trigger selection" )

triggerEff          = triggerEfficiency(options.year)
isr                 = ISRweight()

#Samples: combine if more than one
if len(samples)>1:
    sample_name =  samples[0].name+"_comb"
    logger.info( "Combining samples %s to %s.", ",".join(s.name for s in samples), sample_name )
    sample      = Sample.combine(sample_name, samples, maxN = maxN)
    sampleForPU = Sample.combine(sample_name, samples, maxN = -1)
elif len(samples)==1:
    sample      = samples[0]
    sampleForPU = samples[0]
else:
    raise ValueError( "Need at least one sample. Got %r",samples )

# Add scale etc. friends
has_susy_weight_friend = False
if options.susySignal and options.fastSim:
    # Make friend sample
    friend_dir = "/afs/hephy.at/data/cms05/nanoTuples/signalWeights/%s/%s"% (options.year, sample.name )
    if os.path.exists( friend_dir ):
        weight_friend = Sample.fromDirectory( "weight_friend", directory = [friend_dir] ) 
        if weight_friend.chain.BuildIndex("luminosityBlock", "event")>0:
            has_susy_weight_friend = True
    else:
        raise RuntimeError( "We need the LHE weight friend tries. Not found in: %s" % friend_dir )

if options.reduceSizeBy > 1:
    logger.info("Sample size will be reduced by a factor of %s", options.reduceSizeBy)
    logger.info("Recalculating the normalization of the sample. Before: %s", sample.normalization)
    if isData:
        NotImplementedError ( "Data samples shouldn't be reduced in size!!" )
    sample.reduceFiles( factor = options.reduceSizeBy )
    # recompute the normalization
    sample.clear()
    sample.name += "_redBy%s"%options.reduceSizeBy
    sample.normalization = sample.getYieldFromDraw(weightString="genWeight")['val']
    logger.info("New normalization: %s", sample.normalization)

nameForISR = copy.deepcopy(sample.name)
if options.massSkim:
    sample.name += "_%s_%s"%(stop,lsp)

if isMC:
    from Analysis.Tools.puReweighting import getReweightingFunction
    if options.year == 2016:
        nTrueInt_puRW       = getReweightingFunction(data="PU_2016_35920_XSecCentral", mc="Summer16")
        nTrueInt_puRWDown   = getReweightingFunction(data="PU_2016_35920_XSecDown",    mc="Summer16")
        nTrueInt_puRWVDown  = getReweightingFunction(data="PU_2016_35920_XSecVDown",   mc="Summer16")
        nTrueInt_puRWUp     = getReweightingFunction(data="PU_2016_35920_XSecUp",      mc="Summer16")
        nTrueInt_puRWVUp    = getReweightingFunction(data="PU_2016_35920_XSecVUp",     mc="Summer16")
        nTrueInt_puRWVVUp   = getReweightingFunction(data="PU_2016_35920_XSecVVUp",    mc="Summer16")
    elif options.year == 2017:
        # keep the weight name for now. Should we update to a more general one?
        puProfiles = puProfile( source_sample = sampleForPU )
        mcHist = puProfiles.cachedTemplate( selection="( 1 )", weight='genWeight', save=False, overwrite=False ) # use genWeight for amc@NLO samples. No problems encountered so far
        print "1", mcHist, mcHist.Integral()
        nTrueInt_puRW       = getReweightingFunction(data="PU_2017_41860_XSecCentral",  mc=mcHist)
        nTrueInt_puRWDown   = getReweightingFunction(data="PU_2017_41860_XSecDown",     mc=mcHist)
        nTrueInt_puRWVDown  = getReweightingFunction(data="PU_2017_41860_XSecVDown",    mc=mcHist)
        nTrueInt_puRWUp     = getReweightingFunction(data="PU_2017_41860_XSecUp",       mc=mcHist)
        nTrueInt_puRWVUp    = getReweightingFunction(data="PU_2017_41860_XSecVUp",      mc=mcHist)
        nTrueInt_puRWVVUp   = getReweightingFunction(data="PU_2017_41860_XSecVVUp",     mc=mcHist)
    elif options.year == 2018:
        # keep the weight name for now. Should we update to a more general one?
        nTrueInt_puRW       = getReweightingFunction(data="PU_2018_58830_XSecCentral",  mc="Autumn18")
        nTrueInt_puRWDown   = getReweightingFunction(data="PU_2018_58830_XSecDown",     mc="Autumn18")
        nTrueInt_puRWVDown  = getReweightingFunction(data="PU_2018_58830_XSecVDown",    mc="Autumn18")
        nTrueInt_puRWUp     = getReweightingFunction(data="PU_2018_58830_XSecUp",       mc="Autumn18")
        nTrueInt_puRWVUp    = getReweightingFunction(data="PU_2018_58830_XSecVUp",      mc="Autumn18")
        nTrueInt_puRWVVUp   = getReweightingFunction(data="PU_2018_58830_XSecVVUp",     mc="Autumn18")

## lepton SFs

leptonTrackingSF    = LeptonTrackingEfficiency(options.year)
leptonSF            = leptonSF_(options.year)
leptonHit0SF        = leptonHit0SF_(options.year)
leptonSip3dSF       = leptonSip3dSF_(options.year)

if options.fastSim:
   leptonFastSimSF  = leptonFastSimSF_(options.year)

options.skim = options.skim + '_small' if options.small else options.skim

# LHE cut (DY samples)
if options.LHEHTCut>0:
    sample.name+="_lheHT"+str(options.LHEHTCut)
    logger.info( "Adding upper LHE cut at %f", options.LHEHTCut )
    skimConds.append( "LHE_HTIncoming<%f"%options.LHEHTCut )

sampleName = sample.name

# output directory (store temporarily when running on dpm)
if options.writeToDPM:
    from StopsDilepton.tools.user import dpm_directory  as user_directory
    from Samples.Tools.config  import redirector as redirector_hephy
    # Allow parallel processing of N threads on one worker
    output_directory = os.path.join( '/tmp/%s'%os.environ['USER'], str(uuid.uuid4()) )
    targetPath       = redirector_hephy + os.path.join( user_directory, 'postprocessed',  options.processingEra, options.skim, sampleName )
else:
    # User specific
    from StopsDilepton.tools.user import postprocessing_output_directory as user_directory
    directory = os.path.join( options.targetDir, options.processingEra ) 
    output_directory = os.path.join( directory, options.skim, sampleName )

renormISR = False
if options.susySignal:
    from StopsDilepton.samples.helpers import getT2ttSignalWeight, getT2ttISRNorm
    logger.info( "SUSY signal samples to be processed: %s", ",".join(s.name for s in samples) )
    assert len(samples)==1, "Can only process one SUSY sample at a time."
    logger.info( "Signal weights will be drawn from %s files. If that's not the whole sample, stuff will be wrong.", len(samples[0].files))
    logger.info( "Fetching signal weights..." )
    logger.info( "Weights will be stored in %s for future use.", output_directory)
    print "loading weights from here:", os.path.join( user.cache_dir, str(options.year)), samples[0], samples[0].name
    signalWeight = getT2ttSignalWeight( samples[0], lumi = targetLumi, cacheDir = os.path.join( user.cache_dir, str(options.year)) ) #Can use same x-sec/weight for T8bbllnunu as for T2tt
    logger.info("Done fetching signal weights.")

    logger.info("Fetching the normalization for the ISR weights.")
    masspoints = signalWeight.keys()
    if getT2ttISRNorm(sample, masspoints[0][0], masspoints[0][1], masspoints, options.year, signal=nameForISR, cacheDir = os.path.join( user.cache_dir, str(options.year))):
        renormISR = True
        logger.info("Successfully loaded ISR normalzations.")
    else:
        logger.info("!!WARNING!! No ISR normaliztion factors found. Using the ISR weights will therefore change the normalization. Be careful!")
        #raise NotImplementedError ("Couldn't load ISR normalization factors.")

if options.TTDM:
    signalWeight = getTTDMSignalWeight( samples[0], lumi = targetLumi, year=options.year)

#raise NotImplementedError

len_orig = len(sample.files)
## sort the list of files?
sample = sample.split( n=options.nJobs, nSub=options.job)
logger.info( "fileBasedSplitting: Run over %i/%i files for job %i/%i."%(len(sample.files), len_orig, options.job, options.nJobs))
logger.debug( "fileBasedSplitting: Files to be run over:\n%s", "\n".join(sample.files) )

# top pt reweighting
from StopsDilepton.tools.topPtReweighting import getUnscaledTopPairPtReweightungFunction, getTopPtDrawString, getTopPtsForReweighting
# Decision based on sample name -> whether TTJets or TTLep is in the sample name
isTT = sample.name.startswith("TTJets") or sample.name.startswith("TTLep") or sample.name.startswith("TT_pow")
doTopPtReweighting = isTT and not options.noTopPtReweighting

if sample.name.startswith("TTLep"):
    sample.topScaleF = 1.002 ## found to be universal for years 2016-2018, and in principle negligible

if doTopPtReweighting:
    logger.info( "Sample will have top pt reweighting." )
    topPtReweightingFunc = getUnscaledTopPairPtReweightungFunction(selection = "dilep")
    # Compute x-sec scale factor on unweighted events
    selectionString = "&&".join(skimConds)
    if hasattr(sample, "topScaleF"):
        # If you don't want to get the SF for each subjob run the script and add the topScaleF to the sample
        topScaleF = sample.topScaleF
    else:
        reweighted  = sample.getYieldFromDraw( selectionString = selectionString, weightString = getTopPtDrawString(selection = "dilep") + '*genWeight')
        central     = sample.getYieldFromDraw( selectionString = selectionString, weightString = 'genWeight')

        topScaleF = central['val']/reweighted['val']

    logger.info( "Found topScaleF %f", topScaleF )
else:
    topScaleF = 1
    logger.info( "Sample will NOT have top pt reweighting. topScaleF=%f",topScaleF )

if options.susySignal:
    from Analysis.Tools.isrWeight import *
    isr = ISRweight()
    logger.info("Sample will have ISR reweighting called reweight_nISR.")

# systematic variations
logger.info("Adding systematic variations")
addSystematicVariations = (not isData) and (not options.skipSystematicVariations)

# B tagging SF
from Analysis.Tools.BTagEfficiency import BTagEfficiency
btagEff = BTagEfficiency( fastSim = options.fastSim, year=options.year, tagger='DeepCSV' )

# Directory for individual signal files
if options.susySignal:
    signalSubDir = options.samples[0].split('_')[1]
    signalDir = os.path.join(directory, options.skim, signalSubDir)
    logger.info("Separate files for each mass point will be located in %s"%signalDir)
    if not os.path.exists(signalDir): 
        try:
            os.makedirs(signalDir)
        except OSError:
            pass

if os.path.exists(output_directory) and options.overwrite:
    if options.nJobs > 1:
        logger.warning( "NOT removing directory %s because nJobs = %i", output_directory, options.nJobs )
    else:
        logger.info( "Output directory %s exists. Deleting.", output_directory )
        shutil.rmtree(output_directory)

try:    #Avoid trouble with race conditions in multithreading
    os.makedirs(output_directory)
    logger.info( "Created output directory %s.", output_directory )
except:
    pass

logger.info("Checking existing root files.")
filename, ext = os.path.splitext( os.path.join(output_directory, sample.name + '.root') )
fileNumber = options.job if options.job is not None else 0
outfilename = filename+ext
# checking overwrite or file exists
if not options.overwrite and options.writeToDPM:
    try:
        # ls the directory on DPM
        checkFile = "/cms" + targetPath.split("/cms")[1] + "/"
        cmd = [ "xrdfs", redirector_hephy, "ls", checkFile ]
        fileList = subprocess.check_output( cmd ).split("\n")[:-1]
        fileList = [ line.split(checkFile)[1].split(".root")[0] for line in fileList ]
    except:
        # Not even the directory exists on dpm
        fileList = []

    if sample.name in fileList:
        # Sample found on dpm, check if it is ok
        target  = os.path.join( targetPath, sample.name+".root" )
        if checkRootFile( target, checkForObjects=["Events"] ) and deepCheckRootFile( target ) and deepCheckWeight( target ):
            logger.info( "File already processed. Source: File check ok! Skipping." ) # Everything is fine, no overwriting
            sys.exit(0)
        else:
            logger.info( "File corrupt. Removing file from target." )
            cmd = [ "xrdfs", redirector_hephy, "rm", "/cms" + target.split("/cms")[1] ]
            subprocess.call( cmd )
            logger.info( "Reprocessing." )
    else:
        logger.info( "Sample not processed yet." )
        logger.info( "Processing." )

elif not options.overwrite and not options.writeToDPM:
    if os.path.isfile(outfilename):
        logger.info( "Output file %s found.", outfilename)
        if checkRootFile( outfilename, checkForObjects=["Events"] ) and deepCheckRootFile( outfilename ) and deepCheckWeight( outfilename ):
            logger.info( "File already processed. Source: File check ok! Skipping." ) # Everything is fine, no overwriting
            sys.exit(0)
        else:
            logger.info( "File corrupt. Removing file from target." )
            os.remove( outfilename )
            logger.info( "Reprocessing." )
    else:
        logger.info( "Sample not processed yet." )
        logger.info( "Processing." )

else:
    logger.info( "Overwriting.")

logger.info("Defining all the branches.")
#branches to be kept for data and MC
branchKeepStrings_DATAMC = [\
    "run", "luminosityBlock", "event", "fixedGridRhoFastjetAll", "PV_npvs", "PV_npvsGood",
    "MET_*",
    "CaloMET_*",
    "RawMET_phi", "RawMET_pt", "RawMET_sumEt",
    "Flag_*",
    "nJet", "Jet_*",
    "nElectron", "Electron_*",
    "nMuon", "Muon_*",
]
if not options.fastSim:
    branchKeepStrings_DATAMC += ["HLT_*"]

if options.year == 2017:
    branchKeepStrings_DATAMC += [\
        "METFixEE2017_*",
    ]

#branches to be kept for MC samples only
branchKeepStrings_MC = [ "Generator_*", "GenPart_*", "nGenPart", "genWeight", "Pileup_nTrueInt","GenMET_*", "nISR", "nGenJet", "GenJet_*"]
if not options.fastSim:
    branchKeepStrings_MC.extend([ "*LHEScaleWeight", "*LHEPdfWeight"])
#branches to be kept for data only
branchKeepStrings_DATA = [ ]

# Jet variables to be read from chain
jetCorrInfo = []
jetMCInfo   = ['genJetIdx/I','hadronFlavour/I']

if isData:
    lumiScaleFactor=None
    branchKeepStrings = branchKeepStrings_DATAMC + branchKeepStrings_DATA
    from FWCore.PythonUtilities.LumiList import LumiList
    # Apply golden JSON
    lumiList = LumiList(os.path.expandvars(sample.json))
    logger.info( "Loaded json %s", sample.json )
else:
    for x in xSecModifications.keys():
        if sample.name.startswith(x):
            xSectionOld = xSection
            if xSecModifications[x]['xsec']:
                xSection = xSecModifications[x]['xsec']
            if xSecModifications[x]['kFactor']:
                xSection = xSection * xSecModifications[x]['kFactor']
            logger.info("CAREFUL! Sample x-sec will be modified. Old x-sec: %s. New x-sec: %s.", xSectionOld, xSection)
    #if sample.name.startswith('DYJetsToLL_M10to') or sample.name.startswith('DYJetsToLL_M4to') or sample.name.startswith('DYJetsToLL_M5to'):
    #    logger.info("CAREFUL! Detected low mass DY sample, will apply an additional k-factor of 1.23 to x-sec: %s!"%xSection)
    #    lumiScaleFactor = lumiScaleFactor*1.23
    lumiScaleFactor = xSection*targetLumi/float(sample.normalization) if xSection is not None else None
    branchKeepStrings = branchKeepStrings_DATAMC + branchKeepStrings_MC

if options.TTDM:
    branchKeepStrings += ["GenModel*"]

jetVars         = ['pt/F', 'chEmEF/F', 'chHEF/F', 'neEmEF/F', 'neHEF/F', 'rawFactor/F', 'eta/F', 'phi/F', 'jetId/I', 'btagDeepB/F', 'btagCSVV2/F', 'area/F', 'pt_nom/F', 'corr_JER/F'] + jetCorrInfo
if isMC:
    jetVars     += jetMCInfo
    jetVars     += ['pt_jesTotalUp/F', 'pt_jesTotalDown/F', 'pt_jerUp/F', 'pt_jerDown/F', 'corr_JER/F', 'corr_JEC/F']
jetVarNames     = [x.split('/')[0] for x in jetVars]
genLepVars      = ['pt/F', 'phi/F', 'eta/F', 'pdgId/I', 'genPartIdxMother/I', 'status/I', 'statusFlags/I'] # some might have different types
genLepVarNames  = [x.split('/')[0] for x in genLepVars]
# those are for writing leptons
lepVars         = ['pt/F','eta/F','phi/F','pdgId/I','cutBased/I','miniPFRelIso_all/F','pfRelIso03_all/F','sip3d/F','lostHits/I','convVeto/I','dxy/F','dz/F','charge/I','deltaEtaSC/F','mediumId/I','eleIndex/I','muIndex/I']
lepVarNames     = [x.split('/')[0] for x in lepVars]

read_variables = map(TreeVariable.fromString, [ 'MET_pt/F', 'MET_phi/F', 'run/I', 'luminosityBlock/I', 'event/l', 'PV_npvs/I', 'PV_npvsGood/I'] )
if options.year == 2017:
    read_variables += map(TreeVariable.fromString, [ 'METFixEE2017_pt/F', 'METFixEE2017_phi/F', 'METFixEE2017_pt_nom/F', 'METFixEE2017_phi_nom/F'])
    if isMC:
        read_variables += map(TreeVariable.fromString, [ 'METFixEE2017_pt_jesTotalUp/F', 'METFixEE2017_pt_jesTotalDown/F', 'METFixEE2017_pt_jerUp/F', 'METFixEE2017_pt_jerDown/F', 'METFixEE2017_pt_unclustEnDown/F', 'METFixEE2017_pt_unclustEnUp/F', 'METFixEE2017_phi_jesTotalUp/F', 'METFixEE2017_phi_jesTotalDown/F', 'METFixEE2017_phi_jerUp/F', 'METFixEE2017_phi_jerDown/F', 'METFixEE2017_phi_unclustEnDown/F', 'METFixEE2017_phi_unclustEnUp/F', 'METFixEE2017_pt_jer/F', 'METFixEE2017_phi_jer/F'])
else:
    read_variables += map(TreeVariable.fromString, [ 'MET_pt_nom/F', 'MET_phi_nom/F' ])
    if isMC:
        read_variables += map(TreeVariable.fromString, [ 'MET_pt_jesTotalUp/F', 'MET_pt_jesTotalDown/F', 'MET_pt_jerUp/F', 'MET_pt_jerDown/F', 'MET_pt_unclustEnDown/F', 'MET_pt_unclustEnUp/F', 'MET_phi_jesTotalUp/F', 'MET_phi_jesTotalDown/F', 'MET_phi_jerUp/F', 'MET_phi_jerDown/F', 'MET_phi_unclustEnDown/F', 'MET_phi_unclustEnUp/F', 'MET_pt_jer/F', 'MET_phi_jer/F'])
if isMC:
    read_variables += map(TreeVariable.fromString, [ 'GenMET_pt/F', 'GenMET_phi/F' ])

read_variables += [ TreeVariable.fromString('nPhoton/I'),
                    VectorTreeVariable.fromString('Photon[pt/F,eta/F,phi/F,mass/F,cutBased/I,pdgId/I]') if (options.year == 2016) else VectorTreeVariable.fromString('Photon[pt/F,eta/F,phi/F,mass/F,cutBasedBitmap/I,pdgId/I]') ]

new_variables = [ 'weight/F']
if isMC:
    read_variables += [ TreeVariable.fromString('Pileup_nTrueInt/F') ]
    # reading gen particles for top pt reweighting
    read_variables.append( TreeVariable.fromString('nGenPart/I') )
    read_variables.append( TreeVariable.fromString('nISR/I') )
    read_variables.append( VectorTreeVariable.fromString('GenPart[pt/F,mass/F,phi/F,eta/F,pdgId/I,genPartIdxMother/I,status/I,statusFlags/I]', nMax=200 )) # default nMax is 100, which would lead to corrupt values in this case
    read_variables.append( TreeVariable.fromString('genWeight/F') )
    read_variables.append( TreeVariable.fromString('nGenJet/I') )
    read_variables.append( VectorTreeVariable.fromString('GenJet[pt/F,eta/F,phi/F]' ) )
    new_variables.extend([ 'reweightTopPt/F', 'reweight_nISR/F', 'reweight_nISRUp/F', 'reweight_nISRDown/F', 'reweightPU/F','reweightPUUp/F','reweightPUDown/F', 'reweightPUVUp/F','reweightPUVVUp/F', 'reweightPUVDown/F', 'reweightL1Prefire/F', 'reweightL1PrefireUp/F', 'reweightL1PrefireDown/F'])
    if not options.skipGenLepMatching:
        TreeVariable.fromString( 'nGenLep/I' ),
        new_variables.append( 'GenLep[%s]'% ( ','.join(genLepVars) ) )

read_variables += [\
    TreeVariable.fromString('nElectron/I'),
    VectorTreeVariable.fromString('Electron[pt/F,eta/F,phi/F,pdgId/I,cutBased/I,miniPFRelIso_all/F,pfRelIso03_all/F,sip3d/F,lostHits/b,convVeto/O,dxy/F,dz/F,charge/I,deltaEtaSC/F,vidNestedWPBitmap/I{mcString}]'.format(mcString=",genPartFlav/B,genPartIdx/I" if isMC else "")),
    TreeVariable.fromString('nMuon/I'),
    VectorTreeVariable.fromString('Muon[pt/F,eta/F,phi/F,pdgId/I,mediumId/O,miniPFRelIso_all/F,pfRelIso03_all/F,sip3d/F,dxy/F,dz/F,charge/I{mcString}]'.format(mcString=",genPartFlav/B,genPartIdx/I" if isMC else "")),
    TreeVariable.fromString('nJet/I'),
    VectorTreeVariable.fromString('Jet[%s]'% ( ','.join(jetVars) ) ),
]

new_variables += [\
    'overlapRemoval/I','nlep/I', "min_dl_mass/F",
    'JetGood[%s]'% ( ','.join(jetVars+['index/I']) + ',genPt/F' ),
    'met_pt/F', 'met_phi/F', 'met_pt_min/F'
]

# Add weight branches for susy signal samples from friend tree
if has_susy_weight_friend:
    new_variables.extend([ "LHE[weight/F]", "LHE_weight_original/F"] )

if sample.isData: new_variables.extend( ['jsonPassed/I','isData/I'] )
new_variables.extend( ['nBTag/I', 'ht/F', 'metSig/F'] )
new_variables += ["reweightHEM/F"]

new_variables.append( 'lep[%s]'% ( ','.join(lepVars) ) )

if isSingleLep:
    new_variables.extend( ['m3/F', 'm3_ind1/I', 'm3_ind2/I', 'm3_ind3/I'] )
if isTriLep or isDiLep or isSingleLep:
    new_variables.extend( ['nGoodMuons/I', 'nGoodElectrons/I', 'nGoodLeptons/I' ] )
    new_variables.extend( ['l1_pt/F', 'l1_eta/F', 'l1_phi/F', 'l1_pdgId/I', 'l1_index/I', 'l1_jetPtRelv2/F', 'l1_jetPtRatiov2/F', 'l1_miniRelIso/F', 'l1_relIso03/F', 'l1_dxy/F', 'l1_dz/F', 'l1_mIsoWP/I', 'l1_eleIndex/I', 'l1_muIndex/I' ] )
    new_variables.extend( ['mlmZ_mass/F'])
    if isMC:
            new_variables.extend(['l1_genPartFlav/B', 'l1_genPartIdx/I']) 
            new_variables.extend(['reweightLeptonSF/F', 'reweightLeptonSFUp/F', 'reweightLeptonSFDown/F'])
            new_variables.extend(['reweightLeptonHit0SF/F', 'reweightLeptonHit0SFUp/F', 'reweightLeptonHit0SFDown/F'])
            new_variables.extend(['reweightLeptonSip3dSF/F', 'reweightLeptonSip3dSFUp/F', 'reweightLeptonSip3dSFDown/F'])
if isTriLep or isDiLep:
    new_variables.extend( ['l2_pt/F', 'l2_eta/F', 'l2_phi/F', 'l2_pdgId/I', 'l2_index/I', 'l2_jetPtRelv2/F', 'l2_jetPtRatiov2/F', 'l2_miniRelIso/F', 'l2_relIso03/F', 'l2_dxy/F', 'l2_dz/F', 'l2_mIsoWP/I', 'l2_eleIndex/I', 'l2_muIndex/I' ] )
    new_variables.extend( ['isEE/I', 'isMuMu/I', 'isEMu/I', 'isOS/I' ] )
    new_variables.extend( ['dl_pt/F', 'dl_eta/F', 'dl_phi/F', 'dl_mass/F'])
    new_variables.extend( ['dl_mt2ll/F', 'dl_mt2bb/F', 'dl_mt2blbl/F' ] )
    if isMC: new_variables.extend( \
        [   'l2_genPartFlav/B', 'l2_genPartIdx/I', 'zBoson_genPt/F', 'zBoson_genEta/F', 
            'reweightDilepTrigger/F', 'reweightDilepTriggerUp/F', 'reweightDilepTriggerDown/F',
            'reweightLeptonTrackingSF/F',
         ] )
    if options.susySignal and options.fastSim:
        new_variables.extend( ['dl_mt2ll_gen/F', 'dl_mt2bb_gen/F', 'dl_mt2blbl_gen/F' ] )
new_variables.extend( ['nPhotonGood/I','photon_pt/F','photon_eta/F','photon_phi/F','photon_idCutBased/I'] )
if isMC: new_variables.extend( ['photon_genPt/F', 'photon_genEta/F', 'genZ_mass/F', 'isOnShellTTZ/I'] )
new_variables.extend( ['photonJetdR/F','photonLepdR/F'] )
if isTriLep or isDiLep:
  new_variables.extend( ['dlg_mass/F' ] )

## ttZ related variables
new_variables.extend( ['Z1_l1_index/I', 'Z1_l2_index/I', 'Z2_l1_index/I', 'Z2_l2_index/I', 'nonZ1_l1_index/I', 'nonZ1_l2_index/I'] )
for i in [1,2]:
    new_variables.extend( ['Z%i_pt/F'%i, 'Z%i_eta/F'%i, 'Z%i_phi/F'%i, 'Z%i_lldPhi/F'%i, 'Z%i_lldR/F'%i,  'Z%i_mass/F'%i] )

if options.checkTTGJetsOverlap:
    new_variables.extend( ['TTGJetsEventType/I'] )

if addSystematicVariations:

    for var in ['jesTotalUp', 'jesTotalDown', 'jerUp', 'jer', 'jerDown', 'unclustEnUp', 'unclustEnDown']:
        if not var.startswith('unclust'):
            new_variables.extend( ['nJetGood_'+var+'/I', 'nBTag_'+var+'/I'] )
        new_variables.extend( ['met_pt_'+var+'/F', 'met_phi_'+var+'/F'] )
        if isTriLep or isDiLep:
            new_variables.extend( ['dl_mt2ll_'+var+'/F', 'dl_mt2bb_'+var+'/F', 'dl_mt2blbl_'+var+'/F'] )

# Btag weights Method 1a
for var in btagEff.btagWeightNames:
    if var!='MC':
        new_variables.append('reweightBTag_'+var+'/F')

#if options.susySignal or options.TTDM:
#    read_variables += map(TreeVariable.fromString, ['met_genPt/F', 'met_genPhi/F'] )
if options.susySignal:
    new_variables  += ['reweightXSecUp/F', 'reweightXSecDown/F', 'mStop/I', 'mNeu/I']
    if  'T8bbllnunu' in options.samples[0]:
        new_variables  += ['mCha/I', 'mSlep/I', 'sleptonPdg/I']
    if 'T2tt' in options.samples[0]:
        new_variables  += ['weight_pol_L/F', 'weight_pol_R/F']
if options.TTDM:
    new_variables  += ['reweightXSecUp/F', 'reweightXSecDown/F', 'mPhi/I', 'mChi/I', 'scalar/I', 'pseudoscalar/I']

#myBranchNames = [\
#'GenModel__TTbarDMJets_Dilepton_pseudoscalar_LO_Mchi_1_Mphi_100_TuneCP5_13TeV_madgraph_mcatnlo_pythia8/O',\
#'GenModel__TTbarDMJets_Dilepton_pseudoscalar_LO_Mchi_1_Mphi_150_TuneCP5_13TeV_madgraph_mcatnlo_pythia8/O',\
#'GenModel__TTbarDMJets_Dilepton_pseudoscalar_LO_Mchi_1_Mphi_200_TuneCP5_13TeV_madgraph_mcatnlo_pythia8/O',\
#'GenModel__TTbarDMJets_Dilepton_pseudoscalar_LO_Mchi_1_Mphi_250_TuneCP5_13TeV_madgraph_mcatnlo_pythia8/O',\
#'GenModel__TTbarDMJets_Dilepton_pseudoscalar_LO_Mchi_1_Mphi_300_TuneCP5_13TeV_madgraph_mcatnlo_pythia8/O',\
#'GenModel__TTbarDMJets_Dilepton_pseudoscalar_LO_Mchi_1_Mphi_350_TuneCP5_13TeV_madgraph_mcatnlo_pythia8/O',\
#'GenModel__TTbarDMJets_Dilepton_pseudoscalar_LO_Mchi_1_Mphi_400_TuneCP5_13TeV_madgraph_mcatnlo_pythia8/O',\
#'GenModel__TTbarDMJets_Dilepton_pseudoscalar_LO_Mchi_1_Mphi_450_TuneCP5_13TeV_madgraph_mcatnlo_pythia8/O',\
#'GenModel__TTbarDMJets_Dilepton_pseudoscalar_LO_Mchi_1_Mphi_500_TuneCP5_13TeV_madgraph_mcatnlo_pythia8/O',\
#'GenModel__TTbarDMJets_Dilepton_pseudoscalar_LO_Mchi_1_Mphi_50_TuneCP5_13TeV_madgraph_mcatnlo_pythia8/O',\
#'GenModel__TTbarDMJets_Dilepton_pseudoscalar_LO_Mchi_20_Mphi_100_TuneCP5_13TeV_madgraph_mcatnlo_pythia8/O',\
#'GenModel__TTbarDMJets_Dilepton_pseudoscalar_LO_Mchi_30_Mphi_100_TuneCP5_13TeV_madgraph_mcatnlo_pythia8/O',\
#'GenModel__TTbarDMJets_Dilepton_pseudoscalar_LO_Mchi_40_Mphi_100_TuneCP5_13TeV_madgraph_mcatnlo_pythia8/O',\
#'GenModel__TTbarDMJets_Dilepton_pseudoscalar_LO_Mchi_45_Mphi_100_TuneCP5_13TeV_madgraph_mcatnlo_pythia8/O',\
#'GenModel__TTbarDMJets_Dilepton_pseudoscalar_LO_Mchi_49_Mphi_100_TuneCP5_13TeV_madgraph_mcatnlo_pythia8/O',\
#'GenModel__TTbarDMJets_Dilepton_pseudoscalar_LO_Mchi_51_Mphi_100_TuneCP5_13TeV_madgraph_mcatnlo_pythia8/O',\
#'GenModel__TTbarDMJets_Dilepton_pseudoscalar_LO_Mchi_55_Mphi_100_TuneCP5_13TeV_madgraph_mcatnlo_pythia8/O',\
#]

if options.samples[0].count('_pseudoscalar'):
    print getTTDMBranchNames('pseudoscalar', year=options.year)
    read_variables += map(TreeVariable.fromString, ['%s/O'%x for x in getTTDMBranchNames('pseudoscalar', year=options.year) ] ) # try b or O
    #read_variables = map(TreeVariable.fromString, [ 'MET_pt/F', 'MET_phi/F', 'run/I', 'luminosityBlock/I', 'event/l', 'PV_npvs/I', 'PV_npvsGood/I'] )
if options.samples[0].count('_scalar'):
    print getTTDMBranchNames('scalar', year=options.year)
    read_variables += map(TreeVariable.fromString, ['%s/O'%x for x in getTTDMBranchNames('scalar', year=options.year) ] )


if options.fastSim and (isTriLep or isDiLep):
    new_variables  += ['reweightLeptonFastSimSF/F', 'reweightLeptonFastSimSFUp/F', 'reweightLeptonFastSimSFDown/F']

logger.info("Next up: nanoAOD-tools.")
if not options.skipNanoTools:
    ### nanoAOD postprocessor
    from importlib import import_module
    from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor   import PostProcessor
    from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel       import Collection
    from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop       import Module
    
    ## modules for nanoAOD postprocessor
    #from PhysicsTools.NanoAODTools.postprocessing.modules.jme.jetmetUncertainties   import jetmetUncertaintiesProducer
    #from PhysicsTools.NanoAODTools.postprocessing.modules.jme.jetRecalib            import jetRecalib
    from PhysicsTools.NanoAODTools.postprocessing.modules.jme.METSigProducer        import METSigProducer 
    from PhysicsTools.NanoAODTools.postprocessing.modules.common.ISRcounter        import ISRcounter
    
    logger.info("Preparing nanoAOD postprocessing")
    logger.info("Will put files into directory %s", output_directory)
    cut = '&&'.join(skimConds)
    # different different METSig parameters for data/MC and years
    if options.year == 2016:
        metSigParamsMC      = [1.617529475909303, 1.617529475909303, 1.4505983036866312, 1.4505983036866312, 1.411498565372343, 1.411498565372343, 1.4087559908291813, 1.4087559908291813, 1.3633674107893856, 1.3633674107893856, 0.0019861227075085516, 0.6539410816436597]
        metSigParamsData    = [1.843242937068234, 1.843242937068234, 1.64107911184195,   1.64107911184195,   1.567040591823117, 1.567040591823117, 1.5077143780804294, 1.5077143780804294, 1.614014783345394,  1.614014783345394, -0.0005986196920895609, 0.6071479349467596]
        JER                 = "Summer16_25nsV1_MC"          if not sample.isData else "Summer16_25nsV1_DATA"
        JERera              = "Summer16_25nsV1"

    elif options.year == 2017:
        metSigParamsMC      = [1.9648214119268503, 1.5343086462230238, 1.9167197601498538, 1.5145044341064964, 1.8069380221985405, 1.3217263662622654, 1.5506294867561126, 1.272977540964842,  1.50742322311234,   1.6542883449796797, -0.0017865650107230548,  0.6593106706741719]
        metSigParamsData    = [2.228118299837604,  1.2420725475347338, 2.227630982417529,  1.256752205787215,  2.0215250734187853, 1.1557507029911258, 1.7350536144535336, 1.1587692458345757, 1.9385081854607988, 1.8726188460472792, -2.6697894266706265e-05, 0.646984812801919]
        JER                 = "Fall17_V3_MC"                if not sample.isData else "Fall17_V3_DATA"
        JERera              = "Fall17_V3"

    elif options.year == 2018:
        metSigParamsMC      = [1.9033100259447273, 1.7374326087706358, 1.6957838005481387, 1.7283187962364615, 1.5868244614361302, 1.526252837049335, 1.3744055574137417, 1.4500298644941831, 1.4796204632654997, 1.4481227819959115, 0.0019899110367503207, 0.6927496536100137]
        metSigParamsData    = [1.7492714572981323, 1.3430198915313956, 1.911348554103867, 1.3718438058490257, 1.5885672661901442, 1.4385903478138795, 1.521901070409261, 1.4522895772008289, 1.8870084799263003, 1.7138750357657668, -1.359837542505224e-06, 0.7434240965656385]
        JER                 = "Autumn18_V1_MC"                if not sample.isData else "Autumn18_V1_DATA"
        JERera              = "Autumn18_V1"

    # set the params for MET Significance calculation
    metSigParams            = metSigParamsMC                if not sample.isData else metSigParamsData

    if options.overwriteJEC is not None:
        JEC = options.overwriteJEC

    logger.info("Using JERs for MET significance: %s", JER)
    logger.info("Will use the following parameters for MET significance: %s", metSigParams)
    
    from PhysicsTools.NanoAODTools.postprocessing.modules.jme.jetmetHelperRun2 import *
    METBranchName = 'MET' if not options.year == 2017 else 'METFixEE2017'

    # check if files are available (e.g. if dpm is broken this should result in an error)
    for f in sample.files:
        if not checkRootFile(f):
            raise IOError ("File %s not available"%f)

    # remove empty files. this is necessary in 2018 because empty miniAOD files exist.
    sample.files = [ f for f in sample.files if nonEmptyFile(f) ]
    newFileList = []
    logger.info("Starting nanoAOD postprocessing")
    for f in sample.files:
        JMECorrector = createJMECorrector(isMC=(not sample.isData), dataYear=options.year, runPeriod=era, jesUncert="Total", jetType = "AK4PFchs", metBranchName=METBranchName, isFastSim=options.fastSim, applySmearing=False)
        modules = [
            JMECorrector()
        ]
        
        if not sample.isData:
            modules.append( ISRcounter() )

        if options.year == 2016:
            modules.append( METSigProducer(JER, metSigParams, METCollection="MET", useRecorr=True, calcVariations=(not isData), jetThreshold=15.) )
            if not isData: modules.append( METSigProducer(JER, metSigParams, METCollection="GenMET", useRecorr=False, calcVariations=False, jetThreshold=15., outputBranch="GenMET") )
        elif options.year == 2017:
            modules.append( METSigProducer(JER, metSigParams, METCollection="METFixEE2017", useRecorr=True, calcVariations=(not isData), jetThreshold=15., vetoEtaRegion=(2.65,3.14)) )
            if not isData: modules.append( METSigProducer(JER, metSigParams, METCollection="GenMET", useRecorr=False, calcVariations=False, jetThreshold=15., vetoEtaRegion=(2.65,3.14),  outputBranch="GenMET") )
        elif options.year == 2018:
            modules.append( METSigProducer(JER, metSigParams, METCollection="MET", useRecorr=True, calcVariations=(not isData), jetThreshold=25.) )
            if not isData: modules.append( METSigProducer(JER, metSigParams, METCollection="GenMET", useRecorr=False, calcVariations=False, jetThreshold=25., outputBranch="GenMET") )
        # need a hash to avoid data loss
        file_hash = str(hash(f))
        p = PostProcessor(output_directory, [f], cut=cut, modules=modules, postfix="_for_%s_%s"%(sample.name, file_hash))
        if not options.reuseNanoAOD:
            p.run()
        newFileList += [output_directory + '/' + f.split('/')[-1].replace('.root', '_for_%s_%s.root'%(sample.name, file_hash))]
    logger.info("Done. Replacing input files for further processing.")
    
    sample.files = newFileList

# Define a reader

reader = sample.treeReader( \
    variables = read_variables ,
    selectionString = "&&".join(skimConds)
    )

# using miniRelIso 0.2 as baseline 
eleSelector_ = eleSelector( "tightMiniIso02", year = options.year )
muSelector_  = muonSelector("tightMiniIso02", year = options.year )

genPhotonSel_TTG_OR = genPhotonSelector( 'overlapTTGamma' )

mothers = {"D":0, "B":0}
grannies_D = {}
grannies_B = {}

def filler( event ):
    # shortcut
    r = reader.event
    #workaround  = (r.run, r.luminosityBlock, r.event) # some fastsim files seem to have issues, apparently solved by this.
    event.isData = s.isData
    event.overlapRemoval = 1 
    
    if isMC:

        ## overlap removal by lukas ##
        # GEN Particles
        logger.debug("Getting GenParts")
        gPart = getGenPartsAll(r)
        # GEN Jets
        logger.debug("Getting GenJets")
        gJets = getJets( r, jetVars=['pt','eta','phi','mass','partonFlavour','hadronFlavour','index'], jetColl="GenJet" )

        # Overlap removal flags for ttgamma/ttbar and Zgamma/DY
        logger.debug("GenPhotons")
        GenPhoton                  = filterGenPhotons( gPart, status='last' )

        # OR ttgamma/tt, DY/ZG, WG/WJets
        GenIsoPhoton               = filter( lambda g: isIsolatedPhoton( g, gPart, coneSize=0.2,  ptCut=5, excludedPdgIds=[12,-12,14,-14,16,-16] ), GenPhoton    )
        GenIsoPhotonNoMeson        = filter( lambda g: not hasMesonMother( getParentIds( g, gPart ) ), GenIsoPhoton )

        # OR singleT/tG
        GenIsoPhotonSingleT        = filter( lambda g: isIsolatedPhoton( g, gPart, coneSize=0.05, ptCut=5, excludedPdgIds=[12,-12,14,-14,16,-16] ), GenPhoton    )
        GenIsoPhotonNoMesonSingleT = filter( lambda g: not hasMesonMother( getParentIds( g, gPart ) ), GenIsoPhotonSingleT )
        GenIsoPhotonNoMesonSingleT = filter( lambda g: not photonFromTopDecay( getParentIds( g, gPart ) ), GenIsoPhotonNoMesonSingleT )

        event.isTTGamma = len( filter( lambda g: genPhotonSel_TTG_OR(g), GenIsoPhotonNoMeson        ) ) > 0
        #event.isZWGamma = len( filter( lambda g: genPhotonSel_ZG_OR(g),  GenIsoPhotonNoMeson        ) ) > 0
        #event.isTGamma = len( filter( lambda g: genPhotonSel_T_OR(g), GenIsoPhotonNoMesonSingleT ) ) > 0 

        # new OR flag: Apply overlap removal directly in pp to better handle the plots
        if options.flagTTGamma:
            event.overlapRemoval = event.isTTGamma #good TTgamma event
        elif options.flagTTBar:
            event.overlapRemoval = not event.isTTGamma #good TTbar event


        genLepsFromZ    = genLepFromZ(gPart)
        genZs           = getSortedZCandidates(genLepsFromZ)
        if len(genZs)>0:
            event.genZ_mass = genZs[0][0]
            event.isOnShellTTZ = event.genZ_mass > 70
        
    # weight
    if options.susySignal:
        if has_susy_weight_friend:
            if weight_friend.chain.GetEntryWithIndex(r.luminosityBlock, r.event)>0:
                event.LHE_weight_original =  weight_friend.chain.GetLeaf("LHE_weight_original").GetValue()
                event.nLHE = int(weight_friend.chain.GetLeaf("nLHE").GetValue())
                for nEvt in range(event.nLHE):
                    event.LHE_weight[nEvt] = weight_friend.chain.GetLeaf("LHE_weight").GetValue(nEvt)

        r.GenSusyMStop = max([p['mass']*(abs(p['pdgId']==1000006)) for p in gPart])
        r.GenSusyMNeutralino = max([p['mass']*(abs(p['pdgId']==1000022)) for p in gPart])
        if 'T8bbllnunu' in options.samples[0]:
            r.GenSusyMChargino = max([p['mass']*(abs(p['pdgId']==1000024)) for p in gPart])
            r.GenSusyMSlepton = max([p['mass']*(abs(p['pdgId']==1000011)) for p in gPart]) #FIXME check PDG ID of slepton in sample
            logger.debug("Slepton is selectron with mass %i", r.GenSusyMSlepton)
            event.sleptonPdg = 1000011
            if not r.GenSusyMSlepton > 0:
                r.GenSusyMSlepton = max([p['mass']*(abs(p['pdgId']==1000013)) for p in gPart])
                logger.debug("Slepton is smuon with mass %i", r.GenSusyMSlepton)
                event.sleptonPdg = 1000013
            if not r.GenSusyMSlepton > 0:
                r.GenSusyMSlepton = max([p['mass']*(abs(p['pdgId']==1000015)) for p in gPart])
                logger.debug("Slepton is stau with mass %i", r.GenSusyMSlepton)
                event.sleptonPdg = 1000015
            event.mCha  = int(round(r.GenSusyMChargino,0))
            event.mSlep = int(round(r.GenSusyMSlepton,0))
        #if 'T2tt' in options.samples[0]:
        #    pol_weights = getPolWeights(r)
        #    event.weight_pol_L = pol_weights[0]
        #    event.weight_pol_R = pol_weights[1]

        try:
            event.weight=signalWeight[(int(r.GenSusyMStop), int(r.GenSusyMNeutralino))]['weight'] #* r.genWeight
            #print 'weight', signalWeight[(int(r.GenSusyMStop), int(r.GenSusyMNeutralino))]['weight']
        except KeyError:
            logger.info("Couldn't find weight for %s, %s. Setting weight to 0.", r.GenSusyMStop, r.GenSusyMNeutralino)
            event.weight = 0.
        event.mStop = int(r.GenSusyMStop)
        event.mNeu  = int(r.GenSusyMNeutralino)
        try:
            event.reweightXSecUp    = signalWeight[(r.GenSusyMStop, r.GenSusyMNeutralino)]['xSecFacUp']
            event.reweightXSecDown  = signalWeight[(r.GenSusyMStop, r.GenSusyMNeutralino)]['xSecFacDown']
        except KeyError:
            logger.info("Couldn't find weight for %s, %s. Setting weight to 0.", r.GenSusyMStop, r.GenSusyMNeutralino)
            event.reweightXSecUp    = 0.
            event.reweightXSecDown  = 0.
    elif options.TTDM:
        ttdmWeights             = getTTDMSignalWeightForEvent( sample, r, signalWeight )
        event.weight            = ttdmWeights['weight'] * getattr(r, 'genWeight')
        event.reweightXSecUp    = ttdmWeights['xSecFacUp']
        event.reweightXSecDown  = ttdmWeights['xSecFacDown']
        event.mChi              = ttdmWeights['mChi']
        event.mPhi              = ttdmWeights['mPhi']
        event.scalar            = ttdmWeights['spin'] == 'scalar'
        event.pseudoscalar      = ttdmWeights['spin'] == 'pseudoscalar'

    elif isMC:
        if hasattr(r, "genWeight"):
            event.weight = lumiScaleFactor*r.genWeight if lumiScaleFactor is not None else 1
        else:
            event.weight = lumiScaleFactor if lumiScaleFactor is not None else 1
    elif sample.isData:
        event.weight = 1
    else:
        raise NotImplementedError( "isMC %r isData %r susySignal? %r TTDM? %r" % (isMC, isData, options.susySignal, options.TTDM) )


    # lumi lists and vetos
    if sample.isData:
        #event.vetoPassed  = vetoList.passesVeto(r.run, r.lumi, r.evt)
        event.jsonPassed  = lumiList.contains(r.run, r.luminosityBlock)
        # store decision to use after filler has been executed
        event.jsonPassed_ = event.jsonPassed

    if isMC and hasattr(r, "Pileup_nTrueInt"):
        event.reweightPU     = nTrueInt_puRW       ( r.Pileup_nTrueInt ) # is this correct?
        event.reweightPUDown = nTrueInt_puRWDown   ( r.Pileup_nTrueInt )
        event.reweightPUVDown= nTrueInt_puRWVDown  ( r.Pileup_nTrueInt )
        event.reweightPUUp   = nTrueInt_puRWUp     ( r.Pileup_nTrueInt )
        event.reweightPUVUp  = nTrueInt_puRWVUp    ( r.Pileup_nTrueInt )
        event.reweightPUVVUp = nTrueInt_puRWVVUp   ( r.Pileup_nTrueInt )

    # top pt reweighting
    if isMC:
        event.reweightTopPt     = topPtReweightingFunc(getTopPtsForReweighting(r)) * topScaleF if doTopPtReweighting else 1.
        ISRnorm = getT2ttISRNorm(samples[0], r.GenSusyMStop, r.GenSusyMNeutralino, masspoints, options.year, signal=nameForISR, cacheDir=os.path.join( user.cache_dir, str(options.year))) if renormISR else 1
        event.reweight_nISR     = isr.getWeight(r, norm=ISRnorm )             if options.susySignal else 1
        event.reweight_nISRUp   = isr.getWeight(r, norm=ISRnorm, sigma=1)     if options.susySignal else 1
        event.reweight_nISRDown = isr.getWeight(r, norm=ISRnorm, sigma=-1)    if options.susySignal else 1

    if options.keepAllJets:
        jetAbsEtaCut = 99.
    else:
        jetAbsEtaCut = 2.4
    
    allSlimmedJets      = getJets(r)
    allSlimmedPhotons   = getPhotons(r, year=options.year)
    if options.year == 2018:
        event.reweightL1Prefire, event.reweightL1PrefireUp, event.reweightL1PrefireDown = 1., 1., 1.
    else:
        event.reweightL1Prefire, event.reweightL1PrefireUp, event.reweightL1PrefireDown = L1PW.getWeight(allSlimmedPhotons, allSlimmedJets)

    # get leptons before jets in order to clean jets
    electrons_pt10  = getGoodElectrons(r, ele_selector = eleSelector_)
    muons_pt10      = getGoodMuons(r,     mu_selector = muSelector_ )

    for e in electrons_pt10:
        e['pdgId']      = int( -11*e['charge'] )
        e['eleIndex']   = e['index']
        e['muIndex']    = -1
    for m in muons_pt10:
        m['pdgId']      = int( -13*m['charge'] )
        m['muIndex']    = m['index']
        m['eleIndex']   = -1

    leptons_pt10 = electrons_pt10+muons_pt10

    leptons_pt10.sort(key = lambda p:-p['pt'])

    for iLep, lep in enumerate(leptons_pt10):
        lep['index'] = iLep

    fill_vector_collection( event, "lep", lepVarNames, leptons_pt10)
    event.nlep = len(leptons_pt10)


    leptons      = filter(lambda l:l['pt']>20, leptons_pt10)
    leptons.sort(key = lambda p:-p['pt'])


    # now get jets, cleaned against good leptons

    jetPtVar = 'pt_nom' # see comment below

    # with the latest change, getAllJets calculates the correct jet pt (removing JER) and stores it as Jet_pt again. No need for Jet_pt_nom anymore
    allJetsNotClean = getAllJets(r, [], ptCut=0, absEtaCut=99, jetVars=jetVarNames, jetCollections=["Jet"], idVar=None)
    reallyAllJets= getAllJets(r, leptons, ptCut=0, absEtaCut=99, jetVars=jetVarNames, jetCollections=["Jet"], idVar='jetId') # keeping robert's comment: ... yeah, I know.
    allJets      = filter(lambda j:abs(j['eta'])<jetAbsEtaCut, reallyAllJets)
    jets         = filter(lambda j:jetId(j, ptCut=30, absEtaCut=jetAbsEtaCut, ptVar=jetPtVar), allJets)
    soft_jets    = filter(lambda j:jetId(j, ptCut=0,  absEtaCut=jetAbsEtaCut) and j['pt']<30., allJets) if options.keepAllJets else []
    bJets        = filter(lambda j:isBJet(j, tagger="DeepCSV", year=options.year) and abs(j['eta'])<=2.4    , jets)
    nonBJets     = filter(lambda j:not ( isBJet(j, tagger="DeepCSV", year=options.year) and abs(j['eta'])<=2.4 ), jets)

    nHEMJets = len(filter( lambda j:j['pt']>20 and j['eta']>-3.2 and j['eta']<-1.0 and j['phi']>-2.0 and j['phi']<-0.5, allJets ))

    if isData:
        event.reweightHEM = (r.run>=319077 and nHEMJets==0) or r.run<319077
    else:
        event.reweightHEM = 1 if (nHEMJets==0 or options.year<2018 ) else 0.3518 # 0.2% of Run2018B are HEM affected. Ignore that piece. Thus, if there is a HEM jet, scale the MC to 35.2% which is AB/ABCD=(14.00+7.10)/59.97

    # store the correct MET (EE Fix for 2017, MET_min as backup in 2017)
    
    if options.year == 2017:# and not options.fastSim:
        # v2 recipe. Could also use our own recipe
        event.met_pt    = r.METFixEE2017_pt_nom
        event.met_phi   = r.METFixEE2017_phi_nom
        #event.met_pt_min = r.MET_pt_min not done anymore
    else:
        event.met_pt    = r.MET_pt_nom 
        event.met_phi   = r.MET_phi_nom

        event.met_pt_min = 0

    # Filling jets
    maxNJet = 100
    store_jets = jets if not options.keepAllJets else soft_jets + jets
    store_jets = store_jets[:maxNJet]
    store_jets.sort( key = lambda j:-j[jetPtVar])
    event.nJetGood   = len(store_jets)
    for iJet, jet in enumerate(store_jets):
        event.JetGood_index[iJet] = jet['index']
        for b in jetVarNames:
            getattr(event, "JetGood_"+b)[iJet] = jet[b]
        if isMC:
            if store_jets[iJet]['genJetIdx'] >= 0:
                if r.nGenJet<maxNJet:
                    try:
                        event.JetGood_genPt[iJet] = r.GenJet_pt[store_jets[iJet]['genJetIdx']]
                    except IndexError:
                        event.JetGood_genPt[iJet] = -1
                else:
                    event.JetGood_genPt[iJet] = -1
        getattr(event, "JetGood_pt")[iJet] = jet['pt']

    if isSingleLep:
        # Compute M3 and the three indiced of the jets entering m3
        event.m3, event.m3_ind1, event.m3_ind2, event.m3_ind3 = m3( jets )

    event.ht         = sum([j[jetPtVar] for j in jets])
    event.metSig     = event.met_pt/sqrt(event.ht) if event.ht>0 else float('nan')
    event.nBTag      = len(bJets)

    alljets_sys   = {}
    jets_sys      = {}
    bjets_sys     = {}
    nonBjets_sys  = {}


    # Keep photons and estimate met including (leading pt) photon
    photons = getGoodPhotons(r, ptCut=20, idLevel="tight", isData=isData, year=options.year)
    event.nPhotonGood = len(photons)
    if event.nPhotonGood > 0:
      event.photon_pt         = photons[0]['pt']
      event.photon_eta        = photons[0]['eta']
      event.photon_phi        = photons[0]['phi']
      event.photon_idCutBased = photons[0]['cutBased'] if (options.year == 2016) else photons[0]['cutBasedBitmap']
      if isMC:
        genPhoton       = getGenPhoton(gPart)
        event.photon_genPt  = genPhoton['pt']  if genPhoton is not None else float('nan')
        event.photon_genEta = genPhoton['eta'] if genPhoton is not None else float('nan')

      event.photonJetdR = min(deltaR(photons[0], j) for j in jets) if len(jets) > 0 else 999
      event.photonLepdR = min(deltaR(photons[0], l) for l in leptons_pt10) if len(leptons_pt10) > 0 else 999

    if addSystematicVariations:
        for var in ['jesTotalUp', 'jesTotalDown', 'jerUp', 'jerDown', 'unclustEnUp', 'unclustEnDown']: # don't use 'jer' as of now
            setattr(event, 'met_pt_'+var,  getattr(r, 'METFixEE2017_pt_'+var)  if options.year == 2017 else getattr(r, 'MET_pt_'+var) )
            setattr(event, 'met_phi_'+var, getattr(r, 'METFixEE2017_phi_'+var) if options.year == 2017 else getattr(r, 'MET_phi_'+var) )
            if not var.startswith('unclust'):
                corrFactor = 'corr_JER' if var == 'jer' else None
                alljets_sys[var]    = allJetsNotClean
                jets_sys[var]       = filter(lambda j: jetId(j, ptCut=30, absEtaCut=jetAbsEtaCut, ptVar='pt_'+var if not var=='jer' else 'pt_nom', corrFactor=corrFactor), allJets)
                bjets_sys[var]      = filter(lambda j: isBJet(j) and abs(j['eta'])<2.4, jets_sys[var])
                nonBjets_sys[var]   = filter(lambda j: not ( isBJet(j) and abs(j['eta'])<2.4), jets_sys[var])
                
                # calculate ht
                ht = sum([j['pt_nom']*j['corr_JER'] for j in jets_sys[var]]) if var == 'jer' else sum([j['pt_'+var] for j in jets_sys[var]])

                setattr(event, "nJetGood_"+var, len(jets_sys[var]))
                setattr(event, "ht_"+var,       ht)
                setattr(event, "nBTag_"+var,    len(bjets_sys[var]))

    if isSingleLep or isTriLep or isDiLep:
        event.nGoodMuons      = len(filter( lambda l:abs(l['pdgId'])==13, leptons))
        event.nGoodElectrons  = len(filter( lambda l:abs(l['pdgId'])==11, leptons))
        event.nGoodLeptons    = len(leptons)

        if len(leptons)>=1:
            event.l1_pt         = leptons[0]['pt']
            event.l1_eta        = leptons[0]['eta']
            event.l1_phi        = leptons[0]['phi']
            event.l1_pdgId      = leptons[0]['pdgId']
            event.l1_index      = leptons[0]['index']
            event.l1_miniRelIso = leptons[0]['miniPFRelIso_all']
            event.l1_relIso03   = leptons[0]['pfRelIso03_all']
            event.l1_dxy        = leptons[0]['dxy']
            event.l1_dz         = leptons[0]['dz']
            event.l1_eleIndex   = leptons[0]['eleIndex']
            event.l1_muIndex    = leptons[0]['muIndex']
            if isMC:
               event.l1_genPartFlav = leptons[0]['genPartFlav']
               event.l1_genPartIdx  = leptons[0]['genPartIdx']
 


        # For TTZ studies: find Z boson candidate, and use third lepton to calculate mt
        (event.mlmZ_mass, zl1, zl2) = closestOSDLMassToMZ(leptons_pt10)

        #if fastSim:
        #    event.reweightLeptonFastSimSF     = reduce(mul, [leptonFastSimSF.get2DSF(pdgId=l['pdgId'], pt=l['pt'], eta=l['eta'] , nvtx = r.Pileup_nTrueInt) for l in leptons], 1)
        #    event.reweightLeptonFastSimSFUp   = reduce(mul, [leptonFastSimSF.get2DSF(pdgId=l['pdgId'], pt=l['pt'], eta=l['eta'] , nvtx = r.Pileup_nTrueInt, sigma = +1) for l in leptons], 1)
        #    event.reweightLeptonFastSimSFDown = reduce(mul, [leptonFastSimSF.get2DSF(pdgId=l['pdgId'], pt=l['pt'], eta=l['eta'] , nvtx = r.Pileup_nTrueInt, sigma = -1) for l in leptons], 1)

        if isMC:
            event.reweightDilepTrigger       = 0 
            event.reweightDilepTriggerUp     = 0 
            event.reweightDilepTriggerDown   = 0 

            leptonsForSF   = ( leptons[:2] if isDiLep else (leptons[:3] if isTriLep else leptons[:1]) )
            leptonSFValues = [ leptonSF.getSF(pdgId=l['pdgId'], pt=l['pt'], eta=((l['eta'] + l['deltaEtaSC']) if abs(l['pdgId'])==11 else l['eta'])) for l in leptonsForSF ]
            event.reweightLeptonSF     = reduce(mul, [sf[0] for sf in leptonSFValues], 1)
            event.reweightLeptonSFDown = reduce(mul, [sf[1] for sf in leptonSFValues], 1)
            event.reweightLeptonSFUp   = reduce(mul, [sf[2] for sf in leptonSFValues], 1)  
            if event.reweightLeptonSF ==0:
                logger.error( "reweightLeptonSF is zero!")
            leptonHit0SFValues = [ leptonHit0SF.getSF(pdgId=l['pdgId'], pt=l['pt'], eta=((l['eta'] + l['deltaEtaSC']) if abs(l['pdgId'])==11 else l['eta'])) for l in leptonsForSF ]
            event.reweightLeptonHit0SF     = reduce(mul, [sf[0] for sf in leptonHit0SFValues], 1)
            event.reweightLeptonHit0SFDown = reduce(mul, [sf[1] for sf in leptonHit0SFValues], 1)
            event.reweightLeptonHit0SFUp   = reduce(mul, [sf[2] for sf in leptonHit0SFValues], 1)  
            if event.reweightLeptonHit0SF ==0:
                logger.error( "reweightLeptonHit0SF is zero!")
            leptonSip3dSFValues = [ leptonSip3dSF.getSF(pdgId=l['pdgId'], pt=l['pt'], eta=((l['eta'] + l['deltaEtaSC']) if abs(l['pdgId'])==11 else l['eta'])) for l in leptonsForSF ]
            event.reweightLeptonSip3dSF     = reduce(mul, [sf[0] for sf in leptonSip3dSFValues], 1)
            event.reweightLeptonSip3dSFDown = reduce(mul, [sf[1] for sf in leptonSip3dSFValues], 1)
            event.reweightLeptonSip3dSFUp   = reduce(mul, [sf[2] for sf in leptonSip3dSFValues], 1)  
            if event.reweightLeptonSip3dSF ==0:
                logger.error( "reweightLeptonSip3dSF is zero!")

            if options.fastSim:
                leptonFastSimSFValues = [ leptonFastSimSF.getSF(pdgId=l['pdgId'], pt=l['pt'], eta=((l['eta'] + l['deltaEtaSC']) if abs(l['pdgId'])==11 else l['eta'])) for l in leptonsForSF ]
                event.reweightLeptonFastSimSF     = reduce(mul, [sf[0] for sf in leptonFastSimSFValues], 1)
                event.reweightLeptonFastSimSFDown = reduce(mul, [sf[1] for sf in leptonFastSimSFValues], 1)
                event.reweightLeptonFastSimSFUp   = reduce(mul, [sf[2] for sf in leptonFastSimSFValues], 1)  
                if event.reweightLeptonFastSimSF ==0:
                    logger.error( "reweightLeptonFastSimSF is zero!")

            event.reweightLeptonTrackingSF   = reduce(mul, [leptonTrackingSF.getSF(pdgId = l['pdgId'], pt = l['pt'], eta = ((l['eta'] + l['deltaEtaSC']) if abs(l['pdgId'])==11 else l['eta']))  for l in leptonsForSF], 1)

    if isTriLep or isDiLep:
        if len(leptons)>=2:
            mt2Calculator.reset()
            event.l2_pt         = leptons[1]['pt']
            event.l2_eta        = leptons[1]['eta']
            event.l2_phi        = leptons[1]['phi']
            event.l2_pdgId      = leptons[1]['pdgId']
            event.l2_index      = leptons[1]['index']
            event.l2_miniRelIso = leptons[1]['miniPFRelIso_all']
            event.l2_relIso03   = leptons[1]['pfRelIso03_all']
            event.l2_dxy        = leptons[1]['dxy']
            event.l2_dz         = leptons[1]['dz']
            event.l2_eleIndex   = leptons[1]['eleIndex']
            event.l2_muIndex    = leptons[1]['muIndex']
            if isMC:
               event.l2_genPartFlav = leptons[1]['genPartFlav']
               event.l2_genPartIdx  = leptons[1]['genPartIdx']
            

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
            mt2Calculator.setLeptons(event.l1_pt, event.l1_eta, event.l1_phi, event.l2_pt, event.l2_eta, event.l2_phi)

            # To check MC truth when looking at the TTZToLLNuNu sample
            if isMC:
              trig_eff, trig_eff_err =  triggerEff.getSF(event.l1_pt, event.l1_eta, event.l1_pdgId, event.l2_pt, event.l2_eta, event.l2_pdgId)
              event.reweightDilepTrigger       = trig_eff 
              event.reweightDilepTriggerUp     = trig_eff + trig_eff_err
              event.reweightDilepTriggerDown   = trig_eff - trig_eff_err

              zBoson          = getGenZ(gPart)
              event.zBoson_genPt  = zBoson['pt']  if zBoson is not None else float('nan')
              event.zBoson_genEta = zBoson['eta'] if zBoson is not None else float('nan')

            if event.nPhotonGood > 0:
              gamma = ROOT.TLorentzVector()
              gamma.SetPtEtaPhiM(photons[0]['pt'], photons[0]['eta'], photons[0]['phi'], photons[0]['mass'] )
              dlg = dl + gamma
              event.dlg_mass = dlg.M()

            if options.susySignal and options.fastSim:
                mt2Calculator.setMet(getattr(r, 'GenMET_pt'), getattr(r, 'GenMET_phi'))
                setattr(event, "dl_mt2ll_gen", mt2Calculator.mt2ll())
                if len(jets)>=2:
                    bj0, bj1 = (bJets+nonBJets)[:2]
                    mt2Calculator.setBJets(bj0['pt'], bj0['eta'], bj0['phi'], bj1['pt'], bj1['eta'], bj1['phi'])
                    setattr(event, "dl_mt2bb_gen",   mt2Calculator.mt2bb())
                    setattr(event, "dl_mt2blbl_gen", mt2Calculator.mt2blbl())
                
            mt2Calculator.setMet(getattr(event, 'met_pt'), getattr(event, 'met_phi'))
            setattr(event, "dl_mt2ll", mt2Calculator.mt2ll())

            bj0, bj1 = None, None
            if len(jets)>=2:
                bj0, bj1 = (bJets+nonBJets)[:2]
                mt2Calculator.setBJets(bj0['pt'], bj0['eta'], bj0['phi'], bj1['pt'], bj1['eta'], bj1['phi'])
                setattr(event, "dl_mt2bb",   mt2Calculator.mt2bb())
                setattr(event, "dl_mt2blbl", mt2Calculator.mt2blbl())

            if addSystematicVariations:
                for var in ['jesTotalUp', 'jesTotalDown', 'jerUp', 'jerDown', 'unclustEnUp', 'unclustEnDown']: # don't use 'jer' as of now
                    mt2Calculator.setMet( getattr(event, "met_pt_"+var), getattr(event, "met_phi_"+var) )
                    setattr(event, "dl_mt2ll_"+var,  mt2Calculator.mt2ll())
                    if not var.startswith( 'unclust'):
                        if len(jets_sys[var])>=2:
                            bj0_, bj1_ = (bjets_sys[var]+nonBjets_sys[var])[:2]
                        else: 
                            bj0_, bj1_ = None, None
                    else:
                        bj0_, bj1_ = bj0, bj1
                    if bj0_ and bj1_:
                        mt2Calculator.setBJets(bj0_['pt'], bj0_['eta'], bj0_['phi'], bj1_['pt'], bj1_['eta'], bj1_['phi'])
                        setattr(event, 'dl_mt2bb_'+var,   mt2Calculator.mt2bb())
                        setattr(event, 'dl_mt2blbl_'+var, mt2Calculator.mt2blbl())

            # for quadlep stuff
            allZCands = getSortedZCandidates(leptons_pt10)
            Z_vectors = []
            for i in [0,1]:
                if len(allZCands) > i:
                    (Z_mass, Z_l1_tightLepton_index, Z_l2_tightLepton_index) = allZCands[i]
                    Z_l1_index = Z_l1_tightLepton_index if Z_l1_tightLepton_index>=0 else -1
                    Z_l2_index = Z_l2_tightLepton_index if Z_l2_tightLepton_index>=0 else -1
                    setattr(event, "Z%s_mass"%(i+1),       Z_mass)
                    setattr(event, "Z%s_l1_index"%(i+1),   Z_l1_index)
                    setattr(event, "Z%s_l2_index"%(i+1),   Z_l2_index)
                    Z_l1 = ROOT.TLorentzVector()
                    Z_l1.SetPtEtaPhiM(leptons_pt10[Z_l1_index]['pt'], leptons_pt10[Z_l1_index]['eta'], leptons_pt10[Z_l1_index]['phi'], 0 )
                    Z_l2 = ROOT.TLorentzVector()
                    Z_l2.SetPtEtaPhiM(leptons_pt10[Z_l2_index]['pt'], leptons_pt10[Z_l2_index]['eta'], leptons_pt10[Z_l2_index]['phi'], 0 )
                    Z = Z_l1 + Z_l2
                    setattr(event, "Z%s_pt"%(i+1),         Z.Pt())
                    setattr(event, "Z%s_eta"%(i+1),        Z.Eta())
                    setattr(event, "Z%s_phi"%(i+1),        Z.Phi())
                    setattr(event, "Z%s_lldPhi"%(i+1),     deltaPhi(Z_l1.Phi(), Z_l2.Phi()))
                    Z_vectors.append(Z)

            # take the leptons that are not from the leading Z candidate and assign them as nonZ, ignorant about if they actually form a Z candidate
            # As a start, take the leading two leptons as non-Z. To be overwritten as soon as we have a Z candidate, otherwise one lepton can be both from Z and non-Z
            if len(leptons_pt10)>0:
                event.nonZ1_l1_index = leptons_pt10[0]['index']
            if len(leptons)>1:
                event.nonZ1_l2_index = leptons_pt10[1]['index']
            if len(allZCands)>0:
                # reset nonZ1_leptons
                event.nonZ1_l1_index = -1
                event.nonZ1_l2_index = -1
                nonZ_tightLepton_indices = [ i for i in range(len(leptons_pt10)) if i not in [allZCands[0][1], allZCands[0][2]] ]

                event.nonZ1_l1_index = nonZ_tightLepton_indices[0] if len(nonZ_tightLepton_indices)>0 else -1
                event.nonZ1_l2_index = nonZ_tightLepton_indices[1] if len(nonZ_tightLepton_indices)>1 else -1

            event.min_dl_mass = getMinDLMass(leptons_pt10)[0][0]


    #if addSystematicVariations:
    # B tagging weights method 1a
    if isMC:
        for j in jets:
            btagEff.addBTagEffToJet(j)
        for var in btagEff.btagWeightNames:
            if var!='MC':
                setattr(event, 'reweightBTag_'+var, btagEff.getBTagSF_1a( var, bJets, filter( lambda j: abs(j['eta'])<2.4, nonBJets ) ) )
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
        gLep = gLep[:100]
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
eventRanges = reader.getEventRanges( maxNEvents = options.eventsPerJob, minJobs = options.minNJobs )

logger.info( "Splitting into %i ranges of %i events on average. FileBasedSplitting: %s. Job number %s",  
        len(eventRanges), 
        (eventRanges[-1][1] - eventRanges[0][0])/len(eventRanges), 
        'Yes',
        options.job)

#Define all jobs
jobs = [(i, eventRanges[i]) for i in range(len(eventRanges))]

filename, ext = os.path.splitext( os.path.join(output_directory, sample.name + '.root') )

if len(eventRanges)>1:
    raise RuntimeError("Using fileBasedSplitting but have more than one event range!")

clonedEvents = 0
convertedEvents = 0
outputLumiList = {}
for ievtRange, eventRange in enumerate( eventRanges ):

    logger.info( "Processing range %i/%i from %i to %i which are %i events.",  ievtRange, len(eventRanges), eventRange[0], eventRange[1], eventRange[1]-eventRange[0] )

    _logger.   add_fileHandler( outfilename.replace('.root', '.log'), options.logLevel )
    _logger_rt.add_fileHandler( outfilename.replace('.root', '_rt.log'), options.logLevel )
    
    tmp_directory = ROOT.gDirectory
    outputfile = ROOT.TFile.Open(outfilename, 'recreate')
    tmp_directory.cd()

    if options.small: 
        logger.info("Running 'small'. Not more than 10000 events") 
        nMaxEvents = eventRange[1]-eventRange[0]
        eventRange = ( eventRange[0], eventRange[0] +  min( [nMaxEvents, 10000] ) )

    # Set the reader to the event range
    reader.setEventRange( eventRange )

    logger.info("Cloning tree.")
    clonedTree = reader.cloneTree( branchKeepStrings, newTreename = "Events", rootfile = outputfile )

    clonedEvents += clonedTree.GetEntries()
    # Clone the empty maker in order to avoid recompilation at every loop iteration
    maker = treeMaker_parent.cloneWithoutCompile( externalTree = clonedTree )

    maker.start()
    # Do the thing
    reader.start()

    while reader.run():
        maker.run()
        if sample.isData:
            if maker.event.jsonPassed_:
                if reader.event.run not in outputLumiList.keys():
                    outputLumiList[reader.event.run] = set([reader.event.luminosityBlock])
                else:
                    if reader.event.luminosityBlock not in outputLumiList[reader.event.run]:
                        outputLumiList[reader.event.run].add(reader.event.luminosityBlock)

    convertedEvents += maker.tree.GetEntries()
    logger.info("Writing tree")
    maker.tree.Write()
    outputfile.Close()
    logger.info( "Written %s", outfilename)

  # Destroy the TTree
    maker.clear()
    sample.clear()


logger.info( "Converted %i events of %i, cloned %i",  convertedEvents, reader.nEvents , clonedEvents )

# Storing JSON file of processed events
if sample.isData and convertedEvents>0: # avoid json to be overwritten in cases where a root file was found already
    jsonFile = filename+'_%s.json'%(0 if options.nJobs==1 else options.job)
    LumiList( runsAndLumis = outputLumiList ).writeJSON(jsonFile)
    logger.info( "Written JSON file %s", jsonFile )

if not options.keepNanoAOD and not options.skipNanoTools:
    for f in sample.files:
        try:
            os.remove(f)
            logger.info("Removed nanoAOD file: %s", f)
        except OSError:
            logger.info("nanoAOD file %s seems to be not there", f)

logger.info("Copying log file to %s", output_directory )
copyLog = subprocess.call(['cp', logFile, output_directory] )
if copyLog:
    logger.info( "Copying log from %s to %s failed", logFile, output_directory)
else:
    logger.info( "Successfully copied log file" )
    os.remove(logFile)
    logger.info( "Removed temporary log file" )



# Copying output to DPM or AFS and check the files
if options.writeToDPM:

    for dirname, subdirs, files in os.walk( output_directory ):
        logger.debug( 'Found directory: %s',  dirname )

        for fname in files:

            if fname.startswith("nanoAOD_") or "_for_" in fname: continue # do not copy the nanoAODTools files
            if not fname.endswith(".root"): continue # remove that for copying log files

            source  = os.path.abspath( os.path.join( dirname, fname ) )
            target  = os.path.join( targetPath, fname )

            if fname.endswith(".root"): # redundant, but keep it in case you want to copy log files
                if checkRootFile( source, checkForObjects=["Events"] ) and deepCheckRootFile( source ) and deepCheckWeight( source ):
                    logger.info( "Source: File check ok!" )
                else:
                    raise Exception("Corrupt rootfile at source! File not copied: %s"%source )

            cmd = [ 'xrdcp', '-f',  source, target ]
            logger.info( "Issue copy command: %s", " ".join( cmd ) )
            subprocess.call( cmd )

            if fname.endswith(".root"):
                if checkRootFile( target, checkForObjects=["Events"] ) and deepCheckRootFile( target ) and deepCheckWeight( target ):
                    logger.info( "Target: File check ok!" )
                else:
                    logger.info( "Corrupt rootfile at target! Trying again: %s"%target )
                    logger.info( "2nd try: Issue copy command: %s", " ".join( cmd ) )
                    subprocess.call( cmd )

                    # Many files are corrupt after copying, a 2nd try fixes that
                    if checkRootFile( target, checkForObjects=["Events"] ) and deepCheckRootFile( target ) and deepCheckWeight( target ):
                        logger.info( "2nd try successfull!" )
                    else:
                        # if not successful, the corrupt root file needs to be deleted from DPM
                        logger.info( "2nd try: No success, removing file: %s"%target )
                        logger.info( "Issue rm command: %s", " ".join( cmd ) )
                        cmd = [ "xrdfs", redirector_hephy, "rm", "/cms" + target.split("/cms")[1] ]
                        subprocess.call( cmd )
                        raise Exception("Corrupt rootfile at target! File not copied: %s"%source )

    # Clean up.
    if not options.runOnLxPlus:
        # not needed on condor, container will be removed automatically
        subprocess.call( [ 'rm', '-rf', output_directory ] ) # Let's risk it.

else:
    if checkRootFile( outfilename, checkForObjects=["Events"] ) and deepCheckRootFile( outfilename ) and deepCheckWeight( outfilename ):
        logger.info( "Target: File check ok!" )
    else:
        logger.info( "Corrupt rootfile! Removing file: %s"%outfilename )
        #os.remove( outfilename )
        #raise Exception("Corrupt rootfile! File not copied: %s"%source )

# There is a double free corruption due to stupid ROOT memory management which leads to a non-zero exit code
# Thus the job is resubmitted on condor even if the output is ok
# Current idea is that the problem is with xrootd having a non-closed root file
sample.clear()


