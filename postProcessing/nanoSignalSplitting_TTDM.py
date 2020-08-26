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
from StopsDilepton.tools.helpers import closestOSDLMassToMZ, writeObjToFile, m3, deltaR, bestDRMatchInCollection
from StopsDilepton.samples.helpers           import getTTDMSignalWeightForEvent, getTTDMSignalWeight, getTTDMBranchNames

from Analysis.Tools.helpers import deepCheckRootFile


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
        nargs='?',
        type=int,
        default=0,
        help="Run only job i"
        )

    argParser.add_argument('--targetDir',
        action='store',
        nargs='?',
        type=str,
        default=user.data_output_directory,
        help="Name of the directory the post-processed files will be saved"
        )

    argParser.add_argument('--inputDir',
        action='store',
        nargs='?',
        type=str,
        default=user.data_output_directory,
        help="Name of the directory the post-processed files are read"
        )


    argParser.add_argument('--processingEra',
        action='store',
        nargs='?',
        type=str,
        default='postProcessed_80X_v21',
        help="Name of the processing era"
        )

    argParser.add_argument('--skim',
        action='store',
        nargs='?',
        type=str,
        default='dilepTiny',
        help="Skim conditions to be applied for post-processing"
        )

    argParser.add_argument('--T2bW',
        action='store_true',
        help="Is T2tt signal?"
        )
    argParser.add_argument('--T2bt',
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
    
    argParser.add_argument('--T8bbllnunu',
        action='store_true',
        help="Is T8bbllnunu signal?"
        )

    argParser.add_argument('--T8bbstausnu',
        action='store_true',
        help="Is T8bbstausnu signal?"
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

    argParser.add_argument('--year',
        action='store',
        type=int,
        help="Which year?"
        )


    return argParser

options = get_parser().parse_args()

# Logging
import StopsDilepton.tools.logger as logger
logger = logger.get_logger(options.logLevel, logFile ='%s_%s.txt'%(options.skim, '_'.join(options.samples) ) )
logFileLocation = '%s_%s.txt'%(options.skim, '_'.join(options.samples) )

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

print [ s.name for s in samples ]
directory  = os.path.join(options.targetDir, options.processingEra)
output_directory = os.path.join( directory, options.skim, samples[0].name )
print output_directory

#Samples: Load samples
maxN = 2 if options.small else None

if len(samples)==0:
    logger.info( "No samples found. Was looking for %s. Exiting" % options.samples )
    sys.exit(-1)

isData = False not in [s.isData for s in samples]
isMC   =  True not in [s.isData for s in samples]

sample_name_postFix = ""

inDir = os.path.join(options.inputDir, options.processingEra, options.skim, samples[0].name)
outDir = os.path.join(options.targetDir, options.processingEra, options.skim, samples[0].name)


# Directory for individual signal files
signalDir = os.path.join(options.targetDir, options.processingEra, options.skim, "TTDM")
if not os.path.exists(signalDir): os.makedirs(signalDir)

try:    #Avoid trouble with race conditions in multithreading
    os.makedirs(outDir)
    logger.info( "Created output directory %s.", outDir )
except:
    pass

def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i + n]

signalWeight = getTTDMSignalWeight( samples[0], lumi = 1000) # lumi doesn't matter here

nJobs = options.nJobs
chunkSize = len(signalWeight.keys())/nJobs
if not len(signalWeight.keys())%nJobs == 0: chunkSize += 1


if options.samples[0].count('_pseudoscalar'):
    branchNames = getTTDMBranchNames('pseudoscalar')
if options.samples[0].count('_scalar'):
    branchNames = getTTDMBranchNames('scalar')

masspoints = list(chunks(branchNames, chunkSize))

job = options.job

print "All masspoints:"
print masspoints

print "Running over:"
print masspoints[job]

# Write one file per mass point for T2tt
output = Sample.fromDirectory("TTDM_output", inDir)

print "Initialising chain, otherwise first mass point is empty"
print output.chain
if options.small: output.reduceFiles( to = 1 )
for i,s in enumerate(masspoints[job]):
    sampleName = "TTDM_%s_%s_%s"%(signalWeight[s]['spin'], signalWeight[s]['mPhi'], signalWeight[s]['mChi'])
    #cut = "GenSusyMStop=="+str(s[0])+"&&GenSusyMNeutralino=="+str(s[1]) #FIXME
    logger.info("Going to write: %s", sampleName)
    print s, type(s)
    cut = "(1)&&(%s>0)"%s
    logger.debug("Using cut %s", cut)
    signalFile = os.path.join(signalDir, sampleName + '.root' )
    #signalFile = os.path.join(signalDir, 'T2tt_'+str(s[0])+'_'+str(s[1])+'.root' )
    logger.debug("Ouput file will be %s", signalFile)
    if os.path.exists(signalFile) and deepCheckRootFile(signalFile):
        c = ROOT.TChain("Events")
        c.Add(signalFile)
        if c.GetEntries()==0:
            options.overwrite = True # :-)

    print signalFile
    check = deepCheckRootFile(signalFile)
    if not (os.path.exists(signalFile) and deepCheckRootFile(signalFile)) or options.overwrite:
        outF = ROOT.TFile.Open(signalFile, "RECREATE")
        print cut
        t = output.chain.CopyTree(cut)
        nEvents = t.GetEntries()
        outF.Write()
        outF.Close()
        logger.info( "Number of events %i", nEvents)
        inF = ROOT.TFile.Open(signalFile, "READ")
        try:
            u = inF.Get("Events")
        except ReferenceError:
            logger.info( "Found null pointerfor mStop %i mNeu %i. Continue. ", s[0],s[1]) 
            u = None

        if u:
            nnEvents = u.GetEntries()
            logger.debug("Number of events in tree %i and in file %i", nEvents, nnEvents)
            if nEvents == nnEvents: logger.debug("All events written")
            else: logger.debug("Something went wrong, discrepancy between file and tree")
        inF.Close()
        logger.info( "Written signal file for %s to %s", s, signalFile)
    else:
        logger.info( "Found file %s -> Skipping"%(signalFile) )
    logger.info("Done with %s/%s", i+1, len(masspoints[job])) 

output.clear()

