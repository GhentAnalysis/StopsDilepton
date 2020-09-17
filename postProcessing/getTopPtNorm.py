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

from array import array
from operator import mul
from math import sqrt, atan2, sin, cos

# RootTools
from RootTools.core.standard import *

# User specific
import StopsDilepton.tools.user as user

# top pt reweighting
from StopsDilepton.tools.topPtReweighting import getUnscaledTopPairPtReweightungFunction, getTopPtDrawString, getTopPtsForReweighting

def get_parser():
    ''' Argument parser for post-processing module.
    '''
    import argparse
    argParser = argparse.ArgumentParser(description = "Argument parser for cmgPostProcessing")

    argParser.add_argument('--logLevel',    action='store',         nargs='?',  choices=['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'TRACE', 'NOTSET'],   default='INFO', help="Log level for logging" )
    argParser.add_argument('--samples',     action='store',         nargs='*',  type=str, default=['TTZToLLNuNu_ext'],                  help="List of samples to be post-processed, given as CMG component name" )
    argParser.add_argument('--skim',        action='store',         nargs='?',  type=str, default='dilepTiny',                          help="Skim conditions to be applied for post-processing" )
    argParser.add_argument('--small',       action='store_true',                                                                        help="Run the file on a small sample (for test purpose), bool flag set to True if used" )
    argParser.add_argument('--year',        action='store',                     type=int,                                               help="Which year?" )

    return argParser

options = get_parser().parse_args()

# Logging
import StopsDilepton.tools.logger as _logger
logger  = _logger.get_logger(options.logLevel, logFile = None)

import RootTools.core.logger as _logger_rt
logger_rt = _logger_rt.get_logger(options.logLevel, logFile = None )

# Flags 
isDiLep         = options.skim.lower().startswith('dilep')
isTriLep        = options.skim.lower().startswith('trilep')
isSingleLep     = options.skim.lower().startswith('singlelep')
isTiny          = options.skim.lower().count('tiny')
isSmall         = options.skim.lower().count('small')
isInclusive     = options.skim.lower().count('inclusive')

# Skim condition
skimConds = []

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

#Samples: Load samples
maxN = 1 if options.small else None

if options.year == 2016:
    from Samples.nanoAOD.Summer16_private_legacy_v1 import allSamples as bkgSamples
    from Samples.nanoAOD.Spring16_private           import allSamples as signalSamples
    from Samples.nanoAOD.Run2016_17Jul2018_private  import allSamples as dataSamples
    allSamples = bkgSamples + signalSamples + dataSamples
elif options.year == 2017:
    from Samples.nanoAOD.Fall17_private_legacy_v1   import allSamples as bkgSamples
    from Samples.nanoAOD.Run2017_31Mar2018_private  import allSamples as dataSamples
    allSamples = bkgSamples + dataSamples
elif options.year == 2018:
    from Samples.nanoAOD.Spring18_private           import allSamples as HEMSamples
    from Samples.nanoAOD.Run2018_26Sep2018_private  import allSamples as HEMDataSamples
    from Samples.nanoAOD.Autumn18_private_legacy_v1 import allSamples as bkgSamples
    from Samples.nanoAOD.Run2018_17Sep2018_private  import allSamples as dataSamples
    allSamples = HEMSamples + HEMDataSamples + bkgSamples + dataSamples
else:
    raise NotImplementedError

samples = []
for selectedSamples in options.samples:
    for sample in allSamples:
        if selectedSamples == sample.name:
            samples.append(sample)

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


if options.small:
    sample.reduceFiles(to=10)

logger.info("Obtaining top pt scale factor for sample: %s", sample.name)

selectionString = "&&".join(skimConds)

logger.info("Using selection: %s", selectionString)

reweighted  = sample.getYieldFromDraw( selectionString = selectionString, weightString = getTopPtDrawString(selection = "dilep") + '*genWeight')
central     = sample.getYieldFromDraw( selectionString = selectionString, weightString = 'genWeight')

topScaleF = central['val']/reweighted['val']

logger.info("Got the following SF: %s", topScaleF)
