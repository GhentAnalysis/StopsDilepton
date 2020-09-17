#!/usr/bin/env python
''' Analysis script for standard plots
'''
#
# Standard imports and batch mode
#
import ROOT, os
ROOT.gROOT.SetBatch(True)
import itertools

from math                                import sqrt, cos, sin, pi, atan2
from RootTools.core.standard             import *
from StopsDilepton.tools.user            import plot_directory, analysis_results
from Samples.Tools.metFilters            import getFilterCut
from StopsDilepton.tools.cutInterpreter  import cutInterpreter

from StopsDilepton.tools.resultsDB      import resultsDB

#
# Arguments
# 
import argparse
argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--logLevel',           action='store',      default='INFO',          nargs='?', choices=['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'TRACE', 'NOTSET'], help="Log level for logging")
argParser.add_argument('--small',                                   action='store_true',     help='Run only on a small subset of the data?', )
argParser.add_argument('--plot_directory',     action='store',      default='v0p3')
argParser.add_argument('--year',               action='store', type=int,      default=2016)
args = argParser.parse_args()

#
# Logger
#
import StopsDilepton.tools.logger as logger
import RootTools.core.logger as logger_rt
logger    = logger.get_logger(   args.logLevel, logFile = None)
logger_rt = logger_rt.get_logger(args.logLevel, logFile = None)

if args.year == 2016:
    data_directory = "/afs/hephy.at/data/dspitzbart01/nanoTuples/"
    postProcessing_directory = "stops_2016_nano_v0p7/dilep/"
    from StopsDilepton.samples.nanoTuples_Summer16_postProcessed import *
    data_directory = "/afs/hephy.at/data/dspitzbart03/nanoTuples/"
    postProcessing_directory = "stops_2016_nano_v0p5/dilep/"
    from StopsDilepton.samples.nanoTuples_Run2016_17Jul2018_postProcessed import *
    mc          = [ Top_pow_16, TTXNoZ_16, TTZ_16, multiBoson_16, DY_HT_LO_16, ZZ4l_16]
    data_sample = Run2016

elif args.year == 2017:
    data_directory = "/afs/hephy.at/data/dspitzbart01/nanoTuples/"
    postProcessing_directory = "stops_2017_nano_v0p7/dilep/"
    from StopsDilepton.samples.nanoTuples_Fall17_postProcessed import *
    data_directory = "/afs/hephy.at/data/dspitzbart03/nanoTuples/"
    postProcessing_directory = "stops_2017_nano_v0p6/dilep/"
    from StopsDilepton.samples.nanoTuples_Run2017_31Mar2018_postProcessed import *
    mc          = [ Top_pow_17, TTXNoZ_17, TTZ_17, multiBoson_17, DY_HT_LO_17, ZZ4l_17]
    data_sample = Run2017

elif args.year == 2018:
    data_directory = "/afs/hephy.at/data/dspitzbart01/nanoTuples/"
    postProcessing_directory = "stops_2018_nano_v0p7/dilep/"
    from StopsDilepton.samples.nanoTuples_Autumn18_postProcessed import *
    data_directory = "/afs/hephy.at/data/dspitzbart03/nanoTuples/"
    postProcessing_directory = "stops_2018_nano_v0p5/dilep/"
    from StopsDilepton.samples.nanoTuples_Run2018_PromptReco_postProcessed import *
    mc          = [ Top_pow_18, TTXNoZ_18, TTZ_18, multiBoson_18, DY_HT_LO_18, ZZ4l_18]
    data_sample = Run2018


## cache results
res = resultsDB(analysis_results+'ttZ_results.sq', "ttZ", ["year", "sample", "channel", "region"])

## import regions
from StopsDilepton.analysis.regions         import *
from StopsDilepton.tools.trilepSelection    import *

channels = ["3mu", "2mu1e", "2e1mu", "3e"]
#regions  = getRegions2D("nJetGood", [2,3], "nBTag", [2,-1]) + getRegions2D("nJetGood", [3,4,-1], "nBTag", [1,2,-1])
lumi     = data_sample.lumi/1000.
logger.info("Using lumi: %s", lumi)
presel   = "abs(Z1_mass-91.2)<20 && nJetGood>=2 && nBTag>=1 " #  && min_dl_mass>12" # not in data samples yet
#presel   = "Sum$(lep_pt>20&&abs(lep_eta)<2.4&&lep_pfRelIso03_all<0.12)>2&&l1_pt>40&&l2_pt>20&&nJetGood>=3&&nBTag>=1&&abs(Z1_mass-91.1876)<20"

lepSel = 'Sum$(lep_pt>20&&abs(lep_eta)<2.4&&lep_pfRelIso03_all<0.12)>2&&l1_pt>40&&l2_pt>20'

histos = {}
yields = {}

## MC yields ##
filterCut   = getFilterCut(isData=False, year=args.year)
selection   = '&&'.join([lepSel, filterCut, presel])
for sample in mc:
    # X: nJet, Y: nBTag
    tmp = sample.get2DHistoFromDraw("nBTag:nJetGood", [10, 0, 10, 10, 0, 10], selectionString=selection, weightString="weight * %s * reweightPU*reweightDilepTrigger*reweightLeptonSF*reweightBTag_SF*reweightLeptonTrackingSF"%lumi)
    histos[sample.name] = tmp
    yields[sample.name] = { '2j2b': tmp.Integral(3,3,3,3),
                            '3j1b': tmp.Integral(4,4,2,2),
                            '3j2b': tmp.Integral(4,4,3,10),
                            '4j1b': tmp.Integral(5,10,2,2),
                            '4j2b': tmp.Integral(5,10,3,10)}

## Data yields ##
filterCut   = getFilterCut(isData=False, year=args.year)
selection   = '&&'.join([lepSel, filterCut, presel])
tmp = data_sample.get2DHistoFromDraw("nBTag:nJetGood", [10, 0, 10, 10, 0, 10], selectionString=selection, weightString="(1)")
histos['data'] = tmp
yields['data'] = { '2j2b': tmp.Integral(3,3,3,3),
                   '3j1b': tmp.Integral(4,4,2,2),
                   '3j2b': tmp.Integral(4,4,3,10),
                   '4j1b': tmp.Integral(5,10,2,2),
                   '4j2b': tmp.Integral(5,10,3,10)}

import pickle
pickle.dump(yields, file('ttZ_yields_%s.pkl'%args.year, 'w'))


#for channel in channels:
#    logger.info("Channel: %s", channel)
#    for region in regions:
#        logger.info("Region: %s", region.cutString())
#        for sample in mc:
#            logger.info("Sample: %s", sample.name)
#            # should switch to setup/estimated based workflow once this is setup
#            filterCut   = getFilterCut(isData=False, year=args.year)
#            selection   = '&&'.join([presel, filterCut, region.cutString(), getTrilepSelection(channel)])
#            
#            if not res.contains({'year':args.year, 'sample':sample.name, 'channel':channel, 'region':region.cutString()}):
#                pred = sample.getYieldFromDraw( selection, "weight * %s * reweightPU36fb*reweightDilepTrigger*reweightLeptonSF*reweightBTag_SF*reweightLeptonTrackingSF"%lumi)['val']
#                res.add({'year':args.year, 'sample':sample.name, 'channel':channel, 'region':region.cutString()}, pred, overwrite=False)
#            else:
#                pred = res.get({'year':args.year, 'sample':sample.name, 'channel':channel, 'region':region.cutString()})
#
#            logger.info("Results: %s", pred)
#
#        logger.info("Data")
#        sample = data_sample
#        filterCut   = getFilterCut(isData=True, year=args.year)
#        selection   = '&&'.join([presel, filterCut, region.cutString(), getTrilepSelection(channel)])
#        if not res.contains({'year':args.year, 'sample':sample.name, 'channel':channel, 'region':region.cutString()}):
#            pred = sample.getYieldFromDraw( selection, "(1)")['val']
#            res.add({'year':args.year, 'sample':sample.name, 'channel':channel, 'region':region.cutString()}, pred, overwrite=False)
#        else:
#            pred = res.get({'year':args.year, 'sample':sample.name, 'channel':channel, 'region':region.cutString()})
#
#        logger.info("Results: %s", pred)
#
            
    


## combine channels (and years)



## plot



