#!/usr/bin/env python
''' Analysis script for standard plots
'''
#
# Standard imports and batch mode
#
import ROOT, os
ROOT.gROOT.SetBatch(True)
import itertools
import copy

# for smearing
import numpy as np

from math                         import sqrt, cos, sin, pi, atan2
from RootTools.core.standard      import *
from StopsDilepton.tools.user            import plot_directory
from StopsDilepton.tools.helpers         import deltaR, deltaPhi, getObjDict, getVarValue
from Samples.Tools.metFilters import getFilterCut
from StopsDilepton.tools.cutInterpreter  import cutInterpreter


import argparse
argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--logLevel',           action='store',      default='INFO',          nargs='?', choices=['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'TRACE', 'NOTSET'], help="Log level for logging")
argParser.add_argument('--small',                                   action='store_true',     help='Run only on a small subset of the data?', )
argParser.add_argument('--noPUReweighting',    action='store_true',     help='Dont do PU reweighting', )
argParser.add_argument('--plot_directory',     action='store',      default='dataToData_HEM_newApproach_v0p2')
argParser.add_argument('--selection',          action='store',      default='njet2p-btag0-relIso0.12-looseLeptonVeto-mll20-dPhiJet0-dPhiJet1-onZ') #njet2p-btag0-relIso0.12-looseLeptonVeto-mll20-met80-metSig5-dPhiJet0-dPhiJet1 # lepSel-OS-looseLeptonVeto-njet2p-btag0-relIso0.12-mll20
args = argParser.parse_args()

#
# Logger
#
import StopsDilepton.tools.logger as logger
import RootTools.core.logger as logger_rt
logger    = logger.get_logger(   args.logLevel, logFile = None)
logger_rt = logger_rt.get_logger(args.logLevel, logFile = None)

if args.small:                        args.plot_directory += "_small"
if args.noPUReweighting:              args.plot_directory += "_noPUReweighting"

data_directory           = "/afs/hephy.at/data/dspitzbart03/nanoTuples/"
postProcessing_directory = "stops_2018_nano_v0p2/dilep/"
from StopsDilepton.samples.nanoTuples_Run2018_PromptReco_postProcessed import *


# define 2l selections
offZ = "&&abs(dl_mass-91.1876)>15" if not (args.selection.count("onZ") or args.selection.count("allZ") or args.selection.count("offZ")) else ""
def getLeptonSelection( mode ):
  if   mode=="mumu": return "nGoodMuons==2&&nGoodElectrons==0&&isOS&&isMuMu" + offZ
  elif mode=="mue":  return "nGoodMuons==1&&nGoodElectrons==1&&isOS&&isEMu"
  elif mode=="ee":   return "nGoodMuons==0&&nGoodElectrons==2&&isOS&&isEE" + offZ
  elif mode=='all':  return "nGoodMuons+nGoodElectrons==2"


mode = 'mumu'

data_2018_sample            = DoubleMuon_Run2018
data_2018_sample_preHEM     = copy.deepcopy(data_2018_sample)
data_2018_sample_postHEM    = copy.deepcopy(data_2018_sample)

data_2018_sample_preHEM.setSelectionString(["run < 319077"])
data_2018_sample_postHEM.setSelectionString(["run >= 319077"])

data_2018_sample_preHEM.texName    = "data 2018 (2#mu) (run < 319077)"
data_2018_sample_postHEM.texName   = "data 2018 (2#mu) (run #geq 319077)"
data_2018_samples           = [data_2018_sample, data_2018_sample_preHEM, data_2018_sample_postHEM]

colors2018 = [ROOT.kGreen+1, ROOT.kBlue+1, ROOT.kRed+1]
for i_s, data_2018_sample in enumerate(data_2018_samples):
    data_2018_sample.addSelectionString([getFilterCut(isData=True, year=2018), getLeptonSelection(mode)])
    data_2018_sample.read_variables = ["event/I","run/I"]
    data_2018_sample.style          = styles.errorStyle(colors2018[i_s], markerStyle=20+i_s, width=2)
    lumi_2018_scale                 = data_2018_sample.lumi/1000

# reweight to the proper pile-up distribution
def get_nVtx_reweight( histo ):
    def get_histo_reweight( event, sample) :
        return histo.GetBinContent(sample.nVert_histo.FindBin( event.PV_npvsGood ))/sample.nVert_histo.GetBinContent(sample.nVert_histo.FindBin( event.PV_npvsGood ) )
    return get_histo_reweight

reweight_binning = [3*i for i in range(10)]+[30,35,40,50,100]
for i_s, data_2018_sample in enumerate(data_2018_samples):

    print i_s
    logger.info('nVert Histo for %s', data_2018_sample.name)
    data_2018_sample.nVert_histo     = data_2018_sample.get1DHistoFromDraw("PV_npvsGood", reweight_binning, binningIsExplicit = True)
    data_2018_sample.nVert_histo.Scale(1./data_2018_sample.nVert_histo.Integral())

    data_2018_samples[i_s].weight = get_nVtx_reweight(data_2018_samples[0].nVert_histo)

selection_string = cutInterpreter.cutString(args.selection)

totalEvents = data_2018_sample.getYieldFromDraw(selection_string)
preHEMEvents = data_2018_sample_preHEM.getYieldFromDraw(selection_string)
postHEMEvents = data_2018_sample_postHEM.getYieldFromDraw(selection_string)

hist_preHEM = data_2018_sample_preHEM.get1DHistFromDraw("dl_mt2ll", binning=[0,25,50,75,100,140,200,240,500], selectionString=selection_string, binningIsExplicit = True, addOverflow='upper')
hist_preHEM.Scale(1/hist_preHEM.Integral())
hist_postHEM = data_2018_sample_preHEM.get1DHistFromDraw("dl_mt2ll", binning=[0,25,50,75,100,140,200,240,500], selectionString=selection_string, binningIsExplicit = True, addOverflow='upper')
hist_postHEM.Scale(1/hist_postHEM.Integral())

hist_weight = hist_postHEM.Divide(hist_preHEM)
hist_weight.Scale(postHEMEvents/totalEvents)




