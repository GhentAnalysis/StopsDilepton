'''
Produces cutflow tables for reinterpretation
Quick hack from Priya's version - closes with card files up to sub-percent levels.
'''

import ROOT
import os
import pickle
from  math import sqrt

# RootTools
from RootTools.core.standard import *
from Analysis.Tools.metFilters            import getFilterCut

# argParser
import argparse
argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--logLevel', 
      action='store',
      nargs='?',
      choices=['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'TRACE', 'NOTSET'],
      default='INFO',
      help="Log level for logging"
)
args = argParser.parse_args()

# Logging
import StopsDilepton.tools.logger as logger
from StopsDilepton.tools.user import plot_directory
logger = logger.get_logger(args.logLevel, logFile = None )
import RootTools.core.logger as logger_rt
logger_rt = logger_rt.get_logger(args.logLevel, logFile = None )

from StopsDilepton.analysis.Setup import Setup

#Samples
allsamples  = {}
allweights  = {}
allcuts     = {}
filters     = {}

setup16 = Setup(year=2016)
allweights[2016]    = setup16.preselection('MC', channel='all', isFastSim=True)['weightStr'] + '*%s'%(setup16.lumi/1000.)
filters[2016]       = getFilterCut(isData=False, year=2016, isFastSim=True)

setup17 = Setup(year=2017, puWeight='VUp')
allweights[2017]    = setup17.preselection('MC', channel='all', isFastSim=True)['weightStr'] + '*%s'%(setup17.lumi/1000.)
filters[2017]       = getFilterCut(isData=False, year=2017, isFastSim=True) + '&&Sum$((2.6<abs(Jet_eta)&&abs(Jet_eta)<3&&Jet_pt>30))==0' ## EE veto

setup18 = Setup(year=2018)
allweights[2018]    = setup18.preselection('MC', channel='all', isFastSim=True)['weightStr'] + '*%s'%(setup18.lumi/1000.)
filters[2018]       = getFilterCut(isData=False, year=2018, isFastSim=True)

for year in [2016,2017,2018]:
    allweights[year] += '*reweightLeptonFastSimSF*reweight_nISR'

data_directory              = '/afs/hephy.at/data/cms06/nanoTuples/'
postProcessing_directory    = 'stops_2016_nano_v0p23/dilep/'
from StopsDilepton.samples.nanoTuples_FastSim_Summer16_postProcessed import *
allsamples[2016] = [T2tt_350_150, T2tt_800_100, T2tt_600_300]

data_directory              = '/afs/hephy.at/data/cms06/nanoTuples/'
postProcessing_directory    = 'stops_2017_nano_v0p23/dilep/'
from StopsDilepton.samples.nanoTuples_FastSim_Fall17_postProcessed import *
allsamples[2017] = [T2tt_350_150, T2tt_800_100, T2tt_600_300]

data_directory              = '/afs/hephy.at/data/cms06/nanoTuples/'
postProcessing_directory    = 'stops_2018_nano_v0p23/dilep/'
from StopsDilepton.samples.nanoTuples_FastSim_Autumn18_postProcessed import *
allsamples[2018]    = [T2tt_350_150, T2tt_800_100, T2tt_600_300]


cuts=[
  ("==2 leptons,l1pt>30,l2pt>20",           "$N_{\\textrm{leptons}}=2$, opposite charge",                  "(nGoodMuons+nGoodElectrons)==2&&(isEMu==1||isMuMu==1||isEE==1)&&l1_pt>30&&l2_pt>20&&isOS==1"),
  ("miniIso <=0.2",                         "miniIso <= 0.2",                           "l1_miniRelIso < 0.2 && l2_miniRelIso < 0.2"),
  ("looseLeptonMiniIsoVeto",                "loose lepton veto",               "(Sum$(Electron_pt>15&&abs(Electron_eta)<2.4&&Electron_miniPFRelIso_all<0.4) + Sum$(Muon_pt>15&&abs(Muon_eta)<2.4&&Muon_miniPFRelIso_all<0.4) )==2"),
  ("m(ll)>20",                              "$m(\ell\ell)>20$ GeV",                           "dl_mass>=20"),
  ("|m(ll) - mZ|>15 for SF",                "$m_{Z}-m(\ell\ell) > 15$ GeV (SF)",            "( (isMuMu==1||isEE==1)&&abs(dl_mass-91.1876)>=15 || isEMu==1 )"),
  (">=2 jets",                              "$N_{\\textrm{jets}}\geq2$",                             "nJetGood>=2"),
  (">=1 b-tags",                            "$N_{\\textrm{b}}\geq1$",                           "nBTag>=1"),
  ("MET_significance >= 12",                "$\pazocal{S}>12$",              "MET_significance>=12"), 
  ("dPhiJetMET",                            "$\Delta \phi(p_{\\textrm{T}}^{\\textrm{miss}}, \\textrm{jets})$ veto",             "cos(met_phi-JetGood_phi[0])<0.8&&cos(met_phi-JetGood_phi[1])<cos(0.25)"),
  ("MT2(ll) >= 100",                        "$M_{\\textrm{T2}}(\ell\ell)>100$ GeV",                    "dl_mt2ll>=100"),
  ("MT2(ll) >= 140",                        "$M_{\\textrm{T2}}(\ell\ell)>140$ GeV",                    "dl_mt2ll>=140"),
  ("MT2(ll) >= 240",                        "$M_{\\textrm{T2}}(\ell\ell)>240$ GeV",                    "dl_mt2ll>=240"),
    ]

allcuts[2016] = [('filters', 'filters', filters[2016])] + cuts
allcuts[2017] = [('filters', 'filters', filters[2017])] + cuts
allcuts[2018] = [('filters', 'filters', filters[2018])] + cuts



import copy
template = { name:0 for name, texname, sel in allcuts[2016] }
table_row = { texname:0 for name, texname, sel in allcuts[2016] }
results = { sample.name:copy.deepcopy(template) for sample in allsamples[2016] }
table = { sample.name:copy.deepcopy(table_row) for sample in allsamples[2016] }

for i_cut, cut in enumerate(allcuts[2016]):
    name, texname, sel = cut
    for year in [2016,2017,2018]:
        samples = allsamples[year]
        cut = allcuts[year]
        weight_string = allweights[year]
        for i_sample, sample in enumerate(samples):
            res = sample.getYieldFromDraw(selectionString = "&&".join([ '('+c[2]+')' for c in cut[:i_cut+1]]), weightString = weight_string)['val']
            results[sample.name][name] += res
            table[sample.name][texname] += res

from StopsDilepton.tools.xSecSusy import xSecSusy
xSecSusy_ = xSecSusy()
for sample in allsamples[2016]:
    mStop = int(sample.mStop)
    table[sample.name]['$\sigma \\times \pazocal{B}$'] = xSecSusy_.getXSec(channel='stop13TeV',mass=mStop,sigma=0) * (3*0.108)**2 * (setup16.lumi+setup17.lumi+setup18.lumi)

import pandas as pd
df = pd.DataFrame(results)
tex_table = pd.DataFrame(table)
print
print "### Cut-flow ###"
print df.sort_values('T2tt_800_100', ascending=False)


print
print "### Latex table ###"
print tex_table.round(2).sort_values('T2tt_800_100', ascending=False).to_latex(escape=False)


