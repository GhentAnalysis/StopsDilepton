''' Analysis script for 1D 2l plots (RootTools)
'''

#Standard imports
import ROOT
from math import sqrt, cos, sin, pi, acos
import itertools,os
import copy

import argparse
argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--logLevel',           action='store',      default='INFO',          nargs='?', choices=['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'TRACE', 'NOTSET'], help="Log level for logging")
args = argParser.parse_args()


#RootTools
from RootTools.core.standard import *
from StopsDilepton.tools.user import data_directory as user_data_directory

from StopsDilepton.samples.color import color

import StopsDilepton.tools.logger as logger
import RootTools.core.logger as logger_rt
logger    = logger.get_logger(   args.logLevel, logFile = None)
logger_rt = logger_rt.get_logger(args.logLevel, logFile = None)

data_directory              = '/afs/hephy.at/data/cms06/nanoTuples/'
postProcessing_directory    = 'stops_2018_nano_v0p23/dilep/'
from StopsDilepton.samples.nanoTuples_FastSim_Autumn18_postProcessed import *

data_directory              = '/afs/hephy.at/data/cms03/nanoTuples/'
postProcessing_directory    = 'stops_2018_nano_v0p25/dilep/'
from StopsDilepton.samples.nanoTuples_Autumn18_TTDM_postProcessed import *

presel  = '((isMuMu==1&&nGoodMuons==2&&nGoodElectrons==0) || (isEMu==1&&nGoodMuons==1&&nGoodElectrons==1) || (isEE==1&&nGoodMuons==0&&nGoodElectrons==2))'
presel += '&&(Sum$(Electron_pt>15&&abs(Electron_eta)<2.4&&Electron_pfRelIso03_all<0.4) + Sum$(Muon_pt>15&&abs(Muon_eta)<2.4&&Muon_pfRelIso03_all<0.4) )==2'
presel += '&&dl_mass>=20&&l1_pt>30&&nJetGood>=2&&nBTag>=1'

#h_800_central   = T2tt_800_100.get1DHistoFromDraw( "dl_mt2ll", selectionString=presel, binning = [0, 25,50,75,100,140,240,340], binningIsExplicit=True, addOverFlowBin = 'upper', weightString = "60 * weight * reweightDilepTrigger * reweightLeptonSF * reweightBTag_SF * reweightPU" )
#h_800_vup       = T2tt_800_100.get1DHistoFromDraw( "dl_mt2ll", selectionString=presel, binning = [0, 25,50,75,100,140,240,340], binningIsExplicit=True, addOverFlowBin = 'upper', weightString = "60 * weight * reweightDilepTrigger * reweightLeptonSF * reweightBTag_SF * reweightPUVUp" )
#
#h_350_central   = T2tt_350_150.get1DHistoFromDraw( "dl_mt2ll", selectionString=presel, binning = [0, 25,50,75,100,140,240,340], binningIsExplicit=True, addOverFlowBin = 'upper', weightString = "60 * weight * reweightDilepTrigger * reweightLeptonSF * reweightBTag_SF * reweightPU" )
#h_350_vup       = T2tt_350_150.get1DHistoFromDraw( "dl_mt2ll", selectionString=presel, binning = [0, 25,50,75,100,140,240,340], binningIsExplicit=True, addOverFlowBin = 'upper', weightString = "60 * weight * reweightDilepTrigger * reweightLeptonSF * reweightBTag_SF * reweightPUVUp" )
#
#h_200_central   = T2tt_200_50.get1DHistoFromDraw( "dl_mt2ll", selectionString=presel, binning = [0, 25,50,75,100,140,240,340], binningIsExplicit=True, addOverFlowBin = 'upper', weightString  = "60 * weight * reweightDilepTrigger * reweightLeptonSF * reweightBTag_SF * reweightPU" )
#h_200_vup       = T2tt_200_50.get1DHistoFromDraw( "dl_mt2ll", selectionString=presel, binning = [0, 25,50,75,100,140,240,340], binningIsExplicit=True, addOverFlowBin = 'upper', weightString  = "60 * weight * reweightDilepTrigger * reweightLeptonSF * reweightBTag_SF * reweightPUVUp" )
#
#h_800_central.style = styles.lineStyle( ROOT.kRed+2,  width=2, errors=True )
#h_800_vup.style     = styles.lineStyle( ROOT.kBlue+1, width=2, errors=True )
#
#h_350_central.style = styles.lineStyle( ROOT.kOrange+1, width=2, errors=True )
#h_350_vup.style     = styles.lineStyle( ROOT.kGreen+1,  width=2, errors=True )
#
#h_200_central.style = styles.lineStyle( ROOT.kCyan+1, width=2, errors=True )
#h_200_vup.style     = styles.lineStyle( ROOT.kMagenta+1,  width=2, errors=True )
#
#h_800_central.legendText = "T2tt(800,100), central"
#h_800_vup.legendText     = "T2tt(800,100), VUp"
#
#h_350_central.legendText = "T2tt(350,150), central"
#h_350_vup.legendText     = "T2tt(350,150), VUp"
#
#h_200_central.legendText = "T2tt(200,50), central"
#h_200_vup.legendText     = "T2tt(200,50), VUp"


h_T2tt_750   = T2tt_750_0.get1DHistoFromDraw( "dl_mt2ll", selectionString=presel, binning = [0, 25,50,75,100,140,240,340], binningIsExplicit=True, addOverFlowBin = 'upper', weightString = "60 * weight * reweightDilepTrigger * reweightLeptonSF * reweightBTag_SF * reweightPU" )
h_TTDM       = TTDM_scalar_400_1.get1DHistoFromDraw( "dl_mt2ll", selectionString=presel, binning = [0, 25,50,75,100,140,240,340], binningIsExplicit=True, addOverFlowBin = 'upper', weightString = "60 * weight * reweightDilepTrigger * reweightLeptonSF * reweightBTag_SF * reweightPU" )

h_T2tt_750.style = styles.lineStyle( ROOT.kRed+2,  width=2, errors=True )
h_T2tt_750.legendText = "T2tt(750,0)"

h_TTDM.style = styles.lineStyle( ROOT.kBlue+1,  width=2, errors=True )
h_TTDM.legendText = "tt+DM, scalar 400,1"


h_met_T2tt_750   = T2tt_750_0.get1DHistoFromDraw( "met_pt", selectionString=presel, binning = [20, 0, 500], addOverFlowBin = 'upper', weightString = "60 * weight * reweightDilepTrigger * reweightLeptonSF * reweightBTag_SF * reweightPU" )
h_met_TTDM       = TTDM_scalar_400_1.get1DHistoFromDraw( "met_pt", selectionString=presel, binning = [20, 0, 500], addOverFlowBin = 'upper', weightString = "60 * weight * reweightDilepTrigger * reweightLeptonSF * reweightBTag_SF * reweightPU" )

h_met_T2tt_750.style = styles.lineStyle( ROOT.kRed+2,  width=2, errors=True )
h_met_T2tt_750.legendText = "T2tt(750,0)"

h_met_TTDM.style = styles.lineStyle( ROOT.kBlue+1,  width=2, errors=True )
h_met_TTDM.legendText = "tt+DM, scalar 400,1"



def drawObjects( lumi ):
    tex = ROOT.TLatex()
    tex.SetNDC()
    tex.SetTextSize(0.04)
    tex.SetTextAlign(11) # align right

    lines = [ (0.15, 0.95, 'CMS Simulation') ]
    lines.append( (0.65, 0.95, 'L=%3.1f fb{}^{-1} (13 TeV)'% ( lumi ) ) )
    return [tex.DrawLatex(*l) for l in lines]

plotting.draw(
    Plot.fromHisto("dl_mt2ll",
                #[[ h_800_central ], [ h_800_vup ], [h_350_central], [h_350_vup], [h_200_central], [h_200_vup]],
                [[ h_T2tt_750], [h_TTDM]],
                texX = "M_{T}^{ll} (GeV)"
            ),
    plot_directory = "/afs/hephy.at/user/d/dspitzbart/www/stopsDileptonLegacy/signal_ttDM/",
    logX = False, logY = True, #sorting = True, 
    yRange = (0.03, 300),
    drawObjects = drawObjects(60.),
    scaling = {0:1}
)

plotting.draw(
    Plot.fromHisto("met_pt_lin",
                #[[ h_800_central ], [ h_800_vup ], [h_350_central], [h_350_vup], [h_200_central], [h_200_vup]],
                [[ h_met_T2tt_750], [h_met_TTDM]],
                texX = "p_{T}^{miss} (GeV)"
            ),
    plot_directory = "/afs/hephy.at/user/d/dspitzbart/www/stopsDileptonLegacy/signal_ttDM/",
    logX = False, logY = False, #sorting = True, 
    #yRange = (0.03, 30),
    yRange = (0.03, 3),
    drawObjects = drawObjects(60.),
    scaling = {0:1}
)



