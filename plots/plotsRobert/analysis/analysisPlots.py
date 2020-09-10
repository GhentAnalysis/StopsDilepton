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
import array
import operator

from math                                import sqrt, cos, sin, pi, atan2, cosh
from RootTools.core.standard             import *
from StopsDilepton.tools.user            import plot_directory
from StopsDilepton.tools.helpers         import deltaPhi, deltaR, add_histos
from Analysis.Tools.metFilters            import getFilterCut
from StopsDilepton.tools.cutInterpreter  import cutInterpreter
from StopsDilepton.tools.mt2Calculator   import mt2Calculator
from Analysis.Tools.puProfileCache import *

#
# Arguments
# 
import argparse
argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--logLevel',           action='store',      default='INFO',          nargs='?', choices=['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'TRACE', 'NOTSET'], help="Log level for logging")
argParser.add_argument('--signal',             action='store',      default=None,            nargs='?', choices=[None, "T2tt"], help="Add signal to plot")
argParser.add_argument('--noData',             action='store_true', default=False,           help='also plot data?')
argParser.add_argument('--small',                                 action='store_true',     help='Run only on a small subset of the data?', )
#argParser.add_argument('--sorting',                               action='store', default=None, choices=[None, "forDYMB"],  help='Sort histos?', )
argParser.add_argument('--dpm',                                   action='store_true',     help='Use dpm?', )
argParser.add_argument('--dataMCScaling',      action='store_true',     help='Data MC scaling?', )
argParser.add_argument('--DYInc',              action='store_true',     help='Use Inclusive DY sample?', )
argParser.add_argument('--plot_directory',     action='store',      default='v0p3')
argParser.add_argument('--era',                action='store', type=str,      default="2016")
argParser.add_argument('--selection',          action='store',      default='lepSel-njet2p-btag0-looseLeptonVeto-mll20-dPhiJet0-dPhiJet1')
argParser.add_argument('--nvtxReweightSelection',          action='store',      default=None)
argParser.add_argument('--unblinded',          action='store_true', default=False)
argParser.add_argument('--blinded',            action='store_true', default=False)
argParser.add_argument('--reweightPU',         action='store', default='Central', choices=['VDown', 'Down', 'Central', 'Up', 'VUp', 'VVUp', 'noPUReweighting', 'nvtx'])
argParser.add_argument('--isr',                action='store_true', default=False)
argParser.add_argument('--plotUPara',          action='store_true',     help='Plot u_para?', )
#argParser.add_argument('--splitDiBoson',       action='store_true',     help='Split diBoson?' )
argParser.add_argument('--splitMET',           action='store_true',     help='Split in MET bins?' )
argParser.add_argument('--splitMETSig',        action='store_true',     help='Split in METSig bins?' )
argParser.add_argument('--splitNvtx',          action='store_true',     help='Split in Nvtx bins?' )
argParser.add_argument('--rwHit0',              action='store_true',     help='reweight Hit0?', )
argParser.add_argument('--rwSip3d',              action='store_true',     help='reweight Sip3d?', )
args = argParser.parse_args()

#
# Logger
#
import StopsDilepton.tools.logger as logger
import RootTools.core.logger as logger_rt
logger    = logger.get_logger(   args.logLevel, logFile = None)
logger_rt = logger_rt.get_logger(args.logLevel, logFile = None)

if args.small:                        args.plot_directory += "_small"
if args.splitMET:                     args.plot_directory += "_splitMET"
if args.splitMETSig:                  args.plot_directory += "_splitMETSig"
if args.splitNvtx:                    args.plot_directory += "_splitNvtx"
if args.DYInc:                        args.plot_directory += "_DYInc"
if args.noData:                       args.plot_directory += "_noData"
if args.rwHit0:                       args.plot_directory += "_rwHit0"
if args.rwSip3d:                      args.plot_directory += "_rwSip3d"
if args.signal == "DM":               args.plot_directory += "_DM"
if args.reweightPU:                   args.plot_directory += "_%s"%args.reweightPU
#
# Make samples, will be searched for in the postProcessing directory
#
from Analysis.Tools.puReweighting import getReweightingFunction

if "2016" in args.era:
    year = 2016
elif "2017" in args.era:
    year = 2017
elif "2018" in args.era:
    year = 2018

logger.info( "Working in year %i", year )

# Load from DPM?
if args.dpm:
    data_directory = "/dpm/oeaw.ac.at/home/cms/store/user/rschoefbeck/Stops2l-postprocessed/"

if year == 2016:
    from StopsDilepton.samples.nanoTuples_Run2016_17Jul2018_postProcessed import *
    #data_directory          = "/afs/hephy.at/data/rschoefbeck02/nanoTuples/"
    #postProcessing_directory= "stops_2016_nano_v0p15/dilep/"
    from StopsDilepton.samples.nanoTuples_Summer16_postProcessed import *
    if args.DYInc:
        mc             = [ Top_pow_16, TTXNoZ_16, TTZ_16, multiBoson_16, DY_LO_16]
    else:
        mc             = [ Top_pow_16, TTXNoZ_16, TTZ_16, multiBoson_16, DY_HT_LO_16]
    #if args.reweightPU and not args.reweightPU in ["noPUReweighting", "nvtx"]:
    #    nTrueInt_puRW = getReweightingFunction(data="PU_2016_35920_XSec%s"%args.reweightPU, mc="Summer16")
elif year == 2017:
    from StopsDilepton.samples.nanoTuples_Fall17_postProcessed import *
    #data_directory           = "/afs/hephy.at/data/rschoefbeck02/nanoTuples/"
    #postProcessing_directory = "stops_2017_nano_v0p12_JECV6/dilep/"
    from StopsDilepton.samples.nanoTuples_Run2017_31Mar2018_postProcessed import *

    if args.DYInc:
        mc             = [ Top_pow_17, TTXNoZ_17, TTZ_17, multiBoson_17, DY_LO_17]
    else:
        mc             = [ Top_pow_17, TTXNoZ_17, TTZ_17, multiBoson_17, DY_HT_LO_17]
    #if args.reweightPU:
    #    # need sample based weights
    #    pass
elif year == 2018:
    from StopsDilepton.samples.nanoTuples_Run2018_PromptReco_postProcessed import *
    #data_directory          = "/afs/hephy.at/data/rschoefbeck02/nanoTuples/"
    #postProcessing_directory= "stops_2018_nano_v0p15/dilep/"
    from StopsDilepton.samples.nanoTuples_Autumn18_postProcessed import *
    if args.DYInc:
        mc             = [ Top_pow_18, TTXNoZ_18, TTZ_18, multiBoson_18, DY_LO_18]
    else:
        mc             = [ Top_pow_18, TTXNoZ_18, TTZ_18, multiBoson_18, DY_HT_LO_18]

    #from StopsDilepton.tools.vetoList import vetoList
    #Run2018D.vetoList = vetoList.fromDirectory('/afs/hephy.at/data/rschoefbeck02/StopsDilepton/splitMuonVeto/')
    #if args.reweightPU and not args.reweightPU in ["noPUReweighting", "nvtx"]:
    #    nTrueInt_puRW = getReweightingFunction(data="PU_2018_58830_XSec%s"%args.reweightPU, mc="Autumn18")

#if args.sorting == "forDYMB":
#    mc = [ mc[4], mc[3], mc[0], mc[2], mc[1] ]

try:
  data_sample = eval(args.era)
except Exception as e:
  logger.error( "Didn't find %s", args.era )
  raise e

lumi_scale                 = data_sample.lumi/1000
data_sample.scale          = 1.
for sample in mc:
    sample.scale          = lumi_scale

if args.small:
    for sample in mc + [data_sample]:
        sample.normalization = 1.
        #sample.reduceFiles( factor = 40 )
        sample.reduceFiles( to=1)
        sample.scale /= sample.normalization

#from Analysis.Tools.RecoilCorrector import RecoilCorrector
#recoilCorrector = RecoilCorrector( os.path.join( "/afs/hephy.at/data/rschoefbeck01/StopsDilepton/results/", "recoil_v0p10_fine", "%s_lepSel-njet1p-btag0-relIso0.12-looseLeptonVeto-mll20-onZ_recoil_fitResults_SF.pkl"%args.era ) )

def get_quantiles( histo, quantiles = [1-0.9545, 1-0.6826, 0.5, 0.6826, 0.9545]):
    thresholds = array.array('d', [ROOT.Double()] * len(quantiles) )
    histo.GetQuantiles( len(quantiles), thresholds, array.array('d', quantiles) )
    return thresholds 

def splitMetMC(mc):
    dy = mc[-1]
    dy_1 = copy.deepcopy( dy )
    dy_1.name += "_1"
    dy_1.addSelectionString( "met_pt<40" )
    dy_1.texName += " (MET<40)"
    dy_1.color   = ROOT.kGreen + 1 
    dy_2 = copy.deepcopy( dy )
    dy_2.name += "_2"
    dy_2.addSelectionString( "met_pt>40&&met_pt<80" )
    dy_2.texName += " (40<MET<80)"
    dy_2.color   = ROOT.kGreen + 2
    dy_3 = copy.deepcopy( dy )
    dy_3.name += "_3"
    dy_3.addSelectionString( "met_pt>80" )
    dy_3.texName += " (80<MET)"
    dy_3.color   = ROOT.kGreen + 3
    tt = mc[0]
    tt_1 = copy.deepcopy( tt )
    tt_1.name += "_1"
    tt_1.addSelectionString( "met_pt<40" )
    tt_1.texName += " (MET<40)"
    tt_1.color   = ROOT.kAzure + 1 
    tt_2 = copy.deepcopy( tt )
    tt_2.name += "_2"
    tt_2.addSelectionString( "met_pt>40&&met_pt<80" )
    tt_2.texName += " (40<MET<80)"
    tt_2.color   = ROOT.kAzure + 2
    tt_3 = copy.deepcopy( tt )
    tt_3.name += "_3"
    tt_3.addSelectionString( "met_pt>80" )
    tt_3.texName += " (80<MET)"
    tt_3.color   = ROOT.kAzure + 3

    return [ dy_1, dy_2, dy_3, tt_1, tt_2, tt_3] + mc[1:-1]

def splitMetSigMC(mc):
    dy = mc[-1]
    dy_1 = copy.deepcopy( dy )
    dy_1.name += "_1"
    dy_1.addSelectionString( "MET_significance<6" )
    dy_1.texName += " (S<6)"
    dy_1.color   = ROOT.kGreen + 1 
    dy_2 = copy.deepcopy( dy )
    dy_2.name += "_2"
    dy_2.addSelectionString( "MET_significance>6&&MET_significance<12" )
    dy_2.texName += " (6<S<12)"
    dy_2.color   = ROOT.kGreen + 2
    dy_3 = copy.deepcopy( dy )
    dy_3.name += "_3"
    dy_3.addSelectionString( "MET_significance>12" )
    dy_3.texName += " (12<S)"
    dy_3.color   = ROOT.kGreen + 3
    tt = mc[0]
    tt_1 = copy.deepcopy( tt )
    tt_1.name += "_1"
    tt_1.addSelectionString( "MET_significance<6" )
    tt_1.texName += " (S<6)"
    tt_1.color   = ROOT.kAzure + 1 
    tt_2 = copy.deepcopy( tt )
    tt_2.name += "_2"
    tt_2.addSelectionString( "MET_significance>6&&MET_significance<12" )
    tt_2.texName += " (6<S<12)"
    tt_2.color   = ROOT.kAzure + 2
    tt_3 = copy.deepcopy( tt )
    tt_3.name += "_3"
    tt_3.addSelectionString( "MET_significance>12" )
    tt_3.texName += " (12<S)"
    tt_3.color   = ROOT.kAzure + 3

    return [ dy_1, dy_2, dy_3, tt_1, tt_2, tt_3] + mc[1:-1]

def splitNvtxMC(mc):
    dy = mc[-1]
    dy_1 = copy.deepcopy( dy )
    dy_1.name += "_1"
    dy_1.addSelectionString( "PV_npvsGood<25" )
    dy_1.texName += " (N_{vtx}<25)"
    dy_1.color   = ROOT.kGreen + 1 
    dy_2 = copy.deepcopy( dy )
    dy_2.name += "_2"
    dy_2.addSelectionString( "PV_npvsGood>=25&&PV_npvsGood<35" )
    dy_2.texName += " (25#leq N_{vtx}<35)"
    dy_2.color   = ROOT.kGreen + 2
    dy_3 = copy.deepcopy( dy )
    dy_3.name += "_3"
    dy_3.addSelectionString( "PV_npvsGood>=35" )
    dy_3.texName += " (35#leq N_{vtx})"
    dy_3.color   = ROOT.kGreen + 3
    tt = mc[0]
    tt_1 = copy.deepcopy( tt )
    tt_1.name += "_1"
    tt_1.addSelectionString( "PV_npvsGood<25" )
    tt_1.texName += " (N_{vtx}<25)"
    tt_1.color   = ROOT.kAzure + 1 
    tt_2 = copy.deepcopy( tt )
    tt_2.name += "_2"
    tt_2.addSelectionString( "PV_npvsGood>=25&&PV_npvsGood<35" )
    tt_2.texName += " (25#leq N_{vtx}<35)"
    tt_2.color   = ROOT.kAzure + 2
    tt_3 = copy.deepcopy( tt )
    tt_3.name += "_3"
    tt_3.addSelectionString( "PV_npvsGood>=35" )
    tt_3.texName += " (35#leq N_{vtx})"
    tt_3.color   = ROOT.kAzure + 3

    return [ dy_1, dy_2, dy_3, tt_1, tt_2, tt_3] + mc[1:-1]

signals = []
if args.signal == "T2tt":
    # Load 2017 signal
    #data_directory           = "/afs/hephy.at/data/dspitzbart03/nanoTuples/"
    #postProcessing_directory = "stops_2017_nano_v0p9/dilep"
    if args.era=="Run2016":
        from StopsDilepton.samples.nanoTuples_FastSim_Summer16_postProcessed import *
    elif args.era == "Run2017":
        from StopsDilepton.samples.nanoTuples_FastSim_Fall17_postProcessed import *
    elif args.era == "Run2018":
        from StopsDilepton.samples.nanoTuples_FastSim_Autumn18_postProcessed import *

    T2tt                    = T2tt_650_0
    T2tt2                   = T2tt_500_250
    T2tt2.style             = styles.lineStyle( ROOT.kBlack, width=3, dotted=True )
    T2tt.style              = styles.lineStyle( ROOT.kBlack, width=3 )
    signals = [ T2tt, T2tt2]
#
# Text on the plots
#
tex = ROOT.TLatex()
tex.SetNDC()
tex.SetTextSize(0.04)
tex.SetTextAlign(11) # align right

def drawObjects( plotData, dataMCScale, lumi_scale ):
    lines = [
      (0.15, 0.95, 'CMS Preliminary' if plotData else 'CMS Simulation'), 
      (0.45, 0.95, 'L=%3.1f fb{}^{-1} (13 TeV) Scale %3.2f'% ( lumi_scale, dataMCScale ) ) if plotData else (0.45, 0.95, 'L=%3.1f fb{}^{-1} (13 TeV)' % lumi_scale)
    ]
    if "mt2ll100" in args.selection and args.noData: lines += [(0.55, 0.5, 'M_{T2}(ll) > 100 GeV')] # Manually put the mt2ll > 100 GeV label
    return [tex.DrawLatex(*l) for l in lines] 

def drawPlots(plots, mode, dataMCScale):
  for log in [False, True]:
    plot_directory_ = os.path.join(plot_directory, 'analysisPlots', args.plot_directory, args.era, mode + ("_log" if log else ""), args.selection)
    for plot in plots:
      if not max(l[0].GetMaximum() for l in plot.histos): 
        logger.warning( "Plot %s apparently empty. Skipping.", plot.name )
        continue # Empty plot
      if not args.noData: 
        if mode == "all": plot.histos[1][0].legendText = "Data"
        if mode == "SF":  plot.histos[1][0].legendText = "Data (SF)"

      _drawObjects = []

      if ("u_para" in plot.name or "u_perp" in plot.name) and not args.noData:
          h_mc   = plot.histos_added[0][0].Clone()
          h_data = plot.histos_added[1][0].Clone()
          if h_mc.Integral()>0:
              h_mc.Scale(h_data.Integral()/h_mc.Integral())
          q_mc   = tuple(get_quantiles( h_mc ))
          q_data = tuple(get_quantiles( h_data ))
          median_shift = q_data[2]-q_mc[2]
          sigma1_ratio = (q_data[3]-q_data[1])/(q_mc[3]-q_mc[1]) if q_mc[3]-q_mc[1]!=0 else 0
          sigma2_ratio = (q_data[4]-q_data[0])/(q_mc[4]-q_mc[0]) if q_mc[4]-q_mc[0]!=0 else 0

          _drawObjects.append( tex.DrawLatex(0.22, 0.62, '#Delta(med): %+3.1f   1#sigma: %4.3f  2#sigma  %4.3f' % ( median_shift, sigma1_ratio, sigma2_ratio) ) )

      if isinstance( plot, Plot):
          plotting.draw(plot,
            plot_directory = plot_directory_,
            ratio = {'yRange':(0.1,1.9)} if not args.noData else None,
            logX = False, logY = log, sorting = not (args.splitMET or args.splitMETSig or args.splitNvtx), # and args.sorting is not None,
            yRange = (0.03, "auto") if log else (0.001, "auto"),
            scaling = {0:1} if args.dataMCScaling else {},
            legend = ( (0.18,0.88-0.03*sum(map(len, plot.histos)),0.9,0.88), 2),
            drawObjects = drawObjects( not args.noData, dataMCScale , lumi_scale ) + _drawObjects,
            copyIndexPHP = True, extensions = ["png", "pdf"],
          )
      elif isinstance( plot, Plot2D ):

          p_mc = Plot2D.fromHisto( plot.name+'_mc', plot.histos[:1], texX = plot.texX, texY = plot.texY )
          plotting.draw2D(p_mc,
            plot_directory = plot_directory_,
            #ratio = {'yRange':(0.1,1.9)},
            logX = False, logY = False, logZ = log, #sorting = True,
            #yRange = (0.03, "auto") if log else (0.001, "auto"),
            #scaling = {},
            #legend = (0.50,0.88-0.04*sum(map(len, plot.histos)),0.9,0.88),
            drawObjects = drawObjects( not args.noData, dataMCScale , lumi_scale ),
            copyIndexPHP = True, extensions = ["png", "pdf"], 
          )
          p_data = Plot2D.fromHisto( plot.name+'_data', plot.histos[1:], texX = plot.texX, texY = plot.texY )
          plotting.draw2D(p_data,
            plot_directory = plot_directory_,
            #ratio = {'yRange':(0.1,1.9)},
            logX = False, logY = False, logZ = log, #sorting = True,
            #yRange = (0.03, "auto") if log else (0.001, "auto"),
            #scaling = {},
            #legend = (0.50,0.88-0.04*sum(map(len, plot.histos)),0.9,0.88),
            drawObjects = drawObjects( not args.noData, dataMCScale , lumi_scale ),
            copyIndexPHP = True, extensions = ["png", "pdf"], 
          )
          # 2D plots where we want to make a ratio start with "forRatio"
          if plot.name.startswith("eff_num"):
            # find the denominator plot
            plot_den = filter( lambda p:p.name == plot.name.replace("eff_num", "eff_den"), plots)[0]
            h_eff_data   = plot.histos[1][0] 
            h_den_data   = plot_den.histos[1][0]
            h_eff_data.Divide(h_den_data)
            h_eff_mc     = add_histos(plot.histos[0])
            h_den_mc     = add_histos(plot_den.histos[0])
            h_eff_mc.Divide(h_den_mc)

            p_eff_mc =  Plot2D.fromHisto( plot.name.replace('eff_num', 'eff_mc'), [[h_eff_mc]], texX = plot.texX, texY = plot.texY )
            plotting.draw2D(p_eff_mc,
              plot_directory = plot_directory_,
              zRange = (0.8,1.2),
              logX = False, logY = False, logZ = log, #sorting = True,
              drawObjects = drawObjects( not args.noData, dataMCScale , lumi_scale ),
              copyIndexPHP = True, extensions = ["png", "pdf", "root"], 
            )
            p_eff_data =  Plot2D.fromHisto( plot.name.replace('eff_num', 'eff_data'), [[h_eff_data]], texX = plot.texX, texY = plot.texY )
            plotting.draw2D(p_eff_data,
              plot_directory = plot_directory_,
              logX = False, logY = False, logZ = log, #sorting = True,
              zRange = (0.8,1.2),
              drawObjects = drawObjects( not args.noData, dataMCScale , lumi_scale ),
              copyIndexPHP = True, extensions = ["png", "pdf", "root"], 
            )
            h_eff_data.Divide(h_eff_mc) 
            p_ratio_datamc =  Plot2D.fromHisto( plot.name.replace('eff_num', 'effRatio_dataMC'), [[h_eff_data]], texX = plot.texX, texY = plot.texY )
            plotting.draw2D(p_ratio_datamc,
              plot_directory = plot_directory_,
              logX = False, logY = False, logZ = log, #sorting = True,
              zRange = (0.8,1.2),
              drawObjects = drawObjects( not args.noData, dataMCScale , lumi_scale ),
              copyIndexPHP = True, extensions = ["png", "pdf", "root"], 
            )
            
#
# Read variables and sequences
#
read_variables  = ["weight/F", "l1_pt/F", "dl_phi/F", "dl_pt/F", "l2_pt/F", "l1_eta/F" , "l1_phi/F", "l2_eta/F", "l2_phi/F", "JetGood[pt/F,eta/F,phi/F]", "dl_mass/F", "dl_eta/F", "dl_mt2ll/F", "dl_mt2bb/F", "dl_mt2blbl/F", "met_pt/F", "met_phi/F", "MET_significance/F", "nBTag/I", "nJetGood/I", "PV_npvsGood/I", "RawMET_pt/F", "RawMET_phi/F"]
read_variables += ["l2_eleIndex/I"]
read_variables += ["event/l", "luminosityBlock/I", "run/I"]
if "2017" in args.era:
    read_variables.append( "MET_pt_min/F" ) 

sequence = []

if True: #"2017" in args.era:
    from StopsDilepton.tools.objectSelection    import muonSelector, eleSelector, getGoodMuons, getGoodElectrons, getGoodJets, getAllJets, isBJet
    ele_selector = eleSelector( "tight", year = year )
    mu_selector = muonSelector( "tight", year = year )

    jetVars         = ['pt/F', 'chEmEF/F', 'chHEF/F', 'neEmEF/F', 'neHEF/F', 'rawFactor/F', 'eta/F', 'phi/F', 'jetId/I', 'btagDeepB/F', 'btagCSVV2/F', 'area/F', 'pt_nom/F', 'corr_JER/F'] 
    jetVarNames     = map( lambda s:s.split('/')[0], jetVars)
    read_variables += [\
        TreeVariable.fromString('nElectron/I'),
        VectorTreeVariable.fromString('Electron[pt/F,eta/F,phi/F,pdgId/I,cutBased/I,miniPFRelIso_all/F,pfRelIso03_all/F,sip3d/F,lostHits/b,convVeto/O,dxy/F,dz/F,charge/I,deltaEtaSC/F]'),
        TreeVariable.fromString('nMuon/I'),
        VectorTreeVariable.fromString('Muon[pt/F,eta/F,phi/F,pdgId/I,mediumId/O,miniPFRelIso_all/F,pfRelIso03_all/F,sip3d/F,dxy/F,dz/F,charge/I]'),
        TreeVariable.fromString('nJet/I'),
        VectorTreeVariable.fromString('Jet[%s]'% ( ','.join(jetVars) ) ),
    ]

    def make_all_jets( event, sample ):

        electrons_pt10  = getGoodElectrons(event, ele_selector = ele_selector)
        muons_pt10      = getGoodMuons(event, mu_selector = mu_selector )
        for e in electrons_pt10:
            e['pdgId']      = int( -11*e['charge'] )
        for m in muons_pt10:
            m['pdgId']      = int( -13*m['charge'] )
        leptons_pt10 = electrons_pt10+muons_pt10
        leptons_pt10.sort(key = lambda p:-p['pt'])

        leptons      = filter(lambda l:l['pt']>20, leptons_pt10)
        leptons.sort(key = lambda p:-p['pt'])
        event.jets  = getAllJets(event, leptons, ptCut=30, absEtaCut=99, jetVars=jetVarNames, jetCollections=["Jet"], idVar='jetId')     
        event.nJet_EE = len(filter(lambda j:abs(j['eta'])>2.6 and abs(j['eta'])<3.1, event.jets))
        event.nJet_EE_pt30To50 = len(filter(lambda j:abs(j['eta'])>2.6 and abs(j['eta'])<3.1 and j['pt']>30 and j['pt']<50, event.jets))
        event.nJet_EE_ptTo50 = len(filter(lambda j:abs(j['eta'])>2.6 and abs(j['eta'])<3.1 and j['pt']<50, event.jets))
        event.nJet_EE_ptTo40 = len(filter(lambda j:abs(j['eta'])>2.6 and abs(j['eta'])<3.1 and j['pt']<40, event.jets))
        event.nJet_EE_ptTo30 = len(filter(lambda j:abs(j['eta'])>2.6 and abs(j['eta'])<3.1 and j['pt']<30, event.jets))
        event.nJet_EE_pt30 = len(filter(lambda j:abs(j['eta'])>2.6 and abs(j['eta'])<3.1 and j['pt']>30, event.jets))
        event.nJet_EE_pt40 = len(filter(lambda j:abs(j['eta'])>2.6 and abs(j['eta'])<3.1 and j['pt']>40, event.jets))
        event.nJet_EE_pt50 = len(filter(lambda j:abs(j['eta'])>2.6 and abs(j['eta'])<3.1 and j['pt']>50, event.jets))

        event.badJetE  = sum( [ j['pt']*j['neEmEF']*cosh(j['eta']) for j in event.jets if abs(j['eta'])>2.6 and abs(j['eta'])<3.1], 0. )
        event.badJetPt = sum( [ j['pt']*j['neEmEF'] for j in event.jets if abs(j['eta'])>2.6 and abs(j['eta'])<3.1], 0. )

        bJets        = filter(lambda j:isBJet(j, tagger="DeepCSV", year=year) and abs(j['eta'])<=2.4    , event.jets)
        nonBJets     = filter(lambda j:not ( isBJet(j, tagger="DeepCSV", year=year) and abs(j['eta'])<=2.4 ), event.jets)

        event.bj0, event.bj1 = None, None
        if len(bJets)+len(nonBJets)>=2:
            event.bj0, event.bj1  = (bJets+nonBJets)[:2]
            event.mt2blbl_bj_dR   = deltaR  ( event.bj0, event.bj1 ) 
            event.mt2blbl_bj_dPhi = deltaPhi( event.bj0['phi'], event.bj1['phi'] ) 
            event.mt2blbl_bj_dEta = event.bj0['eta'] - event.bj1['eta']
            event.mt2blbl_bj_mass = sqrt(2.*event.bj0['pt']*event.bj1['pt']*(cosh(event.bj0['eta']-event.bj1['eta'])-cos(event.bj0['phi']-event.bj1['phi']))) 
        else:
            event.mt2blbl_bj_dR   = float('nan') 
            event.mt2blbl_bj_dPhi = float('nan') 
            event.mt2blbl_bj_dEta = float('nan') 
            event.mt2blbl_bj_mass = float('nan') 
    
    sequence.append( make_all_jets )

## veto list
#def make_veto( event, sample ):
#    if hasattr( sample, "vetoList" ):
#        event.passes_veto = sample.vetoList.passesVeto( event.run, event.luminosityBlock, event.event )
#        #event.passes_veto = event.event%2==0
#
#sequence.append( make_veto )

## recoil
#def recoil_weight( var_bin, qt_bin):
#    #if args.recoil == 'v4':
#    def _weight_( event, sample):
#        return event.weight*(event.dl_phi>var_bin[0])*(event.dl_phi<=var_bin[1])*(event.dl_pt>qt_bin[0])*(event.dl_pt<qt_bin[1]) 
#    #elif args.recoil == 'v5':
#    #    def _weight_( event, sample):
#    #        return event.weight*(event.PV_npvsGood>var_bin[0])*(event.PV_npvsGood<=var_bin[1])*(event.dl_pt>qt_bin[0])*(event.dl_pt<qt_bin[1]) 
#    return _weight_
#
#def corr_recoil( event, sample ):
#    mt2Calculator.reset()
#    if not sample.isData: 
#        # Parametrisation vector - # define qt as GenMET + leptons
#        qt_px = event.l1_pt*cos(event.l1_phi) + event.l2_pt*cos(event.l2_phi) + event.GenMET_pt*cos(event.GenMET_phi)
#        qt_py = event.l1_pt*sin(event.l1_phi) + event.l2_pt*sin(event.l2_phi) + event.GenMET_pt*sin(event.GenMET_phi)
#
#        qt = sqrt( qt_px**2 + qt_py**2 )
#        qt_phi = atan2( qt_py, qt_px )
#
#        #ref_phi = qt_phi
#        ref_phi = event.dl_phi
#
#        # recoil-correct met_pt, met_phi
#        # compute fake MET 
#        fakeMET_x = event.met_pt*cos(event.met_phi) - event.GenMET_pt*cos(event.GenMET_phi)
#        fakeMET_y = event.met_pt*sin(event.met_phi) - event.GenMET_pt*sin(event.GenMET_phi)
#        fakeMET = sqrt( fakeMET_x**2 + fakeMET_y**2 )
#        fakeMET_phi = atan2( fakeMET_y, fakeMET_x )
#        # project fake MET on qT
#        fakeMET_para = fakeMET*cos( fakeMET_phi - ref_phi ) 
#        fakeMET_perp = fakeMET*cos( fakeMET_phi - ( ref_phi - pi/2) ) 
#        fakeMET_para_corr = - recoilCorrector.predict_para( ref_phi, qt, -fakeMET_para ) 
#        fakeMET_perp_corr = - recoilCorrector.predict_perp( ref_phi, qt, -fakeMET_perp )
#        # rebuild fake MET vector
#        fakeMET_px_corr = fakeMET_para_corr*cos(ref_phi) + fakeMET_perp_corr*cos(ref_phi - pi/2) 
#        fakeMET_py_corr = fakeMET_para_corr*sin(ref_phi) + fakeMET_perp_corr*sin(ref_phi - pi/2) 
#        #print "%s qt: %3.2f para %3.2f->%3.2f perp %3.2f->%3.2f fakeMET(%3.2f,%3.2f) -> (%3.2f,%3.2f)" % ( sample.name, qt, fakeMET_para, fakeMET_para_corr, fakeMET_perp, fakeMET_perp_corr, fakeMET, fakeMET_phi, sqrt( fakeMET_px_corr**2+fakeMET_py_corr**2), atan2( fakeMET_py_corr, fakeMET_px_corr) )
#        met_px_corr = event.met_pt*cos(event.met_phi) - fakeMET_x + fakeMET_px_corr 
#        met_py_corr = event.met_pt*sin(event.met_phi) - fakeMET_y + fakeMET_py_corr
#        event.met_pt_corr  = sqrt( met_px_corr**2 + met_py_corr**2 ) 
#        event.met_phi_corr = atan2( met_py_corr, met_px_corr ) 
#
#        ## recoil-correct RawMET_pt, RawMET_phi
#        ## compute fake MET 
#        #rawFakeMET_x = event.RawMET_pt*cos(event.RawMET_phi) - event.GenMET_pt*cos(event.GenMET_phi)
#        #rawFakeMET_y = event.RawMET_pt*sin(event.RawMET_phi) - event.GenMET_pt*sin(event.GenMET_phi)
#        #rawFakeMET = sqrt( rawFakeMET_x**2 + rawFakeMET_y**2 )
#        #rawFakeMET_phi = atan2( rawFakeMET_y, rawFakeMET_x )
#        ## project fake MET on qT
#        #rawFakeMET_para = rawFakeMET*cos( rawFakeMET_phi - ref_phi ) 
#        #rawFakeMET_perp = rawFakeMET*cos( rawFakeMET_phi - ( ref_phi - pi/2) ) 
#        #rawFakeMET_para_corr = - recoilCorrector_raw.predict_para( ref_phi, qt, -rawFakeMET_para ) 
#        #rawFakeMET_perp_corr = - recoilCorrector_raw.predict_perp( ref_phi, qt, -rawFakeMET_perp )
#        ## rebuild fake MET vector
#        #rawFakeMET_px_corr = rawFakeMET_para_corr*cos(ref_phi) + rawFakeMET_perp_corr*cos(ref_phi - pi/2) 
#        #rawFakeMET_py_corr = rawFakeMET_para_corr*sin(ref_phi) + rawFakeMET_perp_corr*sin(ref_phi - pi/2) 
#        #RawMET_px_corr = event.RawMET_pt*cos(event.RawMET_phi) - rawFakeMET_x + rawFakeMET_px_corr 
#        #RawMET_py_corr = event.RawMET_pt*sin(event.RawMET_phi) - rawFakeMET_y + rawFakeMET_py_corr
#        #event.RawMET_pt_corr  = sqrt( RawMET_px_corr**2 + RawMET_py_corr**2 ) 
#        #event.RawMET_phi_corr = atan2( RawMET_py_corr, RawMET_px_corr ) 
#
#    else:
#        event.met_pt_corr     = event.met_pt 
#        event.met_phi_corr    = event.met_phi
#        #event.RawMET_pt_corr  = event.RawMET_pt 
#        #event.RawMET_phi_corr = event.RawMET_phi
#
#    mt2Calculator.setLeptons(event.l1_pt, event.l1_eta, event.l1_phi, event.l2_pt, event.l2_eta, event.l2_phi)
#    mt2Calculator.setMet(event.met_pt_corr, event.met_phi_corr)
#    event.dl_mt2ll_corr     = mt2Calculator.mt2ll()
#
#    #mt2Calculator.setMet(event.RawMET_pt, event.RawMET_phi)
#    #event.dl_mt2ll_raw      =  mt2Calculator.mt2ll()
#
#    #mt2Calculator.setMet(event.RawMET_pt_corr, event.RawMET_phi_corr)
#    #event.dl_mt2ll_raw_corr =  mt2Calculator.mt2ll()
#
#    #print event.dl_mt2ll, event.dl_mt2ll_corr
#
#sequence.append( corr_recoil )
  
#
#
# default offZ for SF
offZ = "&&abs(dl_mass-91.1876)>15" if not (args.selection.count("onZ") or args.selection.count("allZ") or args.selection.count("offZ")) else ""
def getLeptonSelection( mode ):
  if   mode=="mumu": return "nGoodMuons==2&&nGoodElectrons==0&&isOS&&isMuMu" + offZ
  elif mode=="mue":  return "nGoodMuons==1&&nGoodElectrons==1&&isOS&&isEMu"
  elif mode=="ee":   return "nGoodMuons==0&&nGoodElectrons==2&&isOS&&isEE" + offZ
  elif mode=="SF":   return "nGoodMuons+nGoodElectrons==2&&isOS&&(isEE||isMuMu)" + offZ
  elif mode=="all":   return "nGoodMuons+nGoodElectrons==2&&isOS&&(((isEE||isMuMu)" + offZ+")||isEMu)"

# get nvtx reweighting histo
if args.reweightPU =='nvtx':
    logger.info( "Now obtain nvtx reweighting histo" )
    data_selectionString = "&&".join([getFilterCut(isData=True, year=year,skipVertexFilter=True), getLeptonSelection("SF"), cutInterpreter.cutString(args.nvtxReweightSelection)])
    data_nvtx_histo = data_sample.get1DHistoFromDraw( "PV_npvsGood", [100/5, 0, 100], selectionString=data_selectionString, weightString = "weight" )
    data_nvtx_histo.Scale(1./data_nvtx_histo.Integral())

    mc_selectionString = "&&".join([getFilterCut(isData=False, year=year,skipVertexFilter=True), getLeptonSelection("SF"), cutInterpreter.cutString(args.nvtxReweightSelection)])
    mc_histos  = [ s.get1DHistoFromDraw( "PV_npvsGood", [100/5, 0, 100], selectionString=mc_selectionString, weightString = "weight*reweightHEM*reweightDilepTrigger*reweightLeptonSF*reweightBTag_SF*reweightLeptonTrackingSF") for s in mc]
    mc_nvtx_histo     = mc_histos[0]
    for h in mc_histos[1:]:
        mc_nvtx_histo.Add( h )
    mc_nvtx_histo.Scale(1./mc_nvtx_histo.Integral())
    def nvtx_puRW( nvtx ):
        i_bin = mc_nvtx_histo.FindBin( nvtx )
        mc_val = mc_nvtx_histo.GetBinContent( i_bin)
        return data_nvtx_histo.GetBinContent( i_bin )/mc_val if mc_val>0 else 1

#
# Loop over channels
#
yields     = {}
allPlots   = {}
allModes   = ['mumu','mue','ee']
for index, mode in enumerate(allModes):
  yields[mode] = {}

  data_sample.setSelectionString([getFilterCut(isData=True, year=year,skipVertexFilter=True), getLeptonSelection(mode)])
  data_sample.name           = "data"
  data_sample.read_variables = ["event/I","run/I", "reweightHEM/F"]
  data_sample.style          = styles.errorStyle(ROOT.kBlack)
  weight_ = lambda event, sample: event.weight*event.reweightHEM

  #data_sample_filtered = copy.deepcopy( data_sample )
  #data_sample_filtered.style = styles.errorStyle(ROOT.kRed)
  #data_sample_filtered.weight = lambda event, sample: event.weight*event.passes_veto
  #data_sample_filtered.name   += "_filtered"
  #data_sample_filtered.texName+= " (filtered)"

  for sample in mc + signals:
    sample.read_variables = ['reweightPU/F', 'reweightL1Prefire/F', 'Pileup_nTrueInt/F', 'reweightDilepTrigger/F','reweightLeptonSF/F','reweightBTag_SF/F', 'reweightLeptonTrackingSF/F', 'GenMET_pt/F', 'GenMET_phi/F', 'reweightHEM/F', 'reweightLeptonHit0SF/F', 'reweightLeptonSip3dSF/F']
    # Need individual pu reweighting functions for each sample in 2017, so nTrueInt_puRW is only defined here

    if args.rwHit0:
        weight_Hit0  = lambda event: 1 
    else:
        weight_Hit0  = operator.attrgetter( 'reweightLeptonHit0SF' ) 
    if args.rwSip3d:
        weight_sip3d = lambda event: 1 
    else:
        weight_sip3d = operator.attrgetter( 'reweightLeptonSip3dSF' ) 

    if args.reweightPU and args.reweightPU not in ["noPUReweighting", "nvtx"]:
        sample.read_variables.append( 'reweightPU/F' if args.reweightPU=='Central' else 'reweightPU%s/F'%args.reweightPU )

    if args.reweightPU == "noPUReweighting":
        sample.weight         = lambda event, sample: event.reweightDilepTrigger*event.reweightLeptonSF*event.reweightBTag_SF*event.reweightLeptonTrackingSF*event.reweightL1Prefire
    elif args.reweightPU == "nvtx":
        sample.weight         = lambda event, sample: nvtx_puRW(event.PV_npvsGood) * event.reweightDilepTrigger*event.reweightLeptonSF*event.reweightBTag_SF*event.reweightLeptonTrackingSF*event.reweightL1Prefire
    elif args.reweightPU:
        pu_getter = operator.attrgetter( 'reweightPU' if args.reweightPU=='Central' else 'reweightPU%s'%args.reweightPU )
        sample.weight         = lambda event, sample: pu_getter(event) * event.reweightDilepTrigger*event.reweightLeptonSF*event.reweightBTag_SF*event.reweightLeptonTrackingSF*event.reweightL1Prefire*weight_sip3d(event)*weight_Hit0(event)
    else: #default
        sample.weight         = lambda event, sample: event.reweightPU*event.reweightDilepTrigger*event.reweightLeptonSF*event.reweightBTag_SF*event.reweightLeptonTrackingSF*event.reweightL1Prefire*weight_sip3d(event)*weight_Hit0(event)

    sample.setSelectionString([getFilterCut(isData=False, year=year,skipVertexFilter=True), getLeptonSelection(mode)])

  if args.splitMET:
    mc_ = splitMetMC(mc)
  elif args.splitMETSig:
    mc_ = splitMetSigMC(mc)
  elif args.splitNvtx:
    mc_ = splitNvtxMC(mc)
  else:
    mc_ = mc 

  for sample in mc_: sample.style = styles.fillStyle(sample.color)

  if not args.noData:
    stack = Stack(mc_, data_sample)#, data_sample_filtered)
  else:
    stack = Stack(mc_)

  stack.extend( [ [s] for s in signals ] )

  # Use some defaults
  Plot  .setDefaults(stack = stack, weight = staticmethod(weight_), selectionString = cutInterpreter.cutString(args.selection), addOverFlowBin='upper', histo_class=ROOT.TH1D)
  Plot2D.setDefaults(stack = stack, weight = staticmethod(weight_), selectionString = cutInterpreter.cutString(args.selection), histo_class=ROOT.TH2D)
  
  plots   = []
  plots2D = []
  plots.append(Plot(
    name = 'yield', texX = 'yield', texY = 'Number of Events',
    attribute = lambda event, sample: 0.5 + index,
    binning=[3, 0, 3],
  ))

  plots.append(Plot(
    name = 'PV_npvsGood', texX = 'N_{PV} (good)', texY = 'Number of Events',
    attribute = TreeVariable.fromString( "PV_npvsGood/I" ),
    binning=[100,0,100],
  ))
    
  plots.append(Plot(
    name = 'PV_npvs', texX = 'N_{PV} (total)', texY = 'Number of Events',
    attribute = TreeVariable.fromString( "PV_npvs/I" ),
    binning=[100,0,100],
  ))

  plots.append(Plot(
      texX = 'E_{T}^{miss} (GeV)', texY = 'Number of Events / 20 GeV',
      attribute = TreeVariable.fromString( "met_pt/F" ),
      binning=[400/20,0,400],
  ))

  plots2D.append(Plot2D(
    name = 'allJets_occupancy_pt30',
    stack = stack,
    attribute = (
      lambda event, sample: [ j['eta'] for j in event.jets if j['pt']>30 ],
      lambda event, sample: [ j['phi'] for j in event.jets if j['pt']>30 ],
    ),
    texX = '#eta (all jets)', texY = '#phi',
    binning=[52, -5.2, 5.2, 32, -pi, pi],
  ))

  plots.append(Plot(
    name = 'allJets_eta_pt30',
    stack = stack,
    attribute = lambda event, sample: [ j['eta'] for j in event.jets if j['pt']>30 ],
    texX = '#eta (all jets)', 
    binning=[52, -5.2, 5.2],
  ))

  plots.append(Plot(
    name = 'allJets_eta_pt50',
    stack = stack,
    attribute = lambda event, sample: [ j['eta'] for j in event.jets if j['pt']>50 ],
    texX = '#eta (all jets)', 
    binning=[52, -5.2, 5.2],
  ))

  plots.append(Plot(
    name = 'allJets_eta_pt100',
    stack = stack,
    attribute =  lambda event, sample: [ j['eta'] for j in event.jets if j['pt']>100 ],
    texX = '#eta (all jets)', 
    binning=[52, -5.2, 5.2],
  ))

  if "2017" in args.era:
    plots.append(Plot(
        texX = 'min E_{T}^{miss} (GeV)', texY = 'Number of Events / 20 GeV',
        attribute = TreeVariable.fromString( "MET_pt_min/F" ),
        binning=[400/20,0,400],
    ))
    plots.append(Plot(name = "MET_pt_min_delta",
        texX = '#Delta min E_{T}^{miss} (GeV)', texY = 'Number of Events / 10 GeV',
        attribute = lambda event, sample: event.met_pt - event.MET_pt_min,
        binning=[200/10,0,200],
    ))

#Sum$(abs(Jet_eta)>2.6&&abs(Jet_eta)<3.1&&Jet_pt<50)==0

  #plots.append(Plot( name = "met_pt_raw",
  #    texX = 'E_{T}^{miss} (GeV)', texY = 'Number of Events / 20 GeV',
  #    attribute = TreeVariable.fromString( "RawMET_pt/F" ),
  #    binning=[400/20,0,400],
  #))

  plots.append(Plot(
      texX = 'E_{T}^{miss} significance', texY = 'Number of Events',
      attribute = TreeVariable.fromString( "MET_significance/F" ),
      binning=[34,0,102],
  ))

  plots.append(Plot(
      texX = '#phi(E_{T}^{miss})', texY = 'Number of Events / 20 GeV',
      attribute = TreeVariable.fromString( "met_phi/F" ),
      binning=[10,-pi,pi],
  ))

  #plots.append(Plot( name = "met_phi_corr",
  #    texX = '#phi(E_{T}^{miss})', texY = 'Number of Events / 20 GeV',
  #    attribute = lambda event, sample: event.met_phi_corr,
  #    binning=[10,-pi,pi],
  #))

  #plots.append(Plot( name = "met_phi_raw",
  #    texX = 'raw #phi(E_{T}^{miss})', texY = 'Number of Events / 20 GeV',
  #    attribute = TreeVariable.fromString( "RawMET_phi/F" ),
  #    binning=[10,-pi,pi],
  #))

  #plots.append(Plot( name = "met_phi_raw_corr",
  #    texX = 'raw #phi(E_{T}^{miss})', texY = 'Number of Events / 20 GeV',
  #    attribute = lambda event, sample: event.RawMET_phi_corr,
  #    binning=[10,-pi,pi],
  #))

  if not args.blinded:
    plots.append(Plot(
      texX = 'M_{T2}(ll) (GeV)', texY = 'Number of Events / 20 GeV',
      attribute = TreeVariable.fromString( "dl_mt2ll/F" ),
      binning=[300/20,0,300],
    ))

    plots.append(Plot( name = "dl_mt2ll_coarse",
      texX = 'M_{T2}(ll) (GeV)', texY = 'Number of Events / 20 GeV',
      attribute = TreeVariable.fromString( "dl_mt2ll/F" ),
      binning=Binning.fromThresholds([0,20,40,60,80,100,140,240,340]),
    ))
    #plots.append(Plot( name = "dl_mt2ll_raw",
    #  texX = 'M_{T2}(ll) (GeV)', texY = 'Number of Events / 20 GeV',
    #  attribute = lambda event, sample: event.dl_mt2ll_raw,
    #  binning=[300/20,0,300],
    #))

  #plots.append(Plot( name = "met_pt_corr",
  #    texX = 'corr E_{T}^{miss} (GeV)', texY = 'Number of Events / 20 GeV',
  #    attribute = lambda event, sample: event.met_pt_corr,
  #    binning=[400/20,0,400],
  #))

  #plots.append(Plot( name = "met_pt_raw_corr",
  #    texX = 'corr E_{T}^{miss} (GeV)', texY = 'Number of Events / 20 GeV',
  #    attribute = lambda event, sample: event.RawMET_pt_corr,
  #    binning=[400/20,0,400],
  #))
    
  #plots.append(Plot( name = "dl_mt2ll_corr",
  #    texX = 'corr M_{T2}(ll) (GeV)', texY = 'Number of Events / 20 GeV',
  #    attribute = lambda event, sample: event.dl_mt2ll_corr,
  #    binning=[300/20,0,300],
  #))

  #plots.append(Plot( name = "dl_mt2ll_raw_corr",
  #    texX = 'raw corr M_{T2}(ll) (GeV)', texY = 'Number of Events / 20 GeV',
  #    attribute = lambda event, sample: event.dl_mt2ll_raw_corr,
  #    binning=[300/20,0,300],
  #))

  plots.append(Plot( name = "qT",
    texX = 'q_{T} (GeV)', texY = 'Number of Events / 50 GeV',
    attribute = lambda event, sample: sqrt((event.l1_pt*cos(event.l1_phi) + event.l2_pt*cos(event.l2_phi) + event.met_pt*cos(event.met_phi))**2 + (event.l1_pt*sin(event.l1_phi) + event.l2_pt*sin(event.l2_phi) + event.met_pt*sin(event.met_phi))**2),
    binning= [1000/50,0,1000]),
  )

  plots.append(Plot(
    texX = 'number of jets', texY = 'Number of Events',
    attribute = TreeVariable.fromString('nJetGood/I'),
    binning=[14,0,14],
  ))
  if "2017" in args.era:  
    plots.append(Plot( name = "nJet_EE",
      texX = 'number of jets', texY = 'Number of Events',
      attribute = lambda event, sample: event.nJet_EE,
      binning=[14,0,14],
    ))
    plots.append(Plot( name = "nJet_EE_pt30To50",
      texX = 'number of jets', texY = 'Number of Events',
      attribute = lambda event, sample: event.nJet_EE_pt30To50,
      binning=[14,0,14],
    ))
    plots.append(Plot( name = "nJet_EE_ptTo50",
      texX = 'number of jets', texY = 'Number of Events',
      attribute = lambda event, sample: event.nJet_EE_ptTo50,
      binning=[14,0,14],
    ))
    plots.append(Plot( name = "nJet_EE_ptTo40",
      texX = 'number of jets', texY = 'Number of Events',
      attribute = lambda event, sample: event.nJet_EE_ptTo40,
      binning=[14,0,14],
    ))
    plots.append(Plot( name = "nJet_EE_ptTo30",
      texX = 'number of jets', texY = 'Number of Events',
      attribute = lambda event, sample: event.nJet_EE_ptTo30,
      binning=[14,0,14],
    ))
    plots.append(Plot( name = "nJet_EE_pt50",
      texX = 'number of jets', texY = 'Number of Events',
      attribute = lambda event, sample: event.nJet_EE_pt50,
      binning=[14,0,14],
    ))
    plots.append(Plot( name = "nJet_EE_pt40",
      texX = 'number of jets', texY = 'Number of Events',
      attribute = lambda event, sample: event.nJet_EE_pt40,
      binning=[14,0,14],
    ))
    plots.append(Plot( name = "nJet_EE_pt30",
      texX = 'number of jets', texY = 'Number of Events',
      attribute = lambda event, sample: event.nJet_EE_pt30,
      binning=[14,0,14],
    ))
    plots.append(Plot( name = "badJetE",
      texX = 'badEEJetEnergy', texY = 'Number of Events',
      attribute = lambda event, sample: event.badJetE,
      binning=[40,0,400],
    ))
    plots.append(Plot( name = "badJetPt",
      texX = 'badEEJetPt', texY = 'Number of Events',
      attribute = lambda event, sample: event.badJetPt,
      binning=[40,0,200],
    ))

  plots.append(Plot(
    texX = 'number of medium b-tags (CSVM)', texY = 'Number of Events',
    attribute = TreeVariable.fromString('nBTag/I'),
    binning=[8,0,8],
  ))

  #plots.append(Plot(
  #  texX = 'H_{T} (GeV)', texY = 'Number of Events / 25 GeV',
  #  attribute = TreeVariable.fromString( "ht/F" ),
  #  binning=[500/25,0,600],
  #))

  plots.append(Plot(
    texX = 'm(ll) of leading dilepton (GeV)', texY = 'Number of Events / 4 GeV',
    attribute = TreeVariable.fromString( "dl_mass/F" ),
    binning=[200/4,0,200],
  ))

  plots.append(Plot(
    texX = 'p_{T}(ll) (GeV)', texY = 'Number of Events / 10 GeV',
    attribute = TreeVariable.fromString( "dl_pt/F" ),
    binning=[20,0,400],
  ))

  plots.append(Plot(
      texX = '#eta(ll) ', texY = 'Number of Events',
      name = 'dl_eta', attribute = lambda event, sample: abs(event.dl_eta), read_variables = ['dl_eta/F'],
      binning=[10,0,3],
  ))

  plots.append(Plot(
    texX = '#phi(ll)', texY = 'Number of Events',
    attribute = TreeVariable.fromString( "dl_phi/F" ),
    binning=[10,-pi,pi],
  ))

  plots.append(Plot(
    texX = 'Cos(#Delta#phi(ll, E_{T}^{miss}))', texY = 'Number of Events',
    name = 'cosZMetphi',
    attribute = lambda event, sample: cos( event.dl_phi - event.met_phi ), 
    read_variables = ["met_phi/F", "dl_phi/F"],
    binning = [10,-1,1],
  ))

  plots.append(Plot(
    texX = 'p_{T}(l_{1}) (GeV)', texY = 'Number of Events / 15 GeV',
    attribute = TreeVariable.fromString( "l1_pt/F" ),
    binning=[20,0,300],
  ))

  plots.append(Plot(
    texX = '#eta(l_{1})', texY = 'Number of Events',
    name = 'l1_eta', attribute = lambda event, sample: abs(event.l1_eta), read_variables = ['l1_eta/F'],
    binning=[15,0,3],
  ))

  plots.append(Plot(
    texX = '#phi(l_{1})', texY = 'Number of Events',
    attribute = TreeVariable.fromString( "l1_phi/F" ),
    binning=[10,-pi,pi],
  ))

  plots.append(Plot(
    texX = 'I_{mini}(l_{1})', texY = 'Number of Events',
    name = 'l1_miniRelIso', attribute = lambda event, sample: event.l1_miniRelIso, read_variables = ['l1_miniRelIso/F'],
    binning=[20,0,.5],
  ))

  plots.append(Plot(
    texX = 'pdgId(l1)', texY = 'Number of Events',
    attribute = TreeVariable.fromString( "l1_pdgId/I" ),
    binning=[30,-15,15],
  ))

  plots.append(Plot(
    texX = 'p_{T}(l_{2}) (GeV)', texY = 'Number of Events / 15 GeV',
    attribute = TreeVariable.fromString( "l2_pt/F" ),
    binning=[20,0,300],
  ))

  plots.append(Plot(
    texX = '#eta(l_{2})', texY = 'Number of Events',
    name = 'l2_eta', attribute = lambda event, sample: abs(event.l2_eta), read_variables = ['l2_eta/F'],
    binning=[15,0,3],
  ))

  plots.append(Plot(
    texX = '#phi(l_{2})', texY = 'Number of Events',
    attribute = TreeVariable.fromString( "l2_phi/F" ),
    binning=[10,-pi,pi],
  ))
  
  plots.append(Plot(
    texX = 'I_{mini}(l_{2})', texY = 'Number of Events',
    name = 'l2_miniRelIso', attribute = lambda event, sample: event.l2_miniRelIso, read_variables = ['l2_miniRelIso/F'],
    binning=[20,0,.5],
  ))

  plots.append(Plot(
    texX = 'pdgId(l2)', texY = 'Number of Events',
    attribute = TreeVariable.fromString( "l2_pdgId/I" ),
    binning=[30,-15,15],
  ))

  # plot trailing lepton quantities
  plots2D.append(Plot2D(
    name = 'eff_num_trailing_ele_sip3Dlt4',
    stack = stack,
    attribute = (
      lambda event, sample: event.l2_eta if event.l2_eleIndex>0 and event.Electron_sip3d[event.l2_eleIndex]<4 else float('nan'),
      lambda event, sample: event.l2_pt  if event.l2_eleIndex>0 and event.Electron_sip3d[event.l2_eleIndex]<4 else float('nan'),
    ),
    texX = 'trailing lepton #eta', texY = 'trailing lepton p_{T} (GeV)',
    binning=[Binning.fromThresholds([-2.5, -2, -1.566, -1.444, -0.8, 0, 0.8, 1.444, 1.566, 2.0, 2.5]), Binning.fromThresholds([20, 35, 50, 100, 200])],
  ))
  plots2D.append(Plot2D(
    name = 'eff_den_trailing_ele_sip3Dlt4', #keep the cut in the name, but this is the denominator, so remove it from the lambda function
    stack = stack,
    attribute = (
      lambda event, sample: event.l2_eta if event.l2_eleIndex>0 else float('nan'),
      lambda event, sample: event.l2_pt  if event.l2_eleIndex>0 else float('nan'),
    ),
    texX = 'trailing lepton #eta', texY = 'trailing lepton p_{T} (GeV)',
    binning=[Binning.fromThresholds([-2.5, -2, -1.566, -1.444, -0.8, 0, 0.8, 1.444, 1.566, 2.0, 2.5]), Binning.fromThresholds([20, 35, 50, 100, 200])],
  ))

  # Plots only when at least one jet:
  if args.selection.count('njet2') or args.selection.count('njet1') or args.selection.count('njet01'):

    plots2D.append(Plot2D(
      name = 'leading_jet_occ',
      stack = stack,
      attribute = (
        lambda event, sample: event.JetGood_eta[0],
        lambda event, sample: event.JetGood_phi[0],
      ),
      texX = '#eta(leading jet) (GeV)', texY = '#phi(leading jet) (GeV)',
      binning=[16, -3.0, 3.0, 10, -pi, pi],
    ))

    plots.append(Plot(
      texX = 'p_{T}(leading jet) (GeV)', texY = 'Number of Events / 30 GeV',
      name = 'jet1_pt', attribute = lambda event, sample: event.JetGood_pt[0],
      binning=[600/30,0,600],
    ))

    plots.append(Plot(
      texX = '#eta(leading jet) (GeV)', texY = 'Number of Events',
      name = 'jet1_eta', attribute = lambda event, sample: event.JetGood_eta[0],
      binning=[20,-3,3],
    ))

    plots.append(Plot(
      texX = '#phi(leading jet) (GeV)', texY = 'Number of Events',
      name = 'jet1_phi', attribute = lambda event, sample: event.JetGood_phi[0],
      binning=[20,-pi,pi],
    ))

    plots.append(Plot(
      name = 'cosMetJet1phi',
      texX = 'Cos(#Delta#phi(E_{T}^{miss}, leading jet))', texY = 'Number of Events',
      attribute = lambda event, sample: cos( event.met_phi - event.JetGood_phi[0]), 
      read_variables = ["met_phi/F", "JetGood[phi/F]"],
      binning = [10,-1,1],
    ))
    
    plots.append(Plot(
      name = 'cosMetJet1phi_fine',
      texX = 'Cos(#Delta#phi(E_{T}^{miss}, leading jet))', texY = 'Number of Events',
      attribute = lambda event, sample: cos( event.met_phi - event.JetGood_phi[0] ) , 
      read_variables = ["met_phi/F", "JetGood[phi/F]"],
      binning = [100,-1,1],
    ))

    plots.append(Plot(
      name = 'cosZJet1phi',
      texX = 'Cos(#Delta#phi(Z, leading jet))', texY = 'Number of Events',
      attribute = lambda event, sample: cos( event.dl_phi - event.JetGood_phi[0] ),
      read_variables =  ["dl_phi/F", "JetGood[phi/F]"],
      binning = [10,-1,1],
    ))
    plots.append(Plot( name = "u_para", 
      texX = "u_{#parallel} (GeV)", texY = 'Number of Events / 30 GeV',
      attribute = lambda event, sample: - event.met_pt*cos(event.met_phi-event.dl_phi),
      binning=[80, -200,200],
    ))
    plots.append(Plot( name = "u_perp", 
      texX = "u_{#perp} (GeV)", texY = 'Number of Events / 30 GeV',
      attribute = lambda event, sample: - event.met_pt*cos(event.met_phi-(event.dl_phi-pi/2)),
      binning=[80, -200,200],
    ))
    #plots.append(Plot( name = "u_para_raw", 
    #  texX = "u_{#parallel} (GeV)", texY = 'Number of Events / 30 GeV',
    #  attribute = lambda event, sample: - event.RawMET_pt*cos(event.RawMET_phi-event.dl_phi),
    #  binning=[80, -200,200],
    #))
    #plots.append(Plot( name = "u_perp_raw", 
    #  texX = "u_{#perp} (GeV)", texY = 'Number of Events / 30 GeV',
    #  attribute = lambda event, sample: - event.RawMET_pt*cos(event.RawMET_phi-(event.dl_phi-pi/2)),
    #  binning=[80, -200,200],
    #))
    #plots.append(Plot( name = "u_para_corr", 
    #  texX = "u_{#parallel} corr. (GeV)", texY = 'Number of Events / 30 GeV',
    #  attribute = lambda event, sample: - event.met_pt_corr*cos(event.met_phi_corr-event.dl_phi),
    #  binning=[80, -200,200],
    #))
    #plots.append(Plot( name = "u_perp_corr", 
    #  texX = "u_{#perp} corr. (GeV)", texY = 'Number of Events / 30 GeV',
    #  attribute = lambda event, sample: - event.met_pt_corr*cos(event.met_phi_corr-(event.dl_phi-pi/2)),
    #  binning=[80, -200,200],
    #))
    #plots.append(Plot( name = "u_para_raw_corr", 
    #  texX = "u_{#parallel} corr. (GeV)", texY = 'Number of Events / 30 GeV',
    #  attribute = lambda event, sample: - event.RawMET_pt_corr*cos(event.RawMET_phi_corr-event.dl_phi),
    #  binning=[80, -200,200],
    #))
    #plots.append(Plot( name = "u_perp_raw_corr", 
    #  texX = "u_{#perp} corr. (GeV)", texY = 'Number of Events / 30 GeV',
    #  attribute = lambda event, sample: - event.RawMET_pt_corr*cos(event.RawMET_phi_corr-(event.dl_phi-pi/2)),
    #  binning=[80, -200,200],
    #))

    if args.plotUPara:
        # u_para u_perp closure plots
        u_para_binning   =  [ i*20 for i in range(-10, 11) ]
        qt_binning    = [0, 50, 100, 150, 200, 300 ]
        qt_bins = [ (qt_binning[i],qt_binning[i+1]) for i in range(len(qt_binning)-1) ]
        var_binning   = [ pi*(i-5)/5. for i in range(0,11) ]
        var_bins      = [ (var_binning[i],var_binning[i+1]) for i in range(len(var_binning)-1) ]
        #for var_bin in var_bins:
        #    for qt_bin in qt_bins:
        #        postfix = "phill_%3.2f_%3.2f_qt_%i_%i"%( var_bin[0], var_bin[1], qt_bin[0], qt_bin[1] )
        #        plots.append(Plot( name = "u_para_" + postfix, 
        #          texX = "u_{#parallel} (GeV)", texY = 'Number of Events / 30 GeV',
        #          attribute = lambda event, sample: - event.met_pt*cos(event.met_phi-event.dl_phi),
        #          weight = recoil_weight(var_bin, qt_bin),
        #          binning=[80, -200,200],
        #        ))
        #        plots.append(Plot( name = "u_perp_" + postfix, 
        #          texX = "u_{#perp} (GeV)", texY = 'Number of Events / 30 GeV',
        #          attribute = lambda event, sample: - event.met_pt*cos(event.met_phi-(event.dl_phi-pi/2)),
        #          weight = recoil_weight(var_bin, qt_bin),
        #          binning=[80, -200,200],
        #        ))
        #        plots.append(Plot( name = "u_para_corr_" + postfix, 
        #          texX = "u_{#parallel} corr. (GeV)", texY = 'Number of Events / 30 GeV',
        #          attribute = lambda event, sample: - event.met_pt_corr*cos(event.met_phi_corr-event.dl_phi),
        #          weight = recoil_weight(var_bin, qt_bin),
        #          binning=[80, -200,200],
        #        ))
        #        plots.append(Plot( name = "u_perp_corr_" + postfix, 
        #          texX = "u_{#perp} corr. (GeV)", texY = 'Number of Events / 30 GeV',
        #          attribute = lambda event, sample: - event.met_pt_corr*cos(event.met_phi_corr-(event.dl_phi-pi/2)),
        #          weight = recoil_weight(var_bin, qt_bin),
        #          binning=[80, -200,200],
        #        ))

#  # Plots only when at least two jets:
  if args.selection.count('njet2'):
    plots.append(Plot(
      texX = 'p_{T}(2nd leading jet) (GeV)', texY = 'Number of Events / 30 GeV',
      name = 'jet2_pt', attribute = lambda event, sample: event.JetGood_pt[1],
      binning=[600/30,0,600],
    ))

    plots.append(Plot(
      texX = '#eta(2nd leading jet) (GeV)', texY = 'Number of Events',
      name = 'jet2_eta', attribute = lambda event, sample: event.JetGood_eta[1],
      binning=[20,-3,3],
    ))

    plots.append(Plot(
      texX = '#phi(2nd leading jet) (GeV)', texY = 'Number of Events',
      name = 'jet2_phi', attribute = lambda event, sample: event.JetGood_phi[1],
      binning=[20,-pi,pi],
    ))

    plots.append(Plot(
      name = 'cosMetJet2phi',
      texX = 'Cos(#Delta#phi(E_{T}^{miss}, second jet))', texY = 'Number of Events',
      attribute = lambda event, sample: cos( event.met_phi - event.JetGood_phi[1] ) , 
      read_variables = ["met_phi/F", "JetGood[phi/F]"],
      binning = [10,-1,1],
    ))
    
    plots.append(Plot(
      name = 'cosMetJet2phi_fine',
      texX = 'Cos(#Delta#phi(E_{T}^{miss}, second jet))', texY = 'Number of Events',
      attribute = lambda event, sample: cos( event.met_phi - event.JetGood_phi[1] ) , 
      read_variables = ["met_phi/F", "JetGood[phi/F]"],
      binning = [20,-1,1],
    ))

    plots.append(Plot(
      name = 'cosZJet2phi',
      texX = 'Cos(#Delta#phi(Z, 2nd leading jet))', texY = 'Number of Events',
      attribute = lambda event, sample: cos( event.dl_phi - event.JetGood_phi[0] ),
      read_variables = ["dl_phi/F", "JetGood[phi/F]"],
      binning = [10,-1,1],
    ))

    plots.append(Plot(
      name = 'cosJet1Jet2phi',
      texX = 'Cos(#Delta#phi(leading jet, 2nd leading jet))', texY = 'Number of Events',
      attribute = lambda event, sample: cos( event.JetGood_phi[1] - event.JetGood_phi[0] ) ,
      read_variables =  ["JetGood[phi/F]"],
      binning = [10,-1,1],
    ))

    plots.append(Plot(
      texX = 'M_{T2}(bb) (GeV)', texY = 'Number of Events / 30 GeV',
      attribute = TreeVariable.fromString( "dl_mt2bb/F" ),
      binning=[420/30,70,470],
    ))

    plots.append(Plot(
      texX = 'M_{T2}(blbl) (GeV)', texY = 'Number of Events / 30 GeV',
      attribute = TreeVariable.fromString( "dl_mt2blbl/F" ),
      binning=[420/30,0,400],

    ))

    plots.append(Plot( name = "dl_mt2blbl_coarse",       # SR binning of MT2ll
      texX = 'M_{T2}(blbl) (GeV)', texY = 'Number of Events / 30 GeV',
      attribute = TreeVariable.fromString( "dl_mt2blbl/F" ),
      binning=Binning.fromThresholds([0,20,40,60,80,100,120,140,160,200,250,300,350]),
    ))

    plots.append(Plot( name = "dl_mt2blbl_bj_dR",       # SR binning of MT2ll
      texX = 'dR(bj0,bj1))', texY = 'Number of Events / 30 GeV',
      attribute = lambda event, sample: event.mt2blbl_bj_dR,
      binning=[12,0,6],
    ))

    plots.append(Plot( name = "dl_mt2blbl_bj_dPhi",       # SR binning of MT2ll
      texX = 'dPhi(bj0,bj1))', texY = 'Number of Events / 30 GeV',
      attribute = lambda event, sample: event.mt2blbl_bj_dPhi,
      binning=[12,0,pi],
    ))

    plots.append(Plot( name = "dl_mt2blbl_bj_dEta",       # SR binning of MT2ll
      texX = 'dEta(bj0,bj1))', texY = 'Number of Events / 30 GeV',
      attribute = lambda event, sample: event.mt2blbl_bj_dEta,
      binning=[12,-6,6],
    ))

    plots.append(Plot( name = "dl_mt2blbl_bj_mass",       # SR binning of MT2ll
      texX = 'mass(bj0,bj1))', texY = 'Number of Events / 30 GeV',
      attribute = lambda event, sample: event.mt2blbl_bj_mass,
      binning=[12,0,600],
    ))

    plots.append(Plot( name = "dl_mt2blbl_bj0_pt",       # SR binning of MT2ll
      texX = 'bj0(pt)', texY = 'Number of Events / 30 GeV',
      attribute = lambda event, sample: event.bj0['pt'] if event.bj0 is not None else float('nan'),
      binning=[12,0,300],
    ))
    plots.append(Plot( name = "dl_mt2blbl_bj1_pt",       # SR binning of MT2ll
      texX = 'bj1(pt)', texY = 'Number of Events / 30 GeV',
      attribute = lambda event, sample: event.bj1['pt'] if event.bj1 is not None else float('nan'),
      binning=[12,0,300],
    ))
    plots.append(Plot( name = "dl_mt2blbl_bj0_phi",       # SR binning of MT2ll
      texX = 'bj0(phi)', texY = 'Number of Events / 30 GeV',
      attribute = lambda event, sample: event.bj0['phi'] if event.bj0 is not None else float('nan'),
      binning=[12,-pi,pi],
    ))
    plots.append(Plot( name = "dl_mt2blbl_bj1_phi",       # SR binning of MT2ll
      texX = 'bj1(phi)', texY = 'Number of Events / 30 GeV',
      attribute = lambda event, sample: event.bj1['phi'] if event.bj1 is not None else float('nan'),
      binning=[12,-pi,pi],
    ))
    plots.append(Plot( name = "dl_mt2blbl_bj0_eta",       # SR binning of MT2ll
      texX = 'bj0(eta)', texY = 'Number of Events / 30 GeV',
      attribute = lambda event, sample: event.bj0['eta'] if event.bj0 is not None else float('nan'),
      binning=[12,-3,3],
    ))
    plots.append(Plot( name = "dl_mt2blbl_bj1_eta",       # SR binning of MT2ll
      texX = 'bj1(eta)', texY = 'Number of Events / 30 GeV',
      attribute = lambda event, sample: event.bj1['eta'] if event.bj1 is not None else float('nan'),
      binning=[12,-3,3],
    ))
    plots.append(Plot( name = "dl_mt2blbl_bj0_DeepCSV",       # SR binning of MT2ll
      texX = 'bj0(DeepCSV)', texY = 'Number of Events / 30 GeV',
      attribute = lambda event, sample: event.bj0['btagDeepB'] if event.bj0 is not None else float('nan'),
      binning=[10,0,1],
    ))
    plots.append(Plot( name = "dl_mt2blbl_bj1_DeepCSV",       # SR binning of MT2ll
      texX = 'bj1(DeepCSV)', texY = 'Number of Events / 30 GeV',
      attribute = lambda event, sample: event.bj1['btagDeepB'] if event.bj1 is not None else float('nan'),
      binning=[10,0,1],
    ))

   
  plotting.fill(plots + plots2D, read_variables = read_variables, sequence = sequence)

  # Get normalization yields from yield histogram
  for plot in plots:
    if plot.name == "yield":
      for i, l in enumerate(plot.histos):
        for j, h in enumerate(l):
          yields[mode][plot.stack[i][j].name] = h.GetBinContent(h.FindBin(0.5+index))
          h.GetXaxis().SetBinLabel(1, "#mu#mu")
          h.GetXaxis().SetBinLabel(2, "e#mu")
          h.GetXaxis().SetBinLabel(3, "ee")
  if args.noData: yields[mode]["data"] = 0

  yields[mode]["MC"] = sum(yields[mode][s.name] for s in mc_)
  dataMCScale        = yields[mode]["data"]/yields[mode]["MC"] if yields[mode]["MC"] != 0 else float('nan')

  drawPlots(plots + plots2D, mode, dataMCScale)
  allPlots[mode] = plots + plots2D

# Add the different channels into SF and all
for mode in ["SF","all"]:
  yields[mode] = {}
  for y in yields[allModes[0]]:
    try:    yields[mode][y] = sum(yields[c][y] for c in (['ee','mumu'] if mode=="SF" else ['ee','mumu','mue']))
    except: yields[mode][y] = 0
  dataMCScale = yields[mode]["data"]/yields[mode]["MC"] if yields[mode]["MC"] != 0 else float('nan')

  for plot in allPlots['mumu']:
    for plot2 in (p for p in (allPlots['ee'] if mode=="SF" else allPlots["mue"]) if p.name == plot.name):  #For SF add EE, second round add EMu for all
      for i, j in enumerate(list(itertools.chain.from_iterable(plot.histos))):
	for k, l in enumerate(list(itertools.chain.from_iterable(plot2.histos))):
	  if i==k:
	    j.Add(l)

  drawPlots(allPlots['mumu'], mode, dataMCScale)


logger.info( "Done with prefix %s and selectionString %s", args.selection, cutInterpreter.cutString(args.selection) )

