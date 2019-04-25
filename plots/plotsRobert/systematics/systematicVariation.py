#!/usr/bin/env python
''' analysis script for standard plots with systematic errors
'''

# Standard imports and batch mode
import ROOT
ROOT.gROOT.SetBatch(True)
import operator
import pickle, os, time
from math                                import sqrt, cos, sin, pi, atan2

# RootTools
from RootTools.core.standard             import *

#Analysis / StopsDilepton / Samples
from StopsDilepton.tools.user            import plot_directory
from StopsDilepton.tools.helpers         import deltaPhi
from Samples.Tools.metFilters            import getFilterCut
from StopsDilepton.tools.cutInterpreter  import cutInterpreter
from StopsDilepton.tools.RecoilCorrector import RecoilCorrector
from StopsDilepton.tools.mt2Calculator   import mt2Calculator
from Analysis.Tools.puReweighting        import getReweightingFunction
from Analysis.Tools.DirDB                import DirDB


# Arguments
import argparse
argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--logLevel',          action='store',      default='INFO',     nargs='?', choices=['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'TRACE', 'NOTSET'], help="Log level for logging")
argParser.add_argument('--signal',            action='store',      default=None,        nargs='?', choices=['None', "T2tt",'DM'], help="Add signal to plot")
argParser.add_argument('--noData',            action='store_true', default=False,       help='also plot data?')
argParser.add_argument('--plot_directory',    action='store',      default='v1')
argParser.add_argument('--selection',         action='store',            default='njet2p-btag1p-relIso0.12-looseLeptonVeto-mll20-met80-metSig5-dPhiJet0-dPhiJet1')
#argParser.add_argument('--normalizationSelection',  action='store',      default='njet2p-btag1p-relIso0.12-looseLeptonVeto-mll20-met80-metSig5-dPhiJet0-dPhiJet1-mt2llTo100')
argParser.add_argument('--variation',         action='store',      default='central', help="Which systematic variation to run.")
argParser.add_argument('--showOnly',          action='store',      default=None)
argParser.add_argument('--small',             action='store_true',     help='Run only on a small subset of the data?')
argParser.add_argument('--mode',              action='store',       choices = ['mumu','ee','mue'],   help='Which mode?')
argParser.add_argument('--normalizeBinWidth', action='store_true', default=False,       help='normalize wider bins?')
argParser.add_argument('--reweightPU',         action='store', default='Central', choices=[ 'Central', 'VUp'] )
#argParser.add_argument('--recoil',             action='store', type=str,      default="Central", choices = ["nvtx", "VUp", "Central"])
argParser.add_argument('--era',                action='store', type=str,      default="Run2016")

args = argParser.parse_args()

# Logger
import StopsDilepton.tools.logger as logger
import RootTools.core.logger as logger_rt
import Analysis.Tools.logger as logger_an
logger    = logger.get_logger(   args.logLevel, logFile = None)
logger_rt = logger_rt.get_logger(args.logLevel, logFile = None)
logger_an = logger_an.get_logger(args.logLevel, logFile = None)

# Year
if "2016" in args.era:
    year = 2016
elif "2017" in args.era:
    year = 2017
elif "2018" in args.era:
    year = 2018

logger.info( "Working in year %i", year )

if args.reweightPU == 'Central':
    nominalPuWeight, upPUWeight, downPUWeight = "reweightPU", "reweightPUUp", "reweightPUDown"
elif args.reweightPU == 'VUp':
    nominalPuWeight, upPUWeight, downPUWeight = "reweightPUVUp", "reweightPUVVUp", "reweightPUUp"

def jetSelectionModifier( sys, returntype = "func"):
    #Need to make sure all jet variations of the following observables are in the ntuple
    variiedJetObservables = ['nJetGood', 'nBTag', 'dl_mt2ll', 'dl_mt2blbl', 'MET_significance', 'met_pt', 'metSig']
    if returntype == "func":
        def changeCut_( string ):
            for s in variiedJetObservables:
                string = string.replace(s, s+'_'+sys)
            return string
        return changeCut_
    elif returntype == "list":
        return [ v+'_'+sys for v in variiedJetObservables ]

def metSelectionModifier( sys ):
    #Need to make sure all MET variations of the following observables are in the ntuple
    variiedMetObservables = ['dl_mt2ll', 'dl_mt2blbl', 'MET_significance', 'met_pt', 'metSig']
    if returntype == "func":
        def changeCut_( string ):
            for s in variiedMetObservables:
                string = string.replace(s, s+'_'+sys)
            return string
        return changeCut_
    elif returntype == "list":
        return [ v+'_'+sys for v in variiedMetObservables ]

def metSelectionModifier( sys ):
    def changeCut_( string ):
        return string.replace('met_pt', 'met_pt_'+sys).replace('metSig', 'metSig_'+sys).replace('dl_mt2ll', 'dl_mt2ll_'+sys).replace('dl_mt2bb', 'dl_mt2bb_'+sys).replace('dl_mt2blbl', 'dl_mt2blbl_'+sys).replace('MET_significance', 'MET_significance_'+sys)
    return changeCut_

# these are the nominal MC weights we always apply
nominalMCWeights = ["weight", "reweightLeptonSF", "reweightPU", "reweightDilepTrigger", "reweightBTag_SF"]

# weight the MC according to a variation
def MC_WEIGHT( variation, returntype = "string"):
    variiedMCWeights = list(nominalMCWeights)   # deep copy
    if variation.has_key('replaceWeight'):
        for i_w, w in enumerate(variiedMCWeights):
            if w == variation['replaceWeight'][0]:
                variiedMCWeights[i_w] = variation['replaceWeight'][1]
                break
        # Let's make sure we don't screw it up ... because we mostly do.
        if variiedMCWeights==nominalMCWeights:
            raise RuntimeError( "Tried to change weight %s to %s but didn't find it in list %r" % ( variation['replaceWeight'][0], variation['replaceWeight'][1], variiedMCWeights ))
    # multiply strings for ROOT weights
    if returntype == "string":
        return "*".join(variiedMCWeights)
    # create a function that multiplies the attributes of the event
    elif returntype == "func":
        getters = map( operator.attrgetter, variiedMCWeights)
        def weight_( event, sample):
            return reduce(operator.mul, [g(event) for g in getters], 1)
        return weight_
    elif returntype == "list":
        return variiedMCWeights

def data_weight( event, sample ):
    return event.weight

data_weight_string = "weight"

# Define all systematic variations
variations = {
    'central'           : {},
    'jesTotalUp'        : {'selectionModifier':jetSelectionModifier('jesTotalUp'),      'read_variables' : [ '%s/F'%v for v in nominalMCWeights + jetSelectionModifier('jesTotalUp','list')]},
    'jesTotalDown'      : {'selectionModifier':jetSelectionModifier('jesTotalDown'),    'read_variables' : [ '%s/F'%v for v in nominalMCWeights + jetSelectionModifier('jesTotalDown','list')]},
    'unclustEnUp'       : {'selectionModifier':metSelectionModifier('unclustEnUp'),     'read_variables' : [ '%s/F'%v for v in nominalMCWeights + jetSelectionModifier('unclustEnUp','list')]},
    'unclustEnDown'     : {'selectionModifier':metSelectionModifier('unclustEnDown'),   'read_variables' : [ '%s/F'%v for v in nominalMCWeights + jetSelectionModifier('unclustEnDown','list')]},
    'PUUp'              : {'replaceWeight':(nominalPuWeight,upPUWeight),                'read_variables' : [ '%s/F'%v for v in nominalMCWeights + [upPUWeight] ]},
    'PUDown'            : {'replaceWeight':(nominalPuWeight,downPUWeight),              'read_variables' : [ '%s/F'%v for v in nominalMCWeights + [downPUWeight] ]},
    'BTag_SF_b_Down'    : {'replaceWeight':('reweightBTag_SF','reweight_BTag_SF_b_Down'),        'read_variables' : [ '%s/F'%v for v in nominalMCWeights + ['reweight_BTag_SF_b_Down']]},  
    'BTag_SF_b_Up'      : {'replaceWeight':('reweightBTag_SF','reweight_BTag_SF_b_Up'),          'read_variables' : [ '%s/F'%v for v in nominalMCWeights + ['reweight_BTag_SF_b_Up'] ]},
    'BTag_SF_l_Down'    : {'replaceWeight':('reweightBTag_SF','reweight_BTag_SF_l_Down'),        'read_variables' : [ '%s/F'%v for v in nominalMCWeights + ['reweight_BTag_SF_l_Down']]},
    'BTag_SF_l_Up'      : {'replaceWeight':('reweightBTag_SF','reweight_BTag_SF_l_Up'),          'read_variables' : [ '%s/F'%v for v in nominalMCWeights + ['reweight_BTag_SF_l_Up'] ]},
    'DilepTriggerDown'  : {'replaceWeight':('reweightDilepTrigger','reweightDilepTriggerDown'),  'read_variables' : [ '%s/F'%v for v in nominalMCWeights + ['reweightDilepTriggerDown']]},
    'DilepTriggerUp'    : {'replaceWeight':('reweightDilepTrigger','reweightDilepTriggerUp'),    'read_variables' : [ '%s/F'%v for v in nominalMCWeights + ['reweightDilepTriggerUp']]},
    'LeptonSFDown'      : {'replaceWeight':('reweightLeptonSF','reweightLeptonSFDown'),          'read_variables' : [ '%s/F'%v for v in nominalMCWeights + ['reweightLeptonSFDown']]},
    'LeptonSFUp'        : {'replaceWeight':('reweightLeptonSF','reweightLeptonSFUp'),            'read_variables' : [ '%s/F'%v for v in nominalMCWeights + ['reweightLeptonSFUp']]},
#    'TopPt':{},
#   'JERUp':{},
#   'JERDown':{},
}

# Add a default selection modifier that does nothing
for key, variation in variations.iteritems():
    if not variation.has_key('selectionModifier'):
        variation['selectionModifier'] = lambda string:string
    if not variation.has_key('read_variables'):
        variation['read_variables'] = [] 

# Check if we know the variation
if args.variation not in variations.keys():
    raise RuntimeError( "Variation %s not among the known: %s", args.variation, ",".join( variation.keys() ) )

# Systematic pairs:( 'name', 'up', 'down' )
systematics = [\
    ('JEC',         'jesTotalUp', 'jesTotalDown'),
    ('Unclustered', 'unclustEnUp', 'unclustEnDown'), 
    ('PU',          'PUUp', 'PUDown'),
    ('BTag_b',      'BTag_SF_b_Down', 'BTag_SF_b_Up' ),
    ('BTag_l',      'BTag_SF_l_Down', 'BTag_SF_l_Up'),
    ('trigger',     'DilepTriggerDown', 'DilepTriggerUp'),
    ('leptonSF',    'LeptonSFDown', 'LeptonSFUp'),
    # ('TopPt',       'TopPt', 'central'),
    # ('JER',         'JERUp', 'JERDown'),
]

#if args.selectSys != "all" and args.selectSys != "combine": all_systematics = [args.selectSys if args.selectSys != 'None' else None]
#else:                                                       all_systematics = [None] + weight_systematics + jme_systematics
#
if args.noData:                   args.plot_directory += "_noData"
if args.signal == "DM":           args.plot_directory += "_DM"
if args.signal == "T2tt":         args.plot_directory += "_T2tt"
if args.small:                    args.plot_directory += "_small"
if args.reweightPU:               args.plot_directory += "_reweightPU%s"%args.reweightPU
#if args.recoil:                   args.plot_directory += '_recoil_'+args.recoil
    
if year == 2016:
    from StopsDilepton.samples.nanoTuples_Summer16_postProcessed import *
    from StopsDilepton.samples.nanoTuples_Run2016_17Jul2018_postProcessed import *
    Top_pow, TTXNoZ, TTZ_LO, multiBoson, DY_HT_LO = Top_pow_16, TTXNoZ_16, TTZ_16, multiBoson_16, DY_HT_LO_16
    mc              = [ Top_pow_16, TTXNoZ_16, TTZ_16, multiBoson_16, DY_HT_LO_16]
elif year == 2017:
    from StopsDilepton.samples.nanoTuples_Fall17_postProcessed import *
    from StopsDilepton.samples.nanoTuples_Run2017_31Mar2018_postProcessed import *
    Top_pow, TTXNoZ, TTZ_LO, multiBoson, DY_HT_LO = Top_pow_17, TTXNoZ_17, TTZ_17, multiBoson_17, DY_LO_17
    mc              = [ Top_pow_17, TTXNoZ_17, TTZ_17, multiBoson_17, DY_LO_17]
elif year == 2018:
    from StopsDilepton.samples.nanoTuples_Run2018_PromptReco_postProcessed import *
    from StopsDilepton.samples.nanoTuples_Autumn18_postProcessed import *
    Top_pow, TTXNoZ, TTZ_LO, multiBoson, DY_HT_LO = Top_pow_18, TTXNoZ_18, TTZ_18, multiBoson_18, DY_LO_18
    mc              = [ Top_pow_18, TTXNoZ_18, TTZ_18, multiBoson_18, DY_LO_18]

try:
  data_sample = eval(args.era)
except Exception as e:
  logger.error( "Didn't find %s", args.era )
  raise e

#if args.recoil:
#    from Analysis.Tools.RecoilCorrector import RecoilCorrector
#    if args.recoil == "nvtx":
#        recoilCorrector = RecoilCorrector( os.path.join( "/afs/hephy.at/data/rschoefbeck01/StopsDilepton/results/", "recoil_v4.3_fine_nvtx_loop", "%s_lepSel-njet1p-btag0-relIso0.12-looseLeptonVeto-mll20-onZ_recoil_fitResults_SF.pkl"%args.era ) )
#    elif args.recoil == "VUp":
#        recoilCorrector = RecoilCorrector( os.path.join( "/afs/hephy.at/data/rschoefbeck01/StopsDilepton/results/", "recoil_v4.3_fine_VUp_loop", "%s_lepSel-njet1p-btag0-relIso0.12-looseLeptonVeto-mll20-onZ_recoil_fitResults_SF.pkl"%args.era ) )
#    elif args.recoil is "Central":
#        recoilCorrector = RecoilCorrector( os.path.join( "/afs/hephy.at/data/rschoefbeck01/StopsDilepton/results/", "recoil_v4.3_fine", "%s_lepSel-njet1p-btag0-relIso0.12-looseLeptonVeto-mll20-onZ_recoil_fitResults_SF.pkl"%args.era ) )
#

# Text on the plots
def drawObjects( plotData, dataMCScale, lumi_scale ):
    tex = ROOT.TLatex()
    tex.SetNDC()
    tex.SetTextSize(0.04)
    tex.SetTextAlign(11) # align right
    lines = [
      (0.15, 0.95, 'CMS Preliminary' if plotData else 'CMS Simulation'), 
      (0.45, 0.95, 'L=%3.1f fb{}^{-1} (13 TeV) Scale %3.2f'% ( lumi_scale, dataMCScale ) ) if False else (0.45, 0.95, 'L=%3.1f fb{}^{-1} (13 TeV)' % lumi_scale)
    ]
    return [tex.DrawLatex(*l) for l in lines] 

# Read variables and sequences
read_variables = ["weight/F", "l1_pt/F", "l2_pt/F", "l1_eta/F" , "l1_phi/F", "l2_eta/F", "l2_phi/F", "JetGood[pt/F,eta/F,phi/F]", "dl_mass/F", "dl_eta/F", "dl_mt2ll/F", "dl_mt2bb/F", "dl_mt2blbl/F",
                  "met_pt/F", "met_phi/F", "dl_pt/F", "dl_phi/F",
                  #"LepGood[pt/F,eta/F,miniRelIso/F]", "nGoodMuons/F", "nGoodElectrons/F", "l1_mIsoWP/F", "l2_mIsoWP/F",
                  "metSig/F", "ht/F", "nBTag/I", "nJetGood/I","run/I","event/l"]

sequence = []
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
#        # compute fake MET 
#        fakeMET_x = event.met_pt*cos(event.met_phi) - event.GenMET_pt*cos(event.GenMET_phi)
#        fakeMET_y = event.met_pt*sin(event.met_phi) - event.GenMET_pt*sin(event.GenMET_phi)
#
#        fakeMET = sqrt( fakeMET_x**2 + fakeMET_y**2 )
#        fakeMET_phi = atan2( fakeMET_y, fakeMET_x )
#
#        # project fake MET on qT
#        fakeMET_para = fakeMET*cos( fakeMET_phi - ref_phi )
#        fakeMET_perp = fakeMET*cos( fakeMET_phi - ( ref_phi - pi/2) )
#
#        fakeMET_para_corr = - recoilCorrector.predict_para( ref_phi, qt, -fakeMET_para )
#        fakeMET_perp_corr = - recoilCorrector.predict_perp( ref_phi, qt, -fakeMET_perp )
#
#        # rebuild fake MET vector
#        fakeMET_px_corr = fakeMET_para_corr*cos(ref_phi) + fakeMET_perp_corr*cos(ref_phi - pi/2)
#        fakeMET_py_corr = fakeMET_para_corr*sin(ref_phi) + fakeMET_perp_corr*sin(ref_phi - pi/2)
#        #print "%s qt: %3.2f para %3.2f->%3.2f perp %3.2f->%3.2f fakeMET(%3.2f,%3.2f) -> (%3.2f,%3.2f)" % ( sample.name, qt, fakeMET_para, fakeMET_para_corr, fakeMET_perp, fakeMET_perp_corr, fakeMET, fakeMET_phi, sqrt( fakeMET_px_corr**2+fakeMET_py_corr**2), atan2( fakeMET_py_corr, fakeMET_px_corr) )
#   
#        for var in [""] + jme_systematics:
#            if var: var = "_"+var
#            met_px_corr = getattr(event, "met_pt"+var)*cos(getattr(event, "met_phi"+var)) - fakeMET_x + fakeMET_px_corr 
#            met_py_corr = getattr(event, "met_pt"+var)*sin(getattr(event, "met_phi"+var)) - fakeMET_y + fakeMET_py_corr
#    
#            setattr(event, "met_pt_corr"+var, sqrt( met_px_corr**2 + met_py_corr**2 ) )
#            setattr(event, "met_phi_corr"+var, atan2( met_py_corr, met_px_corr ) )
#            
#            mt2Calculator.setLeptons(event.l1_pt, event.l1_eta, event.l1_phi, event.l2_pt, event.l2_eta, event.l2_phi)
#            mt2Calculator.setMet(getattr(event,"met_pt_corr"+var), getattr(event,"met_phi_corr"+var))
#            setattr(event, "dl_mt2ll_corr"+var, mt2Calculator.mt2ll() )
#
#    else:
#        event.met_pt_corr  = event.met_pt 
#        event.met_phi_corr = event.met_phi
#
#        mt2Calculator.setLeptons(event.l1_pt, event.l1_eta, event.l1_phi, event.l2_pt, event.l2_eta, event.l2_phi)
#        mt2Calculator.setMet(event.met_pt_corr, event.met_phi_corr)
#        event.dl_mt2ll_corr =  mt2Calculator.mt2ll()
#
#    #print event.dl_mt2ll, event.dl_mt2ll_corr
#
#sequence.append( corr_recoil )

signals = []

# selection
offZ = "&&abs(dl_mass-91.1876)>15" if not (args.selection.count("onZ") or args.selection.count("allZ") or args.selection.count("offZ")) else ""
def getLeptonSelection( mode ):
  if   mode=="mumu": return "nGoodMuons==2&&nGoodElectrons==0&&isOS&&isMuMu" + offZ
  elif mode=="mue":  return "nGoodMuons==1&&nGoodElectrons==1&&isOS&&isEMu"
  elif mode=="ee":   return "nGoodMuons==0&&nGoodElectrons==2&&isOS&&isEE" + offZ

allPlots   = {}

logger.info('Working on mode ' + str(args.mode))

if year == 2016:
    data_sample = Run2016
    data_sample.texName = "data (2016)"
elif year == 2017:
    data_sample = Run2017
    data_sample.texName = "data (2017)"
elif year == 2018:
    data_sample = Run2018
    data_sample.texName = "data (2018)"

data_sample.setSelectionString([getFilterCut(isData=True, year=year), getLeptonSelection(args.mode)])
data_sample.name           = "data"
data_sample.read_variables = ["event/I","run/I"]
data_sample.style          = styles.errorStyle(ROOT.kBlack)
data_sample.scale          = 1.
lumi_scale                 = data_sample.lumi/1000

logger.info('Lumi scale is ' + str(lumi_scale))

for sample in mc:
    sample.scale           = lumi_scale
    sample.style           = styles.fillStyle(sample.color, lineColor = sample.color)
    sample.read_variables  = ['Pileup_nTrueInt/F', 'GenMET_pt/F', 'GenMET_phi/F'] + variations[args.variation]['read_variables']
    sample.setSelectionString([getFilterCut(isData=False, year=year), getLeptonSelection(args.mode)])

if args.small:
  for sample in mc:
    sample.normalization = 1.
    #sample.reduceFiles( factor = 40 )
    sample.reduceFiles( to = 1 )
    sample.scale /= sample.normalization

  # data
  data_sample.normalization = 1.
  data_sample.reduceFiles( to = 1 )
  data_sample.scale /= data_sample.normalization

#    # Apply scale factors in the mt2ll > 100 GeV signal region (except Top which will be already scaled anyway)
#    if args.selection.count('njet2-btagM-looseLeptonVeto-mll20-met80-metSig5-dPhiJet0-dPhiJet1-mt2ll100') and False: # Turn on when scalefactors are rederived
#      if sample == DY_HT_LO:   sample.scale = lumi_scale*1.30
#      if sample == multiBoson: sample.scale = lumi_scale*1.45
#      if sample == TTZ_LO:     sample.scale = lumi_scale*0.89

  #if args.signal == "T2tt":
  #  for s in signals:
  #    s.scale          = lumi_scale
  #    s.read_variables = ['reweightDilepTrigger/F','reweightLeptonSF/F','reweightLeptonFastSimSF/F','reweightBTag_SF/F','reweightPU36fb/F','Pileup_nTrueInt/F']
  #    s.weight         = lambda event, sample: event.reweightLeptonFastSimSF
  #    s.setSelectionString([getFilterCut(isData=False, year=year), getLeptonSelection(args.mode)])

  #if args.signal == "DM":
  #  for s in signals:
  #    s.scale          = lumi_scale
  #    s.read_variables = ['reweightDilepTrigger/F','reweightLeptonSF/F','reweightBTag_SF/F','reweightPU36fb/F','Pileup_nTrueInt/F']
  #    s.setSelectionString([getFilterCut(isData=False, year=year), getLeptonSelection(args.mode)])

# Use some defaults
Plot.setDefaults( selectionString = cutInterpreter.cutString(args.selection) )

selectionModifier = variations[args.variation]['selectionModifier']

# Stack
stack_mc   = Stack( mc )
stack_data = Stack( data_sample )

plots = []

mt2llBinning = [0,20,40,60,80,100,140,240,340]

if args.variation == 'central':
    dl_mt2ll_data  = Plot(
        name = "dl_mt2ll_data",
        texX = 'M_{T2}(ll) (GeV)', texY = 'Number of Events / 20 GeV' if args.normalizeBinWidth else "Number of Events",
        binning=Binning.fromThresholds(mt2llBinning),
        stack = stack_data,
        attribute = TreeVariable.fromString( "dl_mt2ll/F" ),
        weight = data_weight )
    plots.append( dl_mt2ll_data )

dl_mt2ll_mc  = Plot(\
    name            = "dl_mt2ll" if sys is None else "dl_mt2ll_mc_%s" % sys,
    texX            = 'M_{T2}(ll) (GeV)', texY = 'Number of Events / 20 GeV' if args.normalizeBinWidth else "Number of Events",
    binning         = Binning.fromThresholds(mt2llBinning),
    stack           = stack_mc,
    attribute       = TreeVariable.fromString( selectionModifier("dl_mt2ll/F")),
    selectionString = selectionModifier(cutInterpreter.cutString(args.selection)),
    weight          = MC_WEIGHT( variation = variations[args.variation], returntype='func') )
plots.append( dl_mt2ll_mc )

#  #mt2llCorrBinning = [0,20,40,60,80,100,140,240,340]
#  #
#  #dl_mt2llCorr_data  = Plot(
#  #    name = "dl_mt2ll_corr_data",
#  #    texX = 'M_{T2}(ll) (GeV)', texY = 'Number of Events / 20 GeV' if args.normalizeBinWidth else "Number of Events",
#  #    binning=Binning.fromThresholds(mt2llBinning),
#  #    stack = stack_data,
#  #    attribute = TreeVariable.fromString( "dl_mt2ll/F" ),
## #?     attribute = lambda event, sample: event.dl_mt2ll_corr,
#  #    weight = data_weight,
#  #    )
#  #plots.append( dl_mt2llCorr_data )
#
#  #dl_mt2llCorr_mc  = { sys:Plot(\
#  #    name            = "dl_mt2ll_corr" if sys is None else "dl_mt2ll_corr_mc_%s" % sys,
#  #    texX            = 'corrected M_{T2}(ll) (GeV)', texY = 'Number of Events / 20 GeV' if args.normalizeBinWidth else "Number of Events",
#  #    binning         = Binning.fromThresholds(mt2llBinning),
#  #    stack           = sys_stacks[sys],
#  #    attribute       = lambda event, sample: event.dl_mt2ll_corr if sys is None or sys in weight_systematics else getattr(event, "dl_mt2ll_corr_"+sys) ,
## #     attribute        = TreeVariable.fromString( "dl_mt2ll/F" ) if sys is None or sys in weight_systematics else TreeVariable.fromString( "dl_mt2ll_%s/F" % sys ),
#  #    selectionString = addSys(cutInterpreter.cutString(args.selection), sys),
#  #    weight          = MC_WEIGHT( sys = sys )[0],
#  #    ) for sys in all_systematics }
#  #plots.extend( dl_mt2llCorr_mc.values() )
#
#  if args.selection.count('njet2'):
#
#    dl_mt2blbl_data  = Plot( 
#        name = "dl_mt2blbl_data",
#        texX = 'M_{T2}(blbl) (GeV)', texY = 'Number of Events / 20 GeV' if args.normalizeBinWidth else "Number of Events",
#        stack = stack_data,
#        attribute = TreeVariable.fromString( "dl_mt2blbl/F" ),
#        binning=Binning.fromThresholds([0,20,40,60,80,100,120,140,160,200,250,300,350]),
#        weight = data_weight,
#        ) 
#    plots.append( dl_mt2blbl_data )
#
#    dl_mt2blbl_mc  = {sys: Plot(
#        name = "dl_mt2blbl" if sys is None else "dl_mt2blbl_mc_%s" % sys,
#        texX = 'M_{T2}(blbl) (GeV)', texY = 'Number of Events / 20 GeV' if args.normalizeBinWidth else "Number of Events",
#        stack = sys_stacks[sys],
#        attribute = TreeVariable.fromString( "dl_mt2blbl/F" ) if sys is None or sys in weight_systematics else TreeVariable.fromString( "dl_mt2blbl_%s/F" % sys ),
#        binning=Binning.fromThresholds([0,20,40,60,80,100,120,140,160,200,250,300,350]),
#        selectionString = addSys(cutInterpreter.cutString(args.selection), sys),
#        weight = MC_WEIGHT( sys = sys )[0],
#        ) for sys in all_systematics }
#    plots.extend( dl_mt2blbl_mc.values() )
#
#  nBtagBinning = [5, 1, 6] if args.selection.count('btag1p') else [1,0,1]
#
#  nbtags_data  = Plot( 
#      name = "nbtags_data",
#      texX = 'number of b-tags (CSVM)', texY = 'Number of Events',
#      stack = stack_data,
#      attribute = TreeVariable.fromString('nBTag/I'),
#      binning=nBtagBinning,
#      weight = data_weight,
#      ) 
#  plots.append( nbtags_data )
#
#  nbtags_mc  = {sys: Plot(
#      name = "nbtags" if sys is None else "nbtags_mc_%s" % sys,
#      texX = 'number of b-tags (CSVM)', texY = 'Number of Events',
#      stack = sys_stacks[sys],
#      attribute = TreeVariable.fromString('nBTag/I') if sys is None or sys in weight_systematics or sys in met_systematics else TreeVariable.fromString( "nBTag_%s/I" % sys ),
#      binning=nBtagBinning,
#      selectionString = addSys(cutInterpreter.cutString(args.selection), sys),
#      weight = MC_WEIGHT( sys = sys )[0],
#      ) for sys in all_systematics }
#  plots.extend( nbtags_mc.values() )
#
#  jetBinning = [8,2,10] if args.selection.count('njet2') else [2,0,2]
#
#  njets_data  = Plot( 
#      name = "njets_data",
#      texX = 'number of jets', texY = 'Number of Events',
#      stack = stack_data,
#      attribute = TreeVariable.fromString('nJetGood/I'),
#      binning=jetBinning,
#      weight = data_weight,
#      )
#  plots.append( njets_data )
#
#  njets_mc  = {sys: Plot(
#      name = "njets" if sys is None else "njets_mc_%s" % sys,
#      texX = 'number of jets', texY = 'Number of Events',
#      stack = sys_stacks[sys],
#      attribute = TreeVariable.fromString('nJetGood/I') if sys is None or sys in weight_systematics or sys in met_systematics else TreeVariable.fromString( "nJetGood_%s/I" % sys ),
#      binning= jetBinning,
#      selectionString = addSys(cutInterpreter.cutString(args.selection), sys),
#      weight = MC_WEIGHT( sys = sys )[0],
#      ) for sys in all_systematics }
#  plots.extend( njets_mc.values() )
#
#  metBinning = [0,20,40,60,80] if args.selection.count('metInv') else [80,130,180,230,280,320,420,520,800] if args.selection.count('met80') else [0,80,130,180,230,280,320,420,520,800]
#
#  met_data  = Plot( 
#      name = "met_data",
#      texX = 'E_{T}^{miss} (GeV)', texY = 'Number of Events / 50 GeV' if args.normalizeBinWidth else "Number of Events",
#      stack = stack_data, 
#      attribute = TreeVariable.fromString( "met_pt/F" ),
#      binning=Binning.fromThresholds( metBinning ),
#      weight = data_weight,
#      )
#  plots.append( met_data )
#
#  met_mc  = {sys: Plot(
#      name = "met_pt" if sys is None else "met_pt_mc_%s" % sys,
#      texX = 'E_{T}^{miss} (GeV)', texY = 'Number of Events / 50 GeV' if args.normalizeBinWidth else "Number of Events",
#      stack = sys_stacks[sys],
#      attribute = TreeVariable.fromString('met_pt/F') if sys not in met_systematics else TreeVariable.fromString( "met_pt_%s/F" % sys ),
#      binning=Binning.fromThresholds( metBinning ),
#      selectionString = addSys(cutInterpreter.cutString(args.selection), sys),
#      weight = MC_WEIGHT( sys = sys )[0],
#      ) for sys in all_systematics }
#  plots.extend( met_mc.values() )
#
#  metBinning2 = [0,20,40,60,80] if args.selection.count('metInv') else [80,100,120,140,160,200,500] if args.selection.count('met80') else [0,80,100,120,140,160,200,500]
#
#  met2_data  = Plot(
#      name = "met2_data",
#      texX = 'E_{T}^{miss} (GeV)', texY = 'Number of Events / 20 GeV' if args.normalizeBinWidth else "Number of Events",
#      stack = stack_data,
#      attribute = TreeVariable.fromString( "met_pt/F" ),
#      binning=Binning.fromThresholds( metBinning2 ),
#      weight = data_weight,
#      )
#  plots.append( met2_data )
#
#  met2_mc  = {sys: Plot(
#      name = "met2_pt" if sys is None else "met2_pt_mc_%s" % sys,
#      texX = 'E_{T}^{miss} (GeV)', texY = 'Number of Events / 20 GeV' if args.normalizeBinWidth else "Number of Events",
#      stack = sys_stacks[sys],
#      attribute = TreeVariable.fromString('met_pt/F') if sys not in met_systematics else TreeVariable.fromString( "met_pt_%s/F" % sys ),
#      binning=Binning.fromThresholds( metBinning2 ),
#      selectionString = addSys(cutInterpreter.cutString(args.selection), sys),
#      weight = MC_WEIGHT( sys = sys )[0],
#      ) for sys in all_systematics }
#  plots.extend( met2_mc.values() )
#
#  metSigBinning = [0,2,4,6,8,10,12] if args.selection.count('POGMetSig0To12') else [12,16,20,24,28,32,36,40,45,50,55,60,65,70,75,80,85,90,95,100] if args.selection.count('POGMetSig12') else [0,4,8,12,16,20,24,28,32,36,40,45,50,55,60,65,70,75,80,85,90,95,100]
#
#  metSig_data  = Plot( 
#      name = "MET_significance_data",
#      texX = 'E_{T}^{miss} significance (GeV)', texY = 'Number of Events / 5 GeV' if args.normalizeBinWidth else "Number of Events",
#      stack = stack_data, 
#      attribute = TreeVariable.fromString( "MET_significance/F" ),
#      binning=Binning.fromThresholds( metSigBinning ),
#      weight = data_weight,
#      )
#  plots.append( metSig_data )
#
#  metSig_mc  = {sys: Plot(
#      name = "MET_Significance" if sys is None else "MET_significance_mc_%s" % sys,
#      texX = 'E_{T}^{miss} significance (GeV)', texY = 'Number of Events / 5 GeV' if args.normalizeBinWidth else "Number of Events",
#      stack = sys_stacks[sys],
#      attribute = TreeVariable.fromString('MET_significance/F') if sys not in jme_systematics else TreeVariable.fromString( "MET_significance_%s/F" % sys ),
#      binning=Binning.fromThresholds( metSigBinning ),
#      selectionString = addSys(cutInterpreter.cutString(args.selection), sys),
#      weight = MC_WEIGHT( sys = sys )[0],
#      ) for sys in all_systematics }
#  plots.extend( metSig_mc.values() )

# Make plot directory
plot_directory_ = os.path.join(plot_directory, 'systematicPlots', args.plot_directory, args.selection, args.era, args.mode)
result_file = os.path.join(plot_directory_, 'results.pkl')

try: 
    os.makedirs(plot_directory_)
except: 
    pass

# Fire up the cache
dirDB = DirDB(os.path.join(plot_directory, 'systematicPlots', args.plot_directory, args.selection, 'cache'))

key  = ( args.era, args. mode )

if dirDB.contains(key):
    normalisation_mc, normalisation_data, histos = dirDB.get( key )
    for i_p, h_s in enumerate(histos):
        plots[i_p].histos = h_s
else:
    logger.info( "Calculating normalisations.")
    # Calculate the normalisation yield
    normalization_selection_string = selectionModifier(cutInterpreter.cutString(args.selection + '-mt2llTo100'))
    mc_normalization_weight_string    = MC_WEIGHT(variations[args.variation], returntype='string')
    normalisation_mc = {s.name :s.scale*s.getYieldFromDraw(selectionString = normalization_selection_string, weightString = mc_normalization_weight_string)['val'] for s in mc}

    if args.variation == 'central':
        normalisation_data = data_sample.getYieldFromDraw( selectionString = normalization_selection_string, weightString = data_weight_string)['val']
    else:
        normalisation_data = -1

    logger.info( "Making plots.")
    plotting.fill(plots, read_variables = read_variables, sequence = sequence)

    #    allPlots.update({p.name : p.histos for p in plots})
    #    yields.update(yield_mc)

    logger.info( "Done with %s in channel %s.", args.variation, args.mode)
    # Delete lambda because we can't serialize it
    for plot in plots:
        del plot.weight
    dirDB.add( key, (normalisation_mc, normalisation_data, [plot.histos for plot in plots]) )

#plotConfigs = [\
#       [ dl_mt2ll_mc, dl_mt2ll_data, 20],
#       [ dl_mt2llCorr_mc, dl_mt2llCorr_data, 20],
#       [ nbtags_mc, nbtags_data, -1],
#       [ njets_mc, njets_data, -1],
#       [ met_mc, met_data, 50],
#       [ met2_mc, met2_data, 20],
#       [ metSig_mc, metSig_data, 5],
#  ]
#if args.selection.count('njet2'):
#    plotConfigs.append([ dl_mt2blbl_mc, dl_mt2blbl_data, 20])
#
#

#else:
#  (allPlots, yields) = pickle.load(file( result_file ))
#
#  from RootTools.plot.Plot import addOverFlowBin1D
#  for p in plots:
#    p.histos = allPlots[p.name]
#    for s in p.histos:
#      for h in s:
#        addOverFlowBin1D(h, "upper")
#        if h.Integral()==0: logger.warning( "Found empty histogram %s in results file %s", h.GetName(), result_file )
#
#  topName = Top_pow.name
#  top_sf = {}
#
#  dataMCScaling = True
#  if dataMCScaling:
#    yield_data    = yields['data']
#    yield_non_top = sum(yields[s.name + 'None'] for s in mc if s.name != topName)
#    top_sf[None]  = (yield_data - yield_non_top)/yields[topName+'None']
#    total         = yield_data
#    logger.info( "Data: %i MC TT %3.2f MC other %3.2f SF %3.2f", yield_data, yields[topName+'None'], yield_non_top, top_sf[None] )
#    if args.selection.count('njet01-btag0-looseLeptonVeto-mll20-metInv') and mode != "mue":
#      top_sf[None] = 1
#  else:
#    top_sf[None] = 1
#    total        = sum(yield_mc.values())
#
#  #Scaling systematic shapes to MT2ll<100 region
#  for sys_pair in sys_pairs:
#    for sys in sys_pair[1:]:
#      if not top_sf.has_key( sys ):
#          mc_sys_weight_func, mc_sys_weight_string = MC_WEIGHT( sys = sys )
#          non_top                                  = sum(yields[s.name+sys] for s in mc if s.name != topName)
#          top_sf[sys]                              = (total - non_top)/yields[topName+sys]
#          logger.info( "Total: %i sys %s MC TT %3.2f MC other %3.2f SF %3.2f", total, sys, yields[topName+sys], non_top, top_sf[sys] )
#
#          if args.selection.count('njet01-btag0-looseLeptonVeto-mll20-metInv') and mode != "mue":
#            top_sf[sys] = 1
#            logger.info( "NOT scaling top for " + args.selection + " (mode " + mode + ")" )
#
#  for plot_mc, plot_data, bin_width in plotConfigs:
#    if args.normalizeBinWidth and bin_width>0:
#      for p in plot_mc.values() + [plot_data]:
#        for histo in sum(p.histos, []): 
#          for ib in range(histo.GetXaxis().GetNbins()+1):
#            val = histo.GetBinContent( ib )
#            err = histo.GetBinError( ib )
#            width = histo.GetBinWidth( ib )
#            histo.SetBinContent(ib, val / (width / bin_width)) 
#            histo.SetBinError(ib, err / (width / bin_width)) 
#    topHist = None
#    ttzHist = None
#    ttxHist = None
#    mbHist  = None
#    dyHist  = None
#
#    # Scaling Top
#    for k in plot_mc.keys():
#      for s in plot_mc[k].histos:
#    #    for h in s:
#    #      h.Scale(lumi_scale)
#        pos_top = [i for i,x in enumerate(mc) if x == Top_pow][0]
#        pos_ttz = [i for i,x in enumerate(mc) if x == TTZ_LO][0]
#        pos_ttx = [i for i,x in enumerate(mc) if x == TTXNoZ][0]
#        pos_dy  = [i for i,x in enumerate(mc) if x == DY_HT_LO][0]
#        pos_mb  = [i for i,x in enumerate(mc) if x == multiBoson][0]
#        plot_mc[k].histos[0][pos_top].Scale(top_sf[k])
#        topHist = plot_mc[k].histos[0][pos_top]
#        ttzHist = plot_mc[k].histos[0][pos_ttz]
#        ttxHist = plot_mc[k].histos[0][pos_ttx]
#        mbHist  = plot_mc[k].histos[0][pos_mb]
#        dyHist  = plot_mc[k].histos[0][pos_dy]
#        
#    #Calculating systematics
#    h_summed = {k: plot_mc[k].histos_added[0][0].Clone() for k in plot_mc.keys()}
#
#    ##Normalize systematic shapes
#    #if args.sysScaling:
#    #    for k in h_summed.keys():
#    #        if k is None: continue
#    #        h_summed[k].Scale( top_sf[ k ] )
#
#    h_rel_err = h_summed[None].Clone()
#    h_rel_err.Reset()
#
#    #MC statistical error
#    for ib in range( 1 + h_rel_err.GetNbinsX() ):
#        h_rel_err.SetBinContent(ib, h_summed[None].GetBinError(ib)**2 )
#
#    h_sys = {}
#    goOn = False
#    for k, s1, s2 in ([s for s in sys_pairs if s[0] == args.showOnly] if args.showOnly else sys_pairs):
#      goOn = True
#      h_sys[k] = h_summed[s1].Clone()
#      h_sys[k].Scale(-1)
#      h_sys[k].Add(h_summed[s2])
#    if not goOn: continue
#
#    # Adding in quadrature
#    for k in h_sys.keys():
#        for ib in range( 1 + h_rel_err.GetNbinsX() ):
#          h_rel_err.SetBinContent(ib, h_rel_err.GetBinContent(ib) + h_sys[k].GetBinContent(ib)**2 )
#
#    # When making plots with mt2ll > 100 GeV, include also our background shape uncertainties
#    if args.selection.count('mt2ll100') or plot_mc == dl_mt2ll_mc and False:
#     for ib in range(1 + h_rel_err.GetNbinsX() ):
#       if plot_mc == dl_mt2ll_mc and h_rel_err.GetBinCenter(ib) < 100: continue
#       topUnc = 1 if (plot_mc == dl_mt2ll_mc and h_rel_err.GetBinCenter(ib) > 240) else 0.5
#       h_rel_err.SetBinContent(ib, h_rel_err.GetBinContent(ib) + (topUnc*topHist.GetBinContent(ib))**2 )
#       h_rel_err.SetBinContent(ib, h_rel_err.GetBinContent(ib) + (0.2*ttzHist.GetBinContent(ib))**2 )
#       h_rel_err.SetBinContent(ib, h_rel_err.GetBinContent(ib) + (0.25*ttxHist.GetBinContent(ib))**2 )
#       h_rel_err.SetBinContent(ib, h_rel_err.GetBinContent(ib) + (0.25*dyHist.GetBinContent(ib))**2 )
#       h_rel_err.SetBinContent(ib, h_rel_err.GetBinContent(ib) + (0.25*mbHist.GetBinContent(ib))**2 )
#
#    # take sqrt
#    for ib in range( 1 + h_rel_err.GetNbinsX() ):
#        h_rel_err.SetBinContent(ib, sqrt( h_rel_err.GetBinContent(ib) ) )
#
#    # Divide
#    h_rel_err.Divide(h_summed[None])
#
#    plot = plot_mc[None]
#    if args.normalizeBinWidth: plot.name += "_normalizeBinWidth"
#    signal_histos = plot_data.histos[1:]
#    data_histo    = plot_data.histos[0][0]
#    for h in plot_data.histos[0][1:]:
#      data_histo.Add(h)
#
#    data_histo.style = styles.errorStyle( ROOT.kBlack )
#    plot.histos += [[ data_histo ]]
#    for h in signal_histos: plot.histos += [h]
#    plot_data.stack[0][0].texName = data_sample.texName
#    plot.stack += [[ plot_data.stack[0][0] ]]
#    for i, signal in enumerate(signals):
#      plot_data.stack[i+1][0].texName = signal.texName
#      plot_data.stack[i+1][0].style   = signal.style
#      plot.stack += [[ plot_data.stack[i+1][0] ]]
#
#    boxes = []
#    ratio_boxes = []
#    for ib in range(1, 1 + h_rel_err.GetNbinsX() ):
#        val = h_summed[None].GetBinContent(ib)
#        if val<0: continue
#        sys = h_rel_err.GetBinContent(ib)
#        box = ROOT.TBox( h_rel_err.GetXaxis().GetBinLowEdge(ib),  max([0.03, (1-sys)*val]), h_rel_err.GetXaxis().GetBinUpEdge(ib), max([0.03, (1+sys)*val]) )
#        box.SetLineColor(ROOT.kBlack)
#        box.SetFillStyle(3444)
#        box.SetFillColor(ROOT.kBlack)
#        r_box = ROOT.TBox( h_rel_err.GetXaxis().GetBinLowEdge(ib),  max(0.1, 1-sys), h_rel_err.GetXaxis().GetBinUpEdge(ib), min(1.9, 1+sys) )
#        r_box.SetLineColor(ROOT.kBlack)
#        r_box.SetFillStyle(3444)
#        r_box.SetFillColor(ROOT.kBlack)
#
#        boxes.append( box )
#        ratio_boxes.append( r_box )
#
#        ratio = {'yRange':(0.1,1.9), 'drawObjects':ratio_boxes}
#           
#
#    for log in [False, True]:
#      plotDir = os.path.join(plot_directory, 'systematicPlots', args.plot_directory, args.selection, str(year), mode + ("_log" if log else "") + "_scaled")
#      #plotDir = os.path.join(plot_directory, args.plot_directory,  mode + ("_log" if log else "") + "_scaled", args.selection)
#      if args.showOnly: plotDir = os.path.join(plotDir, "only_" + args.showOnly)
#      plotting.draw(plot,
#          plot_directory = plotDir,
#          ratio = ratio,
#          legend = (0.50,0.88-0.04*sum(map(len, plot.histos)),0.95,0.88),
#          logX = False, logY = log, #sorting = True,
#          yRange = (0.03, "auto"),
#          drawObjects = drawObjects( True, top_sf[None], lumi_scale ) + boxes,
#          copyIndexPHP = True
#      )
