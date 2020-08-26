#!/usr/bin/env python
''' analysis script for standard plots with systematic errors
'''

# quick marks
#NN     implementation of input uncertainties
#

# Standard imports and batch mode
import ROOT
ROOT.gROOT.SetBatch(True)
import operator
import pickle, os, time, sys, copy
from math                                import sqrt, cos, sin, pi, atan2, isnan

# RootTools
from RootTools.core.standard             import *

#Analysis / StopsDilepton / Samples
from StopsDilepton.tools.user            import plot_directory
from StopsDilepton.tools.helpers         import deltaPhi, add_histos
from Analysis.Tools.metFilters           import getFilterCut
from StopsDilepton.tools.cutInterpreter  import cutInterpreter
#from StopsDilepton.tools.RecoilCorrector import RecoilCorrector
from StopsDilepton.tools.mt2Calculator   import mt2Calculator
from Analysis.Tools.puReweighting        import getReweightingFunction
from Analysis.Tools.DirDB                import DirDB

# JEC corrector
from JetMET.JetCorrector.JetCorrector    import JetCorrector, correction_levels_data, correction_levels_mc
corrector_data     = JetCorrector.fromTarBalls( [(1, 'Autumn18_RunD_V8_DATA') ], correctionLevels = correction_levels_data )
corrector_mc       = JetCorrector.fromTarBalls( [(1, 'Autumn18_RunD_V8_DATA') ], correctionLevels = correction_levels_mc )


# Arguments
import argparse
argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--logLevel',          action='store',      default='INFO',     nargs='?', choices=['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'TRACE', 'NOTSET'], help="Log level for logging")
argParser.add_argument('--signal',            action='store',      default=None,        nargs='?', choices=['None', "T2tt",'DM'], help="Add signal to plot")
argParser.add_argument('--add',               action='store_true', help="Append missing commands to file (instead of overwriting)")
# define plot 
argParser.add_argument('--era',               action='store', type=str,      default="Run2016")
argParser.add_argument('--plot_directory',    action='store',      default='v1')
argParser.add_argument('--reweightPU',        action='store',      default='Central', choices=[ 'Central', 'VUp'] )
argParser.add_argument('--selection',         action='store',            default='njet2p-btag1p-miniIso0.2-looseLeptonVeto-mll20-met80-metSig5-dPhiJet0-dPhiJet1')
argParser.add_argument('--variation',         action='store',      default=None, help="Which systematic variation to run. Don't specify for producing plots.")
argParser.add_argument('--mode',              action='store',      default = 'all', choices = ['mumu', 'ee', 'mue', 'SF', 'all'],   help='Which mode?')
argParser.add_argument('--normalizeBinWidth', action='store_true', default=False,       help='normalize wider bins?')
argParser.add_argument('--small',             action='store_true',     help='Run only on a small subset of the data?')
argParser.add_argument('--private',           action='store_true',     help='Produce private plots')
# loading samples
argParser.add_argument('--dpm',               action='store_true',     help='Use dpm?', )
argParser.add_argument('--scratch',           action='store_true',     help='Use scratch?', )
# write caches
argParser.add_argument('--overwrite',         action='store_true',     help='Overwrite?')
# Scalings
argParser.add_argument('--scaling',           action='store',      default=None, choices = [None, 'mc', 'top'],     help='Scale top to data in mt2ll<100?')
argParser.add_argument('--variation_scaling', action='store_true', help='Scale the variations individually to mimick bkg estimation?')
argParser.add_argument('--normalize',         action='store_true', help='Perform area normalization mc to data?')
# add extra label to path
argParser.add_argument('--beta',              action='store',      default=None, help="Add an additional directory label for minor changes to the plots")
argParser.add_argument('--newMetSigPlots',    action='store_true', help='')

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

# Load from DPM?
if args.dpm:
    data_directory          = "/dpm/oeaw.ac.at/home/cms/store/user/rschoefbeck/Stops2l-postprocessed/"
if args.scratch:
    data_directory          = "/scratch/robert.schoefbeck/nanoTuples/"
    
if year == 2016:
    from StopsDilepton.samples.nanoTuples_Summer16_postProcessed import *
    from StopsDilepton.samples.nanoTuples_Run2016_17Jul2018_postProcessed import *
    Top_pow, TTXNoZ, TTZ_LO, multiBoson, DY_HT_LO = Top_pow_16, TTXNoZ_16, TTZ_16, multiBoson_16, DY_HT_LO_16
    #Top_pow, TTXNoZ, TTZ_LO, multiBoson = Top_pow_16, TTXNoZ_16, TTZ_16, multiBoson_16
    if args.selection.count("onZ") > 0:
        mc          = [ DY_HT_LO_16, Top_pow_16, multiBoson_16, TTXNoZ_16, TTZ_16 ]
    else:
        mc          = [ Top_pow_16, multiBoson_16, DY_HT_LO_16, TTXNoZ_16, TTZ_16 ]
elif year == 2017:
    from StopsDilepton.samples.nanoTuples_Fall17_postProcessed import *
    from StopsDilepton.samples.nanoTuples_Run2017_31Mar2018_postProcessed import *
    Top_pow, TTXNoZ, TTZ_LO, multiBoson, DY_HT_LO = Top_pow_17, TTXNoZ_17, TTZ_17, multiBoson_17, DY_HT_LO_17
    #Top_pow, TTXNoZ, TTZ_LO, multiBoson = Top_pow_17, TTXNoZ_17, TTZ_17, multiBoson_17
    if args.selection.count("onZ") > 0:
        mc          = [ DY_HT_LO_17, Top_pow_17, multiBoson_17, TTXNoZ_17, TTZ_17 ]
    else:
        mc          = [ Top_pow_17, multiBoson_17, DY_HT_LO_17, TTXNoZ_17, TTZ_17 ]
elif year == 2018:
    from StopsDilepton.samples.nanoTuples_Run2018_PromptReco_postProcessed import *
    from StopsDilepton.samples.nanoTuples_Autumn18_postProcessed import *
    Top_pow, TTXNoZ, TTZ_LO, multiBoson, DY_HT_LO = Top_pow_18, TTXNoZ_18, TTZ_18, multiBoson_18, DY_HT_LO_18
    #Top_pow, TTXNoZ, TTZ_LO, multiBoson = Top_pow_18, TTXNoZ_18, TTZ_18, multiBoson_18
    if args.selection.count("onZ") > 0:
        mc          = [ DY_HT_LO_18, Top_pow_18, multiBoson_18, TTXNoZ_18, TTZ_18 ]
    else:
        mc          = [ Top_pow_18, multiBoson_18, DY_HT_LO_18, TTXNoZ_18, TTZ_18 ]

# postions of MC components in list
position = {s.name:i_s for i_s,s in enumerate(mc)}


# signals
signals = []
if args.signal == "T2tt":
    if year == 2016:
        #data_directory           = "/afs/hephy.at/data/cms06/nanoTuples/"
        #postProcessing_directory = "stops_2016_nano_v0p23/dilep/"
        from StopsDilepton.samples.nanoTuples_FastSim_Summer16_postProcessed import signals_T2tt as jobs
    elif year == 2017:
        #data_directory           = "/afs/hephy.at/data/cms06/nanoTuples/"
        #postProcessing_directory = "stops_2017_nano_v0p23/dilep/"
        from StopsDilepton.samples.nanoTuples_FastSim_Fall17_postProcessed import signals_T2tt as jobs
    elif year == 2018:
        #data_directory           = "/afs/hephy.at/data/cms06/nanoTuples/"
        #postProcessing_directory = "stops_2018_nano_v0p23/dilep/"   
        from StopsDilepton.samples.nanoTuples_FastSim_Autumn18_postProcessed import signals_T2tt as jobs

    jobNames = [ x.name for x in jobs ]
    #print jobNames
    T2tt_800_100 = jobs[jobNames.index("T2tt_800_100")]
    T2tt_350_150 = jobs[jobNames.index("T2tt_350_150")]
    T2tt_750_0   = jobs[jobNames.index("T2tt_750_0")]
    T2tt_600_300 = jobs[jobNames.index("T2tt_600_300")]
    signals = [T2tt_800_100, T2tt_350_150]



def jetSelectionModifier( sys, returntype = "func"):
    #Need to make sure all jet variations of the following observables are in the ntuple
    variiedJetObservables = ['nJetGood', 'nBTag', 'dl_mt2ll', 'dl_mt2blbl', 'MET_significance', 'met_pt'] 
    if returntype == "func":
        def changeCut_( string ):
            for s in variiedJetObservables:
                string = string.replace(s, s+'_'+sys)
            return string
        return changeCut_
    elif returntype == "list":
        return [ v+'_'+sys for v in variiedJetObservables ]

def metSelectionModifier( sys, returntype = 'func'):
    #Need to make sure all MET variations of the following observables are in the ntuple
    variiedMetObservables = ['dl_mt2ll', 'dl_mt2blbl', 'MET_significance', 'met_pt']
    if returntype == "func":
        def changeCut_( string ):
            for s in variiedMetObservables:
                string = string.replace(s, s+'_'+sys)
            return string
        return changeCut_
    elif returntype == "list":
        return [ v+'_'+sys for v in variiedMetObservables ]

# these are the nominal MC weights we always apply
if args.reweightPU == 'Central': 
    nominalMCWeights = ["weight", "reweightLeptonSF", "reweightLeptonHit0SF", "reweightLeptonSip3dSF", "reweightPU", "reweightDilepTrigger", "reweightBTag_SF", "reweightLeptonTrackingSF", "reweightL1Prefire", "reweightHEM"]
if args.reweightPU == 'VUp':
    nominalMCWeights = ["weight", "reweightLeptonSF", "reweightLeptonHit0SF", "reweightLeptonSip3dSF", "reweightPUVUp", "reweightDilepTrigger", "reweightBTag_SF", "reweightLeptonTrackingSF", "reweightL1Prefire", "reweightHEM"]

# weights to use for PU variation
if args.reweightPU == 'Central':
    nominalPuWeight, upPUWeight, downPUWeight = "reweightPU", "reweightPUUp", "reweightPUDown"
elif args.reweightPU == 'VUp':
    nominalPuWeight, upPUWeight, downPUWeight = "reweightPUVUp", "reweightPUVVUp", "reweightPUUp"

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
#NN ----------------------------------------------------------------
    # input uncertainties
    # list
#    elif returntype == "inputList":
#        inputSample, inputWeight, extraSelection = variation['addSampleWeight']
#        weightList = []
#        for sample in mc:
#            sampleWeights = variiedMCWeights
#            if sample.name == inputSample.name:
#                sampleWeights.append(inputWeight)
#            weightList.append(sampleWeights)
#        return [weightList]

    # string
    elif returntype == "inputString":
        inputSample, inputWeight, extraSelection = variation['addSampleWeight']
        print "---- extraSelection", extraSelection

        sampleWeights = variiedMCWeights
        if sample.name == inputSample.name:
            sampleWeights.append(inputWeight)
        return "*".join(sampleWeights)

    # weight function
    elif returntype == "input":
        inputSample, inputWeight, extraSelection = variation['addSampleWeight']
        def weight_( event, sample ):
            ww = 1
            for iw, w in enumerate(variiedMCWeights):
                ww *= getattr(event, w)
            if sample.name == inputSample.name:
 #               print sample.name
                ww *= float(inputWeight)
            return ww
        return weight_
# ------------------------------------------------------------------

def data_weight( event, sample ):
    return event.weight*event.reweightHEM

data_weight_string = "weight*reweightHEM"

def signal_weight( event, sample ):
    #SW
    #return event.weight*event.reweightHEM*getattr(event, "reweightPU"+args.reweightPU if args.reweightPU != "Central" else "reweightPU")*event.reweightDilepTrigger*event.reweightLeptonSF*event.reweightLeptonHit0SF*event.reweightLeptonSip3dSF*event.reweightBTag_SF*event.reweightL1Prefire*event.reweightLeptonTrackingSF*event.reweight_nISR*event.reweightLeptonFastSimSF
    return event.weight*event.reweightHEM*getattr(event, "reweightPU"+args.reweightPU if args.reweightPU != "Central" else "reweightPU")*event.reweightDilepTrigger*event.reweightLeptonSF*event.reweightBTag_SF*event.reweightL1Prefire*event.reweightLeptonTrackingSF*event.reweight_nISR*event.reweightLeptonFastSimSF

# Define all systematic variations
variations = {
    'central'           : {'read_variables': [ '%s/F'%v for v in nominalMCWeights ]},
    'jesTotalUp'        : {'selectionModifier':jetSelectionModifier('jesTotalUp'),               'read_variables' : [ '%s/F'%v for v in nominalMCWeights + jetSelectionModifier('jesTotalUp','list')]},
    'jesTotalDown'      : {'selectionModifier':jetSelectionModifier('jesTotalDown'),             'read_variables' : [ '%s/F'%v for v in nominalMCWeights + jetSelectionModifier('jesTotalDown','list')]},
    'jerUp'             : {'selectionModifier':jetSelectionModifier('jerUp'),                    'read_variables' : [ '%s/F'%v for v in nominalMCWeights + jetSelectionModifier('jerUp','list')]},
    'jerDown'           : {'selectionModifier':jetSelectionModifier('jerDown'),                  'read_variables' : [ '%s/F'%v for v in nominalMCWeights + jetSelectionModifier('jerDown','list')]},
    'unclustEnUp'       : {'selectionModifier':metSelectionModifier('unclustEnUp'),              'read_variables' : [ '%s/F'%v for v in nominalMCWeights + metSelectionModifier('unclustEnUp','list')]},
    'unclustEnDown'     : {'selectionModifier':metSelectionModifier('unclustEnDown'),            'read_variables' : [ '%s/F'%v for v in nominalMCWeights + metSelectionModifier('unclustEnDown','list')]},
    'PUUp'              : {'replaceWeight':(nominalPuWeight,upPUWeight),                         'read_variables' : [ '%s/F'%v for v in nominalMCWeights + [upPUWeight] ]},
    'PUDown'            : {'replaceWeight':(nominalPuWeight,downPUWeight),                       'read_variables' : [ '%s/F'%v for v in nominalMCWeights + [downPUWeight] ]},
    'BTag_SF_b_Down'    : {'replaceWeight':('reweightBTag_SF','reweightBTag_SF_b_Down'),         'read_variables' : [ '%s/F'%v for v in nominalMCWeights + ['reweightBTag_SF_b_Down']]},  
    'BTag_SF_b_Up'      : {'replaceWeight':('reweightBTag_SF','reweightBTag_SF_b_Up'),           'read_variables' : [ '%s/F'%v for v in nominalMCWeights + ['reweightBTag_SF_b_Up'] ]},
    'BTag_SF_l_Down'    : {'replaceWeight':('reweightBTag_SF','reweightBTag_SF_l_Down'),         'read_variables' : [ '%s/F'%v for v in nominalMCWeights + ['reweightBTag_SF_l_Down']]},
    'BTag_SF_l_Up'      : {'replaceWeight':('reweightBTag_SF','reweightBTag_SF_l_Up'),           'read_variables' : [ '%s/F'%v for v in nominalMCWeights + ['reweightBTag_SF_l_Up'] ]},
    'DilepTriggerDown'  : {'replaceWeight':('reweightDilepTrigger','reweightDilepTriggerDown'),  'read_variables' : [ '%s/F'%v for v in nominalMCWeights + ['reweightDilepTriggerDown']]},
    'DilepTriggerUp'    : {'replaceWeight':('reweightDilepTrigger','reweightDilepTriggerUp'),    'read_variables' : [ '%s/F'%v for v in nominalMCWeights + ['reweightDilepTriggerUp']]},
    'LeptonSFDown'      : {'replaceWeight':('reweightLeptonSF','reweightLeptonSFDown'),          'read_variables' : [ '%s/F'%v for v in nominalMCWeights + ['reweightLeptonSFDown']]},
    'LeptonSFUp'        : {'replaceWeight':('reweightLeptonSF','reweightLeptonSFUp'),            'read_variables' : [ '%s/F'%v for v in nominalMCWeights + ['reweightLeptonSFUp']]},
    'LeptonHit0SFDown'  : {'replaceWeight':('reweightLeptonHit0SF','reweightLeptonHit0SFDown'),  'read_variables' : [ '%s/F'%v for v in nominalMCWeights + ['reweightLeptonHit0SFDown']]},
    'LeptonHit0SFUp'    : {'replaceWeight':('reweightLeptonHit0SF','reweightLeptonHit0SFUp'),    'read_variables' : [ '%s/F'%v for v in nominalMCWeights + ['reweightLeptonHit0SFUp']]},
    'LeptonSip3dSFDown' : {'replaceWeight':('reweightLeptonSip3dSF','reweightLeptonSip3dSFDown'),'read_variables' : [ '%s/F'%v for v in nominalMCWeights + ['reweightLeptonSip3dSFDown']]},
    'LeptonSip3dSFUp'   : {'replaceWeight':('reweightLeptonSip3dSF','reweightLeptonSip3dSFUp'),  'read_variables' : [ '%s/F'%v for v in nominalMCWeights + ['reweightLeptonSip3dSFUp']]},
    'L1PrefireDown'     : {'replaceWeight':('reweightL1Prefire','reweightL1PrefireDown'),        'read_variables' : [ '%s/F'%v for v in nominalMCWeights + ['reweightL1PrefireDown']]},
    'L1PrefireUp'       : {'replaceWeight':('reweightL1Prefire','reweightL1PrefireUp'),          'read_variables' : [ '%s/F'%v for v in nominalMCWeights + ['reweightL1PrefireUp']]},
#NN
#    'DYInputUp'             : {'addSampleWeight':(DY_HT_LO, '1.25', '(1)'),                                               'read_variables' : ['%s/F'%v for v in nominalMCWeights] },
#    'DYInputDown'           : {'addSampleWeight':(DY_HT_LO, '0.75', '(1)'),                                               'read_variables' : ['%s/F'%v for v in nominalMCWeights] },
#    'TTInputUp'             : {'addSampleWeight':(Top_pow, '1.25', '(1)'),                                                'read_variables' : ['%s/F'%v for v in nominalMCWeights] },
#    'TTInputDown'           : {'addSampleWeight':(Top_pow, '0.75', '(1)'),                                                'read_variables' : ['%s/F'%v for v in nominalMCWeights] },
#    'MBInputUp'             : {'addSampleWeight':(multiBoson, '1.25', '(1)'),                                             'read_variables' : ['%s/F'%v for v in nominalMCWeights] },
#    'MBInputDown'           : {'addSampleWeight':(multiBoson, '0.75', '(1)'),                                             'read_variables' : ['%s/F'%v for v in nominalMCWeights] },
#    'TTZInputUp'            : {'addSampleWeight':(TTZ_LO, '1.2', '(1)'),                                                 'read_variables' : ['%s/F'%v for v in nominalMCWeights] },
#    'TTZInputDown'          : {'addSampleWeight':(TTZ_LO, '0.8', '(1)'),                                                 'read_variables' : ['%s/F'%v for v in nominalMCWeights] },
#    'OtherInputUp'          : {'addSampleWeight':(TTXNoZ, '1.25', '(1)'),                                                'read_variables' : ['%s/F'%v for v in nominalMCWeights] },
#    'OtherInputDown'        : {'addSampleWeight':(TTXNoZ, '0.75', '(1)'),                                                'read_variables' : ['%s/F'%v for v in nominalMCWeights] },

#
#    'TT1JetMismInputUp'     : {'addSampleWeight':(Top_pow, '1.3', '(Sum$( abs(JetGood_pt - JetGood_genPt) >= 40) >=1)'), 'read_variables' : ['%s/F'%v for v in nominalMCWeights] },
#    'TT1JetMismInputDown'   : {'addSampleWeight':(Top_pow, '0.7', '(Sum$( abs(JetGood_pt - JetGood_genPt) >= 40) >=1)'), 'read_variables' : ['%s/F'%v for v in nominalMCWeights] },
#    'TTTotJetMismInputUp'   : {'addSampleWeight':(Top_pow, '1.15', 
#                              '(Sum$( abs(JetGood_pt - JetGood_genPt) >= 40) ==0 && Sum$(abs(JetGood_pt - JetGood_genPt)) >= 40)'), 
#                                                                                                                         'read_variables' : ['%s/F'%v for v in nominalMCWeights] },
#    'TTTotJetMismInputDown' : {'addSampleWeight':(Top_pow, '0.85',
#                              '(Sum$( abs(JetGood_pt - JetGood_genPt) >= 40) ==0 && Sum$(abs(JetGood_pt - JetGood_genPt)) >= 40)'), 
#                                                                                                                         'read_variables' : ['%s/F'%v for v in nominalMCWeights] },
#    'TTNonPromptInputUp'    : {'addSampleWeight':(Top_pow, '1.5',
#                              '(Sum$( abs(JetGood_pt - JetGood_genPt) >= 40) ==0 && Sum$(abs(JetGood_pt - JetGood_genPt)) < 40 && ((l1_muIndex>=0 && (Muon_genPartFlav[l1_muIndex])!=1) || (l2_muIndex>=0 && (Muon_genPartFlav[l2_muIndex])!=1)))'), 
#                                                                                                                         'read_variables' : ['%s/F'%v for v in nominalMCWeights] },
#    'TTNonPromptInputDown'  : {'addSampleWeight':(Top_pow, '0.5',
#                              '(Sum$( abs(JetGood_pt - JetGood_genPt) >= 40) ==0 && Sum$(abs(JetGood_pt - JetGood_genPt)) < 40 && ((l1_muIndex>=0 && (Muon_genPartFlav[l1_muIndex])!=1) || (l2_muIndex>=0 && (Muon_genPartFlav[l2_muIndex])!=1)))'),
#
#    'TopPt':{},
}

# Add a default selection modifier that does nothing
for key, variation in variations.iteritems():
    if not variation.has_key('selectionModifier'):
        variation['selectionModifier'] = lambda string:string
    if not variation.has_key('read_variables'):
        variation['read_variables'] = [] 

# Check if we know the variation
if args.variation is not None and args.variation not in variations.keys():
    raise RuntimeError( "Variation %s not among the known: %s", args.variation, ",".join( variation.keys() ) )

# arguments & directory
plot_subdirectory = args.plot_directory
if args.small:                    plot_subdirectory += "_small"
if args.signal == "DM":           plot_subdirectory += "_DM"
if args.signal == "T2tt":         plot_subdirectory += "_T2tt"
if args.reweightPU:               plot_subdirectory += "_reweightPU%s"%args.reweightPU
#if args.recoil:                  plot_subdirectory  += '_recoil_'+args.recoil

#if args.recoil:
#    from Analysis.Tools.RecoilCorrector import RecoilCorrector
#    if args.recoil == "nvtx":
#        recoilCorrector = RecoilCorrector( os.path.join( "/afs/hephy.at/data/rschoefbeck01/StopsDilepton/results/", "recoil_v4.3_fine_nvtx_loop", "%s_lepSel-njet1p-btag0-relIso0.12-looseLeptonVeto-mll20-onZ_recoil_fitResults_SF.pkl"%args.era ) )
#    elif args.recoil == "VUp":
#        recoilCorrector = RecoilCorrector( os.path.join( "/afs/hephy.at/data/rschoefbeck01/StopsDilepton/results/", "recoil_v4.3_fine_VUp_loop", "%s_lepSel-njet1p-btag0-relIso0.12-looseLeptonVeto-mll20-onZ_recoil_fitResults_SF.pkl"%args.era ) )
#    elif args.recoil is "Central":
#        recoilCorrector = RecoilCorrector( os.path.join( "/afs/hephy.at/data/rschoefbeck01/StopsDilepton/results/", "recoil_v4.3_fine", "%s_lepSel-njet1p-btag0-relIso0.12-looseLeptonVeto-mll20-onZ_recoil_fitResults_SF.pkl"%args.era ) )
#

# Read variables and sequences
read_variables = ["weight/F", "l1_pt/F", "l2_pt/F", "l1_eta/F" , "l1_phi/F", "l2_eta/F", "l2_phi/F", "JetGood[pt/F,eta/F,phi/F]", "dl_mass/F", "dl_eta/F", "dl_mt2ll/F", "dl_mt2bb/F", "dl_mt2blbl/F",
                  "met_pt/F", "met_phi/F", "dl_pt/F", "dl_phi/F",
#                  "l1_pdgId/I", "l2_pdgId/I",
#                  "Jet[pt/F,rawFactor/F,pt_nom/F,eta/F,area/F]", "run/I", "fixedGridRhoFastjetAll/F",
#                  "nMuon/I", "Muon[dxy/F,dxyErr/F,dz/F,dzErr/F,eta/F,ip3d/F,jetRelIso/F,mass/F,miniPFRelIso_all/F,miniPFRelIso_chg/F,pfRelIso03_all/F,pfRelIso03_chg/F,pfRelIso04_all/F,phi/F,pt/F,ptErr/F,segmentComp/F,sip3d/F,mvaTTH/F,charge/I,jetIdx/I,nStations/I,nTrackerLayers/I,pdgId/I,tightCharge/I,highPtId/b,inTimeMuon/O,isGlobal/O,isPFcand/O,isTracker/O,mediumId/O,mediumPromptId/O,miniIsoId/b,multiIsoId/b,mvaId/b,pfIsoId/b,softId/O,softMvaId/O,tightId/O,tkIsoId/b,triggerIdLoose/O,cleanmask/b]",
                  #"LepGood[pt/F,eta/F,miniRelIso/F]", "nGoodMuons/F", "nGoodElectrons/F", "l1_mIsoWP/F", "l2_mIsoWP/F",
                  "metSig/F", "ht/F", "nBTag/I", "nJetGood/I","run/I","event/l","reweightHEM/F"]

sequence = []

# selection
offZ = "&&abs(dl_mass-91.1876)>15" if not (args.selection.count("onZ") or args.selection.count("allZ") or args.selection.count("offZ")) else ""
def getLeptonSelection( mode ):
  if   mode=="mumu": return "nGoodMuons==2&&nGoodElectrons==0&&isOS&&isMuMu" + offZ
  elif mode=="mue":  return "nGoodMuons==1&&nGoodElectrons==1&&isOS&&isEMu"
  elif mode=="ee":   return "nGoodMuons==0&&nGoodElectrons==2&&isOS&&isEE" + offZ

modes = ['mumu', 'mue', 'ee'] if args.mode=='all' else ['mumu', 'ee'] if args.mode=='SF' else [ args.mode ]
print "running modes: ", modes

logger.info('Working on modes: %s', ','.join(modes))

try:
  data_sample = eval(args.era)
except Exception as e:
  logger.error( "Didn't find %s", args.era )
  raise e

# Define samples
data_sample.name           = "data"
data_sample.read_variables = ["event/I","run/I"]
data_sample.style          = styles.errorStyle(ROOT.kBlack)
data_sample.scale          = 1.
lumi_scale                 = data_sample.lumi/1000
logger.info('Lumi scale is ' + str(lumi_scale))
for sample in mc:
    sample.scale           = lumi_scale
    sample.style           = styles.fillStyle(sample.color, lineColor = sample.color)
    sample.read_variables  = ['Pileup_nTrueInt/F', 'GenMET_pt/F', 'GenMET_phi/F', "l1_muIndex/I", "l2_muIndex/I"]
    # append variables for systematics
    if args.variation is not None:
        sample.read_variables+=list(set(variations[args.variation]['read_variables']))

for i_s, sample in enumerate(signals):
    sample.scale           = lumi_scale
    sample.style           = styles.lineStyle( ROOT.kBlack, width=2, dotted=not bool(i_s-2 % 3), dashed=not bool(i_s-1 % 3) )
    sample.read_variables  = ['Pileup_nTrueInt/F', 'GenMET_pt/F', 'GenMET_phi/F', "l1_muIndex/I", "l2_muIndex/I"]
    sample.read_variables += ['reweightPU/F', 'Pileup_nTrueInt/F', 'reweightDilepTrigger/F','reweightLeptonSF/F',
    #SW
    #'reweightLeptonHit0SF/F','reweightLeptonSip3dSF/F',
    'reweightLeptonFastSimSF/F','reweight_nISR/F','reweightBTag_SF/F', 'reweightLeptonTrackingSF/F', 'GenMET_pt/F', 'GenMET_phi/F', "l1_muIndex/I", "l2_muIndex/I", "l1_eleIndex/I", "l2_eleIndex/I", "reweightHEM/F", "reweightL1Prefire/F"]
    sample.read_variables += ['reweightPU%s/F'%args.reweightPU if args.reweightPU != "Central" else "reweightPU/F"]
    # weight has its own function!
    #sample.weight         = lambda event, sample: getattr(event, "reweightPU"+args.reweightPU if args.reweightPU != "Central" else "reweightPU")*event.reweightDilepTrigger*event.reweightLeptonSF*event.reweightBTag_SF*event.reweightL1Prefire*event.reweightLeptonTrackingSF*event.reweight_nISR*event.reweightLeptonFastSimSF

# reduce if small
if args.small:
  data_sample.normalization = 1.
  data_sample.reduceFiles( factor = 40 )
  #data_sample.reduceFiles( to = 1 )
  data_sample.scale /= data_sample.normalization
  for sample in mc:
    sample.normalization = 1.
    #sample.reduceFiles( factor = 40 )
    sample.reduceFiles( to = 1 )
    sample.scale /= sample.normalization

# Fire up the cache
dirDB = DirDB(os.path.join(plot_directory, 'systematicPlots', args.era, plot_subdirectory, args.selection, 'cache'))

#print "weight"
#print "&&".join([getFilterCut(isData=True, year=year), getLeptonSelection("mumu")])
#print cutInterpreter.cutString(args.selection)

# loop over modes
for mode in modes:
    logger.info('Working on mode: %s', mode)

    # set selection
    data_sample.setSelectionString([getFilterCut(isData=True, year=year,skipVertexFilter=True), getLeptonSelection(mode)])
    for sample in mc + signals:                                                              
        sample.setSelectionString([getFilterCut(isData=False, year=year,skipVertexFilter=True), getLeptonSelection(mode)])


    # Use some defaults
    Plot.setDefaults( selectionString = cutInterpreter.cutString(args.selection), addOverFlowBin = 'upper' )

    # if we're running a variation specify
#NN
    #print "Is input uncertainty (i.e. has key 'addSampleWeight'): ", variations[args.variation].has_key('addSampleWeight')
    if args.variation is not None and  "Input" in args.variation: # input variation
        addSampleWeight   = variations[args.variation]['addSampleWeight']
        #selectionModifier = None
        selectionModifier = variations[args.variation]['selectionModifier']
        mc_weight         = MC_WEIGHT( variation = variations[args.variation], returntype='input')
# ------------------------------------------------------------------------------------------------
    elif args.variation is not None:
        selectionModifier = variations[args.variation]['selectionModifier']
        mc_weight         = MC_WEIGHT( variation = variations[args.variation], returntype='func')
    else:
        selectionModifier = None 
        mc_weight         = None 

    # Stack
    stack_mc   = Stack( mc )
    if args.signal:
        stack_signal = Stack()
        stack_signal.extend([ [s] for s in signals ])
    stack_data = Stack( data_sample )

    plots      = []

    mt2llBinning = [0,20,40,60,80,100,140,240,340]
    if args.variation == 'central':
        dl_mt2ll_data   = Plot(
            name        = "dl_mt2ll_data",
            texX        = 'M_{T2}(ll) (GeV)', texY = 'Number of Events / 20 GeV' if args.normalizeBinWidth else "Number of Events",
            binning     = Binning.fromThresholds(mt2llBinning),
            stack       = stack_data,
            attribute   = TreeVariable.fromString( "dl_mt2ll/F" ),
            weight      = data_weight )
        plots.append( dl_mt2ll_data )
        if args.signal:
            dl_mt2ll_signal   = Plot(
                name        = "dl_mt2ll_signal",
                texX        = 'M_{T2}(ll) (GeV)', texY = 'Number of Events / 20 GeV' if args.normalizeBinWidth else "Number of Events",
                binning     = Binning.fromThresholds(mt2llBinning),
                stack       = stack_signal,
                attribute   = TreeVariable.fromString( "dl_mt2ll/F" ),
                weight      = signal_weight )
            plots.append( dl_mt2ll_signal )
        
    dl_mt2ll_mc  = Plot(\
        name            = "dl_mt2ll_mc",
        texX            = 'M_{T2}(ll) (GeV)', texY = 'Number of Events / 20 GeV' if args.normalizeBinWidth else "Number of Events",
        binning         = Binning.fromThresholds(mt2llBinning),
        stack           = stack_mc,
        attribute       = TreeVariable.fromString( selectionModifier("dl_mt2ll/F"))   if selectionModifier is not None else None,
        selectionString = selectionModifier(cutInterpreter.cutString(args.selection)) if selectionModifier is not None else None,
        weight          = mc_weight )
    plots.append( dl_mt2ll_mc )

    if args.selection.count('njet2'):
        if args.variation == 'central':
            dl_mt2blbl_fine_data  = Plot( 
                name = "dl_mt2blbl_fine_data",
                texX = 'M_{T2}(blbl) (GeV)', texY = 'Number of Events / 30 GeV' if args.normalizeBinWidth else "Number of Events",
                binning=[420/30,0,400],
                stack = stack_data,
                attribute = TreeVariable.fromString( "dl_mt2blbl/F" ),
                weight = data_weight,
                ) 
            plots.append( dl_mt2blbl_fine_data )
            if args.signal:
                dl_mt2blbl_fine_signal  = Plot( 
                    name = "dl_mt2blbl_fine_signal",
                    texX = 'M_{T2}(blbl) (GeV)', texY = 'Number of Events / 30 GeV' if args.normalizeBinWidth else "Number of Events",
                    binning=[420/30,0,400],
                    stack = stack_signal,
                    attribute = TreeVariable.fromString( "dl_mt2blbl/F" ),
                    weight      = signal_weight )
                plots.append( dl_mt2blbl_fine_signal )
    
        dl_mt2blbl_fine_mc  = Plot(
            name = "dl_mt2blbl_fine_mc",
            texX = 'M_{T2}(blbl) (GeV)', texY = 'Number of Events / 30 GeV' if args.normalizeBinWidth else "Number of Events",
            binning=[420/30,0,400],
            stack = stack_mc,
            attribute = TreeVariable.fromString( selectionModifier("dl_mt2blbl/F") )      if selectionModifier is not None else None,
            selectionString = selectionModifier(cutInterpreter.cutString(args.selection)) if selectionModifier is not None else None,
            weight          = mc_weight )
        plots.append( dl_mt2blbl_fine_mc )

        if args.variation == 'central':
            dl_mt2blbl_data  = Plot( 
                name = "dl_mt2blbl_data",
                texX = 'M_{T2}(blbl) (GeV)', texY = 'Number of Events / 20 GeV' if args.normalizeBinWidth else "Number of Events",
                binning=Binning.fromThresholds([0,20,40,60,80,100,120,140,160,200,250,300,350]),
                stack = stack_data,
                attribute = TreeVariable.fromString( "dl_mt2blbl/F" ),
                weight = data_weight,
                ) 
            plots.append( dl_mt2blbl_data )
            if args.signal:
                dl_mt2blbl_signal  = Plot( 
                    name = "dl_mt2blbl_signal",
                    texX = 'M_{T2}(blbl) (GeV)', texY = 'Number of Events / 20 GeV' if args.normalizeBinWidth else "Number of Events",
                    binning=Binning.fromThresholds([0,20,40,60,80,100,120,140,160,200,250,300,350]),
                    stack = stack_signal,
                    attribute = TreeVariable.fromString( "dl_mt2blbl/F" ),
                    weight      = signal_weight )
                plots.append( dl_mt2blbl_signal )
    
        dl_mt2blbl_mc  = Plot(
            name = "dl_mt2blbl_mc",
            texX = 'M_{T2}(blbl) (GeV)', texY = 'Number of Events / 20 GeV' if args.normalizeBinWidth else "Number of Events",
            binning=Binning.fromThresholds([0,20,40,60,80,100,120,140,160,200,250,300,350]),
            stack = stack_mc,
            attribute = TreeVariable.fromString( selectionModifier("dl_mt2blbl/F") )      if selectionModifier is not None else None,
            selectionString = selectionModifier(cutInterpreter.cutString(args.selection)) if selectionModifier is not None else None,
            weight          = mc_weight )
        plots.append( dl_mt2blbl_mc )
    
    nBtagBinning = [6, 0, 6] 
    if args.variation == 'central':
        nbtags_data   = Plot(
            name        = "nbtags_data",
            texX        = 'number of b-tags (DeepCSV)', texY = 'Number of Events' if args.normalizeBinWidth else "Number of Events",
            binning     = nBtagBinning,
            stack       = stack_data,
            attribute   = TreeVariable.fromString( "nBTag/I" ),
            weight      = data_weight )
        plots.append( nbtags_data )
        if args.signal:
            nbtags_signal   = Plot(
                name        = "nbtags_signal",
                texX        = 'number of b-tags (DeepCSV)', texY = 'Number of Events' if args.normalizeBinWidth else "Number of Events",
                binning     = nBtagBinning,
                stack       = stack_signal,
                attribute   = TreeVariable.fromString( "nBTag/I" ),
                weight      = signal_weight )
            plots.append( nbtags_signal )

    nbtags_mc  = Plot(\
        name            = "nbtags_mc",
        texX            = 'number of b-tags (DeepCSV)', texY = 'Number of Events' if args.normalizeBinWidth else "Number of Events",
        binning         = nBtagBinning,
        stack           = stack_mc,
        attribute       = TreeVariable.fromString( selectionModifier("nBTag/I"))   if selectionModifier is not None else None,
        selectionString = selectionModifier(cutInterpreter.cutString(args.selection)) if selectionModifier is not None else None,
        weight          = mc_weight )
    plots.append( nbtags_mc )
    

    jetBinning = [8,2,10] if args.selection.count('njet2') else [2,0,2]
    if args.variation == 'central':
        njets_data   = Plot(
            name        = "njets_data",
            texX        = 'number of jets', texY = 'Number of Events' if args.normalizeBinWidth else "Number of Events",
            binning     = jetBinning,
            stack       = stack_data,
            attribute   = TreeVariable.fromString( "nJetGood/I" ),
            weight      = data_weight )
        plots.append( njets_data )
        if args.signal:
            njets_signal   = Plot(
                name        = "njets_signal",
                texX        = 'number of jets', texY = 'Number of Events' if args.normalizeBinWidth else "Number of Events",
                binning     = jetBinning,
                stack       = stack_signal,
                attribute   = TreeVariable.fromString( "nJetGood/I" ),
                weight      = signal_weight )
            plots.append( njets_signal )

    njets_mc  = Plot(\
        name            = "njets_mc",
        texX            = 'number of jets', texY = 'Number of Events' if args.normalizeBinWidth else "Number of Events",
        binning         = jetBinning,
        stack           = stack_mc,
        attribute       = TreeVariable.fromString( selectionModifier("nJetGood/I"))   if selectionModifier is not None else None,
        selectionString = selectionModifier(cutInterpreter.cutString(args.selection)) if selectionModifier is not None else None,
        weight          = mc_weight )
    plots.append( njets_mc )
    
    
    metBinning = [0,20,40,60,80] if args.selection.count('metInv') else [80,130,180,230,280,320,420,520,800] if args.selection.count('met80') else [0,80,130,180,230,280,320,420,520,800]
    if args.variation == 'central':
        met_data = Plot( 
            name        = "met_data",
            texX        = 'E_{T}^{miss} (GeV)', texY = 'Number of Events / 50 GeV' if args.normalizeBinWidth else "Number of Events",
            binning     = Binning.fromThresholds( metBinning ),
            stack       = stack_data, 
            attribute   = TreeVariable.fromString( "met_pt/F" ),
            weight      = data_weight,
            )
        plots.append( met_data )
        if args.signal:
            met_signal = Plot( 
                name        = "met_signal",
                texX        = 'E_{T}^{miss} (GeV)', texY = 'Number of Events / 50 GeV' if args.normalizeBinWidth else "Number of Events",
                binning     = Binning.fromThresholds( metBinning ),
                stack       = stack_signal, 
                attribute   = TreeVariable.fromString( "met_pt/F" ),
                weight      = signal_weight,
                )
            plots.append( met_signal )
  
    met_mc  = Plot(\
        name            = "met_pt_mc",
        texX            = 'E_{T}^{miss} (GeV)', texY = 'Number of Events / 50 GeV' if args.normalizeBinWidth else "Number of Events",
        binning         = Binning.fromThresholds( metBinning ),
        stack           = stack_mc,
        attribute       = TreeVariable.fromString( selectionModifier("met_pt/F") )      if selectionModifier is not None else None,
        selectionString = selectionModifier( cutInterpreter.cutString(args.selection) ) if selectionModifier is not None else None,
        weight          = mc_weight )
    plots.append( met_mc )
    

    metBinning2 = [0,20,40,60,80] if args.selection.count('metInv') else [80,100,120,140,160,200,500] if args.selection.count('met80') else [0,80,100,120,140,160,200,500]
    if args.variation == 'central':
        met2_data   = Plot(
            name        = "met2_data",
            texX        = 'E_{T}^{miss} (GeV)', texY = 'Number of Events / 20 GeV' if args.normalizeBinWidth else "Number of Events",
            binning     = Binning.fromThresholds(metBinning2),
            stack       = stack_data,
            attribute   = TreeVariable.fromString( "met_pt/F" ),
            weight      = data_weight )
        plots.append( met2_data )
        if args.signal:
            met2_signal   = Plot(
                name        = "met2_signal",
                texX        = 'E_{T}^{miss} (GeV)', texY = 'Number of Events / 20 GeV' if args.normalizeBinWidth else "Number of Events",
                binning     = Binning.fromThresholds(metBinning2),
                stack       = stack_signal,
                attribute   = TreeVariable.fromString( "met_pt/F" ),
                weight      = signal_weight )
            plots.append( met2_signal )

    met2_mc  = Plot(\
        name            = "met2_mc",
        texX            = 'E_{T}^{miss} (GeV)', texY = 'Number of Events / 20 GeV' if args.normalizeBinWidth else "Number of Events",
        binning         = Binning.fromThresholds(metBinning2),
        stack           = stack_mc,
        attribute       = TreeVariable.fromString( selectionModifier("met_pt/F"))   if selectionModifier is not None else None,
        selectionString = selectionModifier(cutInterpreter.cutString(args.selection)) if selectionModifier is not None else None,
        weight          = mc_weight )
    plots.append( met2_mc )


    metSigBinning = [0,2,4,6,8,10,12] if args.selection.count('POGMetSig0To12') else [12,16,20,24,28,32,36,40,45,50,55,60,65,70,75,80,85,90,95,100] if args.selection.count('POGMetSig12') else [0,4,8,12,16,20,24,28,32,36,40,45,50,55,60,65,70,75,80,85,90,95,100]
    if args.variation == 'central': 
        metSig_data  = Plot( 
            name        = "MET_significance_data",
            texX        = 'p_{T}^{miss} significance', texY = 'Number of Events' if args.normalizeBinWidth else "Number of Events",
            binning     = Binning.fromThresholds( metSigBinning ),
            stack       = stack_data, 
            attribute   = TreeVariable.fromString( "MET_significance/F" ),
            weight      = data_weight,
            )
        plots.append( metSig_data )
        if args.signal:
            metSig_signal  = Plot( 
                name        = "MET_significance_signal",
                texX        = 'p_{T}^{miss} significance', texY = 'Number of Events' if args.normalizeBinWidth else "Number of Events",
                binning     = Binning.fromThresholds( metSigBinning ),
                stack       = stack_signal, 
                attribute   = TreeVariable.fromString( "MET_significance/F" ),
                weight      = signal_weight,
                )
            plots.append( metSig_signal )
    
    metSig_mc  = Plot(\
        name = "MET_significance_mc",
        texX = 'p_{T}^{miss} significance', texY = 'Number of Events' if args.normalizeBinWidth else "Number of Events",
        stack = stack_mc,
        attribute = TreeVariable.fromString( selectionModifier("MET_significance/F") )  if selectionModifier is not None else None,
        binning=Binning.fromThresholds( metSigBinning ),
        selectionString = selectionModifier( cutInterpreter.cutString(args.selection) ) if selectionModifier is not None else None,
        weight          = mc_weight )
    plots.append( metSig_mc )

    if args.newMetSigPlots:
        metSigBinning400 = [0,2,4,6,8,10,12] if args.selection.count('POGMetSig0To12') else [12,32,52,72,92,112,132,200,300,400] if args.selection.count('POGMetSig12') else [0,12,32,52,72,92,112,132,166,200,250,300,350,400]
        if args.variation == 'central': 
            metSig_400_data  = Plot( 
                name        = "MET_significance_400_data",
                texX        = 'p_{T}^{miss} significance', texY = 'Number of Events' if args.normalizeBinWidth else "Number of Events",
                binning     = Binning.fromThresholds( metSigBinning400 ),
                stack       = stack_data, 
                attribute   = TreeVariable.fromString( "MET_significance/F" ),
                weight      = data_weight,
                )
            plots.append( metSig_400_data )
            if args.signal:
                metSig_400_signal  = Plot( 
                    name        = "MET_significance_400_signal",
                    texX        = 'p_{T}^{miss} significance', texY = 'Number of Events' if args.normalizeBinWidth else "Number of Events",
                    binning     = Binning.fromThresholds( metSigBinning400 ),
                    stack       = stack_signal, 
                    attribute   = TreeVariable.fromString( "MET_significance/F" ),
                    weight      = signal_weight,
                    )
                plots.append( metSig_400_signal )
        
        metSig_400_mc  = Plot(\
            name = "MET_significance_400_mc",
            texX = 'p_{T}^{miss} significance', texY = 'Number of Events' if args.normalizeBinWidth else "Number of Events",
            stack = stack_mc,
            attribute = TreeVariable.fromString( selectionModifier("MET_significance/F") )  if selectionModifier is not None else None,
            binning=Binning.fromThresholds( metSigBinning400 ),
            selectionString = selectionModifier( cutInterpreter.cutString(args.selection) ) if selectionModifier is not None else None,
            weight          = mc_weight )
        plots.append( metSig_400_mc )

        #sqrt MET_significance
        sqrtMetSigBinning = [3.4+i*1 for i in range(27)]
        if args.variation == 'central': 
            sqrtMetSig_data  = Plot( 
                name        = "sqrt_MET_significance_data",
                texX        = '#sqrt{p_{T}^{miss} significance}', texY = 'Number of Events' if args.normalizeBinWidth else "Number of Events",
                binning     = Binning.fromThresholds( sqrtMetSigBinning ),
                stack       = stack_data, 
                attribute   = lambda event, sample: sqrt(getattr(event, "MET_significance")),
                weight      = data_weight,
                )
            plots.append( sqrtMetSig_data )
            if args.signal:
                sqrtMetSig_signal  = Plot( 
                    name        = "sqrt_MET_significance_signal",
                    texX        = '#sqrt{p_{T}^{miss} significance}', texY = 'Number of Events' if args.normalizeBinWidth else "Number of Events",
                    binning     = Binning.fromThresholds( sqrtMetSigBinning ),
                    stack       = stack_signal, 
                    attribute   = lambda event, sample: sqrt(getattr(event, "MET_significance")),
                    weight      = signal_weight,
                    )
                plots.append( sqrtMetSig_signal )
        
        sqrtMetSig_mc  = Plot(\
            name = "sqrt_MET_significance_mc",
            texX = '#sqrt{p_{T}^{miss} significance}', texY = 'Number of Events' if args.normalizeBinWidth else "Number of Events",
            stack = stack_mc,
            #attribute = lambda event, sample: sqrt(getattr(event, selectionModifier("MET_significance")))  if selectionModifier is not None else None,
            attribute = TreeVariable.fromString( selectionModifier("MET_significance/F") )  if selectionModifier is not None else None,
            binning=Binning.fromThresholds( sqrtMetSigBinning ),
            selectionString = selectionModifier( cutInterpreter.cutString(args.selection) ) if selectionModifier is not None else None,
            weight          = mc_weight )
        plots.append( sqrtMetSig_mc )


#    if args.variation == 'central': 
#        jetRelIsoRecorr_data  = Plot( 
#            name        = "l1_jetRelIsoRecorrHad_data",
#            texX        = 'relIso(l_{1}) (GeV) with recorrected jet pt', texY = "Number of Events",
#            binning     = [50,0,.5], 
#            stack       = stack_data, 
#            attribute   = lambda event, sample: event.l1_jetRelIsoRecorrHad, 
#            weight      = data_weight,
#            )
#        plots.append( jetRelIsoRecorr_data )
#    
#    jetRelIsoRecorr_mc  = Plot(\
#        name = "l1_jetRelIsoRecorrHad_mc",
#        texX = 'relIso(l_{1}) (GeV) with recorrected jet pt', texY = "Number of Events",
#        stack = stack_mc,
#        binning     = [50,0,.5], 
#        attribute   = lambda event, sample: event.l1_jetRelIsoRecorrHad, 
#        #attribute   = lambda event, sample: getattr(event, selectionModifier("l1_jetRelIsoRecorrHad/F") if selectionModifier is not None else None, 
#        selectionString = selectionModifier( cutInterpreter.cutString(args.selection) ) if selectionModifier is not None else None,
#        weight          = mc_weight )
#    plots.append( jetRelIsoRecorr_mc )

    ## Make plot directory
    #plot_directory_ = os.path.join(plot_directory, 'systematicPlots', args.plot_directory, args.selection, args.era, mode)
    #try: 
    #    os.makedirs(plot_directory_)
    #except: 
    #    pass

    if args.variation is not None:
        key  = (args.era, mode, args.variation)

        success = False
        if dirDB.contains(key) and not args.overwrite:
            normalisation_mc, normalisation_data, histos = dirDB.get( key )
            #FIXME: remove next line
            #normalisation_mc = {s.name: 0 for s in mc}
            for i_p, h_s in enumerate(histos):
                plots[i_p].histos = h_s
            logger.info( "Loaded normalisations and histograms for %s in mode %s from cache.", args.era, mode)
            logger.debug("Loaded normalisation_mc: %r normalisation_data: %r", normalisation_mc, normalisation_data )
            if normalisation_mc['Top_pow']<=0:
                if args.scaling is not None:
                    success = False
                    logger.info( "Top_pow histo is zero! Maybe the '-mt2llTo100' selection interferes with yours" )
                else:
                    success = True
            else: 
                success = True
        if not success:
            logger.info( "Obtain normalisations and histograms for %s in mode %s.", args.era, mode)
            # Calculate the normalisation yield for mt2ll<100
            normalization_selection_string = selectionModifier(cutInterpreter.cutString(args.selection + '-mt2llTo100')) 
#NN ------------------------------------------------------------------------------------------------------------
####            mc_normalization_weight_string    = MC_WEIGHT(variations[args.variation], returntype='string')
            if args.variation is not None and "Input" in args.variation: # input variation
                mc_normalization_weight_string    = MC_WEIGHT(variations[args.variation], returntype='inputString')
            else:
                mc_normalization_weight_string    = MC_WEIGHT(variations[args.variation], returntype='string')
# --------------------------------------------------------------------------------------------------------------
            if args.scaling is not None:
                normalisation_mc = {s.name :s.scale*s.getYieldFromDraw(selectionString = normalization_selection_string, weightString = mc_normalization_weight_string)['val'] for s in mc}
                print "~~~ normalization mc: ", normalisation_mc
            else:
                normalisation_mc = {s.name: 0 for s in mc}
            #print normalization_selection_string, mc_normalization_weight_string

            if args.variation == 'central':
                if args.scaling is not None:
                    normalisation_data = data_sample.scale*data_sample.getYieldFromDraw( selectionString = normalization_selection_string, weightString = data_weight_string)['val']
                else:
                    normalisation_data = -5
            else:
                normalisation_data = -1

            logger.info( "Making plots.")
            plotting.fill(plots, read_variables = read_variables, sequence = sequence)

            # Delete lambda because we can't serialize it
            for plot in plots:
                del plot.weight

            # save
            #print "normalisation_mc %f"%(normalisation_mc)
            dirDB.add( key, (normalisation_mc, normalisation_data, [plot.histos for plot in plots]), overwrite = args.overwrite)

            logger.info( "Done with %s in channel %s.", args.variation, mode)

if args.variation is not None:
    logger.info( "Done with modes %s and variation %s of selection %s. Quit now.", ",".join( modes ), args.variation, args.selection )
    sys.exit(0)

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#                               COMBINE CACHED HISTOS
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# Systematic pairs:( 'name', 'up', 'down' )
systematics = [\
    {'name':'JEC',          'pair':('jesTotalUp', 'jesTotalDown')},
    {'name':'Unclustered',  'pair':('unclustEnUp', 'unclustEnDown')},
    {'name':'PU',           'pair':('PUUp', 'PUDown')},
    {'name':'BTag_b',       'pair':('BTag_SF_b_Down', 'BTag_SF_b_Up' )},
    {'name':'BTag_l',       'pair':('BTag_SF_l_Down', 'BTag_SF_l_Up')},
    {'name':'trigger',      'pair':('DilepTriggerDown', 'DilepTriggerUp')},
    {'name':'leptonSF',     'pair':('LeptonSFDown', 'LeptonSFUp')},
    {'name':'leptonHit0SF', 'pair':('LeptonHit0SFDown', 'LeptonHit0SFUp')},
    {'name':'leptonSip3dSF','pair':('LeptonSip3dSFDown', 'LeptonSip3dSFUp')},
    #{'name': 'TopPt',      'pair':(  'TopPt', 'central')},
    {'name': 'JER',         'pair':('jerUp', 'jerDown')},
    {'name': 'L1Prefire',   'pair':('L1PrefireUp', 'L1PrefireDown')},
#NN
#    {'name': 'DYInput',           'pair':('DYInputUp', 'DYInputDown')},
#    {'name': 'TTInput',           'pair':('TTInputUp', 'TTInputDown')},
#    {'name': 'MBInput',           'pair':('MBInputUp', 'MBInputDown')},
#    {'name': 'TTZInput',          'pair':('TTZInputUp', 'TTZInputDown')},
#    {'name': 'OtherInput',        'pair':('OtherInputUp', 'OtherInputDown')},
#
#    {'name': 'TT1JetMismInput',   'pair':('TT1JetMismUp', 'TT1JetMismDown')},
#    {'name': 'TTTotJetMismInput', 'pair':('TTTotJetMismUp', 'TTTotJetMismDown')},
#    {'name': 'TTNonPromptInput',  'pair':('TTNonPromptInputUp', 'TTNonPromptInputDown')},
]

# loop over modes
missing_cmds   = []
variation_data = {}
for mode in modes:
    logger.info('Working on mode: %s', mode)
    logger.info('Now attempting to load all variations from dirDB %s', dirDB.directory)
   
    for variation in variations.keys():
        key  = (args.era, mode, variation)
        success = False
        if dirDB.contains(key) and not args.overwrite:
            normalisation_mc, normalisation_data, histos = dirDB.get(key)
            #FIXME: remove next line
            #normalisation_mc = {s.name: 0 for s in mc}
            variation_data[(mode, variation)] = {'histos':histos, 'normalisation_mc':normalisation_mc, 'normalisation_data':normalisation_data}
            logger.info( "Loaded normalisations and histograms for variation %s, era %s in mode %s from cache.", variation, args.era, mode)
            if normalisation_mc['Top_pow']<=0:
                if args.scaling is not None:
                    success = False
                    logger.info( "Top_pow histo is zero! Maybe the '-mt2llTo100' selection for scaling interferes with yours" )
                else:
                    success = True
            else: 
                success = True
        if not success:
            # prepare sub variation command
            cmd = ['python', 'systematicVariation.py']
            cmd.append('--logLevel=%s'%args.logLevel)
            if args.signal is not None: cmd.append( '--signal=%s'%args.signal )
            cmd.append('--era=%s'%args.era)
            cmd.append('--plot_directory=%s'%args.plot_directory)
            cmd.append('--reweightPU=%s'%args.reweightPU)
            cmd.append('--selection=%s'%args.selection)
            cmd.append('--mode=%s'%args.mode)
            cmd.append('--variation=%s'%variation)
            if args.scaling is not None: cmd.append('--scaling=%s'%args.scaling)
            if args.normalizeBinWidth: cmd.append('--normalizeBinWidth')
            if args.small: cmd.append('--small')
            if args.newMetSigPlots: cmd.append('--newMetSigPlots')
            if args.dpm: cmd.append('--dpm')
            if args.scratch: cmd.append('--scratch')
            if args.overwrite: cmd.append('--overwrite')

            cmd_string = ' '.join( cmd )
            missing_cmds.append( cmd_string )
            logger.info("Missing variation %s, era %s in mode %s in cache. Need to run: \n%s", variation, args.era, mode, cmd_string)

# write missing cmds
filename = 'missing.sh'
if os.path.exists(filename) and args.add:
    append_write = 'a' # append if already exists
else:
    append_write = 'w' # make a new file if not
missing_cmds = list(set(missing_cmds))
if len(missing_cmds)>0:
    with file( filename, append_write ) as f:
        # if we start the file:
        if append_write == 'w':
            f.write("#!/bin/sh\n")
        for cmd in missing_cmds:
            f.write( cmd + '\n')
    logger.info( "Written %i variation commands to ./missing.sh. Now I quit!", len(missing_cmds) )
    sys.exit(0)

# make 'all' and 'SF' from ee/mumu/mue
new_modes = []
all_modes = list(modes)
if 'mumu' in modes and 'ee' in modes:
    new_modes.append( ('SF', ('mumu', 'ee')) )
    all_modes.append( 'SF' )
if 'mumu' in modes and 'ee' in modes and 'mue' in modes:
    new_modes.append( ('all', ('mue', 'SF')) )
    all_modes.append( 'all' )
for variation in variations:
    for new_mode, old_modes in new_modes:
        new_key = ( new_mode, variation )
        variation_data[new_key] = {}
        # Adding up data_normalisation 
        if variation == 'central':
            variation_data[new_key]['normalisation_data'] = sum( variation_data[( old_mode, variation )]['normalisation_data'] for old_mode in old_modes )
        else:
            variation_data[new_key]['normalisation_data'] = -1 

        # Adding up mc normalisation
        sample_keys = variation_data[( old_modes[0], variation )]['normalisation_mc'].keys()
        variation_data[new_key]['normalisation_mc'] = {}
        for sample_key in sample_keys: 
            variation_data[new_key]['normalisation_mc'][sample_key] = variation_data[( old_modes[0], variation )]['normalisation_mc'][sample_key]
            for mode in old_modes[1:]:
                variation_data[new_key]['normalisation_mc'][sample_key] += variation_data[( mode, variation )]['normalisation_mc'][sample_key]

        # Adding up histos (clone old_modes[0] at 3rd level, then add)
        variation_data[new_key]['histos'] = [[[ h.Clone() for h in hs ] for hs in plot_histos ] for plot_histos in variation_data[( old_modes[0], variation )]['histos']]
        for mode in old_modes[1:]:
            for i_plot_histos, plot_histos in  enumerate(variation_data[( mode, variation )]['histos']):
                for i_hs, hs in enumerate(plot_histos):
                    for i_h, h in enumerate(hs):
                        variation_data[new_key]['histos'][i_plot_histos][i_hs][i_h].Add(h)

# SF for top central such that we get area normalisation 
dataMC_SF = {}
for mode in all_modes:
    # All SF to 1
    dataMC_SF[mode] = {variation:{s.name:1 for s in mc} for variation in variations} 
    yield_data = variation_data[(mode,'central')]['normalisation_data'] 
    if args.scaling == 'top': 
        # scale variations individually
        if args.variation_scaling:
            logger.info( "Scaling top yield to data for mt2ll<100 individually for all variations." )
            for variation in variations.keys():
                #print ""%()
                yield_non_top = sum( val for name, val in variation_data[(mode,variation)]['normalisation_mc'].iteritems() if name != Top_pow.name)
                yield_top     = variation_data[(mode,variation)]['normalisation_mc'][Top_pow.name]
                #print "mode %s yield_data %f yield_non_top %f yield_top %f"%(mode, yield_data, yield_non_top, yield_top)
                dataMC_SF[mode][variation][Top_pow.name] = (yield_data - yield_non_top)/yield_top
                #if mode=='mumu' and variation=='central': assert False, ''
        # scale all variations with the central factor
        else:
            logger.info( "Scaling top yield to data for mt2ll<100 ( all variations are scaled by central SF)" )
            yield_non_top = sum( val for name, val in variation_data[(mode,'central')]['normalisation_mc'].iteritems() if name != Top_pow.name)
            yield_top     = variation_data[(mode,'central')]['normalisation_mc'][Top_pow.name]
            sf = (yield_data - yield_non_top)/yield_top
            for variation in variations.keys():
                dataMC_SF[mode][variation][Top_pow.name] = sf 
    elif args.scaling == 'mc':
        # scale variations individually
        if args.variation_scaling:
            logger.info( "Scaling MC yield to data for mt2ll<100 individually for all variations." )
            for variation in variations.keys():
                yield_mc = sum( val for name, val in variation_data[(mode,variation)]['normalisation_mc'].iteritems())
                for s in mc:
                    dataMC_SF[mode][variation][s.name] = yield_data/yield_mc
        # scale all variations with the central factor
        else:
            logger.info( "Scaling MC yield to data for mt2ll<100 ( all variations are scaled by central SF)" )
            yield_mc = sum( val for name, val in variation_data[(mode,'central')]['normalisation_mc'].iteritems())
            sf = yield_data/yield_mc
            for variation in variations.keys():
                for s in mc:
                    dataMC_SF[mode][variation][s.name] = sf 

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # area normalization scale factor
    if args.scaling is None and args.normalize:
        #assert False, "Shouldn't be here!"
        # scale variations individually
        logger.info( "Scaling MC yield to data ( all variations are scaled by central SF)" )
        for i_plot, plot in enumerate(plots):
            if plot.name=="sqrt_MET_significance_mc":
                continue
            for variation in variations.keys():
                if args.signal:
                    yield_data = variation_data[(mode,'central')]['histos'][3*i_plot][0][0].Integral()
                    yield_mc = sum(variation_data[(mode,'central')]['histos'][3*i_plot+2][0][i].Integral() for i, s in enumerate(mc))
                else:
                    yield_data = variation_data[(mode,'central')]['histos'][2*i_plot][0][0].Integral()
                    yield_mc = sum(variation_data[(mode,'central')]['histos'][2*i_plot+1][0][i].Integral() for i, s in enumerate(mc))
                sf = yield_data/yield_mc
                #if variation == "central" and plot.name == "MET_significance_400_mc": 
                print "{}: {} / {} = {}".format(plot.name, yield_data, yield_mc, sf)
                for s in mc:
                    if variation == "central" and plot.name == "MET_significance_400_mc": print "central sf: {}".format(sf)
                    dataMC_SF[mode][variation][s.name] = sf
        print "{}: {}".format(mode, dataMC_SF[mode])
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    # perform the scaling & store styles and texName
    for i_plot, plot in enumerate(plots):
        for variation in variations.keys():
            for s in mc:
                pos_plot = i_plot if variation != 'central' else ( 3*i_plot+2 if args.signal else 2*i_plot+1 )
                variation_data[(mode, variation)]['histos'][pos_plot][0][position[s.name]].Scale( dataMC_SF[mode][variation][s.name] ) 
                variation_data[(mode, variation)]['histos'][pos_plot][0][position[s.name]].style = s.style 
                variation_data[(mode, variation)]['histos'][pos_plot][0][position[s.name]].legendText = s.texName
            for i_s, s in enumerate(signals):
                if variation == 'central' and args.signal: 
                    variation_data[(mode, 'central')]['histos'][3*i_plot+1][i_s][0].style = signals[i_s].style 
                    variation_data[(mode, 'central')]['histos'][3*i_plot+1][i_s][0].legendText = signals[i_s].texName
            if variation == 'central': 
                pos_plot =  3*i_plot if args.signal else 2*i_plot 
                variation_data[(mode, 'central')]['histos'][pos_plot][0][0].style = data_sample.style 
                variation_data[(mode, 'central')]['histos'][pos_plot][0][0].legendText = data_sample.texName
                # debugging:
                yield_mc = sum(variation_data[(mode,'central')]['histos'][3*i_plot+2][0][i].Integral() for i, s in enumerate(mc))
                scale_mc = [dataMC_SF[mode][variation][s.name] for s in mc]
                if plot.name == "MET_significance_400_mc": print "{}: {} scale: {}".format(plot.name, yield_mc, scale_mc)

# store everything in the dir_db
dirdb_key =   'variation_data_scaling_%s'%(args.scaling if args.scaling is not None else "None")
dirdb_key += "_variation_scaling_%s"%bool(args.variation_scaling)
dirdb_key += "_normalize_%s"%bool(args.normalize)
print "DB directory:", os.path.join(plot_directory, 'systematicPlots', args.era, plot_subdirectory, args.selection, 'cache')
print "dirdb_key:", dirdb_key
save_plots = []
for plot in plots:
    p = copy.deepcopy(plot)
    p.stack = None
    save_plots.append( p )
stack_mc_  = copy.deepcopy(stack_mc)
stack_data_  = copy.deepcopy(stack_data)
if args.signal:
    stack_signal_  = copy.deepcopy(stack_signal)
else:
    stack_signal_=Stack([])
for stack in [stack_mc_, stack_data_, stack_signal_]:
    for _s in stack:
        for s in _s:
            del s.style 
#FIXME: add back
dirDB.add(dirdb_key, ( variation_data, save_plots, stack_mc_, stack_data_, stack_signal_), overwrite = True)

#pickle.dump( ( variation_data, save_plots), file(os.path.join(plot_directory, 'systematicPlots', args.era, plot_subdirectory, args.selection, 'variation_data.pkl'), 'w' ) )


def drawObjects( scaling ):
    tex = ROOT.TLatex()
    tex.SetNDC()
    tex.SetTextSize(0.04)
    tex.SetTextAlign(11) # align right
    lines = [
      (0.15, 0.95, 'CMS Private') if args.private else (0.15, 0.95, 'CMS Preliminary'),
      ]
    #if scaling == 'mc':
    #  lines += [(0.45, 0.95, 'L=%3.1f fb{}^{-1} (13 TeV) SF(mc)=%3.2f'% ( lumi_scale, scaleFactor ) )]
    #elif scaling == 'top':
    #  lines += [(0.45, 0.95, 'L=%3.1f fb{}^{-1} (13 TeV) SF(top)=%3.2f'% ( lumi_scale, scaleFactor ) )]
    #elif scaling is None and args.normalize:
    #  lines += [(0.45, 0.95, 'L=%3.1f fb{}^{-1} (13 TeV) scale=%3.2f'% ( lumi_scale, scaleFactor ) )]
    #else:
    #  lines += [(0.45, 0.95, 'L=%3.1f fb{}^{-1} (13 TeV)'% ( lumi_scale) )]
    lines += [(0.45, 0.95, 'L=%3.1f fb{}^{-1} (13 TeV)'% ( lumi_scale) )]
    if "mt2ll100" in args.selection: 
        if args.signal:
            lines += [(0.55, 0.6, 'M_{T2}(ll) > 100 GeV')] # Manually put the mt2ll > 100 GeV label
        else:
            lines += [(0.55, 0.65, 'M_{T2}(ll) > 100 GeV')] # Manually put the mt2ll > 100 GeV label
    return [tex.DrawLatex(*l) for l in lines]

# We plot now. 
if args.normalize: plot_subdirectory += "_normalized"
if args.beta:      plot_subdirectory += "_%s"%args.beta
if args.private:   plot_subdirectory += "_private"
for mode in all_modes:
    for i_plot, plot in enumerate(plots):
        
        # for central (=no variation), we store plot_data_1, plot_mc_1, plot_data_2, plot_mc_2, ...
        if args.signal:
            data_histo_list     = variation_data[(mode, 'central')]['histos'][3*i_plot]
            signal_histo_list   = variation_data[(mode, 'central')]['histos'][3*i_plot+1]
            mc_histo_list       = {'central': variation_data[(mode, 'central')]['histos'][3*i_plot+2] }
        else:
            data_histo_list = variation_data[(mode, 'central')]['histos'][2*i_plot]
            mc_histo_list   = {'central': variation_data[(mode, 'central')]['histos'][2*i_plot+1] }
        # for the other variations, there is no data
        for variation in variations.keys():
            if variation=='central': continue
            mc_histo_list[variation] = variation_data[(mode, variation)]['histos'][i_plot]

        # Add histos, del the stack (which refers to MC only )
        plot.histos =  mc_histo_list['central'] + data_histo_list
        plot.stack  = Stack( mc, [data_sample] )
        if args.signal != None: 
            plot.histos += signal_histo_list
            plot.stack.extend( [ [s] for s in signals ] ) 
        
        # Make boxes and ratio boxes
        boxes           = []
        ratio_boxes     = []
        # Compute all variied MC sums
        total_mc_histo   = {variation:add_histos( mc_histo_list[variation][0]) for variation in variations.keys() }

        # loop over bins & compute shaded uncertainty boxes
        boxes   = []
        r_boxes = []
        for i_b in range(1, 1 + total_mc_histo['central'].GetNbinsX() ):
            # Only positive yields
            total_central_mc_yield = total_mc_histo['central'].GetBinContent(i_b)
            # include overflow bin for the last bin
            overflowBin = True
            if i_b==total_mc_histo['central'].GetNbinsX() and overflowBin:
                total_central_mc_yield += total_mc_histo['central'].GetBinContent(i_b+1)
            if total_central_mc_yield<=0: continue
            variance = 0.
            for systematic in systematics:
                # Use 'central-variation' (factor 1) and 0.5*(varUp-varDown)
                if 'central' in systematic['pair']: 
                    factor = 1
                else:
                    factor = 0.5
                # sum in quadrature
                if i_b==total_mc_histo['central'].GetNbinsX() and overflowBin: # add overflow bin
                    up = total_mc_histo[systematic['pair'][0]].GetBinContent(i_b) + total_mc_histo[systematic['pair'][0]].GetBinContent(i_b+1)
                    down = total_mc_histo[systematic['pair'][1]].GetBinContent(i_b) + total_mc_histo[systematic['pair'][1]].GetBinContent(i_b+1)
                    if isnan(up) or isnan(down):
                        print "ISNAN", "up", up, "down", down
                        variance += ( factor*(total_mc_histo[systematic['pair'][0]].GetBinContent(i_b) - total_mc_histo[systematic['pair'][1]].GetBinContent(i_b)) )**2
                    else: 
                        variance += ( factor*(up - down) )**2
                else:
                    variance += ( factor*(total_mc_histo[systematic['pair'][0]].GetBinContent(i_b) - total_mc_histo[systematic['pair'][1]].GetBinContent(i_b)) )**2

            sigma     = sqrt(variance)
            sigma_rel = sigma/total_central_mc_yield 

            box = ROOT.TBox( 
                    total_mc_histo['central'].GetXaxis().GetBinLowEdge(i_b),
                    max([0.03, (1-sigma_rel)*total_central_mc_yield]),
                    total_mc_histo['central'].GetXaxis().GetBinUpEdge(i_b), 
                    max([0.03, (1+sigma_rel)*total_central_mc_yield]) )
            box.SetLineColor(ROOT.kBlack)
            box.SetFillStyle(3444)
            box.SetFillColor(ROOT.kBlack)
            boxes.append(box)
            
            r_box = ROOT.TBox( 
                total_mc_histo['central'].GetXaxis().GetBinLowEdge(i_b),  
                max(0.1, 1-sigma_rel), 
                total_mc_histo['central'].GetXaxis().GetBinUpEdge(i_b), 
                min(1.9, 1+sigma_rel) )
            r_box.SetLineColor(ROOT.kBlack)
            r_box.SetFillStyle(3444)
            r_box.SetFillColor(ROOT.kBlack)
            ratio_boxes.append(r_box)

        for log in [False, True]:
            plot_directory_ = os.path.join(plot_directory, 'systematicPlots', args.era, plot_subdirectory, args.selection, mode + ("_log" if log else ""))
            #if not max(l[0].GetMaximum() for l in plot.histos): continue # Empty plot
            texMode = "#mu#mu" if mode == "mumu" else "#mue" if mode == "mue" else mode
            if    mode == "all": plot.histos[1][0].legendText = "data (%s)"%args.era
            else:                plot.histos[1][0].legendText = "data (%s, %s)"%(args.era, texMode)

            _drawObjects = []

            plotting.draw(plot,
              plot_directory = plot_directory_,
              ratio = {'yRange':(0.1,1.9), 'drawObjects':ratio_boxes},
              logX = False, logY = log, sorting = False,
              yRange = (0.03, "auto") if log else (0.001, "auto"),
              scaling = {0:1} if args.normalize else {},
              legend = ( (0.18,0.88-0.03*sum(map(len, plot.histos)),0.9,0.88), 2),
              drawObjects = drawObjects( args.scaling ) + boxes,
              copyIndexPHP = True, extensions = ["png", "pdf", "root"],
            )
