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
import pickle
import os, time, sys, copy
from math                                import sqrt, cos, sin, pi, atan2, isnan

# RootTools
from RootTools.core.standard             import *

#Analysis / StopsDilepton / Samples
from StopsDilepton.tools.user            import plot_directory, analysis_results
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
argParser.add_argument('--plot_directory',    action='store',      default='v1')
argParser.add_argument('--selection',         action='store',      default='lepSel-POGMetSig12-njet2p-btag1p-miniIso0.2-looseLeptonMiniIsoVeto-mll20-dPhiJet0-dPhiJet1')
# Scalings
argParser.add_argument('--scaling',           action='store',      default=None, choices = [None, 'mc', 'top'],     help='Scale top to data in mt2ll<100?')
argParser.add_argument('--variation_scaling', action='store_true', help='Scale the variations individually to mimick bkg estimation?')
argParser.add_argument('--normalize',         action='store_true', help='Perform area normalization mc to data?')
argParser.add_argument('--beta',              action='store',      default=None, help="Add an additional directory label for minor changes to the plots")
argParser.add_argument('--small',             action='store_true',     help='Run only on a small subset of the data?')
argParser.add_argument('--private',           action='store_true',     help='Produce private plots')
argParser.add_argument('--directories',       action='store',         nargs='*',  type=str, default=[],                  help="Input" )
argParser.add_argument('--fitChi2',           action='store_true', help='')
argParser.add_argument('--SFUnc',             action='store',      default=None, help="SF pkl file, e.g 'SF', searched for in user results directory")


args = argParser.parse_args()

# Logger
import StopsDilepton.tools.logger as logger
import RootTools.core.logger as logger_rt
import Analysis.Tools.logger as logger_an
logger    = logger.get_logger(   args.logLevel, logFile = None)
logger_rt = logger_rt.get_logger(args.logLevel, logFile = None)
logger_an = logger_an.get_logger(args.logLevel, logFile = None)

# Define all systematic variations
variations = {
    'central'           : {},
    'jesTotalUp'        : {}, 
    'jesTotalDown'      : {}, 
    'jerUp'             : {}, 
    'jerDown'           : {}, 
    'unclustEnUp'       : {}, 
    'unclustEnDown'     : {}, 
    'PUUp'              : {}, 
    'PUDown'            : {}, 
    'BTag_SF_b_Down'    : {}, 
    'BTag_SF_b_Up'      : {}, 
    'BTag_SF_l_Down'    : {}, 
    'BTag_SF_l_Up'      : {}, 
    'DilepTriggerDown'  : {}, 
    'DilepTriggerUp'    : {}, 
    'LeptonSFDown'      : {}, 
    'LeptonSFUp'        : {}, 
    'LeptonHit0SFDown'  : {}, 
    'LeptonHit0SFUp'    : {}, 
    'LeptonSip3dSFDown' : {},
    'LeptonSip3dSFUp'   : {}, 
    'L1PrefireDown'     : {}, 
    'L1PrefireUp'       : {}, 
}

# Systematic pairs:( 'name', 'up', 'down' )
systematics = [\
    {'name':'JEC',          'correlated':True, 'pair':('jesTotalUp', 'jesTotalDown')},
    {'name':'Unclustered',  'correlated':True, 'pair':('unclustEnUp', 'unclustEnDown')},
    {'name':'PU',           'correlated':True, 'pair':('PUUp', 'PUDown')},
    {'name':'BTag_b',       'correlated':True, 'pair':('BTag_SF_b_Down', 'BTag_SF_b_Up' )},
    {'name':'BTag_l',       'correlated':True, 'pair':('BTag_SF_l_Down', 'BTag_SF_l_Up')},
    {'name':'trigger',      'correlated':True, 'pair':('DilepTriggerDown', 'DilepTriggerUp')},
    {'name':'leptonSF',     'correlated':True, 'pair':('LeptonSFDown', 'LeptonSFUp')},
    {'name':'leptonHit0SF', 'correlated':True, 'pair':('LeptonHit0SFDown', 'LeptonHit0SFUp')},
    {'name':'leptonSip3dSF','correlated':True, 'pair':('LeptonSip3dSFDown', 'LeptonSip3dSFUp')},
    #{'name': 'TopPt',      'correlated':True, 'pair':(  'TopPt', 'central')},
    {'name': 'JER',         'correlated':True, 'pair':('jerUp', 'jerDown')},
    {'name': 'L1Prefire',   'correlated':True, 'pair':('L1PrefireUp', 'L1PrefireDown')},
]


# store everything in the dir_db
lumi_scale = 137
#args.directories=\
#    [
#        "/afs/hephy.at/user/r/rschoefbeck/www/StopsDilepton/systematicPlots/Run2016/v0p19_small_T2tt_reweightPUCentral/lepSel-POGMetSig12-njet2p-btag1p-miniIso0.2-looseLeptonMiniIsoVeto-mll20-dPhiJet0-dPhiJet1",
#        "/afs/hephy.at/user/r/rschoefbeck/www/StopsDilepton/systematicPlots/Run2016/v0p19_small_T2tt_reweightPUCentral/lepSel-POGMetSig12-njet2p-btag1p-miniIso0.2-looseLeptonMiniIsoVeto-mll20-dPhiJet0-dPhiJet1",
#        "/afs/hephy.at/user/r/rschoefbeck/www/StopsDilepton/systematicPlots/Run2016/v0p19_small_T2tt_reweightPUCentral/lepSel-POGMetSig12-njet2p-btag1p-miniIso0.2-looseLeptonMiniIsoVeto-mll20-dPhiJet0-dPhiJet1",
#        "/afs/hephy.at/user/r/rschoefbeck/www/StopsDilepton/systematicPlots/Run2017/v0p19_small_T2tt_reweightPUCentral/lepSel-badEEJetVeto-POGMetSig12-njet2p-btag1p-miniIso0.2-looseLeptonMiniIsoVeto-mll20-dPhiJet0-dPhiJet1",
#        "/afs/hephy.at/user/r/rschoefbeck/www/StopsDilepton/systematicPlots/Run2018/v0p19_small_T2tt_reweightPUVUp/lepSel-POGMetSig12-njet2p-btag1p-miniIso0.2-looseLeptonMiniIsoVeto-mll20-dPhiJet0-dPhiJet1",
#    ]

dirdb_key =   'variation_data_scaling_%s'%(args.scaling if args.scaling is not None else "None")
dirdb_key += "_variation_scaling_%s"%bool(args.variation_scaling)
dirdb_key += "_normalize_%s"%bool(args.normalize)

variation_data = {}
plots          = {}
for i_directory, directory in enumerate(args.directories):

    dirDB = DirDB(os.path.join(directory, 'cache'))
    try:
        v_data, v_plots, stack_mc, stack_data, stack_signal = dirDB.get(dirdb_key)
    except TypeError as e:
        logger.error("Something wrong with directory %i %s key: %s", i_directory, directory, dirdb_key)
        raise e
    variation_data[directory] = v_data
    plots[directory] = v_plots 
    if variation_data[directory] is None: 
        logger.error( "Did not find %s in %s", dirdb_key, directory )
    else:
        logger.info( "Found %s in %s", dirdb_key, directory )
# just take the first one for the plots. They must all be the same
plots = plots.values()[0]     

def drawObjects( scaling ):
    tex = ROOT.TLatex()
    tex.SetNDC()
    tex.SetTextSize(0.04)
    tex.SetTextAlign(11) # align right
    lines = [
      (0.15, 0.95, 'CMS Private') if args.private else (0.15, 0.95, 'CMS #it{Preliminary}'),
      ]
    lines += [(0.45, 0.95, 'L=%i fb^{-1} (13 TeV)'% ( int(lumi_scale) ) )]
    if "mt2ll100" in args.selection:
        if args.signal:
            lines += [(0.55, 0.58, 'M_{T2}(ll) > 100 GeV')] # Manually put the mt2ll > 100 GeV label
        else:
            lines += [(0.55, 0.65, 'M_{T2}(ll) > 100 GeV')] # Manually put the mt2ll > 100 GeV label
    return [tex.DrawLatex(*l) for l in lines]

def accumulate_level_2( lists_of_lists ):
    #res = copy.deepcopy( lists_of_lists[0] )
    res = [ map(lambda h:h.Clone(), l) for l in lists_of_lists[0] ]
    for stack_histos in lists_of_lists[1:]:
        for i_histos, histos in enumerate(stack_histos):
            for i_histo,histo in enumerate(histos):
                res[i_histos][i_histo].Add( histo )
    return res

# arguments & directory
plot_subdirectory = args.plot_directory

# Loading & scaling + SF uncertainties
if args.SFUnc is not None:
    SFUnc = pickle.load( file(os.path.join( analysis_results, args.SFUnc )) )
    # restructure 
    for key in SFUnc.keys():
        mc_key, year = key.split("_")
        if not SFUnc.has_key( int(year) ):
            SFUnc[int(year)] = {}
        SFUnc[int(year)][mc_key] = SFUnc[key] 
        systematics.append( {'name': 'ScaleFactor'+mc_key,   'correlated':True, 'pair':('ScaleFactor'+mc_key+'Up', 'ScaleFactor'+mc_key+'Down')} )
    daniel_translation_dict = {'TTZ':'_TTZ_', 'TTJets':'_Top_pow_', 'TTXNoZ':'_TTXNoZ_', 'multiBoson':'_multiBoson_', 'DY':"_DY_HT_LO_" }

    # scale all MC histograms by the SF
    for year, directory in zip( [2016, 2017, 2018], args.directories ):
        for key in variation_data[directory].keys():
            # naming conventions (central variation):
            # [[<ROOT.TH1F object ("nbtags_data_data_c8eeb1ca_36fb_4e72_a3e3_8218fabfaed0") at 0xcf546b0>]],
            # [[<ROOT.TH1F object ("nbtags_signal_T2tt_800_100_187e383a_4a46_49ab_bfcc_cbe8f11cabe8") at 0xcf55a70>],
            #  [<ROOT.TH1F object ("nbtags_signal_T2tt_350_150_3d22a083_1d90_4db2_b5b1_ed88f876677e") at 0xcf56380>]],
            # [[<ROOT.TH1F object ("nbtags_mc_Top_pow_5065a108_9539_4a3d_bcc7_25c8aefec150") at 0xcf56f80>,
            #   <ROOT.TH1F object ("nbtags_mc_multiBoson_8c5e291e_37fe_47af_b019_a5e58e2b30c0") at 0xcf57b70>,
            #   <ROOT.TH1F object ("nbtags_mc_DY_HT_LO_04a40d9f_b4ac_42b8_a578_9a03e26041b2") at 0xcf58780>,
            #   <ROOT.TH1F object ("nbtags_mc_TTXNoZ_ed2cc2e1_1e40_4728_ac35_3f24b99c00c2") at 0xcf593d0>,
            #   <ROOT.TH1F object ("nbtags_mc_TTZ_55917593_e110_437a_b183_45947c2f3006") at 0xcf59fa0>]],
            for p_histo in variation_data[directory][key]['histos']:
               # p_histo is a list of plot_histos, for 'central' it is data/mc/signal, for variations it is just mc
                for s_h in p_histo:
                    for histo in s_h:
                        # If the histo name contains a string that we have in SFUnc[year], we scale the MC histogram
                        if '_mc_' not in histo.GetName(): continue
                        for var_key in SFUnc[year].keys():
                            if daniel_translation_dict[var_key] in histo.GetName():
                                histo.Scale(SFUnc[year][var_key]['val'])
                                #print "Scaled", key, var_key, histo.GetName(), SFUnc[year][var_key]['val']
                                break

            # Now we add the keys for the scale factor uncertainties! I'm so glad I had this idea. I hope it's a good one.
            # let's copy the central MC histos
            if 'central' in key:
                for var_key in SFUnc[year].keys():
                    # Define the new ScaleFactor systematic. It's correlated.
                    for mod, sigma in [("Down", -1), ("Up", +1)]:
                        # need to add the variations
                        variations["ScaleFactor"+var_key+mod] = {}
                        # make a new key with the same syntax we have in variation_data, e.g. (mode, "ScaleFactor"+XX+Up/Down)
                        new_key = (key[0], "ScaleFactor"+var_key+mod)

                        variation_data[directory][new_key] = {} 
                        variation_data[directory][new_key]['histos'] = []
                        for i_plot, plot in enumerate(plots):
                            # Now we find the central MC hists. They are already scaled.
                            mc_pos   = 3*i_plot+2 if args.signal else 2*i_plot+1
                            # Copy the central MC histos 
                            mc_histos = [ [ h.Clone() for h in s_h] for s_h in variation_data[directory][key]['histos'][mc_pos] ]
                            # Now we scale them according to the variation
                            for s_h in mc_histos:
                                for histo in s_h:
                                    if '_mc_' not in histo.GetName(): raise RuntimeError("You should have a MC histogram here, but the name doesn't match: %s" %  histo.GetName())
                                    if daniel_translation_dict[var_key] in histo.GetName():
                                        histo.Scale((1+sigma*SFUnc[year][var_key]['sigma']))
                                        #print "Scaled", key, var_key, histo.GetName(), SFUnc[year][var_key]['val']
                                        break
        
                            variation_data[directory][new_key]['histos'].append( mc_histos )
else:
    SFUnc = None

def niceTexName(texName):
    if "Run" in texName and "SF" in texName:
        return "Observed (SF)"
    elif "multi" in texName:
        return "Multiboson"
    elif "DY" in texName:
        return "Drell-Yan"
    else:
        return texName

# We plot now. 
if args.normalize: plot_subdirectory += "_normalized"
if args.beta:      plot_subdirectory += "_%s"%args.beta
if args.fitChi2:   plot_subdirectory += "_chi2"
if args.private:   plot_subdirectory += "_private"
for mode in ['mumu', 'ee', 'mue', 'SF', 'all']:
    for i_plot, plot in enumerate(plots):
        
        # for central (=no variation), we store plot_data_1, plot_mc_1, plot_data_2, plot_mc_2, ...
        if args.signal:
            data_histo_list     = accumulate_level_2( [variation_data[directory][(mode, 'central')]['histos'][3*i_plot] for directory in args.directories] )
            signal_histo_list   = accumulate_level_2( [variation_data[directory][(mode, 'central')]['histos'][3*i_plot+1]  for directory in args.directories] )
            mc_histo_list       = {'central': accumulate_level_2( [variation_data[directory][(mode, 'central')]['histos'][3*i_plot+2] for directory in args.directories]) }
        else:
            data_histo_list     = accumulate_level_2( [variation_data[directory][(mode, 'central')]['histos'][2*i_plot]  for directory in args.directories] )
            signal_histo_list   = []
            mc_histo_list       = {'central': accumulate_level_2( [variation_data[directory][(mode, 'central')]['histos'][2*i_plot+1] for directory in args.directories]) }

        data_histo_list[0][0].style = styles.errorStyle(ROOT.kBlack)
        for i_mc, mc in enumerate(stack_mc[0]):
            mc_histo_list['central'][0][i_mc].legendText = niceTexName(mc.texName)
            mc_histo_list['central'][0][i_mc].style = styles.fillStyle(mc.color) 
        if args.signal:
            for i_sig, sig in enumerate(stack_signal):
                signal_histo_list[i_sig][0].legendText = sig[0].texName

        # for the other variations, there is no data
        for variation in variations.keys():
            if variation=='central': continue
            mc_histo_list[variation] = accumulate_level_2( [variation_data[directory][(mode, variation)]['histos'][i_plot] for directory in args.directories ] )
        for directory in args.directories:
            mc_histo_list[directory] = {variation: variation_data[directory][(mode, variation)]['histos'][i_plot] for variation in variations.keys()}

        # Add histos, del the stack (which refers to MC only )
        if args.signal:
            plot.histos =  mc_histo_list['central'] + data_histo_list + signal_histo_list
        else:
            plot.histos =  mc_histo_list['central'] + data_histo_list
        #plot.stack  = Stack( mc, [data_sample] )
        #if args.signal != None: 
        #    plot.histos += signal_histo_list
        #    plot.stack.extend( [ [s] for s in signals ] ) 
        
        # Make boxes and ratio boxes
        boxes           = []
        ratio_boxes     = []
        # Compute all variied MC sums
        total_mc_histo  = {variation:add_histos( mc_histo_list[variation][0]) for variation in variations.keys() }
        # store also per directory
        for directory in args.directories:
            total_mc_histo[directory]  = {variation:add_histos( mc_histo_list[directory][variation][0]) for variation in variations.keys() }

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
            #print "bin {} | total_central_mc_yield {}".format(i_b, total_central_mc_yield)
            variance = 0.
            for systematic in systematics:
                # Use 'central-variation' (factor 1) and 0.5*(varUp-varDown)
                if 'central' in systematic['pair']: 
                    factor = 1
                else:
                    factor = 0.5
                # sum in quadrature
                if systematic['correlated']:
                    if i_b==total_mc_histo['central'].GetNbinsX() and overflowBin: # add overflow bin
                        up = factor*(total_mc_histo[systematic['pair'][0]].GetBinContent(i_b) + total_mc_histo[systematic['pair'][0]].GetBinContent(i_b+1))
                        down = factor*(total_mc_histo[systematic['pair'][1]].GetBinContent(i_b) + total_mc_histo[systematic['pair'][1]].GetBinContent(i_b+1))
                        # overflow bin can contain NaN
                        if not isnan(up) and not isnan(down):
                            variance += ( up - down )**2
                        else: 
                            print total_mc_histo[systematic['pair'][0]].GetBinContent(i_b), total_mc_histo[systematic['pair'][1]].GetBinContent(i_b), variance
                            variance += ( factor*(total_mc_histo[systematic['pair'][0]].GetBinContent(i_b) - total_mc_histo[systematic['pair'][1]].GetBinContent(i_b)) )**2
                    else:
                        up = total_mc_histo[systematic['pair'][0]].GetBinContent(i_b)
                        down = total_mc_histo[systematic['pair'][1]].GetBinContent(i_b)
                        # bin can contain NaN !? #FIXME leptonSFDown does contain NaN every once in a while
                        if not isnan(up) and not isnan(down):
                            variance += ( factor*( up - down ) )**2
                        else:
                            print systematic['name'], up, down, variance
                else:
                    if i_b==total_mc_histo['central'].GetNbinsX() and overflowBin: # add overflow bin
                        var = 0.
                        for directory in args.directories:
                            up = total_mc_histo[directory][systematic['pair'][0]].GetBinContent(i_b) + total_mc_histo[directory][systematic['pair'][0]].GetBinContent(i_b+1)
                            down = total_mc_histo[directory][systematic['pair'][1]].GetBinContent(i_b) + total_mc_histo[directory][systematic['pair'][1]].GetBinContent(i_b+1)
                            variance += ( factor*( up - down ) )**2
                    else:
                        for directory in args.directories:
                            variance += ( factor*(total_mc_histo[directory][systematic['pair'][0]].GetBinContent(i_b) - total_mc_histo[directory][systematic['pair'][1]].GetBinContent(i_b)) )**2

            if isnan(sqrt(variance)): assert False
            #print "{}|plot {} bin {}/{} variance: {:.3f}".format(mode, i_plot, i_b, total_mc_histo['central'].GetNbinsX(), sqrt(variance))
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

        # -------------------------------------------------------
        # DY pTmiss plot with chi-squared
        if args.fitChi2 and plot.name == "MET_significance_mc":
            chi2=ROOT.TF1("chi2","2.25*10**7*0.5*exp(-0.5*x)",0,100)
            chi2.SetNpx(1000)
            chi2.SetLineColor(ROOT.kRed)

            chi2_histo = chi2.GetHistogram()

            # Legend
            chi2_histo.legendText = "#chi^{2}(N_{dof}=2)"
            chi2_histo.legendOption = "l"
            chi2_histo.SetMarkerStyle(0)
            chi2_histo.SetMarkerSize(0)
            #chi2_histo.SetBorderSize(0)
            #ROOT.gStyle.SetLegendBorderSize(0)

            plot.histos[1].append(chi2_histo)
            plot.texX = 'p_{T}^{miss} significance'
            #chiSquared = []
        # -------------------------------------------------------

        for log in [True]:
            plot_directory_ = os.path.join(plot_directory, 'systematicPlots', 'combined', plot_subdirectory, args.selection, mode + ("_log" if log else ""))
            #if not max(l[0].GetMaximum() for l in plot.histos): continue # Empty plot
            texMode = "#mu#mu" if mode == "mumu" else "#mue" if mode == "mue" else mode
            if    mode == "all": plot.histos[1][0].legendText = "Observed"
            else:                plot.histos[1][0].legendText = "Observed (%s)"%(texMode)

            _drawObjects = []

            plotting.draw(plot,
              plot_directory = plot_directory_,
              ratio = {'yRange':(0.1,1.9), 'drawObjects':ratio_boxes},
              logX = False, logY = log, sorting = False,
              yRange = (0.03, "auto") if log else (0.6, "auto") if args.fitChi2 else (0.001, "auto"),
              #yRange = (0.6, 10**10) if log else (0.6, "auto") if args.fitChi2 else (0.001, "auto"),
              scaling = {0:1} if args.normalize else {},
              legend = ( (0.18,0.88-0.03*sum(map(len, plot.histos)),0.9,0.88), 2),
              drawObjects = drawObjects( args.scaling ) + boxes,
              copyIndexPHP = True, extensions = ["png", "pdf", "root"],
            )
