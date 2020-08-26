#!/usr/bin/env python
'''
Combination on card file level
'''

import ROOT
import os
import argparse
from RootTools.core.Sample import Sample
argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--logLevel',       action='store',        default='INFO',         nargs='?', choices=['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'TRACE', 'NOTSET'],             help="Log level for logging")
argParser.add_argument("--signal",         action='store',        default='T2tt',         nargs='?', choices=["T2tt","T8bbllnunu_XCha0p5_XSlep0p5","T8bbllnunu_XCha0p5_XSlep0p05","T8bbllnunu_XCha0p5_XSlep0p95","T2bW", "ttHinv", "TTbarDM"], help="which signal scan?")
argParser.add_argument("--version",        action='store',        default='v8',           nargs='?',                                                                          help="which signal scan?")
argParser.add_argument("--overwrite",      action = "store_true", default = False,                                                                                                             help="Overwrite existing output files")
argParser.add_argument("--controlRegions", action='store',        default='signalOnly',   nargs='?', choices=["control2016","controlAll","signalOnly","controlDYVV","controlTTZ","controlTT","fitAll"],                help="which signal scan?")
argParser.add_argument("--model",          action='store',        default='dim6top_LO',   nargs='?', choices=["dim6top_LO", "ewkDM"],                                                          help="which signal model?")
argParser.add_argument("--only",           action='store',        default=None,           nargs='?',                                                                                           help="pick only one signal point?")
argParser.add_argument("--skipYear",       action='store',        default=None, type=int, nargs='?', choices=[2016,2017,2018],                                                                 help="pick only one signal point?")
argParser.add_argument("--includeCR",      action='store_true',                                                                                                                                help="Do simultaneous SR and CR fit")
argParser.add_argument("--expected",       action='store_true',                                                                                                                                help="Do simultaneous SR and CR fit")
argParser.add_argument("--calcNuisances",  action='store_true',                                                                         help="Extract the nuisances and store them in text files?")
argParser.add_argument("--signalInjection",action='store_true',                                                                         help="Would you like some signal with your background?")
argParser.add_argument("--clean",          action='store_true',                                                                         help="Remove potentially failed fits?")
argParser.add_argument("--useTxt",         action='store_true',                                                                         help="Use txt files?")
argParser.add_argument("--significance",   action='store_true',                                                                         help="Calculate significances instead of limits?")
argParser.add_argument("--dryRun",   action='store_true',                                                                         help="Only write card file")


args = argParser.parse_args()

# Logging
import StopsDilepton.tools.logger as logger
logger = logger.get_logger(args.logLevel, logFile = None )
import RootTools.core.logger as logger_rt
logger_rt = logger_rt.get_logger(args.logLevel, logFile = None )

from math                               import sqrt
from copy                               import deepcopy


from StopsDilepton.analysis.Setup               import Setup

#from StopsDilepton.tools.resultsDB             import resultsDB
from StopsDilepton.tools.u_float                import u_float
from StopsDilepton.tools.user                   import analysis_results,  plot_directory
from StopsDilepton.tools.cardFileWriter         import cardFileWriter
from StopsDilepton.analysis.Cache               import Cache

from StopsDilepton.samples.nanoTuples_Summer16_postProcessed            import *
from StopsDilepton.samples.nanoTuples_Run2016_17Jul2018_postProcessed   import *
from StopsDilepton.samples.nanoTuples_Fall17_postProcessed              import *
from StopsDilepton.samples.nanoTuples_Run2017_31Mar2018_postProcessed   import *
from StopsDilepton.samples.nanoTuples_Autumn18_postProcessed            import *
from StopsDilepton.samples.nanoTuples_Run2018_PromptReco_postProcessed  import *

# some fake setup
setup = Setup(2016)

years = [2016,2017,2018]
if args.skipYear:
    years.remove(args.skipYear)
#analysis_results = '/afs/hephy.at/work/p/phussain/StopsDileptonLegacy/results/v3/'

analysis_results = analysis_results.replace('v8', args.version)

vetoList = ["T2tt_275_25"]

overWrite = args.overwrite
def wrapper(s):

    logger.info("Now working on %s", s.name)

    c = cardFileWriter.cardFileWriter()
    c.releaseLocation = os.path.abspath('.')
    cards = {}
    
    # get the seperated cards
    for year in years:
        sSubDir = 'expected' if args.expected else 'observed'
        if args.signalInjection: sSubDir += '_signalInjected'
        baseDir  = analysis_results+"/%s/%s/"%(year,args.controlRegions)
        limitDir = baseDir+"/cardFiles/%s/%s/"%(args.signal, sSubDir)
        cardFileName = os.path.join(limitDir, s.name+'_shapeCard.txt') if not args.useTxt else os.path.join(limitDir, s.name+'.txt')

        #print cardFileName
        if not os.path.isfile(cardFileName):
            #raise IOError("File %s doesn't exist!"%cardFileName)
            print "File %s doesn't exist!"%cardFileName
            return

        cards[year] = cardFileName
    
    baseDir  = baseDir.replace('2018','comb')
    limitDir = limitDir.replace('2018','comb')
    
    if not args.useTxt:
        cacheFileName = os.path.join(limitDir, 'calculatedLimits')
    else:
        cacheFileName = os.path.join(limitDir, 'calculatedLimits_txt')
    limitCache    = Cache(cacheFileName, verbosity=2)


    # run combine and store results in sqlite database
    if not os.path.isdir(limitDir):
        os.makedirs(limitDir)
    #resDB = resultsDB(limitDir+'/results.sq', "results", setup.resultsColumns)

    if   args.signal == "TTbarDM":                      sConfig = s.mChi, s.mPhi, s.type
    elif args.signal == "T2tt":                         sConfig = s.mStop, s.mNeu
    elif args.signal == "T2bt":                         sConfig = s.mStop, s.mNeu
    elif args.signal == "T2bW":                         sConfig = s.mStop, s.mNeu
    elif args.signal == "T8bbllnunu_XCha0p5_XSlep0p05": sConfig = s.mStop, s.mNeu
    elif args.signal == "T8bbllnunu_XCha0p5_XSlep0p09": sConfig = s.mStop, s.mNeu
    elif args.signal == "T8bbllnunu_XCha0p5_XSlep0p5":  sConfig = s.mStop, s.mNeu
    elif args.signal == "T8bbllnunu_XCha0p5_XSlep0p95": sConfig = s.mStop, s.mNeu
    elif args.signal == "ttHinv":                       sConfig = ("ttHinv", "2l")

    #overWrite = False

    if not overWrite and limitCache.contains(sConfig):
        
        res = limitCache.get(sConfig)
        logger.info("Found result for %s, reusing", s.name)


    else:
        print "*" *100
        combinedCard = c.combineCards( cards )
        print combinedCard
        if args.dryRun:
            return None
        try:
            res = c.calcLimit(combinedCard)
        except:
            print "Limit failed"
            res = {}
        print res 
        print "+"*10
        if args.significance:
            sig = c.calcSignif(combinedCard)
        try:
            significance = sig['-1.000']
        except:
            significance = -1
        res['significance'] = significance
        if args.calcNuisances:
            c.calcNuisances(combinedCard)

            ###################
            # extract the SFs #
            ###################
            #if not args.useTxt 
            if True:
                # Would be a bit more complicated with the classical txt files, so only automatically extract the SF when using shape based datacards
                from StopsDilepton.tools.getPostFit import getPrePostFitFromMLF
                print combinedCard
                print cardFileName
                combineWorkspace = combinedCard.replace('.txt','_FD.root')
                print "Extracting fit results from %s"%combineWorkspace

                postFitResults = getPrePostFitFromMLF(combineWorkspace)

                top_prefit  = 0
                top_postfit = 0
                ttZ_prefit  = 0
                ttZ_postfit = 0
                DY_prefit   = 0
                DY_postfit  = 0
                MB_prefit   = 0
                MB_postfit   = 0
                other_prefit  = 0
                other_postfit  = 0

                for year in years:
                    binName = "dc_%s"%year

                    top_prefit  += postFitResults['results']['shapes_prefit'][binName]['TTJets']# + postFitResults['results']['shapes_prefit'][binName]['TTJetsG'] + postFitResults['results']['shapes_prefit'][binName]['TTJetsNG']
                    top_postfit += postFitResults['results']['shapes_fit_b'][binName]['TTJets'] #+ postFitResults['results']['shapes_fit_b'][binName]['TTJetsG'] + postFitResults['results']['shapes_fit_b'][binName]['TTJetsNG']

                    ttZ_prefit  += postFitResults['results']['shapes_prefit'][binName]['TTZ']
                    ttZ_postfit += postFitResults['results']['shapes_fit_b'][binName]['TTZ']

                    DY_prefit  += postFitResults['results']['shapes_prefit'][binName]['DY']
                    DY_postfit += postFitResults['results']['shapes_fit_b'][binName]['DY']

                    MB_prefit  += postFitResults['results']['shapes_prefit'][binName]['multiBoson']
                    MB_postfit += postFitResults['results']['shapes_fit_b'][binName]['multiBoson']

                    other_prefit  += postFitResults['results']['shapes_prefit'][binName]['TTXNoZ']
                    other_postfit += postFitResults['results']['shapes_fit_b'][binName]['TTXNoZ']

                print
                print "## Scale Factors for backgrounds: ##"
                print "{:20}{:4.2f}{:3}{:4.2f}".format('top:',          (top_postfit/top_prefit).val, '+/-',  top_postfit.sigma/top_postfit.val)
                print "{:20}{:4.2f}{:3}{:4.2f}".format('ttZ:',          (ttZ_postfit/ttZ_prefit).val, '+/-',  ttZ_postfit.sigma/ttZ_postfit.val)
                print "{:20}{:4.2f}{:3}{:4.2f}".format('Drell-Yan:',    (DY_postfit/DY_prefit).val,   '+/-',  DY_postfit.sigma/DY_postfit.val)
                print "{:20}{:4.2f}{:3}{:4.2f}".format('multiBoson:',   (MB_postfit/MB_prefit).val,   '+/-',  MB_postfit.sigma/MB_postfit.val)
                print "{:20}{:4.2f}{:3}{:4.2f}".format('other:',        (other_postfit/other_prefit).val,   '+/-',  other_postfit.sigma/other_postfit.val)




        #print "Result: %r obs %5.3f exp %5.3f -1sigma %5.3f +1sigma %5.3f"%(sString, res['-1.000'], res['0.500'], res['0.160'], res['0.840'])
        ## this is dangerous, but no way around it atm ##
        if (args.signal.startswith("T8") and s.mStop<301):# or (args.signal.startswith("T2bW") and s.mStop<201):
            for k in res:
                res[k] *= 0.01
        limitCache.add(sConfig, res)
        logger.info("Adding results to database")
        logger.info("Results stored in %s", limitDir )


    if res:
        if   args.signal == "TTbarDM":                        sString = "mChi %i mPhi %i type %s" % sConfig
        elif args.signal == "T2tt":                           sString = "mStop %i mNeu %i" % sConfig
        elif args.signal == "T2bt":                           sString = "mStop %i mNeu %i" % sConfig
        elif args.signal == "T2bW":                           sString = "mStop %i mNeu %i" % sConfig
        elif args.signal == "T8bbllnunu_XCha0p5_XSlep0p05":   sString = "mStop %i mNeu %i" % sConfig
        elif args.signal == "T8bbllnunu_XCha0p5_XSlep0p09":   sString = "mStop %i mNeu %i" % sConfig
        elif args.signal == "T8bbllnunu_XCha0p5_XSlep0p5":    sString = "mStop %i mNeu %i" % sConfig
        elif args.signal == "T8bbllnunu_XCha0p5_XSlep0p95":   sString = "mStop %i mNeu %i" % sConfig
        elif args.signal == "ttHinv":                         sString = "ttH->inv"
    #print sString, res
    try:
        print "Result: {:30} obs {:<10.2f} exp {:<10.2f} -1sigma {:<10.2f} +1sigma {:<10.2f}".format(sString, res['-1.000'], res['0.500'], res['0.160'], res['0.840'])
        if res['-1.000']>4*res['0.840']: # was 4
            print "WARNING: This point could be problematic!"
            if args.clean:
                vetoList += [ s.name ]
        return sConfig, res
    except:
        print "Problem with limit: %r"%str(res)
        return None
    try:
        print "Significance: %.2f"%res['significance']
    except:
        print "No significance calculated"


if args.signal == "T2tt":
    data_directory              = '/afs/hephy.at/data/cms09/nanoTuples/'
    postProcessing_directory    = 'stops_2016_nano_v0p22/dilep/'
    from StopsDilepton.samples.nanoTuples_FastSim_Fall17_postProcessed import signals_T2tt as jobs
elif args.signal == "T2bW":
    data_directory              = '/afs/hephy.at/data/cms09/nanoTuples/'
    postProcessing_directory    = 'stops_2016_nano_v0p22/dilep/'
    from StopsDilepton.samples.nanoTuples_FastSim_Summer16_postProcessed import signals_T2bW as jobs
elif args.signal == "T8bbllnunu_XCha0p5_XSlep0p05":
    data_directory              = '/afs/hephy.at/data/cms09/nanoTuples/'
    postProcessing_directory    = 'stops_2016_nano_v0p22/dilep/'
    from StopsDilepton.samples.nanoTuples_FastSim_Summer16_postProcessed import signals_T8bbllnunu_XCha0p5_XSlep0p05 as jobs
elif args.signal == "T8bbllnunu_XCha0p5_XSlep0p5":
    data_directory              = '/afs/hephy.at/data/cms09/nanoTuples/'
    postProcessing_directory    = 'stops_2016_nano_v0p22/dilep/'
    from StopsDilepton.samples.nanoTuples_FastSim_Summer16_postProcessed import signals_T8bbllnunu_XCha0p5_XSlep0p5 as jobs
elif args.signal == "T8bbllnunu_XCha0p5_XSlep0p95":
    data_directory              = '/afs/hephy.at/data/cms09/nanoTuples/'
    postProcessing_directory    = 'stops_2016_nano_v0p22/dilep/'
    from StopsDilepton.samples.nanoTuples_FastSim_Summer16_postProcessed import signals_T8bbllnunu_XCha0p5_XSlep0p95 as jobs
elif args.signal == "ttHinv":
    data_directory              = '/afs/hephy.at/data/cms09/nanoTuples/'
    postProcessing_directory    = 'stops_2018_nano_v0p22/dilep/'
    ttH_HToInvisible_M125 = Sample.fromDirectory(name="ttH_HToInvisible_M125", treeName="Events", isData=False, color=1, texName="ttH(125)", directory=os.path.join(data_directory,postProcessing_directory,'ttH_HToInvisible'))
    jobs = [ttH_HToInvisible_M125]
elif args.signal == "TTbarDM":
    data_directory              = '/afs/hephy.at/data/cms03/nanoTuples/'
    postProcessing_directory    = 'stops_2016_nano_v0p25/dilep/'
    from StopsDilepton.samples.nanoTuples_Summer16_TTDM_postProcessed import signals as jobs


# FIXME: removing 1052_0 from list
for i, j in enumerate(jobs):
    if j.name == "T8bbllnunu_XCha0p5_XSlep0p05_1052_0":
        print "~removing ", j.name
        del jobs[i]

#vetoList = ['T2tt_150_63', 'T2tt_200_100', 'T2tt_200_113', 'T2tt_250_150', 'T2tt_200_0', 'T2tt_526_438', 'T2tt_550_450', 'T2tt_576_475', 'T2tt_600_475', 'T2tt_600_514', 'T2tt_626_526', 'T2tt_626_538', 'T2tt_650_475']
#vetoList += ['T2tt_150_63', 'T2tt_200_100', 'T2tt_200_113', 'T2tt_250_150', 'T2tt_200_0']

allJobs = [j for j in jobs if (j.name not in vetoList)]
if args.only is not None:
    if args.only.isdigit():
        wrapper(jobs[int(args.only)])
    else:
        jobNames = [ x.name for x in jobs ]
        wrapper(jobs[jobNames.index(args.only)])
    exit(0)
#i= [j for j in jobs if j.name != 'T2tt_150_63']

allJobs.sort(key=lambda x: x.name, reverse=False)
results = map(wrapper, allJobs)
results = [r for r in results if r]


## do something nice with pandas
import pandas as pd
import pickle
 
# start by getting a dataframe with the results
r_list = []
for r in results:
    tmp = { 'stop':r[0][0], 'lsp':r[0][1] }
    tmp.update(r[1])
    r_list.append(tmp)
results_df = pd.DataFrame(r_list)


#results_df[results_df['0.500']<1.].sort_values('stop')[['stop','lsp']]

if args.signal=='ttHinv' or args.signal == 'TTbarDM': exit()

exculuded_exp_stop = results_df[results_df['0.500']<1.]['stop'].tolist()
exculuded_exp_lsp  = results_df[results_df['0.500']<1.]['lsp'].tolist()

nexculuded_exp_stop = results_df[results_df['0.500']>1.]['stop'].tolist()
nexculuded_exp_lsp  = results_df[results_df['0.500']>1.]['lsp'].tolist()

exculuded_obs_stop = results_df[results_df['-1.000']<1.]['stop'].tolist()
exculuded_obs_lsp  = results_df[results_df['-1.000']<1.]['lsp'].tolist()

nexculuded_obs_stop = results_df[results_df['-1.000']>1.]['stop'].tolist()
nexculuded_obs_lsp  = results_df[results_df['-1.000']>1.]['lsp'].tolist()

#bulk_df = results_df[results_df['stop']-results_df['lsp']>170]
bulk_df = results_df[results_df['stop']%5==0][results_df['lsp']%5==0] #results_df[abs(results_df['stop']-results_df['lsp']-175)>15]
comp_df = results_df[results_df['stop']-results_df['lsp']<175]

#########################################################################################
# Process the results. Make 2D hists for SUSY scans, or table for the DM interpretation #
#########################################################################################

# Make histograms for T2tt
baseDir  = analysis_results+"/comb/%s/"%(args.controlRegions)
limitPrefix = args.signal
if not os.path.isdir(os.path.join(baseDir, 'limits', args.signal, limitPrefix)):
    os.makedirs(os.path.join(baseDir, 'limits', args.signal, limitPrefix))
limitResultsFilename = os.path.join(baseDir, 'limits', args.signal, limitPrefix,'limitResults.root')

print "Root file is here:", limitResultsFilename

## new try, other thing is buggy
def toGraph2D(name,title,length,x,y,z):
    result = ROOT.TGraph2D(length)
    result.SetName(name)
    result.SetTitle(title)
    for i in range(length):
        result.SetPoint(i,x[i],y[i],z[i])
    h = result.GetHistogram()
    h.SetMinimum(min(z))
    h.SetMaximum(max(z))
    c = ROOT.TCanvas()
    result.Draw()
    del c
    #res = ROOT.TGraphDelaunay(result)
    return result

mStop_list = []
mLSP_list  = []
exp_list   = []
obs_list   = []
signif_list   = []
exp_up_list   = []
exp_down_list   = []

for r in results:
    s, res = r
    mStop, mNeu = s
    #if mStop%50>0: continue
    #if mNeu%50>0 and not mNeu>(mStop-125): continue
    mStop_list.append(mStop)
    mLSP_list.append(mNeu)
    exp_list.append(res['0.500'])
    exp_up_list.append(res['0.160'])
    exp_down_list.append(res['0.840'])
    obs_list.append(res['-1.000'])
    try:
        signif_list.append(res['significance'])
    except:
        signif_list.append(-1)

scatter         = ROOT.TGraph(len(mStop_list))
scatter.SetName('scatter')
for i in range(len(mStop_list)):
    scatter.SetPoint(i,mStop_list[i],mLSP_list[i])

scatter_excl_exp = ROOT.TGraph(len(exculuded_exp_stop))
scatter_excl_exp.SetName("scatter_excl_exp")
for i in range(len(exculuded_exp_stop)):
    scatter_excl_exp.SetPoint(i, exculuded_exp_stop[i], exculuded_exp_lsp[i])


scatter_nexcl_exp = ROOT.TGraph(len(nexculuded_exp_stop))
scatter_nexcl_exp.SetName("scatter_nexcl_exp")
for i in range(len(nexculuded_exp_stop)):
    scatter_nexcl_exp.SetPoint(i, nexculuded_exp_stop[i], nexculuded_exp_lsp[i])

scatter_excl_obs = ROOT.TGraph(len(exculuded_obs_stop))
scatter_excl_obs.SetName("scatter_excl_obs")
for i in range(len(exculuded_obs_stop)):
    scatter_excl_obs.SetPoint(i, exculuded_obs_stop[i], exculuded_obs_lsp[i])


scatter_nexcl_obs = ROOT.TGraph(len(nexculuded_obs_stop))
scatter_nexcl_obs.SetName("scatter_nexcl_obs")
for i in range(len(nexculuded_obs_stop)):
    scatter_nexcl_obs.SetPoint(i, nexculuded_obs_stop[i], nexculuded_obs_lsp[i])


exp_graph       = toGraph2D('exp','exp',len(mStop_list),mStop_list,mLSP_list,exp_list)
exp_up_graph    = toGraph2D('exp_up','exp_up',len(mStop_list),mStop_list,mLSP_list,exp_up_list)
exp_down_graph  = toGraph2D('exp_down','exp_down',len(mStop_list),mStop_list,mLSP_list,exp_down_list)
#obs_graph       = toGraph2D('obs','obs',len(mStop_list),mStop_list,mLSP_list,obs_list)
signif_graph    = toGraph2D('signif','signif',len(mStop_list),mStop_list,mLSP_list,signif_list)

obs_bulk = toGraph2D('obs_bulk','obs_bulk',len(bulk_df['stop'].tolist()),bulk_df['stop'].tolist(),bulk_df['lsp'].tolist(),bulk_df['-1.000'].tolist())
obs_graph = toGraph2D('obs','obs',len(bulk_df['stop'].tolist()),bulk_df['stop'].tolist(),bulk_df['lsp'].tolist(),bulk_df['-1.000'].tolist())
obs_comp = toGraph2D('obs_comp','obs_comp',len(comp_df['stop'].tolist()),comp_df['stop'].tolist(),comp_df['lsp'].tolist(),comp_df['-1.000'].tolist())

pickle.dump(results_df, file(limitResultsFilename.replace('root', 'pkl'), 'w'))

outfile = ROOT.TFile(limitResultsFilename, "recreate")
scatter        .Write()
scatter_excl_exp.Write()
scatter_nexcl_exp.Write()
scatter_excl_obs.Write()
scatter_nexcl_obs.Write()
exp_graph      .Write()
exp_down_graph .Write()
exp_up_graph   .Write()
obs_graph      .Write()
obs_bulk      .Write()
obs_comp      .Write()
signif_graph      .Write()
outfile.Close()
print "Written %s"%limitResultsFilename

