#!/usr/bin/env python
import ROOT
import os
import argparse
import shutil

from StopsDilepton.tools.user           import analysis_results, plot_directory


argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--logLevel',       action='store', default='INFO',          nargs='?', choices=['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'TRACE', 'NOTSET'],    help="Log level for logging")
argParser.add_argument("--signal",         action='store', default='T2tt',          nargs='?', choices=["T2tt","TTbarDM","ttHinv"],  help="Which signal?")
argParser.add_argument("--removeDir",      action='store_true',                                                             help="Remove the directory in the combine release after study is done?")
argParser.add_argument("--expected",       action='store_true',                                                             help="Use expected results?")
argParser.add_argument("--observed",       action='store_true',                                                             help="Use expected results?")
argParser.add_argument("--combined",       action='store_true',                                                             help="Use expected results?")
argParser.add_argument("--signalInjection",action='store_true',                                                             help="Inject some signal?")
argParser.add_argument("--useTxt",         action='store_true',                                                             help="Use txt based card files?")
argParser.add_argument("--cores",          action='store', default=8,               nargs='?',                              help="Run on n cores in parallel")
argParser.add_argument("--year",           action='store', default=2017,               nargs='?',                           help="Which year?")
argParser.add_argument("--only",           action='store', default=None,            nargs='?',                              help="pick only one masspoint?")
argParser.add_argument("--bkgOnly",        action='store_true',                                                             help="Allow no signal?")
args = argParser.parse_args()


# Logging
import StopsDilepton.tools.logger as logger
logger = logger.get_logger(args.logLevel, logFile = None )
import RootTools.core.logger as logger_rt
logger_rt = logger_rt.get_logger(args.logLevel, logFile = None )


year = int(args.year)

def wrapper(s):
    logger.info("Processing mass point %s"%s.name)
    cardFile = "%s_combination_shapeCard.txt"%s.name if not args.useTxt else "%s.txt"%s.name
    #analysis_results = '/afs/hephy.at/work/p/phussain/StopsDileptonLegacy/results/v2/'
    sSubDir = 'expected' if args.expected else 'observed'
    if args.signalInjection: sSubDir += '_signalInjected'
    cardFilePath = "%s/%s/fitAll/cardFiles/%s/%s/%s"%(analysis_results, args.year if not args.combined else 'COMBINED', args.signal, sSubDir, cardFile)
    #cardFilePath = "%s/%s/controlAll/cardFiles/%s/%s/%s"%(analysis_results, args.year if not args.combined else 'COMBINED', args.signal, 'expected' if args.expected else 'observed', cardFile)
    #cardFilePath = "%s/%s/signalOnly/cardFiles/%s/%s/%s"%(analysis_results, args.year if not args.combined else 'COMBINED', args.signal, 'expected' if args.expected else 'observed', cardFile)
    #cardFilePath = "%s/%s/controlDYVV/cardFiles/%s/%s/%s"%(analysis_results, args.year if not args.combined else 'COMBINED', args.signal, 'expected' if args.expected else 'observed', cardFile)
    combineDirname = os.path.join(os.path.abspath('.'), s.name)
    print cardFilePath
    logger.info("Creating %s"%combineDirname)
    if not os.path.isdir(combineDirname): os.makedirs(combineDirname)
    shutil.copyfile(cardFilePath,combineDirname+'/'+cardFile)
    if not args.combined: shutil.copyfile(cardFilePath.replace('shapeCard.txt','shape.root'),combineDirname+'/'+cardFile.replace('shapeCard.txt','shape.root'))
    shutil.copyfile(cardFilePath.replace('shapeCard.txt', 'shapeCard.root'),combineDirname+'/'+cardFile.replace('shapeCard.txt', 'shapeCard.root'))
    prepWorkspace   = "text2workspace.py %s -m 125"%cardFile
    if args.bkgOnly:
        robustFit       = "combineTool.py -M Impacts -d %s_shapeCard.root -m 125 -t -1 --expectSignal 0 --doInitialFit --robustFit 1 --rMin -10 --rMax 10"%s.name
        impactFits      = "combineTool.py -M Impacts -d %s_shapeCard.root -m 125 -t -1 --expectSignal 0 --robustFit 1 --doFits --parallel %s --rMin -10 --rMax 10"%(s.name,str(args.cores))
    elif args.observed:
        robustFit       = "combineTool.py -M Impacts -d %s_shapeCard.root -m 125 --doInitialFit --robustFit 1 --rMin -10 --rMax 10"%s.name
        impactFits      = "combineTool.py -M Impacts -d %s_shapeCard.root -m 125 --robustFit 1 --doFits --parallel %s --rMin -10 --rMax 10"%(s.name,str(args.cores))
    else:
        robustFit       = "combineTool.py -M Impacts -d %s_shapeCard.root -m 125 -t -1 --expectSignal 1 --doInitialFit --robustFit 1 --rMin -10 --rMax 10"%s.name
        impactFits      = "combineTool.py -M Impacts -d %s_shapeCard.root -m 125 -t -1 --expectSignal 1 --robustFit 1 --doFits --parallel %s --rMin -10 --rMax 10"%(s.name,str(args.cores))
    extractImpact   = "combineTool.py -M Impacts -d %s_shapeCard.root -m 125 -o impacts.json"%s.name
    plotImpacts     = "plotImpacts.py -i impacts.json -o impacts"
    combineCommand  = "cd %s;%s;%s;%s;%s;%s"%(combineDirname,prepWorkspace,robustFit,impactFits,extractImpact,plotImpacts)
    logger.info("Will run the following command, might take a few hours:\n%s"%combineCommand)
    
    os.system(combineCommand)
    
    plotDir = plot_directory + "/impacts_combination/"
    if args.expected:
        s.name += '_expected'
    if args.bkgOnly:
        s.name += '_bkgOnly'
    if args.observed:
        s.name += '_observed'
    if not os.path.isdir(plotDir): os.makedirs(plotDir)
    elif args.combined:
        shutil.copyfile(combineDirname+'/impacts.pdf', "%s/%s_combined%s.pdf"%(plotDir,s.name,'_signalInjected' if args.signalInjection else ''))
    elif args.year:
        shutil.copyfile(combineDirname+'/impacts.pdf', "%s/%s_%s%s.pdf"%(plotDir,s.name,args.year,'_signalInjected' if args.signalInjection else ''))
    else:
        shutil.copyfile(combineDirname+'/impacts.pdf', "%s/%s.pdf"%(plotDir,s.name))
    logger.info("Copied result to %s"%plotDir)

    if args.removeDir:
        logger.info("Removing directory in release location")
        rmtree(combineDirname)


if args.signal == "T2tt":
    data_directory              = '/afs/hephy.at/data/cms07/nanoTuples/'
    postProcessing_directory    = 'stops_2016_nano_v0p22/dilep/'
    from StopsDilepton.samples.nanoTuples_FastSim_Fall17_postProcessed import signals_T2tt as jobs

allJobs = [j.name for j in jobs]

if args.only is not None:
    if args.only.isdigit():
        wrapper(jobs[int(args.only)])
    else:
        jobNames = [ x.name for x in jobs ]
        wrapper(jobs[jobNames.index(args.only)])
    exit(0)


results = map(wrapper, jobs)
results = [r for r in results if r]
