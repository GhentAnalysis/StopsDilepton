#!/usr/bin/env python
import ROOT
import os
import argparse
import shutil

from StopsDilepton.tools.user           import combineReleaseLocation, analysis_results, plot_directory


releaseLocation = combineReleaseLocation+"/HiggsAnalysis/CombinedLimit/"

argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--logLevel',       action='store', default='INFO',          nargs='?', choices=['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'TRACE', 'NOTSET'],    help="Log level for logging")
argParser.add_argument("--signal",         action='store', default='TTbarDM',       nargs='?', choices=["T2tt","TTbarDM","ttHinv"],  help="Which signal?")
argParser.add_argument("--removeDir",      action='store_true',                                                             help="Remove the directory in the combine release after study is done?")
argParser.add_argument("--cores",          action='store', default=8,               nargs='?',                              help="Run on n cores in parallel")
argParser.add_argument("--only",           action='store', default=None,            nargs='?',                              help="pick only one masspoint?")
argParser.add_argument("--bkgOnly",        action='store_true',                                                             help="Allow no signal?")
args = argParser.parse_args()


# Logging
import StopsDilepton.tools.logger as logger
logger = logger.get_logger(args.logLevel, logFile = None )
import RootTools.core.logger as logger_rt
logger_rt = logger_rt.get_logger(args.logLevel, logFile = None )


def wrapper(s):
    logger.info("Processing mass point %s"%s.name)
    cardFile = "%s.txt"%s.name
    cardFilePath = "%s/fitAll/cardFiles/%s/%s"%(analysis_results, args.signal, cardFile)
    combineDirname = os.path.join(releaseLocation, s.name)
    logger.info("Creating %s"%combineDirname)
    if not os.path.isdir(combineDirname): os.makedirs(combineDirname)
    shutil.copyfile(cardFilePath,combineDirname+'/'+cardFile)
    if args.bkgOnly:
        prepWorkspace   = "text2workspace.py %s --X-allow-no-signal -m 125"%cardFile
    else:
        prepWorkspace   = "text2workspace.py %s -m 125"%cardFile
    if args.bkgOnly:
        robustFit       = "combineTool.py -M Impacts -d %s.root -m 125 --doInitialFit --robustFit 1 --rMin -0.01 --rMax 0.0"%s.name
        impactFits      = "combineTool.py -M Impacts -d %s.root -m 125 --robustFit 1 --doFits --parallel %s --rMin -0.01 --rMax 0.0"%(s.name,str(args.cores))
    else:
        robustFit       = "combineTool.py -M Impacts -d %s.root -m 125 --doInitialFit --robustFit 1 --rMin -10 --rMax 10"%s.name
        impactFits      = "combineTool.py -M Impacts -d %s.root -m 125 --robustFit 1 --doFits --parallel %s --rMin -10 --rMax 10"%(s.name,str(args.cores))
    extractImpact   = "combineTool.py -M Impacts -d %s.root -m 125 -o impacts.json"%s.name
    plotImpacts     = "plotImpacts.py -i impacts.json -o impacts"
    combineCommand  = "cd %s;eval `scramv1 runtime -sh`;%s;%s;%s;%s;%s"%(combineDirname,prepWorkspace,robustFit,impactFits,extractImpact,plotImpacts)
    logger.info("Will run the following command, might take a few hours:\n%s"%combineCommand)
    
    os.system(combineCommand)
    
    plotDir = plot_directory + "/impactsV4/"
    if not os.path.isdir(plotDir): os.makedirs(plotDir)
    if args.bkgOnly:
        shutil.copyfile(combineDirname+'/impacts.pdf', "%s/%s_bkgOnly.pdf"%(plotDir,s.name))
    else:
        shutil.copyfile(combineDirname+'/impacts.pdf', "%s/%s.pdf"%(plotDir,s.name))
    logger.info("Copied result to %s"%plotDir)

    if args.removeDir:
        logger.info("Removing directory in release location")
        rmtree(combineDirname)


if   args.signal == "T2tt":
    postProcessing_directory = "postProcessed_80X_v40/dilepTiny"
    from StopsDilepton.samples.cmgTuples_FastSimT2tt_mAODv2_25ns_postProcessed import signals_T2tt as jobs
elif args.signal == "TTbarDM":
    postProcessing_directory = "postProcessed_80X_v35/dilepTiny"
    from StopsDilepton.samples.cmgTuples_FullSimTTbarDM_mAODv2_25ns_postProcessed import signals_TTbarDM as jobs
elif args.signal == "ttHinv":
    postProcessing_directory = "postProcessed_80X_v35/dilepTiny"
    from StopsDilepton.samples.cmgTuples_Higgs_mAODv2_25ns_postProcessed import *
    jobs = [ttH_HToInvisible_M125]
else: raise Exception("%s not implemented!"%args.signal)


if args.only is not None:
  wrapper(jobs[int(args.only)])
  exit(0)

results = map(wrapper, jobs)
results = [r for r in results if r]
