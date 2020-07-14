'''
/afs/hephy.at/work/p/phussain/StopsDileptonLegacy/results/v3//COMBINED/controlAll//cardFiles/T2tt/observed/
Get a signal region plot from the cardfiles
'''

#!/usr/bin/env python
from optparse import OptionParser
parser = OptionParser()
parser.add_option("--signal",               dest='signal',  action='store', default='T2tt',    choices=["T2tt", "T2bW", "ttHinv", "TTbarDM"], help="which signal?")
parser.add_option("--only",                 dest='only',  action='store', default='T2tt_800_100', help="which masspoints??")
parser.add_option("--small",                action='store_true', help="small?")
parser.add_option('--logLevel',             dest="logLevel",              default='INFO',              action='store',      help="log level?", choices=['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'TRACE', 'NOTSET'])
parser.add_option('--blinded',              action="store_true")
parser.add_option('--overwrite',            dest="overwrite", default = False, action = "store_true", help="Overwrite existing output files, bool flag set to True  if used")
parser.add_option('--expected',             action = "store_true", help="Run expected?")
parser.add_option('--preliminary',          action = "store_true", help="Run expected?")
parser.add_option('--testGrayscale',        action = "store_true", help="Do the most important test for this collaboration?")
parser.add_option('--signalCard',           action = "store_true", help="Make card for signal?")
parser.add_option("--year",                 action='store',      default=0, type="int", help='Which year?')
parser.add_option("--postFix",              action='store',      default="", help='Add sth?')
parser.add_option("--version",              action='store', default='v8',  help='which version to use?')
(options, args) = parser.parse_args()

# Standard imports
import ROOT
import os
import sys
import pickle
import math

# Analysis
from StopsDilepton.tools.u_float           import u_float
from StopsDilepton.tools.user              import plot_directory, analysis_results
from StopsDilepton.samples.color           import color
from StopsDilepton.tools.getPostFit        import *
from StopsDilepton.analysis.regions        import *
from StopsDilepton.tools.cardFileWriter import cardFileWriter

from RootTools.core.standard import *
# logger
import StopsDilepton.tools.logger as logger
import RootTools.core.logger as logger_rt
logger    = logger.get_logger(   options.logLevel, logFile = None)
logger_rt = logger_rt.get_logger(options.logLevel, logFile = None)

#lumiStr = setup.dataLumi/1000
years=[2016,2017,2018]

## ttH example: /afs/hephy.at/data/cms05/StopsDileptonLegacy/results/v7/2018/fitAll/cardFiles/ttHinv/observed/ttH_HToInvisible_M125_shapeCard.txt
cardName = "%s"%(options.only) if not options.signal == 'ttHinv' else "ttH_HToInvisible_M125"

if options.signal == 'TTbarDM':
    spin, mStop, mLSP = options.only.split('_')[1:4]
else:
    spin = None
    mStop, mLSP = options.only.split('_')[1:3]

mStop, mLSP = int(mStop), int(mLSP)

analysis_results = analysis_results.replace('v8', options.version)

cardDir = analysis_results+"/COMBINED/fitAll/cardFiles/%s/%s/"%(options.signal,'expected' if options.expected else 'observed')

cardFile = "%s/%s.txt"%(cardDir, cardName)


logger.info("Starting to work from cardfile %s"%cardFile)

# This is ugly.
import subprocess
if not os.path.isfile(cardFile):
    for year in years:
        cmd = "python ../run/run_limit.py --signal %s --year %s --only=%s --dryRun --useTxt --version %s --skipFitDiagnostics --unblind --fitAll"%(options.signal, year, options.only, options.version)
        logger.info("Running cmd: %s", cmd)
        subprocess.call(cmd, shell=True)
    cmd = "python ../run/run_combination.py --signal %s --only=%s --dryRun --useTxt --version %s --controlRegions fitAll"%(options.signal, options.only, options.version)
    subprocess.call(cmd, shell=True)

preFitHist={}
postFitHist={}
bhistos=[]
hists={}
histos={}
bkgHist=[]
processes = [   ('TTJets',     'TTJets',     't#bar{t}/t'),
                ('DY',         'DY',         'Drell-Yan'),
                ('multiBoson', 'multiBoson', 'VV/VVV'),
                ('TTZ',        'TTZDL',      't#bar{t}Z'),
                ('TTXNoZ',     'TTXNoZ',     't#bar{t}X, rare')]

# tex, name, list of nuisances, combine nuisances?

shape = 'shape'

systematics = [ 
    ("Jet energy scale",                        "JES",      ["JEC_2016", "JEC_2017", "JEC_2018"],                    shape,  True), 
    ("Integrated luminosity",                   "lumi",     ["Lumi_2016", "Lumi_2017", "Lumi_2018"],                'lnN',  True),
    ("Pileup modeling",                         "PU_weight",       ["PU_2016", "PU_2017", "PU_2018"],                shape,  True),
    ("Jet energy resolution",                   "JER",      ["JER_2018", "JER_2017",  "JER_2016"],                   shape,  True),
    ("Modeling of unclust. en.",                "unclEn",   ["unclEn_2016", "unclEn_2017", "unclEn_2018"],           shape,  True), 
    ("Trigger efficiency",                      "trigger2l",  ["trigger_2018", "trigger_2016", "trigger_2017"],        shape,  True),
    ("b tagging light flavor",                  "SFl",      ["SFl_2018", "SFl_2017", "SFl_2016"],                    shape,  True),
    ("b tagging heavy flavor",                  "b",        ["SFb_2017", "SFb_2016", "SFb_2018"],                    shape,  True),
    ("Lepton scale factors",                    "leptonSF", ["leptonSF"],                                            shape,  False),
    ("L1 prefire correction",                   "Prefire_weight", ["L1prefire"],                                     shape,  False),
    ("0 missing hit scale factor",              "leptonHit0SF", ["leptonHit0SF"],                                    shape,  False),
    ("Impact parameter scale factor",           "leptonSIP3DSF", ["leptonSIP3DSF"],                                  shape,  False),
    ("PDF choice",                              "PDF_Weight",      ["PDF"],                                                 shape,  False),
    ("$p_{T}(\\textrm{top})$",                  "toppt",    ["topPt"],                                               shape,  False),
    ("t#bar{t} cross section",                  "topXSec",  ["topXSec"],                                             shape,  False),
    ("t#bar{t}Z background",                    "ttZ_SR",   ["ttZ_SR"],                                              shape,  False),
    ("fake/non-prompt leptons",                 "topFakes", ["topFakes"],                                            shape,  False),
    ("Drell-Yan background",                    "DY_SR",    ["DY_SR"],                                               shape,  False),
    ("non-gaussian jet mismeasurements",        "topNonGauss", ["topNonGauss"],                                      shape,  False),
    ("multiboson background",                   "MB_SR",    ["MB_SR"],                                               shape,  False),
    ("$\mu_{R}$ and $\mu_{F}$ choice $(t\bar{t})$", "scaleTT", ["scaleTT"],                                          shape,  False),
    ("rare background",                         "Rare2l",     ["rare"],                                                'lnN',  False),
    ("$\mu_{R}$ and $\mu_{F}$ choice $(t\bar{t}Z)$", "scaleTTZ", ["scaleTTZ"],                                       shape,  False),
    ("Drell-Yan tail",                          "DY_hMT2ll", ["DY_hMT2ll"],                                          shape,  False), 
]

#from infoFromCards import *
#getPUFC = getPreFitUncFromCard

from makeDataFrameFromCard import getDict

result = getDict(cardFile)

if options.year >0 : years = [options.year]

bins    = range(46)
years   = [2016,2017,2018]

c = cardFileWriter.cardFileWriter()
c.releaseLocation = os.path.abspath('.')
cardFileName = cardFile.replace('.txt', '_combination.txt') # name of the new card file
c.reset()
c.setPrecision(3)

# nuisances
for sys_tex, sys_name, sys_corr, shape, merge in systematics:
    if merge or len(sys_corr)==1:
        c.addUncertainty(sys_name, shape)
    else:
        for sys_c in sys_corr:
            c.addUncertainty(sys_c, shape)

if options.signal == "ttHinv":
    c.addUncertainty('xsec_QCD',    'lnN')
    c.addUncertainty('xsec_PDF',    'lnN')
else:
    c.addUncertainty('leptonFS',    shape)
    #c.addUncertainty('btagFS',      shape)
    #c.addUncertainty('isr',         shape)
    #c.addUncertainty('FSmet',       shape)
    c.addUncertainty('b_fast',      shape)
    c.addUncertainty('ISR_Weight',  shape)
    c.addUncertainty('MET_Unc',     shape)
    c.addUncertainty('LHESigScale', shape)

# rate parameters
c.addRateParameter('DY',            1, '[0,10]')
#c.addRateParameter('TTZ',         1, '[0,10]')
c.addRateParameter('TTZDL',         1, '[0,10]')
c.addRateParameter('TTJets',        1, '[0,10]')
c.addRateParameter('multiBoson',    1, '[0.6,1.4]')

print cardFile

for b in bins:
    binname = "Bin%s"%b
    #obs     = sum( [ getObservationFromCard(cardFile, "dc_%s_Bin%s"%(year, b)) for year in years ])
    obs     = sum( [ result["dc_%s_Bin%s"%(year, b)]["obs"]["yield"] for year in years ])
    #signal   = sum( [ getEstimateFromCard(cardFile, "signal", "dc_%s_Bin%s"%(year,b)) for year in years ])
    signal   = sum( [ result["dc_%s_Bin%s"%(year,b)]["signal"]["yield"] for year in years ])
    c.addBin(binname, [ np for p,np,t in processes ], binname)
    #c.addBin(binname, [ p for p,t in processes ], binname, noSignal=True)
    #c.addBin(binname, [], binname)
    for proc, newproc, tex in processes + [('signal', 'signal', '')]:
    #for proc, tex in processes:
    #for proc, tex in [('signal','')]:
        est = sum( [ result["dc_%s_Bin%s"%(year,b)][proc]["yield"] for year in years ])
        ests = { year: result["dc_%s_Bin%s"%(year,b)][proc]["yield"] for year in years }
        c.specifyExpectation(binname, newproc, est )
        for sys_tex, sys_name, sys_corr, shape, merge in systematics:
            if len(sys_corr) == 1:
                sigma = sum( [ result["dc_%s_%s"%(year, binname)][proc][sys_corr[0]]*ests[year] for year in years ])
                if est > 0 and sigma > 0:
                    c.specifyUncertainty(sys_name,    binname, newproc, 1 + sigma/est)
            else:
                if merge:
                    listOfSigmas = [ (result["dc_%s_%s"%(year, binname)][proc][sorted(sys_corr)[i]]*ests[year])**2 for i,year in enumerate(years) ]
                    sigma = math.sqrt(sum( listOfSigmas ))
                    if est > 0 and sigma > 0:
                        c.specifyUncertainty(sys_name,    binname, newproc, 1 + sigma/est)
                else:
                    for i, year in enumerate(years):
                        sigma = result["dc_%s_%s"%(year, binname)][proc][sorted(sys_corr)[i]]*ests[year]
                        if est > 0 and sigma > 0:
                            c.specifyUncertainty(sorted(sys_corr)[i],    binname, newproc, 1 + sigma/est)
                    

        # stat uncertainty
        listOfSigmas = [ (result["dc_%s_%s"%(year,binname)][proc]["stat"]*ests[year])**2 for i,year in enumerate(years) ]
        sigma = math.sqrt(sum( listOfSigmas ))
        uname = 'Stat_'+binname+'_'+newproc
        c.addUncertainty(uname, 'lnN')
        if est>0:
            c.specifyUncertainty(uname,    binname, newproc, 1+sigma/est)

    ## OBSERVATION
    c.specifyObservation(binname, int(obs))

    ### SIGNAL
    if options.signal == "ttHinv":
        # x-sec uncertainties for ttH: https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CERNYellowReportPageBSMAt13TeV#ttH_Process
        c.specifyUncertainty('xsec_QCD',      binname, 'signal', 1.092)
        c.specifyUncertainty('xsec_PDF',      binname, 'signal', 1.036)
    elif options.signal == "T2tt":
        if est>0:
            c.specifyUncertainty('LHESigScale',    binname, 'signal', 1 +  sum( [ (result["dc_%s_%s"%(year,binname)][proc]["scale"]*ests[year]) for year in years ])/est )
            c.specifyUncertainty('leptonFS', binname, 'signal', 1 +  sum( [ (result["dc_%s_%s"%(year,binname)][proc]["leptonFS"]*ests[year]) for year in years ])/est )
            c.specifyUncertainty('b_fast',   binname, 'signal', 1 +  sum( [ (result["dc_%s_%s"%(year,binname)][proc]["btagFS"]*ests[year]) for year in years ])/est )
            c.specifyUncertainty('ISR_Weight',      binname, 'signal', 1 +  sum( [ (result["dc_%s_%s"%(year,binname)][proc]["isr"]*ests[year]) for year in years ])/est )
            c.specifyUncertainty('MET_Unc',    binname, 'signal', 1 +  sum( [ (result["dc_%s_%s"%(year,binname)][proc]["FSmet"]*ests[year]) for year in years ])/est )

#c.writeToFile(cardFileName)
c.writeToShapeFile(cardFileName.replace('.txt', '_shape.root'))

calcLimit = True

if calcLimit:
    import pandas
    resFile = "calculatedLimits.pkl"
    if not os.path.isfile(resFile):
        pickle.dump([], file(resFile, 'w'))
    
    results = pickle.load(file(resFile))
    found = False
    if len(results)>0:
        df = pandas.DataFrame(results)
        selection = ((df['spin']==spin) & (df['mStop']==mStop) & (df['mLSP']==mLSP))
        if not df[selection].empty:
            print df[selection]
            found = True

    if not found:
        res = c.calcLimit(cardFileName.replace('.txt', '_shapeCard.txt'))
        res['spin']     = spin
        res['mStop']    = mStop
        res['mLSP']     = mLSP
        results.append(res)
        pickle.dump(results, file(resFile, 'w'))

print "Card file:", cardFileName.replace('.txt', '_shapeCard.txt')

