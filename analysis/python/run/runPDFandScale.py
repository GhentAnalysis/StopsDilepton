'''
Extraction of PDF and scale uncertainties in the SRs
'''

#!/usr/bin/env python
from optparse import OptionParser
parser = OptionParser()
#parser.add_option("--noMultiThreading",     dest="noMultiThreading",      default = False,             action="store_true", help="noMultiThreading?")
parser.add_option("--selectWeight",         dest="selectWeight",       default=None,                action="store",      help="select weight?")
parser.add_option("--selectRegion",         dest="selectRegion",          default=None, type="int",    action="store",      help="select region?")
parser.add_option("--sample",               dest='sample',  action='store', default='TTZ',    choices=["TTZ", "DY", "multiBoson", "TTJets", "other"], help="which sample?")
parser.add_option("--year",                 dest='year',  action='store', default='2016',    choices=["2016", "2017", "2018"], help="which year?")
parser.add_option("--small",                action='store_true', help="small?")
parser.add_option("--reducedPDF",           action='store_true', help="Don't use all PDF variations for tests?")
parser.add_option("--combine",              action='store_true', help="Combine results?")
parser.add_option("--noKeepNorm",           action='store_true', help="Keep the normalization = acceptance uncertainty only?")
parser.add_option("--signal",               action='store', default=None, help="Which signal, if any?")
parser.add_option("--only",                 action='store', default=None, help="pick only one masspoint?")
parser.add_option('--logLevel',             dest="logLevel",              default='INFO',              action='store',      help="log level?", choices=['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'TRACE', 'NOTSET'])
parser.add_option('--overwrite',            dest="overwrite", default = False, action = "store_true", help="Overwrite existing output files, bool flag set to True  if used")
parser.add_option('--skipCentral',          dest="skipCentral", default = False, action = "store_true", help="Skip central weights")
parser.add_option('--nJobs',                dest="nJobs", default=1, type="int", action = "store", help="How many jobs?")
parser.add_option('--job',                  dest="job", default=0, type="int", action = "store", help="Which job?")
parser.add_option('--dpm',                  dest='dpm',         default=False,      action='store_true', help='Use dpm?')
(options, args) = parser.parse_args()

# Logging
import StopsDilepton.tools.logger as logger
logger = logger.get_logger(options.logLevel, logFile = None )
import RootTools.core.logger as logger_rt
logger_rt = logger_rt.get_logger(options.logLevel, logFile = None )

# Standard imports
import ROOT
import os
import sys
import pickle
import math

#RootTools
from RootTools.core.standard import *

from StopsDilepton.analysis.SetupHelpers import channels, allChannels, trilepChannels
from StopsDilepton.analysis.estimators   import *
from StopsDilepton.analysis.regions      import regionsLegacy, noRegions
from StopsDilepton.tools.resultsDB       import resultsDB
from Analysis.Tools.u_float              import u_float

# use this for job splitting
from RootTools.core.helpers import partition

if options.dpm:
    data_directory = "/dpm/oeaw.ac.at/home/cms/store/user/rschoefbeck/Stops2l-postprocessed/"


from StopsDilepton.analysis.Setup import Setup

year        = int(options.year)

# dummy values for now
PDFset = 'NNPDF30' # keep this for now, it has no physical impact although 2017 and 2018 LHAPDFs are used
PDFType = "hessian" if (year == 2017 or year == 2018 or options.signal=='TTbarDM') else "replicas"
#PDFType = "hessian"
PSweights = False

print "PDF type:", PDFType

#from StopsDilepton.samples.nanoTuples_Summer16_postProcessed            import Top_pow_16, DY_HT_LO_16, TTZ_16, multiBoson_16,TTXNoZ_16
#from StopsDilepton.samples.nanoTuples_Fall17_postProcessed              import Top_pow_17, DY_HT_LO_17, TTZ_17, multiBoson_17, TTXNoZ_17
#from StopsDilepton.samples.nanoTuples_Autumn18_postProcessed            import Top_pow_18, DY_HT_LO_18, TTZ_18, multiBoson_18, TTXNoZ_18

setupSR     = Setup(year=year)

setupDYVV   = setupSR.sysClone(parameters={'nBTags':(0,0 ), 'dPhi': False, 'dPhiInv': False,  'zWindow': 'onZ', 'metSigMin' : 12})
setupTTZ1   = setupSR.sysClone(parameters={'triLep': True, 'zWindow' : 'onZ', 'mllMin': 0, 'metMin' : 0, 'metSigMin' : 0, 'nJets':(2,2),  'nBTags':(2,-1), 'dPhi': False, 'dPhiInv': False})
setupTTZ2   = setupSR.sysClone(parameters={'triLep': True, 'zWindow' : 'onZ', 'mllMin': 0, 'metMin' : 0, 'metSigMin' : 0, 'nJets':(3,3),  'nBTags':(1,1),  'dPhi': False, 'dPhiInv': False})
setupTTZ3   = setupSR.sysClone(parameters={'triLep': True, 'zWindow' : 'onZ', 'mllMin': 0, 'metMin' : 0, 'metSigMin' : 0, 'nJets':(3,3),  'nBTags':(2,-1), 'dPhi': False, 'dPhiInv': False})
setupTTZ4   = setupSR.sysClone(parameters={'triLep': True, 'zWindow' : 'onZ', 'mllMin': 0, 'metMin' : 0, 'metSigMin' : 0, 'nJets':(4,-1), 'nBTags':(1,1),  'dPhi': False, 'dPhiInv': False})
setupTTZ5   = setupSR.sysClone(parameters={'triLep': True, 'zWindow' : 'onZ', 'mllMin': 0, 'metMin' : 0, 'metSigMin' : 0, 'nJets':(4,-1), 'nBTags':(2,-1), 'dPhi': False, 'dPhiInv': False})
setupTT     = setupSR.sysClone()

setupIncl   = setupSR.sysClone(parameters={'mllMin':0, 'nJets':(0,-1), 'nBTags':(0,-1), 'zWindow':'allZ', 'metSigMin':0})

if options.signal:
    if options.signal == 'T2tt':
        if year == 2016:
            data_directory              = '/afs/hephy.at/data/cms06/nanoTuples/'
            postProcessing_directory    = 'stops_2016_nano_v0p23/dilep/'
            from StopsDilepton.samples.nanoTuples_FastSim_Summer16_postProcessed import signals_T2tt as jobs
        elif year == 2017:
            data_directory              = '/afs/hephy.at/data/cms07/nanoTuples/'
            postProcessing_directory    = 'stops_2017_nano_v0p22/dilep/'
            from StopsDilepton.samples.nanoTuples_FastSim_Fall17_postProcessed import signals_T2tt as jobs
        if year == 2018:
            data_directory              = '/afs/hephy.at/data/cms07/nanoTuples/'
            postProcessing_directory    = 'stops_2018_nano_v0p21/dilep/'
            from StopsDilepton.samples.nanoTuples_FastSim_Autumn18_postProcessed import signals_T2tt as jobs

    elif options.signal == 'T2bW':
        if year == 2016:
            data_directory              = '/afs/hephy.at/data/cms02/nanoTuples/'
            postProcessing_directory    = 'stops_2016_nano_v0p19/dilep/'
            from StopsDilepton.samples.nanoTuples_FastSim_Summer16_postProcessed import signals_T2bW as jobs
        elif year == 2017:
            data_directory              = '/afs/hephy.at/data/cms01/nanoTuples/'
            postProcessing_directory    = 'stops_2017_nano_v0p19/dilep/'
            from StopsDilepton.samples.nanoTuples_FastSim_Fall17_postProcessed import signals_T2bW as jobs
        if year == 2018:
            data_directory              = '/afs/hephy.at/data/cms02/nanoTuples/'
            postProcessing_directory    = 'stops_2018_nano_v0p19/dilep/'
            from StopsDilepton.samples.nanoTuples_FastSim_Autumn18_postProcessed import signals_T2bW as jobs

    elif options.signal == 'T8bbllnunu_XCha0p5_XSlep0p05':
        if year == 2016:
            data_directory              = '/afs/hephy.at/data/cms02/nanoTuples/'
            postProcessing_directory    = 'stops_2016_nano_v0p19/dilep/'
            from StopsDilepton.samples.nanoTuples_FastSim_Summer16_postProcessed import signals_T8bbllnunu_XCha0p5_XSlep0p05 as jobs
        elif year == 2017:
            data_directory              = '/afs/hephy.at/data/cms01/nanoTuples/'
            postProcessing_directory    = 'stops_2017_nano_v0p19/dilep/'
            from StopsDilepton.samples.nanoTuples_FastSim_Fall17_postProcessed import signals_T8bbllnunu_XCha0p5_XSlep0p05 as jobs
        if year == 2018:
            data_directory              = '/afs/hephy.at/data/cms02/nanoTuples/'
            postProcessing_directory    = 'stops_2018_nano_v0p19/dilep/'
            from StopsDilepton.samples.nanoTuples_FastSim_Autumn18_postProcessed import signals_T8bbllnunu_XCha0p5_XSlep0p05 as jobs

    elif options.signal == 'T8bbllnunu_XCha0p5_XSlep0p5':
        if year == 2016:
            data_directory              = '/afs/hephy.at/data/cms02/nanoTuples/'
            postProcessing_directory    = 'stops_2016_nano_v0p19/dilep/'
            from StopsDilepton.samples.nanoTuples_FastSim_Summer16_postProcessed import signals_T8bbllnunu_XCha0p5_XSlep0p5 as jobs
        elif year == 2017:
            data_directory              = '/afs/hephy.at/data/cms01/nanoTuples/'
            postProcessing_directory    = 'stops_2017_nano_v0p19/dilep/'
            from StopsDilepton.samples.nanoTuples_FastSim_Fall17_postProcessed import signals_T8bbllnunu_XCha0p5_XSlep0p5 as jobs
        if year == 2018:
            data_directory              = '/afs/hephy.at/data/cms02/nanoTuples/'
            postProcessing_directory    = 'stops_2018_nano_v0p19/dilep/'
            from StopsDilepton.samples.nanoTuples_FastSim_Autumn18_postProcessed import signals_T8bbllnunu_XCha0p5_XSlep0p5 as jobs

    elif options.signal == 'T8bbllnunu_XCha0p5_XSlep0p95':
        if year == 2016:
            data_directory              = '/afs/hephy.at/data/cms02/nanoTuples/'
            postProcessing_directory    = 'stops_2016_nano_v0p19/dilep/'
            from StopsDilepton.samples.nanoTuples_FastSim_Summer16_postProcessed import signals_T8bbllnunu_XCha0p5_XSlep0p95 as jobs
        elif year == 2017:
            data_directory              = '/afs/hephy.at/data/cms01/nanoTuples/'
            postProcessing_directory    = 'stops_2017_nano_v0p19/dilep/'
            from StopsDilepton.samples.nanoTuples_FastSim_Fall17_postProcessed import signals_T8bbllnunu_XCha0p5_XSlep0p95 as jobs
        if year == 2018:
            data_directory              = '/afs/hephy.at/data/cms02/nanoTuples/'
            postProcessing_directory    = 'stops_2018_nano_v0p19/dilep/'
            from StopsDilepton.samples.nanoTuples_FastSim_Autumn18_postProcessed import signals_T8bbllnunu_XCha0p5_XSlep0p95 as jobs

    elif options.signal == 'ttHinv':
        if year == 2016:
            data_directory              = '/afs/hephy.at/data/cms03/nanoTuples/'
            postProcessing_directory    = 'stops_2016_nano_v0p26/dilep/'
        elif year == 2017:
            data_directory              = '/afs/hephy.at/data/cms09/nanoTuples/'
            postProcessing_directory    = 'stops_2017_nano_v0p22/dilep/'
        elif year == 2018:
            data_directory              = '/afs/hephy.at/data/cms09/nanoTuples/'
            postProcessing_directory    = 'stops_2018_nano_v0p22/dilep/'
        ttH_HToInvisible_M125 = Sample.fromDirectory(name="ttH_HToInvisible_M125", treeName="Events", isData=False, color=1, texName="ttH(125)", directory=os.path.join(data_directory,postProcessing_directory,'ttH_HToInvisible'))
        jobs = [ttH_HToInvisible_M125]

    elif options.signal == 'TTbarDM':
        if year == 2016:
            data_directory              = '/afs/hephy.at/data/cms03/nanoTuples/'
            postProcessing_directory    = 'stops_2016_nano_v0p26/dilep/'
            from StopsDilepton.samples.nanoTuples_Summer16_TTDM_postProcessed import signals as jobs
        elif year == 2017:
            data_directory              = '/afs/hephy.at/data/cms03/nanoTuples/'
            postProcessing_directory    = 'stops_2017_nano_v0p26/dilep/'
            from StopsDilepton.samples.nanoTuples_Fall17_TTDM_postProcessed import signals as jobs
        elif year == 2018:
            data_directory              = '/afs/hephy.at/data/cms03/nanoTuples/'
            postProcessing_directory    = 'stops_2018_nano_v0p26/dilep/'
            from StopsDilepton.samples.nanoTuples_Autumn18_TTDM_postProcessed import signals as jobs


    if options.only.isdigit():
        sample = jobs[int(options.only)]
    else:
        jobNames = [ x.name for x in jobs ]
        #print jobNames
        sample = jobs[jobNames.index(options.only)]

else:
    sample = setupSR.samples[options.sample]

logger.info("Sample: %s", sample.name)

allRegions = regionsLegacy
#if options.aggregate: allRegions = regionsAgg

from StopsDilepton.analysis.MCBasedEstimate import MCBasedEstimate
from StopsDilepton.analysis.DataObservation import DataObservation


'''
PDF4LHC variations are hessian eigenvectors
NNPDF30 variations are usually mc replicas (unless otherwise stated)

## Summer16 ##
LHE scale variation weights (w_var / w_nominal); [0] is muR=0.50000E+00 muF=0.50000E+00 ; [1] is muR=0.50000E+00 muF=0.10000E+01 ; [2] is muR=0.50000E+00 muF=0.20000E+01 ; [3] is muR=0.10000E+01 muF=0.50000E+00 ; [4] is muR=0.10000E+01 muF=0.10000E+01 ; [5] is muR=0.10000E+01 muF=0.20000E+01 ; [6] is muR=0.20000E+01 muF=0.50000E+00 ; [7] is muR=0.20000E+01 muF=0.10000E+01 ; [8] is muR=0.20000E+01 muF=0.20000E+01
Therefore, use 0, 1,3,4 (central?), 5, 7, 8
Used PDF sets
TTLep_pow: 260000 (NNPDF30_nlo_as_0118)
TTZToLLNuNu: 292200 (NNPDF30_nlo_nf_5_pdfas)
LowMass TTZ: 262000 (NNPDF30_lo_as_0118) -> sample could go to ttX anyway (nunu not included)
SUSY: no weights stored. damn

## Fall17 ##
LHE scale variation weights (w_var / w_nominal); [0] is renscfact=0.5d0 facscfact=0.5d0 ; [1] is renscfact=0.5d0 facscfact=1d0 ; [2] is renscfact=0.5d0 facscfact=2d0 ; [3] is renscfact=1d0 facscfact=0.5d0 ; [4] is renscfact=1d0 facscfact=1d0 ; [5] is renscfact=1d0 facscfact=2d0 ; [6] is renscfact=2d0 facscfact=0.5d0 ; [7] is renscfact=2d0 facscfact=1d0 ; [8] is renscfact=2d0 facscfact=2d0
Therefore, use 0, 1,3,4 (central?), 5, 7, 8
Used PDF sets
TTLep_pow: 91400 (PDF4LHC15_nnlo_30_pdfas)
TTZToLLNuNu: 91400 (PDF4LHC15_nnlo_30_pdfas)
LowMass TTZ: 91400
SUSY: no weights stored atm.

Autumn18:
LHE scale variation weights (w_var / w_nominal); [0] is renscfact=0.5d0 facscfact=0.5d0 ; [1] is renscfact=0.5d0 facscfact=1d0 ; [2] is renscfact=0.5d0 facscfact=2d0 ; [3] is renscfact=1d0 facscfact=0.5d0 ; [4] is renscfact=1d0 facscfact=1d0 ; [5] is renscfact=1d0 facscfact=2d0 ; [6] is renscfact=2d0 facscfact=0.5d0 ; [7] is renscfact=2d0 facscfact=1d0 ; [8] is renscfact=2d0 facscfact=2d0 
Therefore, use 0, 1,3,4 (central?), 5, 7, 8
Used PDF sets
TTLep_pow: 91400 (PDF4LHC15_nnlo_30_pdfas)
TTZToLLNuNu: 91400 (PDF4LHC15_nnlo_30_pdfas)
LowMass TTZ: 91400
SUSY: no weights stored atm.
'''

scale_indices = [0,1,3,4,5,7,8]
if options.signal == 'ttHinv':
    LHEweight_original = 'abs(LHEScaleWeight[4])'
elif options.signal == 'TTbarDM':
    scale_indices = [0,1,2,3,4,6,8] # skip extremes (indices 5 and 7)
    LHEweight_original = 'abs(LHE_weight[0])'
elif options.signal.startswith('T'):
    LHEweight_original = 'abs(LHE_weight[4])'
else:
    LHEweight_original = 'abs(LHEScaleWeight[4])'

if options.signal.startswith('T'):
    LHEweight_original_PDF = LHEweight_original
else:
    LHEweight_original_PDF = 'abs(LHEPdfWeight[0])' if (year == 2017 or year == 2018) else LHEweight_original

centralWeight = LHEweight_original

pdf_indices = range(100) if year == 2016 else range(30)

if options.signal == 'TTbarDM':
    pdf_indices = range(9,110) # 9 is central
elif options.signal.startswith('T'):
    pdf_indices = []

'''
alpha_S indices TTDM: 110, 111
PS indices in TTDM samples:
- 112 central
- 113,114,115,116 reduced: isrRedHi,fsrRedHi,isrRedLo,fsrRedLo
- 117,118,119,120 default
- 121,122,123,124 conservative
'''

#if year == 2016:
#    if options.sample == 'TTLep_pow': #only use ttbar sample, no single-t
#        PDF_indices = range(100) # no central weight stored
#        aS_indices = []
#    elif options.sample == 'TTZToLLNuNu':
#        PDF_indices = range(100)
#        aS_indices = [100,101]
#    #elif options.sample.contains('T2tt'):
#    #    raise NotImplementedError

if not options.selectWeight:
    centralPDFIndex = 0 if  not options.signal == 'TTbarDM' else 9
    if options.signal == 'ttHinv':
        scaleWeightString   = 'LHEScaleWeight'
    elif options.signal.startswith('T'):
        scaleWeightString   ='LHE_weight' # also applies to TTbarDM
    else:
        scaleWeightString   = 'LHEScaleWeight'

    scale_variations    = [ "abs(%s[%s])"%(scaleWeightString, str(i)) for i in scale_indices ]
    pdfWeightString     = 'LHEPdfWeight' if not options.signal == 'TTbarDM' else 'LHE_weight'
    if year == 2016 and not options.signal == 'TTbarDM':
        PDF_variations      = [ "abs(%s[%s])"%(pdfWeightString, str(i)) for i in pdf_indices ]
    else:
        if not options.signal == 'TTbarDM':
            PDF_variations      = [ "(abs(%s[%s])/abs(%s[%s]))"%(pdfWeightString, str(i), pdfWeightString, str(centralPDFIndex)) for i in pdf_indices ]
        else:
            PDF_variations      = [ "(%s[%s])"%(pdfWeightString, str(i)) for i in pdf_indices ]

    aS_variations       = [] #[ "abs(LHEPdfWeight[100])", "abs(LHEPdfWeight[101])"] if year == 2016 else [ "abs(LHEPdfWeight[31])", "abs(LHEPdfWeight[32])"]
    if options.signal == 'TTbarDM':
        aS_variations       = [ "(%s[%s])"%(pdfWeightString, str(i)) for i in [110,111] ]

    if options.signal == 'ttHinv':
        variations          = scale_variations + PDF_variations + ['(1)'] 
    elif options.signal.startswith('T2') or options.signal.startswith('T8'):
        variations          = scale_variations + ['(1)']
    else:
        variations          = scale_variations + PDF_variations + ['(1)'] + aS_variations

#print variations

# only properly works for selectRegion>0
selectRegion = True if options.selectRegion >= 0 else False
regions = allRegions if not selectRegion else [allRegions[options.selectRegion]]


results = {}

scale_systematics = {}

cacheDir = "/afs/hephy.at/data/cms05/StopsDileptonLegacy/results/PDF_v2_%s/%s/"%(PDFset,year)
#cacheDir = "/afs/hephy.at/data/cms05/StopsDileptonLegacy/results/PDF_v3_%s/%s/"%(PDFset,year) # v2 is used for pre-approval results. some caches got lost, so rerunning in v3

estimate = MCBasedEstimate(name=sample.name, sample=sample )
estimate.initCache(cacheDir)

## Results DB for scale and PDF uncertainties

PDF_cache = resultsDB(cacheDir+'PDFandScale_unc.sq', "PDF", ["name", "region", "CR", "channel", "PDFset"])
scale_cache = resultsDB(cacheDir+'PDFandScale_unc.sq', "scale", ["name", "region", "CR", "channel", "PDFset"])
PS_cache = resultsDB(cacheDir+'PDFandScale_unc.sq', "PSscale", ["name", "region", "CR", "channel", "PDFset"])

print cacheDir+'PDFandScale_unc.sq'

def wrapper(args):
        r, c, s = args
        res = estimate.cachedEstimate(r, c, s)
        logger.debug("Done with one of the jobs in region %s and channel %s", r, c)
        return (estimate.uniqueKey(r, c, s), res )


setupSR.regions     = regionsLegacy[1:]
setupDYVV.regions   = regionsLegacy[1:]
setupTTZ1.regions   = noRegions
setupTTZ2.regions   = noRegions
setupTTZ3.regions   = noRegions
setupTTZ4.regions   = noRegions
setupTTZ5.regions   = noRegions
setupTT.regions     = [regionsLegacy[0]]

setups = [setupSR,setupTTZ1,setupTTZ2,setupTTZ3,setupTTZ4,setupTTZ5,setupTT,setupDYVV]

jobs=[]

if not options.skipCentral:
    # First run over seperate channels
    jobs.append((noRegions[0], 'all', setupIncl))
    jobs.append((noRegions[0], 'all', setupIncl.sysClone(sys={'reweight':[LHEweight_original]})))
    jobs.append((noRegions[0], 'all', setupIncl.sysClone(sys={'reweight':[PDF_variations[0]]})))
    jobs.append((noRegions[0], 'all', setupIncl.sysClone(sys={'reweight':[LHEweight_original_PDF]})))
    #jobs.append((noRegions[0], 'all', setupIncl.sysClone(sys={'reweight':['(1)']})))
    for var in variations:
        #for c in ['EE', 'MuMu', 'EMu']:
        jobs.append((noRegions[0], 'all', setupIncl.sysClone(sys={'reweight':[var]})))


if not options.combine:
    for region in regions:
        logger.info("Queuing jobs for region %s", region)
        for c in ['EE', 'MuMu', 'EMu']:
            logger.info("Queuing jobs for channel %s", c)
            jobs.append((region, c, setupSR))
            jobs.append((region, c, setupTT))
            jobs.append((region, c, setupSR.sysClone(sys={'reweight':[LHEweight_original]})))
            jobs.append((region, c, setupSR.sysClone(sys={'reweight':[PDF_variations[0]]})))
            jobs.append((region, c, setupTT.sysClone(sys={'reweight':[LHEweight_original]})))
            jobs.append((region, c, setupTT.sysClone(sys={'reweight':[PDF_variations[0]]})))
            #jobs.append((region, c, setupTT.sysClone(sys={'reweight':['(1)']})))
            if not c == 'EMu':
                jobs.append((region, c, setupDYVV))
                jobs.append((region, c, setupDYVV.sysClone(sys={'reweight':[LHEweight_original]})))
                jobs.append((region, c, setupDYVV.sysClone(sys={'reweight':[PDF_variations[0]]})))
            for var in variations:
                jobs.append((region, c, setupSR.sysClone(sys={'reweight':[var]})))
                jobs.append((region, c, setupTT.sysClone(sys={'reweight':[var]})))
                if not c == 'MuMu': # use a rendom c
                    jobs.append((region, 'SF', setupDYVV.sysClone(sys={'reweight':[var]})))
                #sigma_reweight  = estimate.cachedEstimate(region, c, setup.sysClone(sys={'reweight':[var]}))
    
    # how to not run this for every --selectRegion?
    if ( selectRegion and options.selectRegion == 0 ) or len(regions)>1:
        for c in ['3mu', '2mu1e', '2e1mu', '3e']:
            for setup in [setupTTZ1,setupTTZ2,setupTTZ3,setupTTZ4,setupTTZ5]:
                jobs.append((noRegions[0], c, setup))
                jobs.append((noRegions[0], c, setup.sysClone(sys={'reweight':[LHEweight_original]})))
                jobs.append((noRegions[0], c, setup.sysClone(sys={'reweight':[PDF_variations[0]]})))
                #jobs.append((noRegions[0], c, setup.sysClone(sys={'reweight':['(1)']})))
                for var in variations:
                    jobs.append((noRegions[0], c, setup.sysClone(sys={'reweight':[var]})))

    ## PDF central weights (do last)
    for var in variations:
        #jobs.append((noRegions[0], 'all', setupSR.sysClone(sys={'reweight':[var]})))
        #jobs.append((noRegions[0], 'all', setupTT.sysClone(sys={'reweight':[var]})))
        for setup in [setupTTZ1,setupTTZ2,setupTTZ3,setupTTZ4,setupTTZ5, setupSR, setupTT]:
            jobs.append((noRegions[0], 'all', setup.sysClone(sys={'reweight':[var]})))
            

    for region in regions:
        logger.info("Queuing PDF jobs for region %s", region)
        for var in variations:
            # combine the individual channels here and only now to avoid running the same job twice at the same time
            jobs.append((region, 'all', setupSR.sysClone(sys={'reweight':[var]})))
        for c in ['EE', 'MuMu', 'EMu']:
            logger.info("Queuing jobs for channel %s", c)
            jobs.append((region, c, setupSR.sysClone(sys={'reweight':[LHEweight_original_PDF]})))
            jobs.append((region, c, setupTT.sysClone(sys={'reweight':[LHEweight_original_PDF]})))
            #jobs.append((region, c, setupSR.sysClone(sys={'reweight':['(1)']})))
            #jobs.append((region, c, setupTT.sysClone(sys={'reweight':['(1)']})))
            if not c == 'EMu':
                jobs.append((region, c, setupDYVV.sysClone(sys={'reweight':[LHEweight_original_PDF]})))
                #jobs.append((region, c, setupDYVV.sysClone(sys={'reweight':['(1)']})))
    
    if ( selectRegion and options.selectRegion == 0 ) or len(regions)>1:
        for c in ['3mu', '2mu1e', '2e1mu', '3e']:
            for setup in [setupTTZ1,setupTTZ2,setupTTZ3,setupTTZ4,setupTTZ5]:
                jobs.append((noRegions[0], c, setup.sysClone(sys={'reweight':[LHEweight_original_PDF]})))
                #jobs.append((noRegions[0], c, setup.sysClone(sys={'reweight':['(1)']})))
    

    logger.info("Created %s jobs",len(jobs))

    jobs = partition(jobs, options.nJobs)[options.job]

    logger.info("Running over %s jobs", len(jobs))

    results = map(wrapper, jobs)

    ## multithreading doesn't work with DirDB/MergeDirDB ...
    #if options.noMultiThreading: 
    #    results = map(wrapper, jobs)
    #else:
    #    from multiprocessing import Pool
    #    pool = Pool(processes=1)
    #    results = pool.map(wrapper, jobs)
    #    pool.close()
    #    pool.join()
    
    logger.info("All done.")

setupSR.channels   = ['SF','EMu']
setupDYVV.channels = ['SF']
setupTTZ1.channels = ['all']
setupTTZ2.channels = ['all']
setupTTZ3.channels = ['all']
setupTTZ4.channels = ['all']
setupTTZ5.channels = ['all']
setupTT.channels = ['SF','EMu']


PDF_unc     = []
Scale_unc   = []
PS_unc      = []

if options.combine:
    for setup in setups:
        for c in setup.channels:#allChannels:
        
            for region in setup.regions:
                logger.info("Region: %s", region)
                
                scales = []
                showerScales = []
                deltas = []
                delta_squared = 0
                # central yield inclusive and in region
                logger.debug("Getting inclusive (noRegions) yield")
                sigma_incl_central  = estimate.cachedEstimate(noRegions[0], 'all', setupIncl.sysClone(sys={'reweight':[LHEweight_original]}))
                logger.debug("Getting yield for region with LHEweight_original")
                sigma_central       = estimate.cachedEstimate(region, c, setup.sysClone(sys={'reweight':[LHEweight_original]}))
                logger.debug("Getting yield for region with centralWeight")
                sigma_centralWeight = estimate.cachedEstimate(region, c, setup.sysClone(sys={'reweight':[centralWeight]}))

                for var in scale_variations:
                    #print var
                    logger.debug("Getting inclusive yield with varied weight")
                    simga_incl_reweight = estimate.cachedEstimate(noRegions[0], 'all', setupIncl.sysClone(sys={'reweight':[var]}))
                    norm = sigma_incl_central/simga_incl_reweight if not options.noKeepNorm else 1
                    
                    logger.debug("Getting yield for region with varied weight")
                    sigma_reweight  = estimate.cachedEstimate(region, c, setup.sysClone(sys={'reweight':[var]}))
                    sigma_reweight_acc = sigma_reweight * norm
                    
                    #logger.info("Using norm of %s", norm)
                    unc = abs( ( sigma_reweight_acc - sigma_central) / sigma_central ) if sigma_central > 0 else u_float(1)
                    scales.append(unc.val)
                
                scale_rel = max(scales)

                ## PDF stuff
                #sigma_incl_central_PDF  = estimate.cachedEstimate(noRegions[0], 'all', setupIncl.sysClone(sys={'reweight':[LHEweight_original_PDF]}))
                logger.debug("sigma_incl_central_PDF")
                sigma_incl_central_PDF  = estimate.cachedEstimate(noRegions[0], 'all', setupIncl.sysClone(sys={'reweight':[PDF_variations[0]]})) # changed from (1)
                logger.debug("sigma_centralWeight_PDF")
                sigma_centralWeight_PDF = estimate.cachedEstimate(region, c, setup.sysClone(sys={'reweight':[PDF_variations[0]]}))
                for i, var in enumerate(PDF_variations):
                    # calculate x-sec noramlization
                    logger.debug("simga_incl_reweight")
                    simga_incl_reweight = estimate.cachedEstimate(noRegions[0], 'all', setupIncl.sysClone(sys={'reweight':[var]}))
                    norm = sigma_incl_central_PDF/simga_incl_reweight if not options.noKeepNorm else 1

                    logger.debug("simga_reweight")
                    sigma_reweight  = estimate.cachedEstimate(region, c, setup.sysClone(sys={'reweight':[var]}))
                    sigma_reweight_acc = sigma_reweight * norm

                    ## For replicas, just get a list of all sigmas, sort it and then get the 68% interval
                    deltas.append(sigma_reweight_acc.val)
                    ## recommendation for hessian is to have delta_sigma = sum_k=1_N( (sigma_k - sigma_0)**2 )
                    ## so I keep the norm for both sigma_k and sigma_0 to obtain the acceptance uncertainty. Correct?
                    delta_squared += ( sigma_reweight.val - sigma_centralWeight_PDF.val )**2
                    #print i, delta_squared
                    #print sigma_reweight.val, sigma_centralWeight_PDF.val
                
                deltas = sorted(deltas)

                # calculate uncertainty
                delta_sigma = 0
                if PDFType == "replicas" and len(deltas)>0:
                    # get the 68% interval
                    upper = len(deltas)*84/100-1
                    lower = len(deltas)*16/100 - 1
                    delta_sigma = (deltas[upper]-deltas[lower])/2
                elif PDFType == "hessian":
                    delta_sigma = math.sqrt(delta_squared)

                # recommendation is to multiply uncertainty by 1.5
                deltas_as = []
                for var in aS_variations:
                    simga_incl_reweight = estimate.cachedEstimate(noRegions[0], 'all', setupIncl.sysClone(sys={'reweight':[var]}))
                    norm = sigma_incl_central_PDF/simga_incl_reweight if not options.noKeepNorm else 1
                    
                    sigma_reweight  = estimate.cachedEstimate(region, c, setup.sysClone(sys={'reweight':[var]}))
                    sigma_reweight_acc = sigma_reweight * norm

                    deltas_as.append(sigma_reweight_acc.val)

                if len(deltas_as)>0:
                    scale = 1.5 if PDFset.count("NNPDF") else 1.0
                    delta_sigma_alphaS = scale * ( deltas_as[0] - deltas_as[1] ) / 2.

                    # add alpha_s and PDF in quadrature
                    delta_sigma_total = math.sqrt( delta_sigma_alphaS**2 + delta_sigma**2 )
                    try:
                        delta_sigma_alphaS_rel = delta_sigma_alphaS/sigma_central.val
                    except:
                        0.001
                else:
                    delta_sigma_total = delta_sigma

                # make it relative wrt central value in region
                try:
                    delta_sigma_rel = delta_sigma/sigma_central.val
                except:
                    delta_sigma_rel = 0.001 # eh wurscht

                if delta_sigma_rel > 1: print "############# ALERTA #################", delta_sigma_rel

                # calculate the PS uncertainties
                if PSweights:
                    sigma_incl_central  = estimate.cachedEstimate(noRegions[0], 'all', setupIncl.sysClone(sys={'reweight':[PSweight_original]}))
                    sigma_central       = estimate.cachedEstimate(region, c, setup.sysClone(sys={'reweight':[PSweight_original]}))
                    print "Count:", sigma_central.val/0.5148500
                    shower_scales  = []
                    for var in PS_variations:
                        simga_incl_reweight = estimate.cachedEstimate(noRegions[0], channel(-1,-1), setupIncl.sysClone(sys={'reweight':[var]}))
                        norm = sigma_incl_central/simga_incl_reweight
                        
                        sigma_reweight  = estimate.cachedEstimate(region, c, setup.sysClone(sys={'reweight':[var]}))
                        sigma_reweight_acc = sigma_reweight #* norm

                        #unc = ( ( sigma_reweight_acc - sigma_central) / sigma_central ) # no abs atm
                        unc = sigma_reweight_acc / sigma_central
                        #print ( sigma_reweight_acc - sigma_central) / sigma_central
                        showerScales.append(unc.val)

                    print "ISR up/down", round(showerScales[0], 3), round( showerScales[2], 3)
                    print "FSR up/down", round(showerScales[1], 3), round( showerScales[3], 3)

                    PS_scale_rel = max(showerScales)
                else:
                    PS_scale_rel = 0.

                niceName = ' '.join([c, region.__str__()])
                if setup == setupDYVV: niceName += "_controlDYVV"
                if setup == setupTTZ1: niceName += "_controlTTZ1"
                if setup == setupTTZ2: niceName += "_controlTTZ2"
                if setup == setupTTZ3: niceName += "_controlTTZ3"
                if setup == setupTTZ4: niceName += "_controlTTZ4"
                if setup == setupTTZ5: niceName += "_controlTTZ5"
                if setup == setupTT:   niceName += "_controlTTBar"

                logger.info("Calculated PDF and alpha_s uncertainties for region %s in channel %s"%(region, c))
                logger.info("Central x-sec: %s", sigma_central)
                logger.info("Delta x-sec using PDF variations: %s", delta_sigma)
                logger.info("Delta x-sec using alpha_S variations: %s", delta_sigma_alphaS)
                logger.info("Delta x-sec total: %s", delta_sigma_total)
                logger.info("Relative PDF uncertainty: %s", delta_sigma_rel)
                logger.info("Relative alpha_S uncertainty: %s", delta_sigma_alphaS_rel)
                logger.info("Relative scale uncertainty: %s", scale_rel)
                #logger.info("Relative shower scale uncertainty: %s", PS_scale_rel)
                
                if sigma_central.val>0:
                    #print "stat uncertainty", sigma_central.sigma/sigma_central.val
                    #print "scale uncertainty", scale_rel
                    if sigma_central.sigma/sigma_central.val < scale_rel or True:
                        Scale_unc.append(scale_rel)
                    #Scale_unc.append(scale_rel)
                    if sigma_central.sigma/sigma_central.val < 0.50: # 0.15
                        PDF_unc.append(delta_sigma_rel)
                        #if scale_rel < 1: Scale_unc.append(scale_rel) # only append here if we have enough stats
                #PS_unc.append(PS_scale_rel)
                
                # Store results

                if not options.reducedPDF:
                    PDF_cache.add({"name": sample.name, "region":region, "CR":niceName, "channel":c, "PDFset":PDFset}, delta_sigma_rel, overwrite=True)
                scale_cache.add({"name": sample.name, "region":region, "CR":niceName, "channel":c, "PDFset":'scale'}, scale_rel, overwrite=True)
                #PS_cache.add({"region":region, "channel":c, "PDFset":'PSscale'}, PS_scale_rel, overwrite=True)

                if not options.reducedPDF:
                    PDF_cache.get({"region":region, "channel":c, "PDFset":PDFset})
                scale_cache.get({"name": sample.name, "region":region, "CR":niceName, "channel":c, "PDFset":'scale'})
                #PS_cache.get({"region":region, "channel":c, "PDFset":'PSscale'})

    cleanPDF = PDF_unc #[ x for x in PDF_unc if x<1 ]

    logger.info('Min. PDF uncertainty: %.3f', min(cleanPDF))
    logger.info('Max. PDF uncertainty: %.3f', max(cleanPDF))
    #logger.info('Min. PDF uncertainty: %.3f', min(PDF_unc))
    #logger.info('Max. PDF uncertainty: %.3f', max(PDF_unc))
    
    logger.info('Min. scale uncertainty: %.3f', min(Scale_unc))
    logger.info('Max. scale uncertainty: %.3f', max(Scale_unc))
    
    #logger.info('Min. PS scale uncertainty: %.3f', min(PS_unc))
    #logger.info('Max. PS scale uncertainty: %.3f', max(PS_unc))




