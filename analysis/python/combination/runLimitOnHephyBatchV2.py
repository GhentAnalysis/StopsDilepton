#!/usr/bin/env python

#!/usr/bin/env python
from optparse import OptionParser
parser = OptionParser()
parser.add_option("--year",                 dest='year',  action='store', default='2016',    choices=["2016", "2017", "2018"], help="which year?")
parser.add_option("--submit",               action='store_true', help="Really submit?")
parser.add_option("--combine",              action='store_true', help="Combine years?")
parser.add_option("--significance",         action='store_true', help="Calculate significance?")
parser.add_option("--useTxt",               action='store_true', help="Calculate significance?")
parser.add_option("--dryRun",               action='store_true', help="Calculate significance?")
parser.add_option("--subset",               action='store_true', help="Don't run all the points")
parser.add_option("--overwrite",            action='store_true', help="Overwrite")
parser.add_option("--signal",               action='store', default=None, help="Which signal, if any?")
parser.add_option("--jobTitle",             action='store', default=None, help="How to call the job?")
(options, args) = parser.parse_args()

import os

year = int(options.year)

if options.signal:
    if options.signal == 'T2tt':
        if year == 2016:
            data_directory              = '/afs/hephy.at/data/cms06/nanoTuples/'
            postProcessing_directory    = 'stops_2016_nano_v0p23/dilep/'
            from StopsDilepton.samples.nanoTuples_FastSim_Summer16_postProcessed import signals_T2tt as jobs
        elif year == 2017:
            data_directory              = '/afs/hephy.at/data/cms06/nanoTuples/'
            postProcessing_directory    = 'stops_2017_nano_v0p23/dilep/'
            from StopsDilepton.samples.nanoTuples_FastSim_Fall17_postProcessed import signals_T2tt as jobs
        if year == 2018:
            data_directory              = '/afs/hephy.at/data/cms06/nanoTuples/'
            postProcessing_directory    = 'stops_2018_nano_v0p23/dilep/'
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

    elif options.signal == 'TTbarDM':
        if year == 2016:
            data_directory              = '/afs/hephy.at/data/cms03/nanoTuples/'
            postProcessing_directory    = 'stops_2016_nano_v0p25/dilep/'
            from StopsDilepton.samples.nanoTuples_Summer16_TTDM_postProcessed import signals as jobs
        elif year == 2017:
            data_directory              = '/afs/hephy.at/data/cms03/nanoTuples/'
            postProcessing_directory    = 'stops_2017_nano_v0p25/dilep/'
            from StopsDilepton.samples.nanoTuples_Fall17_TTDM_postProcessed import signals as jobs
        elif year == 2018:
            data_directory              = '/afs/hephy.at/data/cms03/nanoTuples/'
            postProcessing_directory    = 'stops_2018_nano_v0p25/dilep/'
            from StopsDilepton.samples.nanoTuples_Autumn18_TTDM_postProcessed import signals as jobs

signalEstimators = [ s.name for s in jobs ]


import time

jobTitle = 'li'+str(year-2000)+options.signal if not options.jobTitle else options.jobTitle

#print jobTitle

if options.submit:
    cmd = "submitBatch.py --title='%s'"%jobTitle
else:
    cmd = "echo"

#print len(signalEstimators)
for i, estimator in enumerate(signalEstimators):
    mass = estimator.split('_')[1:]
    os.system(cmd+" 'python run_limit_combined.py --signal %s --only=%s  --version v9 --scaleTTZ --overwrite'"%(options.signal, estimator))
    #    else:
    #        os.system(cmd+" 'python run_limit.py --signal %s --year %s --only=%s --version v9 --skipFitDiagnostics --unblind --fitAll'"%(options.signal, str(year), estimator))

