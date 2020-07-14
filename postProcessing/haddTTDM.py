import ROOT
import os
import subprocess

def get_parser():
    ''' Argument parser for post-processing module.
    '''
    import argparse
    argParser = argparse.ArgumentParser(description = "Argument parser for cmgPostProcessing")

    argParser.add_argument('--logLevel',    action='store',         nargs='?',  choices=['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'TRACE', 'NOTSET'],   default='INFO', help="Log level for logging" )
    argParser.add_argument('--samples',     action='store',         nargs='*',  type=str, default=['TTZToLLNuNu_ext'],                  help="List of samples to be post-processed, given as CMG component name" )
    argParser.add_argument('--year',        action='store',                     type=int,                                               help="Which year?" )
#    argParser.add_argument('--nJobs',       action='store',         nargs='?',  type=int, default=1,                                    help="Maximum number of simultaneous jobs." )
#    argParser.add_argument('--job',         action='store',                     type=int, default=0,                                    help="Run only jobs i" )
    return argParser

options = get_parser().parse_args()

# Logging
import StopsDilepton.tools.logger as _logger
logger  = _logger.get_logger(options.logLevel, logFile = None)

import Analysis.Tools.logger as _logger_an
logger_an = _logger_an.get_logger('DEBUG', logFile = None )

import RootTools.core.logger as _logger_rt
logger_rt = _logger_rt.get_logger(options.logLevel, logFile = None )

if options.year == 2016:
    from Samples.nanoAOD.Summer16_private_legacy_v1 import allSamples as mcSamples
    allSamples = mcSamples
elif options.year == 2017:
    from Samples.nanoAOD.Fall17_private_legacy_v1   import allSamples as mcSamples
    allSamples = mcSamples
elif options.year == 2018:
    from Samples.nanoAOD.Autumn18_private_legacy_v1 import allSamples as mcSamples
    allSamples = mcSamples
else:
    raise NotImplementedError

samples = []
for selectedSamples in options.samples:
    for sample in allSamples:
        if selectedSamples == sample.name:
            samples.append(sample)


if len(samples)==0:
    logger.info( "No samples found. Was looking for %s. Exiting" % options.samples )
    sys.exit(-1)

#Samples: combine if more than one
if len(samples)>1:
    sample_name =  samples[0].name+"_comb"
    logger.info( "Combining samples %s to %s.", ",".join(s.name for s in samples), sample_name )
    sample      = Sample.combine(sample_name, samples, maxN = maxN)
elif len(samples)==1:
    sample      = samples[0]
else:
    raise ValueError( "Need at least one sample. Got %r",samples )

sampleName = sample.name

len_orig = len(sample.files)
## sort the list of files?

outDir = '/afs/hephy.at/data/cms05/nanoAOD/TTDM/%s/%s/'%(options.year, sampleName)
if not os.path.isdir(outDir):
    os.makedirs(outDir)

nJobs = len(sample.files)/20 # merge max 20 nanoAOD files

for j in  range(nJobs):

    splitSample = sample.split( n=nJobs, nSub=j)
    files = ' '.join(splitSample.files)
    # use the local redirector, some files are not found when using the global ones
    files = files.replace('root://cms-xrd-global.cern.ch///store','root://hephyse.oeaw.ac.at//dpm/oeaw.ac.at/home/cms/store')
    # use haddnano.py and not the processor to speed up things a bit
    command = "python ../../PhysicsTools/NanoAODTools/scripts/haddnano.py %s/nanoAOD_%s.root %s"%(outDir, j, files)
    logger.info("Working on job %s of %s", j, nJobs)
    subprocess.call(command, shell=True)


