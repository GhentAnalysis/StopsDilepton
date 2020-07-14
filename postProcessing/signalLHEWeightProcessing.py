#!/usr/bin/env python
''' Make flat ntuple with LHE weights 
'''
#
# Standard imports and batch mode
#
import ROOT
import os, sys
ROOT.gROOT.SetBatch(True)
from math                                import sqrt, cos, sin, pi, acos
import imp

# RootTools
from RootTools.core.standard             import *

# Analysis
from Analysis.Tools.helpers              import checkRootFile, deepCheckRootFile

# User specific
import StopsDilepton.tools.user as user

# FWLite
from DataFormats.FWLite import Handle,Lumis

#
# Arguments
# 
import argparse
argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--logLevel',  action='store',      default='INFO',          nargs='?', choices=['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'TRACE', 'NOTSET'], help="Log level for logging")
argParser.add_argument('--small',     action='store_true', help='Run only on a small subset of the data?')#, default = True)
argParser.add_argument('--overwrite', action='store_true', help='Overwrite?')#, default = True)
argParser.add_argument('--TTDM',      action='store_true', help='Overwrite?')#, default = True)
argParser.add_argument('--targetDir', action='store',      nargs='?',  type=str, default=user.postprocessing_output_directory, help="Name of the directory the post-processed files will be saved" )
argParser.add_argument('--sample',    action='store',      default='SMS_T2tt_mStop_150to250', help="Name of the sample loaded from fwlite_benchmarks. Only if no inputFiles are specified")
argParser.add_argument('--year',      action='store',      type=int,                        help="Which year?" )
argParser.add_argument('--nJobs',     action='store',      nargs='?', type=int, default=1,  help="Maximum number of simultaneous jobs.")
argParser.add_argument('--job',       action='store',      nargs='?', type=int, default=0,  help="Run only job i")
argParser.add_argument('--fullsim',   action='store_true', help='You want to load fullsim samples?')#, default = True)
args = argParser.parse_args()

#
# Logger
#
import StopsDilepton.tools.logger as _logger
import RootTools.core.logger as _logger_rt
logger    = _logger.get_logger(   args.logLevel, logFile = None)
logger_rt = _logger_rt.get_logger(args.logLevel, logFile = None)

if args.fullsim:
    if   args.year == 2016:
        from Samples.miniAOD.Summer16_miniAODv3 import samples
    elif args.year == 2017:
        from Samples.miniAOD.Fall17_miniAODv2 import samples
    elif args.year == 2018:
        from Samples.miniAOD.Autumn18_miniAODv1 import samples
else:
    if   args.year == 2016:
        from Samples.miniAOD.Summer16_Fast_miniAODv3 import samples
    elif args.year == 2017:
        from Samples.miniAOD.Fall17_Fast_miniAODv2 import samples
    elif args.year == 2018:
        from Samples.miniAOD.Autumn18_Fast_miniAODv1 import samples

sample = getattr(samples, args.sample)

sub_directory = '%s'%args.year
maxEvents = -1
if args.small: 
    sub_directory += "_small"
    #maxEvents=100 
    sample.files=sample.files[:1]

# output directory
directory = os.path.join( args.targetDir, 'signalWeights', sub_directory ) 
output_directory = os.path.join( directory, sample.name )

if not os.path.exists( output_directory ): 
    os.makedirs( output_directory )
    logger.info( "Created output directory %s", output_directory )

variables = ["run/I", "luminosityBlock/I", "event/l"]
variables += ["LHE_weight_original/F", "LHE[weight/F]"]

# read weight names
lumis           = Lumis(sample.files)
genLumiHandle   = Handle('GenLumiInfoHeader')
LHEWeightsNames = []
for lumi in lumis:
    if lumi.getByLabel('generator',genLumiHandle):
        weightNames = genLumiHandle.product().weightNames()
        for wn in weightNames:  #direct cast is not working properly, copy of elements is needed
            LHEWeightsNames.append(wn)
        break

# Run only job number "args.job" from total of "args.nJobs"
if args.nJobs>1:
    n_files_before = len(sample.files)
    sample = sample.split(args.nJobs)[args.job]
    n_files_after  = len(sample.files)
    logger.info( "Running job %i/%i over %i files from a total of %i.", args.job, args.nJobs, n_files_after, n_files_before)

# products to read
products = { 'generator':{'type':'GenEventInfoProduct', 'label':('generator') } }

def fill_vector_collection( event, collection_name, collection_varnames, objects):
    setattr( event, "n"+collection_name, len(objects) )
    for i_obj, obj in enumerate(objects):
        for var in collection_varnames:
            getattr(event, collection_name+"_"+var)[i_obj] = obj[var]
def fill_vector( event, collection_name, collection_varnames, obj):
    for var in collection_varnames:
        try:
            setattr(event, collection_name+"_"+var, obj[var] )
        except TypeError as e:
            logger.error( "collection_name %s var %s obj[var] %r", collection_name, var,  obj[var] )
            raise e
        except KeyError as e:
            logger.error( "collection_name %s var %s obj[var] %r", collection_name, var,  obj[var] )
            raise e

logger.info( "Running over files: %s", ", ".join(sample.files ) )


TTDM_scale_indices  = [6,11,16,21,26,31,36,41] 
TTDM_PDF_indices    = range(46,147)
TTDM_aS_indices     = [147,148]


def filler( event ):
    event.run, event.luminosityBlock, event.event = reader.evt
    LHE_weights                 = reader.products['generator'].weights()
    default_indices = LHE_weights[1:10] if not args.TTDM else TTDM_scale_indices+TTDM_PDF_indices+TTDM_aS_indices
    event.LHE_weight_original   = LHE_weights[0]
    fill_vector_collection( event, "LHE", ["weight"], [ {'weight':w} for w in default_indices ] )

tmp_dir         = ROOT.gDirectory
output_filename = os.path.join(output_directory, sample.name + '.root')

_logger.   add_fileHandler( output_filename.replace('.root', '.log'), args.logLevel )
_logger_rt.add_fileHandler( output_filename.replace('.root', '_rt.log'), args.logLevel )

if os.path.exists( output_filename ) and checkRootFile( output_filename, checkForObjects=["Events"] ) and deepCheckRootFile( output_filename ) and not args.overwrite:
    logger.info( "File %s found. Quit.", output_filename )
    sys.exit(0)

# FWLite reader if this is an EDM file
reader = sample.fwliteReader( products = products )

output_file = ROOT.TFile( output_filename, 'recreate')
output_file.cd()
maker = TreeMaker(
    sequence  = [ filler ],
    variables = [ TreeVariable.fromString(x) for x in variables ],
    treeName = "Events"
    )

tmp_dir.cd()

counter = 0
reader.start()
maker.start()

while reader.run():
    maker.run()
    counter += 1
    if counter == maxEvents:  break

logger.info( "Done with running over %i events.", reader.nEvents )

output_file.cd()
maker.tree.Write()
output_file.Close()

logger.info( "Written output file %s", output_filename )
