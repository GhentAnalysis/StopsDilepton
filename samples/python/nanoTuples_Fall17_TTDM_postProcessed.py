import copy, os, sys
from RootTools.core.Sample import Sample
import ROOT

# Logging
import logging
logger = logging.getLogger(__name__)

from StopsDilepton.samples.color import color

# Data directory
try:
    data_directory_ = sys.modules['__main__'].data_directory
except:
    from StopsDilepton.samples.default_locations import default_locations
    data_directory_ = default_locations.mc_2017_data_directory 

# Take post processing directory if defined in main module
try:
  import sys
  postProcessing_directory_ = sys.modules['__main__'].postProcessing_directory
except:
  from StopsDilepton.samples.default_locations import default_locations
  postProcessing_directory_ = default_locations.mc_2017_postProcessing_directory 

logger.info("Loading MC samples from directory %s", os.path.join(data_directory_, postProcessing_directory_))

directory = os.path.join( data_directory_, postProcessing_directory_)

from Samples.nanoAOD.Fall17_private_legacy_v1 import TTbarDMJets_Dilepton_pseudoscalar, TTbarDMJets_Dilepton_scalar
from helpers import getTTDMSignalWeight

pseudo = getTTDMSignalWeight(TTbarDMJets_Dilepton_pseudoscalar, 1, year=2017)
scalar = getTTDMSignalWeight(TTbarDMJets_Dilepton_scalar, 1, year=2017)

signals = []


for key in pseudo.keys():
    #pseudo[s]
    name = "TTDM_%s_%s_%s"%(pseudo[key]['spin'], pseudo[key]['mPhi'], pseudo[key]['mChi'])
    tmp = Sample.fromFiles(\
        name = name,
        files = [os.path.join(data_directory_, postProcessing_directory_,'TTDM',name+'.root')],
        treeName = "Events",
        isData = False,
        color = 8 ,
        texName = "PS (%s,%s)"%(pseudo[key]['mPhi'],pseudo[key]['mChi'])
    )

    tmp.mChi = int(pseudo[key]['mChi'])
    tmp.mPhi = int(pseudo[key]['mPhi'])
    tmp.type = 'PS'
    tmp.isFastSim = False

    exec("%s=tmp"%name)
    exec("signals.append(%s)"%name)


for key in scalar.keys():
    #pseudo[s]
    name = "TTDM_%s_%s_%s"%(scalar[key]['spin'], scalar[key]['mPhi'], scalar[key]['mChi'])
    tmp = Sample.fromFiles(\
        name = name,
        files = [os.path.join(data_directory_, postProcessing_directory_,'TTDM',name+'.root')],
        treeName = "Events",
        isData = False,
        color = 8 ,
        texName = "S (%s,%s)"%(scalar[key]['mPhi'],scalar[key]['mChi'])
    )

    tmp.mChi = int(scalar[key]['mChi'])
    tmp.mPhi = int(scalar[key]['mPhi'])
    tmp.type = 'S'
    tmp.isFastSim = False

    exec("%s=tmp"%name)
    exec("signals.append(%s)"%name)

logger.info("Loaded %i TTDM signals", len(signals))
logger.debug("Loaded TTDM signals: %s", ",".join([s.name for s in signals]))



## Pseudoscalar
#TTDM_pseudoscalar_100_1     = Sample.fromFiles(name="TTDM_pseudoscalar_100_1",  treeName="Events", isData=False, texName="PS(100,1)", files=[directory+"TTDM_pseudoscalar_100_1.root"])
#TTDM_pseudoscalar_100_20     = Sample.fromFiles(name="TTDM_pseudoscalar_100_20",  treeName="Events", isData=False, texName="PS(100,1)", files=[directory+"TTDM_pseudoscalar_100_20.root"])
#TTDM_pseudoscalar_100_30     = Sample.fromFiles(name="TTDM_pseudoscalar_100_30",  treeName="Events", isData=False, texName="PS(100,1)", files=[directory+"TTDM_pseudoscalar_100_30.root"])
#TTDM_pseudoscalar_100_40     = Sample.fromFiles(name="TTDM_pseudoscalar_100_40",  treeName="Events", isData=False, texName="PS(100,1)", files=[directory+"TTDM_pseudoscalar_100_40.root"])
#TTDM_pseudoscalar_100_45     = Sample.fromFiles(name="TTDM_pseudoscalar_100_45",  treeName="Events", isData=False, texName="PS(100,1)", files=[directory+"TTDM_pseudoscalar_100_45.root"])
#TTDM_pseudoscalar_100_49     = Sample.fromFiles(name="TTDM_pseudoscalar_100_49",  treeName="Events", isData=False, texName="PS(100,1)", files=[directory+"TTDM_pseudoscalar_100_49.root"])
#TTDM_pseudoscalar_100_51     = Sample.fromFiles(name="TTDM_pseudoscalar_100_51",  treeName="Events", isData=False, texName="PS(100,1)", files=[directory+"TTDM_pseudoscalar_100_51.root"])
#TTDM_pseudoscalar_100_55     = Sample.fromFiles(name="TTDM_pseudoscalar_100_55",  treeName="Events", isData=False, texName="PS(100,1)", files=[directory+"TTDM_pseudoscalar_100_55.root"])
#TTDM_pseudoscalar_150_1     = Sample.fromFiles(name="TTDM_pseudoscalar_150_1",  treeName="Events", isData=False, texName="PS(100,1)", files=[directory+"TTDM_pseudoscalar_150_1.root"])
#TTDM_pseudoscalar_200_1     = Sample.fromFiles(name="TTDM_pseudoscalar_200_1",  treeName="Events", isData=False, texName="PS(100,1)", files=[directory+"TTDM_pseudoscalar_200_1.root"])
#TTDM_pseudoscalar_250_1     = Sample.fromFiles(name="TTDM_pseudoscalar_250_1",  treeName="Events", isData=False, texName="PS(100,1)", files=[directory+"TTDM_pseudoscalar_250_1.root"])
#TTDM_pseudoscalar_300_1     = Sample.fromFiles(name="TTDM_pseudoscalar_300_1",  treeName="Events", isData=False, texName="PS(100,1)", files=[directory+"TTDM_pseudoscalar_300_1.root"])
#TTDM_pseudoscalar_350_1     = Sample.fromFiles(name="TTDM_pseudoscalar_350_1",  treeName="Events", isData=False, texName="PS(100,1)", files=[directory+"TTDM_pseudoscalar_350_1.root"])
#TTDM_pseudoscalar_400_1     = Sample.fromFiles(name="TTDM_pseudoscalar_400_1",  treeName="Events", isData=False, texName="PS(100,1)", files=[directory+"TTDM_pseudoscalar_400_1.root"])
#TTDM_pseudoscalar_450_1     = Sample.fromFiles(name="TTDM_pseudoscalar_450_1",  treeName="Events", isData=False, texName="PS(100,1)", files=[directory+"TTDM_pseudoscalar_450_1.root"])
#TTDM_pseudoscalar_500_1     = Sample.fromFiles(name="TTDM_pseudoscalar_500_1",  treeName="Events", isData=False, texName="PS(100,1)", files=[directory+"TTDM_pseudoscalar_500_1.root"])
#TTDM_pseudoscalar_50_1     = Sample.fromFiles(name="TTDM_pseudoscalar_50_1",  treeName="Events", isData=False, texName="PS(100,1)", files=[directory+"TTDM_pseudoscalar_50_1.root"])



