# manually hadded nanoAODs
from RootTools.core.standard import *

from Samples.nanoAOD.Fall17_private_legacy_v1 import TTbarDMJets_Dilepton_scalar as TTbarDMJets_Dilepton_scalar_nanoAOD
from Samples.nanoAOD.Fall17_private_legacy_v1 import TTbarDMJets_Dilepton_pseudoscalar as TTbarDMJets_Dilepton_pseudoscalar_nanoAOD

TTDM_Dilepton_scalar = Sample.fromDirectory('TTDM_Dilepton_scalar', '/afs/hephy.at/data/cms05/nanoAOD/TTDM/2017/TTbarDMJets_Dilepton_scalar/')
TTDM_Dilepton_scalar.normalization, TTDM_Dilepton_scalar.xSection = TTbarDMJets_Dilepton_scalar_nanoAOD.normalization, TTbarDMJets_Dilepton_scalar_nanoAOD.xSection

TTDM_Dilepton_pseudoscalar = Sample.fromDirectory('TTDM_Dilepton_pseudoscalar', '/afs/hephy.at/data/cms05/nanoAOD/TTDM/2017/TTbarDMJets_Dilepton_pseudoscalar/')
TTDM_Dilepton_pseudoscalar.normalization, TTDM_Dilepton_pseudoscalar.xSection = TTbarDMJets_Dilepton_pseudoscalar_nanoAOD.normalization, TTbarDMJets_Dilepton_pseudoscalar_nanoAOD.xSection

allSamples = [TTDM_Dilepton_scalar, TTDM_Dilepton_pseudoscalar]
