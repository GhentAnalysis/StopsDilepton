jmeVariations = ["JER", "JERUp", "JERDown", "JECUp", "JECDown"]

btagVariationWeights = [
'reweightBTag_central',
'reweightBTag_up_hfstats1',
'reweightBTag_down_hfstats1',
'reweightBTag_up_hfstats2',
'reweightBTag_down_hfstats2',
'reweightBTag_up_hf',
'reweightBTag_down_hf',
'reweightBTag_up_lfstats1',
'reweightBTag_down_lfstats1',
'reweightBTag_up_lfstats2',
'reweightBTag_down_lfstats2',
'reweightBTag_up_lf',
'reweightBTag_down_lf',
'reweightBTag_up_cferr1',
'reweightBTag_down_cferr1',
'reweightBTag_up_cferr2',
'reweightBTag_down_cferr2',
#'reweightBTag_up_jes',
#'reweightBTag_down_jes',
]

#Abstract base class for systematics
import os
import abc
from math import sqrt
from cache import cache
import json
class SystematicBaseClass:
  __metaclass__ = abc.ABCMeta

  def __init__(self, name, cacheDir=None):
    self.name = name
    if cacheDir:
      self.cacheDir=cacheDir
      cacheFileName = os.path.join(cacheDir, self.name+'.pkl')
      if not os.path.exists(os.path.dirname(cacheFileName)):
        os.makedirs(os.path.dirname(cacheFileName))
      self.cache = cache(cacheFileName, verbosity=2)
    else:
      self.cache=None

  def uniqueKey(self, region, channel, setup):
    return region, channel, json.dumps(setup.lumi, sort_keys=True), json.dumps(setup.sys, sort_keys=True)

  def cachedEstimate(self, region, channel, setup, save=True):
    key =  self.uniqueKey(region, channel, setup)
    if self.cache and self.cache.contains(key):
      res = self.cache.get(key)
      if setup.verbose: print "Loading cached %s result for %r : %r"%(self.name, key, res)
      return res
    elif self.cache:
      return self.cache.add( key, self._estimate( region, channel, setup), save=save)
    else:
      return self._estimate( region, channel, setup)

  @abc.abstractmethod
  def _estimate(self, region, channel, setup):
    '''Estimate yield in 'region' using setup'''
    return

  def PUSystematic(self, region, channel, setup):
    ref  = self.cachedEstimate(region, channel, setup)
    up   = self.cachedEstimate(region, channel, setup.sysClone({'weight':'weightPUUp'}))
    down = self.cachedEstimate(region, channel, setup.sysClone({'weight':'weightPUDown'}))
    return 0.5*(up-down)/ref

  def topPtSystematic(self, region, channel, setup):
    ref   = self.cachedEstimate(region, channel, setup)
    up   = self.cachedEstimate(region, channel, setup.sysClone({'reweight':['reweightTopPt']}))
    return 0.5*(up-ref)/ref
 
  def JERSystematic(self, region, channel, setup):
    ref   = self.cachedEstimate(region, channel, setup)
    up   = self.cachedEstimate(region, channel, setup.sysClone({'selectionModifier':'JERUp'}))
    down = self.cachedEstimate(region, channel, setup.sysClone({'selectionModifier':'JERDown'}))
    return 0.5*(up-down)/ref

  def JECSystematic(self, region, channel, setup):
    ref   = self.cachedEstimate(region, channel, setup)
    up   = self.cachedEstimate(region, channel, setup.sysClone({'selectionModifier':'JECUp'}))
    down = self.cachedEstimate(region, channel, setup.sysClone({'selectionModifier':'JECDown'}))
    return 0.5*(up-down)/ref 

  def leptonFSSystematic(self, region, channel, setup):
    ref  = self.cachedEstimate(region, channel, setup.sysClone({'reweight':['reweightLeptonFastSimSF']}))
    up   = self.cachedEstimate(region, channel, setup.sysClone({'reweight':['reweightLeptonFastSimSFUp']}))
    down = self.cachedEstimate(region, channel, setup.sysClone({'reweight':['reweightLeptonFastSimSFDown']}))
    return 0.5*(up-down)/ref 

  def btaggingSFbSystematic(self, region, channel, setup):
    ref     = self.cachedEstimate(region, channel, setup.sysClone({'useBTagWeights':'SF'}))     
    up      = self.cachedEstimate(region, channel, setup.sysClone({'useBTagWeights':'SF_b_Up'}))
    down    = self.cachedEstimate(region, channel, setup.sysClone({'useBTagWeights':'SF_b_Down'}))
    return 0.5*(up-down)/ref 

  def btaggingSFlSystematic(self, region, channel, setup):
    ref     = self.cachedEstimate(region, channel, setup.sysClone({'useBTagWeights':'SF'}))
    up      = self.cachedEstimate(region, channel, setup.sysClone({'useBTagWeights':'SF_l_Up'}))
    down    = self.cachedEstimate(region, channel, setup.sysClone({'useBTagWeights':'SF_l_Down'}))
    return 0.5*(up-down)/ref 

  def btaggingSFFSSystematic(self, region, channel, setup):
    ref     = self.cachedEstimate(region, channel, setup.sysClone({'useBTagWeights':'SF'}))
    up      = self.cachedEstimate(region, channel, setup.sysClone({'useBTagWeights':'SF_FS_Up'}))
    down    = self.cachedEstimate(region, channel, setup.sysClone({'useBTagWeights':'SF_FS_Down'}))
    return 0.5*(up-down)/ref 

  def getAllSysJobs(self, region, channel, setup):
    return [
      (region, channel, setup.sysClone({'weight':'weightPUUp'})),
      (region, channel, setup.sysClone({'weight':'weightPUDown'})),

      (region, channel, setup.sysClone({'reweight':['reweightTopPt']})),
   
      (region, channel, setup.sysClone({'selectionModifier':'JERUp'})),
      (region, channel, setup.sysClone({'selectionModifier':'JERDown'})),

      (region, channel, setup.sysClone({'selectionModifier':'JECUp'})),
      (region, channel, setup.sysClone({'selectionModifier':'JECDown'})),

#      (region, channel, setup.sysClone({'reweight':['reweightLeptonFastSimSF']})),
#      (region, channel, setup.sysClone({'reweight':['reweightLeptonFastSimSFUp']})),
#      (region, channel, setup.sysClone({'reweight':['reweightLeptonFastSimSFDown']})),

      (region, channel, setup.sysClone({'useBTagWeights':'SF'})),
      (region, channel, setup.sysClone({'useBTagWeights':'SF_b_Up'})),
      (region, channel, setup.sysClone({'useBTagWeights':'SF_b_Down'})),

      (region, channel, setup.sysClone({'useBTagWeights':'SF'})),
      (region, channel, setup.sysClone({'useBTagWeights':'SF_l_Up'})),
      (region, channel, setup.sysClone({'useBTagWeights':'SF_l_Down'})),

      (region, channel, setup.sysClone({'useBTagWeights':'SF'})),
#      (region, channel, setup.sysClone({'useBTagWeights':'SF_FS_Up'})),
#      (region, channel, setup.sysClone({'useBTagWeights':'SF_FS_Down'})),
    ]
