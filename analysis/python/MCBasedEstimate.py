from StopsDilepton.tools.helpers import getYieldFromChain
from setupHelpers import loadChain
from u_float import u_float 
from systematics import SystematicBaseClass
from StopsDilepton.tools.helpers import printHeader

class MCBasedEstimate(SystematicBaseClass):
  def __init__(self, name, sample, cacheDir=None):
    super(MCBasedEstimate, self).__init__(name, cacheDir=cacheDir)
    self.sample=sample
#Concrete implementation of abstract method 'estimate' as defined in Systematic
  def _estimate(self, region, channel, setup):
    printHeader("MC prediction for %s channel %s" %(self.name, channel))

    if channel=='all':
      return sum( [ self.cachedEstimate(region, c, setup) for c in ['MuMu', 'EE', 'EMu'] ], u_float(0., 0.) )
    else:
      preSelection = setup.preselection('MC', channel=channel)
      cut = "&&".join([region.cutString(setup.sys['selectionModifier']), preSelection['cut']])
      weight = preSelection['weightStr']

      if setup.verbose: 
        print "Using cut %s and weight %s"%(cut, weight)
      if not self.sample[channel].has_key('chain'):
        loadChain(self.sample[channel])
      return setup.lumi[channel]/1000.*u_float(getYieldFromChain(self.sample[channel]['chain'], cutString = cut, weight=weight, returnError = True) )
