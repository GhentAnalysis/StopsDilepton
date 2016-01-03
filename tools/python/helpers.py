import ROOT
from math import pi, sqrt, cos, sin, sinh, log, cosh
from array import array
ROOT.gROOT.LoadMacro("$CMSSW_BASE/src/StopsDilepton/tools/scripts/tdrstyle.C")
ROOT.setTDRStyle()


def deltaPhi(phi1, phi2):
  dphi = phi2-phi1
  if  dphi > pi:
    dphi -= 2.0*pi
  if dphi <= -pi:
    dphi += 2.0*pi
  return abs(dphi)

def deltaR2(l1, l2):
  return deltaPhi(l1['phi'], l2['phi'])**2 + (l1['eta'] - l2['eta'])**2

def deltaR(l1, l2):
  return sqrt(deltaR2(l1,l2))

def getFileList(dir, histname='histo', maxN=-1):
  import os
  filelist = os.listdir(os.path.expanduser(dir))
  filelist = [dir+'/'+f for f in filelist if histname in f and f.endswith(".root")]
  if maxN>=0:
    filelist = filelist[:maxN]
  return filelist

def getChain(sampleList, histname='', maxN=-1, treeName="Events"):
  if not type(sampleList)==type([]):
    sampleList_ = [sampleList]
  else:
    sampleList_= sampleList 
  c = ROOT.TChain(treeName)
  i=0
  for s in sampleList_:
    if type(s)==type(""):
      for f in getFileList(s, histname, maxN):
        if checkRootFile(f, checkForObjects=[treeName]): 
          i+=1
          c.Add(f)
        else:
          print "File %s looks broken."%f
      print "Added ",i,'files from samples %s' %(", ".join([s['name'] for s in sampleList_]))
    elif type(s)==type({}):
      if s.has_key('file'):
        c.Add(s['file'])
#        print "Added file %s"%s['file']
        i+=1
      if s.has_key('bins'):
        for b in s['bins']:
          dir = s['dirname'] if s.has_key('dirname') else s['dir']
          for f in getFileList(dir+'/'+b, histname, maxN):
            if checkRootFile(f, checkForObjects=[treeName]): 
              i+=1
              c.Add(f)
            else:
              print "File %s looks broken."%f
#      print 'Added %i files from %i elements' %(i, len(sampleList))
    else: 
#      print sampleList
      print "Could not load chain from sampleList %s"%repr(sampleList)
  return c

#def getChunks(sample,  maxN=-1):
##  print "sample" , sample , maxN
#  import os, subprocess
#  #print "sample dir:" , sample['dir']
#  chunks = [{'name':x} for x in os.listdir(sample['dir']) if x.startswith(sample['chunkString']+'_Chunk') or x==sample['name']]
#  #print chunks
#  chunks=chunks[:maxN] if maxN>0 else chunks
#  sumWeights=0
#  failedChunks=[]
#  goodChunks  =[] 
#  const = 'All Events' if sample['isData'] else 'Sum Weights'
#  for i, s in enumerate(chunks):
#      if not sample.has_key("skimAnalyzerDir"):
#        logfile = sample['dir']+'/'+s['name']+'/SkimReport.txt'
#      else:
#        logfile = sample['dir']+'/'+s['name']+"/"+sample["skimAnalyzerDir"]+'/SkimReport.txt'
#      if os.path.isfile(logfile):
#        line = [x for x in subprocess.check_output(["cat", logfile]).split('\n') if x.count(const)]
#        assert len(line)==1,"Didn't find normalization constant '%s' in  number in file %s"%(const, logfile)
#        #n = int(float(line[0].split()[2]))
#        sumW = float(line[0].split()[2])
#        inputFilename = sample['dir']+'/'+s['name']+'/'+sample['rootFileLocation']
#        #print sumW, inputFilename
#        if os.path.isfile(inputFilename):
#          sumWeights+=sumW
#          s['file']=inputFilename
#          goodChunks.append(s)
#        else:
#          failedChunks.append(chunks[i])
#      else:
#        print "log file not found:  ", logfile
#        failedChunks.append(chunks[i])
##    except: print "Chunk",s,"could not be added"
#  eff = round(100*len(failedChunks)/float(len(chunks)),3)
#  print "Chunks: %i total, %i good (normalization constant %f), %i bad. Inefficiency: %f"%(len(chunks),len(goodChunks),sumWeights,len(failedChunks), eff)
#  for s in failedChunks: 
#    print "Failed:",s
#  return goodChunks, sumWeights

def checkRootFile(f, checkForObjects=[]):
  rf = ROOT.TFile.Open(f)
  try: 
    good = (not rf.IsZombie()) and (not rf.TestBit(ROOT.TFile.kRecovered))
  except:
    return False
  for o in checkForObjects:
    if not rf.GetListOfKeys().Contains(o):
      print "[checkRootFile] Failed to find object %s in file %s"%(o, f) 
      rf.Close()
      return False
#    print "Keys recoveredd %i zombie %i tb %i"%(rf.Recover(), rf.IsZombie(), rf.TestBit(ROOT.TFile.kRecovered))
  rf.Close()
  return good

def getChunks(sample,  maxN=-1):
  import os, subprocess
  chunks = [{'name':x} for x in os.listdir(sample.path) if x.startswith(sample.chunkString+'_Chunk') or x==sample.chunkString]
  chunks=chunks[:maxN] if maxN>0 else chunks
  sumWeights=0
  failedChunks=[]
  goodChunks  =[] 
  const = 'All Events' if sample.isData else 'Sum Weights'
  for i, s in enumerate(chunks):
      logfile = "/".join([sample.path, s['name'], sample.skimAnalyzerDir,'SkimReport.txt'])
#      print logfile 
      if os.path.isfile(logfile):
        line = [x for x in subprocess.check_output(["cat", logfile]).split('\n') if x.count(const)]
        assert len(line)==1,"Didn't find normalization constant '%s' in  number in file %s"%(const, logfile)
        n = int(float(line[0].split()[2]))
        sumW = float(line[0].split()[2])
        inputFilename = '/'.join([sample.path, s['name'], sample.rootFileLocation])
#        print sumW, inputFilename
        if os.path.isfile(inputFilename) and checkRootFile(inputFilename):
          sumWeights+=sumW
          s['file']=inputFilename
          goodChunks.append(s)
        else:
          failedChunks.append(chunks[i])
      else:
        print "log file not found:  ", logfile
        failedChunks.append(chunks[i])
#    except: print "Chunk",s,"could not be added"
  eff = round(100*len(failedChunks)/float(len(chunks)),3)
  print "Chunks: %i total, %i good (normalization constant %f), %i bad. Inefficiency: %f"%(len(chunks),len(goodChunks),sumWeights,len(failedChunks), eff)
  for s in failedChunks: 
    print "Failed:",s
  return goodChunks, sumWeights

def getObjFromFile(fname, hname):
  f = ROOT.TFile(fname)
  assert not f.IsZombie()
  f.cd()
  htmp = f.Get(hname)
  if not htmp:  return htmp
  ROOT.gDirectory.cd('PyROOT:/')
  res = htmp.Clone()
  f.Close()
  return res

def writeObjToFile(fname, obj):
  gDir = ROOT.gDirectory.GetName()
  f = ROOT.TFile(fname, 'recreate')
  objw = obj.Clone()
  objw.Write()
  f.Close()
  ROOT.gDirectory.cd(gDir+':/')
  return

def getVarValue(c, var, n=-1):
  att = getattr(c, var)
  if n>=0:
#    print "getVarValue %s %i"%(var,n)
    if n<att.__len__():
      return att[n]
    else:
      return float('nan')
  return att

  
#def getVarValue(c, var, n=-1):
#  try:
#    att = getattr(c, var)
#    if n>=0:return att[n]
#    return att
#  except:  
#    l = c.GetLeaf(var)
#    if l:return l.GetValue(n)
#    return float('nan')

def getEList(chain, cut, newname='eListTMP'):
  chain.Draw('>>eListTMP_t', cut)
  #elistTMP_t = ROOT.gROOT.Get('eListTMP_t')
  elistTMP_t = ROOT.gDirectory.Get('eListTMP_t')
  elistTMP = elistTMP_t.Clone(newname)
  del elistTMP_t
  return elistTMP

def getObjDict(c, prefix, variables, i):
  res={var: getVarValue(c, prefix+var, i) for var in variables}
  res['index']=i
  return res
#  return {var: c.GetLeaf(prefix+var).GetValue(i) for var in variables}

## There are several problems with the function below and it seems not to be used
#def getWeight(c,sample,lumi,n=0):
#  genweight_value    = c.GetLeaf("genWeight").GetValue(n)
#  lumi_value         = lumi
#  xsec_value         = c.GetLeaf("xsec").GetValue(n)
#  sumofweights_value = sum(sample['totalweight'])
#  return (genweight_value*lumi_value*xsec_value)/sumofweights_value

def getCutYieldFromChain(c, cutString = "(1)", cutFunc = None, weight = "weight", weightFunc = None, returnVar=False):
  c.Draw(">>eList", cutString)
  elist = ROOT.gDirectory.Get("eList")
  number_events = elist.GetN()
  res = 0.
  resVar=0.
  for i in range(number_events): 
    c.GetEntry(elist.GetEntry(i))
    if (not cutFunc) or cutFunc(c):
      if weight:
        w = c.GetLeaf(weight).GetValue()
      else:
        w=1.
      if weightFunc:
        w*=weightFunc(c)
      res += w
      resVar += w**2
  del elist
  if returnVar:
    return res, resVar
  return res

def getYieldFromChain(c, cutString = "(1)", weight = "weight", returnError=False):
  h = ROOT.TH1D('h_tmp', 'h_tmp', 1,0,2)
  h.Sumw2()
  c.Draw("1>>h_tmp", "("+weight+")*("+cutString+")", 'goff')
  res = h.GetBinContent(1)
  resErr = h.GetBinError(1)
#  print "1>>h_tmp", weight+"*("+cutString+")",res,resErr
  del h
  if returnError:
    return res, resErr
  return res 

def getPlotFromChain(c, var, binning, cutString = "(1)", weight = "weight", binningIsExplicit=False, addOverFlowBin=''):
  if binningIsExplicit:
    h = ROOT.TH1D('h_tmp', 'h_tmp', len(binning)-1, array('d', binning))
#    h.SetBins(len(binning), array('d', binning))
  else:
    if len(binning)==6:
      h = ROOT.TH2D('h_tmp', 'h_tmp', *binning)
    else:
      h = ROOT.TH1D('h_tmp', 'h_tmp', *binning)
  c.Draw(var+">>h_tmp", weight+"*("+cutString+")", 'goff')
  res = h.Clone()
  h.Delete()
  del h
  if addOverFlowBin.lower() == "upper" or addOverFlowBin.lower() == "both":
    nbins = res.GetNbinsX()
#    print "Adding", res.GetBinContent(nbins + 1), res.GetBinError(nbins + 1)
    res.SetBinContent(nbins , res.GetBinContent(nbins) + res.GetBinContent(nbins + 1))
    res.SetBinError(nbins , sqrt(res.GetBinError(nbins)**2 + res.GetBinError(nbins + 1)**2))
  if addOverFlowBin.lower() == "lower" or addOverFlowBin.lower() == "both":
    res.SetBinContent(1 , res.GetBinContent(0) + res.GetBinContent(1))
    res.SetBinError(1 , sqrt(res.GetBinError(0)**2 + res.GetBinError(1)**2))
  return res

## FIXME: Doesn't do what the name says ... plz move to private if really needed
#def genmatching(lepton,genparticles):
#  for gen in genparticles:
#      deltaphi = abs(lepton['phi'] - gen['phi'])
#      if (deltaphi > pi): deltaphi = 2*pi - deltaphi
#      deltaeta = abs(lepton['eta'] - gen['eta'])
#      deltar = sqrt(deltaphi**2 + deltaeta**2)
#      if deltar<0.01:
#        print deltar
#        print gen['motherId']
