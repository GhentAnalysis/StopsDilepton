#Standard import
import copy
import os

# RootTools
from RootTools.core.standard import *

# Logging
import logging
logger = logging.getLogger(__name__)

#user specific
from StopsDilepton.tools.user       import analysis_results
from StopsDilepton.tools.helpers    import getObjFromFile

from StopsDilepton.analysis.SystematicEstimator         import jmeVariations, metVariations
from StopsDilepton.analysis.SetupHelpers                import getZCut, channels, allChannels, trilepChannels, allTrilepChannels
from StopsDilepton.analysis.fastSimGenMetReplacements   import fastSimGenMetReplacements
from Analysis.Tools.metFilters                          import getFilterCut

##define samples
#from StopsDilepton.samples.nanoTuples_Summer16_postProcessed            import *
#from StopsDilepton.samples.nanoTuples_Run2016_17Jul2018_postProcessed   import *
#from StopsDilepton.samples.nanoTuples_Fall17_postProcessed              import *
#from StopsDilepton.samples.nanoTuples_Run2017_31Mar2018_postProcessed   import *
#from StopsDilepton.samples.nanoTuples_Autumn18_postProcessed            import *
#from StopsDilepton.samples.nanoTuples_Run2018_PromptReco_postProcessed  import *
#Define defaults here
zMassRange            = 15
default_mllMin        = 20
default_metMin        = 0
default_metSigMin     = 12
default_zWindow       = "offZ"
default_dPhi          = True
default_triLep        = False
default_dPhiInv       = False
default_dPhiJetMet    = 0.25 # To fix: this only applies to the 2nd jet, the first jet is fixes see below
default_nJets         = (2, -1)   # written as (min, max)
default_nBTags        = (1, -1)
default_leptonCharges = "isOS"


class Setup:
    def __init__(self, year=2016, puWeight=False, puUpOrDown=False):
        self.analysis_results = analysis_results
        self.zMassRange       = zMassRange
        self.prefixes         = []
        self.externalCuts     = []
        self.year = year

        #Default cuts and requirements. Those three things below are used to determine the key in the cache!
        self.parameters   = {
            'mllMin':        default_mllMin,
            'metMin':        default_metMin,
            'metSigMin':     default_metSigMin,
            'zWindow':       default_zWindow,
            'dPhi':          default_dPhi,
            'dPhiInv':       default_dPhiInv,
            'nJets':         default_nJets,
            'nBTags':        default_nBTags,
            'leptonCharges': default_leptonCharges,
            'triLep':        default_triLep,
        }

        # can use arbitrary combinations of PU weights (minBias x-secs) for central value and uncertainties
        minBiasXSecs = ['VVUp','VUp','Up','','Down']
        if puWeight==False:
            puWeight = 'VUp' if self.year == 2018 else ''

        # default: 1 sigma up/down
        if puUpOrDown == False:
            self.puUpOrDown = [ minBiasXSecs[minBiasXSecs.index(puWeight)-1], minBiasXSecs[minBiasXSecs.index(puWeight)+1] ]
        else:
            self.puUpOrDown = puUpOrDown
        self.puWeight = 'reweightPU'+puWeight

        self.sys = {'weight':'weight', 'reweight':[self.puWeight, 'reweightLeptonSF', 'reweightBTag_SF','reweightLeptonTrackingSF', 'reweightDilepTrigger', 'reweightL1Prefire','reweightHEM', 'reweightLeptonHit0SF', 'reweightLeptonSip3dSF'], 'selectionModifier':None} 

        if year == 2016:
            #define samples
            from StopsDilepton.samples.nanoTuples_Summer16_postProcessed            import Top_pow_16, DY_HT_LO_16, TTZ_16, multiBoson_16,TTXNoZ_16, TZX_16, diBoson_ll_16, triboson_16, WWTo2L2Nu_16, ZZTo2L2Nu_16
            top         = Top_pow_16
            DY          = DY_HT_LO_16
            TTZ         = TTZ_16
            multiBoson  = multiBoson_16
            TTXNoZ      = TTXNoZ_16
            TZX         = TZX_16
            # for splitting multiboson
            diBoson     = diBoson_ll_16
            triBoson    = triboson_16
            WW          = WWTo2L2Nu_16
            ZZ          = ZZTo2L2Nu_16
            from StopsDilepton.samples.nanoTuples_Run2016_17Jul2018_postProcessed   import Run2016
            data        = Run2016
        elif year == 2017:
            from StopsDilepton.samples.nanoTuples_Fall17_postProcessed              import Top_pow_17, DY_HT_LO_17, TTZ_17, multiBoson_17, TTXNoZ_17, TZX_17, diBoson_ll_17, triboson_17,WWTo2L2Nu_17, ZZTo2L2Nu_17
            top         = Top_pow_17
            DY          = DY_HT_LO_17
            TTZ         = TTZ_17
            multiBoson  = multiBoson_17
            TTXNoZ      = TTXNoZ_17
            TZX         = TZX_17
            # for splitting multiboson
            diBoson     = diBoson_ll_17
            triBoson    = triboson_17
            WW          = WWTo2L2Nu_17
            ZZ          = ZZTo2L2Nu_17
            from StopsDilepton.samples.nanoTuples_Run2017_31Mar2018_postProcessed   import Run2017 
            data        = Run2017  
            #data        = Run2017BCDE
        elif year == 2018:
            from StopsDilepton.samples.nanoTuples_Autumn18_postProcessed            import Top_pow_18, DY_HT_LO_18, TTZ_18, multiBoson_18, TTXNoZ_18, TZX_18, diBoson_ll_18, triboson_18, WWTo2L2Nu_18, ZZTo2L2Nu_18
            top         = Top_pow_18
            DY          = DY_HT_LO_18
            TTZ         = TTZ_18
            multiBoson  = multiBoson_18
            TTXNoZ      = TTXNoZ_18
            TZX         = TZX_18
            # for splitting multiboson
            diBoson     = diBoson_ll_18
            triBoson    = triboson_18
            WW          = WWTo2L2Nu_18
            ZZ          = ZZTo2L2Nu_18
            from StopsDilepton.samples.nanoTuples_Run2018_PromptReco_postProcessed  import Run2018 
            data        = Run2018

        self.samples = {
            'Top_gaussian' :    top,
            'Top_nongaussian':  top,
            'Top_fakes' :       top,
            'DY':               DY,
            'TTJets' :          top,
            'TTZ' :             TTZ,
            'multiBoson' :      multiBoson,
            'TTXNoZ' :          TTXNoZ,
            'TZX' :             TZX,
            'diBoson' :         diBoson,
            'triBoson' :        triBoson,
            'WW' :              WW,
            'ZZ' :              ZZ,
            'other'  :          TTXNoZ,
            'Data'   :          data,
        }

        self.lumi     = data.lumi
        self.dataLumi = data.lumi # get from data samples later
        
        dataPUHistForSignalPath = "$CMSSW_BASE/src/StopsDilepton/tools/data/puFastSimUncertainty/dataPU.root"
        self.dataPUHistForSignal = getObjFromFile(os.path.expandvars(dataPUHistForSignalPath), "data")

    def prefix(self):
        return '_'.join(self.prefixes+[self.preselection('MC')['prefix']])

    def defaultCacheDir(self):
        cacheDir = os.path.join(self.analysis_results, str(self.year), self.prefix(), 'cacheFiles')
        logger.info('Default cache dir is: %s', cacheDir)
        return cacheDir

    #Clone the setup and optinally modify the systematic variation
    def sysClone(self, sys=None, parameters=None):
        '''Clone setup and change systematic if provided'''

        res            = copy.copy(self)
        res.sys        = copy.deepcopy(self.sys)
        res.parameters = copy.deepcopy(self.parameters)

        if sys:
            for k in sys.keys():
                if k=='remove':
                    for i in sys[k]:
                      res.sys['reweight'].remove(i)
                elif k=='reweight':
                    res.sys[k] = list(set(res.sys[k]+sys[k])) #Add with unique elements
                    
                    # need to treat PU seperately
                    #puUpOrDown = ['VVUp','VUp','Up','','Down'] # we can just loop over every possible case  # if self.year == 2018 else ['Up','Down']
                    #puUpOrDown = ['VVUp','Up'] # we can just loop over every possible case  # if self.year == 2018 else ['Up','Down']
                    for upOrDown in self.puUpOrDown:
                      #print upOrDown, res.sys[k]
                      if 'reweightPU'+upOrDown in res.sys[k]:
                        res.sys[k].remove(self.puWeight)
                        #break
                    # all the other weight systematics
                    for upOrDown in ['Up','Down']:
                      if 'reweightDilepTrigger'+upOrDown        in res.sys[k]: res.sys[k].remove('reweightDilepTrigger')
                      if 'reweightL1Prefire'+upOrDown           in res.sys[k]: res.sys[k].remove('reweightL1Prefire')
                      if 'reweightBTag_SF_b_'+upOrDown          in res.sys[k]: res.sys[k].remove('reweightBTag_SF')
                      if 'reweightBTag_SF_l_'+upOrDown          in res.sys[k]: res.sys[k].remove('reweightBTag_SF')
                      if 'reweightBTag_SF_FS_'+upOrDown         in res.sys[k]: res.sys[k].remove('reweightBTag_SF')
                      if 'reweightLeptonSF'+upOrDown            in res.sys[k]: res.sys[k].remove('reweightLeptonSF')
                      if 'reweightLeptonHit0SF'+upOrDown        in res.sys[k]: res.sys[k].remove('reweightLeptonHit0SF')
                      if 'reweightLeptonSip3dSF'+upOrDown       in res.sys[k]: res.sys[k].remove('reweightLeptonSip3dSF')
                      if 'reweightLeptonFastSimSF'+upOrDown     in res.sys[k]: res.sys[k].remove('reweightLeptonFastSimSF')
                      if 'reweight_nISR'+upOrDown               in res.sys[k]: res.sys[k].remove('reweight_nISR')
                else:
                    res.sys[k] = sys[k] # if sys[k] else res.sys[k]

        if parameters:
            for k in parameters.keys():
                res.parameters[k] = parameters[k]


        #print res.puWeight
        return res

    def defaultParameters(self, update={}):
        assert type(update)==type({}), "Update arguments with key arg dictionary. Got this: %r"%update
        res = copy.deepcopy(self.parameters)
        res.update(update)
        return res

    def weightString(self):
        return "*".join([self.sys['weight']] + (self.sys['reweight'] if self.sys['reweight'] else []))

    def preselection(self, dataMC , channel='all', isFastSim = False):
        '''Get preselection  cutstring.'''
        return self.selection(dataMC, channel = channel, isFastSim = isFastSim, hadronicSelection = False, **self.parameters)

    def selection(self, dataMC,
                        mllMin, metMin, metSigMin, zWindow, dPhi, dPhiInv,
                        nJets, nBTags, leptonCharges, triLep,
                        channel = 'all', hadronicSelection = False,  isFastSim = False):
        '''Define full selection
           dataMC: 'Data' or 'MC'
           channel: all, EE, MuMu or EMu
           zWindow: offZ, onZ, or allZ
           hadronicSelection: whether to return only the hadronic selection
       isFastSim: adjust filter cut etc. for fastsim
        '''
        #Consistency checks
        if self.sys['selectionModifier']:
          assert self.sys['selectionModifier'] in jmeVariations+metVariations+['GenMET'] or 'nVert' in self.sys['selectionModifier'] or 'MVA' in self.sys['selectionModifier'], "Don't know about systematic variation %r, take one of %s"%(self.sys['selectionModifier'], ",".join(jmeVariations + ['GenMET']))
        assert dataMC in ['Data','MC'],                                                   "dataMC = Data or MC, got %r."%dataMC
        assert not leptonCharges or leptonCharges in ["isOS", "isSS"],                    "Don't understand leptonCharges %r. Should take isOS or isSS."%leptonCharges

        #Postfix for variables (only for MC and if we have a jme variation)
        sysStr = ""
        metStr = ""
        if dataMC == "MC" and self.sys['selectionModifier'] in jmeVariations: sysStr = "_" + self.sys['selectionModifier']
        if dataMC == "MC" and self.sys['selectionModifier'] in metVariations: metStr = "_" + self.sys['selectionModifier']

        res={'cuts':[], 'prefixes':[]}

        if leptonCharges and not hadronicSelection and not triLep:
            res['cuts'].append(leptonCharges)
            res['prefixes'].append(leptonCharges)

        if nJets and not (nJets[0]==0 and nJets[1]<0):
            assert nJets[0]>=0 and (nJets[1]>=nJets[0] or nJets[1]<0), "Not a good nJets selection: %r"%nJets
            njetsstr = "nJetGood"+sysStr+">="+str(nJets[0])
            prefix   = "nJets"+str(nJets[0])
            if nJets[1]>=0:
                njetsstr+= "&&"+"nJetGood"+sysStr+"<="+str(nJets[1])
                if nJets[1]!=nJets[0]: prefix+=str(nJets[1])
            else:
                prefix+='p'
            res['cuts'].append(njetsstr)
            res['prefixes'].append(prefix)

        if nBTags and not (nBTags[0]==0 and nBTags[1]<0):
            assert nBTags[0]>=0 and (nBTags[1]>=nBTags[0] or nBTags[1]<0), "Not a good nBTags selection: %r"% nBTags
            nbtstr = "nBTag"+sysStr+">="+str(nBTags[0])
            prefix = "nbtag"+str(nBTags[0])
            if nBTags[1]>=0:
                nbtstr+= "&&nBTag"+sysStr+"<="+str(nBTags[1])
                if nBTags[1]!=nBTags[0]: prefix+=str(nBTags[1])
            else:
                prefix+='p'
            res['cuts'].append(nbtstr)
            res['prefixes'].append(prefix)

        if metMin and metMin>0:
          res['cuts'].append('met_pt'+sysStr+metStr+'>='+str(metMin))
          res['prefixes'].append('met'+str(metMin))
        if metSigMin and metSigMin>0:
          res['cuts'].append('MET_significance'+sysStr+metStr+'>='+str(metSigMin))
          #res['cuts'].append('metSig'+sysStr+metStr+'>='+str(metSigMin))
          res['prefixes'].append('METsig'+str(metSigMin))

        if dPhi:
          res['cuts'].append('cos(met_phi'+sysStr+metStr+'-JetGood_phi[0])<0.8&&cos(met_phi'+sysStr+metStr+'-JetGood_phi[1])<cos(0.25)')
          res['prefixes'].append('dPhiJet0-dPhiJet')
        elif dPhiInv:
          res['cuts'].append('!(cos(met_phi'+sysStr+metStr+'-JetGood_phi[0])<0.8&&cos(met_phi'+sysStr+metStr+'-JetGood_phi[1])<cos(0.25))')
          res['prefixes'].append('dPhiInv')

        if not hadronicSelection:
            if mllMin and mllMin>0:
              res['cuts'].append('dl_mass>='+str(mllMin))
              res['prefixes'].append('mll'+str(mllMin))

            if triLep:
              assert channel in allTrilepChannels, "channel must be one of "+",".join(allTrilepChannels)+". Got %r."%channel
              from StopsDilepton.tools.trilepSelection import getTrilepSelection
              if channel =="all": chStr = '(' + '||'.join([getTrilepSelection(c) for c in ['3mu','2mu1e','2e1mu','3e']]) + ')'
              else:               chStr = getTrilepSelection(channel)
              res['cuts'].append(chStr)

              assert zWindow == 'onZ', "zWindow must be onZ for trilep selection"
              res['cuts'].append('abs(Z1_mass-91.1876)<10')

            else:

              preselMuMu = "isMuMu==1&&nGoodMuons==2&&nGoodElectrons==0"
              preselEE   = "isEE==1&&nGoodMuons==0&&nGoodElectrons==2"
              preselEMu  = "isEMu==1&&nGoodMuons==1&&nGoodElectrons==1"

              #Z window
              assert zWindow in ['offZ', 'onZ', 'allZ'], "zWindow must be one of onZ, offZ, allZ. Got %r"%zWindow
              if zWindow == 'onZ':                     res['cuts'].append(getZCut(zWindow, self.zMassRange))
              if zWindow == 'offZ' and channel!="EMu": res['cuts'].append(getZCut(zWindow, self.zMassRange))  # Never use offZ when in emu channel, use allZ instead

              #lepton channel
              assert channel in allChannels, "channel must be one of "+",".join(allChannels)+". Got %r."%channel

              if channel=="MuMu":  chStr = preselMuMu
              elif channel=="EE":  chStr = preselEE
              elif channel=="EMu": chStr = preselEMu
              elif channel=="all": chStr = "("+preselMuMu+'||'+preselEE+'||'+preselEMu+')'
              elif channel=="SF":  chStr = "("+preselMuMu+'||'+preselEE+')'

              res['cuts'].append(chStr)

              res['prefixes'].append('looseLeptonMiniIsoVeto')
              #res['cuts'].append('Sum$(LepGood_pt>15&&LepGood_relIso03<0.4)==2')
              res['cuts'].append('(Sum$(Electron_pt>15&&abs(Electron_eta)<2.4&&Electron_miniPFRelIso_all<0.4) + Sum$(Muon_pt>15&&abs(Muon_eta)<2.4&&Muon_miniPFRelIso_all<0.4) )==2')
 
              res['prefixes'].append('miniIso0.2')
              res['cuts'].append("l1_miniRelIso<0.2&&l2_miniRelIso<0.2")

              # add HEMJetVetoWide in 2018
             # if self.year == 2018:
             #   res['prefixes'].append('HEMJetVetoWide')
             #   res['cuts'].append("Sum$(Jet_pt>20&&Jet_eta<-1.0&&Jet_eta>-3.2&&Jet_phi<-0.5&&Jet_phi>-2.0)==0")
              # add BadeeJetVeto in 2017
              if self.year == 2017:
                res['prefixes'].append('BadEEJetVeto')
                res['cuts'].append("Sum$((2.6<abs(Jet_eta)&&abs(Jet_eta)<3&&Jet_pt>30))==0")


              res['prefixes'].append('lepSel')
              res['cuts'].append("l1_pt>30&&l2_pt>20")
              
        #res['cuts'].append(getFilterCut(isData=(dataMC=='Data'), year=self.year, isFastSim=isFastSim, skipBadPFMuon=False, skipVertexFilter=True))
        res['cuts'].append(getFilterCut(isData=(dataMC=='Data'), year=self.year, isFastSim=isFastSim, skipBadPFMuon=False))
        if dataMC=='Data' and self.year == 2018:
            res['cuts'].append('reweightHEM>0')
        res['cuts'].extend(self.externalCuts)


         
        # for SUSY PU uncertainty of 2016 samples
        if self.sys['selectionModifier']:
            if "nVert" in self.sys['selectionModifier']:
                res['cuts'].append(self.sys['selectionModifier'])

        # for MVA studies
        if self.sys['selectionModifier']:
            if self.sys['selectionModifier'].count('MVA'):
                res['cuts'].append(self.sys['selectionModifier'])

        # for SUSY fast sim MET uncertainty
        if self.sys['selectionModifier'] == 'GenMET':
            res['cuts'] = [ fastSimGenMetReplacements(r) for r in res['cuts'] ]
        return {'cut':"&&".join(res['cuts']), 'prefix':'-'.join(res['prefixes']), 'weightStr': self.weightString()}
