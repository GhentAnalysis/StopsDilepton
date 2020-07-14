'''
/afs/hephy.at/work/p/phussain/StopsDileptonLegacy/results/v3//COMBINED/controlAll//cardFiles/T2tt/observed/
Get a signal region plot from the cardfiles
'''

#!/usr/bin/env python
from optparse import OptionParser
parser = OptionParser()
parser.add_option("--signal",               dest='signal',  action='store', default='T2tt',    choices=["T2tt", "T2bW", "ttHinv"], help="which signal?")
parser.add_option("--massPoints",           dest='massPoints',  action='store', default='800_100,350_150', help="which masspoints??")
parser.add_option("--version",              dest='version', action='store', default='v8',    help="Which version of estimates should be used?")
parser.add_option("--channel",              dest='channel',  action='store', default='all', choices=['all','OF','SF'], help="which channel??")
parser.add_option("--small",                action='store_true', help="small?")
parser.add_option('--logLevel',             dest="logLevel",              default='INFO',              action='store',      help="log level?", choices=['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'TRACE', 'NOTSET'])
parser.add_option('--blinded',              action="store_true")
parser.add_option('--overwrite',            dest="overwrite", default = False, action = "store_true", help="Overwrite existing output files, bool flag set to True  if used")
parser.add_option('--postFit',              dest="postFit", default = False, action = "store_true", help="Apply pulls?")
parser.add_option('--signalPostFit',        dest="signalPostFit", default = False, action = "store_true", help="Apply pulls to signal too?")
parser.add_option('--expected',             action = "store_true", help="Run expected?")
parser.add_option('--preliminary',          action = "store_true", help="Run expected?")
parser.add_option('--combined',             action = "store_true", help="combined fit for all years?")
parser.add_option('--extractSF',            action = "store_true", help="extract the SFs for backgrounds?")
parser.add_option('--testGrayscale',        action = "store_true", help="Do the most important test for this collaboration?")
parser.add_option('--splitBosons',          action = "store_true", help="Split multiboson component?")
parser.add_option('--signalOnly',           action = "store_true", help="Show only signals?")
parser.add_option("--year",                 action='store',      default=0, type="int", help='Which year?')
parser.add_option("--region",               action='store',      default="controlAll", choices=['fitAll', 'controlAll', 'signalOnly', 'controlDYVV'], help='Which year?')
parser.add_option("--postFix",              action='store',      default="", help='Add sth?')
(options, args) = parser.parse_args()

# Standard imports
import ROOT
import os
import sys
import pickle
import math
import yaml

# Analysis
from StopsDilepton.tools.u_float           import u_float
from StopsDilepton.tools.user              import plot_directory, analysis_results
from StopsDilepton.samples.color           import color
from StopsDilepton.tools.getPostFit        import *
from StopsDilepton.analysis.regions        import *

if options.region == 'controlAll':
    pass
elif options.region == 'signalOnly':
    regions = regionsLegacy[1:] + regionsLegacy[1:]
elif options.region == 'fitAll':
    regions =  [0]*45 #noRegions + noRegions +regionsLegacy + regionsLegacy

from RootTools.core.standard import *
# logger
import StopsDilepton.tools.logger as logger
import RootTools.core.logger as logger_rt
logger    = logger.get_logger(   options.logLevel, logFile = None)
logger_rt = logger_rt.get_logger(options.logLevel, logFile = None)

# get a setup
from StopsDilepton.analysis.Setup import Setup
if options.combined and options.year==0:
    setup=Setup(2016)
    lumiStr = 137
else:
    setup = Setup(options.year)
    lumiStr = setup.dataLumi/1000
#analysis_results = '/afs/hephy.at/work/p/phussain/StopsDileptonLegacy/results/v3/'
isData = True if not options.expected else False
#lumiStr = setup.dataLumi/1000
years=[2016,2017,2018]
controlRegions = options.region
massPoints = options.massPoints.split(',')

## ttH example: /afs/hephy.at/data/cms05/StopsDileptonLegacy/results/v7/2018/fitAll/cardFiles/ttHinv/observed/ttH_HToInvisible_M125_shapeCard.txt
#cardName = "%s_%s_shapeCard"%(options.signal,massPoints[0]) if not options.signal == 'ttHinv' else "ttH_HToInvisible_M125_shapeCard"
cardName = "%s_%s_shapeCard"%(options.signal,massPoints[0]) if not options.signal == 'ttHinv' else "ttH_HToInvisible_M125_shapeCard"

analysis_results = analysis_results.replace('v8', options.version)

inSignalRegions = not options.region.count('control')>0
if inSignalRegions and len(massPoints)>1:
    cardName2 = "%s_%s_shapeCard"%(options.signal,massPoints[1])
if options.combined:
    cardDir = analysis_results+"/COMBINED/%s/cardFiles/%s/%s/"%(controlRegions,options.signal,'expected' if options.expected else 'observed')
else:
    cardDir = analysis_results+"/%s/%s/cardFiles/%s/%s/"%(options.year,controlRegions,options.signal,'expected' if options.expected else 'observed')

cardFile = "%s/%s.txt"%(cardDir, cardName)
if inSignalRegions and len(massPoints)>1:
    cardFile2 = "%s/%s.txt"%(cardDir, cardName2)


logger.info("Plotting from cardfile %s"%cardFile)

# get the results
postFitResults = getPrePostFitFromMLF(cardFile.replace('.txt','_FD.root'))
print cardFile
if inSignalRegions and len(massPoints)>1:
    print cardFile2
    postFitResults2 = getPrePostFitFromMLF(cardFile2.replace('.txt','_FD.root'))
else: postFitResults2 = False

print cardFile.replace('.txt','_FD.root')
covariance = getCovarianceFromMLF(cardFile.replace('.txt','_FD.root'), postFit=options.postFit)

'''
Hard-code uncertainties on rate parameters now. Unclear how to extract them from combine workspace. Can be read when running FitDiagnostics with verbosity -v 3
'''
rateParams = yaml.load(file('rateParams.yaml','r'))

preFitHist={}
postFitHist={}
bhistos=[]
hists={}
histos={}
bkgHist=[]
if not options.splitBosons:
    processes = [   ('TTJets', 't#bar{t}/single t'),
                ('DY', 'Drell-Yan'),
                ('multiBoson', 'Multiboson'),
                ('TTZ', 't#bar{t}Z'),
                ('TTXNoZ', 'Rare')]
else:
    processes = [   ('TTJets', 't#bar{t}/t'),
                ('DY', 'Drell-Yan'),
                ('diBoson', 'WZ+others'),
                ('WW', 'WW(ll#nu#nu)'),
                ('ZZ', 'ZZ(ll#nu#nu)'),
                ('triBoson', 'VVV'),
                ('TTZ', 't#bar{t}Z'),
                ('TTXNoZ', 't#bar{t}X, rare')]

if options.year >0 : years = [options.year]

SFs = {}

if options.combined:
    for year in years:
        preFitHist[year]        = postFitResults['hists']['shapes_prefit']['dc_%s'%year]
        postFitHist[year]       = postFitResults['hists']['shapes_fit_s']['dc_%s'%year]
        hists[year] = preFitHist[year] if not options.postFit else postFitHist[year]

        # signal is always prefit for the plots
        if options.postFit and options.signalPostFit:
            hists[year]['signal1'] = postFitResults['hists']['shapes_fit_s']['dc_%s'%year]['signal']
            if inSignalRegions and postFitResults2: hists[year]['signal2'] = postFitResults2['hists']['shapes_fit_s']['dc_%s'%year]['signal']
        else:
            hists[year]['signal1'] = postFitResults['hists']['shapes_prefit']['dc_%s'%year]['signal']
            if inSignalRegions and postFitResults2: hists[year]['signal2'] = postFitResults2['hists']['shapes_prefit']['dc_%s'%year]['signal']

    for i,(p,tex) in enumerate(processes):
        bhistos.append( hists[years[0]][p])
    dataHist = hists[years[0]]['DY'].Clone()
    dataHist.Reset()
    dataHist.SetName('data')
    dataHist.legendText = 'Observed'
    #dataHist.drawOption = 'e0'
    for n,(p,tex) in enumerate(processes):
        for i in range(bhistos[n].GetNbinsX()):
            v = 0
            v = sum([hists[x][p].GetBinContent(i+1) for x in years])
            bhistos[n].SetBinContent(i+1, v)
        if tex:
            bhistos[n].legendText = tex
        #SFs[p] = {'val': sum([postFitHist[x][p].Integral()/preFitHist[x][p].Integral() for x in years]), 'sigma': math.sqrt(covariance['%s_norm'%p]['%s_norm'%p])}
        try:
            SFs[p+'_%s'%year] = {'val': sum([postFitHist[x][p].Integral()/preFitHist[x][p].Integral() for x in years]), 'sigma': rateParams[p+'_%s'%year]['sigma']}
        except KeyError:
            print "No rate parameter found for %s, %s"%(p, year)
            SFs[p+'_%s'%year] = {'val': 1, 'sigma': 0}

        histos[p]=bhistos[n]
        bhistos[n].style = styles.fillStyle( getattr(color, p), lineColor=getattr(color,p), errors=False )
        histos[p].style = styles.fillStyle( getattr(color, p), lineColor=getattr(color,p), errors=False )
        bkgHist.append( bhistos[n])

    dataHist.SetBinErrorOption(ROOT.TH1F.kPoisson)
    for i in range(dataHist.GetNbinsX()):
        dataHist.SetBinContent(i+1, sum([hists[x]['data'].Eval(i+0.5) for x in years]))

    histos['data'] = dataHist
    histos['data'].style = styles.errorStyle( ROOT.kBlack, markerSize = 1., drawOption='e0' )
    histos['data'].legendOption = 'p'

    signalHist = hists[years[0]]['DY'].Clone()
    signalHist.Reset()
    signalHist.SetName('signal')
    massPoint1 = [ x for x in massPoints[0].split('_') if not x.isalpha() ]
    signalHist.legendText = options.signal+' (%s,%s)'%(tuple(massPoint1))
    if options.signal == 'ttHinv': signalHist.legendText = "ttH, B(H#rightarrow inv)=100%"
    for year in years:
        print year, hists[year]['signal1'].Integral()
    for i in range(signalHist.GetNbinsX()):
        signalHist.SetBinContent(i+1, sum([hists[x]['signal1'].GetBinContent(i+1) for x in years]))
    histos['signal1'] = signalHist
    histos['signal1'].style = styles.lineStyle( ROOT.kBlack, width=2 )

    if inSignalRegions and postFitResults2:
        signalHist2 = hists[years[0]]['DY'].Clone()
        signalHist2.Reset()
        signalHist2.SetName('signal')
        massPoint2 = [ x for x in massPoints[1].split('_') if not x.isalpha() ]
        signalHist2.legendText = options.signal+' (%s,%s)'%(tuple(massPoint2))
        for year in years:
            print year, hists[year]['signal2'].Integral()
        for i in range(signalHist2.GetNbinsX()):
            signalHist2.SetBinContent(i+1, sum([hists[x]['signal2'].GetBinContent(i+1) for x in years]))
        histos['signal2'] = signalHist2
        histos['signal2'].style = styles.lineStyle( ROOT.kBlack, width=2, dashed=True )

else:
    preFitHists     = postFitResults['hists']['shapes_prefit']['Bin0']
    postFitHists    = postFitResults['hists']['shapes_fit_b']['Bin0']
    hists           = preFitHists if not options.postFit else postFitHists

    bkgHists = []

    dataHist = hists['DY'].Clone()
    dataHist.Reset()
    dataHist.SetName('data')
    dataHist.legendText = 'Data'

    for p,tex in processes:
        hists[p].style = styles.fillStyle( getattr(color, p), lineColor=getattr(color,p), errors=False )
        bkgHists.append(hists[p])
        if tex:
            hists[p].legendText = tex
        
    dataHist.SetBinErrorOption(ROOT.TH1F.kPoisson)
    for i in range(dataHist.GetNbinsX()):
        dataHist.SetBinContent(i+1, hists['data'].Eval(i+0.5))
        dataHist.SetBinError(i+1, math.sqrt(hists['data'].Eval(i+0.5)))

    hists['data'] = dataHist
    hists['data'].style = styles.errorStyle( ROOT.kBlack, markerSize = 1., drawOption='e0' )
    hists['data'].legendOption = 'p'

    if options.signalOnly:
        hists['total_background'].style = styles.lineStyle( ROOT.kOrange+1, width=2, errors=True)
        hists['total_background'].legendText = "Total SM"

    print postFitResults['hists']['shapes_prefit']['Bin0'].keys()
    if options.postFit and options.signalPostFit:
        hists['signal1'] = postFitResults['hists']['shapes_fit_s']['Bin0']['signal']
    else:
        hists['signal1'] = postFitResults['hists']['shapes_prefit']['Bin0']['signal']
    hists['signal1'].style = styles.lineStyle( ROOT.kBlack, width=2 )
    massPoint1 = [ x for x in massPoints[0].split('_') if not x.isalpha() ]
    hists['signal1'].legendText = options.signal+' (%s,%s)'%(tuple(massPoint1))
    if options.signal == 'ttHinv': hists['signal1'].legendText = "ttH, B(H#rightarrow inv)=100%"

    if inSignalRegions and postFitResults2:
        if options.postFit and options.signalPostFit:
            hists['signal2'] = postFitResults2['hists']['shapes_fit_s']['Bin0']['signal']
        else:
            hists['signal2'] = postFitResults2['hists']['shapes_prefit']['Bin0']['signal']
        hists['signal2'].style = styles.lineStyle( ROOT.kBlack, width=2, dashed=True )
        massPoint2 = [ x for x in massPoints[1].split('_') if not x.isalpha() ]
        hists['signal2'].legendText = options.signal+' (%s,%s)'%(tuple(massPoint2))    

if options.extractSF:
    SF_file = analysis_results + 'SF.pkl'
    if os.path.isfile(SF_file):
        SF_dict = pickle.load(file(SF_file, 'r'))
    else:
        SF_dict = {}
    
    print SFs
    
    SF_dict.update(SFs)
    pickle.dump(SF_dict, file(SF_file, 'w'))
    yaml.dump(SF_dict, file('SFs.yaml', 'w'))

boxes = []
ratio_boxes = []

signal_boxes = []
signal_ratio_boxes = []

if options.combined:
    for ib in range(1, 1 + hists[years[0]]['total_background'].GetNbinsX() ):
        val = sum([hists[x]['total_background'].GetBinContent(ib) for x in years])
        if val<0: continue
        variance = sum( [ covariance['dc_%s_%s'%(comb[0], ib-1)]['dc_%s_%s'%(comb[1], ib-1)] for comb in itertools.combinations_with_replacement(years,2) ] )
        sys = math.sqrt(variance)
        if val > 0:
            sys_rel = sys/val
        else:
            sys_rel = 1.
        
        # uncertainty box in main histogram
        box = ROOT.TBox( hists[years[0]]['total_background'].GetXaxis().GetBinLowEdge(ib),  max([0.006, val-sys]), hists[years[0]]['total_background'].GetXaxis().GetBinUpEdge(ib), max([0.006, val+sys]) )
        box.SetLineColor(ROOT.kGray+1)
        box.SetFillStyle(3244)
        box.SetFillColor(ROOT.kGray+1)
        
        # uncertainty box in ratio histogram
        r_box = ROOT.TBox( hists[years[0]]['total_background'].GetXaxis().GetBinLowEdge(ib),  max(0.11, 1-sys_rel), hists[years[0]]['total_background'].GetXaxis().GetBinUpEdge(ib), min(1.9, 1+sys_rel) )
        r_box.SetLineColor(ROOT.kGray+1)
        r_box.SetFillStyle(3244)
        r_box.SetFillColor(ROOT.kGray+1)

        boxes.append( box )
        hists[years[0]]['total_background'].SetBinError(ib, 0)
        ratio_boxes.append( r_box )
        print "Bin %s: %s +/- %s"%(ib, val, sys)
else:

    for ib in range(1, 1 + hists['total_background'].GetNbinsX() ):
        val = hists['total_background'].GetBinContent(ib)
        if val<0: continue
        sys = hists['total_background'].GetBinError(ib)
        if val > 0:
            sys_rel = sys/val
        else:
            sys_rel = 1.
        
        # uncertainty box in main histogram
        box = ROOT.TBox( hists['total_background'].GetXaxis().GetBinLowEdge(ib),  max([0.006, val-sys]), hists['total_background'].GetXaxis().GetBinUpEdge(ib), max([0.006, val+sys]) )
        box.SetLineColor(ROOT.kGray+1)
        box.SetFillStyle(3244)
        box.SetFillColor(ROOT.kGray+1)
        
        # uncertainty box in ratio histogram
        r_box = ROOT.TBox( hists['total_background'].GetXaxis().GetBinLowEdge(ib),  max(0.11, 1-sys_rel), hists['total_background'].GetXaxis().GetBinUpEdge(ib), min(1.9, 1+sys_rel) )
        r_box.SetLineColor(ROOT.kGray+1)
        r_box.SetFillStyle(3244)
        r_box.SetFillColor(ROOT.kGray+1)

        boxes.append( box )
        if not options.signalOnly:
            hists['total_background'].SetBinError(ib, 0)
        ratio_boxes.append( r_box )
        
        # signal uncertainties
        val = hists['signal1'].GetBinContent(ib)
        if val<0: continue
        sys = hists['signal1'].GetBinError(ib)
        if val > 0:
            sys_rel = sys/val
        else:
            sys_rel = 1.
        
        # uncertainty box in main histogram
        signal_box = ROOT.TBox( hists['signal1'].GetXaxis().GetBinLowEdge(ib),  max([0.006, val-sys]), hists['signal1'].GetXaxis().GetBinUpEdge(ib), max([0.006, val+sys]) )
        signal_box.SetLineColor(ROOT.kGray+1)
        signal_box.SetFillStyle(3244)
        signal_box.SetFillColor(ROOT.kGray+1)
        
        signal_box.SetFillColor(ROOT.kGray+1)
        
        # uncertainty box in ratio histogram
        signal_r_box = ROOT.TBox( hists['signal1'].GetXaxis().GetBinLowEdge(ib),  max(0.11, 1-sys_rel), hists['signal1'].GetXaxis().GetBinUpEdge(ib), min(1.9, 1+sys_rel) )
        signal_r_box.SetLineColor(ROOT.kGray+1)
        signal_r_box.SetFillStyle(3244)
        signal_r_box.SetFillColor(ROOT.kGray+1)

        signal_boxes.append( signal_box )
        signal_ratio_boxes.append( signal_r_box )


def drawObjects( isData=False, lumi=36. ):
    tex = ROOT.TLatex()
    tex.SetNDC()
    tex.SetTextSize(0.05)
    tex.SetTextAlign(11) # align right
    lines = [
      (0.08, 0.945, 'CMS Simulation') if not isData else ( (0.08, 0.945, 'CMS') if not options.preliminary else (0.08, 0.945, 'CMS #bf{#it{Preliminary}}')),
      (0.80, 0.945, '#bf{%s fb^{-1} (13 TeV)}'%lumi )
    ]
    return [tex.DrawLatex(*l) for l in lines]

binLabels = [
    'TT CR SF',
    'TT CR DF',
    'TTZ 2j2b',
    'TTZ 3j1b',
    'TTZ 3j2b',
    'TTZ 4j1b',
    'TTZ 4j2b',]

binLabels += [ 'CR%s'%i for i in range(13) ]
for i in range(13):
    binLabels.append('SR%s SF'%i)
    binLabels.append('SR%s DF'%i)
#binLables += [ 'SR%s SF'%i for i in range(13) ]

def setBinLabels( hist ):
    for i in range(1, hist.GetNbinsX()+1):
        hist.GetXaxis().SetBinLabel(i, binLabels[i-1])
#        hist.GetXaxis().SetBinLabel(i, "   %s"%((i+1)/2 if i%2==1 else ''))

def drawDivisions(regions):
    #print len(regions)
    min = 0.08
    max = 0.95
    diff = (max-min) / len(regions)
    lines = []
    lines2 = []
    line = ROOT.TLine()
#   line.SetLineColor(38)
    line.SetLineWidth(1)
    line.SetLineStyle(2)
    lines  = [ (min+2*diff,  0.015, min+2*diff, 0.90) ]
    lines += [ (min+7*diff,  0.015, min+7*diff, 0.90) ]
    lines += [ (min+20*diff,  0.015, min+20*diff, 0.90) ]
    lines += [ (min+32*diff,  0.015, min+32*diff, 0.90) ]
    lines += [ (min+44*diff,  0.015, min+44*diff, 0.90) ]
    #lines = [ (min+4*i*diff,  0.015, min+4*i*diff, 0.85) if min+4*i*diff<0.74 else (min+4*i*diff,  0.013, min+4*i*diff, 0.54) for i in range(1,len(regions)/4 + 1) ]
    #lines += [ (min+4*3*diff, 0.85, min+4*3*diff, 0.93 ), (min+4*6*diff, 0.85, min+4*6*diff, 0.93 )]
    return [line.DrawLineNDC(*l) for l in lines] + [tex.DrawLatex(*l) for l in []] + [tex2.DrawLatex(*l) for l in lines2]

def drawLabels( regions ):
    tex = ROOT.TLatex()
    tex.SetNDC()
    tex.SetTextSize(0.042)
    tex.SetTextAngle(0)
    tex.SetTextAlign(12) # align right
    min = 0.08
    max = 0.95
    diff = (max-min) / len(regions)
    lines  = [ (min + 0.41, 0.85, "#bf{100 < M_{T2}(ll) < 140 GeV}"), (min + 0.65, 0.85, "#bf{140 < M_{T2}(ll) < 240 GeV}") ] 
    lines += [ (min + 0.04, 0.85, "#bf{N_{l} = 3, on-Z}"), (min + 0.18, 0.85, "#bf{N_{l} = 2, N_{b} = 0, on-Z}") ] 
    #lines =  [(min+(8*i+0.90)*diff, 0.850,  "M_{T2}(ll)=3")   for i, r in enumerate(regions[:-4][::4])]
    return [tex.DrawLatex(*l) for l in lines] if len(regions)>12 else []


def drawTexLabels( regions ):
    tex = ROOT.TMathText()
    tex.SetNDC()
    tex.SetTextSize(0.042)
    tex.SetTextAngle(0)
    tex.SetTextAlign(12) # align right
    min = 0.08
    max = 0.95
    diff = (max-min) / len(regions)
    lines  = [ (min + 0.41, 0.85, "\\mathrm{100 < M_{T2}(ll) < 140 GeV}"), (min + 0.65, 0.85, "\\mathrm{140 < M_{T2}(ll) < 240 GeV}") ] 
    lines += [ (min + 0.045, 0.85, "\\mathrm{N_{l} = 3, on-Z}"), (min + 0.20, 0.85, "\\mathrm{N_{l} = 2, N_{b} = 0, on-Z}") ] 
    #lines =  [(min+(8*i+0.90)*diff, 0.850,  "M_{T2}(ll)=3")   for i, r in enumerate(regions[:-4][::4])]
    return [tex.DrawMathText(*l) for l in lines] if len(regions)>12 else []

def drawLabelsRot( regions ):
    tex = ROOT.TLatex()
    tex.SetNDC()
    tex.SetTextSize(0.035)
    tex.SetTextAngle(90)
    tex.SetTextAlign(12) # align right
    min = 0.08
    max = 0.95
    diff = (max-min) / len(regions)
    #lines = [(min+(i*4+2.7)*diff, 0.545 if i<3 else 0.285,  r.texStringForVar('dl_mt2blbl')) for i, r in enumerate(regions[::2]) if i < 6]
    lines  = [(min+(44.8)*diff, 0.66, "#bf{M_{T2}(ll) > 240 GeV}")]
    lines += [(min+(1.5)*diff-0.002, 0.66, "#bf{M_{T2}(ll) < 100 GeV}")]
    return [tex.DrawLatex(*l) for l in lines] 

lumiStr = round(lumiStr,1) if not options.combined else int(lumiStr)

def getLegendRight():
    leg = ROOT.TLegend(0.71,0.78-9*0.052, 0.92, 0.78)
    leg.SetFillColor(ROOT.kWhite)
    leg.SetShadowColor(ROOT.kWhite)
    leg.SetBorderSize(0)
    leg.SetTextSize(0.045)
    for proc in processes:
        histos[proc[0]].SetFillColor(getattr(color, proc[0]))
        histos[proc[0]].SetLineColor(getattr(color, proc[0]))
        leg.AddEntry(histos[proc[0]], proc[1], 'f' )
    histos['signal1'].SetLineColor(1)
    histos['signal1'].SetLineWidth(2)
    histos['signal1'].SetLineStyle(1)
    leg.AddEntry(histos['signal1'], options.signal+' (%s,%s)'%(tuple(massPoint1)), 'l')
    histos['signal2'].SetLineColor(1)
    histos['signal2'].SetLineWidth(2)
    histos['signal2'].SetLineStyle(2)
    leg.AddEntry(histos['signal2'], options.signal+' (%s,%s)'%(tuple(massPoint2)), 'l')
    histos['data'].SetLineColor(1)
    histos['data'].SetLineWidth(1)
    histos['data'].SetMarkerStyle(8)
    leg.AddEntry(histos['data'], 'Observed', 'e1p')
    boxes[0].SetLineWidth(0)
    leg.AddEntry(boxes[0], 'Uncertainty', 'f')
    return [leg]

if options.signalOnly:
    drawObjects = drawObjects( isData=isData, lumi=lumiStr ) + signal_boxes + drawDivisions( regions ) + drawLabelsRot( regions ) + drawLabels( regions )
else:
    drawObjects = drawObjects( isData=isData, lumi=lumiStr ) + boxes + drawDivisions( regions ) + drawLabelsRot( regions ) + drawLabels( regions ) + getLegendRight()

if options.region == 'signalOnly':
    drawObjects += drawDivisions( regions ) + drawLabels( regions ) + drawLabelsRot( regions )

if options.combined:
    if inSignalRegions and len(massPoints)>1:
        plots = [ bkgHist, [histos['data']], [histos['signal1']], [histos['signal2']]]
    else:
        plots = [ bkgHist, [histos['data']], [histos['signal1']]]
    print "Total signal yield:", histos['signal1'].Integral(), histos['signal2'].Integral()
        
else:
    if options.signalOnly:
        plots = [ [hists['signal1']], [hists['signal2']], [hists['data']], [hists['total_background']] ]
    else:
        plots = [ bkgHists, [hists['data']], [hists['signal1']],[hists['signal2'] ]] if postFitResults2 else [ bkgHists, [hists['data']], [hists['signal1']] ]
    print "Total signal yield:", hists['signal1'].Integral(), hists['signal2'].Integral()
    



plotName = options.region
if options.combined: plotName += "_COMBINED"
if options.year > 0: plotName += "_%s"%options.year

if options.postFit:
    plotName += '_postFit'

if options.signal is not "T2tt":
    plotName += '_%s'%options.signal

if options.postFix:
    plotName += '_%s'%options.postFix

if options.signalOnly:
    plotName += '_signalOnly'

if options.combined and options.year==0:
    yMax = 90000000.
else:
    yMax = 900000.

canvasModifications = [ 
    lambda c : c.SetLeftMargin(0.08),
    lambda c : c.GetPad(2).SetLeftMargin(0.08),
    lambda c : c.GetPad(1).SetLeftMargin(0.08),
    lambda c : c.GetPad(2).SetBottomMargin(0.60),
    lambda c : c.GetPad(1).SetRightMargin(0.03), 
    lambda c : c.GetPad(2).SetRightMargin(0.03) 
    ]

if options.testGrayscale:
    canvasModifications += [ lambda c : c.SetGrayscale() ]


plotting.draw(
    Plot.fromHisto(plotName,
                plots,
                texX = "",
                texY = 'Number of events',
            ),
    plot_directory = os.path.join(plot_directory, "controlRegions", options.version),
    logX = False, logY = True, sorting = False, 
    legend = None,
    #legend = (0.75,0.80-0.010*32, 0.95, 0.80),
    #legend = (0.71,0.38, 0.92, 0.78) if not options.signalOnly else (0.70,0.55, 0.92, 0.75),
    widths = {'x_width':1300, 'y_width':600, 'y_ratio_width':250},
    yRange = (0.2,yMax),
    #yRange = (0.03, [0.001,0.5]), ### need to go up to 3.19 if we want to show full excess!! ##
    ratio = {'yRange': (0.11, 2.19), 'texY':'Data/Pred.', 'histos':[(1,0)], 'drawObjects':ratio_boxes if not options.signalOnly else signal_ratio_boxes, #+ drawLabelsLower( regions ) +drawHeadlineLower( regions ) + drawDivisionsLower(regions),
            'histModifications': [lambda h: setBinLabels(h), lambda h: h.GetYaxis().SetTitleSize(32), lambda h: h.GetYaxis().SetLabelSize(28), lambda h: h.GetYaxis().SetTitleOffset(1.0), lambda h: h.GetXaxis().SetTitleSize(32), lambda h: h.GetXaxis().SetLabelSize(27), lambda h: h.GetXaxis().SetLabelOffset(0.035), lambda h: h.LabelsOption('v'), lambda h: h.GetYaxis().SetTickLength(0.035), lambda h: h.GetXaxis().SetTickLength(0.02)]} ,
    drawObjects = drawObjects,
    histModifications = [lambda h: h.GetYaxis().SetTitleSize(32), lambda h: h.GetYaxis().SetLabelSize(28), lambda h: h.GetYaxis().SetTitleOffset(1.0), lambda h: h.GetYaxis().SetTickLength(0.015)],
    canvasModifications = canvasModifications,
    copyIndexPHP = True,
)

