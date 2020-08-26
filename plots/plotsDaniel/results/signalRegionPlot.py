'''
/afs/hephy.at/work/p/phussain/StopsDileptonLegacy/results/v3//COMBINED/controlAll//cardFiles/T2tt/observed/
Get a signal region plot from the cardfiles
'''

#!/usr/bin/env python
from optparse import OptionParser
parser = OptionParser()
parser.add_option("--signal",               dest='signal',  action='store', default='T2tt',    choices=["T2tt"], help="which signal?")
parser.add_option("--massPoints",           dest='massPoints',  action='store', default='800_100,350_150', help="which masspoints??")
parser.add_option("--version",              action='store', default='v8',            nargs='?',      help="Which version of estimates should be used?")
parser.add_option("--channel",              dest='channel',  action='store', default='all', choices=['all','OF','SF'], help="which channel??")
parser.add_option("--small",                action='store_true', help="small?")
parser.add_option('--logLevel',             dest="logLevel",              default='INFO',              action='store',      help="log level?", choices=['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'TRACE', 'NOTSET'])
parser.add_option('--blinded',              action="store_true")
parser.add_option('--overwrite',            dest="overwrite", default = False, action = "store_true", help="Overwrite existing output files, bool flag set to True  if used")
parser.add_option('--postFit',              dest="postFit", default = False, action = "store_true", help="Apply pulls?")
parser.add_option('--expected',             action = "store_true", help="Run expected?")
parser.add_option('--preliminary',             action = "store_true", help="Run expected?")
parser.add_option('--combined',             action = "store_true", help="combined fit for all years?")
parser.add_option("--year",                 action='store',      default=2017, type="int", help='Which year?')
parser.add_option("--region",               action='store',      default="controlAll", choices=['fitAll', 'controlAll', 'signalOnly', 'controlDYVV'], help='Which year?')
(options, args) = parser.parse_args()

# Standard imports
import ROOT
import os
import sys
import pickle
import math

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
if options.combined:
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
cardName = "T2tt_%s_shapeCard"%massPoints[0]

inSignalRegions = not options.region.count('control')>0
if inSignalRegions:
    cardName2 = "T2tt_%s_shapeCard"%massPoints[1]
if options.combined:
    cardDir = analysis_results.replace('v8','') +'/' + str(options.version)+"/COMBINED/%s/cardFiles/%s/%s/"%(controlRegions,options.signal,'expected' if options.expected else 'observed')
else:
    cardDir = analysis_results.replace('v8','') +'/' + str(options.version)+"/%s/%s/cardFiles/%s/%s/"%(options.year,controlRegions,options.signal,'expected' if options.expected else 'observed')

cardFile = "%s/%s.txt"%(cardDir, cardName)
if inSignalRegions:
    cardFile2 = "%s/%s.txt"%(cardDir, cardName2)


logger.info("Plotting from cardfile %s"%cardFile)

# get the results
postFitResults = getPrePostFitFromMLF(cardFile.replace('.txt','_FD.root'))
if inSignalRegions: postFitResults2 = getPrePostFitFromMLF(cardFile2.replace('.txt','_FD.root'))

covariance = getCovarianceFromMLF(cardFile.replace('.txt','_FD.root'), postFit=options.postFit)

#raise NotImplementedError()

preFitHist={}
postFitHist={}
bhistos=[]
hists={}
histos={}
bkgHist=[]
processes = [   ('TTJets', 't#bar{t}/t'),
                ('DY', 'Drell-Yan'),
                ('multiBoson', 'VV/VVV'),
                ('TTZ', 't#bar{t}Z'),
                ('TTXNoZ', 't#bar{t}X, rare')]

if options.combined:
    for year in years:
        preFitHist[year]        = postFitResults['hists']['shapes_prefit']['dc_%s'%year]
        postFitHist[year]       = postFitResults['hists']['shapes_fit_s']['dc_%s'%year]
        hists[year] = preFitHist[year] if not options.postFit else postFitHist[year]

        # signal is always prefit for the plots
        hists[year]['signal1'] = postFitResults['hists']['shapes_prefit']['dc_%s'%year]['signal']
        if inSignalRegions: hists[year]['signal2'] = postFitResults2['hists']['shapes_prefit']['dc_%s'%year]['signal']

    for i,(p,tex) in enumerate(processes):
        bhistos.append( hists[2016][p])
    dataHist = hists[2016]['DY'].Clone()
    dataHist.Reset()
    dataHist.SetName('data')
    dataHist.legendText = 'Data'
    for n,(p,tex) in enumerate(processes):
        for i in range(bhistos[n].GetNbinsX()):
            v=0
            v=hists[2016][p].GetBinContent(i+1) + hists[2017][p].GetBinContent(i+1) + hists[2018][p].GetBinContent(i+1)
            bhistos[n].SetBinContent(i+1, v)
        if tex:
            bhistos[n].legendText = tex
        histos[p]=bhistos[n]
        bhistos[n].style = styles.fillStyle( getattr(color, p), lineColor=getattr(color,p), errors=False )
        bkgHist.append( bhistos[n])

    dataHist.SetBinErrorOption(ROOT.TH1F.kPoisson)
    for i in range(dataHist.GetNbinsX()):
        dataHist.SetBinContent(i+1, (hists[2016]['data'].Eval(i+0.5) + hists[2017]['data'].Eval(i+0.5) + hists[2018]['data'].Eval(i+0.5)))

    histos['data'] = dataHist
    histos['data'].style = styles.errorStyle( ROOT.kBlack, markerSize = 1. )
    histos['data'].legendOption = 'p'

    signalHist = hists[2016]['DY'].Clone()
    signalHist.Reset()
    signalHist.SetName('signal')
    signalHist.legendText = 'T2tt (800,100)'
    for i in range(signalHist.GetNbinsX()):
        signalHist.SetBinContent(i+1, (hists[2016]['signal1'].GetBinContent(i+1) + hists[2017]['signal1'].GetBinContent(i+1) + hists[2018]['signal1'].GetBinContent(i+1)))
    histos['signal1'] = signalHist
    histos['signal1'].style = styles.lineStyle( ROOT.kBlack, width=2 )

    if inSignalRegions:
        signalHist2 = hists[2016]['DY'].Clone()
        signalHist2.Reset()
        signalHist2.SetName('signal')
        signalHist2.legendText = 'T2tt (350,150)'
        for i in range(signalHist2.GetNbinsX()):
            signalHist2.SetBinContent(i+1, (hists[2016]['signal2'].GetBinContent(i+1) + hists[2017]['signal2'].GetBinContent(i+1) + hists[2018]['signal2'].GetBinContent(i+1)))
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
        
    for i in range(dataHist.GetNbinsX()):
        dataHist.SetBinContent(i+1, hists['data'].Eval(i+0.5))
        dataHist.SetBinError(i+1, math.sqrt(hists['data'].Eval(i+0.5)))

    hists['data'] = dataHist
    hists['data'].style = styles.errorStyle( ROOT.kBlack, markerSize = 1. )
    hists['data'].legendOption = 'p'

    print postFitResults['hists']['shapes_prefit']['Bin0'].keys()
    hists['signal1'] = postFitResults['hists']['shapes_prefit']['Bin0']['signal']
    hists['signal1'].style = styles.lineStyle( ROOT.kBlack, width=2 )
    if inSignalRegions:
        hists['signal2'] = postFitResults2['hists']['shapes_prefit']['Bin0']['signal']

#hists['BSM'].legendOption = 'l'

boxes = []
ratio_boxes = []
if options.combined:
    for ib in range(1, 1 + hists[2016]['total_background'].GetNbinsX() ):
        val = hists[2016]['total_background'].GetBinContent(ib) + hists[2017]['total_background'].GetBinContent(ib) + hists[2018]['total_background'].GetBinContent(ib)
        if val<0: continue
        #sys = math.sqrt((hists[2016]['total_background'].GetBinError(ib) * hists[2016]['total_background'].GetBinError(ib))+ (hists[2017]['total_background'].GetBinError(ib) * hists[2017]['total_background'].GetBinError(ib)) +( hists[2018]['total_background'].GetBinError(ib) * hists[2018]['total_background'].GetBinError(ib)))
        variance = covariance['dc_2016_%s'%(ib-1)]['dc_2016_%s'%(ib-1)] + covariance['dc_2017_%s'%(ib-1)]['dc_2017_%s'%(ib-1)] + covariance['dc_2018_%s'%(ib-1)]['dc_2018_%s'%(ib-1)]
        variance += covariance['dc_2016_%s'%(ib-1)]['dc_2017_%s'%(ib-1)] + covariance['dc_2016_%s'%(ib-1)]['dc_2018_%s'%(ib-1)] + covariance['dc_2017_%s'%(ib-1)]['dc_2018_%s'%(ib-1)]
        sys = math.sqrt(variance)
        if val > 0:
            sys_rel = sys/val
        else:
            sys_rel = 1.
        
        # uncertainty box in main histogram
        box = ROOT.TBox( hists[2016]['total_background'].GetXaxis().GetBinLowEdge(ib),  max([0.006, val-sys]), hists[2016]['total_background'].GetXaxis().GetBinUpEdge(ib), max([0.006, val+sys]) )
        box.SetLineColor(ROOT.kGray+1)
        box.SetFillStyle(3244)
        box.SetFillColor(ROOT.kGray+1)
        
        # uncertainty box in ratio histogram
        r_box = ROOT.TBox( hists[2016]['total_background'].GetXaxis().GetBinLowEdge(ib),  max(0.11, 1-sys_rel), hists[2016]['total_background'].GetXaxis().GetBinUpEdge(ib), min(1.9, 1+sys_rel) )
        r_box.SetLineColor(ROOT.kGray+1)
        r_box.SetFillStyle(3244)
        r_box.SetFillColor(ROOT.kGray+1)

        boxes.append( box )
        hists[2016]['total_background'].SetBinError(ib, 0)
        ratio_boxes.append( r_box )
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
        hists['total_background'].SetBinError(ib, 0)
        ratio_boxes.append( r_box )

def drawObjects( isData=False, lumi=36. ):
    tex = ROOT.TLatex()
    tex.SetNDC()
    tex.SetTextSize(0.05)
    tex.SetTextAlign(11) # align right
    lines = [
      (0.15, 0.945, 'CMS Simulation') if not isData else ( (0.15, 0.945, 'CMS') if not options.preliminary else (0.15, 0.945, 'CMS #bf{#it{Preliminary}}')),
      (0.70, 0.945, '#bf{%s fb^{-1} (13 TeV)}'%lumi )
    ]
    return [tex.DrawLatex(*l) for l in lines]

def setBinLabels( hist ):
    for i in range(1, hist.GetNbinsX()+1):
        hist.GetXaxis().SetBinLabel(i, "   %s"%((i+1)/2 if i%2==1 else ''))

def drawDivisions(regions):
    #print len(regions)
    min = 0.15
    max = 0.95
    diff = (max-min) / len(regions)
    lines = []
    lines2 = []
    line = ROOT.TLine()
#   line.SetLineColor(38)
    line.SetLineWidth(1)
    line.SetLineStyle(2)
    lines = [ (min+4*i*diff,  0.015, min+4*i*diff, 0.85) if min+4*i*diff<0.74 else (min+4*i*diff,  0.013, min+4*i*diff, 0.54) for i in range(1,len(regions)/4 + 1) ]
    lines += [ (min+4*3*diff, 0.85, min+4*3*diff, 0.93 ), (min+4*6*diff, 0.85, min+4*6*diff, 0.93 )]
    return [line.DrawLineNDC(*l) for l in lines] + [tex.DrawLatex(*l) for l in []] + [tex2.DrawLatex(*l) for l in lines2]

def drawLabels( regions ):
    tex = ROOT.TLatex()
    tex.SetNDC()
    tex.SetTextSize(0.028)
    tex.SetTextAngle(0)
    tex.SetTextAlign(12) # align right
    min = 0.15
    max = 0.95
    diff = (max-min) / len(regions)
    lines = [ (0.15 + 0.10, 0.88, "100 < M_{T2}(ll) < 140 GeV"), (0.15 + 0.45, 0.88, "140 < M_{T2}(ll) < 240 GeV") ] 
    #lines =  [(min+(8*i+0.90)*diff, 0.850,  "M_{T2}(ll)=3")   for i, r in enumerate(regions[:-4][::4])]
    return [tex.DrawLatex(*l) for l in lines] if len(regions)>12 else []

def drawLabelsRot( regions ):
    tex = ROOT.TLatex()
    tex.SetNDC()
    tex.SetTextSize(0.028)
    tex.SetTextAngle(90)
    tex.SetTextAlign(12) # align right
    min = 0.15
    max = 0.95
    diff = (max-min) / len(regions)
    lines = [(min+(i*4+2.7)*diff, 0.545 if i<3 else 0.285,  r.texStringForVar('dl_mt2blbl')) for i, r in enumerate(regions[::2]) if i < 6]
    lines += [(min+(24.8)*diff, 0.7, "M_{T2}(ll) > 240 GeV")]
    return [tex.DrawLatex(*l) for l in lines] 

lumiStr = round(lumiStr,1) if not options.combined else int(lumiStr)
drawObjects = drawObjects( isData=isData, lumi=lumiStr) + boxes
if options.region == 'signalOnly':
    drawObjects += drawDivisions( regions ) + drawLabels( regions ) + drawLabelsRot( regions )
if options.combined:
    if inSignalRegions:
        plots = [ bkgHist, [histos['data']], [histos['signal1']], [histos['signal2']]]
    else:
        plots = [ bkgHist, [histos['data']], [histos['signal1']]]
        
    plotName = "%s_COMBINED"%options.region
else:
    plots = [ bkgHists, [hists['data']], [hists['signal1']],[hists['signal2'] ]]
    plotName = "%s_%s"%(options.region, options.year)
if options.postFit:
    plotName += '_postFit'


yMax = 90000. if not options.combined else 900000.

plotting.draw(
    Plot.fromHisto(plotName,
                plots,
                texX = "",
                texY = 'Number of events',
            ),
    plot_directory = os.path.join(plot_directory, "controlRegions", 'v7'),
    logX = False, logY = True, sorting = False, 
    #legend = (0.75,0.80-0.010*32, 0.95, 0.80),
    legend = (0.70,0.55, 0.95, 0.85),
    widths = {'x_width':900, 'y_width':600, 'y_ratio_width':250},
    yRange = (0.2,yMax),
    #yRange = (0.03, [0.001,0.5]),
    ratio = {'yRange': (0.11, 1.89), 'texY':'Data/pred', 'histos':[(1,0)], 'drawObjects':ratio_boxes, #+ drawLabelsLower( regions ) +drawHeadlineLower( regions ) + drawDivisionsLower(regions),
            'histModifications': [lambda h: setBinLabels(h), lambda h: h.GetYaxis().SetTitleSize(32), lambda h: h.GetYaxis().SetLabelSize(28), lambda h: h.GetYaxis().SetTitleOffset(1.2), lambda h: h.GetXaxis().SetTitleSize(32), lambda h: h.GetXaxis().SetLabelSize(27), lambda h: h.GetXaxis().SetLabelOffset(0.035)]} ,
    drawObjects = drawObjects,
    histModifications = [lambda h: h.GetYaxis().SetTitleSize(32), lambda h: h.GetYaxis().SetLabelSize(28), lambda h: h.GetYaxis().SetTitleOffset(1.2)],
    #canvasModifications = [ lambda c : c.SetLeftMargin(0.08), lambda c : c.GetPad(2).SetLeftMargin(0.08), lambda c : c.GetPad(1).SetLeftMargin(0.08), lambda c: c.GetPad(2).SetBottomMargin(0.60), lambda c : c.GetPad(1).SetRightMargin(0.03), lambda c: c.GetPad(2).SetRightMargin(0.03) ],
    copyIndexPHP = True,
)

