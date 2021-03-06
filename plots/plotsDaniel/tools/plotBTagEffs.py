import ROOT, pickle, itertools
import os

#from Workspace.HEPHYPythonTools.helpers import *
from StopsDilepton.tools.helpers import *

ROOT.gROOT.LoadMacro('$CMSSW_BASE/src/StopsDilepton/tools/scripts/tdrstyle.C')
ROOT.setTDRStyle()

def varBinHalfOpen(vb):
  if vb[0] < vb[1] : return '[' + str(vb[0]) + ',' +str(vb[1]) + ')'
  if vb[1]==-1 : return '#geq'+ str(vb[0])
  if vb[0]==vb[1] : return str(vb[0])


def getHistMCTruthEfficiencies(MCEff, histname, etaBin = (0,0.8), hadron='b'):
  nBins = len(MCEff)
  hist = ROOT.TH1F(histname,'MC truth b-tag efficiency',nBins,0,nBins)
  effs = []
  for a in sorted(MCEff):
    effs.append(MCEff[a][etaBin][hadron])
  for b in range(1,nBins+1):
    hist.SetBinContent(b,effs[b-1])
    hist.GetXaxis().SetBinLabel(b,varBinHalfOpen(sorted(MCEff.keys())[b-1]))
  return hist

for year in [2016,2017,2018]:

    can = ROOT.TCanvas('can','can',600,600)
    can.SetBottomMargin(0.22)
    
    bTagEffFile = '$CMSSW_BASE/src/Analysis/Tools/data/btagEfficiencyData/TTLep_pow_%s_2j_2l_DeepB_eta.pkl'%year
    effs = pickle.load(file(os.path.expandvars(bTagEffFile)))
 
    etaBin = (0,2.4) if year == 2016 else (0,2.5)
   
    h_b_1 = getHistMCTruthEfficiencies(effs, 'h_b_1', etaBin = etaBin, hadron='b')
    
    h_c_1 = getHistMCTruthEfficiencies(effs, 'h_c_1', etaBin = etaBin, hadron='c')
    
    h_l_1 = getHistMCTruthEfficiencies(effs, 'h_l_1', etaBin = etaBin, hadron='other')
    
    h_b_1.SetLineColor(ROOT.kAzure+12)
    
    h_c_1.SetLineColor(ROOT.kGreen+3)
    
    h_l_1.SetLineColor(ROOT.kRed+1)
    
    hists = [h_c_1,h_l_1]
    
    for h in hists+[h_b_1]:
      h.SetLineWidth(2)
      h.SetMarkerSize(0)
      h.SetMarkerColor(h.GetLineColor())
    
    h_b_1.GetYaxis().SetTitle('Efficiency')
    h_b_1.GetXaxis().SetTitle('p_{T} [GeV]')
    h_b_1.GetXaxis().SetTitleSize(0.05)
    h_b_1.GetXaxis().SetTitleOffset(2.3)
    
    h_b_1.SetMaximum(1)
    h_b_1.SetMinimum(0.00)
    h_b_1.LabelsOption("v")
    
    
    
    leg = ROOT.TLegend(0.16,0.78,0.44,0.95)
    leg.SetFillColor(ROOT.kWhite)
    leg.SetShadowColor(ROOT.kWhite)
    leg.SetBorderSize(1)
    leg.SetTextSize(0.035)
    
    leg.AddEntry(None,'b hadrons','')
    leg.AddEntry(h_b_1,'#bf{0.0 #leq |#eta| < 2.4}')
    
    leg_2 = ROOT.TLegend(0.44,0.78,0.71,0.95)
    leg_2.SetFillColor(ROOT.kWhite)
    leg_2.SetShadowColor(ROOT.kWhite)
    leg_2.SetBorderSize(1)
    leg_2.SetTextSize(0.035)
    
    leg_2.AddEntry(None,'c hadrons','')
    leg_2.AddEntry(h_c_1,'#bf{0.0 #leq |#eta| < 2.4}')
    
    leg_3 = ROOT.TLegend(0.71,0.78,0.98,0.95)
    leg_3.SetFillColor(ROOT.kWhite)
    leg_3.SetShadowColor(ROOT.kWhite)
    leg_3.SetBorderSize(1)
    leg_3.SetTextSize(0.035)
    
    leg_3.AddEntry(None,'light/gluon','')
    leg_3.AddEntry(h_l_1,'#bf{0.0 #leq |#eta| < 2.4}')
    
    h_b_1.Draw('hist')
    for h in hists:
      h.Draw('hist same')
    
    leg.Draw()
    leg_2.Draw()
    leg_3.Draw()
    
    
    latex1 = ROOT.TLatex()
    latex1.SetNDC()
    latex1.SetTextSize(0.04)
    latex1.SetTextAlign(11)
    
    latex1.DrawLatex(0.16,0.96,'CMS #bf{#it{Simulation}}')
    latex1.DrawLatex(0.85,0.96,'#bf{(13TeV)}')
    
    #can.SetLogy()
     
    can.Print('/afs/hephy.at/user/d/dspitzbart/www/stopsDileptonLegacy/btagEfficiency/TTLep_pow_%s_old.png'%year)
    can.Print('/afs/hephy.at/user/d/dspitzbart/www/stopsDileptonLegacy/btagEfficiency/TTLep_pow_%s_old.pdf'%year)
    can.Print('/afs/hephy.at/user/d/dspitzbart/www/stopsDileptonLegacy/btagEfficiency/TTLep_pow_%s_old.root'%year)
