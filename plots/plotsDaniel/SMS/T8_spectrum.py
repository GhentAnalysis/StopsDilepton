import ROOT
from RootTools.core.standard    import *
from RootTools.core.helpers     import *
from StopsDilepton.tools.user   import plot_directory

stop1       = ROOT.TLine(0.02, 0.95,0.18, 0.95)
chargino1   = ROOT.TLine(0.02, 0.50,0.18, 0.50)
slepton1    = ROOT.TLine(0.02, 0.48,0.18, 0.48)
lsp1        = ROOT.TLine(0.02, 0.05,0.18, 0.05)

stop2       = ROOT.TLine(0.22, 0.95,0.38, 0.95)
chargino2   = ROOT.TLine(0.22, 0.50,0.38, 0.50)
slepton2    = ROOT.TLine(0.22, 0.275,0.38, 0.275)
lsp2        = ROOT.TLine(0.22, 0.05,0.38, 0.05)

stop3       = ROOT.TLine(0.42, 0.95,0.58, 0.95)
chargino3   = ROOT.TLine(0.42, 0.50,0.58, 0.50)
slepton3    = ROOT.TLine(0.42, 0.07,0.58, 0.07)
lsp3        = ROOT.TLine(0.42, 0.05,0.58, 0.05)

stops       = [stop1, stop2, stop3]
charginos   = [chargino1, chargino2, chargino3]
sleptons    = [slepton1, slepton2, slepton3]
lsps        = [lsp1, lsp2, lsp3]

for h in stops:
    h.SetLineColor(ROOT.kBlack)
    h.SetLineWidth(3)

for h in charginos:
    h.SetLineColor(ROOT.kRed+1)
    h.SetLineWidth(3)

for h in sleptons:
    h.SetLineColor(ROOT.kBlue+1)
    h.SetLineWidth(3)

for h in lsps:
    h.SetLineColor(ROOT.kGreen+1)
    h.SetLineWidth(3)

leg = ROOT.TLegend(0.65,0.60,0.90,0.85)
leg.SetFillColor(ROOT.kWhite)
leg.SetShadowColor(ROOT.kWhite)
leg.SetBorderSize(0) 
leg.SetTextSize(0.045)
leg.AddEntry(stop1,     "top squark", 'l')
leg.AddEntry(chargino1, "chargino", 'l')
leg.AddEntry(slepton1,  "slepton", 'l')
leg.AddEntry(lsp1,      "neutralino", 'l')

dummyHist = ROOT.TH1F("h", "", 10,0,1) 
dummyHist.drawOption = "AH"


def drawObjects():
    tex = ROOT.TLatex()
    tex.SetNDC()
    tex.SetTextSize(0.05)
    tex.SetTextAlign(11) # align right
    lines = [
      (0.16, 0.07, '#bf{x=0.95}' ),
      (0.32, 0.07, '#bf{x=0.50}' ),
      (0.48, 0.07, '#bf{x=0.05}' )

    ]
    return [tex.DrawLatex(*l) for l in lines]


def drawYLabel():
    tex = ROOT.TLatex()
    tex.SetNDC()
    tex.SetTextSize(0.05)
    tex.SetTextAlign(11) # align right
    tex.SetTextAngle(90)
    lines = [
      (0.10, 0.50, '#bf{sparticle mass (a.u.)}' ),

    ]
    return [tex.DrawLatex(*l) for l in lines]

plot_path = plot_directory

plotting.draw(
    Plot.fromHisto(name = 'T8_spectrum', histos = [[ dummyHist ]], texX = "", texY = "Mass"),
    drawObjects = stops + charginos + sleptons + lsps + [leg] + drawObjects() + drawYLabel(),
    plot_directory = plot_path, #ratio = ratio, 
    #drawOptions = "AH",
    logX = False, logY = False, sorting = False,
    yRange = (0,1),
    legend = None ,
)
