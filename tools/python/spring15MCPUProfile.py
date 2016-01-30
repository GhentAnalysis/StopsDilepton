vals = [\
4.8551E-07,
1.74806E-06,
3.30868E-06,
1.62972E-05,
4.95667E-05,
0.000606966,
0.003307249,
0.010340741,
0.022852296,
0.041948781,
0.058609363,
0.067475755,
0.072817826,
0.075931405,
0.076782504,
0.076202319,
0.074502547,
0.072355135,
0.069642102,
0.064920999,
0.05725576,
0.047289348,
0.036528446,
0.026376131,
0.017806872,
0.011249422,
0.006643385,
0.003662904,
0.001899681,
0.00095614,
0.00050028,
0.000297353,
0.000208717,
0.000165856,
0.000139974,
0.000120481,
0.000103826,
8.88868E-05,
7.53323E-05,
6.30863E-05,
5.21356E-05,
4.24754E-05,
3.40876E-05,
2.69282E-05,
2.09267E-05,
1.5989E-05,
4.8551E-06,
2.42755E-06,
4.8551E-07,
2.42755E-07,
1.21378E-07,
4.8551E-08,
]
import ROOT
mcProfile = ROOT.TH1D('pileup','pileup',50,0,50)
for i in range(51):
    mcProfile.Fill(i, vals[i])

mcProfile.Scale(1./mcProfile.Integral())
