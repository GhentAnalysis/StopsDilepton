inFile = file('stops_13TeV_NNLL_approx.txt','r')

lines = inFile.readlines()

for l in lines:
    x = l.split()
    if len(x)>0:
        print "{:5}: ({:9}, {:9}),".format(int(float(x[0])), float(x[2]), float(x[4])/100.)


