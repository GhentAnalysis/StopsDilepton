'''
Read a combine data card and create a pandas DataFrame from it for further use.

First, create a python dictionary from the data card.
Bin, process, count, stat,

https://pandas.pydata.org/pandas-docs/stable/user_guide/advanced.html

'''

import pandas

#cardFile = '/afs/hephy.at/data/cms05/StopsDileptonLegacy/results/v8//COMBINED/fitAll/cardFiles/T2tt/observed//T2tt_800_100_combination.txt'



def getDict(cardFile):
    #result = {}
    
    '''
    Get the meta information in the first go, then fill everything
    '''
    systematics = ['yield']
    binList = False
    estimateList = False
    with open(cardFile) as f:
        for line in f:
            if len(line.split())<2: continue
            if line.split()[0] == "bin":
                if not binList:
                    binList = True
                    bins = line.split()[1:]
                else:
                    binList = line.split()[1:]
            if line.split()[1] == 'lnN' and not line.split()[0].startswith('Stat'):
                systematics.append(line.split()[0])
            if line.split()[0] == "process":
                if not estimateList: estimateList = line.split()[1:]
    
    # trick to get unique entries of estimate list
    uniqueEstmates = list(set(estimateList))
    processes   = uniqueEstmates + ['obs']
    systematics += ['stat']
    
    # Build the skelleton
    result = { b: { proc: { sys:0 for sys in systematics} for proc in processes } for b in bins }

    with open(cardFile) as f:
        for line in f:
            if len(line.split())<20: continue
            # fill observation
            if line.split()[0] == 'observation':
                for i,b in enumerate(bins):
                    result[b]['obs']['yield'] = int(line.split()[i+1])
            
            # fill estimates
            if line.split()[0] == 'rate':
                for i,b in enumerate(binList):
                    #print b, uniqueEstmates[(i+1)%len(uniqueEstmates)], line.split()[i+1]
                    #print b, estimateList[i], float(line.split()[i+1])
                    result[b][estimateList[i]]['yield'] = float(line.split()[i+1])
            
            # now put the uncertainties
            if line.split()[1] == 'lnN' and not line.split()[0].startswith('Stat'):
                for i,b in enumerate(binList):
                    try:
                        result[b][estimateList[i]][line.split()[0]] = float(line.split()[i+2]) - 1
                        #print b, line.split()[0], uniqueEstmates[(i)%len(uniqueEstmates)], float(line.split()[i+2]) - 1
                        #print b, line.split()[0], estimateList[i], float(line.split()[i+2]) - 1
                    except:
                        result[b][estimateList[i]][line.split()[0]] = 0

            # stat uncertainties. this needed some level of hardcoding, unfortunately.
            if line.split()[1] == 'lnN' and line.split()[0].startswith('Stat'):
                if len(line.split()[0].split('_'))>3:
                    # examples: dc_2016_Bin0, Stat_Bin0_DY_2016
                    binname = "dc_%s_%s"%(line.split()[0].split('_')[3], line.split()[0].split('_')[1])
                    result[binname][line.split()[0].split('_')[2]]['stat'] = float([ x for x in line.split()[2:] if x != '-' ][0]) - 1

    return result



def getDictFromShape(cardFile):
    return
    '''
    Get the meta information in the first go, then fill everything
    '''
    binList = False
    estimateList = False
    with open(cardFile) as f:
        for line in f:
            if len(line.split())<2: continue
            if line.split()[0] == "bin":
                if not binList:
                    binList = True
                    bins = line.split()[1:]
                else:
                    binList = line.split()[1:]
            if line.split()[1] == 'lnN' or line.split()[1] == 'shape':
                systematics.append(line.split()[0])
            if line.split()[0] == "process":
                if not estimateList: estimateList = line.split()[1:]
    
    # trick to get unique entries of estimate list
    uniqueEstmates = list(set(estimateList))
    processes   = uniqueEstmates + ['obs']
    
    # Build the skelleton
    result = { b: { proc: { sys:0 for sys in systematics} for proc in processes } for b in bins }

    # Not ready yet
    return systematics, processes


'''
first       bin0                bin1          
second     proc1     proc2     proc1     proc2
sys1    0.895717  0.805244 -1.206412  2.565646
sys2    0.410835  0.813850  0.132003 -0.827317
sys3   -1.413681  1.607920  1.024180  0.569605
'''


#arrays = [['bar', 'bar', 'baz', 'baz', 'foo', 'foo', 'qux', 'qux']


