__author__='pbmanis'

"""
RIFunc.py does a rate-intensity function for the selected model
Uses multiprocessing to take advantage of the processors in a given system

"""
import numpy as np
import pyqtgraph.multiprocess as mproc
import endbulb1

RIlevels = np.arange(0, 101, 50)
nLevels = RIlevels.shape[0]
CF0 = 5000.
nWorkers = 8

TASKS = [s for s in range(nLevels)]
res2 = [None]*len(TASKS)
endbulb = [None]*nLevels
print endbulb
with mproc.Parallelize(enumerate(TASKS), results=res2, workers=nWorkers) as tasker:
    for i, x in tasker:
        print i
        endbulb[i] = endbulb1.Endbulb()
        endbulb[i].CF = CF0
        endbulb[i].dBSPL = RIlevels[i]
        res = endbulb[i].run()
        tasker.results[i] = {'r': res, 'pars': endbulb[i].info()}

for i, result in enumerate(res2):
    if result is None:
        continue
    print 'i: %d, result: ' %( i), result['pars']
