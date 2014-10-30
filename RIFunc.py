__author__='pbmanis'

"""
RIFunc.py does a rate-intensity function for the selected model
Uses multiprocessing to take advantage of the processors in a given system

"""
import numpy as np
import pyqtgraph.multiprocess as mproc
import endbulb1
import pyqtgraph as pg
import pylibrary.pyqtgraphPlotHelpers as pgh
import faulthandler

faulthandler.enable()
def plotResults(results):

    win = pgh.figure(title='RIFunc')
    layout = pgh.LayoutMaker(rows=len(results), cols=2, win=win, labelEdges=True, ticks='talbot')

    for i, result in enumerate(results):
        if result is None:
            continue
        layout.plot((i,0), result['r']['t'], result['r']['vm'])
        layout.title((i,0), '%d dBSPL'%result['pars']['SPL'])
        layout.title((i,1), title='AN spikes at %d dBSPL' % result['pars']['SPL'])
        nconverge = result['pars']['nConverge']
        dist = 1./nconverge
        ytick = np.linspace(0, dist, nconverge+1)
        for k in range(nconverge):
            spiketrain = result['r']['preCell%03d'%k]
            nsp = len(spiketrain)
            vt = pg.ScatterPlotItem(spiketrain, [ytick[k+1]]*nsp, symbol='+', pen=(k,nconverge))
            layout.getPlot((i,1)).addItem(vt)
    pgh.show()



RIlevels = np.arange(0, 101, 50)
nLevels = RIlevels.shape[0]
CF0 = 5000.
nWorkers = 8

TASKS = [s for s in range(nLevels)]
results = [None]*len(TASKS)
endbulb = [None]*nLevels

with mproc.Parallelize(enumerate(TASKS), results=results, workers=nWorkers) as tasker:
    for i, x in tasker:
        print i
        endbulb[i] = endbulb1.Endbulb()
        endbulb[i].CF = CF0
        endbulb[i].dBSPL = RIlevels[i]
        endbulb[i].run()
        #tres = {'index': i}
        #print tres
        #tasker.results[i] = {'r': tres, 'pars': {'db': RIlevels[i], 'cf': CF0}}
        tasker.results[i] = {'r': endbulb[i].getResults(), 'pars': endbulb[i].info()}

plotResults(results)


