# encoding: utf-8
"""
pyqtgraphPlotHelpers.py

Routines to help use pyqtgraph and make cleaner plots
as well as get plots read for publication. 

Intially copied from PlotHelpers.py for matplotlib.
Modified to allow us to use a list of axes, and operate on all of those,
or to use just one axis if that's all that is passed.
Therefore, the first argument to these calls can either be a pyqtgraph axis object,
or a list of axes objects.  2/10/2012 pbm.

Created by Paul Manis on 2010-03-09.
"""


import string

stdFont = 'Arial'

import scipy.stats
import numpy as np
try:
    import pyqtgraph as pg
    from PyQt4 import QtCore, QtGui
except:
    pass
import talbotetalTicks as ticks # logical tick formatting... 

"""
Basic functions:
"""


def nice_plot(plotlist, spines=['left', 'bottom'], position=10, direction='inward', axesoff = False):
    """ Adjust a plot so that it looks nicer than the default matplotlib plot
        Also allow quickaccess to things we like to do for publication plots, including:
           using a calbar instead of an axes: calbar = [x0, y0, xs, ys]
           inserting a reference line (grey, 3pt dashed, 0.5pt, at refline = y position)
        
        Paramaetrs
        ----------
        plotlist : list
            a plot handle or list of plot handles to which the "niceplot" will be applied
        spines : list
            a list of which axes should have spines. Not relevant for pyqtgraph
        position : int
            not relevant for pyqtgraph
        direction : string
            need to implement for pyqtgraph
        axesoff : boolean
            flag that forces plots to turn axes off
        
    """
    if type(plotlist) is not list:
        plotlist = [plotlist]
    for pl in plotlist:
        if axesoff is True:
            pl.hideAxis('bottom')
            pl.hideAxis('left')


def noaxes(plotlist, whichaxes = 'xy'):
    """ take away all the axis ticks and the lines
    
    Parameters
    ----------
        plotlist : list
            list of plot handles
        whichaxes : string
            string describing which axes to remove: 'x', 'y', or 'xy' for both
    
    """
    if type(plotlist) is not list:
        plotlist = [plotlist]
    for pl in plotlist:
        if 'x' in whichaxes:
            pl.hideAxis('bottom')
        if 'y' in whichaxes:
            pl.hideAxis('left')


def setY(ax1, ax2):
    """
    Set the Y axis of all the plots in ax2 to be like ax1
    
    Parameters
    ----------
        ax1 : pyqtgraph plot instance
        ax2 : list
            list of target plots that will have the axes properties of ax1
    """
    if type(ax1) is list:
        print 'PlotHelpers: cannot use list as source to set Y axis'
        return
    if type(ax2) is not list:
        ax2 = [ax2]
    y = ax1.getAxis('left')
    refy = y.range # return the current range
    for ax in ax2:
        ax.setRange(refy)


def setX(ax1, ax2):
    """
    Set the X axis of all the plots in ax2 to be like ax1

    Parameters
    ----------
        ax1 : pyqtgraph plot instance
        ax2 : list
            list of target plots that will have the axes properties of ax1
    """
    if type(ax1) is list:
        print 'PlotHelpers: cannot use list as source to set X axis'
        return
    if type(ax2) is not list:
        ax2 = [ax2]
    x = ax1.getAxis('bottom')
    refx = x.range
    for ax in ax2:
        ax.setrange(refx)


def labelPanels(axl, axlist=None, font='Arial', fontsize=18, weight = 'normal'):
    """
    Label the panels like a specific panel
    
    Parameters
    ----------
    axl : dict or list
    axlist : list, optional
         list of labels to use for the axes, defaults to None
    font : str, optional
        Font to use for the labels, defaults to Arial
    fontsize : int, optional
        Font size in points for the labels, defaults to 18
    weight : str, optional
        Font weight to use, defaults to 'normal'
    
    """
    if type(axl) is dict:
        axt = [axl[x] for x in axl]
        axlist = axl.keys()
        axl = axt
    if type(axl) is not list:
        axl = [axl]
    if axlist is None:
        axlist = string.uppercase(1,len(axl)) # assume we wish to go in sequence

    for i, ax in enumerate(axl):
        labelText = pg.TextItem(axlist[i])
        y = ax.getAxis('left').range
        x = ax.getAxis('bottom').range
        ax.addItem(labelText)
        labelText.setPos(x[0], y[1])


def listAxes(axd):
    """
    make a list of the axes from the dictionary of axes
    
    Parameters
    ----------
    axd : dict
        a dict of axes, whose values are returned in a list
    
    Returns
    -------
    list : a list of the axes
    
    """
    if type(axd) is not dict:
        if type(axd) is list:
            return axd
        else:
            print 'listAxes expects dictionary or list; type not known (fix the code)'
            raise
    axl = [axd[x] for x in axd]
    return axl


def cleanAxes(axl):
    """
    """
    if type(axl) is not list:
        axl = [axl]
    # does nothing at the moment, as axes are already "clean"
    # for ax in axl:
    #
    #    update_font(ax)


def formatTicks(axl, axis='xy', fmt='%d', font='Arial'):
    """
    Convert tick labels to intergers
    to do just one axis, set axis = 'x' or 'y'
    control the format with the formatting string
    """
    if type(axl) is not list:
        axl = [axl]

def autoFormatTicks(axl, axis='xy', font='Arial'):
    if type(axl) is not list:
        axl = [axl]
    for ax in axl:
        if 'x' in axis:
            b = ax.getAxis('bottom')
            x0 = b.range
 #           setFormatter(ax,  x0, x1, axis = 'x')
        if 'y' in axis:
            l = ax.getAxis('left')
            y0= l.range
 #           setFormatter(ax, y0, y1, axis = 'y')


def setFormatter(ax, x0, x1, axis='x'):
    datarange = np.abs(x0-x1)
    mdata = np.ceil(np.log10(datarange))
    # if mdata > 0 and mdata <= 4:
    #     majorFormatter = FormatStrFormatter('%d')
    # elif mdata > 4:
    #     majorFormatter = FormatStrFormatter('%e')
    # elif mdata <= 0 and mdata > -1:
    #     majorFormatter = FormatStrFormatter('%5.1f')
    # elif mdata < -1 and mdata > -3:
    #     majorFormatatter = FormatStrFormatter('%6.3f')
    # else:
    #     majorFormatter = FormatStrFormatter('%e')
    # if axis == 'x':
    #     ax.xaxis.set_major_formatter(majorFormatter)
    # else:
    #     ax.yaxis.set_major_formatter(majorFormatter)


def update_font(axl, size=6, font=stdFont):
    pass
    # if type(axl) is not list:
    #     axl = [axl]
    # fontProperties = {'family':'sans-serif','sans-serif':[font],
    #         'weight' : 'normal', 'size' : size}
    # for ax in axl:
    #     for tick in ax.xaxis.get_major_ticks():
    #           tick.label1.set_family('sans-serif')
    #           tick.label1.set_fontname(stdFont)
    #           tick.label1.set_size(size)
    #
    #     for tick in ax.yaxis.get_major_ticks():
    #           tick.label1.set_family('sans-serif')
    #           tick.label1.set_fontname(stdFont)
    #           tick.label1.set_size(size)
    #     ax.set_xticklabels(ax.get_xticks(), fontProperties)
    #     ax.set_yticklabels(ax.get_yticks(), fontProperties)
    #     ax.xaxis.set_smart_bounds(True)
    #     ax.yaxis.set_smart_bounds(True)
    #     ax.tick_params(axis = 'both', labelsize = 9)


def lockPlot(axl, lims, ticks=None):
    """ 
        This routine forces the plot of invisible data to force the axes to take certain
        limits and to force the tick marks to appear. 
        call with the axis and lims = [x0, x1, y0, y1]
    """ 
    if type(axl) is not list:
        axl = [axl]
    plist = []
    for ax in axl:
        y = ax.getAxis('left')
        x = ax.getAxis('bottom')
        x.setRange(lims[0], lims[1])
        y.setRange(lims[2], lims[3])


def adjust_spines(axl, spines = ('left', 'bottom'), direction = 'outward', distance=5, smart=True):
    pass
    # if type(axl) is not list:
    #     axl = [axl]
    # for ax in axl:
    #     # turn off ticks where there is no spine
    #     if 'left' in spines:
    #         ax.yaxis.set_ticks_position('left')
    #     else:
    #         # no yaxis ticks
    #         ax.yaxis.set_ticks([])
    #
    #     if 'bottom' in spines:
    #         ax.xaxis.set_ticks_position('bottom')
    #     else:
    #         # no xaxis ticks
    #         ax.xaxis.set_ticks([])
    #     for loc, spine in ax.spines.iteritems():
    #         if loc in spines:
    #             spine.set_position((direction,distance)) # outward by 10 points
    #             if smart is True:
    #                 spine.set_smart_bounds(True)
    #             else:
    #                 spine.set_smart_bounds(False)
    #         else:
    #             spine.set_color('none') # don't draw spine
    #     return
    #
def calbar(plotlist, calbar=None, axesoff=True, orient='left', unitNames=None):
    """ draw a calibration bar and label it up. The calibration bar is defined as:
        [x0, y0, xlen, ylen]
        
        Parameters
        ----------
        plotlist : list
            a plot item or a list of plot items for which a calbar will be applied
        calbar : list, optional
            a list with 4 elements, describing the calibration bar
                [xposition, yposition, xlength, ylength] in units of the data inthe plot
                defaults to None
        axesoff : boolean, optional
                Set true to turn off the standard axes, defaults to True
        orient : text, optional
                'left': put vertical part of the bar on the left
                'right': put the vertical part of the bar on the right
                defaults to 'left'
        unitNames: str, optional
                a dictionary with the names of the units to append to the calibration bar
                lengths. Example: {'x': 'ms', 'y': 'nA'}
                defaults to None
        Returns
        -------
        Nothing
    """
    
    if type(plotlist) is not list:
        plotlist = [plotlist]
    for pl in plotlist:
        if axesoff is True:
            noaxes(pl)
        Vfmt = '%.0f'
        if calbar[2] < 1.0:
            Vfmt = '%.1f'
        Hfmt = '%.0f'
        if calbar[3] < 1.0:
            Hfmt = '%.1f'
        if unitNames is not None:
            Vfmt = Vfmt + ' ' + unitNames['x']
            Hfmt = Hfmt + ' ' + unitNames['y']
        Vtxt = pg.TextItem(Vfmt % calbar[2], anchor=(0.5, 0.5), color=pg.mkColor('k'))
        Htxt = pg.TextItem(Hfmt % calbar[3], anchor=(0.5, 0.5), color=pg.mkColor('k'))
       # print pl
        if calbar is not None:
            if orient == 'left': # vertical part is on the left
                pl.plot([calbar[0], calbar[0], calbar[0]+calbar[2]],
                    [calbar[1]+calbar[3], calbar[1], calbar[1]],
                    pen=pg.mkPen('k'), linestyle = '-', linewidth = 1.5)
                ht = Htxt.setPos(calbar[0]+0.05*calbar[2], calbar[1]+0.5*calbar[3])
            elif orient == 'right': # vertical part goes on the right
                pl.plot([calbar[0] + calbar[2], calbar[0]+calbar[2], calbar[0]],
                    [calbar[1]+calbar[3], calbar[1], calbar[1]],
                    pen=pg.mkPen('k'), linestyle = '-', linewidth = 1.5)
                ht = Htxt.setPos(calbar[0]+calbar[2]-0.05*calbar[2], calbar[1]+0.5*calbar[3])
            else:
                print "PlotHelpers.py: I did not understand orientation: %s" % (orient)
                print "plotting as if set to left... "
                pl.plot([calbar[0], calbar[0], calbar[0]+calbar[2]],
                    [calbar[1]+calbar[3], calbar[1], calbar[1]],
                    pen=pg.mkPen('k'), linestyle = '-', linewidth = 1.5)
                ht = Htxt.setPos(calbar[0]+0.05*calbar[2], calbar[1]+0.5*calbar[3])
                Htxt.setText(Hfmt % calbar[3])
            xc = float(calbar[0]+calbar[2]*0.5)  # always centered, below the line
            yc = float(calbar[1]-0.1*calbar[3])
            vt = Vtxt.setPos(xc, yc)
            Vtxt.setText(Vfmt % calbar[2])
            pl.addItem(Htxt)
            pl.addItem(Vtxt)


def refline(axl, refline=None, color=[64, 64, 64], linestyle='--' ,linewidth=0.5, orient='horizontal'):
    """
    Draw a reference line at a particular level of the data on the y axis
    
    Parameters
    ----------
    axl : list
        axis handle or list of axis handles
    refline : float, optional
            the position of the reference line, defaults to None
    color : list, optional
            the RGB color list for the line, in format [r,g,b], defaults to [64, 64, 64] (faint grey line)
    linestyle : str, optional
            defines the linestyle to be used: -- for dash, . for doted, - for solid, -. for dash-dot,
            -.. for -.., etc.
            defaults to '--' (dashed)
    linewidth : float, optional
            width of the line, defaults to 0.5
    """
    if type(axl) is not list:
        axl = [axl]
    if linestyle == '--':
        style = QtCore.Qt.DashLine
    elif linestyle == '.':
        style=QtCore.Qt.DotLine
    elif linestyle == '-':
        style=QtCore.Qt.SolidLine
    elif linestyle == '-.':
        style = QtCore.Qt.DsahDotLine
    elif linestyle == '-..':
        style = QtCore.Qt.DashDotDotLine
    else:
        style = QtCore.Qt.SolidLine # default is solid
    if orient is 'horizontal':
        for ax in axl:
            if refline is not None:
                x = ax.getAxis('bottom')
                xlims = x.range
                ax.plot(xlims, [refline, refline], pen=pg.mkPen(color, width=linewidth, style=style))
    if orient is 'vertical':
        for ax in axl:
            if refline is not None:
                y = ax.getAxis('left')
                ylims = y.range
                ax.plot([refline, refline], [ylims[0]+0.5, ylims[1]-0.5], pen=pg.mkPen(color, width=linewidth, style=style))


def tickStrings(values, scale=1, spacing=None, tickPlacesAdd=1):
    """Return the strings that should be placed next to ticks. This method is called 
    when redrawing the axis and is a good method to override in subclasses.
    
    Parameters
    ----------
    values : array or list
         An array or list of tick values
    scale : float, optional
        a scaling factor (see below), defaults to 1
    spacing : float, optional
        spaceing between ticks (this is required since, in some instances, there may be only 
    one tick and thus no other way to determine the tick spacing). Defaults to None
    tickPlacesToAdd : int, optional
        the number of decimal places to add to the ticks, default is 1
    
    Returns
    -------
    list : a list containing the tick strings
    
    The scale argument is used when the axis label is displaying units which may have an SI scaling prefix.
    When determining the text to display, use value*scale to correctly account for this prefix.
    For example, if the axis label's units are set to 'V', then a tick value of 0.001 might
    be accompanied by a scale value of 1000. This indicates that the label is displaying 'mV', and 
    thus the tick should display 0.001 * 1000 = 1.
    Copied rom pyqtgraph; we needed it here.
    """
    if spacing is None:
        spacing = np.mean(np.diff(values))
    places = max(0, np.ceil(-np.log10(spacing*scale))) + tickPlacesAdd
    strings = []
    for v in values:
        vs = v * scale
        if abs(vs) < .001 or abs(vs) >= 10000:
            vstr = "%g" % vs
        else:
            vstr = ("%%0.%df" % places) % vs
        strings.append(vstr)
    return strings


def crossAxes(axl, xyzero=[0., 0.], limits=[None, None, None, None], **kwds):
    """
    Make the plot(s) have crossed axes at the data points set by xyzero, and optionally
    set axes limits
    
    Parameters
    ----------
    axl : pyqtgraph plot/axes instance or list
        the plot to modify
    xyzero : list
        A 2-element list for the placement of x=0 and y=0, defaults to [0., 0.]
    limits : list
        A 4-element list with the min and max limits of the axes, defaults to all None
    **kwds : keyword arguments to pass to make_crossedAxes
    
    """
    if type(axl) is not list:
        axl = [axl]
    for ax in axl:
        make_crossedAxes(ax, xyzero, limits, **kwds)

 
def make_crossedAxes(ax, xyzero=[0., 0.], limits=[None, None, None, None], ndec=3,
            density=(1.0, 1.0), tickl = 0.0125, insideMargin=0.05, pointSize=12, tickPlacesAdd=(0,0)):
    """
    Parameters
    ----------
    axl : pyqtgraph plot/axes instance or list
        the plot to modify
    xyzero : list
        A 2-element list for the placement of x=0 and y=0, defaults to [0., 0.]
    limits : list
        A 4-element list with the min and max limits of the axes, defaults to all None
    ndec : int
        Number of decimals (would be passed to talbotTicks if that was being called)
    density : tuple
        tick density (for talbotTicks), defaults to (1.0, 1.0)
    tickl : float
        Tick length, defaults to 0.0125
    insideMargin : float
        Inside margin space for plot, defaults to 0.05 (5%)
    pointSize : int
         point size for tick text, defaults to 12
    tickPlacesAdd : tuple
        number of decimal places to add in tickstrings for the ticks, pair for x and y axes, defaults to (0,0)
    
    Returns
    -------
    Nothing
    
    
    """
    # get axis limits
    aleft = ax.getAxis('left')
    abottom = ax.getAxis('bottom')
    aleft.setPos(pg.Point(3., 0.))
    yRange = aleft.range
    xRange =  abottom.range
    hl = pg.InfiniteLine(pos=xyzero[0], angle=90, pen=pg.mkPen('k'))
    ax.addItem(hl)
    vl = pg.InfiniteLine(pos=xyzero[1], angle=0, pen=pg.mkPen('k'))
    ax.addItem(vl)
    ax.hideAxis('bottom')
    ax.hideAxis('left')
    # now create substitue tick marks and labels, using Talbot et al algorithm
    xr = np.diff(xRange)[0]
    yr = np.diff(yRange)[0]
    xmin, xmax = (np.min(xRange) - xr * insideMargin, np.max(xRange) + xr * insideMargin)
    ymin, ymax = (np.min(yRange) - yr * insideMargin, np.max(yRange) + yr * insideMargin)
    xtick = ticks.Extended(density=density[0], figure=None, range=(xmin, xmax), axis='x')
    ytick = ticks.Extended(density=density[1], figure=None, range=(ymin, ymax), axis='y')
    xt = xtick()
    yt = ytick()
    ytk = yr*tickl
    xtk = xr*tickl
    y0 = xyzero[1]
    x0 = xyzero[0]
    tsx = tickStrings(xt, tickPlacesAdd=tickPlacesAdd[0])
    tsy = tickStrings(yt, tickPlacesAdd=tickPlacesAdd[1])
    for i, x in enumerate(xt):
        t = pg.PlotDataItem(x=x*np.ones(2), y=[y0-ytk, y0+ytk], pen=pg.mkPen('k'))
        ax.addItem(t)  # tick mark
        # put text in only if it does not overlap the opposite line
        if x == y0:
            continue
        txt = pg.TextItem(tsx[i], anchor=(0.5, 0), color=pg.mkColor('k')) #, size='10pt')
        txt.setFont(pg.QtGui.QFont('Arial', pointSize=pointSize))
        txt.setPos(pg.Point(x, y0-ytk))
        ax.addItem(txt) #, pos=pg.Point(x, y0-ytk))
    for i, y in enumerate(yt):
        t = pg.PlotDataItem(x=np.array([x0-xtk, x0+xtk]), y=np.ones(2)*y, pen=pg.mkPen('k'))
        ax.addItem(t)
        if y == x0:
            continue
        txt = pg.TextItem(tsy[i], anchor=(1, 0.5), color=pg.mkColor('k')) # , size='10pt')
        txt.setFont(pg.QtGui.QFont('Arial', pointSize=pointSize))
        txt.setPos(pg.Point(x0-xtk, y))
        ax.addItem(txt) #, pos=pg.Point(x, y0-ytk))


class polarPlot():
    """
    Create a polar plot, as a PlotItem for pyqtgraph.
    
    
    """
    def __init__(self, plot=None):
        """
        Instantiate a plot as a polar plot
        
        Parmeters
        ---------
        plot : pyqtgraph plotItem
            the plot that will be converted to a polar plot, defaults to None
            if None, then a new PlotItem will be created, accessible
            as polarPlot.plotItem
        """
        if plot is None:
            self.plotItem = pg.PlotItem()  # create a plot item for the plot
        else:
            self.plotItem = plot
        self.plotItem.setAspectLocked()
        self.plotItem.hideAxis('bottom')
        self.plotItem.hideAxis('left')
        self.gridSet = False
        self.data = None
        self.rMax = None

    def setAxes(self, steps=4, rMax=None, makeGrid=True):
        """
        Make the polar plot axes
        
        Parameters
        ----------
        steps : int, optional
            The number of radial steps for the grid, defaults to 4
        rMax : float, optional
            The maximum radius of the plot, defaults to None (the rMax is 1)
        makeGrid : boolean, optional
            Whether the grid will actually be plotted or not, defaults to True
            
        """
        if makeGrid is False or self.gridSet:
            return
        if rMax is None:
            if self.data is None:
                rMax = 1.0
            else:
                rMax = np.max(self.data['y'])
        self.rMax = rMax
        # Add radial grid lines (theta markers)
        gridPen = pg.mkPen(width=0.55, color='k',  style=QtCore.Qt.DotLine)
        ringPen = pg.mkPen(width=0.75, color='k',  style=QtCore.Qt.SolidLine)  
        for th in np.linspace(0., np.pi*2, 8, endpoint=False):
            rx = np.cos(th)*rMax
            ry = np.sin(th)*rMax
            self.plotItem.plot(x=[0, rx], y=[0., ry], pen=gridPen)
            ang = th*360./(np.pi*2)
            # anchor is odd: 0,0 is upper left corner, 1,1 is lower right corner
            if ang < 90.:
                x=0.
                y=0.5
            elif ang == 90.:
                x=0.5
                y=1
            elif ang < 180:
                x=1.0
                y=0.5
            elif ang == 180.:
                x=1
                y=0.5
            elif ang < 270:
                x=1
                y=0
            elif ang== 270.:
                x=0.5
                y=0
            elif ang < 360:
                x=0
                y=0
            ti = pg.TextItem("%d" % (int(ang)), color=pg.mkColor('k'), anchor=(x,y))
            self.plotItem.addItem(ti)
            ti.setPos(rx, ry)
        # add polar grid lines (r)
        for gr in np.linspace(rMax/steps, rMax, steps):
            circle = pg.QtGui.QGraphicsEllipseItem(-gr, -gr, gr*2, gr*2)
            if gr < rMax:
                circle.setPen(gridPen)
            else:
                circle.setPen(ringPen)
            self.plotItem.addItem(circle)
            ti = pg.TextItem("%d" % (int(gr)), color=pg.mkColor('k'), anchor=(1, 1))
            ti.setPos(gr, 0.)
            self.plotItem.addItem(ti)
        self.gridSet = True

    
    def plot(self, r, theta, vectors=False, arrowhead=True, normalize=False, sort=False, **kwds):
        """
        plot puts the data into a polar plot.
        the plot will be converted to a polar graph
        
        Parameters
        ----------
        r : list or numpy array
            a list or array of radii
        theta : list or numpy array
            a list or array of angles (in radians) corresponding to the values in r
        vectors : boolean, optional
            vectors True means that plot is composed of vectors to each point radiating from the origin, defaults to False
        arrowhead : boolean, optional
            arrowhead True plots arrowheads at the end of the vectors, defaults to True
        normalize : boolean, optional
            normalize forces the plot to be scaled to the max values in r, defaults to False
        sort : boolean, optional
            causes data r, theta to be sorted by theta, defaults to False
        **kwds are passed to the data plot call.

        """

        # sort r, theta by r

        rs = np.array(r)
        thetas = np.array(theta)
        
        if sort:
            indx = np.argsort(thetas)
            theta = thetas
            if not isinstance(indx, np.int64):
                for i, j in enumerate(indx):
                    rs[i] = r[j]
                    thetas[i] = theta[j]


        # Transform to cartesian and plot
        if normalize:
            rs = rs/np.max(rs)
        x = rs * np.cos(thetas)
        y = rs * np.sin(thetas)
        try:
            len(x)
        except:
            x = [x]
            y = [y]
        if vectors:  # plot r,theta as lines from origin
            for i, xi in enumerate(x):
    #            print x[i], y[i]
                if arrowhead:
                    arrowAngle = -(thetas[i]*360/(2*np.pi)+180) # convert to degrees, and correct orientation
                    arrow = pg.ArrowItem(angle=arrowAngle, tailLen=0, tailWidth=1.5, **kwds)
                    arrow.setPos(x[i], y[i])
                    self.plotItem.addItem(arrow)
                self.plotItem.plot([0., x[i]], [0., y[i]], **kwds)
                    
        else:
            self.plotItem.plot(x, y, **kwds)
        self.rMax = np.max(y)
        self.data = {'x': x, 'y': y}

    def hist(self, r, theta, binwidth=np.pi/6., normalize=False, density=False, mode='straight', **kwds):
        """
        plot puts the data into a polar plot as a histogram of the number of observations
        within a wedge
        the plot will be converted to a polar graph
        
        Parameters
        ----------
        r : list or numpy array
            a list or array of radii
        theta : list or numpy array
            a list or array of angles (in radians) corresponding to the values in r
        binwidth : bin width, in radians optional
            vectors True means that plot is composed of vectors to each point radiating from the origin, defaults to 30 degrees (np.pi/6)
        normalize : boolean, optional
            normalize forces the plot to be scaled to the max values in r, defaults to False
        density : boolean, optional
            plot a count histogram, or a density histogram weighted by r values, defaults to False
        mode : str, optional
            'straight' selects straight line between bars. 'arc' makes the end of the bar an arc (truer representation), defaults to 'straight'
        **kwds are passed to the data plot call.

        Returns
        -------
        tuple : (list of rHist, list of bins)
            The histogram that was plotted (use for statistical comparisions)
        """

        rs = np.array(r)
        thetas = np.array(theta)
        twopi = np.pi*2.0
        for i, t in enumerate(thetas):  # restrict to positive half plane [0....2*pi]
            while t < 0.0:
                t += twopi
            while t > twopi:
                t -= twopi
            thetas[i] = t
        bins = np.arange(0, np.pi*2+1e-12, binwidth)
        # compute histogram
        (rhist, rbins) = np.histogram(thetas, bins=bins, weights=rs, density=density)
        # Transform to cartesian and plot
        if normalize:
            rhist = rhist/np.max(rhist)
        xo = rhist * np.cos(bins[:-1])  # get cartesian form
        xp = rhist * np.cos(bins[:-1] + binwidth)
        yo = rhist * np.sin(bins[:-1])
        yp = rhist * np.sin(bins[:-1] + binwidth)
        arcinc = np.pi/100.  # arc increments
        for i in range(len(xp)):
            if mode is 'arc':
                self.plotItem.plot([xo[i], 0., xp[i]], [yo[i], 0., yp[i]], **kwds) # "v" segement
                arcseg = np.arange(bins[i], bins[i+1], arcinc)
                x = np.array(rhist[i] * np.cos(arcseg))
                y = np.array(rhist[i] * np.sin(arcseg))
                self.plotItem.plot(x, y, **kwds)
                
            else:
                self.plotItem.plot([0., xo[i], xp[i], 0.], [0., yo[i], yp[i], 0.], **kwds)
        self.data = {'x': xo, 'y': yo}
        self.rMax = np.max(yo)
        return (rhist, rbins)


    def circmean(self, alpha, axis=None):
        """
        Compute the circular mean of a set of angles along the axis
        
        Parameters
        ----------
        alpha : numpy array
            the angles to compute the circular mean of
        axis : int
            The axis of alpha for the computatoin, defaults to None
        
        Returns
        -------
        float : the mean angle
        """
        mean_angle = np.arctan2(np.mean(np.sin(alpha),axis),np.mean(np.cos(alpha),axis))
        return mean_angle



def talbotTicks(axl, **kwds):
    """
    Adjust the tick marks using the talbot et al algorithm, on an existing plot.
    """
    if type(axl) is not list:
        axl = [axl]
    for ax in axl:
        do_talbotTicks(ax, **kwds)

def do_talbotTicks(ax, ndec=3,
                   density=(1.0, 1.0), insideMargin=0.05, pointSize=None, tickPlacesAdd=(0,0)):
    """
    Change the axis ticks to use the talbot algorithm for ONE axis
    Paramerters control the ticks
    
    Parameters
    ----------
    ax : pyqtgraph axis instance
        the axis to change the ticks on
    ndec : int
        Number of decimals (would be passed to talbotTicks if that was being called)
    density : tuple
        tick density (for talbotTicks), defaults to (1.0, 1.0)
    insideMargin : float
        Inside margin space for plot, defaults to 0.05 (5%)
    pointSize : int
         point size for tick text, defaults to 12
    tickPlacesAdd : tuple
        number of decimal places to add in tickstrings for the ticks, pair for x and y axes, defaults to (0,0)
    
    """
    # get axis limits
    aleft = ax.getAxis('left')
    abottom = ax.getAxis('bottom')
    yRange = aleft.range
    xRange =  abottom.range
    # now create substitue tick marks and labels, using Talbot et al algorithm
    xr = np.diff(xRange)[0]
    yr = np.diff(yRange)[0]
    xmin, xmax = (np.min(xRange) - xr * insideMargin, np.max(xRange) + xr * insideMargin)
    ymin, ymax = (np.min(yRange) - yr * insideMargin, np.max(yRange) + yr * insideMargin)
    xtick = ticks.Extended(density=density[0], figure=None, range=(xmin, xmax), axis='x')
    ytick = ticks.Extended(density=density[1], figure=None, range=(ymin, ymax), axis='y')
    xt = xtick()
    yt = ytick()
    xts = tickStrings(xt, scale=1, spacing=None, tickPlacesAdd = tickPlacesAdd[0])
    yts = tickStrings(yt, scale=1, spacing=None, tickPlacesAdd = tickPlacesAdd[1])
    xtickl = [[(x, xts[i]) for i, x in enumerate(xt)] , []]  # no minor ticks here
    ytickl = [[(y, yts[i]) for i, y in enumerate(yt)] , []]  # no minor ticks here

    #ticks format: [ (majorTickValue1, majorTickString1), (majorTickValue2, majorTickString2), ... ],
    aleft.setTicks(ytickl)
    abottom.setTicks(xtickl)
    # now set the point size (this may affect spacing from axis, and that would have to be adjusted - see the pyqtgraph google groups)
    if pointSize is not None:
        b=QtGui.QFont()
        b.setPixelSize(pointSize)
        aleft.tickFont = b
        abottom.tickFont = b

def violinPlotScatter(ax, data, symbolColor='k', symbolSize=4, symbol='o'):
    """
    Plot data as violin plot with scatter and error bar
    
    Parameters
    ----------
    ax : pyqtgraph plot instance
        is the axs to plot into
    data : dict
        dictionary containing {pos1: data1, pos2: data2}, where pos is the x position for the data in data. Each data
        set iis plotted as a separate column
    symcolor : string, optional
        color of the symbols, defaults to 'k' (black)
    symbolSize : int, optional
        Size of the symbols in the scatter plot, points, defaults to 4
    symbol : string, optoinal
        The symbol to use, defaults to 'o' (circle)
    """

    y = []
    x = []
    xb=np.arange(0,len(data.keys()), 1)
    ybm = [0]*len(data.keys()) # np.zeros(len(sdat.keys()))
    ybs = [0]*len(data.keys()) # np.zeros(len(sdat.keys()))
    for i, k in enumerate(data.keys()):
        yvals = np.array(data[k])
        xvals = pg.pseudoScatter(yvals, spacing=0.4, bidir=True) * 0.2
        ax.plot(x=xvals+i, y=yvals, pen=None, symbol=symbol, symbolSize=symbolSize,
            symbolBrush=pg.mkBrush(symbolColor))
        y.append(yvals)
        x.append([i]*len(yvals))
        ybm[i] = np.nanmean(yvals)
        ybs[i] = np.nanstd(yvals)
        mbar = pg.PlotDataItem(x=np.array([xb[i]-0.2, xb[i]+0.2]), y=np.array([ybm[i], ybm[i]]), pen={'color':'k', 'width':0.75})
        ax.addItem(mbar)
        bar = pg.ErrorBarItem(x=xb, y=np.array(ybm), height=np.array(ybs), beam=0.2, pen={'color':'k', 'width':0.75})
    violin_plot(ax, y, xb, bp=False)
    ax.addItem(bar)
    ticks = [[(v, k) for v, k in enumerate(data.keys())], []]
    ax.getAxis('bottom').setTicks(ticks)


def violin_plot(ax, data, pos, dist=.0, bp=False):
    '''
    create violin plots on an axis
    '''

    if data is None or len(data) == 0:
        return  # skip trying to do the plot

    dist = max(pos)-min(pos)
    w = min(0.15*max(dist,1.0),0.5)
    for i, d in enumerate(data):
        if d == [] or len(d) == 0:
            continue
        k = scipy.stats.gaussian_kde(d) #calculates the kernel density
        m = k.dataset.min() #lower bound of violin
        M = k.dataset.max() #upper bound of violin
        y = np.arange(m, M, (M-m)/100.) # support for violin
        v = k.evaluate(y) #violin profile (density curve)
        v = v / v.max() * w #scaling the violin to the available space
        c1 = pg.PlotDataItem(y=y, x=pos[i]+v, pen=pg.mkPen('k', width=0.5))
        c2 = pg.PlotDataItem(y=y, x=pos[i]-v, pen=pg.mkPen('k', width=0.5))
        #mean = k.dataset.mean()
        #vm = k.evaluate(mean)
        #vm = vm * w
        #ax.plot(x=np.array([pos[i]-vm[0], pos[i]+vm[0]]), y=np.array([mean, mean]), pen=pg.mkPen('k', width=1.0))
        ax.addItem(c1)
        ax.addItem(c2)
        #ax.addItem(hbar)
        f = pg.FillBetweenItem(curve1=c1, curve2=c2, brush=pg.mkBrush((255, 255, 0, 96)))
        ax.addItem(f)

    if bp:
       pass
       # bpf = ax.boxplot(data, notch=0, positions=pos, vert=1)
       # pylab.setp(bpf['boxes'], color='black')
       # pylab.setp(bpf['whiskers'], color='black', linestyle='-')


def labelAxes(plot, xtext, ytext, **kwargs):
    """
    helper to label up the plot
    
    Parameters
    -----------
    plot : plot item
    xtext : string
        text for x axis
    ytext : string
        text for y axis
    **kwargs : keywords
        additional arguments to pass to pyqtgraph setLabel
    """

    plot.setLabel('bottom', xtext, **kwargs)
    plot.setLabel('left', ytext, **kwargs)


def labelPanels(plot, label=None, **kwargs):
    """
        helper to label up the plot
        Inputs: plot item
                text for x axis
                text for yaxis
                plot title (on top) OR
                plot panel label (for example, "A", "A1")
    """

    if label is not None:
        setPlotLabel(plot, plotlabel="%s" % label, **kwargs)
    else:
        setPlotLabel(plot, plotlabel="")


def labelTitles(plot, title=None, **kwargs):
    """
    Set the title of a plotitem. Basic HTML formatting is allowed, along
    with "size", "bold", "italic", etc..
    If the title is not defined, then a blank label is used
    A title is a text label that appears centered above the plot, in 
    QGridLayout (position 0,2) of the plotitem.
    params
    -------
    :param plotitem: The plot item to label
    :param title: The text string to use for the label
    :kwargs: keywords to pass to the pg.LabelItem
    :return: None
    
    """

    if title is not None:
        plot.setTitle(title="<b><large>%s</large></b>" % title, visible=True, **kwargs)    
    else:  # clear the plot title
        plot.setTitle(title=" ")


def setPlotLabel(plotitem, plotlabel='', **kwargs):
    """
    Set the plotlabel of a plotitem. Basic HTML formatting is allowed, along
    with "size", "bold", "italic", etc..
    If plotlabel is not defined, then a blank label is used
    A plotlabel is a text label that appears the upper left corner of the
    QGridLayout (position 0,0) of the plotitem.
    params
    -------
    :param plotitem: The plot item to label
    :param plotlabel: The text string to use for the label
    :kwargs: keywords to pass to the pg.LabelItem
    :return: None
    
    """
    
    plotitem.LabelItem = pg.LabelItem(plotlabel, **kwargs)
    plotitem.LabelItem.setMaximumHeight(30)
    plotitem.layout.setRowFixedHeight(0, 30)
    plotitem.layout.addItem(plotitem.LabelItem, 0, 0)
    plotitem.LabelItem.setVisible(True)


class LayoutMaker():
    def __init__(self, win=None, cols=1, rows=1, letters=True, titles=False, labelEdges=True, margins=4, spacing=4, ticks='default'):
        self.sequential_letters = string.ascii_uppercase
        self.cols = cols
        self.rows = rows
        self.letters = letters
        self.titles = titles
        self.edges = labelEdges
        self.margins = margins
        self.spacing = spacing
        self.rcmap = [None]*cols*rows
        self.plots = None
        self.win = win
        self.ticks = ticks
        self._makeLayout(letters=letters, titles=titles, margins=margins, spacing=spacing)
        #self.addLayout(win)

    # def addLayout(self, win=None):
    #     if win is not None:
    #         win.setLayout(self.gridLayout)

    def getCols(self):
        return self.cols

    def getRows(self):
        return self.rows

    def mapFromIndex(self, index):
        """
        for a given index, return the row, col tuple associated with the index
        """
        return self.rcmap[index]

    def getPlot(self, index):
        """
        return the plot item in the list corresponding to the index n
        """
        if isinstance(index, tuple):
            r, c = index
        elif isinstance(index, int):
            r, c = self.rcmap[index]
        else:
            raise ValueError ('pyqtgraphPlotHelpers, LayoutMaker plot: index must be int or tuple(r,c)')
        return self.plots[r][c]

    def plot(self, index, *args, **kwargs):
        p = self.getPlot(index).plot(*args, **kwargs)
        if self.ticks == 'talbot':
            talbotTicks(self.getPlot(index))
        return p

    def _makeLayout(self, letters=True, titles=True, margins=4, spacing=4):
        """
        Create a multipanel plot.
        The pyptgraph elements (widget, gridlayout, plots) are stored as class variables.
        The layout is always a rectangular grid with shape (cols, rows)
        if letters is true, then the plot is labeled "A, B, C..." Indices move horizontally first, then vertically
        margins sets the margins around the outside of the plot
        spacing sets the spacing between the elements of the grid
        If a window was specified (self.win is not None) then the grid layout will derive from that window's central 
        item; otherwise we just make a gridLayout that can be put into another container somewhere.
        """
        import string
        if self.win is not None:
            self.gridLayout = self.win.ci.layout  # the window's 'central item' is the main gridlayout.
        else:
            self.gridLayout = QtGui.QGridLayout()  # just create the grid layout to add to another item
        self.gridLayout.setContentsMargins(margins, margins, margins, margins)
        self.gridLayout.setSpacing(spacing)
        self.plots = [[0 for x in xrange(self.cols)] for x in xrange(self.rows)]
        i = 0
        for r in range(self.rows):
            for c in range(self.cols):
                self.plots[r][c] = self.win.addPlot(row=r, col=c) # pg.PlotWidget()
                if letters:
                    labelPanels(self.plots[r][c], label=self.sequential_letters[i], size='14pt', bold=True)
                if titles:
                    labelTitles(self.plots[r][c], title=self.sequential_letters[i], size='14pt', bold=False)

                self.rcmap[i] = (r, c)
                i += 1
                if i > 25:
                    i = 0
        self.labelEdges('T(s)', 'Y', edgeOnly=self.edges)

    def labelEdges(self, xlabel='T(s)', ylabel='Y', edgeOnly=True, **kwargs):
        """
        label the axes on the outer edges of the gridlayout, leaving the interior axes clean
        """
        (lastrow, lastcol) = self.rcmap[-1]
        i = 0
        for (r,c) in self.rcmap:
            if c == 0:
                ylab = ylabel
            elif edgeOnly:
                ylab = ''
            else:
                ylab = ylabel
            if r == self.rows-1:  # only the last row
                xlab = xlabel
            elif edgeOnly:  # but not other rows
                xlab = ''
            else:
                xlab = xlabel  # otherwise, label it
            labelAxes(self.plots[r][c], xlab, ylab, **kwargs)
            i += 1


    def axesEdges(self, edgeOnly=True):
        """
        text labesls only on the axes on the outer edges of the gridlayout, 
        leaving the interior axes clean
        """
        (lastrow, lastcol) = self.rcmap[-1]
        i = 0
        for (r,c) in self.rcmap:
            xshow = True
            yshow = True
            if edgeOnly and c > 0:
                yshow = False
            if edgeOnly and r < self.rows:  # only the last row
                yshow = False
            ax = self.getPlot((r,c))
            leftaxis = ax.getAxis('left')
            bottomaxis = ax.getAxis('bottom')
            #print dir(self.plots[r][c])
            leftaxis.showValues = yshow
            bottomaxis.showValues = xshow
            i += 1


    def columnAutoScale(self, col, axis='left'): 
        """
        autoscale the columns according to the max value in the column.
        Finds outside range of column data, then sets the scale of all plots
        in the column to that range
        """
        atmax = None
        atmin = None
        for (r, c) in self.rcmap:
            if c != col:
                continue
            ax = self.getPlot((r,c))
            thisaxis = ax.getAxis(axis)
            amin, amax = thisaxis.range
            if atmax is None:
                atmax = amax
            else:
                if amax > atmax:
                    atmax = amax
            if atmin is None:
                atmin = amin
            else:
                if amin > atmin:
                    atmin = amin
            
        self.columnSetScale(col, axis=axis, range=(atmin, atmax))
        return(atmin, atmax)
                
    def columnSetScale(self, col, axis='left', range=(0., 1.)):
        """
        Set the column scale
        """
        for (r, c) in self.rcmap:
            if c != col:
                continue
            ax = self.getPlot((r,c))
            if axis == 'left':
                ax.setYRange(range[0], range[1])
            elif axis == 'bottom':
                ax.setXRange(range[0], range[1])

            if self.ticks == 'talbot':
                talbotTicks(ax)

        
    def title(self, index, title='', **kwargs):
        """
        add a title to a specific plot (specified by index) in the layout
        """
        labelTitles(self.getPlot(index), title=title, **kwargs)

def figure(title = None, background='w'):
    if background == 'w':
        pg.setConfigOption('background', 'w')  # set background to white
        pg.setConfigOption('foreground', 'k')
    pg.mkQApp()
    win = pg.GraphicsWindow(title=title)
    return win


def show():
    QtGui.QApplication.instance().exec_()


def test_layout(win):
    """
    Test the various plot types and modifications provided by the helpers above,
    in the context of a layout with various kinds of plots.
    """
    layout = LayoutMaker(cols=4,rows=2, win=win, labelEdges=True, ticks='talbot')
    x=np.arange(0, 10., 0.1)
    y = np.sin(x*3.)  # make an interesting signal
    r = np.random.random(10)  # and a random signal
    theta = np.linspace(0, 2.*np.pi, 10, endpoint=False)  # r, theta for polar plots
    for n in range(4*2):
        if n not in [1,2,3,4]:
            layout.plot(n, x, y)
        p = layout.getPlot(n)
        if n == 0:  # crossed axes plot
            crossAxes(p, xyzero=[5., 0.], density=(0.75, 1.5), tickPlacesAdd=(1, 0), pointSize=12)
            layout.title(n, 'Crossed Axes')
        if n in [1,2,3]:  # two differnt forms of polar plots
            if n == 1:
                po = polarPlot(p)
                po.setAxes(rMax=np.max(r))
                po.plot(r, theta, pen=pg.mkPen('r'))
                layout.title(n, 'Polar Path')
                
            if n == 2: 
                po = polarPlot(p)
                po.plot(r, theta, vectors=True, pen=pg.mkPen('k', width=2.0))
                po.setAxes(rMax=np.max(r))
                po.plot([np.mean(r)], [po.circmean(theta)],
                 vectors=True, pen=pg.mkPen('r', width=2.0))
                layout.title(n, 'Polar Arrows')
            if n == 3:
                po = polarPlot(p)
                po.hist(r, theta, binwidth=np.pi/6., normalize=False, density=False, pen='r')
                po.hist(r, theta, binwidth=np.pi/6., normalize=False, density=False, mode='arc', pen='b')
                po.setAxes(rMax=None)
                layout.title(n, 'Polar Histogram')
                
        if n == 4:  # violin plot with scatter plot data
            data = {2: [3,5,7,9,2,4,6,8,7,2,3,1,2.5], 3: [5, 6, 7, 9, 2, 8, 10, 9.5, 11]}
            violinPlotScatter(p, data, symbolColor='r')
            p.setYRange(0, 12)
            layout.title(n, 'Violin Plots with PseudoScatter')
            
        if n == 5:  # clean plot for physiology with baseline reference and a calibration bar
            calbar(p, calbar=[7.0, -1.5, 2.0, 0.5], axesoff=True, orient='left', unitNames={'x': 'ms', 'y': 'nA'})
            refline(p, refline=0., color = [64, 64, 64], linestyle = '--' ,linewidth = 0.5)
            layout.title(n, 'Calbar and Refline')
    
    
    #talbotTicks(layout.getPlot(1))
    layout.columnAutoScale(col=3, axis='left')
    show()


def test_crossAxes(win):
    layout = LayoutMaker(cols=1,rows=1, win=win, labelEdges=True)
    x=np.arange(-1, 1., 0.01)
    y = np.sin(x*10.)
    layout.plot(0, x, y)
    p = layout.getPlot(0)
    crossAxes(p, xyzero=[0., 0.], limits=[None, None, None, None], density=1.5, tickPlacesAdd=1, pointSize=12)
    show()


def test_polarPlot(win):
    layout = LayoutMaker(cols=1,rows=1, win=win, labelEdges=True)
    po = polarPlot(layout.getPlot((0,0)))  # convert rectangular plot to polar
    po.setAxes(steps=4, rMax=100, makeGrid=True)  # build the axes
    nvecs = 50
    #th = np.linspace(-np.pi*2, np.pi*2-np.pi*2/nvecs, nvecs)
    th = np.linspace(-np.pi*4, 0, nvecs)
    r = np.linspace(10, 100, nvecs)
    po.plot(r, th, vectors=True, arrowhead=True, symbols='o', pen=pg.mkPen('k', width=1.5))  # plot with arrowheads
    nvecs=8
    th = np.linspace(-np.pi*2, np.pi*2-np.pi*2/nvecs, nvecs)
    r = np.linspace(10, 100, nvecs)
    #po.plot(r, th, vectors=True, arrowhead=False, symbols='o', pen=pg.mkPen('r', width=1.5))  # plot with just lines
    
    show()
    

if __name__ == '__main__':
    win = figure(title='testing')
    test_layout(win)
    #test_crossAxes(win)
    #test_polarPlot(win)
