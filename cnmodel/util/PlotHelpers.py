#!/usr/bin/env python
# encoding: utf-8
"""
PlotHelpers.py

Routines to help use matplotlib and make cleaner plots
as well as get plots ready for publication. 

Modified to allow us to use a list of axes, and operate on all of those, 
or to use just one axis if that's all that is passed.
Therefore, the first argument to these calls can either be an axes object,
or a list of axes objects.  2/10/2012 pbm.

Plotter class: a simple class for managing figures with multiple plots.
Uses gridspec to build sets of axes. 

Created by Paul Manis on 2010-03-09.
Copyright 2010-2016  Paul Manis
Distributed under MIT/X11 license. See license.txt for more infofmation.

"""

import sys
import os
import string
from collections import OrderedDict

stdFont = 'Arial'
from matplotlib.ticker import FormatStrFormatter
from matplotlib.font_manager import FontProperties
from matplotlib.offsetbox import AnchoredOffsetbox, TextArea, DrawingArea, HPacker
from scipy.stats import gaussian_kde
import numpy as np
import matplotlib.pyplot as mpl
import matplotlib.gridspec as gridspec
from matplotlib.patches import Circle
from matplotlib.patches import Rectangle
from matplotlib.patches import Ellipse
from matplotlib.collections import PatchCollection
import matplotlib
rcParams = matplotlib.rcParams
rcParams['svg.fonttype'] = 'none' # No text as paths. Assume font installed.
rcParams['pdf.fonttype'] = 42
rcParams['ps.fonttype'] = 42
#rcParams['font.serif'] = ['Times New Roman']
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Arial']})
#rcParams['font.sans-serif'] = ['Arial']
#rcParams['font.family'] = 'sans-serif'
# check for LaTeX install - 
from distutils.spawn import find_executable
latex_avail = False
if find_executable('latex'):
    latex_avail = True
rc('text', usetex=latex_avail)
rcParams['text.latex.unicode'] = latex_avail

def _ax_tolist(ax):
    if isinstance(ax, list):
        return(ax)
    elif isinstance(ax, dict):
        axlist = list(axl.keys())
        return([ax for ax in axl[axlist]])
    else:
        return([ax])
    
    
def nice_plot(axl, spines=['left', 'bottom'], position=10, direction='inward', axesoff=False):
    """ Adjust a plot so that it looks nicer than the default matplotlib plot.
        Also allow quickaccess to things we like to do for publication plots, including:
           using a calbar instead of an axes: calbar = [x0, y0, xs, ys]
           inserting a reference line (grey, 3pt dashed, 0.5pt, at refline = y position)
    
    Parameters
    ----------
    axl : list of axes objects
        If a single axis object is present, it will be converted to a list here.
    
    spines : list of strings (default : ['left', 'bottom'])
        Sets whether spines will occur on particular axes. Choices are 'left', 'right',
        'bottom', and 'top'. Chosen spines will be displayed, others are not
    
    position : float (default : 10)
        Determines position of spines in points, typically outward by x points. The
        spines are the main axes lines, not the tick marks
        if the position is dict, then interpret as such.
        
    direction : string (default : 'inward')
        Sets the direction of spines. Choices are 'inward' and 'outward'

    axesoff : boolean (default : False)
        If true, forces the axes to be turned completely off.
    
    Returns
    -------
        Nothing.
    """
    #print 'NICEPLOT'
    if type(axl) is not list:
        axl = [axl]
    for ax in axl:
        if ax is None:
            continue
        #print 'ax: ', ax
        for loc, spine in ax.spines.items():
            if loc in spines:
                spine.set_color('k')
                #print 'spine color : k'
                if type(position) in [int, float]:
                    spine.set_position(('axes', position))
                elif type(position) is dict:
                    spine.set_position(('axes', position[loc]))
                else:
                    raise ValueError("position must be int, float or dict [ex: ]{'left': -0.05, 'bottom': -0.05}]")
            else:
                spine.set_color('none')
                #print 'spine color : none'
        if axesoff is True:
            noaxes(ax)

        # turn off ticks where there is no spine, if there are axes
        if 'left' in spines and not axesoff:
            ax.yaxis.set_ticks_position('left')
            ax.yaxis.set_tick_params(color='k')
        else:
            ax.yaxis.set_ticks([]) # no yaxis ticks

        if 'bottom' in spines and not axesoff:
            ax.xaxis.set_ticks_position('bottom')
            ax.xaxis.set_tick_params(color='k')
        else:
            ax.xaxis.set_ticks([])  # no xaxis ticks

        if direction == 'inward':
            ax.tick_params(axis='y', direction='in')
            ax.tick_params(axis='x', direction='in')
        else:
            ax.tick_params(axis='y', direction='out')
            ax.tick_params(axis='x', direction='out')


def noaxes(axl, whichaxes = 'xy'):
    """ take away all the axis ticks and the lines
    
    Parameters
    ----------
    
    axl : list of axes objects
        If a single axis object is present, it will be converted to a list here.
    
    whichaxes : string (default : 'xy')
        Sets which axes are turned off. The presence of an 'x' in 
        the string turns off x, the presence of 'y' turns off y.
    
    Returns
    -------
        Nothing
    """
    if type(axl) is not list:
        axl = [axl]
    for ax in axl:
        if ax is None:
            continue
        if 'x' in whichaxes:
            ax.xaxis.set_ticks([])
        if 'y' in whichaxes:
            ax.yaxis.set_ticks([])
        if 'xy' == whichaxes:
            ax.set_axis_off()


def setY(ax1, ax2):
    """
    Set the Y limits for an axes from a source axes to 
    the target axes.
    
    Parameters
    ----------
    
    ax1 : axis object
        The source axis object
    ax2 : list of axes objects
        If a single axis object is present, it will be converted to a list here.
        These are the target axes objects that will take on the limits of the source.
    
    Returns
    -------
        Nothing
    
    """
    if type(ax1) is list:
        print('PlotHelpers: cannot use list as source to set Y axis')
        return
    ax2 = _ax_tolist(ax2)
    # if type(ax2) is not list:
    #     ax2 = [ax2]
    refy = ax1.get_ylim()
    for ax in ax2:
        ax.set_ylim(refy)


def setX(ax1, ax2):
    """
    Set the X limits for an axes from a source axes to 
    the target axes.
    
    Parameters
    ----------
    
    ax1 : axis object
        The source axis object
    ax2 : list of axes objects
        If a single axis object is present, it will be converted to a list here.
        These are the target axes objects that will take on the limits of the source.
    
    Returns
    -------
        Nothing
    
    """
    if type(ax1) is list:
        print('PlotHelpers: cannot use list as source to set Y axis')
        return
    ax2 = _ax_tolist(ax2)
    # if type(ax2) is not list:
    #     ax2 = [ax2]
    refx = ax1.get_xlim()
    for ax in ax2:
        ax.set_xlim(refx)


def labelPanels(axl, axlist=None, font='Arial', fontsize=18, weight='normal', xy=(-0.05, 1.05), 
        horizontalalignment='right', verticalalignment='bottom', rotation=0.):
    """
    Provide labeling of panels in a figure with multiple subplots (axes)
    
    Parameters
    ----------
    axl : list of axes objects
        If a single axis object is present, it will be converted to a list here.
    
    axlist : list of string labels (default : None)
        Contains a list of the string labels. If the default value is provided,
        the axes will be lettered in alphabetical sequence. 
    
    font : string (default : 'Arial')
        Name of a valid font to use for the panel labels
    
    fontsize : float (default : 18, in points)
        Font size to use for axis labeling
    
    weight : string (default : 'normal')
        Font weight to use for labels. 'Bold', 'Italic', and 'Normal' are options
    
    xy : tuple (default : (-0.05, 1.05))
        A tuple (x,y) indicating where the label should go relative to the axis frame.
        Values are normalized as a fraction of the frame size.
    
    Returns
    -------
        list of the annotations

    """
    if isinstance(axl, dict):
        axlist = list(axl.keys())
    axl = _ax_tolist(axl)
    # if isinstance(axl, dict):
    #     axt = [axl[x] for x in axl]
    #     axlist = axl.keys()
    #     axl = axt
    # if not isinstance(axl, list):
    #     axl = [axl]
    if axlist is None:
        axlist = string.ascii_uppercase[0:len(axl)]
    # assume we wish to go in sequence
    if len(axlist) > len(axl):
        raise ValueError ('axl must have more entries than axlist: got axl=%d and axlist=%d for axlist:' % (len(axl), len(axlist)), axlist)
    font = FontProperties()
    font.set_family('sans-serif')
    font.set_weight=weight
    font.set_size=fontsize
    font.set_style('normal')
    labels = []
    for i, ax in enumerate(axl):
        if i >= len(axlist):
            continue
        if ax is None:
            continue
        if isinstance(ax, list):
            ax = ax[0]
        ann = ax.annotate(axlist[i], xy=xy, xycoords='axes fraction',
                annotation_clip=False,
                color="k", verticalalignment=verticalalignment,weight=weight, horizontalalignment=horizontalalignment,
                fontsize=fontsize, family='sans-serif', rotation=rotation
                )
        labels.append(ann)
    return(labels)


def listAxes(axd):
    """
    make a list of the axes from the dictionary
    """
    if type(axd) is not dict:
        if type(axd) is list:
            return axd
        else:
            print('listAxes expects dictionary or list; type not known (fix the code)')
            raise
    axl = [axd[x] for x in axd]
    return axl


def cleanAxes(axl):
    axl = _ax_tolist(axl)
    for ax in axl:
        if ax is None:
            continue
        for loc, spine in ax.spines.items():
            if loc in ['left', 'bottom']:
                spine.set_visible(True)
            elif loc in ['right', 'top']:
                spine.set_visible(False) 
                # spine.set_color('none')
                # do not draw the spine
            else:
                raise ValueError('Unknown spine location: %s' % loc)
            # turn off ticks when there is no spine
            ax.xaxis.set_ticks_position('bottom')
            #pdb.set_trace()
            ax.yaxis.set_ticks_position('left')  # stopped working in matplotlib 1.10
        update_font(ax)


def setTicks(axl, axis='x', ticks=np.arange(0, 1.1, 1.0)):
    axl = _ax_tolist(axl)
    # if type(axl) is dict:
    #     axl = [axl[x] for x in axl.keys()]
    # if type(axl) is not list:
    #     axl = [axl]
    for ax in axl:
        if ax is None:
            continue
        if axis == 'x':
            ax.set_xticks(ticks)
        if axis == 'y':
            ax.set_yticks(ticks)


def formatTicks(axl, axis='xy', fmt='%d', font='Arial'):
    """
    Convert tick labels to integers
    To do just one axis, set axis = 'x' or 'y'
    Control the format with the formatting string
    """
    axl = _ax_tolist(axl)
    # if type(axl) is not list:
    #     axl = [axl]
    majorFormatter = FormatStrFormatter(fmt)
    for ax in axl:
        if ax is None:
            continue
        if 'x' in axis:
            ax.xaxis.set_major_formatter(majorFormatter)
        if 'y' in axis:
            ax.yaxis.set_major_formatter(majorFormatter)


def autoFormatTicks(axl, axis='xy', font='Arial'):
    axl = _ax_tolist(axl)
    # if type(axl) is not list:
    #     axl = [axl]
    for ax in axl:
        if ax is None:
            continue
        if 'x' in axis:
        #    print ax.get_xlim()
            x0, x1= ax.get_xlim()
            setFormatter(ax,  x0, x1, axis = 'x')
        if 'y' in axis:
            y0, y1= ax.get_xlim
            setFormatter(ax, y0, y1, axis = 'y')


def setFormatter(axl, x0, x1, axis='x'):
    axl = _ax_tolist(axl)
    datarange = np.abs(x0-x1)
    mdata = np.ceil(np.log10(datarange))
    if mdata > 0 and mdata <= 4:
        majorFormatter = FormatStrFormatter('%d')
    elif mdata > 4:
        majorFormatter = FormatStrFormatter('%e')
    elif mdata <= 0 and mdata > -1:
        majorFormatter = FormatStrFormatter('%5.1f')
    elif mdata < -1 and mdata > -3:
        majorFormatatter = FormatStrFormatter('%6.3f')
    else:
        majorFormatter = FormatStrFormatter('%e')
    for ax in axl:
        if axis == 'x':
            ax.xaxis.set_major_formatter(majorFormatter)
        elif axis == 'y':
            ax.yaxis.set_major_formatter(majorFormatter)


def update_font(axl, size=9, font=stdFont):
    axl = _ax_tolist(axl)
    # if type(axl) is not list:
    #     axl = [axl]
    fontProperties = {'family':'sans-serif', #'sans-serif': font,
            'weight' : 'normal', 'size' : size}
    for ax in axl:
        if ax is None:
            continue
        for tick in ax.xaxis.get_major_ticks():
              #tick.label1.set_family('sans-serif')
            #  tick.label1.set_fontname(stdFont)
              tick.label1.set_size(size)

        for tick in ax.yaxis.get_major_ticks():
             # tick.label1.set_family('sans-serif')
            #  tick.label1.set_fontname(stdFont)
              tick.label1.set_size(size)
        ax.set_xticklabels(ax.get_xticks(), fontProperties)
        ax.set_yticklabels(ax.get_yticks(), fontProperties)
        ax.xaxis.set_smart_bounds(True)
        ax.yaxis.set_smart_bounds(True) 
        ax.tick_params(axis = 'both', labelsize = size)


def lockPlot(axl, lims, ticks=None):
    """ 
        This routine forces the plot of invisible data to force the axes to take certain
        limits and to force the tick marks to appear. 
        call with the axis and lims (limits) = [x0, x1, y0, y1]
    """
    axl = _ax_tolist(axl)
    # if type(axl) is not list:
    #     axl = [axl]
    plist = []
    for ax in axl:
        if ax is None:
            continue
        lpl = ax.plot([lims[0], lims[0], lims[1], lims[1]], [lims[2], lims[3], lims[2], lims[3]],
            color='none', marker='', linestyle='None')
        plist.extend(lpl)
        ax.axis(lims)
    return(plist)  # just in case you want to modify these plots later.


def adjust_spines(axl, spines = ['left', 'bottom'], direction = 'outward', distance=5, smart=True):
    axl = _ax_tolist(axl)
    # if type(axl) is not list:
    #     axl = [axl]
    for ax in axl:
        if ax is None:
            continue
        # turn off ticks where there is no spine
        if 'left' in spines:
            ax.yaxis.set_ticks_position('left')
        else:
            # no yaxis ticks
            ax.yaxis.set_ticks([])

        if 'bottom' in spines:
            ax.xaxis.set_ticks_position('bottom')
        else:
            # no xaxis ticks
            ax.xaxis.set_ticks([])
        for loc, spine in ax.spines.items():
            if loc in spines:
                spine.set_position((direction,distance)) # outward by 10 points
                if smart is True:
                    spine.set_smart_bounds(True)
                else:
                    spine.set_smart_bounds(False)
            else:
                spine.set_color('none')  # don't draw spine


def getLayoutDimensions(n, pref='height'):
    """
    Return a tuple of optimized layout dimensions for n axes
    
    Parameters
    ----------
    n : int (no default):
        Number of plots needed
    pref : string (default : 'height')
        prefered way to organized the plots (height, or width)
    
    Returns
    -------
    (h, w) : tuple
        height (rows) and width (columns)
    
    """
    nopt = np.sqrt(n)
    inoptw = int(nopt)
    inopth = int(nopt)
    while inoptw*inopth < n:
        if pref == 'width':
            inoptw += 1
            if inoptw * inopth > (n-inopth):
                inoptw -= 1
                inopth += 1
        else:
            inopth += 1
            if inoptw * inopth > (n-inoptw):
                inopth -= 1
                inoptw += 1
            
    return(inopth, inoptw)


def calbar(axl, calbar=None, axesoff=True, orient='left', unitNames=None, fontsize=11, weight='normal', font='Arial'):
    """
        draw a calibration bar and label it. The calibration bar is defined as:
        [x0, y0, xlen, ylen]
    """
    axl = _ax_tolist(axl)
    # if type(axl) is not list:
    #     axl = [axl]
    for ax in axl:
        if ax is None:
            continue
        if axesoff is True:
            noaxes(ax)
        Hfmt = r'{:.0f}'
        if calbar[2] < 1.0:
            Hfmt = r'{:.1f}'
        Vfmt = r' {:.0f}'
        if calbar[3] < 1.0:
            Vfmt = r' {:.1f}'
        if unitNames is not None:
            Vfmt = Vfmt + r' ' + r'{:s}'.format(unitNames['y'])
            Hfmt = Hfmt + r' ' + r'{:s}'.format(unitNames['x'])
        # print(Vfmt, unitNames['y'])
        # print(Vfmt.format(calbar[3]))
        font = FontProperties()
        font.set_family('sans-serif')
        font.set_weight=weight
        font.set_size=fontsize
        font.set_style('normal')
        if calbar is not None:
            if orient == 'left':  # vertical part is on the left
                ax.plot([calbar[0], calbar[0], calbar[0]+calbar[2]], 
                    [calbar[1]+calbar[3], calbar[1], calbar[1]],
                    color = 'k', linestyle = '-', linewidth = 1.5)
                ax.text(calbar[0]+0.05*calbar[2], calbar[1]+0.5*calbar[3], Vfmt.format(calbar[3]), 
                    horizontalalignment = 'left', verticalalignment = 'center',
                    fontsize = fontsize, weight=weight, family='sans-serif',)
            elif orient == 'right':  # vertical part goes on the right
                ax.plot([calbar[0] + calbar[2], calbar[0]+calbar[2], calbar[0]], 
                    [calbar[1]+calbar[3], calbar[1], calbar[1]],
                    color = 'k', linestyle = '-', linewidth = 1.5)
                ax.text(calbar[0]+calbar[2]-0.05*calbar[2], calbar[1]+0.5*calbar[3], Vfmt.format(calbar[3]), 
                    horizontalalignment = 'right', verticalalignment = 'center',
                    fontsize = fontsize, weight=weight, family='sans-serif',)
            else:
                print("PlotHelpers.py: I did not understand orientation: %s" % (orient))
                print("plotting as if set to left... ")
                ax.plot([calbar[0], calbar[0], calbar[0]+calbar[2]], 
                    [calbar[1]+calbar[3], calbar[1], calbar[1]],
                    color = 'k', linestyle = '-', linewidth = 1.5)
                ax.text(calbar[0]+0.05*calbar[2], calbar[1]+0.5*calbar[3],Vfmt.format(calbar[3]), 
                    horizontalalignment = 'left', verticalalignment = 'center',
                    fontsize = fontsize, weight=weight, family='sans-serif',)
            ax.text(calbar[0]+calbar[2]*0.5, calbar[1]-0.1*calbar[3], Hfmt.format(calbar[2]), 
                horizontalalignment = 'center', verticalalignment = 'top',
                fontsize = fontsize, weight=weight, family='sans-serif',)


def referenceline(axl, reference=None, limits=None, color='0.33', linestyle='--' ,linewidth=0.5, dashes=None):
    """
    draw a reference line at a particular level of the data on the y axis
    returns the line object.
    """
    axl = _ax_tolist(axl)
    # if type(axl) is not list:
    #     axl = [axl]
    if reference is None:
        refeference = 0.
    for ax in axl:
        if ax is None:
            continue
        if limits is None or type(limits) is not list or len(limits) != 2:
            xlims = ax.get_xlim()
        else:
            xlims = limits
        rl, = ax.plot([xlims[0], xlims[1]], [reference, reference],
             color=color, linestyle=linestyle, linewidth=linewidth)
        if dashes is not None:
            rl.set_dashes(dashes)
    return rl


def crossAxes(axl, xyzero=[0., 0.], limits=[None, None, None, None]):
    """
    Make plot(s) with crossed axes at the data points set by xyzero, and optionally
    set axes limits
    """
    axl = _ax_tolist(axl)
    # if type(axl) is not list:
    #     axl = [axl]
    for ax in axl:
        if ax is None:
            continue
#        ax.set_title('spines at data (1,2)')
#        ax.plot(x,y)
        ax.spines['left'].set_position(('data',xyzero[0]))
        ax.spines['right'].set_color('none')
        ax.spines['bottom'].set_position(('data',xyzero[1]))
        ax.spines['top'].set_color('none')
        ax.spines['left'].set_smart_bounds(True)
        ax.spines['bottom'].set_smart_bounds(True)
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')
        if limits[0] is not None:
            ax.set_xlim(left=limits[0], right=limits[2])
            ax.set_ylim(bottom=limits[1], top=limits[3])
            
def violin_plot(ax, data, pos, bp=False, median = False):
    '''
    create violin plots on an axis
    '''
    dist = max(pos)-min(pos)
    w = min(0.15*max(dist,1.0),0.5)
    for d,p in zip(data,pos):
        k = gaussian_kde(d)  #calculates the kernel density
        m = k.dataset.min()  #lower bound of violin
        M = k.dataset.max()  #upper bound of violin
        x = np.arange(m, M, (M-m)/100.)  # support for violin
        v = k.evaluate(x)  #violin profile (density curve)
        v = v / v.max() * w  #scaling the violin to the available space
        ax.fill_betweenx(x, p, v+p, facecolor='y', alpha=0.3)
        ax.fill_betweenx(x, p, -v+p, facecolor='y', alpha=0.3)
        if median:
            ax.plot([p-0.5, p+0.5], [np.median(d), np.median(d)], '-')
    if bp:
        bpf = ax.boxplot(data, notch=0, positions=pos, vert=1)
        mpl.setp(bpf['boxes'], color='black')
        mpl.setp(bpf['whiskers'], color='black', linestyle='-')


# # from somewhere on the web:

class NiceScale:
    def __init__(self, minv,maxv):
        self.maxTicks = 6
        self.tickSpacing = 0
        self.lst = 10
        self.niceMin = 0
        self.niceMax = 0
        self.minPoint = minv
        self.maxPoint = maxv
        self.calculate()

    def calculate(self):
        self.lst = self.niceNum(self.maxPoint - self.minPoint, False)
        self.tickSpacing = self.niceNum(self.lst / (self.maxTicks - 1), True)
        self.niceMin = np.floor(self.minPoint / self.tickSpacing) * self.tickSpacing
        self.niceMax = np.ceil(self.maxPoint / self.tickSpacing) * self.tickSpacing

    def niceNum(self, lst, rround):
        self.lst = lst
        exponent = 0 # exponent of range */
        fraction = 0 # fractional part of range */
        niceFraction = 0 # nice, rounded fraction */

        exponent = np.floor(np.log10(self.lst));
        fraction = self.lst / np.power(10, exponent);

        if (self.lst):
            if (fraction < 1.5):
                niceFraction = 1
            elif (fraction < 3):
                niceFraction = 2
            elif (fraction < 7):
                niceFraction = 5;
            else:
                niceFraction = 10;
        else :
            if (fraction <= 1):
                niceFraction = 1
            elif (fraction <= 2):
                niceFraction = 2
            elif (fraction <= 5):
                niceFraction = 5
            else:
                niceFraction = 10

        return niceFraction * np.power(10, exponent)

    def setMinMaxPoints(self, minPoint, maxPoint):
          self.minPoint = minPoint
          self.maxPoint = maxPoint
          self.calculate()

    def setMaxTicks(self, maxTicks):
        self.maxTicks = maxTicks;
        self.calculate()

def circles(x, y, s, c='b', ax=None, vmin=None, vmax=None, **kwargs):
    """
    Make a scatter of circles plot of x vs y, where x and y are sequence 
    like objects of the same lengths. The size of circles are in data scale.

    Parameters
    ----------
    x,y : scalar or array_like, shape (n, )
        Input data
    s : scalar or array_like, shape (n, ) 
        Radius of circle in data scale (ie. in data unit)
    c : color or sequence of color, optional, default : 'b'
        `c` can be a single color format string, or a sequence of color
        specifications of length `N`, or a sequence of `N` numbers to be
        mapped to colors using the `cmap` and `norm` specified via kwargs.
        Note that `c` should not be a single numeric RGB or
        RGBA sequence because that is indistinguishable from an array of
        values to be colormapped.  `c` can be a 2-D array in which the
        rows are RGB or RGBA, however.
    ax : Axes object, optional, default: None
        Parent axes of the plot. It uses gca() if not specified.
    vmin, vmax : scalar, optional, default: None
        `vmin` and `vmax` are used in conjunction with `norm` to normalize
        luminance data.  If either are `None`, the min and max of the
        color array is used.  (Note if you pass a `norm` instance, your
        settings for `vmin` and `vmax` will be ignored.)

    Returns
    -------
    paths : `~matplotlib.collections.PathCollection`

    Other parameters
    ----------------
    kwargs : `~matplotlib.collections.Collection` properties
        eg. alpha, edgecolors, facecolors, linewidths, linestyles, norm, cmap

    Examples
    --------
    a = np.arange(11)
    circles(a, a, a*0.2, c=a, alpha=0.5, edgecolor='none')

    License
    --------
    This code is under [The BSD 3-Clause License]
    (http://opensource.org/licenses/BSD-3-Clause)
    """

    #import matplotlib.colors as colors

    if ax is None:
        ax = mpl.gca()    

    if isinstance(c,str):
        color = c     # ie. use colors.colorConverter.to_rgba_array(c)
    else:
        color = None  # use cmap, norm after collection is created
    kwargs.update(color=color)

    if np.isscalar(x):
        patches = [Circle((x, y), s),]
    elif np.isscalar(s):
        patches = [Circle((x_,y_), s) for x_,y_ in zip(x,y)]
    else:
        patches = [Circle((x_,y_), s_) for x_,y_,s_ in zip(x,y,s)]
    collection = PatchCollection(patches, **kwargs)

    if color is None:
        collection.set_array(np.asarray(c))
        if vmin is not None or vmax is not None:
            collection.set_clim(vmin, vmax)

    ax.add_collection(collection)
    ax.autoscale_view()
    return collection


def rectangles(x, y, sw, sh=None, c='b', ax=None, vmin=None, vmax=None, **kwargs):
    """
    Make a scatter of squares plot of x vs y, where x and y are sequence 
    like objects of the same lengths. The size of sqares are in data scale.

    Parameters
    ----------
    x,y : scalar or array_like, shape (n, )
        Input data
    s : scalar or array_like, shape (n, ) 
        side of square in data scale (ie. in data unit)
    c : color or sequence of color, optional, default : 'b'
        `c` can be a single color format string, or a sequence of color
        specifications of length `N`, or a sequence of `N` numbers to be
        mapped to colors using the `cmap` and `norm` specified via kwargs.
        Note that `c` should not be a single numeric RGB or
        RGBA sequence because that is indistinguishable from an array of
        values to be colormapped.  `c` can be a 2-D array in which the
        rows are RGB or RGBA, however.
    ax : Axes object, optional, default: None
        Parent axes of the plot. It uses gca() if not specified.
    vmin, vmax : scalar, optional, default: None
        `vmin` and `vmax` are used in conjunction with `norm` to normalize
        luminance data.  If either are `None`, the min and max of the
        color array is used.  (Note if you pass a `norm` instance, your
        settings for `vmin` and `vmax` will be ignored.)

    Returns
    -------
    paths : `~matplotlib.collections.PathCollection`

    Other parameters
    ----------------
    kwargs : `~matplotlib.collections.Collection` properties
        eg. alpha, edgecolors, facecolors, linewidths, linestyles, norm, cmap

    Examples
    --------
    a = np.arange(11)
    squaress(a, a, a*0.2, c=a, alpha=0.5, edgecolor='none')

    License
    --------
    This code is under [The BSD 3-Clause License]
    (http://opensource.org/licenses/BSD-3-Clause)
    """
    #import matplotlib.colors as colors

    if ax is None:
        ax = mpl.gca()    

    if isinstance(c,str):
        color = c     # ie. use colors.colorConverter.to_rgba_array(c)
    else:
        color = None  # use cmap, norm after collection is created
    kwargs.update(color=color)
    if sh is None:
        sh = sw
    x = x - sw/2.  # offset as position specified is "lower left corner"
    y = y - sh/2.
    if np.isscalar(x):
        patches = [Rectangle((x, y), sw, sh),]
    elif np.isscalar(sw):
        patches = [Rectangle((x_,y_), sw, sh) for x_,y_ in zip(x,y)]
    else:
        patches = [Rectangle((x_,y_), sw_, sh_) for x_,y_,sw_,sh_ in zip(x,y,sw,sh)]
    collection = PatchCollection(patches, **kwargs)

    if color is None:
        collection.set_array(np.asarray(c))
        if vmin is not None or vmax is not None:
            collection.set_clim(vmin, vmax)

    ax.add_collection(collection)
    ax.autoscale_view()
    return collection


def show_figure_grid(fig, figx=10., figy=10.):
    """
    Create a background grid with major and minor lines like graph paper
    if using default figx and figy, the grid will be in units of the 
    overall figure on a [0,1,0,1] grid
    if figx and figy are in units of inches or cm, then the grid
    will be on that scale.
    
    Figure grid is useful when building figures and placing labels
    at absolute locations on the figure.
    
    Parameters
    ----------
    
    fig : Matplotlib figure handle (no default):
        The figure to which the grid will be applied
    
    figx : float (default: 10.)
        # of major lines along the X dimension
    
    figy : float (default: 10.)
        # of major lines along the Y dimension
    
    """
    backGrid = fig.add_axes([0,0,1,1], frameon=False)
    backGrid.set_ylim(0., figy)
    backGrid.set_xlim(0., figx)
    backGrid.grid(True)

    backGrid.set_yticks(np.arange(0., figy+0.01, 1.))
    backGrid.set_yticks(np.arange(0., figy+0.01, 0.1), minor=True)
    backGrid.set_xticks(np.arange(0., figx+0.01, 1.))
    backGrid.set_xticks(np.arange(0., figx+0.01, 0.1), minor=True)
#   backGrid.get_xaxis().set_minor_locator(matplotlib.ticker.AutoMinorLocator())
#   backGrid.get_yaxis().set_minor_locator(matplotlib.ticker.AutoMinorLocator())
    backGrid.grid(b=True, which='major', color='g', alpha=0.6, linewidth=0.8)
    backGrid.grid(b=True, which='minor', color='g', alpha=0.4, linewidth=0.2)
    return backGrid

def hide_figure_grid(fig, grid):
    grid.grid(False)

def delete_figure_grid(fig, grid):
    mpl.delete(grid)


class Plotter():
    """
    The Plotter class provides a simple convenience for plotting data in 
    an row x column array.
    """
    def __init__(self, rcshape=None, axmap=None, arrangement=None, title=None, label=False, roworder=True, refline=None,
        figsize=(11, 8.5), fontsize=10, position=0, labeloffset=[0., 0.], labelsize=12):
        """
        Create an instance of the plotter. Generates a new matplotlib figure,
        and sets up an array of subplots as defined, initializes the counters
        
        Examples
        --------
        Ex. 1: 
        One way to generate plots on a standard grid, uses gridspec to specify an axis map:
            labels = ['A', 'B1', 'B2', 'C1', 'C2', 'D', 'E', 'F']
            gr = [(0, 4, 0, 1), (0, 3, 1, 2), (3, 4, 1, 2), (0, 3, 2, 3), (3, 4, 2, 3), (5, 8, 0, 1), (5, 8, 1, 2), (5, 8, 2, 3)]
            axmap = OrderedDict(zip(labels, gr))
            P = PH.Plotter((8, 1), axmap=axmap, label=True, figsize=(8., 6.))
            PH.show_figure_grid(P.figure_handle)
        
        Ex. 2:
        Place plots on defined locations on the page - no messing with gridspec or subplots. 
            For this version, we just generate N subplots with labels (used to tag each plot)
            The "sizer" array then maps the tags to specific panel locations
            # define positions for each panel in Figure coordinages (0, 1, 0, 1)
            # you don't have to use an ordered dict for this, I just prefer it when debugging
            sizer = OrderedDict([('A', [0.08, 0.22, 0.55, 0.4]), ('B1', [0.40, 0.25, 0.65, 0.3]), ('B2', [0.40, 0.25, 0.5, 0.1]),
                    ('C1', [0.72, 0.25, 0.65, 0.3]), ('C2', [0.72, 0.25, 0.5, 0.1]),
                    ('D', [0.08, 0.25, 0.1, 0.3]), ('E', [0.40, 0.25, 0.1, 0.3]), ('F', [0.72, 0.25, 0.1, 0.3]),
            ])  # dict elements are [left, width, bottom, height] for the axes in the plot.
            gr = [(a, a+1, 0, 1) for a in range(0, 8)]   # just generate subplots - shape does not matter
            axmap = OrderedDict(zip(sizer.keys(), gr))
            P = PH.Plotter((8, 1), axmap=axmap, label=True, figsize=(8., 6.))
            PH.show_figure_grid(P.figure_handle)
            P.resize(sizer)  # perform positioning magic
            P.axdict['B1'] access the plot associated with panel B1

        Parameters
        ----------
        rcshape : a list or tuple: 2x1 (no default)
                    rcshape is an array [row, col] telling us how many rows and columns to build.
                    default defines a rectangular array r x c of plots
                  a dict : 
                  None: expect axmap to provide the input... 
        
        axmap : 
            list of gridspec slices (default : None)
            define slices for the axes of a gridspec, allowing for non-rectangular arrangements
            The list is defined as:
            [(r1t, r1b, c1l, c1r), slice(r2, c2)]
            where r1t is the top for row 1 in the grid, r1b is the bottom, etc... 
            When using this mode, the axarr returned is a 1-D list, as if r is all plots indexed,
            and the number of columns is 1. The results match in order the list entered in axmap
        

        
        arrangement: Ordered Dict (default: None)
            Arrangement allows the data to be plotted according to a logical arrangement
            The dict keys are the names ("groups") for each column, and the elements are
            string names for the entities in the groups
        
        title : string (default: None)
            Provide a title for the entire plot
        
        label : Boolean (default: False)
            If True, sets labels on panels
        
        roworder : Boolean (default: True)
            Define whether labels run in row order first or column order first
        
        refline : float (default: None)
            Define the position of a reference line to be used in all panels
        
        figsize : tuple (default : (11, 8.5))
            Figure size in inches. Default is for a landscape figure
        
        fontsize : points (default : 10)
            Defines the size of the font to use for panel labels

        position : position of spines (0 means close, 0.05 means break out)
            x, y spines.. 
        Returns
        -------
        Nothing
        """
        self.arrangement = arrangement
        self.fontsize = fontsize
        self.referenceLines = {}
        self.figure_handle = mpl.figure(figsize=figsize) # create the figure
        self.figure_handle.set_size_inches(figsize[0], figsize[1], forward=True)
        self.axlabels = []
        self.axdict = OrderedDict()  # make axis label dictionary for indirect access (better!)
        if isinstance(fontsize, int):
            fontsize = {'tick': fontsize, 'label': fontsize, 'panel': fontsize}
        gridbuilt = False
        # compute label offsets
        p = [0., 0.]
        if label:
            if type(labeloffset) is int:
                p = [labeloffset, labeloffset]
            elif type(labeloffset) is dict:
                p = [position['left'], position['bottom']]
            elif type(labeloffset) in [list, tuple]:
                p = labeloffset
            else:
                p = [0., 0.]
        
        # build axes arrays
        # 1. nxm grid
        if isinstance(rcshape, list) or isinstance(rcshape, tuple):
            rc = rcshape
            gs = gridspec.GridSpec(rc[0], rc[1])  # define a grid using gridspec
            # assign to axarr
            self.axarr = np.empty(shape=(rc[0], rc[1],), dtype=object)  # use a numpy object array, indexing features
            ix = 0
            for r in range(rc[0]):
                for c in range(rc[1]):
                    self.axarr[r,c] = mpl.subplot(gs[ix])
                    ix += 1
            gridbuilt = True
        # 2. specified values - starts with Nx1 subplots, then reorganizes according to shape boxes
        elif isinstance(rcshape, dict):  # true for OrderedDict also
            nplots = len(list(rcshape.keys()))
            gs = gridspec.GridSpec(nplots, 1)
            rc = (nplots, 1)
            self.axarr = np.empty(shape=(rc[0], rc[1],), dtype=object)  # use a numpy object array, indexing features
            ix = 0
            for r in range(rc[0]):
                for c in range(rc[1]):
                    self.axarr[r,c] = mpl.subplot(gs[ix])
                    ix += 1
            gridbuilt = True
            for k, pk in enumerate(rcshape.keys()):
                self.axdict[pk] = self.axarr[k,0]
            plo = labeloffset
            self.axlabels = labelPanels(self.axarr.tolist(), axlist=list(rcshape.keys()), xy=(-0.095+plo[0], 0.95+plo[1]), fontsize=fontsize['panel'])
            self.resize(rcshape)
        else:
            raise ValueError('Input rcshape must be list/tuple or dict')
            
        # create sublots
        if axmap is not None:
            if isinstance(axmap, list) and not gridbuilt:
                self.axarr = np.empty(shape=(len(axmap), 1), dtype=object)
                for k, g in enumerate(axmap):
                    self.axarr[k,] = mpl.subplot(gs[g[0]:g[1], g[2]:g[3]])
            elif isinstance(axmap, dict) or isinstance(axmap, OrderedDict): # keys are panel labels
                if not gridbuilt:
                    self.axarr = np.empty(shape=(len(list(axmap.keys())), 1), dtype=object)
                na = np.prod(self.axarr.shape)  # number of axes
                for k, pk in enumerate(axmap.keys()):
                    g = axmap[pk]  # get the gridspec info
                    if not gridbuilt:
                        self.axarr[k,] = mpl.subplot(gs[g[0]:g[1], g[2]:g[3]])
                    self.axdict[pk] = self.axarr.ravel()[k]
            else:
                raise TypeError('Plotter in PlotHelpers: axmap must be a list or dict')
 
        if len(self.axdict) == 0:
            for i, a in enumerate(self.axarr.flatten()):
                label = string.ascii_uppercase[i]
                self.axdict[label] = a
        
        if title is not None:
            self.figure_handle.canvas.set_window_title(title)
            self.figure_handle.suptitle(title)
        self.nrows = self.axarr.shape[0]
        if len(self.axarr.shape) > 1:
            self.ncolumns = self.axarr.shape[1]
        else:
            self.ncolumns = 1
        self.row_counter = 0
        self.column_counter = 0
        for i in range(self.nrows):
            for j in range(self.ncolumns):
                self.axarr[i, j].spines['top'].set_visible(False)
                self.axarr[i, j].get_xaxis().set_tick_params(direction='out', width=0.8, length=4.)
                self.axarr[i, j].get_yaxis().set_tick_params(direction='out', width=0.8, length=4.)
                self.axarr[i, j].tick_params(axis='both', which='major', labelsize=fontsize['tick'])
#                if i < self.nrows-1:
#                    self.axarr[i, j].xaxis.set_major_formatter(mpl.NullFormatter())
                nice_plot(self.axarr[i, j], position=position)
                if refline is not None:
                    self.referenceLines[self.axarr[i,j]] = referenceline(self.axarr[i,j], reference=refline)

        if label:
            if isinstance(axmap, dict) or isinstance(axmap, OrderedDict):  # in case predefined... 
                self.axlabels = labelPanels(self.axarr.ravel().tolist(), axlist=list(axmap.keys()), xy=(-0.095+p[0], 0.95+p[1]), fontsize=fontsize['panel'])
                return
            self.axlist = []
            if roworder == True:
                for i in range(self.nrows):
                    for j in range(self.ncolumns):
                        self.axlist.append(self.axarr[i, j])
            else:
                for i in range(self.ncolumns):
                    for j in range(self.nrows):
                        self.axlist.append(self.axarr[j, i])
                
            if self.nrows*self.ncolumns > 26:  # handle large plot using "A1..."
                ctxt = string.ascii_uppercase[0:self.ncolumns]  # columns are lettered
                rtxt = [str(x+1) for x in range(self.nrows)] # rows are numbered, starting at 1
                axl = []
                for i in range(self.nrows):
                    for j in range(self.ncolumns):
                        axl.append(ctxt[j]+rtxt[i])
                self.axlabels = labelPanels(self.axlist, axlist=axl, xy=(-0.35+p[0], 0.75))
            else:
                self.axlabels = labelPanels(self.axlist, xy=(-0.095+p[0], 0.95+p[1]))
    
    def _next(self):
        """
        Private function
        _next gets the axis pointer to the next row, column index that is available
        Only sets internal variables
        """
        self.column_counter += 1
        if self.column_counter >= self.ncolumns:
            self.row_counter += 1
            self.column_counter = 0
            if self.row_counter >= self.nrows:
                raise ValueError('Call to get next row exceeds the number of rows requested initially: %d' % self.nrows)
    
    def getaxis(self, group=None):
        """
        getaxis gets the current row, column counter, and calls _next to increment the counter
        (so that the next getaxis returns the next available axis pointer)
        
        Parameters
        ----------
        group : string (default: None)
            forces the current axis to be selected from text name of a "group"
        
        Returns
        -------
        the current axis or the axis associated with a group
        """
        
        if group is None:
            currentaxis = self.axarr[self.row_counter, self.column_counter]
            self._next() # prepare for next call
        else:
            currentaxis = self.getRC(group)
                
        return currentaxis
    
    def getRC(self, group):
        """
        Get the axis associated with a group
        
        Parameters
        ----------
        group : string (default: None)
            returns the matplotlib axis associated with a text name of a "group"
        
        Returns
        -------
        The matplotlib axis associated with the group name, or None if no group by
        that name exists in the arrangement
        """
        
        if self.arrangement is None:
            raise ValueError('specifying a group requires an arrangment dictionary')
        # look for the group label in the arrangement dicts
        for c, colname in enumerate(self.arrangement.keys()):
            if group in self.arrangement[colname]:
                # print ('column name, column: ', colname, self.arrangement[colname])
                # print ('group: ', group)
                r = self.arrangement[colname].index(group)  # get the row position this way
                return(self.axarr[r, c])
        print(('Group {:s} not in the arrangement'.format(group)))
        return None
        
        sizer = {'A': {'pos': [0.08, 0.22, 0.50, 0.4]}, 'B1': {'pos': [0.40, 0.25, 0.60, 0.3]}, 'B2': {'pos': [0.40, 0.25, 0.5, 0.1]},
                'C1': {'pos': [0.72, 0.25, 0.60, 0.3]}, 'C2': {'pos': [0.72, 0.25, 0.5, 0.1]},
                'D': {'pos': [0.08, 0.25, 0.1, 0.3]}, 'E': {'pos': [0.40, 0.25, 0.1, 0.3]}, 'F': {'pos': [0.72, 0.25, 0.1, 0.3]},
        }
    
    def resize(self, sizer):
        """
        Resize the graphs in the array.
        
        Parameters
        ----------
        sizer : dict (no default)
            A dictionary with keys corresponding to the plot labels. 
            The values for each key are a list (or tuple) of [left, width, bottom, height]
            for each panel in units of the graph [0, 1, 0, 1]. 
        
        sizer = {'A': {'pos': [0.08, 0.22, 0.50, 0.4], 'labelpos': (x,y), 'noaxes': True}, 'B1': {'pos': [0.40, 0.25, 0.60, 0.3], 'labelpos': (x,y)},
                 'B2': {'pos': [0.40, 0.25, 0.5, 0.1],, 'labelpos': (x,y), 'noaxes': False},
                'C1': {'pos': [0.72, 0.25, 0.60, 0.3], 'labelpos': (x,y)}, 'C2': {'pos': [0.72, 0.25, 0.5, 0.1], 'labelpos': (x,y)},
                'D': {'pos': [0.08, 0.25, 0.1, 0.3], 'labelpos': (x,y)}, 
                'E': {'pos': [0.40, 0.25, 0.1, 0.3], 'labelpos': (x,y)}, 'F': {'pos': [0.72, 0.25, 0.1, 0.3],, 'labelpos': (x,y)}
                }
        Returns
        -------
        Nothing
        """
        
        for i, s in enumerate(sizer.keys()):
            ax = self.axdict[s]
            bbox = ax.get_position()
            bbox.x0 = sizer[s]['pos'][0]
            bbox.x1 = sizer[s]['pos'][1]+ sizer[s]['pos'][0]
            bbox.y0 = sizer[s]['pos'][2]
            bbox.y1 = sizer[s]['pos'][3] + sizer[s]['pos'][2]  # offsets are in figure fractions
            ax.set_position(bbox)
            if 'labelpos' in list(sizer[s].keys()) and len(sizer[s]['labelpos']) == 2:
                x, y = sizer[s]['labelpos']
                self.axlabels[i].set_x(x)
                self.axlabels[i].set_y(y)
            if 'noaxes' in sizer[s] and sizer[s]['noaxes'] == True:
                noaxes(ax)
                


if __name__ == '__main__':
#    P = Plotter((3,3), axmap=[(0, 1, 0, 3), (1, 2, 0, 2), (2, 1, 2, 3), (2, 3, 0, 1), (2, 3, 1, 2)])
    labels = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I']
    l = [(a, a+2, 0, 1) for a in range(0, 6, 2)]
    r = [(a, a+1, 1, 2) for a in range(0, 6)]
    axmap = OrderedDict(list(zip(labels, l+r)))
    P = Plotter((6,2), axmap=axmap, figsize=(6., 6.), label=True)
#    P = Plotter((2,3), label=True)  # create a figure with plots
    # for a in P.axarr.flatten():
    #     a.plot(np.random.random(10), np.random.random(10))
    
#    hfig, ax = mpl.subplots(2, 3)
    axd = OrderedDict()
    for i, a in enumerate(P.axarr.flatten()):
        label = string.ascii_uppercase[i]
        axd[label] = a
    for a in list(axd.keys()):
        axd[a].plot(np.random.random(10), np.random.random(10))
    nice_plot([axd[a] for a in axd], position=-0.1)
    cleanAxes([axd['B'], axd['C']])
    calbar([axd['B'], axd['C']], calbar=[0.5, 0.5, 0.2, 0.2])
    #labelPanels([axd[a] for a in axd], axd.keys())
    #mpl.tight_layout(pad=2, w_pad=0.5, h_pad=2.0)
    mpl.show()
    
               