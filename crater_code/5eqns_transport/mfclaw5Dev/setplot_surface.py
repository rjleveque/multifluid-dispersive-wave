"""
Set up the plot figures, axes, and items to be done for each frame.

This module is imported by the plotting routines and then the
function setplot is called to set the plot parameters.

"""


import os


#--------------------------
def setplot(plotdata=None):
#--------------------------

    """
    Specify what is to be plotted at each frame.
    Input:  plotdata, an instance of clawpack.visclaw.data.ClawPlotData.
    Output: a modified version of plotdata.

    """


    from clawpack.visclaw import colormaps

    if plotdata is None:
        from clawpack.visclaw.data import ClawPlotData
        plotdata = ClawPlotData()


    plotdata.clearfigures()  # clear any old figures,axes,items data

    plotdata.format = 'ascii'

    # 2D pcolor plot if desired:
    # --------------------------

    plotfigure = plotdata.new_plotfigure(name='2D', figno=0)
    #plotfigure.show = False

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = 'auto'
    plotaxes.title = 'q(6) indicator'
    plotaxes.scaled = True      # so aspect ratio is 1

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = 5
    plotitem.pcolor_cmap = colormaps.blue_yellow_red
    plotitem.add_colorbar = True
    plotitem.show = True       # show on plot?
    plotitem.pcolor_cmin = 0.
    plotitem.pcolor_cmax = 1.0
    #plotitem.amr_data_show = [1,0]  # which levels to plot data
    plotitem.amr_patchedges_show = [1,1,1]
    #plotitem.amr_celledges_show = [1,0,0]



    # Figure for surface plot
    # -----------------------

    plotfigure = plotdata.new_plotfigure(name='surface', figno=2)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    #plotaxes.xlimits = [0,100e3]
    plotaxes.ylimits = [-300, 200]
    plotaxes.title = 'Surface plot'
    plotaxes.grid = True

    # Set up for item on these axes: scatter of 2d data
    plotitem = plotaxes.new_plotitem(plot_type='1d_from_2d_data')

    def surface(current_data):
        # Return radius of each grid cell and surface elevation
        h0 = 4000.  # ocean depth
        r = current_data.x[:,0]
        y = current_data.y
        ylow = y.min()
        dy = y[0,1] - y[0,0]
        q = current_data.q
        zfa = q[5,:,:]  # q(6,:,:) in Fortran, called zfa in Keh-Ming's code
        # compute have as in Keh-Mings out1.f:
        have = (1-zfa).sum(axis=1) * dy   # sum up in y
        for j in range(len(r)):
            if (zfa[j,:].min() < 1e-6) and (zfa[j,:].max() > 0.9):
                # surface seems to exist in jth column of this patch
                assert have[j] > 0, '*** expected have[j]>0'
                have[j] = have[j] + (ylow+h0)  # add depth of water below patch
                have[j] = have[j] - h0
            else:
                have[j] = 0.
        #import pdb; pdb.set_trace()
        if have.max() != 0:
            print(f'have.max() = {have.max()}')
        return r,have

    plotitem.map_2d_to_1d = surface
    #plotitem.plot_var = 0
    plotitem.plotstyle = '-'
    plotitem.color = 'b'
    #plotitem.amr_data_show = [1,1]  # which levels to plot data


    #-----------------------------------------
    # Figures for gauges
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='q', figno=300, \
                    type='each_gauge')
    plotfigure.show = False
    plotfigure.clf_each_gauge = True

    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = 'auto'
    plotaxes.title = 'Pressure'
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = 0
    plotitem.plotstyle = 'b-'


    # Parameters used only when creating html and/or latex hardcopy
    # e.g., via clawpack.visclaw.frametools.printframes:

    plotdata.printfigs = True                # print figures
    plotdata.print_format = 'png'            # file format
    plotdata.print_framenos = 'all'          # list of frames to print
    plotdata.print_fignos = 'all'            # list of figures to print
    plotdata.html = True                     # create html files of plots?
    plotdata.html_homelink = '../README.html'   # pointer for top of index
    plotdata.latex = True                    # create latex file of plots?
    plotdata.latex_figsperline = 2           # layout of plots
    plotdata.latex_framesperline = 1         # layout of plots
    plotdata.latex_makepdf = False           # also run pdflatex?
    plotdata.parallel = True                 # make multiple frame png's at once

    return plotdata
