#
#
# Create 4 panel Held Suarez 94 test case plot
# Contours plots field Z vs X (lat) and Y (pressure)
#

import matplotlib.cm as cm
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
import matplotlib.gridspec as gridspec
import format_zonal_plt as zp

# suggested by: https://stackoverflow.com/questions/23876588/matplotlib-colorbar-in-each-subplot

from mpl_toolkits.axes_grid1 import make_axes_locatable 
### 
def fourpanel(X,Y,Z, u_mn, vptp_mn, tptp_mn, upvp_mn, display):

    test_data=False
    if test_data:
        lat=u_mn.coords["latitude"]
        print(lat)
        press = u_mn.coords["pressure"]
        print(press)
    else:
        lat=u_mn.coords["latitude"]
        press = u_mn.coords["pressure"]
#        Z=u_mn 
#        vptp_mn = Z
#        tptp_mn = Z
#        upvp_mn = Z

    #
    # Info about stacking figures was found here:
    # https://matplotlib.org/stable/gallery/subplots_axes_and_figures/subplots_demo.html
    #

    fig = plt.figure(figsize=(8,8))
    gs = fig.add_gridspec(nrows=2, ncols=2, hspace=0.2)
    axs = gs.subplots()

    #fig, axs = plt.subplots(2,2)
    
    #
    # Define font sizes for various parts of the plots:
    #
    
    fs_title=10
    fs_axlab=8
    fs_major=8
    fs_minor=6
    
    fig.suptitle("Held-Suarez 1994 Test Case Results", fontweight="bold")
    axs[0,0].set_title("Zonal wind                                     m/s", fontdict={'fontsize':fs_title})
    axs[0,1].set_title("Eddy tempurature variance                K$^{2}$", fontdict={'fontsize':fs_title})
    axs[1,0].set_title("Northward eddy momentum flux    m$^{2}$/s$^{2}$", fontdict={'fontsize':fs_title})
    axs[1,1].set_title("Northward eddy heat flux                K m/s", fontdict={'fontsize':fs_title})

    # 
    # Contour plotting function:
    # See documentation at https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.contour.html
    #
    pltno = 0
    for ax in axs.flat:

        if pltno == 0:
            print("\nplotting u_mn\n")
            tmp=u_mn
        elif pltno == 1:
            print("plotting tptp_mn\n")
            tmp=tptp_mn
        elif pltno == 2:
            print("plotting upvp_mn\n")
            tmp=upvp_mn
        elif pltno == 3:
            print("plotting vptp_mn\n")
            tmp=vptp_mn
        else:
            print("ERR: unknown HS94 field, or bad plotting number\n")
            sys.exit()
        pltno=pltno+1
        
        #
        # Tutorial on choices of color map:
        # https://matplotlib.org/stable/tutorials/colors/colormaps.html
        #
        
        CS = ax.contourf(lat, press, tmp, 20, cmap="RdBu_r", alpha=0.8)
        #CS = ax.contour(lat, press, tmp, 20, colors = "k")  # colors = "k" means negative contours default to dashed.
        ax.tick_params(axis='both', which='major', labelsize=fs_major)
        #ax.clabel(CS, fontsize=fs_minor, inline=True)
        ax.set_ylabel('pressure (mbars)', fontsize = fs_axlab)
        
        # x axis
        ax.xaxis.set_major_locator(MultipleLocator(30))
        ax.xaxis.set_major_formatter(plt.FuncFormatter(zp.NS_format_func))

        # For the minor ticks, use no labels; default NullFormatter.
        ax.xaxis.set_minor_locator(MultipleLocator(10))

        # y axis
        # Get the bottom of the atmosphere at the top of the plot!
        ax.invert_yaxis() 
        ax.yaxis.set_major_locator(MultipleLocator(200))
        ax.yaxis.set_major_formatter(plt.FuncFormatter(zp.millibar_format_func))
        #ax.yaxis.set_major_formatter('{x:.0f}')

        # For the minor ticks, use no labels; default NullFormatter.
        ax.yaxis.set_minor_locator(MultipleLocator(50))

        #
        # Add colorbar to each plot
        #
        
        divider = make_axes_locatable(ax)
        cax = divider.append_axes('bottom', size='5%', pad=0.25)
        fig.colorbar(CS, cax=cax, orientation='horizontal')
    
    if display:
        plt.show()
    else:
        plt.savefig('hs94_4panel.png',dpi=600)
    

