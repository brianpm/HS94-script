import math 
import time
import sys
import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt
import pandas as pd
import xarray as xr
import geocat.datafiles as gdf
import held_suarez_plt as hs94               


xr.set_options(keep_attrs=True)

# ====================
# Diagnostic flag
# ====================

verbose = False

# =======================================================
# Set up time management for our analysis
# =======================================================

ntpd = 4         # number of sample times per day (proscribed by HS94 test)

ntstart = 0      # sim time level start (day)
ntstop = 150     # sim time level end (day)
ntwin_start = 60 # sample time level start (time level)
ntwin_stop = 150 # sample time level stop (time level)

# =======================================================
# Initialize the Pressure (mbars) array on an equally-spaced
# set of sigma coordinates. These are the pressures we will
# interpolate to.
# =======================================================

nlev = 32
ntgbl = (ntstop-ntstart)

# =======================================================
# Set pressure levels: note these assume CAM ordering
# i.e. top of atmosphere (TOA) is level 1.
# =======================================================

dsigma = 1.0/nlev
lev = np.arange(0, nlev, 1) # CAM top of atmosphere lies on level 1.                                                                                                                 
sigma_mid = dsigma*(lev+0.5)
Psig = 1000.*sigma_mid

print("pressure on Sigma midpoints:")
print(Psig)

# =========================================================
# We need to pre-build the time-series analysis receptacles
# for this workflow.

#
# The analysis variables (var_gbl), where var=
# * Mean Zonal wind (u),
# * Northward eddy heat flux (vp*tp),
# * Eddy temperature variance(tp*tp), and
# * Northward eddy momentum flux(up*vp)
#
# They will have the shape: (global_time, levels, latitudes)
# ==========================================================

# ===============================================================
# First we need to open one of the datasets to get certain values
# needed by the time-series analysis receptacles
# ===============================================================


# ===============================================================
# path to test data on my (RDL) local laptop:
# ===============================================================

### 
### rho: Density
### Allow for getting rho out of the DataSet (.nc file) or creating it ab initio for testing purposes
### rho_is_present = True means get from file
### rho_is_present = False means create rho and just fill it with Zeros.
###

rho_is_present = False

if rho_is_present:

    # ===============================================================
    # Path on GLADE to real HS94 data:
    # ===============================================================

    hpath = "/glade/scratch/gdicker/val.FHS94.mpasa120.che.gnu/run/convertedOutputs_latlon/"
    
else:
    
    # ===============================================================
    # Use test data on local system
    # ===============================================================
    
    hpath = "/Users/loft/Desktop/Globus-MBP/"
    
hstem = "latlon_val.FHS94.mpasa120.che.gnu.cam.h1."
hextension= "0001-01-01-00000.nc"
hfile = hpath + hstem + hextension
ds = xr.open_dataset(hfile,engine="netcdf4")
ut=ds.U
lat=ut.coords["latitude"]
print(lat)
ds.close()

# Set the axis order of locally created arrays to be either 'F' (FORTRAN) or 'C' (C)

axis_order = 'C'

ds_gbl = xr.Dataset({
    'u_gbl': xr.DataArray(
        data   = np.zeros([ntgbl,32,360],dtype=float,order=axis_order),
        dims   = ['time', 'pressure', 'latitude'],
        coords = {'pressure': Psig, 'latitude': ut.latitude},
        attrs  = {'long_name': "mean Zonal wind (time series)", 'units': "m/s"}
    ),
    
    'vptp_gbl': xr.DataArray(
        data   = np.zeros([ntgbl,32,360],dtype=float,order=axis_order),
        dims   = ['time', 'pressure', 'latitude'],
        coords = {'pressure': Psig, 'latitude': ut.latitude},
        attrs  = {'long_name': "Northward eddy heat flux", 'units': "K m/s"}
    ),

    'tptp_gbl': xr.DataArray(
        data   = np.zeros([ntgbl,32,360],dtype=float,order=axis_order),
        dims   = ['time', 'pressure', 'latitude'],
        coords = {'pressure': Psig, 'latitude': ut.latitude},
        attrs  = {'long_name': "Eddy temperature_variance", 'units': "K^2"}
    ),

    'upvp_gbl': xr.DataArray(
        data   = np.zeros([ntgbl,32,360],dtype=float,order=axis_order),
        dims   = ['time', 'pressure', 'latitude'],
        coords = {'pressure': Psig, 'latitude': ut.latitude},
        attrs  = {'long_name': "Northward eddy momentum flux", 'units': "m^2/s^2"}
    )

    }
                      )

print("\n\n============U_GBL====================")
u_gbl = ds_gbl.u_gbl
print(u_gbl)

print("\n\n============VPTP_GBL=================")
vptp_gbl = ds_gbl.vptp_gbl
print(vptp_gbl)

print("\n\n============TPTP_GBL=================")
tptp_gbl = ds_gbl.tptp_gbl
print(tptp_gbl)

print("\n\n============U_GBL====================")
upvp_gbl = ds_gbl.upvp_gbl
print(upvp_gbl)

print("=====================================\n\n")

print("the size of time-series analysis receptacles in MB:", 4e-6*sys.getsizeof(u_gbl.data))
print("\n\nEntering the history file loop...")

###sys.exit()

#
# =========================================
# Loop over history files
# extracting and processing state variables 
# =========================================

ntime = 0  #starting time level 
nhfiles = 5 
cft=xr.cftime_range(start="2000", periods=nhfiles, freq="180H", calendar="365_day")

tic = time.perf_counter()
for hcount, htime in enumerate(cft):

    print("Begin processing history file #",hcount)
    
    # =========================================
    # Build the history file name
    # =========================================

    year = htime.strftime("%Y")
    year=year.replace("2000","0001")  # correct to year 1
    month = htime.strftime("%m")
    day = htime.strftime("%d")
    hrs = htime.strftime("%H")
    mins = htime.strftime("%M")
    sec = htime.strftime("%S")
    hrs=int(hrs)
    mins = int(mins)
    sec = int(sec)
    spd = 3600*hrs+60*mins+sec
    spd = str(spd)
    if spd == "0":
        spd = "00000"
    hextension= year + "-" + month + "-" + day + "-" + spd + ".nc"
    hfile = hpath + hstem + hextension
    print("opening hfile:", hfile)
    
    ds = xr.open_dataset(hfile,engine="netcdf4")
    if verbose:
        print(ds)

    # ======================================================
    # Determination whether the time level to start the time
    # averaging window has been reached
    # ======================================================

    nt=30 # NB: we should really pull this value from the dataset

    # =======================================================
    #
    #   Following translates lines 106-191 of FHS94_test.ncl
    #
    # =======================================================
    
    ut=ds.U
    vt=ds.V
    tt=ds.T
    if verbose:
        print("\n===================================")
        print("ut Data Array Contents")
        print(ut)
        print("===================================\n")
        print("ut dims = ", ut.dims)
        print("ut coords = ", ut.coords)
        print("ut atttrs = ", ut.attrs)
        print("===================================\n\n")


    if rho_is_present:
        rhot=ds.rho
    else:
        data = np.zeros((30,32,360,720),dtype=float,order=axis_order)
        rho_attrs = ut.attrs
        rhot = xr.DataArray(data, ut.coords, ut.dims, name="rho")
        rhot.attrs["mdims"] = "1"
        rhot.attrs["units"] = "kg/m^3"
        rhot.attrs["long_name"] = "density of dry air"

    if verbose:
        print("\n\n============RHO====================")
        print(rhot)
        print("===================================\n\n")

    ### 
    ### Create the pressure field: Pt 
    ### 
    ### Compute Pt from rhot using the ideal gas law.
    ### if rho_is_present = False Pt, will be zero everywhere.
    ###

    datap = np.zeros((30,32,360,720),dtype=float,order=axis_order)
    Pt_attrs = ut.attrs
    Pt = xr.DataArray(datap, ut.coords, ut.dims, name="P")
    Pt.attrs["mdims"] = "1"
    Pt.attrs["units"] = "J/m^3"
    Pt.attrs["long_name"] = "air pressure"
    if verbose:
        print("\n\n=============P=====================")    
        print(Pt)
        print("===================================\n\n")
    
    Rd = 287.0             # Gas constant in (J kg^-1 K^-1))
    P0 = 100000.0          # Reference pressure (kg m^-1 s^-2)
    Cmb = 1000.0/P0        # Conversion from MKS to millibars
    Pt = Cmb*Rd*rhot*tt    # Get pressure in mb via equation of state
    
    uprime = ut-ut.mean(dim="longitude")
    vprime = vt-vt.mean(dim="longitude")
    tprime = tt-tt.mean(dim="longitude")
        
    vptp = vprime*tprime
    tptp = tprime*tprime
    upvp = uprime*vprime

    # =======================================================
    # Interpolate quantities in the vertical from height
    # to equally-spaced pressure levels.
    # (var -> var_sig)
    # where var is one of {u,vptp,tptp,or upvp}
    # =======================================================

    
    ##=======================================================
    #
    # interpolate.interp1d fails with an error:
    # :
    # "ValueError: x and y arrays must be equal in length along interpolation axis."
    # Comment this out for now, and just copy variables * to  *_sig (11/09/22)
    #
    ##=======================================================
    
    #!! fint = interpolate.interp1d(Psig,ut,kind='linear', axis=2, copy=True, bounds_error=False, fill_value="extrapolate") 
    #!! u_sig    = fint(ut)
    #!! vptp_sig = fint(vptp)
    #!! tptp_sig = fint(tptp)
    #!! upvp_sig = fint(upvp)
    
    ut_sig = ut
    vptp_sig = vptp
    tptp_sig = tptp
    upvp_sig = upvp
    
    # =======================================================
    # Compute zonal means of eddy fluxes and U
    #
    # var_sig -> var_zm
    #
    # where var is one of {u,vptp,tptp,or upvp}
    # =======================================================

    u_zm    = ut_sig.mean(dim="longitude")
    vptp_zm = vptp_sig.mean(dim="longitude")
    tptp_zm = tptp_sig.mean(dim="longitude")
    upvp_zm = upvp_sig.mean(dim="longitude")

    # ===================================================================
    # These {*_zm} variables now have dimension (time, pressure, latitude)
    # Where time is the number of slices in a single history file (e.g. 30)
    # These values now need to be dumped into the global timeseries variables:
    #
    # var_zm -> var_gbl
    # where var is one of {u,vptp,tptp,or upvp} with shape (time_gbl, pressure, latitude)
    #
    # for time averaging over the time analysis window (e.g. days 200-1200).
    # ===================================================================

    ##  ===================================================================
    #
    #!! This doesn't work as a subsetting strategy, but yields the error:
    #!! "ValueError: shape mismatch: value array of shape (30,32,360) could not be broadcast to indexing result of shape (2,32,360)"
    #!! https://docs.xarray.dev/en/stable/user-guide/indexing.html
    #
    ##  ===================================================================
    
##    print("before shape error: nt=",nt," ntime=",ntime)
##    ind_time = xr.DataArray([ntime, ntime+nt-1], dims=["Time"])
    u_gbl[ntime:ntime+nt] = u_zm
    vptp_gbl[ntime:ntime+nt] = vptp_zm
    tptp_gbl[ntime:ntime+nt] = tptp_zm
    upvp_gbl[ntime:ntime+nt] = upvp_zm
        
    ntime+=nt
    print("time level pointer now set to =",ntime)
        
    # =======================================================
    # Close dataset and continue on with the file loop
    # =======================================================
    ds.close()

# =======================================================
# Grab out the days to sample (e.g. 200 to 1200)
# and take the time average
# =======================================================

### ind_gtime = xr.DataArray([ntwin_start, ntwin_stop], dims=["Time"])
time_samp = u_gbl[ntwin_start:ntwin_stop]
u_mn = time_samp.mean(dim="time")
u_mn.assign_coords({"latitude": ut.latitude})
u_mn.assign_coords({"pressure": Psig})
print(u_mn)

time_samp = vptp_gbl[ntwin_start:ntwin_stop]
vptp_mn = time_samp.mean(dim="time")
print(vptp_mn)

time_samp = tptp_gbl[ntwin_start:ntwin_stop]
tptp_mn = time_samp.mean(dim="time")
print(tptp_mn)

time_samp = upvp_gbl[ntwin_start:ntwin_stop]
upvp_mn = time_samp.mean(dim="time")
print(upvp_mn)

# ===========================================================
# These var_mn arrays should have dimensions (nlev, latitude)
# Let's see if this is true...
# ===========================================================

if verbose:
    print("Printing array description for u_mn:")
    print(u_mn)

print("Finished processing files...")
toc = time.perf_counter()
print(f"Time to process files {toc - tic:0.4f} seconds\n")

# =======================================================
# Develop contour plotting capability
# Borrowed from: https://matplotlib.org/stable/gallery/images_contours_and_fields/contour_demo.html
# =======================================================

lat=u_mn.coords["latitude"]
press = u_mn.coords["pressure"]

X, Y = np.meshgrid(lat, press)
scalex = 30.0
scaley = 250.0
scalez = 25.

Z1 = np.exp(-(X/scalex)**2 - (Y/scaley -2.0)**2)
Z2 = np.exp(-((X/scalex) - 1)**2 - ((Y/scaley-2.0) - 1)**2)
Z = scalez*(Z1 - Z2) * 2

hs94.fourpanel(X,Y,Z, u_mn, vptp_mn, tptp_mn, upvp_mn, display=False)

sys.exit()

# ================== DONE FOR THE DAY===============================

