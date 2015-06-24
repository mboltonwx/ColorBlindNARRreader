from netCDF4 import Dataset
import sys
from pylab import *
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.axes_grid1 import make_axes_locatable
from datetime import datetime
import scipy.ndimage
# Requires basemap, netcdf4-python, matplotlib, numpy, scipy
# can be obtained through "conda install <pkgname>"

############################################
#
# NARR Plotter
# Greg Blumberg (OU/CIMMS)
# wblumberg@ou.edu
#
# Modified by Matt Bolton (How The Weatherworks)
# matt.bolton@weatherworks.com
#
#
# This code will make some nice maps from the
# NARR. Good for getting an idea of what 
# happened on a certain day.
#
############################################

def regMap():
    '''
        Define map location. Code borrowed from https://github.com/keltonhalbert/AWIDS	
    '''
    figure(figsize=(10,8))
    m = Basemap(width=1500000,height=1100000,
                  rsphere=(6378137.00,6356752.3142),\
                  resolution='l',area_thresh=1000.,projection='lcc',\
                  lat_1=40,lat_2=30,lat_0=30,lon_0=-87)
    m.drawcoastlines()
    m.drawcountries()
    m.drawstates()
    m.drawcounties()
    return m

def e2td(E):
    '''
        Function to convert vapor pressure to dewpoint.  Needed
        for some of the maps that get plotted.
    
        I think the input units are in mb.
    '''
    B = (np.log(E / 6.108)) / 17.27
    D = (237.3 * B) / (1 - B)
    return D

def cc(temp):
    '''
        Calusius Clapyron equation to convert temperature or dewpoint
        to saturation vapor pressure or vapor pressure.
        
        I think the input units are in Kelvin or Celsius.
        Units of the (saturation) vapor pressure are in mb.
    '''
    e = 6.112 * np.exp((17.67*temp)/(temp + 243.5))
    return e

def q2w(w):
    '''
        I really don't remember what this does.
    
        But it's important.
    '''
    return w/(1.+w) 

def w2e(w, p):
    '''
        Now this converts mixing ratio to vapor pressure.
        w is the mixing ratio
        p is the pressure (mb)
    '''
    eps = 0.622
    e = ((w/eps)*p)/(1 + (w/eps))
    return e



# When running narr_plotter.py, you need to include some command line arguments.
# i.e. python narr_plotter.py 19990503 21 svr
# This example will make a severe weather type map for May 3rd, 1999 at 21 UTC
# Other types of maps can be (here are the arguments):
#   850
#   700
#   500
#   300
#   sfc
#   sfccnt

yyyymmdd = sys.argv[1]
hh = sys.argv[2]
type = sys.argv[3]

# Below modifies the URL to the NARR data using the date/time information supplied by the user.
narr_path = 'http://nomads.ncdc.noaa.gov/thredds/dodsC/narr-a/YYYYMM/YYYYMMDD/narr-a_221_YYYYMMDD_HH00_000.grb'
narr_path = narr_path.replace('YYYYMMDD', yyyymmdd)
narr_path = narr_path.replace('YYYYMM', yyyymmdd[:6])
narr_path = narr_path.replace('HH', hh)

# Print out the path created
print "Here is the path to the NARR data:", narr_path

# Load the data into Python from online using the Dataset() command.
d = Dataset(narr_path)
# d becomes a Python dictionary, where you pass it a key and it returns a value, just like a real dictionary.

#for i in d.variables.keys():

# Print out all the keys (this will show you all the variables the NARR has (i.e. CAPE, winds, temperature, etc.)
print d.variables.keys()

# Load in an extra file that has the latitude, longitude points of the NARR grid.  This isn't online,
# so I have a file contained within this package that has this called 'narr_lat_lon.nc'
ll = Dataset('narr_lat_lon.nc')
lat = ll.variables['lat'][:]
lon = ll.variables['lon'][:]
ll.close()


if type == 'sfc':
    # Make a surface map with color fills

    # Create a datetime object that holds the date and time of the map I'm plotting
    dt = datetime.strptime(yyyymmdd+hh, '%Y%m%d%H')

    # Make a string from the datetime object (this will go in the map title)
    dt_str = datetime.strftime(dt, '%d %b %Y  %H UTC')

    # Pull the different variables we want to plot out from the NARR database
    # Each variable has some dimensions I've had to look up, and it can be imagined like a cube.
    # For example, the dimensions could be east-west, north-south, up-down
    # That's what something like the [0,:,:] means
    mslp = d.variables['Mean_sea_level_pressure_ETA_model'][0,:,:]/100. # convert MSLP to mb from Pascals
    sfc_temp = d.variables['Temperature_height_above_ground'][0,0,:,:] # 2 meter temperature
    u_wind = d.variables['u_wind_height_above_ground'][0,0,:,:] * 1.94384 # 10 meter u wind converted to knots
    v_wind = d.variables['v_wind_height_above_ground'][0,0,:,:] * 1.94384 # 10 meter v wind converted to knots
    
    # Make a regular map and get the map object so we can draw the data on it.
    m = regMap()
    
    # Draw the title on the map.
    title(dt_str + ' ' + 'Surface NARR-A', fontsize=15)

    # Convert the latitude longitude NARR grid into x,y coordinates (Basemap handles this)
    x,y = m(lon, lat)

    # Draw the MSLP pressure contour on the map.  Use an interval of every 4 mb starting at 940 mb to 1104 mb.
    CS = m.contour(x,y, mslp, np.arange(940,1104,4), colors='k', linewidths=2)

    # Draw labels for those contours
    clabel(CS, CS.levels, fmt='%4.0f')

    # We don't want to plot every wind barb, so let's plot every 5.
    stride = 2

    # Convert the surface temperature from Kelvin to Farenheit
    sfc_temp =((sfc_temp - 273.15)*1.8 + 32)

    # Draw a filled contour, where the fill corresponds to the temperature.  Use the colormap "RdYlBu_r"
    # Here is where you'll probably need to play around with different colormaps.
    cb = m.contourf(x,y,sfc_temp, np.arange(-40,132,2), cmap=get_cmap("RdYlBu_r"))

    # Draw the freezing line on the map and label it.
    fz = m.contour(x,y,sfc_temp, np.asarray([32]), colors='m', linestyles='--', linewidths=2)
    clabel(fz, fz.levels, fmt='%4.0f')

    # Draw the surface wind barbs on the map.
    barbs(x[::stride,::stride],y[::stride,::stride],u_wind[::stride,::stride], v_wind[::stride,::stride])
    
    # Stuff to draw the colorbar and label it.
    divider = make_axes_locatable(gca())
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cb = colorbar(cb, cax=cax)
    cb.set_label("Temperature [F]")
    tight_layout()
    #show()

    # Save the map to the disk as a .png
    savefig(yyyymmdd + '.' + hh + '.sfc.png')

if type == 'sfccnt':
    # This will draw surface temperature lines, and a contour fill of dewpoint
    # This won't be commented as througouly as the above stuff.
    dt = datetime.strptime(yyyymmdd+hh, '%Y%m%d%H')
    dt_str = datetime.strftime(dt, '%d %b %Y  %H UTC')
    
    # Download the data we're going to plot.
    mslp = d.variables['Mean_sea_level_pressure_ETA_model'][0,:,:]/100.
    sfc_temp = d.variables['Temperature_height_above_ground'][0,0,:,:]
    u_wind = d.variables['u_wind_height_above_ground'][0,0,:,:] * 1.94384
    v_wind = d.variables['v_wind_height_above_ground'][0,0,:,:] * 1.94384
    sh = d.variables['Specific_humidity_height_above_ground'][0,0,:,:]
    pres = d.variables['Pressure'][0,0,:,:]/100.
    
    # Convert the specific humidity into dewpoint.
    e = w2e(q2w(sh), pres)
    dwpt = e2td(e)
    dwpt = dwpt*1.8 + 32.

    # Convert Sfc temp to F.
    sfc_temp =((sfc_temp - 273.15)*1.8 + 32)
    
    # Draw the map background.
    m = regMap()
    
    title(dt_str + ' ' + 'Surface NARR-A', fontsize=15)
    x,y = m(lon, lat)

    # Draw the dewpoint using a color fill of greens between 50 F to 82 F every 2 F.
    cb = m.contourf(x,y,dwpt, np.arange(50,82,2), cmap=get_cmap('Greens'))

    # Draw the MSLP lines and label them.
    CS = m.contour(x,y, mslp, np.arange(940,1104,4), colors='k', linewidths=2)
    clabel(CS, CS.levels, fmt='%4.0f')
    
    # Run a smoother through the surface temperature field so it looks a little cleaner on the map.
    sfc_temp = scipy.ndimage.gaussian_filter(sfc_temp, 1.5)

    # This is some funky code that will allow us to color the temperature contours based on whether or not
    # they are above, at, or below freezing by using red, purple, and blue colors (all respectively).
    temp_levels = np.arange(-60,150,10)
    zero_level = np.where(temp_levels == 0)[0]
    temp_colors = np.repeat('r', len(temp_levels))
    temp_colors[zero_level] = 'm'
    temp_colors[:zero_level] = 'b'

    # Plot the temperature lines and label them.
    tm = m.contour(x,y,sfc_temp, temp_levels, colors=temp_colors, linewidths=2, linestyles='--')
    clabel(tm, tm.levels, fmt='%4.0f')

    # Plot every 5 surface wind barb
    stride=2 
    barbs(x[::stride,::stride],y[::stride,::stride],u_wind[::stride,::stride], v_wind[::stride,::stride])
    
    # Draw and position the colorbar onto the figure.
    divider = make_axes_locatable(gca())
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cb = colorbar(cb, cax=cax)
    cb.set_label("Dewpoint [F]")
    tight_layout()
    #show()

    # Save the figure
    savefig(yyyymmdd + '.' + hh + '.sfc_dew.png')

if type == 'svr':
    m = regMap()
    mslp = d.variables['Mean_sea_level_pressure_ETA_model'][0,:,:]
    dt = datetime.strptime(yyyymmdd+hh, '%Y%m%d%H')
    dt_str = datetime.strftime(dt, '%d %b %Y  %H UTC')


    v_sfc = d.variables['v_wind_height_above_ground'][0,:,:] 
    u_sfc = d.variables['u_wind_height_above_ground'][0,:,:] 
    
    # Compute the lapse rate (this stuff isn't really used in the program yet.
    # 700 temperature:
    pres_idx = np.where(d.variables['isobaric'][:] == 700)[0]
    temp7 = d.variables['Temperature'][0,pres_idx, :,:]-273.15
    z7 = d.variables['Geopotential_height'][0, pres_idx,:,:]/1000.
    print z7, temp7

    # 500 temp:
    pres_idx = np.where(d.variables['isobaric'][:] == 500)[0]
    temp5 = d.variables['Temperature'][0,pres_idx, :,:]-273.15
    u5 = d.variables['u_wind'][0,pres_idx, :,:]
    v5 = d.variables['v_wind'][0,pres_idx, :,:]
    z5 = d.variables['Geopotential_height'][0, pres_idx,:,:]/1000.
    print z5, temp5
    lr = ((temp7 - temp5) / (z7 - z5))[0]
    # Finished doing the lapse rate stuff

    # Get the lowest 180 mb most unstable CAPE
    cape = d.variables['Convective_available_potential_energy'][0,:,:][0]

    # Find the surface to 500 mb wind shear (convert from m/s to knots)
    u_shear = 1.943 * (u5[0] - u_sfc)
    v_shear = 1.943 * (v5[0] - v_sfc)
    mag = np.sqrt(np.power(u_shear, 2) + np.power(v_shear, 2))
    
    # Mask all of the grid points where the shear is less than 30 knots.
    u_shear = np.ma.masked_where(mag < 30, u_shear)[0]
    v_shear = np.ma.masked_where(mag < 30, v_shear)[0]

    # Get the specific humidity and convert to dewpoint (F)
    sh = d.variables['Specific_humidity_height_above_ground'][0,0,:,:]
    pres = d.variables['Pressure'][0,0,:,:]/100.
    e = w2e(q2w(sh), pres)
    dwpt = e2td(e)
    x,y = m(lon, lat)
    dwpt = dwpt*1.8 + 32.

    # Smooth the dewpoint to make it a little nicer and cleaner
    sfc_dwpt = scipy.ndimage.gaussian_filter(dwpt, 1.5)

    # Draw the dewpoint lines and label them
    cb = m.contour(x,y,sfc_dwpt, np.arange(50,82,2), colors='g')
    clabel(cb, cb.levels, fmt='%4.0f') 

    # Draw the MSLP contours and label them.
    CS = m.contour(x,y, mslp, np.arange(940,1104,4), colors='k', linewidths=3)
    clabel(CS, CS.levels, fmt='%4.0f')
    
    title(dt_str + ' ' + 'Surface NARR-A', fontsize=15)
    
    # Plot every 5 shear vector on the map
    stride=2
    #lr_levels = -1 * np.asarray([5.5, 6, 6.5, 7, 7.5, 8, 8.5, 9, 9.5, 10])
    
    # Plot only CAPE values between 500 and 6500 J/kg at every 500 J/kg
    # Use the reversed spring color map ('spring_r')
    cape_levels = np.arange(500,6500,500)
    cb = m.contourf(x,y,cape, cape_levels, cmap='spring_r')

    # Plot the wind shear vectors
    barbs(x[::stride,::stride],y[::stride,::stride],u_shear[::stride,::stride], v_shear[::stride,::stride])
    
    # Draw the colorbar and position it where we want it.
    divider = make_axes_locatable(gca())
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cb = colorbar(cb, cax=cax)
    cb.set_label("CAPE [J/kg]")
    tight_layout()

    # Save the map as a .png
    savefig(yyyymmdd + '.' + hh + '.svr.png')
    show() # Not 100% needed, but this will show you the map

# At this point of the program, they may want to make a upper air map instead.
# The following code handles that...

def plotUA(level):
    # This is a specific function that runs to create specific upper air charts, like 850, 700, 500 etc.
    # the variable level gets passed to it, and level contains a pressure level (i.e. 850, 700)
    dt = datetime.strptime(yyyymmdd+hh, '%Y%m%d%H')
    dt_str = datetime.strftime(dt, '%d %b %Y  %H UTC')
    
    # Find the array index in the NARR dataset that corresponds to the pressure level requested
    pres_idx = np.where(d.variables['isobaric'][:] == level)[0]
    
    # Get the surface pressure
    sfc_press = d.variables['Pressure'][0,0, :,:]

    # Find all of the gridpoints where the pressure level requested is greater than the surface pressure
    # This finds array indices where the pressure surface is "below ground level."
    mask = np.where(sfc_press/100. < level)
    
    temp = d.variables['Temperature'][0,pres_idx, :,:]-273.15
    rh = d.variables['Specific_humidity'][0, pres_idx, :, :]
    z = d.variables['Geopotential_height'][0, pres_idx,:,:]
    u = d.variables['u_wind'][0,pres_idx, :,:]*1.94384
    v = d.variables['v_wind'][0, pres_idx, :,:]*1.94384
    
    # Convert the specific humidity to relative humidity (%)
    es = cc(temp)
    e = w2e(q2w(rh), level)
    rh = (e/es)*100.
    
    # This part is used to hide the portion of the pressure surface that is "below ground"
    # np.nan means it's not a number, which matplotlib won't plot
    temp[0][mask] = np.nan
    rh[0][mask] = np.nan
    z[0][mask] = np.nan
    u[0][mask] = np.nan
    v[0][mask] = np.nan
    terrain = np.ones(v[0].shape)
    terrain[mask] = 0.
    
    
    m = regMap()
    title(dt_str + ' ' + str(level) + ' mb NARR-A', fontsize=15)
    
    tight_layout()
    x,y = m(lon,lat)
    
    # Draw the terrain in grey ( this is the place where the pressure surface is "below ground level"
    contourf(x,y,terrain, [0,.5,1], colors=['k','#FFFFFF'], alpha=.2)
    
    # This if statement distinguishes what variables ought to be plotted given different pressure levels.
    if level > 500:
        # If the pressure level is below 500 mb, the contour fill should be relative humidity, and should be green
        cb = m.contourf(x,y,rh[0], np.arange(70,105,5), cmap=get_cmap('Greens'))
    else:
        # Instead, the wind speed is probably a more important variable.  Plot that and the wind barbs instead.
        wind_spd = np.sqrt(np.power(u[0],2) + np.power(v[0], 2)) # Pythagorean theorem to get wind speed
        cb = m.contourf(x,y,wind_spd, np.arange(60,240,10), cmap=get_cmap('Reds'), alpha=.8)
        m.barbs(x[::5,::5],y[::5,::5],u[0,::5,::5], v[0,::5,::5])
    
    #CS = m.contour(x,y,z[0], np.arange(5160-(60*10), 5880+(60*10),60), colors='k', linewidths=2)
    
    # Always plot the geopotential height 
    # Need to actually provide specific contour levels for different maps (haven't done that yet)
    CS = m.contour(x,y,z[0], colors='k', linewidths=2) # Plot the geopotential height
    clabel(CS, CS.levels, fmt='%4.0f') # 
    
    # Do the same coloring of different isotherms (temperature contours) 
    # to distinguish whether or not they are above, at, or below freezing
    temp_levels = np.arange(-80,82,2)
    zero_level = np.where(temp_levels == 0)[0]
    temp_colors = np.repeat('r', len(temp_levels)) # Use red (could also use hex colors here?)
    temp_colors[zero_level] = 'm' # use magenta
    temp_colors[:zero_level] = 'b' # use blue 
    CS = m.contour(x,y,temp[0], temp_levels, colors=temp_colors, linestyles='--', linewidths=1.5)
    clabel(CS, CS.levels, fmt='%4.0f')
    
    # Draw the colorbar and position it 
    divider = make_axes_locatable(gca())
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cb = colorbar(cb, cax=cax)
    cb.set_label("Relative Humidity [%]")

    # Save the figure
    savefig(yyyymmdd + '.' + hh + '.ua.' + str(level) + '.png') 
    #show()

print "TYPE OF UA PLOT:", type
try:
    # A really cheap hack to find out whether or not the user requested a upper air map
    type = int(type)
    # If the variable "type" can't be converted to an integer, it crashes the program and the program quits gracefully
    plotUA(type)
    d.close()
except Exception,e:
    print "You haven't selected a valid map type to plot.  Try again!"
    print e
    sys.exit() # Exit the program

    # some left over dummy code...don't remember why
    """
    m = regMap()
    sfc_temp = d.variables['Temperature'][0,0,:,:]
    u_wind = d.variables['u_wind_height_above_ground'][0,0,:,:] * 1.94384
    v_wind = d.variables['v_wind_height_above_ground'][0,0,:,:] * 1.94384
    x,y = m(lon, lat)
    CS = m.contour(x,y, mslp, np.arange(940,1104,4), colors='k', linewidths=2)
    clabel(CS, CS.levels, fmt='%4.0f')
    stride = 2
    sfc_temp =((sfc_temp - 273.15)*1.8 + 32)
    cb = m.contourf(x,y,sfc_temp, np.arange(-40,132,2), cmap=get_cmap("jet"))
    barbs(x[::stride,::stride],y[::stride,::stride],u_wind[::stride,::stride], v_wind[::stride,::stride])
    colorbar(cb)
    tight_layout()
    savefig(yyyymmdd + '.' + hh + '.sfc.png')
    """
    
