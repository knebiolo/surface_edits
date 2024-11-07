# -*- coding: utf-8 -*-
"""
Created on Tue Feb 22 19:13:56 2022

Script Intent: create a wireframe model of a river from 1D HECRAS Cross Sections

@author: KNebiolo
"""
# import modules
import numpy as np
import scipy.interpolate as interpolate
import os
import pandas as pd
import matplotlib.pyplot as plt
from shapely.geometry import LineString
from shapely.geometry import Point
import geopandas

# identify workspace
inputWS = r"C:\Users\knebiolo\OneDrive - Kleinschmidt Associates, Inc\Software\surface_edits\data"

# import data
dat = pd.read_csv(os.path.join(inputWS,'xs_points.csv'))

# create dictionary of bsplines for each cross section
xs_dict = {}
xs_list = dat.XS.unique()

wireframe = pd.DataFrame()

for xs in xs_list:
    # get this cross section's data
    xs_dat = dat[dat.XS == xs]

    # get bspline parameters for cross section
    t, c, k = interpolate.splrep(xs_dat.sta, xs_dat.Z, s= 0, k = 1)


    # identify the xmin and xmax and create a linespace
    xmin = xs_dat.sta.min()
    xmax = xs_dat.sta.max()
    N = 100
    xarr = np.linspace(xmin,xmax,N)

    # create a spline function and write to dictionary
    spline = interpolate.BSpline(t, c, k, extrapolate = False)
    xs_dict[xs] = spline

    # now calculate the X and Y values for our 100 stations across cross section
    xs_pts = tuple(zip(xs_dat.X,xs_dat.Y))
    xs_linestring = LineString(xs_pts)

    counter = 0
    for i in xarr:
        pt = xs_linestring.interpolate(i)
        x = list(pt.coords)[0][0]
        y = list(pt.coords)[0][1]
        z = spline(i)
        newRow = pd.DataFrame(np.array([[xs,i,counter,x,y,z]]), columns = ['XS','sta','id','x','y','z'])
        wireframe = wireframe.append(newRow)
        counter = counter + 1

    # plot and profit
    plt.scatter(xs_dat.sta,xs_dat.Z, label = "Original Points")
    plt.plot(xarr,spline(xarr),'r',label = 'Spline')
    plt.show()

wireframe.dropna(inplace = True)
# for each station across the cross section, create a longitudinal bSpline
long_dict = {}

for sta in np.arange(0,100,1):
    # get longitudinal data at this station
    long_sta_dat = wireframe[wireframe.id == sta]
    long_sta_dat.sort_values(by = 'XS', inplace = True)
    
    # get bspline parameters for x and y coordinates separately
    tx, cx, kx = interpolate.splrep(long_sta_dat.XS, long_sta_dat.x, s= 0, k = 3)
    ty, cy, ky = interpolate.splrep(long_sta_dat.XS, long_sta_dat.y, s= 0, k = 3)
    tz, cz, kz = interpolate.splrep(long_sta_dat.XS, long_sta_dat.z, s= 0, k = 3)


    # identify the xmin and xmax and create a linespace
    xmin = long_sta_dat.XS.min()
    xmax = long_sta_dat.XS.max()
    N = 100
    xarr = np.linspace(xmin,xmax,N)

    # create a longitudinal spline function and write to dictionary
    spline_x = interpolate.BSpline(tx, cx, kx, extrapolate = False)
    spline_y = interpolate.BSpline(ty, cy, ky, extrapolate = False)
    spline_z = interpolate.BSpline(tz, cz, kz, extrapolate = False)

    long_dict[sta] = (spline_x,spline_y,spline_z)

    # now calculate the X and Y values for our 100 stations across cross section
    long_pts = tuple(zip(long_sta_dat.x,long_sta_dat.y))
    long_linestring = LineString(long_pts)

    # plot and profit
    plt.scatter(long_sta_dat.x,long_sta_dat.y, label = "Original Points")
    plt.plot(spline_x(xarr),spline_y(xarr),'r',label = 'Spline')
    plt.show()

# now iterate over cross stations and stations and create a wireframe
xs_min = dat.XS.min()
xs_max = dat.XS.max()

# create empty xyz dataframe
xyz_pts = []

# create empty output arrays
X = np.array([])
Y = np.array([])
Z = np.array([])

# for every station calculate X, Y, and Z
for xs in np.linspace(xs_min,xs_max,1000):
    for sta in np.arange(0,100,1):
        splines = long_dict[sta]
        X = np.append(X,splines[0](xs))
        Y = np.append(Y,splines[1](xs))
        Z = np.append(Z,splines[2](xs))
        xyz_pts.append(Point(splines[0](xs),splines[1](xs),splines[2](xs)))
        
# convert to geo pandas dataframe
xyz_pts = geopandas.GeoDataFrame(xyz_pts, columns = ['geometry'])
xyz_pts.to_file(os.path.join(inputWS,'xyz.shp'))

# plot and profit   
fig = plt.figure()
ax = fig.add_subplot(projection='3d')

# Plot a basic wireframe.
ax.scatter(X, Y, Z)

plt.show()
    
    


