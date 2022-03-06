# -*- coding: utf-8 -*-
"""
Created on Wed Apr  1 16:33:09 2020

@author: Kevin Nebiolo

Script Intent: parse hecras geometry file - pull out state/elevation relationship
by cross section (RM) and normalize for import to SQL
"""

# import modules
import os
import pandas as pd
import numpy as np
from shapely.geometry import LineString
from shapely.geometry import Point
from scipy import interpolate
import fiona

# identify workspace
geomWS = r"J:\3452\001\Calcs\Python\Data\Lower\geom"
shapeWS = r"J:\3452\001\Calcs\Python\Data\Lower\shape"
outputWS = r"J:\3452\001\Calcs\Python\Output"

# find geometry files
geoms = os.listdir(geomWS)

# create XS:xs_geometry dictionary
xs_geom = {}

# create XS:sta/elev dataframe dictionary
xs_sta_elev_dict = {}

# create XS:sta/elev interpolator dictionary
xs_sta_elev_int = {}

for  i in geoms:
    fileName = os.path.join(geomWS,i)
    o_file = open(fileName)
    counter = 0                                                                 # start line counter
    for line in o_file:
        print (line[:-1])

        # we have thalweg geometry
        if line[0:8] == "Reach XY":
            no_verts = int(line.split('=')[1])                                 # get the number of cross section vertices
            nrows = np.ceil(no_verts/2)                                       # calculate the number of rows, we have two XY pairs per row
            xs_geom_import = pd.read_fwf(os.path.join(geomWS,i),
                                         header = None,
                                         names = ['x1','y1','x2','y2'],
                                         colspecs = [(0,16),(16,32),(32,48),(48,64)],
                                         skiprows = counter + 1,
                                         nrows = nrows)                        # import cross section geometry
            x1y1 = xs_geom_import[['x1','y1']]
            x1y1.dropna(inplace = True)
            x1y1['id'] = np.arange(0,no_verts,2)
            x2y2 = xs_geom_import[['x2','y2']]
            x2y2.dropna(inplace = True)
            x2y2['id'] = np.arange(1,no_verts,2)
            x2y2.rename(columns = {'x2':'x1','y2':'y1'},inplace = True)
            xy = x1y1.append(x2y2)
            xy.sort_values(by = 'id', inplace = True)
            xy_tups = list(tuple(zip(xy.x1.values,xy.y1.values)))
            thalweg = LineString(xy_tups)                              # add XS geometry to dictionary as Shapely LineString
            thalweg_XY = xy
            del no_verts, nrows, xs_geom_import, x1y1, x2y2, xy, xy_tups

        # we have a new cross section - I assume that this always come before its attribute data - falls apart if this doesn't hold up
        if line[0:21] == "Type RM Length L Ch R":
            parts = line.split(',')
            new_XS = parts[1]

        # we have cross section geometry
        if line[0:15] == "XS GIS Cut Line":
            no_verts = int(line.split('=')[1])                                 # get the number of cross section vertices
            nrows = np.ceil(no_verts/2)                                       # calculate the number of rows, we have two XY pairs per row
            xs_geom_import = pd.read_fwf(os.path.join(geomWS,i),
                                         header = None,
                                         names = ['x1','y1','x2','y2'],
                                         colspecs = [(0,16),(16,32),(32,48),(48,64)],
                                         skiprows = counter + 1,
                                         nrows = nrows)                        # import cross section geometry
            x1y1 = xs_geom_import[['x1','y1']]
            x1y1.dropna(inplace = True)
            x1y1['id'] = np.arange(0,no_verts,2)
            x2y2 = xs_geom_import[['x2','y2']]
            x2y2.dropna(inplace = True)
            x2y2['id'] = np.arange(1,no_verts,2)
            x2y2.rename(columns = {'x2':'x1','y2':'y1'},inplace = True)
            xy = x1y1.append(x2y2)
            xy.sort_values(by = 'id', inplace = True)
            xy_tups = list(tuple(zip(xy.x1.values,xy.y1.values)))
            xs_geom[new_XS] = LineString(xy_tups)                              # add XS geometry to dictionary as Shapely LineString
            del no_verts, nrows, xs_geom_import, x1y1, x2y2, xy, xy_tups

        # we have a cross section station/elevation table
        if line[0:9] == '#Sta/Elev':
            no_verts = int(line.split('=')[1])
            nrows = np.ceil(no_verts/5)
            xs_sta_elev = pd.read_fwf(os.path.join(geomWS,i),
                                      header = None,
                                      names = ['sta1','elev1','sta2','elev2','sta3','elev3','sta4','elev4','sta5','elev5'],
                                      colspecs = [(0,8),(8,16),(16,24),(24,32),(32,40),(40,48),(48,56),(56,64),(64,72),(72,80)],
                                      skiprows = counter + 1,
                                      nrows = nrows)                        # import cross section geometry
            sta_elev1 = xs_sta_elev[['sta1','elev1']]
            sta_elev1.dropna(inplace = True)
            sta_elev1['id'] = np.arange(0,no_verts,5)
            sta_elev2 = xs_sta_elev[['sta2','elev2']]
            sta_elev2.rename(columns = {'sta2':'sta1','elev2':'elev1'},inplace = True)
            sta_elev2.dropna(inplace = True)
            sta_elev2['id'] = np.arange(1,no_verts,5)
            sta_elev3 = xs_sta_elev[['sta3','elev3']]
            sta_elev3.rename(columns = {'sta3':'sta1','elev3':'elev1'},inplace = True)
            sta_elev3.dropna(inplace = True)
            sta_elev3['id'] = np.arange(2,no_verts,5)
            sta_elev4 = xs_sta_elev[['sta4','elev4']]
            sta_elev4.rename(columns = {'sta4':'sta1','elev4':'elev1'},inplace = True)
            sta_elev4.dropna(inplace = True)
            sta_elev4['id'] = np.arange(3,no_verts,5)
            sta_elev5 = xs_sta_elev[['sta5','elev5']]
            sta_elev5.rename(columns = {'sta5':'sta1','elev5':'elev1'},inplace = True)
            sta_elev5.dropna(inplace = True)
            sta_elev5['id'] = np.arange(4,no_verts,5)
            sta_elev = sta_elev1.append(sta_elev2)
            sta_elev = sta_elev.append(sta_elev3)
            sta_elev = sta_elev.append(sta_elev4)
            sta_elev = sta_elev.append(sta_elev5)
            sta_elev.sort_values(by = 'id', inplace = True)
            sta_elev.reset_index(inplace = True)
            # if there are negative stages - shift things over - it's like these are getting fucked up
            min_sta = sta_elev.sta1.min()

            if min_sta < 0:
                sta_elev['sta1'] = sta_elev.sta1 + np.abs(min_sta)
            # create an interpolator for elevation as a function of stage along the xs
            int_f = interpolate.interp1d(sta_elev.sta1.values,sta_elev.elev1.values,fill_value = 'extrapolate')

            xs_sta_elev_dict[new_XS] = sta_elev
            xs_sta_elev_int[new_XS] = int_f
            del no_verts, nrows, xs_sta_elev, sta_elev1, sta_elev2, sta_elev4, sta_elev5, sta_elev

        counter = counter + 1
    o_file.close()

print ("Geometry file parsed - now perform linear referencing to get list of XYZ coordinates per XS")

xs_points = pd.DataFrame(columns = ['XS','X','Y','Z','sta'])

thalweg_dist_df = pd.DataFrame(columns = ['XS','X','Y','Z','dist'])

# for every cross section, perform linear referencing to get XYZ coordinates along XS
for key in xs_geom.keys():
    print ("Start analyzing XS: %s"%(key))
    # get geometry and station elevation relationship for this cross section
    geom = xs_geom[key]
    sta_elev = xs_sta_elev_dict[key]

    # for every row in the sta/elev - perform interpolation to get XY - then append to dataframe
    for row in sta_elev.iterrows():
        dist = row[1]['sta1']
        elev = row[1]['elev1']
        xy = list(geom.interpolate(dist).coords)
        print ("Station %s is at %s, %s, %s"%(dist,xy[0][0],xy[0][1],elev))
        new_row = np.array([key,xy[0][0],xy[0][1],elev,dist])
        new_row = pd.DataFrame(np.array([new_row]),columns = ['XS','X','Y','Z','sta'])
        xs_points = xs_points.append(new_row)

    # identify where the cross sections intersects the thalweg, ('XS',x, y, z, dist_on_thalweg) and write to dataframe
    intersection = thalweg.intersection(geom)                                  # get position of intersection
    dist_along_thalweg = thalweg.project(intersection)                         # get distance along thalweg
    dist_along_xs = geom.project(intersection)                                 # get distance along cross section
    elev = xs_sta_elev_int[key](dist_along_xs)      # get elevation at distance along cross section - remember, we created that dictionary of interpolator functions?
    print ("The thalweg intersects cross section %s at %s, %s at an elevation of %s"%(key,intersection.coords[0][0],intersection.coords[0][1],elev))
    new_row = np.array(np.array([key,intersection.coords[0][0],intersection.coords[0][1],elev,dist_along_thalweg]))
    new_row = pd.DataFrame(np.array([new_row]), columns = ['XS','X','Y','Z','dist'])
    thalweg_dist_df = thalweg_dist_df.append(new_row)



# create an interpolator for elevation as a function for distance along the stream with knots at each cross section
dist_elev_int = interpolate.interp1d(sta_elev.sta1.values,sta_elev.elev1.values,fill_value = 'extrapolate')

# now let's iterate over every point in the thalweg, project it's distance, and evaluate elevation
thalweg_XY['z'] = np.zeros(len(thalweg_XY))
thalweg_XY.set_index('id')
for row in thalweg_XY.iterrows():
    # get the ID and xy coordinates
    point_id = row[1]['id']
    x = row[1]['x1']
    y = row[1]['y1']
    # create a vertex
    vert = Point(x,y)
    # project the distance
    dist_along_thalweg = thalweg.project(vert)
    # get elevation from our interpolator function
    elev = dist_elev_int(dist_along_thalweg)
    # add elevation to the thalweg XY dataframe
    thalweg_XY.at[point_id,'z'] = elev


# now that we have finished with our geometry file, let's read vector data (banklines)
bank_dist_df = pd.DataFrame(columns = ['bank','xs','X','Y','Z','dist'])
shapes = os.listdir(shapeWS)
for shape in shapes:
    banklines = fiona.open(os.path.join(shapeWS,shape),'r')

    # for each bankline in the banklines shapefile
    for bankline in banklines:
        # create shapely linestring
        bank = LineString(bankline['geometry']['coordinates'])
        bank_id = bankline['id']
        print ("Start analyzing bank")
        # for every cross section, perform linear referencing to get XYZ coordinates along XS
        for key in xs_geom.keys():
            print ("Start analyzing XS: %s"%(key))
            # get geometry and station elevation relationship for this cross section
            geom = xs_geom[key]

            # identify where the cross sections intersects the bank, ('XS',x, y, z, dist_on_thalweg) and write to dataframe
            intersection = bank.intersection(geom)                                  # get position of intersection
            dist_along_bank = bank.project(intersection)                         # get distance along thalweg
            dist_along_xs = geom.project(intersection)                                 # get distance along cross section
            elev = xs_sta_elev_int[key](dist_along_xs)                                 # get elevation at distance along cross section - remember, we created that dictionary of interpolator functions?
            print ("The thalweg intersects cross section %s at %s, %s at an elevation of %s"%(key,intersection.coords[0][0],intersection.coords[0][1],elev))
            new_row = np.array(np.array([bank_id,key,intersection.coords[0][0],intersection.coords[0][1],elev,dist_along_bank]))
            new_row = pd.DataFrame(np.array([new_row]), columns = ['bank','xs','X','Y','Z','dist'])
            bank_dist_df = thalweg_dist_df.append(new_row)


## create an interpolator for elevation as a function for distance along the stream with knots at each cross section
#bank_elev_int = {}
#for
#bank_elev_int = interpolate.interp1d(sta_elev.sta1.values,sta_elev.elev1.values,fill_value = 'extrapolate')
#
#
#
#thalweg_XY['z'] = np.zeros(len(thalweg_XY))
#thalweg_XY.set_index('id')
#for row in thalweg_XY.iterrows():
#    # get the ID and xy coordinates
#    point_id = row[1]['id']
#    x = row[1]['x1']
#    y = row[1]['y1']
#    # create a vertex
#    vert = Point(x,y)
#    # project the distance
#    dist_along_thalweg = thalweg.project(vert)
#    # get elevation from our interpolator function
#    elev = dist_elev_int(dist_along_thalweg)
#    # add elevation to the thalweg XY dataframe
#    thalweg_XY.at[point_id,'z'] = elev
#
#xs_points.to_csv(os.path.join(outputWS,'xs_points.csv'))
#thalweg_XY.to_csv(os.path.join(outputWS,'thalweg_points.csv'))
#print ('XYZ file by cross section exported')
#
#
#
#





























