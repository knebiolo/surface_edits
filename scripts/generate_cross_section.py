# -*- coding: utf-8 -*-
"""
Created on Sun Feb  7 21:17:43 2021

@author: Kevin Nebiolo

Script Intent: Autogenerate points for cross sections along centerline of Wells River,VT
"""

# import modules
import geopandas
import pandas as pd
import fiona
import shapely
import os
import numpy as np
from scipy.interpolate import interp1d

# declare workspace
inputWS = r"J:\534\039\GIS\Langdale_Bahtymetry\Python\langdale_upstream\bathy_v2\Data\Langdale_CrossSection"
outputWS = r"J:\534\039\GIS\Langdale_Bahtymetry\Python\langdale_upstream\bathy_v2\Output\Langdale_CrossSection"
dam_crest_elev = 0.
bank_angle = -75.
thalweg_buffer = 900.

# import data
centerline = geopandas.read_file(os.path.join(inputWS,'centerline_segmented.shp'))
route = geopandas.read_file(os.path.join(inputWS,'centerline.shp'))
bankline = geopandas.read_file(os.path.join(inputWS,'bankline.shp'))
bankline['z_enabled_geom'] = np.empty(len(bankline), dtype = 'object')

# fix elevations to the Z coordinate
for i in bankline.iterrows():
    shore = i[1]['geometry']
    #z = i[1]['Elevation']
    z = dam_crest_elev

    # get list of shoreline geometry - we need to edit in the Z coordinate
    coords = list(shore.coords)

    # for every shoreline coordinate
    idx = 0
    for j in coords:
        # write the elevation to the z coordinate
        coord_list = list(coords[idx])
        coord_list[2] = z
        coords[idx] = tuple(coord_list)
        idx = idx + 1

    # create new linestring
    line = shapely.geometry.LineString(coords)

    # overwrite geometry
    bankline.loc[i[0],'z_enabled_geom'] = line

bankline = bankline.set_geometry('z_enabled_geom')

# create a piecewise linear function of elevation as a function of river mile
knots = pd.read_csv(os.path.join(inputWS,'ref_dep.csv'))
elev_f = interp1d(knots.distance.values,knots.depth.values)

# extract route
route = route.geometry.iloc[0]

# create geopandas dataframe of cross sections
xs_l_gdf = geopandas.GeoDataFrame(columns = ['Distance','geometry'])
xs_r_gdf = geopandas.GeoDataFrame(columns = ['Distance','geometry'])

bank_l_buff_gdf = geopandas.GeoDataFrame(columns = ['ID','geometry'])
bank_r_buff_gdf = geopandas.GeoDataFrame(columns = ['ID','geometry'])

# top of bank
bank_r_tob = []
bank_l_tob = []

# bottom of bank
bank_r_bob = []
bank_l_bob = []

xs_l_list = []
xs_r_list = []

# for a centerline segment
for feature in centerline.iterrows():
    feat_geom = list(feature[1]['geometry'].coords)                             # get the feature geometry
    p0 = np.array([feat_geom[0][0],feat_geom[0][1]])                                        # extract the first xy
    p1 = np.array([feat_geom[1][0],feat_geom[1][1]])                                        # extract the second xy

    p = np.array([p0[0],p0[1]]) - np.array([p1[0],p1[1]])                       # first vertex relative to second
    phat = p/np.linalg.norm(p)                                                  # unit vector

    # create rotation matrices
    rot_l = np.array([[np.cos(np.radians(-90)),
                       -np.sin(np.radians(-90))],
                      [np.sin(np.radians(-90)),
                       np.cos(np.radians(-90))]])
    rot_r = np.array([[np.cos(np.radians(90)),
                       -np.sin(np.radians(90))],
                      [np.sin(np.radians(90)),
                       np.cos(np.radians(90))]])

    # rotate that unit vector
    phat_l = rot_l.dot(phat)
    phat_r = rot_r.dot(phat)

    # calculate that point far off into space in the direction of the cross section that is guranteed to intersect a bankline
    p_l = p0 + thalweg_buffer * phat_l
    p_r = p0 + thalweg_buffer * phat_r

    # create shapely line strings
    xs_l = shapely.geometry.LineString([p0,p_l])
    xs_r = shapely.geometry.LineString([p0,p_r])

    # if we are in the impoundment
    for i in bankline.iterrows():
        #shore = shapely.geometry.LineString(np.array(i[1]['geometry'].coords)[:,:2])
        shore = i[1]['z_enabled_geom']

        if xs_l.intersects(shore):
            # calculate where the intersect, this is now the endpoint
            bank_l = shore.intersection(xs_l)

            try:
                len(bank_l)

            except:
                # calculate length of cross section
                xl_len = shapely.geometry.Point(p0).distance(bank_l)

                # now we worry about our Z elevation
                p0_dist = route.project(shapely.geometry.Point(p0))
                elev = elev_f(p0_dist)

                # convert elevation to depth
                depth = dam_crest_elev - elev

                # how far in from the shore does the bank start to slope up?
                buff_dist = np.tan(np.radians(bank_angle)) * depth
                if buff_dist > xl_len:
                    buff_dist = xl_len

                # find that point along the cross section
                buff_l = shapely.geometry.Point(np.array(list(bank_l.coords)[0][:2]) + buff_dist * phat_r)

                # now make everything 3d
                p0_xyz = np.insert(p0,2,elev)

                # convert to array and insert elevation
                buff_l_arr = np.array(list(buff_l.coords))
                buff_l_arr = np.insert(buff_l_arr,2,elev)

                # create a new XS consisting of these 3 points
                xs_l = shapely.geometry.LineString([p0_xyz,buff_l_arr,list(list(bank_l.coords)[0])])
                #xs_l = shapely.geometry.LineString([list(list(bank_l.coords)[0]),buff_l_arr, p0_xyz])

                del p0_xyz

                # add cross section to output dataframe
                xs_l_list.append(xs_l)

                # add cross section to output dataframe
                # bank_l_bob.append(list(bank_l.coords)[0])
                # bank_l_tob.append(buff_l_arr)
                bank_l_bob.append(buff_l_arr)
                bank_l_tob.append(list(bank_l.coords)[0])

        if xs_r.intersects(shore):
            # calculate where the intersect, this is now the endpoint
            bank_r = shore.intersection(xs_r)

            try:
                len(bank_r)
            except:

                # calculate length of cross section
                xr_len = shapely.geometry.Point(p0).distance(bank_r)

                # now we worry about our Z elevation
                p0_dist = route.project(shapely.geometry.Point(p0))
                elev = elev_f(p0_dist)

                # convert elevation to depth
                depth = dam_crest_elev - elev

                # calculate distance along cross section where bankline starts to slope up
                buff_dist = np.tan(np.radians(bank_angle)) * depth
                if buff_dist > xr_len:
                    buff_dist = xr_len

                # find that point along the cross section
                buff_r = shapely.geometry.Point(np.array(list(bank_r.coords)[0][:2]) + buff_dist * phat_l)

                # now make everything 3d
                p0_xyz = np.insert(p0,2,elev)

                # convert to array and insert elevation
                buff_r_arr = np.array(list(buff_r.coords))
                buff_r_arr = np.insert(buff_r_arr,2,elev)

                # create a new XS consisting of these 3 points
                xs_r = shapely.geometry.LineString([p0_xyz,buff_r_arr,list(list(bank_r.coords)[0])])
                #xs_r = shapely.geometry.LineString([list(list(bank_r.coords)[0]),buff_r_arr,p0_xyz])

                del p0_xyz

                # add cross section to output dataframe
                xs_r_list.append(xs_r)

                # add cross section to output dataframe
                bank_r_tob.append(list(bank_r.coords)[0])
                bank_r_bob.append(buff_r_arr)

# find and remove overlapping cross sections
xs_l_exp = xs_l_list.copy()
rm_list = []
for i in np.arange(0,len(xs_l_list)-2):
    curr_xs = xs_l_list[i]
    next_xs = xs_l_list[i+1]
    if next_xs.intersects(curr_xs):
        rm_list.append(i+1)

xs_l_exp = np.delete(xs_l_exp,rm_list)

xs_r_exp = xs_r_list.copy()
rm_list = []
for i in np.arange(1,len(xs_r_list)-1):
    curr_xs = xs_r_list[i]
    prev_xs = xs_r_list[i-1]
    if curr_xs.intersects(prev_xs):
        rm_list.append(i-1)

xs_r_exp = np.delete(xs_r_exp,rm_list)

# convert to geo pandas dataframe
xs_r_gdf = geopandas.GeoDataFrame(xs_r_exp, columns = ['geometry'])
xs_l_gdf = geopandas.GeoDataFrame(xs_l_exp, columns = ['geometry'])

# save to shapefile
xs_l_gdf.to_file(os.path.join(outputWS,'xs_l.shp'))
xs_r_gdf.to_file(os.path.join(outputWS,'xs_r.shp'))


newRowArr = [1,shapely.geometry.LineString(bank_l_tob)]
imp_l_buff_gdf = geopandas.GeoDataFrame([newRowArr], columns = ['ID','geometry'])

newRowArr = [1,shapely.geometry.LineString(bank_r_tob)]
imp_r_buff_gdf = geopandas.GeoDataFrame([newRowArr], columns = ['ID','geometry'])

imp_l_buff_gdf.to_file(os.path.join(outputWS,'top_of_bank_l.shp'))
imp_r_buff_gdf.to_file(os.path.join(outputWS,'top_of_bank_r.shp'))

newRowArr = [1,shapely.geometry.LineString(bank_l_bob)]
imp_l_buff_gdf = geopandas.GeoDataFrame([newRowArr], columns = ['ID','geometry'])

newRowArr = [1,shapely.geometry.LineString(bank_r_bob)]
imp_r_buff_gdf = geopandas.GeoDataFrame([newRowArr], columns = ['ID','geometry'])

imp_l_buff_gdf.to_file(os.path.join(outputWS,'bot_of_bank_l.shp'))
imp_r_buff_gdf.to_file(os.path.join(outputWS,'bot_of_bank_r.shp'))

new_bank = geopandas.GeoDataFrame(bankline, columns = ['z_enabled_geom'])
new_bank.to_file(os.path.join(outputWS,'new_bankline.shp'))

print ("Script Complete Check Results")





