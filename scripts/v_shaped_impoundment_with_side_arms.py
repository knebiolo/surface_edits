# -*- coding: utf-8 -*-
"""
Created on Tue Feb  1 08:49:31 2022

@author: KNebiolo

script intent: make 2d impoundment centerline a 3d enforced shapefile
"""
#%% Part 1 import modules, import data, set parameters, calculate centerline elevation function
# import modules
import geopandas
import pandas as pd
import fiona
import shapely
import os
import numpy as np
from scipy.interpolate import interp1d


# declare workspace
inputWS = r"J:\534\039\GIS\Langdale_Bahtymetry\Python\langdale_upstream\bathy_v2\Data"
outputWS = r"J:\534\039\GIS\Langdale_Bahtymetry\Python\langdale_upstream\bathy_v2\Output"

# script parameters
dam_crest_elev = 0.01
depth_at_0 = 3.
depth_at_n = 3.7
thalweg_len = 1495.
side_arm_0_len = 1246.
side_arm_1_len = 1007.
bank_angle = 89.

# import data
# centerline_segmented = geopandas.read_file(os.path.join(inputWS,'centerline_segmented.shp'))
centerline_route = geopandas.read_file(os.path.join(inputWS,'centerline.shp'))
side_arm_route = geopandas.read_file(os.path.join(inputWS,'side_arm_dissolve.shp'))
# arm_1 = geopandas.read_file(os.path.join(inputWS,'arm_1.shp'))
# arm_2 = geopandas.read_file(os.path.join(inputWS,'arm_2.shp'))
# arm_3 = geopandas.read_file(os.path.join(inputWS,'arm_3.shp'))
# arms = {1:arm_1,2:arm_2,3:arm_3}
bankline = geopandas.read_file(os.path.join(inputWS,'bankline.shp'))
bankline['z_enabled_geom'] = np.empty(len(bankline), dtype = 'object')

# create a piecewise linear function for centerline
center_knots = [0,thalweg_len]
center_elevs = [dam_crest_elev - depth_at_0, dam_crest_elev- depth_at_n]
center_elev_f = interp1d(center_knots,center_elevs, bounds_error = False, fill_value = 'extrapolate')

# extract route
center_route = centerline_route.geometry.iloc[0]

# create dictionary to hold side arm elevation functions
side_elev_fs = {}
side_routes = {}

for arm in side_arm_route.iterrows():
    # get the feature geometry
    feat_geom = list(arm[1]['geometry'].coords)

    # extract the first point - this intersects centerline
    p0 = shapely.geometry.Point(np.array([feat_geom[0][0],feat_geom[0][1]]))

    # get position of intersection along center point
    p0_dist = center_route.project(shapely.geometry.Point(p0))
    elev = center_elev_f(p0_dist)

    # create piecewise linear function for side arm centerline
    side_knots = [0,arm[1]['Shape_Leng']]
    side_elev = [elev,dam_crest_elev - depth_at_n]
    side_elev_fs[arm[0]] = interp1d(side_knots,side_elev, bounds_error = False, fill_value = 'extrapolate')

    # add geometry to side route
    side_routes[arm[0]] = arm[1]['geometry']

#%% Part 2 Add elevations to Bank lines and create data frames

# create coordinate lists
centerline_list =[]
side_arms = {}
bankline_list =[]

# fix elevations to the Z coordinate
for i in bankline.iterrows():
    shore = i[1]['geometry']
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

    # write to list of banklines
    bankline_list.append(line)

#%% Part 3 Iterate over vertices of centerlines and interpolate elevation
for coord in list(centerline_route.geometry[0].coords):

    # create point from coordinate pair
    pt = shapely.geometry.Point(coord)

    # now we worry about our Z elevation
    pt_dist = center_route.project(pt)
    elev = center_elev_f(pt_dist)

    # 3d point
    
    pt_xyz = list(list(pt.coords)[0])
    pt_xyz[2] = elev
    #pt_xyz = np.append(np.array(list(pt.coords)),elev)

    # create a new Point
    pt_xyz = shapely.geometry.Point(pt_xyz)

    # append to centerline list of points
    centerline_list.append(pt_xyz)
    

for arm in side_arm_route.iterrows():
    # get the feature geometry and ID
    arm_id = arm[0]
    arm_geom = arm[1]['geometry']
    side_arms[arm_id] = []

    for coord in arm_geom.coords:
        # create point from coordinate pair
        pt = shapely.geometry.Point(coord)

        # now we worry about our Z elevation
        pt_dist = arm_geom.project(pt)
        elev = side_elev_fs[arm_id](pt_dist)

        # write elevation to pt
        pt_list = list(list(pt.coords)[0])
        pt_list[2] = elev

        pt_xyz = shapely.geometry.Point(pt_list)

        # append to centerline list of points
        side_arms[arm_id].append(pt_xyz)



#%% Part 4 Convert and Export

# create geopandas dataframe for fixed shapefiles
centerline_line = shapely.geometry.LineString(centerline_list)
centerline_z_enabled = geopandas.GeoDataFrame([centerline_line], columns = ['geometry'])
bankline_z_enabled = geopandas.GeoDataFrame(bankline_list, columns = ['geometry'])
side_arms_z_enabled = geopandas.GeoDataFrame()

for arm in side_arms:
    side_arm_line = shapely.geometry.LineString(side_arms[arm])
    newRowArr = [0,side_arm_line]
    rowDF = geopandas.GeoDataFrame([newRowArr], columns = ['ID','geometry'])
    side_arms_z_enabled = side_arms_z_enabled.append(rowDF)

centerline_z_enabled.to_file(os.path.join(outputWS,'centerline.shp'))
bankline_z_enabled.to_file(os.path.join(outputWS,'bankline.shp'))
side_arms_z_enabled.to_file(os.path.join(outputWS,'side_arms.shp'))


print ("Script Complete Check Results")





