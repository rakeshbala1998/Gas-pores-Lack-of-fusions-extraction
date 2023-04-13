# -*- coding: utf-8 -*-
"""
Created on Sat Mar 18 18:20:27 2023

@author: rakes
"""
import sys
import os
import glob
import numpy as np
import trimesh
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
import networkx as nx
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns 
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Set the directory path
pore_path = "C:/Users/rakes/Dropbox (ASU)/Rakesh AM FEA/pores"
full_path="C:/Users/rakes/Dropbox (ASU)/Rakesh AM FEA/full plots"
# Define an empty list to store the features of each component
features_list = []
for file_name in os.listdir(full_path):  
    if file_name.endswith('.stl'):
# Read STL file and create mesh object
        mesh = trimesh.load(os.path.join(pore_path,file_name))
        #mesh=trimesh.load("C:\\Users\\rakes\\Dropbox (ASU)\\Rakesh AM FEA\\pores\\amga6.stl")
        #print(mesh.bounds)
        z_min=mesh.bounds[0,2];
        z_max=mesh.bounds[1,2];  
        height=z_max-z_min;
        if height < 10.5:
            print(f"{file_name}: STL file loaded is in correct size")
        else:
            print(f"{file_name}: STL file should be resized")
            continue
        # Extract submeshes based on volume
        submeshes = mesh.split(only_watertight=False)
        
        # Remove file extension from file_name
        sample_name = os.path.splitext(file_name)[0]
        #sample_name='amg6a'
        # Save each connected component as a separate STL file
        # Read STL file and create mesh object
        mesh1 = trimesh.load(os.path.join(full_path,file_name))
        #mesh1=trimesh.load("C:\\Users\\rakes\\Dropbox (ASU)\\Rakesh AM FEA\\full plots\\amga6.stl")
        
        print(mesh1.bounds)
        z_min=mesh1.bounds[0,2];
        z_max=mesh1.bounds[1,2];
        
        height1=z_max-z_min;
        
        if height1 < 10.5:
            print(f"{file_name}: STL file loaded is in correct size")
               #print(" STL file loaded is in correct size")
        else:
            print(f"{file_name}: STL file should be resized")
            continue
        # Convert the mesh to a graph
        graph = nx.from_edgelist(mesh.edges_unique)
        
        # Find the connected components
        components = list(nx.connected_components(graph))
        
        # Count the number of components
        num_components = len(components)
        
        print("Number of connectivity components :",num_components)
        # Extract connected components
        full_components = mesh1.split(only_watertight=False)
        
        largest_component = None
        largest_num_vertices = 0
        num_connected_components = len(full_components)
        
        for component in full_components:
            num_vertices = component.vertices.shape[0]
            if num_vertices > largest_num_vertices:
                largest_component = component
                largest_num_vertices = num_vertices
        
        #print("Number of connected components:", num_connected_components)
        #print("Number of vertices in largest component:", largest_num_vertices)
        for i,  submesh in enumerate(submeshes):
           
        
            # Calculate the volume of the component
            component_volume = submesh.volume
            
            #Calculate the equivalent diameter
            r=(((3/4)*(1/np.pi)*component_volume)**(1/3))*2*1000
            
            # Calculate the surface area of the component
            component_area = submesh.area
            
            # Calculate the centroid of the component
            component_centroid = submesh.centroid   
            
            # Calculate the angle from (0,0) with respect to positive x axis in counter clockwise direction
            angle=np.round(np.degrees(np.arctan2(component_centroid[1], component_centroid[0])),0)
            
            if angle<0:
                
                angle=360+angle

            # Find the closest point to the centroid
            valid_vertices = largest_component.vertices
            distances = np.sqrt(np.sum((valid_vertices - component_centroid)**2, axis=1))
            distance = np.min(distances)
            closest_point_idx = np.argmin(distances)
            closest_point = valid_vertices[closest_point_idx]
            
            # Find the unit vector between centroid and closest point
            diff_vector = closest_point - component_centroid
            unit_vector = diff_vector / np.linalg.norm(diff_vector)
            
            #distance from centroid 
            direction = np.array([0, 0, component_centroid[2]]) - component_centroid
            distance_from_sample_centroid = np.linalg.norm(direction)

        
            # Calculate the sphericity of the component
            component_sphericity = (np.pi**(1/3) * (6*component_volume)**(2/3)) / component_area
               
            # Append the features to the features_list
            features_list.append([sample_name,component_volume, r,component_area, component_centroid[0], component_centroid[1], component_centroid[2], angle, component_sphericity,distance, unit_vector, distance_from_sample_centroid ])
            print(features_list)
# Convert the features list to a pandas dataframe
columns = ['Sample name','Volume','equivalent dia(um)','Surface Area', 'Centroid X', 'Centroid Y', 'Centroid Z', 'angle','Sphericity','surface distance', 'unit vector','distance from sample centroid']
features_df = pd.DataFrame(features_list, columns=columns)
filtered_df = features_df.loc[features_df['distance from sample centroid'] > 1.25]
folder_path = "C:/Users/rakes/Dropbox (ASU)/Rakesh AM FEA/Data extraction"
file_name = 'pores_3d_NEW.csv'
file_path = os.path.join(folder_path, file_name)
#filtered_df.to_csv(file_path, index=True)
filtered_list = filtered_df.values.tolist()




# Create a new figure for each sample
fig1 = plt.figure(figsize=(10, 30))
ax1 = fig1.add_subplot(111, projection='3d')

# Set the aspect ratio of the plot box
ax1.set_box_aspect([10,10,30])

# Plot the centroids and direction arrows for each sample
for feature in filtered_list:
    centroid = feature[4:7]
    distance = feature[-3]
    direction = feature[-2]
    ax1.scatter(centroid[0], centroid[1], centroid[2], c='r', marker='o')
    ax1.quiver(centroid[0], centroid[1], centroid[2], direction[0]*distance, direction[1]*distance, direction[2]*distance, color='black')
    ax1.set_xlabel('X Label')
    ax1.set_ylabel('Y Label')
    ax1.set_zlabel('Z Label')
    
 