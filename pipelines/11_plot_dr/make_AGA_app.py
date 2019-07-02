#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 22 11:03:12 2018

@author: doru
"""

# argument variables
import sys
output_folder = sys.argv[1]

from os.path import join

# file names
material_folder      = join(output_folder,   "AGA_folder")
save_to              = join(output_folder,   'AGAlinkage_map_{cat}.html'.format(cat = sys.argv[2]))
colors_fname         = join(material_folder, 'colours.csv')
connectivities_fname = join(material_folder, 'connectivities.csv')
coordinates_fname    = join(material_folder, 'coordinates.csv')

# read data from files in csv formatr
import pandas as pd
connectivities = pd.read_csv(connectivities_fname, index_col = 0, header = 0)
coordinates = pd.read_csv(coordinates_fname, index_col = 0, header = 0)
try:
    colors = pd.read_csv(colors_fname, index_col = 0, header = 0)
except FileNotFoundError:
    cell_types = connectivities.columns
    import random
    cell_types = [f for f in connectivities.columns]
    colours = []
    for cell_type in cell_types:
        r   = lambda: random.randint(0,255)
        col = '#%02X%02X%02X' % (r(),r(),r())
        colours.append({'CellTypes': cell_type, 'Colours': col})
    colors = pd.DataFrame(colours)
    colors = colors.set_index('CellTypes')

scaleScale = 1.4
minX = coordinates.min()[0] * scaleScale
minY = coordinates.min()[1] * scaleScale
maxX = coordinates.max()[0] * scaleScale
maxY = coordinates.max()[1] * scaleScale

# prepare the coordinates and colors data
cell_names = list(coordinates.index)
cell_sizes = coordinates.Size.tolist()
# reorder cell names by population size - so during drawing smaller cell population are not covered by bigger bubbles
cell_names = [cell_name for [cell_size, cell_name] in sorted(zip(cell_sizes, cell_names), reverse = True)]
data_coordinates = []
for cell_name in cell_names:
    row_data = coordinates.loc[cell_name]
    X, Y, R  = row_data.X, row_data.Y, row_data.Size
    X        = (X - minX) / (maxX - minX);
    Y        = (Y - minY) / (maxY - minY);
    color    = colors.loc[cell_name].Colours
    indata   = 'data_coordinates["{cell_name}"] = [{X}, {Y}, {R}, "{C}"]'.format(cell_name = cell_name, 
                                X = X, Y = Y, R = R, C = color)
    data_coordinates.append(indata)
data_coordinates = '\n'.join(data_coordinates)

# prepare edge thickness data
data_edges = []
# rearrange connectivities by order of cell name
for cell_name in cell_names:
    indata = connectivities[cell_name][cell_names].tolist()
    indata = ','.join([str(i) for i in indata])
    indata = 'data_edges["{cell_name}"] = [{indata}]'.format(cell_name = cell_name, indata = indata)
    data_edges.append(indata)
data_edges = '\n'.join(data_edges)

# make cell_names array
cell_names = ['"{cell_name}"'.format(cell_name = cell_name) for cell_name in cell_names]
cell_names = ','.join(cell_names)
cell_names = 'cell_names = [{cell_names}]'.format(cell_names = cell_names)

# prepare all the data
data = '\n'.join([data_coordinates, data_edges, cell_names])

template_fobj = open('template_for_AGA_app.html', 'r')
template = template_fobj.read();
template_fobj.close()

# insert data in template
template = template.replace('// insert data here', data)

# save interactive page
with open(save_to, 'w') as save_fobj:
    save_fobj.write(template)