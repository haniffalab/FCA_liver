#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  7 11:42:51 2018

@author: doru
"""


import sys
args = sys.argv
save_to   = args[1]
expression_data_fname = args[2]
no_of_categories = int(args[3])

import pandas as pd
import numpy as np

data = pd.read_csv(expression_data_fname, index_col = None)

# convert Colours to r, g, b values, then to floats < 1.0
def hexdec_to_1floats(hexdec):
    return np.array([int(hexdec[1:][i:(i+2)], 16) for i in (0, 2, 4)]) / 255.0
    
gene_names = [gene_name for gene_name in data.columns[(2 + 2 * no_of_categories):]]
raw_expression = data.values[:, (2 + 2 * no_of_categories):]
gene_options = []
gene_expression_colour_coded = []
max_expression = raw_expression.max(axis = 1)
raw_expression / max_expression.reshape(max_expression.shape[0], 1)
max_expression_string = []
for index, gene_name in enumerate(gene_names):
    gene_expression = raw_expression[:, index]
    gene_expression = [str(value)[:min(4, len(str(value)))] for value in gene_expression]
    gene_expression = ",".join(gene_expression)
    gene_expression_colour_coded.append("gene_expression['{gn}'] = [{ge}]".format(gn = gene_name, ge = gene_expression))
    gene_options.append("<option value='{gn}'>{gn}</option>".format(gn = gene_name))
    max_expression_string.append("max_expression['{gene}'] = {val}".format(gene = gene_name, val = max_expression[index]))
gene_options = "".join(gene_options)  
gene_expression_colour_coded = ";".join(gene_expression_colour_coded)
max_expression_string = ";".join(max_expression_string)

# make coordinates data string
coordinates = data.values[:, 0:2].astype('float32')
# next few steps are compressing the data into a stadard cube centered at (0,0,0) and L = 200 
Xrange = np.percentile(coordinates[:, 0],  q = [1, 98]) * 1.2
Yrange = np.percentile(coordinates[:, 1],  q = [1, 98]) * 1.2
center = np.array((np.mean(Xrange), np.mean(Yrange)))
coordinates = coordinates - np.tile(center, (coordinates.shape[0], 1))
ratio = max(np.abs(np.percentile(coordinates[:, 0],  q = [1, 98]) * 1.2))
ratio = max(ratio, max(np.abs(np.percentile(coordinates[:, 1],  q = [1, 98]) * 1.2)))
ratio = 1.0 / ratio
coordinates = coordinates * ratio
coordinates = ",".join([str(value)[:min(6, len(str(value)))] for value in coordinates.ravel()])

categories = [str(value).replace(".", " ") for value in data.columns[2:(2 + no_of_categories)]]
categories_options = ["<option value='{cat}'>{cat}</option>".format(cat=cat) for cat in categories]
categories_options = "".join(categories_options)

categories_colours = []
categories_indices = []
for cat_index in range(no_of_categories):
    category_name = data.columns[2 + cat_index]
    category_name = category_name.replace(".", " ")
    category_colours = [hexdec_to_1floats(colour) for colour in data.values[:, 2 + cat_index + no_of_categories]]
    category_colours = [",".join([str(value)[:min(4, len(str(value)))] for value in colour]) for colour in category_colours]
    category_colours = ",".join(category_colours)
    categories_colours.append("categories_colours['{cn}'] = [{cc}]".format(cn = category_name, cc = category_colours))
    types = [value for value in np.unique(data.values[:, 2 + cat_index])]
    cat_indices = []
    categories_indices.append("categories_indices['{cn}'] = []".format(cn = category_name))
    for t_name in types:
        indices = data.values[:, 2 + cat_index] == t_name
        indices = np.where(indices)[0]
        indices = ",".join([str(value) for value in indices])
        cat_indices.append("categories_indices['{cn}']['{tn}'] = [{ind}]".format(cn = category_name, tn = t_name, ind = indices))
    cat_indices = "\n".join(cat_indices)
    categories_indices.append(cat_indices)
categories_indices = "\n".join(categories_indices)
categories_colours = "\n".join(categories_colours)

gene_families_file = open("./gene_families.txt", "r")
gene_families = gene_families_file.read()
gene_families_file.close()
geneFams = [fam.split("=")[0] for fam in gene_families.split("\n") if fam != ""]
geneFams = [fam.split("\'")[1] for fam in geneFams]
geneFams = ["<option value='{cat}'>{cat}</option>".format(cat=cat) for cat in geneFams]
geneFams = "".join(geneFams)

f = open('template.html', "r")
template_str = f.read()
f.close()
        
template_str = template_str.replace('gene_options_here', gene_options)
template_str = template_str.replace('gene_expression_colour_coded', gene_expression_colour_coded)
template_str = template_str.replace('coordinates_data_here', coordinates)
template_str = template_str.replace('category_options_here', categories_options)
template_str = template_str.replace('categories_colours_data_here', categories_colours)
template_str = template_str.replace('categories_indices_data_here', categories_indices)
template_str = template_str.replace('gene_families_options_here', gene_families)
template_str = template_str.replace('feature_family_option_here', geneFams)
template_str = template_str.replace('max_expression_here', max_expression_string)
with open(save_to, 'w') as result:
    result.write(template_str)

     