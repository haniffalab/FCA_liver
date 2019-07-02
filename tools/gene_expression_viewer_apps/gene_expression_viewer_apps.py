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

template_str = """
<!doctype html>

<html lang='en'>
<head>
  <meta charset='utf-8'>

  <title>3D viewer</title>
  <meta name='description' content='The HTML5 Herald'>
  <meta name='author' content='Dorin-Mirel Popescu'>

</head>

<body>
	<table>
		<tr>
			<td align='left'>
				<form>
					<fieldset>
						<legend><b>Visualisation options</b></legend>
						<label for = 'particleSizeBar'>Particle size: </label>
						<input type='range' name = 'particleSizeBar' min = 1 max = 14 step=0.1 oninput='setParticleSize(value)' value = 2 /><br />
				
						<label for = 'alphaInput'>Transparency: </label>
						<input type='range' name = 'alphaInput' min = 0 max = 1000 oninput='setAlpha(value)' value = 1000 /><br />
				
						<label for = 'canvasSizeInput'>Canvas size: </label>
						<input type='range' name = 'canvasSizeInput' min = 200 max = 2000 oninput='setCanvasSize(value)' value = 500 /><br />
				
						<label for = 'bgInput'>Dark background: </label>
						<input type='radio' name = 'bgInput' oninput='setBackground(value)' value = 'dark' />
						<label for = 'bgInput'>White background: </label>
						<input type='radio' name = 'bgInput' oninput='setBackground(value)' value = 'white' checked />
						<br />
					</fieldset>
				</form>
			</td>
			<td style='vertical-align: top' rowspan='2'>
			
				<form>
					<fieldset>
						<legend><b>Colour by:</b></legend>
							<table>
								<tr>
									<td>
										<label for='colourType'><input type='radio' name='colourType' onchange='setColourBy(value)' value='gene_expression' />Gene expression: </label>
									</td>
									<td>
										<label for='geneSelector'><select name='geneSelector' id='geneSelector' onchange='selectFeature()'>gene_options_here</select></label>
									</td>
								</tr>
								<tr>
									<td colspan = '2' align='center'>
										<canvas id='canvasColorScale' width = 200 height=40></canvas>
									</td>
								</tr>
								<tr>
									<td>
										<label for='colourType'><input type='radio' name='colourType' checked onchange='setColourBy(value)' value='category' />Category:</label>
									</td>
									<td>
										<label for='categorySelector'><select name='categorySelector' id='categorySelector' onchange = 'setCategory()'>category_options_here</select></label>
									</td>
								</tr>
							</table>
					</fieldset>
				</form>
				<br />
				<div>
					<fieldset>
						<legend><b>Cell types:</b></legend>
						<label for='toggleRadio'><input type='checkbox' name = 'toggleRadio' id='toggleRadio' onchange='toggleAllTypes()' checked />Show all:</label>
						<form id = 'typesControlPanel'>
							
						</form>
					</fieldset>
				</div>
			</td>
		</tr> 
		<tr>
			<td style='vertical-align: text-top' >
				<canvas id='canvas' width=600 height=600></canvas>
			</td>
		</tr>
	</table>
  <script id='vertex-shader' type='x-shader/x-fragment'>
  	attribute vec4 a_Position;
  	attribute vec3 a_Color;
  	uniform float u_basePointSize;
  	uniform float u_Alpha;
  	uniform int u_PaintFeatureScale;
  	varying vec4 v_Color;
  	void main() {
    	gl_Position = a_Position;
    	gl_PointSize = u_basePointSize;
    	if (u_PaintFeatureScale == 0){
    		v_Color = vec4(a_Color, u_Alpha);
    	}
    	else{
    		float r = 0.0;
    		float g = 0.0;
    		float b = 0.0;
    		r = max(0.0, 2.0 * a_Color.r - 1.0);
    		b = max(0.0, 2.0 * (1.0 - a_Color.r) - 1.0);
    		g = 1.0 - 2.0 * abs(a_Color.r - 0.5);
    		v_Color = vec4(r, g, b, u_Alpha);
    	}
  	}
  </script>
  <script id ='fragment-shader' type='x-shader/x-fragment'>
	precision mediump float;
  	varying vec4 v_Color;
  	void main() {
  		float r = 0.0;
  		vec2 cxy = 2.0 * gl_PointCoord - 1.0;
  		r = dot(cxy, cxy);
  		if (r > 1.0){
  			discard;
  		}
    	gl_FragColor = v_Color;
  	}
  </script>
  <script type = 'text/javascript'>
  	
  	var Matrix4 = function(opt_src) {
  		var i, s, d;
  		if (opt_src && typeof opt_src === 'object' && opt_src.hasOwnProperty('elements')) {
    	s = opt_src.elements;
    	d = new Float32Array(16);
    	for (i = 0; i < 16; ++i) {
      		d[i] = s[i];
    	}
    	this.elements = d;
  		} else {
    	this.elements = new Float32Array([1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1]);
  		}
	};
	
	Matrix4.prototype.setTranslate = function(x, y, z) {
  		var e = this.elements;
  		e[0] = 1;  e[4] = 0;  e[8]  = 0;  e[12] = x;
  		e[1] = 0;  e[5] = 1;  e[9]  = 0;  e[13] = y;
  		e[2] = 0;  e[6] = 0;  e[10] = 1;  e[14] = z;
  		e[3] = 0;  e[7] = 0;  e[11] = 0;  e[15] = 1;
  		return this;
	};
	
	Matrix4.prototype.setLookAt = function(eyeX, eyeY, eyeZ, centerX, centerY, centerZ, upX, upY, upZ) {
  		var e, fx, fy, fz, rlf, sx, sy, sz, rls, ux, uy, uz;

  		fx = centerX - eyeX;
  		fy = centerY - eyeY;
  		fz = centerZ - eyeZ;

  		// Normalize f.
  		rlf = 1 / Math.sqrt(fx*fx + fy*fy + fz*fz);
  		fx *= rlf;
  		fy *= rlf;
  		fz *= rlf;

  		// Calculate cross product of f and up.
  		sx = fy * upZ - fz * upY;
  		sy = fz * upX - fx * upZ;
  		sz = fx * upY - fy * upX;

  		// Normalize s.
  		rls = 1 / Math.sqrt(sx*sx + sy*sy + sz*sz);
  		sx *= rls;
  		sy *= rls;
  		sz *= rls;

  		// Calculate cross product of s and f.
  		ux = sy * fz - sz * fy;
  		uy = sz * fx - sx * fz;
  		uz = sx * fy - sy * fx;

  		// Set to this.
  		e = this.elements;
  		e[0] = sx;
  		e[1] = ux;
  		e[2] = -fx;
  		e[3] = 0;

  		e[4] = sy;
  		e[5] = uy;
  		e[6] = -fy;
  		e[7] = 0;

  		e[8] = sz;
  		e[9] = uz;
  		e[10] = -fz;
  		e[11] = 0;

  		e[12] = 0;
  		e[13] = 0;
  		e[14] = 0;
  		e[15] = 1;

  		// Translate.
  		return this.translate(-eyeX, -eyeY, -eyeZ);
	};
	
	Matrix4.prototype.translate = function(x, y, z) {
  		var e = this.elements;
  		e[12] += e[0] * x + e[4] * y + e[8]  * z;
  		e[13] += e[1] * x + e[5] * y + e[9]  * z;
  		e[14] += e[2] * x + e[6] * y + e[10] * z;
  		e[15] += e[3] * x + e[7] * y + e[11] * z;
  		return this;
	};
	
	Matrix4.prototype.setPerspective = function(fovy, aspect, near, far) {
  		var e, rd, s, ct;

  		if (near === far || aspect === 0) {
    		throw 'null frustum';
  		}
  		if (near <= 0) {
    		throw 'near <= 0';
  		}
  		if (far <= 0) {
    		throw 'far <= 0';
  		}

  		fovy = Math.PI * fovy / 180 / 2;
  		s = Math.sin(fovy);
  		if (s === 0) {
    		throw 'null frustum';
  		}

  		rd = 1 / (far - near);
  		ct = Math.cos(fovy) / s;

  		e = this.elements;

  		e[0]  = ct / aspect;
  		e[1]  = 0;
  		e[2]  = 0;
  		e[3]  = 0;

  		e[4]  = 0;
  		e[5]  = ct;
  		e[6]  = 0;
  		e[7]  = 0;

  		e[8]  = 0;
  		e[9]  = 0;
  		e[10] = -(far + near) * rd;
  		e[11] = -1;

  		e[12] = 0;
  		e[13] = 0;
  		e[14] = -2 * near * far * rd;
  		e[15] = 0;

  		return this;
	};
  </script>
  <script type='text/javascript'>
  
  	function buildCategoryRadioButtons(){
  		category_type = categorySelector.options[categorySelector.selectedIndex].value;
  		current_indices = indices_all;
  		// create radio commands from categories
  		typesControlPanel.innerHTML = "";
  		radio_commands_HTML = "";
  		for(name in categories_indices[category_type]){
  			f_index = categories_indices[category_type][name][0]
  			cols = categories_colours[category_type].slice(3 * f_index, 3 * f_index + 3)
  			col_label = "#";
  			for(k=0;k<cols.length;k++){col_hex = Math.round(255 * cols[k]).toString(16).padStart(2, '0'); col_label = col_label + col_hex}
  			radio_command = "<div style='background-color:" + col_label + "'>";
  			radio_command = radio_command + "<input style='float:left' type='checkbox' id='" + name;
  			radio_command = radio_command + "' checked onchange='toggleCategoryAction()' /><label style='float:left' for='" + name + "'";
  			radio_command = radio_command + ">" + name + ": </label><br/></div>"
  			radio_commands_HTML = radio_commands_HTML + radio_command
  		}
  		typesControlPanel.innerHTML = radio_commands_HTML;
  	}
  	
  	function toggleCategoryAction(){
  		updateBuffer()
  		draw()
  	}
  
  	function setCategory(){
  		buildCategoryRadioButtons()
  		updateBuffer()
  		draw()
  	}
  	
  	function setColourBy(value){
  		colour_by = value;
  		if (colour_by =='category'){
  			PaintFeatureScale = 0;
  		}else{
  			PaintFeatureScale = 1;
  		}
  		gl_context.uniform1i(u_PaintFeatureScale, PaintFeatureScale)
  		updateBuffer()
  		draw()
  	}
  	
  	
  	function toggleAllTypes(){
  		controlRadios = typesControlPanel.elements
  		for(i=0;i<controlRadios.length;i++){
  			controlRadios[i].checked = toggleRadio.checked
  		}
  		updateBuffer()
  		draw()
  	}
  	
  	function selectFeature(){
  		feature = geneSelector.value
  		updateBuffer()
  		draw()
  		console.log('selected features')
  	}
  	
  	function draw(){
  		if(bg_color == "white"){
			gl_context.clearColor(1, 1, 1, 1)
		}else{
			gl_context.clearColor(0, 0, 0, 1)
		}
		gl_context.clear(gl_context.COLOR_BUFFER_BIT);
  		gl_context.bufferData(gl_context.ARRAY_BUFFER, buffer_data_array, gl_context.STATIC_DRAW)
  		gl_context.drawArrays(gl_context.POINTS, 0, n)
  	}
  	
  	function updateBuffer(){
  		var buffer_data = [];
  		// first update indices to be used - for this read the category control panel radio buttons
  		controlRadios = typesControlPanel.elements
  		current_indices = []
  		for(i=0;i<controlRadios.length;i++){
  			if(controlRadios[i].checked){
  				radio_type = controlRadios[i].id
  				current_indices = current_indices.concat(categories_indices[category_type][radio_type])
  			}
  		}
  		// now just populate the buffer_data
  		if(colour_by == 'gene_expression'){
  			current_indices.forEach(function(index, i){
  				buffer_data.push(coordinates_data[2 * index])
  				buffer_data.push(coordinates_data[2 * index + 1])
  				buffer_data.push(gene_expression[feature][index])
  				buffer_data.push(gene_expression[feature][index])
  				buffer_data.push(gene_expression[feature][index])
  			})
  		}else{
  			current_indices.forEach(function(index, i){
  				buffer_data.push(coordinates_data[2 * index])
  				buffer_data.push(coordinates_data[2 * index + 1])
  				buffer_data.push(categories_colours[category_type][3 * index])
  				buffer_data.push(categories_colours[category_type][3 * index + 1])
  				buffer_data.push(categories_colours[category_type][3 * index + 2])	
  			})
  		}
  		buffer_data_array = new Float32Array(buffer_data)
  		n = buffer_data_array.length / 5
  	}
  	
  	function setParticleSize(value){
  		particleSize = parseInt(value)
  		gl_context.uniform1f(u_basePointSize, particleSize)
  		updateBuffer()
  		draw()
  	}
  	
  	function setAlpha(value){
  		alphaValue = parseInt(value) / 1000
  		gl_context.uniform1f(u_Alpha, alphaValue)
  		updateBuffer()
  		draw()
  	}
  	
  	function setCanvasSize(value){
  		value = parseInt(value)
  		canvas.width = value
  		canvas.height = value
  		gl_context = getContext(canvas)
		gl_context = initContext(gl_context)
		gl_context.viewport(0, 0, canvas.width, canvas.height)
		updateBuffer()
		draw()
  	}
  	
  	function setBackground(value){
  		bg_color = value;
  		draw()
  	}
  	
  	function shadersFromScriptElement(gl, ID, type){
  		shaderScript = document.getElementById(ID)
  		var str = ''
  		var k = shaderScript.firstChild;
  		while(k){
  			if (k.nodeType == 3){
  				str += k.textContent;
  			}
  			k = k.nextSibling
  		}
  		var shader = gl.createShader(type)
  		gl.shaderSource(shader, str)
  		gl.compileShader(shader)
  		return shader
  	}
  	
  	function getContext(canvasWidget){
  		var names = ['webgl', 'experimental-webgl', 'webkit-3d', 'moz-webgl'];
  		for(var i=0; i<names.length; i++){
  			try{
  				var gl = canvasWidget.getContext(names[i])
  			}catch(e){}
  			if(gl){i=names.length}
  		}
  	
  		var vshader = shadersFromScriptElement(gl, 'vertex-shader', gl.VERTEX_SHADER),
  			fshader = shadersFromScriptElement(gl, 'fragment-shader', gl.FRAGMENT_SHADER)
  			program = gl.createProgram();
  		gl.attachShader(program, vshader)
  		gl.attachShader(program, fshader)
  		gl.linkProgram(program)
  		gl.useProgram(program)
  		gl.program = program
  		return gl
  	}
  	
  	function initContext(gl){
  		n = buffer_data_array.length / 5
  		var vertexColourBuffer = gl.createBuffer()
  		gl.bindBuffer(gl.ARRAY_BUFFER, vertexColourBuffer)
  	
  		var FSIZE = buffer_data_array.BYTES_PER_ELEMENT;
  	
  		var a_Position = gl.getAttribLocation(gl.program, 'a_Position')
  		gl.vertexAttribPointer(a_Position, 2, gl.FLOAT, false, FSIZE * 5, 0)
  		gl.enableVertexAttribArray(a_Position)
  	
  		var a_Color = gl.getAttribLocation(gl.program, 'a_Color')
  		gl.vertexAttribPointer(a_Color, 3, gl.FLOAT, false, FSIZE * 5, 2 * FSIZE)
  		gl.enableVertexAttribArray(a_Color)
  	
  		u_basePointSize = gl.getUniformLocation(gl.program, 'u_basePointSize')
  		gl.uniform1f(u_basePointSize, particleSize)
  	
  		u_Alpha = gl.getUniformLocation(gl.program, "u_Alpha")
  		gl.uniform1f(u_Alpha, alphaValue)
  		
  		u_PaintFeatureScale = gl.getUniformLocation(gl.program, 'u_PaintFeatureScale')
  		gl.uniform1i(u_PaintFeatureScale, PaintFeatureScale)
  	
  		gl.clearColor(1, 1, 1, 1);
  		if(bg_color == "dark"){
  			gl.clearColor(0, 0, 0, 1)
  		}
  	
  		gl.disable(gl.DEPTH_TEST)
  		gl.enable(gl.BLEND)
  		gl.blendFunc(gl.SRC_ALPHA, gl.ONE_MINUS_SRC_ALPHA)
  	
  		gl.clear(gl.COLOR_BUFFER_BIT);
  		return gl
  	}
  
  	var categorySelector  = document.getElementById('categorySelector'),
  		geneSelector      = document.getElementById('geneSelector'),
  		typesControlPanel = document.getElementById('typesControlPanel'),
  		toggleRadio       = document.getElementById('toggleRadio')
  	
  	var canvas       = document.getElementById('canvas'),
  		particleSize = 5,
  		alphaValue   = 1.0,
  		bg_color     = "white",
  		n            = 0,
  		particleSize = 2,
  		PaintFeatureScale = 0;
  		
  	coordinates_data = [coordinates_data_here]
  	gene_expression  = []; gene_expression_colour_coded;
  	
  	categories_colours = []
  	categories_colours_data_here
  	categories_indices = []
  	categories_indices_data_here
  	
  	// initialize flags
  	// when toggling between gene expression and category, do not slice data i.e. do not recompute index data
  	// when choosing a category always re-initiate index data
  	
  	var colour_by      = 'category', // the other options is can be 'category'
  		category_types = [],
  		category_type  = '',
  		features       = [],
  		feature        = ''; 
  	
  	// set category
  	for(name in categories_colours){category_types.push(name)}
  	category_type = category_types[0]
  	
  	// set feature
  	for(name in gene_expression){features.push(name)}
  	feature = features[0];
  	
  	// create global data holders
  	var indices_all = [],
  		current_indices = [],
  		current_colours = [],
  		buffer_data_array = [];
  		
  	for(j=0;j<categories_colours[category_type].length/3;j++){indices_all.push(j)}
  	
  	// build the categories buttons for the first time
  	buildCategoryRadioButtons()
  	
  	updateBuffer()
  	
  	// create the renderer
  	var gl_context = getContext(canvas);
	gl_context = initContext(gl_context)
	
	// now draw
	draw()
	
	// draw the scale
	var canvasColorScale = document.getElementById('canvasColorScale'),
		canvas_ctx = canvasColorScale.getContext('2d'),
		scale_gradient = canvas_ctx.createLinearGradient(0, 0, 200, 0);
		
	scale_gradient.addColorStop(0, 'blue');
	scale_gradient.addColorStop(0.5, 'green');
	scale_gradient.addColorStop(1, 'red');
	canvas_ctx.fillStyle = scale_gradient;
	canvas_ctx.fillRect(0, 20, canvasColorScale.width, canvasColorScale.height)
	canvas_ctx.fillStyle = 'black'
	canvas_ctx.fillText('0', 10, 10)
	canvas_ctx.fillText('8', 190, 10)
  </script>
</body>
</html>
"""

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
max_expression = round(raw_expression.max())
raw_expression /= max_expression
for index, gene_name in enumerate(gene_names):
    gene_expression = raw_expression[:, index]
    gene_expression = [str(value)[:min(4, len(str(value)))] for value in gene_expression]
    gene_expression = ",".join(gene_expression)
    gene_expression_colour_coded.append("gene_expression['{gn}'] = [{ge}]".format(gn = gene_name, ge = gene_expression))
    gene_options.append("<option value='{gn}'>{gn}</option>".format(gn = gene_name))
gene_options = "".join(gene_options)  
gene_expression_colour_coded = ";".join(gene_expression_colour_coded)

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
        
template_str = template_str.replace('gene_options_here', gene_options)
template_str = template_str.replace('gene_expression_colour_coded', gene_expression_colour_coded)
template_str = template_str.replace('coordinates_data_here', coordinates)
template_str = template_str.replace('category_options_here', categories_options)
template_str = template_str.replace('categories_colours_data_here', categories_colours)
template_str = template_str.replace('categories_indices_data_here', categories_indices)
with open(save_to, 'w') as result:
    result.write(template_str)
