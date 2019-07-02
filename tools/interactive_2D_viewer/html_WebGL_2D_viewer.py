#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  7 13:45:37 2018

@author: Dorin-Mirel Popescu
"""

import sys
args = sys.argv
save_to   = args[1]
data_frame_fname = args[2]

template = """
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
			<th align='left'>
			
				<label for = 'particleSizeBar'>Particle size: </label>
				<input type='range' name = 'particleSizeBar' min = 0 max = 50 onchange='setParticleSize(value)' value = 5 /><br />
				
				<label for = 'alphaInput'>Transparency: </label>
				<input type='range' name = 'alphaInput' min = 0 max = 1000 onchange='setAlpha(value)' value = 1000 /><br />
				
				<label for = 'canvasSizeInput'>Canvas size: </label>
				<input type='range' name = 'canvasSizeInput' min = 200 max = 2000 onchange='setCanvasSize(value)' value = 500 /><br />
				
				<label for = 'bgInput'>Dark background: </label>
				<input type='radio' name = 'bgInput' onchange='setBackground(value)' value = 'dark' />
				<label for = 'bgInput'>White background: </label>
				<input type='radio' name = 'bgInput' onchange='setBackground(value)' value = 'white' checked />
				<br />
				
			</th>
		</tr> 
		<tr>
			<td>
				<canvas id='canvas' width=600 height=600></canvas>
			</td>
			<td style='vertical-align: top'>
				<label>Categories: </label><br />
				<label for='toggleRadio'>Show all:</label>
				<input type='checkbox' name = 'toggleRadio' id='toggleRadio' onchange='toggleShowTypes()' checked /><br />
				<form id = 'ControlPanel'>
					radiocommands
				</form>
			</td>
		</tr>
	</table>
  <script id='vertex-shader' type='x-shader/x-fragment'>
  	attribute vec4 a_Position;
  	attribute vec3 a_Color;
  	uniform float u_basePointSize;
  	uniform float u_Alpha;
  	varying vec4 v_Color;
  	void main() {
    	gl_Position = a_Position;
    	gl_PointSize = u_basePointSize;
    	v_Color = vec4(a_Color, u_Alpha);
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
  <script type='text/javascript'>
  
  	function selectData(){
  		controlPanel = document.getElementById('ControlPanel')
  		controlRadios = controlPanel.elements
  		values = []
  		for(i=0;i<controlRadios.length;i++){
  			if(controlRadios[i].checked){
  				values = values.concat(index_table[controlRadios[i].id])
  			}
  		}
  		new_indices = []
  		for (i=0;i<values.length;i++){
  			v = values[i]
  			new_indices.push(5*v)
  			new_indices.push(5*v+1)
  			new_indices.push(5*v+2)
  			new_indices.push(5*v+3)
  			new_indices.push(5*v+4)
  		}
  		current_data_buffer = []
  		new_indices.forEach(function(val, i){current_data_buffer.push(data_buffer[val])})
  		current_data_buffer = new Float32Array(current_data_buffer)
  		
  		gl_context.bufferData(gl_context.ARRAY_BUFFER, current_data_buffer, gl_context.STATIC_DRAW); // load data to buffer
  		n = current_data_buffer.length / 5
  		gl_context.clear(gl_context.COLOR_BUFFER_BIT)
  		gl_context.drawArrays(gl_context.POINTS, 0, n)
  	}
  	
  	function toggleShowTypes(){
  		toggleRadio = document.getElementById('toggleRadio')
  		controlPanel = document.getElementById('ControlPanel')
  		controlRadios = controlPanel.elements
  		for(i=0;i<controlRadios.length;i++){
  			controlRadios[i].checked = toggleRadio.checked
  		}
  		selectData()
  	}
  
  	function setParticleSize(value){
  		particleSize = parseInt(value)
  		gl_context.uniform1f(u_basePointSize, particleSize)
  		gl_context.clear(gl_context.COLOR_BUFFER_BIT)
  		gl_context.drawArrays(gl_context.POINTS, 0, n)
  	}
  	
  	function setAlpha(value){
  		alphaValue = parseInt(value) / 1000
  		gl_context.uniform1f(u_Alpha, alphaValue)
  		gl_context.clear(gl_context.COLOR_BUFFER_BIT);
  		gl_context.drawArrays(gl_context.POINTS, 0, n)
  	}
  	
  	function setCanvasSize(value){
  		value = parseInt(value)
  		canvas.width = value
  		canvas.height = value
  		gl_context = getContext(canvas)
		gl_context = initContext(gl_context)
		gl_context.viewport(0, 0, canvas.width, canvas.height)
		if(bg_color == "white"){
			gl_context.clearColor(1, 1, 1, 1)
		}else{
			gl_context.clearColor(0, 0, 0, 1)
		}
		gl_context.clear(gl_context.COLOR_BUFFER_BIT)
  		gl_context.drawArrays(gl_context.POINTS, 0, n)
  	}
  	
  	function setBackground(value){
  		if(value == "dark"){
  			gl_context.clearColor(0, 0, 0, 1)
  			bg_color = "dark"
  		}else{
  			gl_context.clearColor(1, 1, 1, 1)
  			bg_color = "white"
  		}
  		gl_context.clear(gl_context.COLOR_BUFFER_BIT)
  		gl_context.drawArrays(gl_context.POINTS, 0, n)
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
  		n = current_data_buffer.length / 5
  		var vertexColourBuffer = gl.createBuffer()
  		gl.bindBuffer(gl.ARRAY_BUFFER, vertexColourBuffer)
  		gl.bufferData(gl.ARRAY_BUFFER, current_data_buffer, gl.STATIC_DRAW)
  	
  		var FSIZE = data_buffer.BYTES_PER_ELEMENT;
  	
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
  	
  		gl.clearColor(1, 1, 1, 1); // add ternary conditional
  		if(bg_color == "dark"){
  			gl.clearColor(0, 0, 0, 1)
  		}
  	
  		gl.disable(gl.DEPTH_TEST)
  		gl.enable(gl.BLEND)
  		gl.blendFunc(gl.SRC_ALPHA, gl.ONE_MINUS_SRC_ALPHA)
  		//gl.blendFunc(gl.ONE, gl.ONE_MINUS_SRC_ALPHA)
  	
  		gl.clear(gl.COLOR_BUFFER_BIT);
  		return gl
  	}
  	
  	var canvas = document.getElementById('canvas'),
  		particleSize = 5,
  		alphaValue = 1.0,
  		bg_color = "white"
  		
  	data_buffer = new Float32Array([
  		datahere
  	])
  	
  	index_table = []
  	indiceshere
  	
  	current_data_buffer = data_buffer
  	
  	gl_context = getContext(canvas)
	gl_context = initContext(gl_context)
  	gl_context.drawArrays(gl_context.POINTS, 0, n)
  	
  </script>
</body>
</html>
"""

import pandas as pd
import numpy as np

data = pd.read_csv(data_frame_fname, index_col = None)

# convert Colours to r, g, b values, then to floats < 1.0
def hexdec_to_1floats(hexdec):
    return np.array([int(hexdec[1:][i:(i+2)], 16) for i in (0, 2, 4)]) / 255.0

# map Labels to colours
labels = sorted(list(data.Labels.unique()))
index_table = []
radio_commands = []
for index, label in enumerate(labels):
    indices = data.Labels == label
    indices = indices.values
    indices = np.where(indices)
    indices = ','.join([str(i) for i in indices[0]])
    indices = "[{indices}]".format(indices = indices)
    index_table.append("index_table['{label}'] = {indices}".format(label = label, indices = indices))
    colour = data.Colours[data.Labels == label].values[0]
    radio_command = "<div style='background-color:{colour}'><input style='float:left' type='checkbox' id='{label}' checked onchange='selectData()' /><label style='float:left' for='{label}'>{label}: </label><br /></div>".format(colour = colour, label = label)
    radio_commands.append(radio_command)
index_table = ';\n  	'.join(index_table)
radio_commands = '\n					'.join(radio_commands)

# make data string
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

# next few steps the buffer data is created as string
colours     = data.Colours.values
buffer_data = []
for index in range(coordinates.shape[0]):
    coordinate = [str(i) for i in coordinates[index, :]]
    colour = [str(i) for i in hexdec_to_1floats(colours[index]).astype('float32')]
    vertex_data = coordinate + colour
    buffer_data.append(",".join(vertex_data))
buffer_data = ",".join(buffer_data)     

template_str = template.replace('datahere', buffer_data)
template_str = template_str.replace('indiceshere', index_table)
template_str = template_str.replace('radiocommands', radio_commands)
with open(save_to, 'w') as result:
    result.write(template_str)
     