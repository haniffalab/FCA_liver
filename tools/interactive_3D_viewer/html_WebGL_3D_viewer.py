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
				<input type='range' name = 'particleSizeBar' min = 10 max = 300 onchange='setParticleSize(value)' value = 150 /><br />
				
				<label for = 'alphaInput'>Transparency: </label>
				<input type='range' name = 'alphaInput' min = 0 max = 1000 onchange='setAlpha(value)' value = 1000 /><br />
				
				<label for = 'canvasSizeInput'>Canvas size: </label>
				<input type='range' name = 'canvasSizeInput' min = 200 max = 2000 onchange='setCanvasSize(value)' value = 500 /><br />
				
				<label for = "zoom">Zoom: </label>
				<input type='range' name = 'zoom' min = 100 max = 1000 onchange='setZoom(value)' value = 400 /><br />
				
				<label for = 'bgInput'>Dark background: </label>
				<input type='radio' name = 'bgInput' onchange='setBackground(value)' value = 'dark' />
				<label for = 'bgInput'>White background: </label>
				<input type='radio' name = 'bgInput' onchange='setBackground(value)' value = 'white' checked />
				<br />
				
			</th>
		</tr> 
		<tr>
			<td style='vertical-align: text-top' >
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
  	uniform mat4 u_ModelMatrix;
  	uniform mat4 u_ViewMatrix;
  	uniform mat4 u_ProjMatrix;
  	uniform float u_basePointSize;
  	uniform float u_Alpha;
  	varying vec4 v_Color;
  	void main() {
  		vec4 cubePos = u_ProjMatrix * u_ViewMatrix * u_ModelMatrix * a_Position;
  		float currentWidth = 0.0;
  		//currentWidth = u_basePointSize;
  		currentWidth = 3.0 + (u_basePointSize - 3.0) * (1.0 - cubePos.z / cubePos.w) / 2.0;
    	gl_Position = cubePos;
    	gl_PointSize = currentWidth;
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
  		vec2 D = vec2(0.0, 0.0), centers = vec2(.65, .35);
  		float light = 0.0;
  		light = length(centers - gl_PointCoord);
  		light = .1 + .9 * (pow(50.0, -light));
    	gl_FragColor = v_Color * light + (1.0 - light) * vec4(0.0, 0.0, 0.0, 1.0);
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
  			new_indices.push(6*v)
  			new_indices.push(6*v+1)
  			new_indices.push(6*v+2)
  			new_indices.push(6*v+3)
  			new_indices.push(6*v+4)
  			new_indices.push(6*v+5)
  		}
  		current_data_buffer = []
  		new_indices.forEach(function(val, i){current_data_buffer.push(data_buffer[val])})
  		current_data_buffer = new Float32Array(current_data_buffer)
  		
  		gl_context.bufferData(gl_context.ARRAY_BUFFER, current_data_buffer, gl_context.STATIC_DRAW); // load data to buffer
  		n = current_data_buffer.length / 6
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
  	
  	function setZoom(value){
  		eyeVN = parseInt(value)
  		farField = eyeVN + 100;
  		rotateData(0, 0)
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
  		n = current_data_buffer.length / 6
  		var vertexColourBuffer = gl.createBuffer()
  		gl.bindBuffer(gl.ARRAY_BUFFER, vertexColourBuffer)
  		gl.bufferData(gl.ARRAY_BUFFER, current_data_buffer, gl.STATIC_DRAW)
  	
  		var FSIZE = data_buffer.BYTES_PER_ELEMENT;
  	
  		var a_Position = gl.getAttribLocation(gl.program, 'a_Position')
  		gl.vertexAttribPointer(a_Position, 3, gl.FLOAT, false, FSIZE * 6, 0)
  		gl.enableVertexAttribArray(a_Position)
  	
  		var a_Color = gl.getAttribLocation(gl.program, 'a_Color')
  		gl.vertexAttribPointer(a_Color, 3, gl.FLOAT, false, FSIZE * 6, 3 * FSIZE)
  		gl.enableVertexAttribArray(a_Color)
  	
  		u_basePointSize = gl.getUniformLocation(gl.program, 'u_basePointSize')
  		gl.uniform1f(u_basePointSize, particleSize)
  	
  		u_Alpha = gl.getUniformLocation(gl.program, "u_Alpha")
  		gl.uniform1f(u_Alpha, alphaValue)
  		
  		u_ModelMatrix  = gl.getUniformLocation(gl.program, 'u_ModelMatrix');
  		u_ViewMatrix   = gl.getUniformLocation(gl.program, 'u_ViewMatrix');
  		u_ProjMatrix   = gl.getUniformLocation(gl.program, 'u_ProjMatrix');
  		
  		modelMatrix = new Matrix4(); // The model matrix
  		viewMatrix  = new Matrix4();  // The view matrix
  		projMatrix  = new Matrix4();  // The projection matrix
  		
  		modelMatrix.setTranslate(0, 0, 0);  // 
  		viewMatrix.setLookAt(eyeX, eyeY, eyeZ, 0, 0, 0, upX, upY, upZ); // eyeX, eyeY, eyeZ, camX, camY, camZ, upX, upY, upY
  		projMatrix.setPerspective(30, canvas.width/canvas.height, 100, farField); // fov, ratio, near, far
  		// Pass the model, view, and projection matrix to the uniform variable respectively
  		gl.uniformMatrix4fv(u_ModelMatrix, false, modelMatrix.elements);
  		gl.uniformMatrix4fv(u_ViewMatrix, false, viewMatrix.elements);
  		gl.uniformMatrix4fv(u_ProjMatrix, false, projMatrix.elements);
  		
  		gl.clearColor(1, 1, 1, 1); // add ternary conditional
  	
  		gl.enable(gl.DEPTH_TEST)
  		gl.enable(gl.BLEND)
  		gl.blendFunc(gl.SRC_ALPHA, gl.ONE_MINUS_SRC_ALPHA)
  		//gl.blendFunc(gl.ONE, gl.ONE_MINUS_SRC_ALPHA)
  	
  		gl.clear(gl.COLOR_BUFFER_BIT);
  		return gl
  	}
  	
  	var canvas = document.getElementById('canvas'),
  		particleSize = 150,
  		alphaValue = 1.0,
  		bg_color = "dark",
  		eyeX    = 0.0,
  		eyeY    = 0.0,
  		eyeZ    = 400.0,
  		upX     = 0.0,
  		upY     = 1.0,
  		upZ     = 0.0,
  		eyeVN   = 400.0,
  		farField = 500.0,
  		previousX = null,
  		previousY = null,
  		currentX  = null,
  		currentY  = null;
  		
  	data_buffer = new Float32Array([
  		datahere
  	])
  	
  	index_table = []
  	indiceshere
  	
  	current_data_buffer = data_buffer
  	
  	gl_context = getContext(canvas)
	gl_context = initContext(gl_context)
  	gl_context.drawArrays(gl_context.POINTS, 0, n)
  	
  	function negCrossProduct(vecA, vecB){
  		crossproduct = [ - vecA[1] * vecB[2] + vecA[2] * vecB[1],
  					     - vecA[2] * vecB[0] + vecA[0] * vecB[2],
  					     - vecA[0] * vecB[1] + vecA[1] * vecB[0]
  		]
  		return(crossproduct)
  	}
  	
  	function vectNorm(vector){
  		return(Math.sqrt((vector[0] * vector[0]) + (vector[1] * vector[1]) + (vector[2] * vector[2])))
  	}
  	
  	function rotateData(hAngle, vAngle){
  		// change vector for very small angles is approximately the cross product of the eye vector and up vector
  		change = negCrossProduct([eyeX, eyeY, eyeZ], [upX, upY, upZ])
  		// normalize the change vector
  		normChange = vectNorm(change)
  		// scale the change vector by the horizontal angle
  		change = [hAngle * change[0]/normChange, hAngle * change[1]/normChange, hAngle * change[2]/normChange]
  		// update the eye vector by adding the change vector
  		eyeX = eyeX - change[0]
  		eyeY = eyeY - change[1]
  		eyeZ = eyeZ - change[2]
  		// renormalize the eye vector, otherwise it will increase with each change (due to approx error)
  		normEye = vectNorm([eyeX, eyeY, eyeZ])
  		eyeX = eyeVN * eyeX /  normEye
  		eyeY = eyeVN * eyeY /  normEye
  		eyeZ = eyeVN * eyeZ /  normEye
  		
  		// get the (eye, up) plane normal
  		planeInvNormal = negCrossProduct([eyeX, eyeY, eyeZ], [upX, upY, upZ])
  		// in the case of vertical angle, the up vector is already the change vector
  		normChange = vectNorm([upX, upY, upZ])
  		change = [vAngle * upX / normChange, vAngle * upY / normChange, vAngle * upZ / normChange]
  		// update the eye vector by adding the change vector
  		eyeX = eyeX + change[0]
  		eyeY = eyeY + change[1]
  		eyeZ = eyeZ + change[2]
  		// renormalize the eye vector, other it will increase with each change (due to approx error)
  		normEye = Math.sqrt((eyeX * eyeX)+(eyeY * eyeY)+(eyeZ * eyeZ))
  		eyeX = eyeVN * eyeX /  normEye
  		eyeY = eyeVN * eyeY /  normEye
  		eyeZ = eyeVN * eyeZ /  normEye
  		// but the up vector needs changing as well
  		newUp = negCrossProduct([eyeX, eyeY, eyeZ], planeInvNormal)
  		newUpNormal = vectNorm(newUp)
  		upX = -newUp[0] / newUpNormal
  		upY = -newUp[1] / newUpNormal
  		upZ = -newUp[2] / newUpNormal
  		
  		gl_context.clear(gl_context.COLOR_BUFFER_BIT);
  		viewMatrix.setLookAt(eyeX, eyeY, eyeZ, 0, 0, 0, upX, upY, upZ);
  		projMatrix.setPerspective(30, canvas.width/canvas.height, 100, farField); 
  		gl_context.uniformMatrix4fv(u_ViewMatrix, false, viewMatrix.elements);
  		gl_context.uniformMatrix4fv(u_ProjMatrix, false, projMatrix.elements);
  		gl_context.drawArrays(gl_context.POINTS, 0, n);
  	}
  	
  	function startRotating(ev){
  		previousX = ev.clientX
  		previousY = ev.clientY
  		canvas.addEventListener('mousemove', rotateEvent)
  		canvas.addEventListener('mouseup', stopRotation)
  		canvas.addEventListener('mouseout', stopRotation)
  	}
  	
  	function stopRotation(ev){
  		canvas.removeEventListener('mousemove', rotateEvent)
  		canvas.removeEventListener('mouseup', stopRotation)
  		canvas.removeEventListener('mouseout', stopRotation)
  	}
  	
  	function rotateEvent(ev){
  		currentX = ev.clientX
  		currentY = ev.clientY
  		var dX = currentX - previousX,
  			dY = currentY - previousY;
  		rotateData(2.0 * dX, 2.0 * dY)
  		previousX = currentX;
  		previousY = currentY;
  	}
  	
  	canvas.addEventListener('mousedown', startRotating)
  	
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
coordinates = data.values[:, 0:3].astype('float32')
# next few steps are compressing the data into a stadard cube centered at (0,0,0) and L = 200 
Xrange = np.percentile(coordinates[:, 0],  q = [1, 99]) * 1.2
Yrange = np.percentile(coordinates[:, 1],  q = [1, 99]) * 1.2
Zrange = np.percentile(coordinates[:, 2],  q = [1, 99]) * 1.2
center = np.tile(np.array([np.mean(Xrange), np.mean(Yrange), np.mean(Zrange)]), 
                 (coordinates.shape[0], 1))
coordinates = coordinates - center
Xrange = Xrange[1] - Xrange[0]
Yrange = Yrange[1] - Yrange[0]
Zrange = Zrange[1] - Zrange[0]
maxRange = max((Xrange, Yrange, Zrange))
ratio = 180.0 / maxRange
coordinates = coordinates * ratio

# next few steps the buffer data is created as string
colours     = data.values[:, 4]
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
     