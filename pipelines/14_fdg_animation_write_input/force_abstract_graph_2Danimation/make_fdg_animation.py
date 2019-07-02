
from os import mkdir
from os.path import exists
from shutil import rmtree
from fa2 import ForceAtlas2
import pandas as pd
from scipy.io import mmread
import numpy as np
import subprocess

# smaller steps by:
# - decrease barnesHutOptimize
# - decrease gravity

# number of frames
frames = 2000
# load pca, SNN and label colours data
# the first 2 PC form PCA are used as initial conditions
# SNN is used for building the force directed graph
pca_data = pd.read_csv("./input/pca.csv", index_col = 0)
labels_col   = pd.read_csv("./input/label_colours.csv", squeeze = True, index_col = 0)
snn      = mmread("./input/SNN.smm")

# set initialposition as the first 2 PCs
positions = pca_data.values[:, 0:2]

# initialize force directed graph class instance
forceatlas2 = ForceAtlas2(outboundAttractionDistribution=False, linLogMode=False,
                          adjustSizes=False, edgeWeightInfluence=1.0,
                          jitterTolerance=1.0, barnesHutTheta = .8,
                          barnesHutOptimize=True, multiThreaded=False,
                          scalingRatio=2.0, strongGravityMode=True, gravity=1, verbose=True)

# run force directed graph; for each iterations generates the coordinates use din each frame
discard = forceatlas2.forceatlas2(G = snn, pos = positions, iterations = frames)

if exists("./input/buffers"):
    rmtree("./input/buffers")
if exists("./input/frames"):
    rmtree("./input/frames")
mkdir("./input/buffers")
mkdir("./input/frames")
for index in range(len(forceatlas2.dataContainer)):
    positions = forceatlas2.dataContainer[index]
    fname = "./input/buffers/{index}.csv".format(index = index)
    np.savetxt(fname, positions, delimiter = ",")
    print("Saving buffer: {index}".format(index = index))
    
# run R
subprocess.call(["Rscript", "make_plots.R"], shell = True)


# assemble the frames into a video
    
import cv2
import os

def sortImages(imgPath):
    return int(os.path.splitext(imgPath)[0])

# Arguments
dir_path = './input/frames'
ext = "png"
output = "fdg.mp4"

images = []
for f in os.listdir(dir_path):
    if f.endswith(ext):
        images.append(f)

images = sorted(images, key = sortImages)

legend = cv2.imread("./input/legend.png")
lH, lW, chs = legend.shape
legend = legend[0:(lH-10), 10:lW]
legend = cv2.resize(legend, (0, 0), fx = .8, fy = .8)
lH, lW, chs = legend.shape

# Determine the width and height from the first image
image_path = os.path.join(dir_path, images[0])
frame = cv2.imread(image_path)
cv2.imshow('video',frame)
height, width, channels = frame.shape

# Define the codec and create VideoWriter object
fourcc = cv2.VideoWriter_fourcc(*'mp4v') # Be sure to use lower case
out = cv2.VideoWriter(output, fourcc, 30.0, (width+792, height))
import numpy as np
for image in images:

    image_path = os.path.join(dir_path, image)
    frame = cv2.imread(image_path)
    frame = cv2.resize(frame, (width, height))
    lh1 = width + lW
    template = np.zeros((height, lW, 3), dtype = frame.dtype)
    frame = np.hstack((frame, template))
    frame[0:lH, width:lh1, :] = legend

    #cv2.putText(frame, "by Dorin-Mirel Popescu", (width - 400, height - 30), cv2.FONT_HERSHEY_SIMPLEX, 0.6, (255, 255, 255), thickness = 2) 

    out.write(frame) # Write out frame to video
    print(image)

# Release everything if job is finished
out.release()
cv2.destroyAllWindows()

print("The output video is {}".format(output))
