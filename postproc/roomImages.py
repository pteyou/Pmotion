#! /usr/bin/python
import numpy as np
from matplotlib import pyplot as plt
import os

basePath = '../build/'
outFname = 'reunion.mp4'


X = np.loadtxt(basePath + "room12X.dat")
Y = np.loadtxt(basePath + "room12Y.dat")
wall = np.loadtxt(basePath + "room12Wall.dat")

nbTimeStep, nbPoints = np.shape(X)
nbPoints -= 1

for i in range(nbTimeStep):
    plt.figure(i+1)
    plt.scatter(X[i, 1:], Y[i, 1:])
    plt.scatter(wall[:, 0], wall[:, 1], c='black', marker='.', s=0.05)
    plt.xlim(-1, 9)
    plt.ylim(-1, 5)
    #plt.title("4 seconds")
    plt.savefig('room12_{}.png'.format(i))
    plt.close(i+1)


# 10 images par secondes avec un bitrate de 1800 kb/s
os.system("ffmpeg -r 10 -f image2 -s 1920*1080 -i room12_%d.png -vcodec libx264 -crf 25 -pix_fmt yuv420p -y " + outFname)
os.system("rm *.png")


print 'done'
