#! /usr/bin/python
import numpy as np
from matplotlib import pyplot as plt
import os
import math

basePath = '../build/'
outFname = 'cinema.mp4'

X = np.loadtxt(basePath + "movie350X.dat")
Y = np.loadtxt(basePath + "movie350Y.dat")
wall = np.loadtxt(basePath + "movie350Wall.dat")

area = math.pi * 0.25**2 / 4
siz = 4*area / math.pi

nbTimeStep, nbPoints = np.shape(X)
nbPoints -= 1

for i in range(nbTimeStep):
    plt.figure(i+1)
    plt.scatter(X[i, 1:], Y[i, 1:], s=siz)
    plt.scatter(wall[:, 0], wall[:, 1], c='black', marker='.', s=0.05)
    plt.xlim(-1, 30)
    plt.ylim(-1, 17)
    #plt.title("19.9 secondes")
    plt.savefig('movie350_{}.png'.format(i))
    plt.close(i+1)


# 10 images par secondes avec un bitrate de 1800 kb/s
os.system("ffmpeg -r 10 -f image2 -s 1920*1080 -i movie350_%d.png -vcodec libx264 -crf 25 -pix_fmt yuv420p -y " + outFname)
os.system("rm *.png")

print 'done'
