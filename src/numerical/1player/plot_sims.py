#!/usr/bin/env python3

import sys, re, os.path, itertools
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import pandas as pd

from pylab import *

params={'axes.linewidth' : .5}
rcParams.update(params)

# get the filename from the command line
filename = sys.argv[1]

# first figure out last line of dataset
f = open(filename)
fl = f.readlines()
f.close()

nrow = len(fl)-8

def find_parline(filename, searchvar):

    f = open(filename)
    fl = f.readlines()
    f.close()

    rx = "^" + searchvar + ".*"

    for row, line in enumerate(fl):

        if re.match("^" + searchvar + ".*", line) is not None:
            return(row)

    return(-1)

parline = find_parline(filename, "mum")

data = None

if parline != -1:

    # read the parameters
    parameters = pd.read_csv(filename, 
                                sep=";",
                                skiprows=parline-1,
                                header=None,
                                names=["var","val"])

    data = pd.read_csv(filename, sep=";", nrows=parline-3)
else:
    data = pd.read_csv(filename, sep=";")



# initialize and specify size 
fig = plt.figure(figsize=(10,10))

num_rows = 5

# names of columns in dataframe
colnames = list(data.columns.values)

# add first subplot depicting switch rate modulation
plt.subplot(num_rows,1,1)

plt.plot(
        data["time"],data["p1imm"],'b',
        data["time"],data["p1phil"],'r',
        data["time"],data["p2imm"],'m',
        data["time"],data["p2phil"],'g',
        linewidth=1)

plt.legend((
                r'$p_{1,\mathrm{imm}}$',
                r'$p_{1,\mathrm{phil}}$',
                r'$p_{2,\mathrm{imm}}$',
                r'$p_{2,\mathrm{phil}}$'
            ))
plt.ylabel(r'Prop $z_{1}$ offspring')
plt.ylim(-0.05,1.05)

# add 2nd subplot depicting patch frequencies for immigrants
plt.subplot(num_rows,1,2)

plt.plot(data["time"],data["f1aimm"],'c',
        data["time"],data["f1mimm"],'m',
        data["time"],data["f2aimm"],'y',
        data["time"],data["f2mimm"],'k',
        linewidth=1)
plt.legend((
                r'$f_{1,a,\mathrm{imm}}$',
                r'$f_{1,m,\mathrm{imm}}$',
                r'$f_{2,a, \mathrm{imm}}$',
                r'$f_{2,m, \mathrm{imm}}$'))
plt.ylabel(r'Patch freq immigrant')
plt.ylim(-0.05,1.05)

# add 3rd subplot depicting patch frequencies for philopatrics
plt.subplot(num_rows,1,3)

plt.plot(data["time"],data["f1aphil"],'c',
        data["time"],data["f1mphil"],'m',
        data["time"],data["f2aphil"],'y',
        data["time"],data["f2mphil"],'k',
        linewidth=1)
plt.legend((
                r'$f_{1,a,\mathrm{phil}}$',
                r'$f_{1,m,\mathrm{phil}}$',
                r'$f_{2,a, \mathrm{phil}}$',
                r'$f_{2,m, \mathrm{phil}}$'))
plt.ylabel(r'Patch freq philigrant')
plt.ylim(-0.05,1.05)

# add 4th subplot depicting reproductive values immigrants
plt.subplot(num_rows,1,4)

plt.plot(data["time"],data["v1aimm"],'c',
        data["time"],data["v1mimm"],'m',
        data["time"],data["v2aimm"],'y',
        data["time"],data["v2mimm"],'k',
        linewidth=1)
plt.ylim(-0.05,2.05)
plt.legend((
                r'$v_{1,a,\mathrm{imm}}$',
                r'$v_{1,m, \mathrm{imm}}$',
                r'$v_{2,a, \mathrm{imm}}$',
                r'$v_{2,m, \mathrm{imm}}$'
                ))
plt.ylabel(r'Reproductive value')

# add 5th subplot depicting reproductive values immigrants
plt.subplot(num_rows,1,5)

plt.plot(
        data["time"],data["v1aphil"],'c',
        data["time"],data["v1mphil"],'m',
        data["time"],data["v2aphil"],'y',
        data["time"],data["v2mphil"],'k',
        linewidth=1)
plt.ylim(-0.05,2.05)
plt.legend((
                r'$v_{1,a,\mathrm{phil}}$',
                r'$v_{1,m, \mathrm{phil}}$',
                r'$v_{2,a, \mathrm{phil}}$',
                r'$v_{2,m, \mathrm{phil}}$'
                ))
plt.ylabel(r'Reproductive value')

graphname = os.path.dirname(filename)
if graphname != '':
    graphname += "/"
graphname += "graph_" + os.path.basename(filename) + ".pdf"

plt.savefig(graphname,format="pdf")
