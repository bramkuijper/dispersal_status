#!/usr/bin/env python

import os, re, sys

first = True


def analyze_parameters(lines,first=False):

    pars = {}

    for line in lines:
        mobj = re.match("(.*);(.*)",line)
        if mobj != None:
            pars[mobj.group(1)] = mobj.group(2)

    return(pars)

def analyze_data(lines):

    data = [ [ float(celli) ] for celli in lines[0].split(";")[0:-1] ]

    # loop through the lines to collect the data
    for line in lines[1:]:
        splitline = line.split(";")[0:-1]

        for i in range(0,len(splitline)):
            data[i].append(float(splitline[i]))

    # now take averages

    avgs = []
    for i in range(0,len(data)):
        avgs.append(sum(data[i])/len(data[i]))

    return(avgs)

def analyze_file(filename):

    global first;

    # open file; read data
    f = open(filename)
    fl = f.readlines()
    f.close

    flhead = fl[0]

    lc = len(fl)
    parline = -1

    linerange = range(0,lc)
    linerange.reverse()

    # search the parameter line
    for lineno in linerange:

        if re.match("^mum;",fl[lineno]) != None:
            parline = lineno
            break

    if parline == -1:
        return

    parameters = analyze_parameters(fl[parline:])

    data = fl[parline-3].strip()

    if first:
        print ";".join(parameters.keys()) + ";" + flhead.strip() + "file"
        first = False

    print ";".join(parameters.values()) + ";" + data + filename


def visit(arg, dirname, names):
    for name in names:
        if re.match("(sim|iter).*",name) != None:
#            print dirname + "/" + name
            data = analyze_file(dirname + "/" + name)



os.path.walk(sys.argv[1], visit, None)
