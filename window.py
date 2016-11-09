#! /usr/bin/env python



import sys
#1 meter per second is 26.8224 miutes per mile
#1609.34 meters per mile

#look at line 971 for time without distance

import numpy as np
import matplotlib.pyplot as plt

METER_PER_MILE = 1609.34


distsel = 1000
intrack = 0
intrackpoint = 0
timestamp = []
diststamp = []
speedstamp = []
allInSec = 0
allInSecOld = 0
dist = 0
distOld = 0
timestr = ""

def speedCalc (tsamp, dsamp, tsampPrev, dsampPrev):
    #use difference in distance and differece in time
    #to calculate the pace between two samples.
    deld = dsamp - dsampPrev
    delt = tsamp - tsampPrev
    mps = float(deld / delt) #meters per second
    #accomodate div/0
    if mps == 0:
        return 0
    else:
        return (26.8224 / mps)
#mps / 1609.34 = mile per second
#*60 mile per mile

def paceFilter (paceAry):
    #moving window filter to remove bogus samples
    i = 0
    paceFiltered = []
    while i < len(paceAry):
        if i < 2:
            paceFiltered.append((paceAry[i] + paceAry[i + 1] + paceAry[i + 2]) / 3)
        elif i > (len(paceAry) - 3):
            paceFiltered.append((paceAry[i] + paceAry[i - 1] + paceAry[i - 2]) / 3)
        else:
            paceFiltered.append((paceAry[i - 1] + paceAry[i] + paceAry[i + 1]) / 3)
#paceFiltered.append((paceAry[i - 2] + paceAry[i - 1] + paceAry[i] + paceAry[i + 1] + paceAry[i + 2]) / 5)
        i = i + 1
    return paceFiltered

def interp_dist(diststampEnd, diststampEndLast, timestampEnd, timestampEndLast, distsel, diststampStart, timestampStart):
    deltat = timestampEnd - timestampEndLast
    deltad = diststampEnd - diststampStart
    deltadLast = diststampEndLast - diststampStart
    dif = deltad - deltadLast
    difLast = distsel - deltadLast
    distcoef = difLast / dif
    timeadj = distcoef * deltat
    return (timeadj + timestampEndLast) - timestampStart


infile = open("samp.tcx", "r")

for line in infile:
	linespl = line.strip()
    #set flags for tcx samples
	if (linespl == "<Track>"):
		intrack = 1
	if (linespl == "</Track>"):
		intrack = 0
	if (linespl == "<Trackpoint>"):
		intrackpoint = 1
	if (linespl == "</Trackpoint>"):
		intrackpoint = 0

#<Time>2016-10-31T21:02:11Z</Time>
	if "<Time>" in linespl and intrackpoint == 1:
		justtime = linespl.split("T")
		justtime = justtime[2].split("Z")
		justtime = justtime[0]
		justtimeAry = justtime.split(":")
		hour =   int(justtimeAry[0])
		minute = int(justtimeAry[1])
		second = int(justtimeAry[2])
		allInSecOld = allInSec
		allInSec = second + (60 * minute) + (60 * 60 * hour)
    
	if "<DistanceMeters>" in linespl and intrackpoint == 1:
		enabletime = 0
		distOld = dist
		dist = linespl.split("<DistanceMeters>")
		dist = dist[1].split("</DistanceMeters>")
		dist = dist[0]
		timestamp.append(float(allInSec))
		diststamp.append(float(dist))
		if len(diststamp) == 1:
			speedstamp.append(speedCalc(float(allInSec), float(dist), float(allInSec) - 0.001, float(dist)))
		else:
			speedstamp.append(speedCalc(float(allInSec), float(dist), float(allInSecOld) - 0.001, float(distOld)))
infile.close()


#Find minimal time for selected distance
mintime = None
mini = None  #minimal start element
minj = None  #minimal finish element
i = 0
while (i < len(timestamp)):
	deltad = 0
	deltat = 0
	j = i
	while j < len(timestamp):
		deltad = float(diststamp[j]) - float(diststamp[i])
		deltat = float(timestamp[j]) - float(timestamp[i])
		if deltad > distsel:
            #This is where intermediate distance function would go to get deltat
            		deltat = interp_dist(diststamp[j], diststamp[j-1], timestamp[j], timestamp[j-1], distsel, diststamp[i], timestamp[i])
			if deltat < mintime or mintime == None:
				mintime = deltat
				mini = i
				minj = j
			j = len(timestamp)
		j = j + 1
	i = i + 1

timestr = str(int(mintime / 60)) + ":" + str((mintime % 60))
print timestr


#Old Velocity Calc, Now in file read loop
'''
i = 0
while i < len(timestamp):
    if i == 0:
        speedstamp.append(speedCalc(timestamp[i], diststamp[i], timestamp[i] - 0.0001, diststamp[i]))
    else:
        speedstamp.append(speedCalc(timestamp[i], diststamp[i], timestamp[i - 1] - 0.0001, diststamp[i - 1]))
    i = i + 1
'''

#matplotlib stuff
X = diststamp

S = paceFilter(speedstamp)

a = diststamp[mini] #hilight start
b = diststamp[minj] #hilight end
                
plt.axvspan(a, b, color='r', alpha=0.1, lw=2)
plt.plot(X,S)
plt.gca().invert_yaxis()
plt.gca().grid(True)
plt.xlabel("Distance - Time for Distance: %s" % timestr)
plt.ylabel('Pace')
plt.yticks(np.arange(min(S), max(S)+1, 0.5))

plt.show()
