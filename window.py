#! /usr/bin/env python


import sys
#1 meter per second is 26.8224 miutes per mile
#1609.34 meters per mile

#look at line 971 for time without distance

import numpy as np
import matplotlib.pyplot as plt


FILENAME = "samp.tcx"



METER_PER_MILE = 1609.34
METER_PER_KILO = 1000



distsel = 1000
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

def meanMax(diststamp, timestamp):
    maxdist = diststamp[len(diststamp) - 1] #distance total
    MMT = []
    MMD = []
    MMP = []
    i = 0
    while i < maxdist:
        MMD.append(i)
        timer = minTimeForDist(diststamp, timestamp, i)[0]
        if timer == None:
            MMT.append(0)
        else:
            MMT.append((timer)/60)
        if timer == 0:
            MMP.append(0)
        else:
            MMP.append(26.8224 / (i / timer))
        i = i + 50
    return [MMT, MMD, MMP]

def interp_dist(diststampEnd, diststampEndLast, timestampEnd, timestampEndLast, distsel, diststampStart, timestampStart):
    deltat = timestampEnd - timestampEndLast
    deltad = diststampEnd - diststampStart
    deltadLast = diststampEndLast - diststampStart
    dif = deltad - deltadLast
    difLast = distsel - deltadLast
    distcoef = difLast / dif
    timeadj = distcoef * deltat
    return (timeadj + timestampEndLast) - timestampStart

def secToMinAry(timeSec):
    minuteAry = []
    for sample in timeSec:
        minuteAry.append((sample / 60))
    return minuteAry

def meterToK(dist):
    kAry = []
    for sample in dist:
        kAry.append((sample / METER_PER_KILO))
    return kAry

def minTimeForDist(diststamp, timestamp, distsel):
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
            #deltat = float(timestamp[j]) - float(timestamp[i])
            if deltad > distsel:
		deltat = interp_dist(diststamp[j], diststamp[j-1], timestamp[j], timestamp[j-1], distsel, diststamp[i], timestamp[i])
		if deltat < mintime or mintime == None:
                    mintime = deltat
                    mini = i
                    minj = j
                j = len(timestamp)
            j = j + 1
        i = i + 1
    result = [mintime, mini, minj]
    return result

infile = open(FILENAME, "r")

for line in infile:
	linespl = line.strip()
    #set flags for tcx samples
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

result = minTimeForDist(diststamp, timestamp, distsel)
#result = [mintime, mini, minj]
#mean maximal pace
#return [MMT, MMD, MMP]
MMR = meanMax(diststamp, timestamp)
mini = result[1]
minj = result[2]
mintime = result[0]
timestr = str(int(mintime / 60)) + ":" + str((mintime % 60))
print timestr

#matplotlib stuff

D = meterToK(diststamp)
S = paceFilter(speedstamp)

a = D[mini]
b = D[minj]
c = timestamp[mini] #hilight start time
d = timestamp[minj] #hilight end time

plt.subplot(3, 1, 1)
plt.title('Pace')
plt.axvspan(a, b, color='r', alpha=0.1, lw=2)
plt.plot(D,S, 'b', lw = 1.5, alpha = 0.8)
plt.gca().invert_yaxis()
plt.gca().grid(True)
plt.xlabel("Distance - Time for Distance: %s" % timestr)
plt.ylabel('Pace')
plt.yticks(np.arange(min(S), max(S)+1, 1))

plt.subplot(3, 1, 2)
plt.title('Maximal Pace')
plt.axvspan(0, distsel, color='b', alpha=0.25, lw=2)
plt.plot(MMR[1],MMR[2], 'r', lw = 1.5, alpha = 0.8)
#plt.xscale('log')
plt.gca().grid(True)
plt.xlabel("Distance")
plt.ylabel('Pace')
plt.yticks(np.arange(min(MMR[2]), max(MMR[2])+1, 1))


plt.subplot(3, 1, 3)
plt.title('Minimal Time')
plt.axvspan(0, distsel, color='b', alpha=0.25, lw=2)
plt.plot(MMR[1],MMR[0], 'g', lw = 1.5, alpha = 0.8)
#plt.xscale('log')
plt.gca().grid(True)
plt.xlabel("Distance")
plt.ylabel('Time')
plt.yticks(np.arange(min(MMR[0]), max(MMR[0])+1, 5))
plt.tight_layout()
plt.show()
