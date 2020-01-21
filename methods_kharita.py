import time, datetime
import numpy as np
import matplotlib
import getopt
import sys
import matplotlib.pyplot as plt
from geopy.distance import vincenty
from sklearn.neighbors import NearestNeighbors
from geojson import MultiLineString

LL = (41, -87);
latconst = vincenty(LL, (LL[0] + 1, LL[1])).meters;
lonconst = vincenty(LL, (LL[0], LL[1] + 1)).meters


def geodist(point1, point2):
#    print(point1, point2, vincenty((point1[1],point1[0]),(point2[1],point2[0])).meters,np.sqrt((lonconst*(point1[0]-point2[0]))**2+(latconst*(point1[1]-point2[1]))**2))
    return(np.sqrt((lonconst*(point1[0]-point2[0]))**2+(latconst*(point1[1]-point2[1]))**2)) #180dg difference equivalent to 80m difference

def taxidist(point1, point2,theta):
    return(lonconst*np.abs(point1[0]-point2[0])+latconst*np.abs(point1[1]-point2[1])+ theta/180*angledist(point2[2],point1[2])) #180dg difference equivalent to 80m difference

def angledist(a1, a2):
    return(min(abs(a1-a2),abs((a1-a2) % 360),abs((a2-a1) % 360),abs(a2-a1)))


def getdata(nsamples, datafile, datestart, datestr):
#3233678911,1080020,83,2015-10-03 06:52:48,57,51.4950963,25.262793500000001,PICKUP,private,100
    datapointwts = [];
    lats = []; lons = []; j = 0;
    print datafile
    with open(datafile,'rb') as f:
        for line in f:
            j = j+1;
            if j>nsamples:
                break;
            line = line[:-1].decode('ascii', 'ignore')
            zz = line.split("\t")
            if zz[6][:10]<datestr and zz[6][:10]>=datestart: 
                ts = time.mktime(datetime.datetime.strptime(zz[6][:-3], "%Y-%m-%d %H:%M:%S").timetuple())
                LL = (float(zz[0][:8]),float(zz[1][:8])); angle = float(zz[-1])-180; speed = float(zz[5])
                if j>1:
                    if oldts<ts and oldts>ts-20:
                        speed = int(geodist(LL,oldLL)/(ts-oldts)*3.6)
                lats.append(LL[0])
                lons.append(LL[1])
                pointwts = (LL[0],LL[1],angle,speed,j,ts);
#                print(pointwts)
                oldLL = LL; oldts = ts;
                datapointwts.append(pointwts)
    return(datapointwts)

def greaterthanangle(alpha,beta):
    if (beta-alpha)%360<180:
        return True
    else:
        return False

def anglebetweentwopoints(LL1, LL2):
    xd = (LL1[0]-LL2[0]); yd =LL1[1]-LL2[1];
#    xd = latconst/lonconst*(LL1[0]-LL2[0]); yd =LL1[1]-LL2[1];
    return(np.arctan2(xd,yd)*180/np.pi)

def is_power2(num):
	return num != 0 and ((num & (num - 1)) == 0)

def getseeds(datapoint,radius,theta):
    chosen = []; seeds = [];
#    random.shuffle(datapoint)
    periodsampl = 500000
    for p in datapoint:
        chosen.append(p);
    for j,p in enumerate(chosen):
        ok = -1;
        if j<periodsampl:
            for q in seeds:
                if taxidist(p,q,theta)<radius:
                    ok = 1
                    break;
            if ok <1:
                seeds.append(p)
        else:
            if j%periodsampl == 0:# and (is_power2(int(j/1000))):
#                print(j,time.time()-start)
                S = [(lonconst * xx[0], latconst * xx[1], theta / 180 * (xx[2]+45)) for xx in seeds];
                nbrs = NearestNeighbors(n_neighbors=1, algorithm='ball_tree').fit(S)
                X = [(lonconst * xx[0], latconst * xx[1], theta / 180 * (xx[2]+45)) for xx in chosen[j:min(len(chosen),j+periodsampl)]];
                distances, indices = nbrs.kneighbors(X)
            if distances[j%periodsampl][0] >radius:
                seeds.append(p)
    print('seeds: ', len(seeds))
    return (seeds)

def avgpoint(cc):
    hh = np.arctan2(sum([np.sin(xx[2] / 360 * 2 * np.pi) for xx in cc]), sum([np.cos(xx[2] / 360 * 2 * np.pi) for xx in cc])) * 180 / np.pi
    return((np.mean([xx[0] for xx in cc]), np.mean([xx[1] for xx in cc]), hh))

def newmeans(datapointwts,seeds,theta):
    newseeds = []; cost = 0; avgspeed = []; pointsperseed = [];
    cluster, p2cluster = point2cluster(datapointwts, seeds,theta);
    for cd in cluster:
        if len(cluster[cd])>0:
            hh = np.arctan2(sum([np.sin(xx[2]/360*2*np.pi) for xx in cluster[cd]]),sum([np.cos(xx[2]/360*2*np.pi) for xx in cluster[cd]]))*180/np.pi
            newseeds.append((np.mean([xx[0] for xx in cluster[cd]]),np.mean([xx[1] for xx in cluster[cd]]),hh))
            hh = [xx[3] for xx in cluster[cd] if xx[3]>0];
            if len(hh)<1:
                hh = [0]
            avgspeed.append(np.mean(hh))
            cost = cost+sum([taxidist(xx,newseeds[-1],theta) for xx in cluster[cd]])
        else:
            newseeds.append(seeds[cd])
            avgspeed.append(0)
        pointsperseed.append(len(cluster[cd]))
    return(newseeds,cost,avgspeed,pointsperseed)

def densify(datapointwts):
    newpoints = [];
    for ii, xx in enumerate(datapointwts):
        if ii>1:
            if datapointwts[ii-1][-1]<datapointwts[ii][-1] and datapointwts[ii-1][-1]>datapointwts[ii][-1]-11 and taxidist(datapointwts[ii-1],datapointwts[ii],theta)<1000:
                delta = int(taxidist(datapointwts[ii][:3],datapointwts[ii-1][:3],theta)/20)+1;
                x1 = datapointwts[ii-1]; x2 = datapointwts[ii];
                if np.abs(datapointwts[ii-1][2]-datapointwts[ii][2])<500:
                    for jj in range(1,delta-1):
                        newpoints.append(tuple([jj/delta*x1[sq]+(delta-jj)/delta*x2[sq] for sq in range(len(x1))]))
    print('original datapoints: ', len(datapointwts), 'densified datapoints:', len(newpoints))
    result = datapointwts+newpoints;
    result.sort(key=lambda x: x[-2],reverse=False)
    return(result)

def getpossibleedges(datapointwts,seeds):
#    datapointwts = densify(datapointwts);
    X = [(xx[0], xx[1]) for xx in datapointwts];    S = [(xx[0], xx[1]) for xx in seeds];cluster = {};p2cluster = []; gedges = {}; gedges1 = {}; nedges = {};
    nbrs = NearestNeighbors(n_neighbors=5, algorithm='ball_tree').fit(S)
    distances, indices = nbrs.kneighbors(X)
    for cd in range(len(seeds)):
        cluster[cd] = []
    for ii, ll in enumerate(indices):
        dd = [taxidist(seeds[xx], datapointwts[ii][:-1],theta) for xx in ll]
        cd = ll[dd.index(min(dd))];
        cluster[cd].append(datapointwts[ii])
        p2cluster.append(cd)
    for ii, xx in enumerate(datapointwts):
        if ii>1:
            if datapointwts[ii-1][-1]<datapointwts[ii][-1] and datapointwts[ii-1][-1]>datapointwts[ii][-1]-11:
                cd1 = p2cluster[ii-1]; cd2 = p2cluster[ii];
            if not cd1== cd2:
                gedges1[(cd1,cd2)] =  gedges1.get((cd1,cd2),0)+1;
    return(gedges1)

def coocurematrix(datapointwts,seeds,theta):
    startcoocurence = time.time();    gedges1 = {}; std = {};
    cluster, p2cluster = point2cluster(datapointwts, seeds,theta);
    for ii, xx in enumerate(datapointwts):
        if ii>1:
            if datapointwts[ii-1][-1]<=datapointwts[ii][-1] and datapointwts[ii-1][-1]>=datapointwts[ii][-1]-121 and taxidist(datapointwts[ii-1],datapointwts[ii],theta)<1000:
                cd1 = p2cluster[ii-1]; cd2 = p2cluster[ii];
                if (not cd1== cd2):
                    gedges1[(cd1, cd2)] =  gedges1.get((cd1,cd2),0)+1;
#                LL = datapointwts[ii]
#                if (LL[0]>-87.657) and (np.abs(LL[0])>81.6568) and (LL[1]>41.8755) and (LL[1]<41.8765) and (np.abs(LL[2]+20)<40):
#                    print(cd1,cd2,LL[:3],seeds[cd1][:2],seeds[cd2][:2])
    gedges2 = {gg: gedges1[gg] for gg in gedges1};
    for gg in gedges2:
        if gg in gedges1 and (gg[1],gg[0]) in gedges1:
            if gedges1[(gg[1],gg[0])]>gedges1[gg]:
                del gedges1[gg]
            elif gedges1[(gg[1],gg[0])]==gedges1[gg]:
                gg0 =(gg[1],gg[0]);   cd1 = gg[0];  cd2 = gg[1];
#                print(gg,cd1,cd2)
                AA = anglebetweentwopoints(seeds[cd1], seeds[cd2]);   AArev = anglebetweentwopoints(seeds[cd2], seeds[cd1]);
#                print(int(angledist(AA, seeds[cd1][2]) + angledist(AA, seeds[cd2][2])),int(angledist(AArev, seeds[cd1][2]) + angledist(AArev, seeds[cd2][2])))
                if (angledist(AA, seeds[cd1][2]) + angledist(AA, seeds[cd2][2]))<(angledist(AArev, seeds[cd1][2]) + angledist(AArev, seeds[cd2][2])):
                    del gedges1[gg0]
                else:
                    del gedges1[gg]
    gedges2 = {gg: gedges1[gg] for gg in gedges1}; neighbors = {}; filneigh = {};
    for gg in gedges2:
        neighbors[gg[0]] = [];  cd1 = gg[0];        cd2 = gg[1];
#        AA = anglebetweentwopoints(seeds[cd1], seeds[cd2]);    AArev = anglebetweentwopoints(seeds[cd2], seeds[cd1]);
#        print(int(angledist(AA, seeds[cd1][2]) + angledist(AA, seeds[cd2][2])), int(angledist(AArev, seeds[cd1][2]) + angledist(AArev, seeds[cd2][2])))
    for gg in gedges2:
        neighbors[gg[0]].append(gedges1[gg])
    for ss in neighbors:
        neighbors[ss] = sorted(neighbors[ss])
#        print(ss,len(neighbors[ss]),sum(neighbors[ss]),neighbors[ss])
    for gg in gedges2:
        hh = min(sum(neighbors.get(gg[1],[0])),sum(neighbors[gg[0]]));
        if gedges1[gg]<np.log(max(1,hh))-1:#filneigh[gg[0]]:
            del gedges1[gg]
    print(len(datapointwts), sum(gedges1.values()), 'coocurence computation time:', time.time() - startcoocurence)
    return (gedges1)

def prunegraph(gedges,seeds):
    neighbors = {}
    for ss in range(len(seeds)):
        neighbors[ss] = [];
        if (ss,ss) in gedges:
            del gedges[(ss,ss)]
    gedges1  = dict(gedges);
    for gg in gedges1:
        neighbors[gg[0]].append(gg[1])
    depth = 5;
    gedges2 = {gg:geodist(seeds[gg[0]],seeds[gg[1]]) for gg in gedges};
    gedges = {gg:geodist(seeds[gg[0]],seeds[gg[1]]) for gg in gedges};
    hopedges = []
    for dd in range(depth):
        gedges1 = dict(gedges2);
        for gg in gedges1:
            for ss in neighbors[gg[1]]:
                if not ss == gg[0]:
                    gedges2[(gg[0], ss)] = min(gedges2[(gg[0],gg[1])] + gedges[(gg[1],ss)],gedges2.get((gg[0], ss),100000))
                    hopedges.append((gg[0],ss))
#        print(len(gedges2))
    for gg in hopedges:
      if gg in gedges and gedges[gg]>0.8*gedges2[gg]:
          del gedges[gg];
    return (gedges)

def point2cluster(datapointwts,seeds,theta):
    cluster = {};p2cluster = []; gedges = {}; gedges1 = {}; nedges = {}; std = {}; seeds1 = []; seedweight = [];
    X = [(lonconst * xx[0], latconst * xx[1], theta / 180 * xx[2]) for xx in datapointwts];    S = [(lonconst * xx[0], latconst * xx[1], theta / 180 * xx[2]) for xx in seeds];
    Xrot = [(lonconst * xx[0], latconst * xx[1], theta / 180 * (xx[2]%360)) for xx in datapointwts];    Srot = [(lonconst * xx[0], latconst * xx[1], theta / 180 * (xx[2]%360)) for xx in seeds];
    for cd in range(len(seeds)):
        cluster[cd] = []
    nbrs = NearestNeighbors(n_neighbors=1, algorithm='ball_tree').fit(S)
    distances, indices = nbrs.kneighbors(X)
    nbrsrot = NearestNeighbors(n_neighbors=1, algorithm='ball_tree').fit(Srot)
    distancesrot, indicesrot = nbrsrot.kneighbors(Xrot)
    for ii, ll in enumerate(indices):
        #        print(distances[ii],distancesrot[ii],ll,indices[ii],indicesrot[ii])
        cd = indicesrot[ii][0]
        if distances[ii][0] < distancesrot[ii][0]:
            cd = indices[ii][0];
            #        print(cd,distances[ii],distancesrot[ii],ll,indices[ii],indicesrot[ii])
        cluster[cd].append(datapointwts[ii])
        p2cluster.append(cd)
    return(cluster,p2cluster)

def splitclusters(datapointwts,seeds,theta):
    std = {}; seeds1 = []; seedweight = [];
    cluster, p2cluster = point2cluster(datapointwts, seeds,theta);
    for cl in cluster:
        mang = seeds[cl][-1];
        if len(cluster[cl]) > 10:
            std[cl] = np.percentile([angledist(xx[2], mang) for xx in cluster[cl]], 90)
            clockwise = [xx for xx in cluster[cl] if greaterthanangle(xx[2], mang)];
            if std[cl]>20 and len(clockwise)>0 and len(clockwise)<len(cluster[cl]):
                seeds1.append(avgpoint(clockwise))
                seeds1.append(avgpoint([xx for xx in cluster[cl] if not greaterthanangle(xx[2], mang)]))
                seedweight.append(len(clockwise))
                seedweight.append(len(cluster[cl]) -len(clockwise))
            else:
                seeds1.append(seeds[cl]); seedweight.append(len(cluster[cl]))
        else:
            seeds1.append(seeds[cl]); seedweight.append(len(cluster[cl]))
    return seeds1, seedweight

def splitclustersparallel(datapointwts,seeds):
    X = [(xx[0], xx[1]) for xx in datapointwts];    S = [(xx[0], xx[1]) for xx in seeds];cluster = {};p2cluster = []; gedges = {}; gedges1 = {}; nedges = {}; std = {}; seeds1 = []; seedweight = []; roadwidth = [];
    nbrs = NearestNeighbors(n_neighbors=20, algorithm='ball_tree').fit(S)
    distances, indices = nbrs.kneighbors(X)
    for cd in range(len(seeds)):
        cluster[cd] = []; roadwidth.append(0);
    for ii, ll in enumerate(indices):
        dd = [taxidist(seeds[xx], datapointwts[ii][:-1],theta) for xx in ll]
        cd = ll[dd.index(min(dd))];
        cluster[cd].append(datapointwts[ii])
        p2cluster.append(cd)
    for cl in cluster:
        mang = seeds[cl][-1];
        scl = seeds[cl]
        if len(cluster[cl]) > 10:
            std[cl] = np.percentile([angledist(xx[2], mang) for xx in cluster[cl]], 90)
            roadwidth[cl] = 1+5*np.std([geodist(scl,xx)*np.sin(anglebetweentwopoints(scl,xx)-scl[-1])  for xx in cluster[cl]])
            print(cl,scl,[(anglebetweentwopoints(scl,xx),scl[-1])  for xx in cluster[cl]])

def printclusters(seeds):
    with open('clusters_uic.txt', 'w') as fdist:
        for pp in seeds:
            fdist.write("%s %s %s\n" % (pp[0],pp[1],pp[2]))

def computeclusters(datapointwts,maxiteration,SEEDRADIUS,theta):
    datapoint = [(x[0], x[1], x[2]) for x in datapointwts];
    seeds = getseeds(datapoint, SEEDRADIUS,theta);
    oldcost = 100000000;
    for ss in range(maxiteration):
        nseeds,cost,avgspeed,pointsperseed = newmeans(datapointwts,seeds,theta)
        print(ss, cost)
        if (oldcost-cost)/cost<0.0001:
            break;
        seeds = nseeds;
        oldcost = cost;
    for ii in range(1):
        seeds, seedweight = splitclusters(datapointwts, seeds,theta);
    return(seeds)

def printedges(gedges, seeds,datapointwts,theta):
    maxspeed = [0 for xx in range(len(seeds))]
    cluster, p2cluster = point2cluster(datapointwts, seeds,theta);
    for cd in cluster:
        maxspeed[cd] = int(np.percentile([0] + [xx[3] for xx in cluster[cd]], 90))
    with open('edgesuic.txt', 'w') as fdist:
        for gg in gedges:
            fdist.write("%s %s %s %s %s %s %s %s\n" % (seeds[gg[0]][0],seeds[gg[0]][1],seeds[gg[0]][2],seeds[gg[1]][0],seeds[gg[1]][1],seeds[gg[1]][2], maxspeed[gg[0]], maxspeed[gg[1]]))

def getgeojson(gedges,seeds):
    inp = []
    for xx in gedges:
        ll1 = seeds[xx[0]]; ll2 = seeds[xx[1]];
        inp.append([(ll1[0],ll1[1]),(ll2[0],ll2[1])])
    with open('map0.geojson', 'w') as fdist:
        fdist.write(MultiLineString(inp))

def plotmap(seeds,gedges,datapointwts):
    plt.figure(figsize=(12, 8))  # in inches!
    ax = plt.gca()
    ax.get_xaxis().get_major_formatter().set_useOffset(False)
    ax.get_yaxis().get_major_formatter().set_useOffset(False)
    seedswedge = set([gg[0] for gg in gedges]+[gg[1] for gg in gedges])
    seedslat, seedslon = [seeds[xx][0] for xx in seedswedge], [seeds[xx][1] for xx in seedswedge]
    plt.scatter(seedslat, seedslon, marker='.', color='g', s=20)  # pointsperseed)
#    seedslat, seedslon = [xx[0] for xx in datapointwts], [xx[1] for xx in datapointwts]
#    plt.scatter(seedslat, seedslon, marker='.', color='k', s=5)  # pointsperseed)
    segs = [];colsegs = []
    for gg in gedges:
        (n1,n2) = gg;
        segs.append(((seeds[n1][0], seeds[n1][1]), (seeds[n2][0], seeds[n2][1])))
        colsegs.append((0,0,1))
#    for xx in seedswedge:
#        gg = seeds[xx]; arl = 0.0001;
#        segs.append(((gg[0], gg[1]), (gg[0] - arl * np.sin(np.pi / 180 * gg[2]),gg[1] - 0.9999 * arl * np.cos(np.pi / 180 * gg[2]) )))
#        colsegs.append((0, 0, 1))
    ln_coll = matplotlib.collections.LineCollection(segs, colors=colsegs)# (rr, 0, 1 - rr))
    ax = plt.gca()
    ax.add_collection(ln_coll)
    plt.draw()
    plt.xlabel('lon')
    plt.ylabel('lat')
    plt.show()
#    print('total seeds:',  len(seeds), 'active seeds: ', len(seedswedge), 'total edges: ',len(gedges), time.time() - start)

def readseeds():
    seeds= [];
    with open('./clusters_uic.txt','rb') as f:
        for line in f:
            line = line[:-1].decode('ascii', 'ignore')
            zz = line.split(" ")
            seeds.append((float(zz[0]),float(zz[1]),float(zz[2])))
    return(seeds)
