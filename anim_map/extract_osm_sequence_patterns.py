"""
author: sofiane abbar
Create the road network by merging trajectories.
Consume edges one at a time.
"""


import geopy
import math
import numpy as np
import sys
import getopt
import datetime
import networkx as nx
from scipy.spatial import cKDTree
from methods import create_trajectories, diffangles, partition_edge, vector_direction_re_north, Cluster, calculate_bearing
from matplotlib import pyplot as plt
from matplotlib import collections as mc
import operator
import copy
import os
from math import cos, sin
from sklearn.cluster import KMeans

def crop_osm_graph(fname):
	max_lat = 25.302769999999999
	min_lat = 25.283760000000001
	max_lon = 51.479749499999997
	min_lon = 51.458219999999997
	# use this awk command: awk 'BEGIN {FS=" ";} {if ($1 < 51.479749499999997 && $1 > 51.458219999999997 && $2 < 25.302769999999999 && $2 > 25.283760000000001 && $4 < 51.479749499999997 && $4 > 51.458219999999997 && $5 < 25.302769999999999 && $5 > 25.283760000000001 ) print}' osmmapclusterangle.txt > osm_bbx.csv


def build_initial_graph_from_osm(fname):
	"""
	Build the OSM graph for a list of edges: source, target.
	:param fname: file generated and provided by Rade
	:return:
	"""
	clusters = []
	clusters_latlon = []
	list_of_edges = []
	edge_weight = []
	osm_roadnet = nx.DiGraph()
	now_ts = datetime.datetime.now()
	with open(fname) as f:
		for line in f:
			slon, slat, sang, tlon, tlat, tang = map(float, line.strip().split(' '))
			if (slat, slon) not in clusters_latlon:
				new_cluster = Cluster(cid=len(clusters), nb_points=1, last_seen=now_ts, lat=slat, lon=slon, angle=sang)
				clusters.append(new_cluster)
				clusters_latlon.append((slat, slon))
			if (tlat, tlon) not in clusters_latlon:
				new_cluster = Cluster(cid=len(clusters), nb_points=1, last_seen=now_ts, lat=tlat, lon=tlon, angle=tang)
				clusters.append(new_cluster)
				clusters_latlon.append((tlat, tlon))
			osm_roadnet.add_edge(clusters_latlon.index((slat, slon)), clusters_latlon.index((tlat, tlon)))
			list_of_edges.append([(slon, slat), (tlon, tlat)])
			edge_weight.append(0)
	clusters_latlon = None
	return clusters, osm_roadnet, list_of_edges, edge_weight


def paths_of_length_3(g):
	'''
	Return all simple paths of length l in a graph g
	:param g: directed graph
	:return: list of paths.
	'''
	paths = set()
	for node in  g.nodes_iter():
		neis = set([node])
		neis_1 = g.neighbors(node)
		for nei1 in neis_1:
			if nei1 == node:
				continue
			neis_2 = g.neighbors(nei1)
			for nei2 in neis_2:
				if nei2 == nei1 or nei2 == node:
					continue
				paths.add(tuple([node, nei1, nei2]))
	return paths


def path_to_vec(clusters, p):
	s, m, e = p
	startpoint = geopy.Point(clusters[s].get_coordinates())
	middlepoint = geopy.Point(clusters[m].get_coordinates())
	endpoint = geopy.Point(clusters[e].get_coordinates())
	d1 = geopy.distance.distance(startpoint, middlepoint).meters
	d2 = geopy.distance.distance(middlepoint, endpoint).meters
	bear_1 = calculate_bearing(clusters[s].lat, clusters[s].lon, clusters[m].lat, clusters[m].lon)
	bear_2 = calculate_bearing(clusters[m].lat, clusters[m].lon, clusters[e].lat, clusters[e].lon)

	# TODO: make sure this diffangles function does exactly what it's meant to do
	ang = diffangles(bear_1, bear_2)
	return [d1, d2, ang]


def path_vec_to_plot(path_vec):
	d1, d2, ang = path_vec
	startpoint = (0,0)
	middlepoint = (d1, 0)
	endpoint = (d1+ (cos(ang)*d2), 0+ sin(ang)*d2)
	return (startpoint, middlepoint, endpoint)


if __name__ == '__main__':

	# Default parameters
	RADIUS_METER = 25
	SAMPLING_DISTANCE = 10 # sparsification of the edges.
	HEADING_ANGLE_TOLERANCE = 10
	# MAX_PATH_LEN = 20
	# MAX_PATH_DISTANCE_FACTOR = 2.7
	FILE_CODE = 'data_bbox'
	DATA_PATH = 'data'

	# Read and prepare the existing map, assume it is coming from OSM.
	clusters, actual_rn, list_of_edges, edge_weight = build_initial_graph_from_osm(fname='../data/osm_bbx.csv')
	original_osm_clusters_lonlats = [c.get_lonlat() for c in clusters]

	print 'generate all path sequences'
	paths_3 = list(paths_of_length_3(actual_rn))
	path_vecs = []
	print 'nb paths:', len(paths_3)
	for path in paths_3:
		path_vecs.append(path_to_vec(clusters, path))

	kmeans = KMeans(n_clusters=9, random_state=0).fit(path_vecs)

	fig = plt.figure(figsize=(18, 20))

	for i, center in enumerate(kmeans.cluster_centers_):
		print center
		s, m, e = path_vec_to_plot(center)
		ax1 = fig.add_subplot(3,3,i+1)
		ax1.plot([s[0], m[0], e[0]], [s[1], m[1], e[1]], marker='o', linewidth=2)
		ax1.set_title('Cluster %s: %s paths' % (i+1, len([c for c in kmeans.labels_ if c == i])))
		ax1.set_ylim(ymin=-1)
	plt.show()