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
from methods import create_trajectories, diffangles, partition_edge, vector_direction_re_north, Cluster
from matplotlib import pyplot as plt
from matplotlib import collections as mc
import operator
import copy
import os

def create_proxy(label, marker='s'):
	line = plt.Line2D((0, 1), (0, 0), color=label, marker=marker, linestyle='')
	return line


def road_color(weight):
	if weight == 0:
		return '0.3'
	if weight == 1:
		return 'green'
	return 'white'

def road_color_regarding_ground(edge, weight, ground_map_edges):
	"""
	If the edge is new (not part of the ground map) paint it in green.
	else, for grounp map edges use two colors: gray for un-used segments, white for used ones.
	:param edge:
	:param weight:
	:param ground_map_edges:
	:return:
	"""

	if edge not in ground_map_edges:
		return 'green'
	if weight == 0:
		return '0.3'
	return 'white'



def generate_image(list_of_edges, edge_weight, nbr, roardnet, clusters, nb_traj, osm, fig, ax, timestamp, lonlat_to_cid,
                   ground_map_edges):

	if os.path.exists('/home/sofiane/projects/2017/kharita/animation_bbx_osm'):
		path_animation = '/home/sofiane/projects/2017/kharita/animation_bbx_osm'
	else:
		path_animation = '/home/sabbar/projects/2017/kharita/animation_bbx_osm'
	print 'generating image:', nbr
	lines = [[s, t] for s, t in list_of_edges]
	# colors based on weight
	# colors = [road_color(edge_weight[i]) for i in range(len(lines))]

	# colors based on whether the edge exists in the ground map or not.
	colors = [road_color_regarding_ground(edge, edge_weight[i], ground_map_edges) for i, edge in enumerate(list_of_edges)]
	for i, edge in enumerate(list_of_edges):
		s_delta = timestamp - clusters[lonlat_to_cid[edge[0]]].last_seen
		t_delta = timestamp - clusters[lonlat_to_cid[edge[1]]].last_seen
		if (s_delta.days*24*3600 + s_delta.seconds) > 3600 and (t_delta.days*24*3600 + t_delta.seconds) > 3600:
			colors[i] = 'red'

	lc = mc.LineCollection(lines, colors=colors, linewidths=2)
	fig, ax = plt.subplots(facecolor='black', figsize=(14, 10))
	# add OSM every 100 frames
	# if nbr % 100 == 0:
	# 	ax.add_collection(copy.copy(osm))
	ax.add_collection(lc)
	#plt.plot(t[0], t[1], marker=(3, 0, 90), markersize=10, linestyle='None')
	plt.plot(t[0], t[1], marker="*", markersize=10, color='red', linestyle='None')



	# # Intersections in different colors?
	# outdegree = roadnet.out_degree()
	# indegree = roadnet.out_degree()
	# intersections = set([n for n in outdegree if outdegree[n] > 1] + [n for n in indegree if indegree[n] > 1])
	# X, Y = [], []
	# for n in intersections:
	# 	X.append(clusters[n].lon)
	# 	Y.append(clusters[n].lat)
	# plt.scatter(X, Y, color='yellow')

	ax.text(0.05, 0.01, 'Time: %s' % (timestamp),
	        verticalalignment='bottom',
	        horizontalalignment='left',
	        transform=ax.transAxes,
	        color='white', fontsize=10)

	ax.text(0.70, 0.01, '# Edges: %s' % len(list_of_edges),
	        verticalalignment='bottom',
	        horizontalalignment='right',
	        transform=ax.transAxes,
	        color='white', fontsize=10)

	ax.text(0.95, 0.01, 'Animation: S. Abbar',
	        verticalalignment='bottom',
	        horizontalalignment='right',
	        transform=ax.transAxes,
	        color='white', fontsize=6)

	ax.autoscale()
	# ax.margins(0.1)
	plt.axis('off')

	# legends
	descriptions = ['Vehicles', 'New Road Seg.', 'Confirmed Road Seg.', 'Unused Road Seg.']
	# descriptions = ['Vehicles', 'New Road Seg.', 'Confirmed Road Seg.']
	labels = ['red', 'green', 'white', 'red']
	pers_markers = ['*', 's', 's', 's']
	proxies = [create_proxy(item, mark) for item, mark in zip(labels, pers_markers)]
	ax.legend(proxies, descriptions, fontsize=6, numpoints=1, markerscale=1, ncol=4, bbox_to_anchor=(0.8, -0.05))

	plt.savefig('%s/frame_%s.png' % (path_animation, nbr), format='PNG',
	            facecolor=fig.get_facecolor(), transparent=True, bbox_inches='tight')

	# ax.clear()
	# fig.clf()
	plt.close()


def crop_osm_graph(fname):
	max_lat = 25.302769999999999
	min_lat = 25.283760000000001
	max_lon = 51.479749499999997
	min_lon = 51.458219999999997
	# use this awk command: awk 'BEGIN {FS=" ";} {if ($1 < 51.479749499999997 && $1 > 51.458219999999997 && $2 < 25.302769999999999 && $2 > 25.283760000000001 && $4 < 51.479749499999997 && $4 > 51.458219999999997 && $5 < 25.302769999999999 && $5 > 25.283760000000001 ) print}' osmmapclusterangle.txt > osm_bbx.csv


def build_initial_graph_from_osm(fname):
	"""
	Build the OSM graph for a list of edges: source, target.
	:param fname:
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


if __name__ == '__main__':
	# For the animation:
	lonlat_to_clusterid = dict()
	list_of_edges = []
	edge_weight = []
	frame_nb = 0
	fig, ax = plt.subplots(facecolor='black')
	osm_rn = nx.read_gpickle('/home/sabbar/PycharmProjects/kharita/anim_map/osm_qmic_roads.gpickle')
	osm_lines = lines = [[s, t] for s, t in osm_rn.edges()]
	osm_lc = mc.LineCollection(osm_lines, colors=['0.2' for _ in range(len(osm_lines))])

	# Default parameters
	RADIUS_METER = 25
	SAMPLING_DISTANCE = 10 # sparsification of the edges.
	HEADING_ANGLE_TOLERANCE = 10
	# MAX_PATH_LEN = 20
	# MAX_PATH_DISTANCE_FACTOR = 2.7
	FILE_CODE = 'data_bbox'
	DATA_PATH = 'data'

	drawmap = False
	(opts, args) = getopt.getopt(sys.argv[1:], "f:m:p:r:s:a:d:h")
	for o, a in opts:
		if o == "-f":
			FILE_CODE = str(a)
		if o == "-p":
			DATA_PATH = str(a)
		if o == "-r":
			RADIUS_METER = int(a)
		if o == "-s":
			SAMPLING_DISTANCE = int(a)
		if o == "-a":
			HEADING_ANGLE_TOLERANCE = int(a)
		if o == "-d":
			drawmap = True
		if o == "-h":
			print "Usage: python sofa_map.py [-f <file_name>] [-p <file repository>] [-r <clustering_radius>] [-s <sampling_rate>] " \
			      "[-a <heading angle tolerance>] [-h <help>]\n"
			exit()

	RADIUS_DEGREE = RADIUS_METER * 10e-6
	total_points = 0
	p_X = []
	p_Y = []
	starting_time = datetime.datetime.now()
	trajectories = create_trajectories(INPUT_FILE_NAME='%s/%s.csv' % (DATA_PATH, FILE_CODE), waiting_threshold=21)
	# push points as they are seen
	print 'creating gps point stream / nb trajectories %s' % len(trajectories)
	building_trajectories = dict()
	gps_point_stream = []
	for i, trajectory in enumerate(trajectories):
		for point in trajectory:
			point.set_traj_id(i)
			gps_point_stream.append(point)
	print 'Sort the gps point stream by timestamp'
	gps_point_stream = sorted(gps_point_stream, key=operator.attrgetter('timestamp'))
	trajectories = None
	update_kd_tree = False
	prev_cluster = -1
	current_cluster = -1
	first_edge = True

	covered_osm_clusters = []
	dead_osm_clusters = []
	new_trajectory_clusters = []

	# ##################### Incrementality starts here! #################################
	# Read and prepare the existing map, assume it is coming from OSM.
	clusters, actual_rn, list_of_edges, edge_weight = build_initial_graph_from_osm(fname='/home/sabbar/PycharmProjects/kharita/data/osm_bbx.csv')
	original_osm_clusters_lonlats = [c.get_lonlat() for c in clusters]

	# X, Y =[], []
	# for clu in clusters:
	# 	X.append(clu.lon)
	# 	Y.append(clu.lat)
	# plt.scatter(X, Y, c='black')
	# plt.show()

	ground_map_edges = copy.copy(list_of_edges)
	cluster_kdtree = cKDTree([c.get_lonlat() for c in clusters])
	lonlat_to_clusterid = {c.get_lonlat():c.cid for c in clusters}
	roadnet = actual_rn
	generate_image(list_of_edges, edge_weight, frame_nb, roadnet, clusters, nb_traj=i, osm=osm_lc, fig=fig, ax=ax,
	               lonlat_to_cid=lonlat_to_clusterid, timestamp=datetime.datetime.now(), ground_map_edges=ground_map_edges)

	for point in gps_point_stream:
		if point.timestamp < datetime.datetime.strptime('2015-11-05 22:00:00', '%Y-%m-%d %H:%M:%S'):
			continue
		traj_id = point.traj_id
		prev_cluster = building_trajectories.get(traj_id, -1)
		p_X.append(point.lon)
		p_Y.append(point.lat)

		# if there's a cluster within x meters and y angle: add to. Else: create new cluster
		nearest_cluster_indices = [clu_index for clu_index in cluster_kdtree.query_ball_point(x=point.get_lonlat(), r=RADIUS_DEGREE, p=2)
		                           if math.fabs(diffangles(point.angle, clusters[clu_index].angle)) <= HEADING_ANGLE_TOLERANCE]

		# *****************
		# Cluster creation
		# *****************
		# TODO: be more conservative in creating clusters! Add something like a threshold, min number of cars, etc.
		if len(nearest_cluster_indices) == 0:
			# create a new cluster
			new_cluster = Cluster(cid=len(clusters), nb_points=1, last_seen=point.timestamp, lat=point.lat, lon=point.lon, angle=point.angle)
			clusters.append(new_cluster)
			lonlat_to_clusterid[new_cluster.get_lonlat()] = new_cluster.cid
			roadnet.add_node(new_cluster.cid)
			current_cluster = new_cluster.cid
			# recompute the cluster index
			update_kd_tree = True
			new_trajectory_clusters.append(new_cluster.get_lonlat())
		else:
			# add the point to the cluster
			pt = geopy.Point(point.get_coordinates())
			close_clusters_distances = [geopy.distance.distance(pt, geopy.Point(clusters[clu_index].get_coordinates())).meters
			                            for clu_index in nearest_cluster_indices]
			closest_cluster_indx = nearest_cluster_indices[close_clusters_distances.index(min(close_clusters_distances))]
			clusters[closest_cluster_indx].add(point)
			current_cluster = closest_cluster_indx
			if clusters[closest_cluster_indx].get_lonlat() in original_osm_clusters_lonlats:
				covered_osm_clusters.append(clusters[closest_cluster_indx].get_lonlat())
		if update_kd_tree:
			cluster_kdtree = cKDTree([c.get_lonlat() for c in clusters])
			update_kd_tree = False
	print 'len new clusters:', len(new_trajectory_clusters)

	print 'covered osm clusters:', len(covered_osm_clusters)

	dead_osm_clusters = [c for c in original_osm_clusters_lonlats if c not in covered_osm_clusters]

	X, Y = [], []
	for pt in new_trajectory_clusters:
		X.append(pt[0])
		Y.append(pt[1])
	plt.scatter(X, Y, c='green')

	X, Y = [], []
	for pt in covered_osm_clusters:
		X.append(pt[0])
		Y.append(pt[1])
	plt.scatter(X, Y, c='0.8')

	# X, Y = [], []
	# for pt in dead_osm_clusters:
	# 	X.append(pt[0])
	# 	Y.append(pt[1])
	# plt.scatter(X, Y, c='red')

	plt.show()
		# TODO: deal with edges later
		# # *****************
		# # Edge creation
		# # *****************
		# # case of very first point in the trajectory (has no previous cluster.)
		# if prev_cluster == -1:
		# 	building_trajectories[traj_id] = current_cluster
		# 	continue
		#
		# edge = [clusters[prev_cluster], clusters[current_cluster]]
		# # TODO: I can add a condition on when to create fictional clusters. E.g., condition on angle diff (prev,curr)
		# intermediate_fictional_clusters = partition_edge(edge, distance_interval=SAMPLING_DISTANCE)
		#
		# # Check if the newly created points belong to any existing cluster:
		# intermediate_fictional_cluster_ids = []
		# for pt in intermediate_fictional_clusters:
		# 	nearest_cluster_indices = [clu_index for clu_index in
		# 	                           cluster_kdtree.query_ball_point(x=pt.get_lonlat(), r=RADIUS_DEGREE, p=2)
		# 	                           if math.fabs(diffangles(pt.angle, clusters[clu_index].angle)) <= HEADING_ANGLE_TOLERANCE]
		# 	if len(nearest_cluster_indices) == 0:
		# 		intermediate_fictional_cluster_ids.append(-1)
		# 	else:
		# 		# identify the cluster to which the intermediate cluster belongs
		# 		PT = geopy.Point(pt.get_coordinates())
		# 		close_clusters_distances = [
		# 			geopy.distance.distance(PT, geopy.Point(clusters[clu_index].get_coordinates())).meters for clu_index
		# 			in nearest_cluster_indices]
		# 		closest_cluster_indx = nearest_cluster_indices[close_clusters_distances.index(min(close_clusters_distances))]
		# 		intermediate_fictional_cluster_ids.append(closest_cluster_indx)
		#
		# # For each fictional point in segment: if ==-1 create new cluster and link to it, else link to the corresponding cluster
		# prev_path_point = prev_cluster
		# for idx, inter_clus_id in enumerate(intermediate_fictional_cluster_ids):
		# 	if inter_clus_id == -1:
		# 		n_cluster_point = intermediate_fictional_clusters[idx]
		# 		# create a new cluster
		# 		new_cluster = Cluster(cid=len(clusters), nb_points=1, last_seen=point.timestamp, lat=n_cluster_point.lat,
		# 		                      lon=n_cluster_point.lon, angle=n_cluster_point.angle)
		# 		clusters.append(new_cluster)
		# 		lonlat_to_clusterid[new_cluster.get_lonlat()] = new_cluster.cid
		# 		roadnet.add_node(new_cluster.cid)
		# 		# recompute the clusters kd-tree index
		# 		update_kd_tree = True
		# 		# create the actual edge: condition on angle differences only.
		# 		if math.fabs(diffangles(clusters[prev_path_point].angle, new_cluster.angle)) > HEADING_ANGLE_TOLERANCE \
		# 			or math.fabs(diffangles(vector_direction_re_north(clusters[prev_path_point], new_cluster),
		# 			                         clusters[prev_path_point].angle )) > HEADING_ANGLE_TOLERANCE:
		# 			prev_path_point = new_cluster.cid
		# 			continue
		# 		# if satisfy_path_condition_distance(prev_path_point, new_cluster.cid, roadnet, clusters, alpha=1.2):
		# 		if (prev_path_point, new_cluster.cid) not in list_of_edges:
		# 			list_of_edges.append([clusters[prev_path_point].get_lonlat(), clusters[new_cluster.cid].get_lonlat()])
		# 			roadnet.add_edge(prev_path_point, new_cluster.cid)
		# 			edge_weight.append(1)
		# 		else:
		# 			edge_weight[list_of_edges.index([clusters[prev_path_point].get_lonlat(), clusters[new_cluster.cid].get_lonlat()])] += 1
		# 		prev_path_point = new_cluster.cid
		# 	else:
		# 		if (prev_path_point, inter_clus_id) not in roadnet.edges():
		# 			list_of_edges.append([clusters[prev_path_point].get_lonlat(), clusters[inter_clus_id].get_lonlat()])
		# 			edge_weight.append(1)
		# 			roadnet.add_edge(prev_path_point, inter_clus_id)
		# 		else:
		# 			edge_weight[list_of_edges.index([clusters[prev_path_point].get_lonlat(), clusters[inter_clus_id].get_lonlat()])] += 1
		# 		prev_path_point = inter_clus_id
		# 		clusters[inter_clus_id].add(intermediate_fictional_clusters[idx])
		# if (len(intermediate_fictional_cluster_ids) == 0 or intermediate_fictional_cluster_ids[-1] != current_cluster) and \
		# 	((prev_path_point, current_cluster) not in roadnet.edges()):
		# 	list_of_edges.append([clusters[prev_path_point].get_lonlat(), clusters[current_cluster].get_lonlat()])
		# 	edge_weight.append(1)
		# 	roadnet.add_edge(prev_path_point, current_cluster)
		# elif (prev_path_point, current_cluster) in roadnet.edges():
		# 	edge_weight[list_of_edges.index([clusters[prev_path_point].get_lonlat(), clusters[current_cluster].get_lonlat()])] += 1
		#
		# building_trajectories[traj_id] = current_cluster
		# # update the index
		# if update_kd_tree:
		# 	cluster_kdtree = cKDTree([c.get_lonlat() for c in clusters])
		#	update_kd_tree = False

		# generate_image(list_of_edges, edge_weight, frame_nb, roadnet, clusters, nb_traj=i, osm=osm_lc, fig=fig,
		#                ax=ax, lonlat_to_cid=lonlat_to_clusterid, timestamp=point.timestamp, ground_map_edges=ground_map_edges)
		# frame_nb += 1
	# exec_time = datetime.datetime.now() - starting_time


	# with open('%s/%s_edges.txt' % (DATA_PATH, FILE_CODE), 'w') as fout:
	# 	for s, t in roadnet.edges():
	# 		fout.write('%s,%s\n%s,%s\n\n' % (clusters[s].lon, clusters[s].lat, clusters[t].lon, clusters[t].lat))
	# print 'Graph generated in %s seconds' % exec_time.seconds
	# if drawmap:
	# 	from matplotlib import collections as mc, pyplot as plt
	# 	lines = [[clusters[s].get_lonlat(), clusters[t].get_lonlat()] for s, t in roadnet.edges()]
	# 	lc = mc.LineCollection(lines)
	# 	fig, ax = plt.subplots()
	# 	ax.autoscale()
	# 	ax.margins(0.1)
	# 	plt.show()
	# 	ax.add_collection(lc)
