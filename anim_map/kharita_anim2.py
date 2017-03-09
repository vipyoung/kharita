"""
author: sofiane
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
import copy
from collections import defaultdict
from random import shuffle
import operator


def create_proxy(label, marker='s'):
	line = plt.Line2D((0, 1), (0, 0), color=label, marker=marker, linestyle='')
	return line

def generate_image(list_of_edges, edge_weight, nbr, roardnet, clusters, nb_traj, osm, fig, ax, timestamp, lonlat_to_cid):
	path_animation = '/home/sofiane/projects/2017/kharita/animation_bbx'
	print 'generating image:', nbr
	lines = [[s, t] for s, t in list_of_edges]
	colors = ['green' if edge_weight[i] == 1 else 'white' for i in range(len(lines))]
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


if __name__ == '__main__':

	# For the animation:
	lonlat_to_clusterid = dict()
	list_of_edges = []
	edge_weight = []
	frame_nb = 0
	fig, ax = plt.subplots(facecolor='black')
	osm_rn = nx.read_gpickle('osm_qmic_roads.gpickle')
	osm_lines = lines = [[s, t] for s, t in osm_rn.edges()]
	osm_lc = mc.LineCollection(osm_lines, colors=['0.2' for _ in range(len(osm_lines))])

	# Default parameters
	RADIUS_METER = 25
	SAMPLING_DISTANCE = 10 # sparsification of the edges.
	HEADING_ANGLE_TOLERANCE = 45
	# MAX_PATH_LEN = 20
	# MAX_PATH_DISTANCE_FACTOR = 2.7
	FILE_CODE = 'data_bbox'
	DATA_PATH = '../data'

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
	clusters = []
	cluster_kdtree = None
	roadnet = nx.DiGraph()
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
	update_cluster_index = False
	prev_cluster = -1
	current_cluster = -1
	first_edge = True

	for point in gps_point_stream:
		if point.timestamp < datetime.datetime.strptime('2015-11-05 22:00:00', '%Y-%m-%d %H:%M:%S'):
			continue
		traj_id = point.traj_id
		prev_cluster = building_trajectories.get(traj_id, -1)
		p_X.append(point.lon)
		p_Y.append(point.lat)
		# very first case: enter only once
		if len(clusters) == 0:
			# create a new cluster
			new_cluster = Cluster(cid=len(clusters), nb_points=1, last_seen=point.timestamp, lat=point.lat,
			                      lon=point.lon, angle=point.angle)
			clusters.append(new_cluster)
			lonlat_to_clusterid[new_cluster.get_lonlat()] = new_cluster.cid
			roadnet.add_node(new_cluster.cid)
			prev_cluster = new_cluster.cid  # all I need is the index of the new cluster
			# recompute the cluster index
			cluster_kdtree = cKDTree([c.get_lonlat() for c in clusters])

			# add the point to its trajectory, to keep track of the last point in trajectory
			building_trajectories[traj_id] = new_cluster.cid
			continue
		# if there's a cluster within x meters and y angle: add to. Else: create new cluster
		close_clusters_indices = [clu_index for clu_index in cluster_kdtree.query_ball_point(x=point.get_lonlat(), r=RADIUS_DEGREE, p=2)
		                          if math.fabs(diffangles(point.angle, clusters[clu_index].angle)) <= HEADING_ANGLE_TOLERANCE ]

		if len(close_clusters_indices) == 0:
			# create a new cluster
			new_cluster = Cluster(cid=len(clusters), nb_points=1, last_seen=point.timestamp, lat=point.lat, lon=point.lon, angle=point.angle)
			clusters.append(new_cluster)
			lonlat_to_clusterid[new_cluster.get_lonlat()] = new_cluster.cid
			roadnet.add_node(new_cluster.cid)
			current_cluster = new_cluster.cid
			# recompute the cluster index
			update_cluster_index = True
		else:
			# add the point to the cluster
			pt = geopy.Point(point.get_coordinates())
			close_clusters_distances = [geopy.distance.distance(pt, geopy.Point(clusters[clu_index].get_coordinates())).meters
			                            for clu_index in close_clusters_indices]
			closest_cluster_indx = close_clusters_indices[close_clusters_distances.index(min(close_clusters_distances))]
			clusters[closest_cluster_indx].add(point)
			current_cluster = closest_cluster_indx
		# Adding the edge:
		if prev_cluster == -1:
			#prev_cluster = current_cluster
			building_trajectories[traj_id] = current_cluster
			continue

		#edge = [clusters[prev_cluster].get_coordinates(), clusters[current_cluster].get_coordinates()]
		edge = [clusters[prev_cluster], clusters[current_cluster]]
		intermediate_clusters = partition_edge(edge, distance_interval=SAMPLING_DISTANCE)

		# Check if the newly created points belong to any existing cluster:
		intermediate_cluster_ids = []
		for pt in intermediate_clusters:
			close_clusters_indices = [clu_index for clu_index in
			                          cluster_kdtree.query_ball_point(x=pt.get_lonlat(), r=RADIUS_DEGREE, p=2)
			                          if math.fabs(diffangles(pt.angle, clusters[clu_index].angle)) <= HEADING_ANGLE_TOLERANCE]

			if len(close_clusters_indices) == 0:
				intermediate_cluster_ids.append(-1)
			else:
				# identify the cluster to which the intermediate cluster belongs
				PT = geopy.Point(pt.get_coordinates())
				close_clusters_distances = [
					geopy.distance.distance(PT, geopy.Point(clusters[clu_index].get_coordinates())).meters for clu_index
					in close_clusters_indices]
				closest_cluster_indx = close_clusters_indices[close_clusters_distances.index(min(close_clusters_distances))]
				intermediate_cluster_ids.append(closest_cluster_indx)

		# For each fictional point in segment: if ==-1 create new cluster and link to it, else link to the corresponding cluster
		prev_path_point = prev_cluster
		for idx, inter_clus_id in enumerate(intermediate_cluster_ids):
			if inter_clus_id == -1:
				n_cluster_point = intermediate_clusters[idx]
				# create a new cluster
				new_cluster = Cluster(cid=len(clusters), nb_points=1, last_seen=point.timestamp, lat=n_cluster_point.lat,
				                      lon=n_cluster_point.lon, angle=n_cluster_point.angle)
				clusters.append(new_cluster)
				lonlat_to_clusterid[new_cluster.get_lonlat()] = new_cluster.cid
				roadnet.add_node(new_cluster.cid)
				# recompute the cluster index
				update_cluster_index = True
				# create the actual edge:
				if math.fabs(diffangles(clusters[prev_path_point].angle, new_cluster.angle)) > HEADING_ANGLE_TOLERANCE \
					or math.fabs(diffangles(vector_direction_re_north(clusters[prev_path_point], new_cluster),
					                         clusters[prev_path_point].angle )) > HEADING_ANGLE_TOLERANCE:
					prev_path_point = new_cluster.cid
					continue
				# if satisfy_path_condition_distance(prev_path_point, new_cluster.cid, roadnet, clusters, alpha=1.2):
				if (prev_path_point, new_cluster.cid) not in list_of_edges:
					list_of_edges.append([clusters[prev_path_point].get_lonlat(), clusters[new_cluster.cid].get_lonlat()])
					roadnet.add_edge(prev_path_point, new_cluster.cid)
					edge_weight.append(1)
				else:
					edge_weight[list_of_edges.index([clusters[prev_path_point].get_lonlat(), clusters[new_cluster.cid].get_lonlat()])] += 1
				prev_path_point = new_cluster.cid
			else:
				if (prev_path_point, inter_clus_id) not in roadnet.edges():
					list_of_edges.append([clusters[prev_path_point].get_lonlat(), clusters[inter_clus_id].get_lonlat()])
					edge_weight.append(1)
					roadnet.add_edge(prev_path_point, inter_clus_id)
				else:
					edge_weight[list_of_edges.index([clusters[prev_path_point].get_lonlat(), clusters[inter_clus_id].get_lonlat()])] += 1
				prev_path_point = inter_clus_id
				clusters[inter_clus_id].add(intermediate_clusters[idx])
		if (len(intermediate_cluster_ids) == 0 or intermediate_cluster_ids[-1] != current_cluster) and \
			((prev_path_point, current_cluster) not in roadnet.edges()):
			list_of_edges.append([clusters[prev_path_point].get_lonlat(), clusters[current_cluster].get_lonlat()])
			edge_weight.append(1)
			roadnet.add_edge(prev_path_point, current_cluster)
		elif (prev_path_point, current_cluster) in roadnet.edges():
			edge_weight[list_of_edges.index([clusters[prev_path_point].get_lonlat(), clusters[current_cluster].get_lonlat()])] += 1

		building_trajectories[traj_id] = current_cluster
		# update the index
		if update_cluster_index:
			cluster_kdtree = cKDTree([c.get_lonlat() for c in clusters])

		generate_image(list_of_edges, edge_weight, frame_nb, roadnet, clusters, nb_traj=i, osm=osm_lc, fig=fig, ax=ax, lonlat_to_cid=lonlat_to_clusterid, timestamp=point.timestamp)
		frame_nb += 1
	exec_time = datetime.datetime.now() - starting_time

	# with open('%s/%s_edges.txt' % (DATA_PATH, FILE_CODE), 'w') as fout:
	# 	for s, t in roadnet.edges():
	# 		fout.write('%s,%s\n%s,%s\n\n' % (clusters[s].lon, clusters[s].lat, clusters[t].lon, clusters[t].lat))
	# print 'Graph generated in %s seconds' % exec_time.seconds
	# if drawmap:
	# 	from matplotlib import collections as mc, pyplot as plt
	# 	lines = [[clusters[s].get_lonlat(), clusters[t].get_lonlat()] for s, t in roadnet.edges()]
	# 	lc = mc.LineCollection(lines)
	# 	fig, ax = plt.subplots()
	# 	ax.add_collection(lc)
	# 	ax.autoscale()
	# 	ax.margins(0.1)
	# 	plt.show()