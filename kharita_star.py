"""
author: sofiane
Create the road network by merging trajectories.

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


if __name__ == '__main__':
	# Default parameters
	RADIUS_METER = 25
	SAMPLING_DISTANCE = 20 # sparsification of the edges.
	HEADING_ANGLE_TOLERANCE = 100
	FILE_CODE = 'data_uic'
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
	#geodist = geopy.distance.VincentyDistance(meters=RADIUS_METER)
	clusters = []
	cluster_kdtree = None
	roadnet = nx.DiGraph()
	total_points = 0
	p_X = []
	p_Y = []
	starting_time = datetime.datetime.now()
	trajectories = create_trajectories(INPUT_FILE_NAME= '%s/%s.csv' % (DATA_PATH, FILE_CODE), waiting_threshold=21)

	starting_time = datetime.datetime.now()
	for i, trajectory in enumerate(trajectories[:-1]):
		sys.stdout.write('\rprocessing trajectory: %s / %s' % (i,len(trajectories)))
		sys.stdout.flush()
		update_cluster_index = False
		prev_cluster = -1
		current_cluster = -1
		first_edge = True
		for point in trajectory:
			p_X.append(point.lon)
			p_Y.append(point.lat)
			# very first case: enter only once
			if len(clusters) == 0:
				# create a new cluster
				new_cluster = Cluster(cid=len(clusters), nb_points=1, last_seen=point.timestamp, lat=point.lat,
									  lon=point.lon, angle=point.angle)
				clusters.append(new_cluster)
				roadnet.add_node(new_cluster.cid)
				prev_cluster = new_cluster.cid  # all I need is the index of the new cluster
				# recompute the cluster index
				cluster_kdtree = cKDTree([c.get_lonlat() for c in clusters])
				continue
			# if there's a cluster within x meters and y angle: add to. Else: create new cluster
			close_clusters_indices = [clu_index for clu_index in cluster_kdtree.query_ball_point(x=point.get_lonlat(), r=RADIUS_DEGREE, p=2)
									  if math.fabs(diffangles(point.angle, clusters[clu_index].angle)) <= HEADING_ANGLE_TOLERANCE ]

			if len(close_clusters_indices) == 0:
				# create a new cluster
				new_cluster = Cluster(cid=len(clusters), nb_points=1, last_seen=point.timestamp, lat=point.lat, lon=point.lon, angle=point.angle)
				clusters.append(new_cluster)
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
				prev_cluster = current_cluster
				continue

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
					continue
				else:
					# identify the cluster to which the intermediate cluster belongs
					PT = geopy.Point(pt.get_coordinates())
					close_clusters_distances = [
						geopy.distance.distance(PT, geopy.Point(clusters[clu_index].get_coordinates())).meters for clu_index
						in close_clusters_indices]
					closest_cluster_indx = close_clusters_indices[close_clusters_distances.index(min(close_clusters_distances))]
					intermediate_cluster_ids.append(closest_cluster_indx)

			# For each element is segment: if ==-1 create new cluster and link to it, else link to the corresponding cluster
			prev_path_point = prev_cluster
			for idx, inter_clus_id in enumerate(intermediate_cluster_ids):
				if inter_clus_id == -1:
					n_cluster_point = intermediate_clusters[idx]
					# create a new cluster
					new_cluster = Cluster(cid=len(clusters), nb_points=1, last_seen=point.timestamp, lat=n_cluster_point.lat,
										  lon=n_cluster_point.lon, angle=n_cluster_point.angle)
					clusters.append(new_cluster)
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
					roadnet.add_edge(prev_path_point, new_cluster.cid)
					prev_path_point = new_cluster.cid
				else:
					roadnet.add_edge(prev_path_point, inter_clus_id)
					prev_path_point = inter_clus_id
					clusters[inter_clus_id].add(intermediate_clusters[idx])
			if len(intermediate_cluster_ids) == 0 or intermediate_cluster_ids[-1] != current_cluster:
				roadnet.add_edge(prev_path_point, current_cluster)
			prev_cluster = current_cluster
		if update_cluster_index:
			cluster_kdtree = cKDTree([c.get_lonlat() for c in clusters])
	exec_time = datetime.datetime.now() - starting_time
	with open('%s/%s_edges.txt' % (DATA_PATH, FILE_CODE), 'w') as fout:
		for s, t in roadnet.edges():
			fout.write('%s,%s\n%s,%s\n\n' % (clusters[s].lon, clusters[s].lat, clusters[t].lon, clusters[t].lat))
	print 'Graph generated in %s seconds' % exec_time.seconds
	if drawmap:
		from matplotlib import collections as mc, pyplot as plt
		lines = [[clusters[s].get_lonlat(), clusters[t].get_lonlat()] for s, t in roadnet.edges()]
		lc = mc.LineCollection(lines)
		fig, ax = plt.subplots()
		ax.add_collection(lc)
		ax.autoscale()
		ax.margins(0.1)
		plt.show()
