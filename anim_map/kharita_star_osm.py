"""
author: sofiane abbar
Create the road network by merging trajectories.
Consume edges one at a time.
"""


import sys
import getopt
from methods import *
from matplotlib import pyplot as plt
from matplotlib import collections as mc
import operator
import copy

if __name__ == '__main__':
	# # For the animation:
	# lonlat_to_clusterid = dict()
	# list_of_edges = []
	# edge_weight = []
	# frame_nb = 0
	# fig, ax = plt.subplots(facecolor='black')
	# osm_rn = nx.read_gpickle('osm_qmic_roads.gpickle')
	# osm_lines = lines = [[s, t] for s, t in osm_rn.edges()]
	# osm_lc = mc.LineCollection(osm_lines, colors=['0.2' for _ in range(len(osm_lines))])

	# Default parameters
	RADIUS_METER = 35
	SAMPLING_DISTANCE = 20 # sparsification of the edges.
	HEADING_ANGLE_TOLERANCE = 25
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

	parameters = {'file_code': FILE_CODE,
	              'data_path': DATA_PATH,
	              'radius_meter': RADIUS_METER,
	              'radius_degree': RADIUS_DEGREE,
	              'sampling_distance': SAMPLING_DISTANCE,
	              'heading_angle': HEADING_ANGLE_TOLERANCE
	              }
	starting_time = datetime.datetime.now()
	# generate_image(list_of_edges, edge_weight, frame_nb, roadnet, clusters, nb_traj=i, osm=osm_lc, fig=fig, ax=ax,
	#               lonlat_to_cid=lonlat_to_clusterid, timestamp=datetime.datetime.now(), ground_map_edges=ground_map_edges)


	# TODO: create a mapmatching function in methods

	# matchedTraj, unmatchedTraj = mapMatching(roadnet, trajectories)
	roadnet, clusters, matched_osm_clusters, new_osm_clusters, dead_osm_clusters = kharitaStar(parameters=parameters)
	print roadnet.number_of_edges(), roadnet.number_of_nodes()
	draw_roadnet_id_colored(roadnet, clusters, matched_osm_clusters, new_osm_clusters, dead_osm_clusters)
	#newMap = mergeMaps(roadnet, inferredMap)




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
