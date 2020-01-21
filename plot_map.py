"""
Author: Sofiane
Date: Jan 21, 2020

Sample code on how to plot the output map of Kharita Star (online).
NOTE: THIS DOESN'T WORK WITH KHARITA OFFLINE.
"""

from matplotlib import collections as mc, pyplot as plt 

edges = []
with open('data/data_uic_edges.txt') as f:
    lines = f.readlines()
    i = 0
    while i  <  len(lines) - 3:
        edges.append([lines[i].strip().split(','),
            lines[i+1].strip().split(',')])
        i += 3

lc = mc.LineCollection(edges)
fig, ax = plt.subplots()
ax.add_collection(lc)
ax.autoscale()
ax.margins(0.1)
plt.savefig('figs/uic_map_new.png', format='PNG')
plt.show()

