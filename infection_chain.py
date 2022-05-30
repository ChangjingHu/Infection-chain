import h5py
import numpy as np
import dynetx as dn
import pandas as pd
import networkx as nx
from matplotlib import pyplot as plt
import csv
from networkx.generators import *

G = nx.Graph()
#construct a graph in order to record the edges between nodes
all_population = [i for i in range(1, 387583)]
#print(all_population)
G.add_nodes_from(all_population)
#all candidates can be considered as isolated nodes
data = pd.read_csv(r"C:\Users\CharlotteHoo\Desktop\internship\data_to_2020_3_8\test.csv", header=None)
#data_1 = pd.read_csv(r"C:\Users\CharlotteHoo\Desktop\internship\data_to_2020_3_8\results.csv", header=None)
#data = data_0.append(data_1)
data.columns = ["timestamp", "location_ids", "infector_ids", "infected_ids", "infection_ids", "location_spaces", "region_names"]
#organize data table
group = data.groupby('timestamp')
#divide infection events by date

dates = []
#list for time series
for key, value in group:
	dates.append(key)
#here "key" is date, and date is put in dates[]

g = dn.DynGraph(edge_removal=True)
#start a dynamic network which can change networks over time, and time step is one day here

chain_record = {}
#define a dictionary to record chain of infection

for i in range(len(dates)):
	subgroup = group.get_group(dates[i])
	#acquire the data of a certain day in chronological order
	seed = subgroup[subgroup["location_spaces"] == "infection_seed"]
	#when initial seed patient
	for j in seed["infected_ids"]:
		chain_record[j] = 1
	subgroup = subgroup[subgroup["location_spaces"] != "infection_seed"]
	#when not initial seed patient
	#g. add_nodes_from(all_population)
	for k in subgroup["infected_ids"]:
		index = subgroup[subgroup["infected_ids"] == k].index.tolist()[0]
		#get the index of infection event
		patient = subgroup[subgroup["infected_ids"] == k]
		#get a complete piece of data at the location of "infected_ids"
		#print(chain_record.get(patient.at[index, "infector_ids"]))
		if chain_record.get(patient.at[index, "infector_ids"]) is not None:
			chain_record[k] = chain_record.get(patient.at[index, "infector_ids"]) + 1
			#add 1 to the infector person's current chain location
		else:
			chain_record[k] = 1
			#if the new patient is not infected by other infected person, it is a seed patient
	#print(chain_record)
	G = nx.from_pandas_edgelist(subgroup, "infector_ids", "infected_ids")
	#add edges to the network graph used to record new rounds of infections, and G is an intermediate container
	#add new edges for each time step
	g.add_interactions_from(G.edges(), t=i)
	#add the newly added edges in G at time t to the dynamic graph g
	degrees = g.degree()
	nodes = g.nodes()
	n_color = np.asarray([degrees[n] for n in nodes])
	#set different colors based on each node's contribution to infection
	print("when ", dates[i])
	print("new addition:", G.edges)
	print("--------------------------------------------------------------------------------------")
#nx.draw(g, cmap=plt.get_cmap('rainbow'), node_color=n_color, node_size=10)
#draw network, node color is related with its degree and node size is 10 mm
#plt.show()
#print the network above
#print(chain_record)

with open("C:\\Users\\CharlotteHoo\\Desktop\\internship\\data_to_2020_3_8\\infection_chain.csv", 'w') as f:
	w = csv.DictWriter(f, chain_record.keys())
	w.writeheader()
	w.writerow(chain_record)
pd.read_csv("C:\\Users\\CharlotteHoo\\Desktop\\internship\\data_to_2020_3_8\\infection_chain.csv", header=None).T.to_csv("C:\\Users\\CharlotteHoo\\Desktop\\internship\\data_to_2020_3_8\\infection_chain.csv", header=False, index=False)
#transpose dateset
#data_1 = pd.read_csv("C:\\Users\\CharlotteHoo\\Desktop\\internship\\data_to_2020_3_8\\infection_chain.csv", header=None)


#Here I want to use different colors and lables to set nodes that have high degree centrality or betweenness centrality
#dynamic g cannot compute centrality and betweeness, there is a new static g graph using edges from dynamic g graph
g_edges = g.edges()
static_g = nx.Graph()
all_population = [i for i in range(1, 387583)]
static_g.add_nodes_from(all_population)
static_g.add_edges_from(g_edges)
#add all the nodes to the network
#when connecting two nodes, we need to ensure all the nodes with the id exist



#******************Computing centrality************************************************

#Degree centrality is a simple count of the total number of connections linked to a vertex. 
#It can be thought of as a kind of popularity measure
degCent = nx.degree_centrality(static_g)

#Descending order sorting centrality
degCent_sorted = dict(sorted(degCent.items(), key=lambda item: item[1], reverse=True))


#*************Computing betweeness*********************************************************

#it represents the degree to which nodes stand between each other. 
#For example, in a telecommunications network, a node with higher betweenness centrality would have more control over the network, because more information will pass through that node.
betCent = nx.betweenness_centrality(static_g, normalized=True, endpoints=True)



#Descending order sorting betweeness
betCent_sorted = dict(sorted(betCent.items(), key=lambda item: item[1], reverse=True))

print(betCent_sorted)

#Color for regular nodes
color_list = static_g.number_of_nodes() * ['lightsteelblue']

#Getting indices on top 10 nodes for each measure
N_top = 10
colors_top_10 = ['tab:orange','tab:blue','tab:green','lightsteelblue']

#Getting the first 10 nodes in centrality and betweeness
keys_deg_top = list(degCent_sorted)[0:N_top]
keys_bet_top = list(betCent_sorted)[0:N_top]


#**********************Computing centrality and betweeness intersection************************

inter_list = list(set(keys_deg_top) & set(keys_bet_top))
print(inter_list)

#Setting up color for all nodes
for i in inter_list:
	color_list.append(colors_top_10[2])
#	color_list[i] = colors_top_10[2]
print(keys_deg_top)
print(keys_bet_top)

#Setting up color for the first 10 nodes in centrality and betweeness
for i in range(N_top):
	if keys_deg_top[i] not in inter_list:
		color_list[keys_deg_top[i]] = colors_top_10[0]
	if keys_bet_top[i] not in inter_list:
		color_list[keys_bet_top[i]] = colors_top_10[1]

#Draw graph
static_g.remove_nodes_from(list(nx.isolates(static_g)))
nx.draw(static_g, with_labels=True, node_color=color_list)
#pos = nx.circular_layout(static_g)

#Setting up legend
labels = ['Top 10 deg cent','Top 10 bet cent','Top 10 deg and bet cent','no top 10']
for i in range(len(labels)):
	plt.scatter([],[],label=labels[i],color=colors_top_10[i])
plt.legend(loc='center')
plt.show()

'''
g_edges = g.edges()
static_g = nx.Graph()
all_population = [i for i in range(1, 387583)]
static_g.add_nodes_from(all_population)
static_g.add_edges_from(g_edges)
#add all the nodes to the network
# when connecting two nodes, we need to ensure all the nodes with the id exist

for v in list(static_g):
	if static_g.degree(v) == 0:
		static_g.remove_node(v)
#remove isolated nodes

mean = lambda l: sum(l)/len(l)
#define "mean"

for vv in list(static_g):
	if len([static_g.degree(x) for x in static_g.neighbors(vv)]) != 0:
		mean_deg_neighbors = mean([static_g.degree(x) for x in static_g.neighbors(vv)])
		if mean_deg_neighbors <= 1.0:
			static_g.remove_node(vv)
	else:
		static_g.remove_node(vv)
#Describe the influence of a node by the average number of people affected by its neighbors
#Here define mean_deg_neighbors should be greater than 1.0
#and remove persons that do not contribute to the infection chain too much
#and remain 15546 persons

for m in list(static_g):
	if static_g.degree(m) == 0:
		static_g.remove_node(m)
print(static_g.number_of_nodes())
#After the processing in the previous step, the isolated nodes are further removed
#and remain 14619 persons

degrees = static_g.degree()
nodes = static_g.nodes()
n_color = np.asarray([degrees[n] for n in nodes])
#decide node color according to its degree
nx.draw(static_g, cmap=plt.get_cmap('rainbow'), node_color=n_color, node_size=10)
#draw network, node color is related with its degree and node size is 10 mm
plt.show()

'''