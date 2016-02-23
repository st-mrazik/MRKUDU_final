import networkx as nx
import numpy as nm
from gurobipy import *
import matplotlib.pyplot as py
import math as math
from random import randint
import re

# Charakteristiky vstupneho grafu:

# Edge:
# ['weight'] vaha hrany

# Node
# ['x'] x-ova suradnica
# ['y'] y-ova suradnica
# ['n'] pocet obsluznych vozidieli <=0 pre obsluhovane miesta
# ['start'] zaciatok casoveho okna vo vrchole
# ['end'] koniec casoveho okna vo vrchole

def readGraph(filePath):
	G = nx.Graph()
	# pridaj atributy vrcholom
	nx.set_node_attributes(G, 'x', 0)
	nx.set_node_attributes(G, 'y', 0)
	nx.set_node_attributes(G, 'n', 0)
	nx.set_node_attributes(G, 'start', 0)
	nx.set_node_attributes(G, 'end', 0)
	# pridaj atributy hranam
	nx.set_edge_attributes(G, 'weight', 0)
	f = open(filePath, 'r')
	s = f.readline()
	regex = re.compile('([0-9]+\s){6}')
	
	# preskoc riadky pred datami vrcholov
	while regex.match(s) is None:
		s = f.readline()
	# nacitaj data vrcholov
	while regex.match(s) is not None:
		data = s.split()
		actNode = int(data[0])
		G.add_node(actNode)
		G.node[actNode]['x'] = int(data[1])
		G.node[actNode]['y'] = int(data[2])
		G.node[actNode]['n'] = int(data[3])
		G.node[actNode]['start'] = int(data[4])
		G.node[actNode]['end'] = int(data[5])
		s = f.readline()
	regex = re.compile('([0-9]+\s){3}')
	# preskoc riadky pred datami hran
	while regex.match(s) is None:
		s = f.readline()
	# nacitaj data hran
	while regex.match(s) is not None:
		data = s.split()
		nodeFrom = int(data[0])
		nodeTo = int(data[1])
		G.add_edge(nodeFrom, nodeTo)
		G.edge[nodeFrom][nodeTo]['weight'] = int(data[2])
		s = f.readline()
	return G
	
def drawGraph(G):
	xcords = nx.get_node_attributes(G,'x')
	ycords = nx.get_node_attributes(G,'y')
	coordinates = dict([(k, [xcords[k], ycords[k]]) for k in xcords])
	ns = nx.get_node_attributes(G, 'n')
	nx.draw_networkx(G,pos=coordinates,node_color=[1 if int(ns[i])>0 else 0 for i in ns],cmap=py.cm.PuBu)
	py.show() 

def drawMultipleGraphs(listG):
	for G in listG:
		xcords = nx.get_node_attributes(G,'x')
		ycords = nx.get_node_attributes(G,'y')
		coordinates = dict([(k, [xcords[k], ycords[k]]) for k in xcords])
		ns = nx.get_node_attributes(G, 'n')
		nx.draw_networkx(G,pos=coordinates,node_color=[1 if int(ns[i])>0 else 0 for i in ns],cmap=py.cm.PuBu)	
	py.show() 

def doLowestWeightNodeClustering(G):
	sourceNodeIndexes = list()
	nodesAssigment = dict()
	for i in G.nodes():
		if G.node[i]['n'] > 0:
			sourceNodeIndexes.append(i)
	clientNodeIndexes = list(set(G.nodes()) - set(sourceNodeIndexes))
	for i in range(0, len(clientNodeIndexes)):
		minWeight = sys.maxsize
		minWeightNodeIndex = -1
		for j in range(0, len(sourceNodeIndexes)):
			if G.has_edge(clientNodeIndexes[i],sourceNodeIndexes[j]) and G.edge[clientNodeIndexes[i]][sourceNodeIndexes[j]]['weight'] < minWeight:
				minWeight = G.edge[clientNodeIndexes[i]][sourceNodeIndexes[j]]['weight']
				minWeightNodeIndex = j
		if minWeightNodeIndex != -1:
			if sourceNodeIndexes[minWeightNodeIndex] not in nodesAssigment:
				nodesAssigment[sourceNodeIndexes[minWeightNodeIndex]] = list()
			nodesAssigment[sourceNodeIndexes[minWeightNodeIndex]].append(clientNodeIndexes[i])
	rG = list()
	for i in nodesAssigment.keys():
		rG.append(nx.Graph())
		workIndex = len(rG)-1
		nodesIncluded = list()
		nodesIncluded.append(i)
		nodesIncluded += nodesAssigment[i]
		for j in nodesIncluded:
			rG[workIndex].add_node(j,G.node[j])
		for j in nodesIncluded:
			for k in nodesIncluded:
				if G.has_edge(j,k):
					rG[workIndex].add_edge(j,k,G.edge[j][k])
	return rG

# parametre: cesta k vstupnemu a vystupnemu suboru
def convertEuc2DtoGraph(tspFile, graphFile):
	f = open(tspFile, 'r')
	s = f.readline()
	isEuc2D = False
	regex = re.compile('EDGE_WEIGHT_TYPE.*EUC_2D')
	while s and regex.match(s) is None:
		s = f.readline()
	if regex.match(s) is not None:
		isEuc2D = True
	if not isEuc2D:
		raise TypeError('Coordinates of file are not Euc2D')
		return
	regex = re.compile('(\s*[0-9]+|([0-9]+e\+[0-9]+)\s*){3}')
	i = 1
	while s and regex.match(s) is None:
		s = f.readline()
	nodes = {}
	while s and regex.match(s) is not None:
		nds = s.split()
		nodes[int(nds[0])] = (int(float(nds[1])), int(float(nds[2])))
		s = f.readline()
	f.close()
	f = open(graphFile, 'w')
	f.write('Nodes:\n')
	f.write('name\tx\ty\tn\tstart\tend\n')
	for nodeKey in nodes.keys():
		f.write(str(nodeKey) + '\t' + str(nodes[nodeKey][0]) + '\t' + str(nodes[nodeKey][1]) + '\t'
		+ '0\t0\t' + str(sys.maxsize) + '\n')
	f.write('Edges:\n')
	f.write('name1\tname2\tweight\n')
	for i in nodes.keys():
		for j in nodes.keys():
			if i < j:
				f.write(str(i) + '\t' + str(j) + '\t' 
				+ str(
				math.floor(
				math.sqrt(
				math.pow(
					nodes[i][0]-nodes[j][0]
				,2)+
				math.pow(
					nodes[i][1]-nodes[j][1]
				,2)
				)
				)
				)
				+ '\n')
	f.close()				

# Parametre	
# G vstupny graf
# max_tau maximalna doba obsluhy (cela okruzna jazda)
# tau doba obsluhy v jednotlivych vrcholoch (doba obsluzenia v kazdom vrchole)
def solveSingle(G,max_tau,tau):
	model = Model("ipamod")
	
	V_Z = list()  		#Zberne miesta
	P = {}      		#Pocet vozidiel (na indexe zberneho miesta
	for i in G.nodes():	
		if int(G.node[i]['n']) > 0:
			V_Z.append(i)
	V = list(G.nodes())   		#vrcholy
	V_K = list(set(V) - set(V_Z)) 	#Klienti
	
	x = {}
	for i in V:
		for j in V:
			if i != j:
				x[i,j] = model.addVar(ub=1, vtype='B', name="x(%s,%s)"%(i,j))
	# slucky pre nevyuzite jazdy - spomaluju vypocet
	#for i in V_Z:
	#	x[i,i] = model.addVar(lb=0, vtype=GRB.INTEGER, name="u(%s)"%j)

	u = {}
	for j in G.nodes():
		u[j] = model.addVar(lb=0, ub=max_tau, vtype=GRB.INTEGER, name="u(%s)"%j)
	model.update()

	#prijazd P[j] aut do miesta v V_Z
	for j in V_Z:
		model.addConstr(quicksum(x[i,j] for i in V if i != j) == int(G.node[j]['n']), "Stlpec_j%s"%(j))
		# verzia zo sluckami
		# model.addConstr(quicksum(x[i,j] for i in V) == int(G.node[j]['n']), "Stlpec_j%s"%(j))

	#ku klinetovi prichadzam prave raz
	for j in V_K:
		model.addConstr(quicksum(x[i,j] for i in V if i != j) == 1, "SumaMiesta_j%s"%(j))
	#vyjazd P[j] aut z miesta v V_Z	
	for i in V_Z:
		model.addConstr(quicksum(x[i,j] for j in V if i != j) == int(G.node[i]['n']), "Stlpec_i%s"%(j))
		# verzia zo sluckami
		# model.addConstr(quicksum(x[i,j] for j in V) == int(G.node[i]['n']), "Stlpec_i%s"%(j))

	#od klienta odchadzam prave raz	
	for i in V_K:
		model.addConstr(quicksum(x[i,j] for j in V if i != j) == 1, "SumaMiestai_i%s"%(j))
	
	for i in V:
		for j in V_K:
			if i != j:
				model.addConstr(u[i] - u[j] + (max_tau + int(G.edge[i][j]['weight'])+ tau)*x[i,j] <= max_tau, "u_%s_%s"%(i,j))	
	for i in V_K:
		for k in V_Z:
			model.addConstr(int(G.edge[i][k]['weight'])*x[i,k] + u[i] + tau <= max_tau, "end_u_%s_%s"%(i,j))	

	for i in V_K:
		model.addConstr(u[i] >= int(G.node[i]['start']))	
	for i in V_K:
		model.addConstr(u[i] <= int(G.node[i]['end']))	
	
	for i in V_Z:
		model.addConstr(u[i] == 0)
	
	model.update()

	model.setObjective(quicksum(int(G.edge[i][j]['weight'])*x[i,j] for i in V for j in V if i != j),GRB.MINIMIZE)

	model.update()

	model.optimize()
	edges = []
	EPS = 1.e-6
	retG = nx.Graph()
	retG.add_nodes_from(G.nodes(data=True))
	nx.set_edge_attributes(retG, 'weight', 0)

	for (i,j) in x:
		if x[i,j].X > EPS:
			if i != j:
				retG.add_edge(i,j)
				retG.edge[i][j]['weight'] = G.edge[i][j]['weight']
				#edges.append((i,j))
	ret_u = []
	for i in V:
		ret_u.append((i,u[i].X))	
	return model.ObjVal,retG, ret_u

# Parametre	
# G vstupny graf
# max_tau maximalna doba obsluhy (cela okruzna jazda)
# tau doba obsluhy v jednotlivych vrcholoch (doba obsluzenia v kazdom vrchole)
def solveMultiple(listG,max_tau,tau):
	results = list()
	for G in listG:
		objVal, retG, ret_u = solveSingle(G,max_tau,tau)
		results.append((objVal, retG, ret_u))
	objVal = 0
	retG = nx.Graph()
	ret_u = list()
	for res in results:
		objVal += res[0]
		retG = nx.compose(retG, res[1])
		ret_u = ret_u + res[2]
	return objVal, retG, ret_u

