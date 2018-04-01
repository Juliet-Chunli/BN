#!/usr/bin/python3.6
from event_tree import *
from node import *
import numpy as np
import os
from computeGraph import diGraph
import matplotlib.pyplot as plt
import warnings
from textwrap import wrap
warnings.filterwarnings('error')


def buildSolveGraph(truncP, ETname):
	# Analysis using ET8_steam_line_break_outside_cont.txt
	# unit 1 top events

	k3_u1  = nodeByDist('k3_u1 ', 'K-3', 1, 'sfn')
	tt1_u1 = nodeByDist('tt1_u1', 'TT-1', 1, 'sfn')
	ms2_u1 = nodeByDist('ms2_u1', 'MS-2', 1, 'sfn')
	l1_u1  = nodeByDist('l1_u1 ', 'L-1', 1, 'sfn')
	op2_u1 = nodeByDist('op2_u1', 'OP-2', 1, 'sfn')
	r3_u1  = nodeByDist('r3_u1 ', 'R-3', 1, 'sfn')
	cf2_u1 = nodeByDist('cf2_u1', 'CF-2', 1, 'sfn')
	cs2_u1 = nodeByDist('cs2_u1', 'CS-2', 1, 'sfn')
	na2_u1 = nodeByDist('na2_u1', 'NA-2', 1, 'sfn')

	TElist_u1  = [k3_u1, tt1_u1, ms2_u1, l1_u1, op2_u1, r3_u1, cf2_u1, cs2_u1, na2_u1]


	# unit 1 buses
	bus2_u1 = nodeByDist('bus2_u1', 'BUS-2', 1, 'sfn')
	bus3_u1 = nodeByDist('bus3_u1', 'BUS-3', 1, 'sfn')
	bus5_u1 = nodeByDist('bus5_u1', 'BUS-5', 1, 'sfn')
	bus6_u1 = nodeByDist('bus6_u1', 'BUS-6', 1, 'sfn')

	busList_u1 = [bus2_u1, bus3_u1, bus5_u1, bus6_u1]

	# unit 1 diesels
	edg21 = nodeByDist('edg21', 'EDG-21', None, 'sfn')
	edg22 = nodeByDist('edg22', 'EDG-22', None, 'sfn')
	edg23 = nodeByDist('edg23', 'EDG-23', None, 'sfn')

	edgList_u1 = [edg21, edg22, edg23]

	# unit 1 plant end state

	plant_u1   = nodeByDist('plant_u1', 'Plant_ES', 1, 'ET')

	cont_u1    = nodeByDist('cont_u1', 'Cont_ES', 1, 'multi', dropZero = False)

	# build TE and bus relationship

	bus5_u1.addParent(edg21)
	bus2_u1.addParent(edg22)
	bus3_u1.addParent(edg22)
	bus6_u1.addParent(edg23)
	for te in TElist_u1:
	    plant_u1.addParent(te)
	    for bus in busList_u1:
	        te.addParent(bus)




	g = diGraph()
	for node in TElist_u1:
		g.addNode(node)
	for node in busList_u1:
		g.addNode(node)
	for node in edgList_u1:
		g.addNode(node)

	g.addNode(plant_u1)
	g.addNode(cont_u1)


	nonCompNode = set([plant_u1, cont_u1])
	ETdataDir = os.path.join(os.getcwd(), 'ET_data_txt')
	for node in g.nodes:
		if node not in nonCompNode:
			# print(node.name)
			fileName = node.name + '.txt'
			file = os.path.join(ETdataDir, fileName)
			node.addCondRelation(file = file)

	plant_u1 = g.nameToNode['Plant_ES'][1]
	cont_u1 = g.nameToNode['Cont_ES'][1]

	

	ETstructDir = os.path.join(os.getcwd(), 'ET_structure_txt')



	ET_u1 = Event_tree(os.path.join(ETstructDir, ETname + '.txt'))

	plant_u1.eventTreeRelationship(ET_u1)







	# for node in busList_u1:
	# 	print(node.parentNameToNode)

	file = os.path.join(os.getcwd(), 'plantEStoRelease.txt')
	cont_u1.addParent(plant_u1)
	cont_u1.addSinglePcondRel(file, plant_u1)
	# print(cont_u1.parentSingleDist)

	# print(plant_u1.stateSpace)
	# 
	g.solve(truncateProb = truncP, normETtrunc = True)






	return cont_u1




def MC(n, truncP, ETname):
	MCresult = {}
	for i in range(n):
		nodeWant = buildSolveGraph(truncP, ETname)
		plant_u1 = nodeWant.parent[0]
		if i % 50 == 0:
			print('plant_u1 truncated ES \n')
			# print(len(plant_u1.dist))
			esName = list(plant_u1.statesNumToStr[key] for key in plant_u1.dist)
			print(esName)
			print(plant_u1.leftSeqNum)


		# esName = list(plant_u1.statesNumToStr[key] for key in plant_u1.dist)
		# if len(esName) > 1:
		# 	print('plant_u1 truncated ES \n')
		# 	print(esName)
		# 	print(plant_u1.leftSeqNum)

		if len(MCresult) == 0:
			for key in nodeWant.dist:
				MCresult[key] = [nodeWant.dist[key]]
		else:
			for key in nodeWant.dist:
				MCresult[key].append(nodeWant.dist[key])
	return (MCresult, nodeWant.statesNumToStr)

def plotMCunsort(MCresult, statNumToStr, title):
    medianList = []
    n = len(MCresult)
    x = np.arange(n)
    xticks = []
    i = 0
    for state, statData in MCresult.items():
        array = np.array(statData)
        median = np.median(array)
        medianList.append(median)
        xticks.append(statNumToStr[state])
        plotX = array.size * [i]
        plt.plot(plotX, array, marker = 'x', linestyle = 'None', markevery = 1)
        try:
        	plt.yscale('log')
        except Warning:
        	plt.yscale('linear')
        i += 1
            

        # plt.plot(x, statData, marker = 'x', linestyle="None", markevery =10)
                  

    plt.plot(x, medianList, color = 'black', linestyle = 'None', marker = 'o', label = 'Median')
    plt.xticks(x, xticks, rotation = 'vertical', fontsize = 6)
    plt.ylabel('Probability')
    plt.title('\n'.join(wrap(title, 60)))
    plt.grid()
    plt.legend(loc = 'best')
    plt.savefig( title + '.png', dpi = 300)
    plt.show()




ETname = 'ET7_loss_of_main_feedwater'

truncP = 1
n = 500

MCresult, statNumToStr = MC(n, truncP, ETname)
plotMCunsort(MCresult, statNumToStr, ETname + '-Containment_End_State_Truncate_Risk_Ratio_' + str(truncP))

