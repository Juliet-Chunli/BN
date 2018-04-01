#!/usr/bin/python3.6
# Needs to include DFS function to get the order of elements.
import numpy as np
from supportFunc import *


class diGraph(object):
	def __init__(self):
		self.nodes = set()
		self.order = []
		self.visited = set()
		self.root = None
		self.nameToNode = {}
		self.timeDepResult = {}   # hold the time dependent results for nodes interested, key node.uniqueName, val List of result
		self.timeSeq = []
		self.topoSorted = False # Use a flag to note whether the graph is sorted or not. If already sorted and no change in structure, no need to resort. 
		# self.timeDep = False  # Need to set self.timeDep to True if the graph is time dependent
		self.timeIterNodes = []   # nodes that are time dependent and need to be updated. Only prev nodes because they are the one need to be updated and we can find corresponding value by node.next pointer.  


	def addNode(self, newNodeInput):
		if type(newNodeInput) is list:
			for newNode in newNodeInput:
				if newNode not in self.nodes:
					self.nodes.add(newNode)
					unit = newNode.unitInfo
					if unit == None:
						self.nameToNode[newNode.name] = newNode
					else:
						if newNode.name not in self.nameToNode:
							self.nameToNode[newNode.name] = {}
							self.nameToNode[newNode.name][unit] = newNode
						else:
							self.nameToNode[newNode.name][unit] = newNode

		else:
			newNode = newNodeInput
			if newNode not in self.nodes:
				self.nodes.add(newNode)
				unit = newNode.unitInfo
				if unit == None:
					self.nameToNode[newNode.name] = newNode
				else:
					if newNode.name not in self.nameToNode:
						self.nameToNode[newNode.name] = {}
						self.nameToNode[newNode.name][unit] = newNode
					else:
						self.nameToNode[newNode.name][unit] = newNode



	def topologicalSort(self):
		def dfsVisit(self, v):
			self.visited.add(v)
			for n in v.child:
				if n not in self.visited:
					dfsVisit(self, n)
			self.order.append(v)
			

		for node in self.nodes:
			if node not in self.visited:
				dfsVisit(self, node)

		self.visited = set()  ## Have to clear previous values here 

		self.order.reverse()

		self.topoSorted = True

		# for node in self.order:
		# 	print(node.uniqueName)



	def solve(self, t = None, deltaT = None, truncateProb = None, truncETmethod = None, reNormalize = False, hashDigit = 2, normETtrunc = False):  ########### add comment on reNormalize

		# node.parentSetDist is updated here. 


		if self.topoSorted == False:
			self.topologicalSort()



		# for node in self.order:
		# 	print(node.uniqueName)



		for node in self.order:
			# print('node.name in solve', node.name)
			if len(node.parent) == 0:  # no parent and also not initialized by dist. 
				if node.timeDep:
					node.dist = node.distT(t, deltaT)
					# print('computeGraph', node.name, node.dist)
					# print('node.distT in computeGraph', node.name, node.distT)
				# elif node.typ == 'tIterPrev':
				# 	if node.dist
				elif node.typ != 'tIterPrev':  # tIterPrev has no parent, but not computed here. They are computed by init and update
					if node.needSample:
						mean = node.mean
						median = node.median
						mu = np.log(median)
						std = (2*(np.log(mean) - np.log(median)))**0.5
						pf = np.random.lognormal(mu, std)
						node.dist[node.stateStrToNum['f']] = pf
						node.dist[node.stateStrToNum['s']] = 1 - pf
					else:
						node.dist[node.stateStrToNum['f']] = node.pFail
						node.dist[node.stateStrToNum['s']] = 1 - node.pFail

			elif node.treeType == 'ET':
				if node.ET == None:
					print('Event tree is not defined yet')
				else:
					if truncETmethod == 'incAllES':
					
						c = next(iter(node.child))   # only applicable to cont only have ET single parent. 
						dic = c.parentSingleDist[node] # dict for ET es to cont ES
						contESset = set()   # set of all cont ES, in num form
						ETcontESdic = {}
						# for key in dic:
						# 	oneContSet = set(dic[key])
						# 	contESset = contESset.union(oneContSet)
						# 	ETcontESdic[node.statesNumToStr[key]] = oneContSet



						contESset = set()   # set of all cont ES, in num form
						ETcontESdic = {}
						for key in dic:
							oneContSet = set(dic[key])
							contESset = contESset.union(oneContSet)
							ETcontESdic[node.statesNumToStr[key]] = oneContSet


						node.solveETrelatinoship(truncETmethod, normET = normETtrunc, ETcontESdic = ETcontESdic, contSet = contESset)
						
					elif truncETmethod == 'threshold':
						node.solveETrelatinoship(truncETmethod, normET = normETtrunc, truncateProb  = truncateProb)
					elif truncETmethod == None:
						node.solveETrelatinoship(truncETmethod)



			elif node.typ == 'sfn' or node.typ == 'multi':        
				# update node.parentSetDist together for two node types


				if len(node.condDist) == 0:
					node.updateCondDist() 
				normConst = 0
				for string in node.condDist:
					prob = 1

					# print(node.uniqueName, 'node.relateOrder', node.relateOrder)
					for (i, pName) in enumerate(node.relateOrder):
						if pName not in node.parentNameToNode:
							print('Parent name ' + pName + ' for node ' + node.uniqueName + ' is not defined correctly.')
						else:
							pNode = node.parentNameToNode[pName]


						# print(pNode.uniqueName, pNode.dist)
						# print(node.uniqueName, node.parentSetDist)

						pNum = string[hashDigit*i : hashDigit*(i + 1)]
						if pNum in pNode.dist:
							prob *= pNode.dist[pNum]
						else:   # if any pNum not in pNode.dist, then this whole string is not a valid string, so prob = 0
							prob = 0
							break

					node.parentSetDist[string] = prob
					normConst += prob
				if reNormalize:
					for key in node.parentSetDist:
						node.parentSetDist[key] /= normConst

				# compute node distribution given parentSetdist and condDist
				# 
				
				if node.typ == 'sfn':
					mean = 0
					Evar = 0
					EcondMean2 = 0
					
					for key in node.condDist:
						pParentSet = node.parentSetDist[key]
						if 'val' in node.condDist[key]:
							mean += node.condDist[key]['val'] * pParentSet
							EcondMean2 += node.condDist[key]['val']**2 * pParentSet
						else:
							mean += node.condDist[key]['mean'] * pParentSet
							EcondMean2 += node.condDist[key]['mean']**2 * pParentSet
							Evar += medianToVar(node.condDist[key]['mean'], node.condDist[key]['median']) * pParentSet
					if node.needSample:
						var = Evar + EcondMean2 - mean**2
						mu, std = meanVarToMuStd(mean, var)
						pFail = np.random.lognormal(mu, std)
					else:
						pFail = mean
					node.dist[node.stateStrToNum['f']] = pFail
					node.dist[node.stateStrToNum['s']] = 1 - pFail
					node.dist[node.stateStrToNum['n']] = 1
				else:   # compute node.dist if typ = 'multi' given node.condDist and node.parentSetDist[key]. This only considers the case when parent set gives a value, uncertainty part to be implemented.
					for key in node.condDist:
						pParentSet = node.parentSetDist[key]
						for selfState in node.condDist[key]:
							if selfState not in node.dist:
								node.dist[selfState] = node.condDist[key][selfState] * pParentSet
							else:
								node.dist[selfState] += node.condDist[key][selfState] * pParentSet

					# delete impossible states to control complexity
					if node.dropZeroInDist:
						node.dist = {x:y for x,y in node.dist.items() if y!=0.0}


			elif node.condDist != {}:
				for pNode in node.parent:
					node.parentSetDist[pNode.name] = pNode.dist
				node.dist = node.condDist(t, deltaT, node.parentSetDist)



			elif node.dist == {}:
				print('Relationship to parent set is not defined correctly.')


	def clear(self):
		for node in self.nodes:
			node.dist = {}
			node.parentSetDist = {}


	# add time dep nodes that need to be updated. Only include prev nodes.  

	def addPrevNode(self, inputIterNode):
		if type(inputIterNode) is list:
			for node in inputIterNode:
				self.timeIterNodes.append(node)
		else:
			self.timeIterNodes.append(inputIterNode)

	# initiate nodes at time 0. dic key node, val initiate dist. 
	def initT0(self, initDic):
		for node in initDic:
			node.dist = initDic[node]


	# Need to initiate first before computeInTime Node of interest is a list with all nodes of interest. Will update self.timeDepResult and self.timeSeq

	def computeInTime(self, startTime, endTime, timeStep, nodeOfInterest):
		if self.topoSorted == False:
			self.topologicalSort()

		# solve and memorize dist of node of interest at startTime

		# self.solve(t = startTime, deltaT = timeStep)
		# self.timeSeq.append(startTime)
		# for node in nodeOfInterest:
		# 	if node.uniqueName not in self.timeDepResult:
		# 		self.timeDepResult[node.uniqueName] = [node.dist]

		# Iteration by time

		for t in np.arange(startTime, endTime, timeStep):
			print('t in computGraph', t)

			self.timeSeq.append(t)
			self.solve(t = t, deltaT = timeStep)

			for node in nodeOfInterest:
				if node.uniqueName not in self.timeDepResult:
					self.timeDepResult[node.uniqueName] = [node.dist]
				else:
					self.timeDepResult[node.uniqueName].append(node.dist)


			for node in self.timeIterNodes:
				if type(node.dist) is dict:
					node.dist = node.next.dist.deepcopy()
				else:
					node.dist = node.next.dist
			






