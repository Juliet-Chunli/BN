#!/usr/bin/python3.6
import numpy as np
import pandas as pd
from supportFunc import *

'''
Node by Distribution in Bayesian Nets


This node represents a distribution. It includes the information of all possible outcome and corresponding probability values for one single event. Advantage is completeness. It will compute next distribution given parent distribution. Disadvantage is complexity. Computational complexity for this method is very high. Since each node is represented by a distribution with full set of possible outcomes, we need to consider everything in order to get to next outcome. It's not too bad if we don't want to consider time dependency of the events, but it will be terrible if the analysis is time dependent. The reason is endstate space is much richer in time dependent situation. It's even worse if the time is continuous because in that case, we need to include all possible continuous time into our variable space.

Attributes:
    name: 
        name of the node, type str. Name is not unique, can't use name to identify node in Bayesian net. The name is used in event tree or fault tree to match the top event or basic event name. When we create a instance of the node type, the name of the instance should be unique to represent the node itself.
    uniqueName: 
        uniqueName of the node, includes information about the unit that this node belongs to. Could use same string as instance name while initiating the instance.
    stateSpace: 
        list of state space of the node. We only have two kinds of state spaces: ['s', 'f', 'n'] for events that can succeed, fail or ignored; state space for end states of event trees which is determined based on specific tree.
    stateSpace:
        state space code for the node itself. Two types: 'sfn' meaning this node is a component or system, only have succ, fail or ignore three values. 
    unitInfo: 
        unit the node belongs to. If it doesn't belong to any specific unit, value is None. Otherwise, it equals to the integer. 
    parent: 
        parent nodes of the node, type list
    parentSet: 
        set of parent node, type set
    parentNameToNode: 
        dictory of parent name to parent node. This is necessary since parent is represented by name in cond. dist, we need to image name to actual node.
    child: 
        child nodes of the node, type set

    relateOrder: 
        list of nodes stores the order of nodes used in parent state vector.
    treeType: 
        the relationship of node and parent is represented by event tree or fault tree or None. Event Tree: 'ET', Fault Tree: 'FT'.
    ET: 
        the event tree that describes the relationship to its parent. self.ET is not None when self.treeType == 'ET'. 
    mean: 
        mean of the node failure prob, if it has uncertainty
    median: 
        median of the node failure prob, if it has uncertainty
    needSample: 
        whether the node prob has uncertainty and needs sample. If any state of its parent has uncertainty, needs sample is true.
    pFail: 
        failure probability of the node. Only applicable to nodes without parent and don't need sample. 
    parentSetDist: 
        For 'sfn' and 'multi' type nodes: PMF of parent set as a whole. This will only be able to solve when len(relateOrder) != 0 and all parents have known state probability. The later one could be checked in computeGraph. Key parent set string in relateOrder, value prob. 

        For 'None' type nodes: key parent node name. The reason why it's not uniqueName is we may use the same logic for multiple units. Value distribution. 
    condDist: 
        conditional dist given parent set status. Key parent set string in relateOrder, value is also a dic. with three possible attributes: 'mean', 'median' and 'val'. When 'val' is not in condDist[key], meaning we need to sample to get cond. dist given mean and median.
    dist: 
        pmf of the node states, type dic, key states, value prob. Only solvable when condDist and parentSetDist are both set. len(condDist) != 0 and len(parentSetDist) != 0


    If treeType not 'ET' or 'FT', the distribution should be solved in computeGraph


    Possible useless attributes:

    relateToP: 
        correspondence relationship from parent state to self state, type dic, key parent states, value self stat. The parent states are represented by a random vector. It is a string in a specific order (relateOrder) that represents the state of each parent node.


Functions:
    addParent(self, newp): 
        add parent node (newp) to parent set, Append newp to parent list
    removeParent(self, delp): 
        if present, remove delp from parent set and parent list. This is expensive because remove from list is expensive, so try to avoid using it.
    addChild(self, newc): 
        add newc to child set
    removeChild(self, delc): 
        delete delc from child set if present.
    eventTreeRelationship(self, ET): 
        update the states, relateToP and relateOrder from an event tree ET. Have to set treeType to be 'ET' before using this function. Only use simple ES in this version because that's sufficient for representing the whole distribution.
    solveETrelatinoship(self, ET): 
        if treeType is 'ET', update top events probability based on parent distriutions and solve the ET. 
    addCondRelation(self, file): 
        input file that indicates the conditional distribution given parent set states and the order of parents. Will update self.relateToP and self.condDist. 

'''


class nodeByDist(object):
    def __init__(self, uniqueName, name, unitInfo, typ = None, dropZero = True):
        self.name = name
        self.uniqueName = uniqueName
        self.typ = typ           ####### 'sfn' means only succ / fail states, 'multi' means multiple states, 'tIterPrev' means time dep iter at previous time. 'timeDep' will set self.timeDep to True
        self.timeDep = False  ########## set timeDep to True if dist of the node is related to t. This only includes the node whose dist or condDist has t as input. Node that is time dependent but don't have t as direct input for dist are not included here. The nodes with timeDep == True will also have distT != None or distT. 
        if self.typ == 'timeDep':
            self.timeDep = True
        self.stateSpace = []   
        self.stateStrToNum = {}     ####################################
        self.statesNumToStr = {}    ##########################################
        self.hashDigits = None     ############################################
        if self.typ == 'sfn':
            self.stateSpace = ['s', 'f', 'n']
            self.hashStateLabel()
        elif self.typ == 'ET':
            self.treeType = 'ET'
        self.unitInfo = unitInfo
        self.dropZeroInDist = dropZero 
        self.parent = [] #list of p nodes
        self.parentSet = set()
        self.parentNameToNode = {}
        self.child = set()
        self.next = None # if the node has time dependent, we need to compute the value iteratively. self.next is the same node next in time. 
        self.prev = None # if the node has time dependent, we need to compute the value iteratively. self.pref is the same node prev in time. 
        # self.relateToP = {}
        self.relateOrder = []   # list of p.name, not node
        self.treeType = None
        self.ET = None
        self.dist = {}
        self.mean = None
        self.median = None
        self.needSample = False
        self.pFail = None
        self.parentSingleDist = {}     ##### key parent node, val cond. dist. given the single parent. 
        self.parentSetDist = {}
        self.condDist = {}             ########################################## redefine key. This is not time dependent. All relationship to parent should be defined in advance.
        self.leftSeqNum = None   #####################
        self.distT = None      # A function. Input time t, output dictionary which is self.dist. 
    def hashStateLabel(self, digits = 2):        ################################# digits should be the same for all nodes !!!
        self.hashDigits = digits
        if len(self.stateSpace) == 0:
            print('No element in state space yet, please add state space for node ' + self.uniqueName)
            return
        for (i, state) in enumerate(self.stateSpace):
            string = str(i).zfill(digits)
            self.stateStrToNum[state] = string
            self.statesNumToStr[string] = state
    def addStateSpace(self, ssList):    # add statespace list and hash at the same time
        self.stateSpace = ssList
        self.hashStateLabel()

    def addParent(self, newpInput):   # input either list of parent nodes or a single node
        if type(newpInput) is list:
            for newp in newpInput:
                if newp not in self.parentSet:
                    self.parent.append(newp)
                    self.parentSet.add(newp)
                    self.parentNameToNode[newp.name] = newp
                    # print('newp.name in node.py', newp.name)
                if self not in newp.child:
                    newp.child.add(self)
        else:
            newp = newpInput
            if newp not in self.parentSet:
                self.parent.append(newp)
                self.parentSet.add(newp)
                self.parentNameToNode[newp.name] = newp
            if self not in newp.child:
                newp.child.add(self)

    def removeParent(self, delpInput):
        if type(delpInput) is list:
            for delp in delpInput:
                if delp in self.parentSet:
                    self.parent.remove(delp)
                    self.parentSet.remove(delp)
                    self.parentNameToNode.pop(delp.name, None)
                if self in delp.child:
                    delp.child.remove(self)
                else:
                    print('The node to be removed is not in parent set')

        else:
            delp = delpInput
            if delp in self.parentSet:
                self.parent.remove(delp)
                self.parentSet.remove(delp)
                self.parentNameToNode.pop(delp.name, None)
            if self in delp.child:
                delp.child.remove(self)
            else:
                print('The node to be removed is not in parent set')


    def addChild(self, newc):
        self.child.add(newc)
    def removeChild(self, delc):
        if delc in self.child:
            self.child.remove(delc)
        else:
            print('The node to be removed is not in child set')

    def addPrev(self, prevNode):
        self.prev = prevNode
        prevNode.next = self

    def setPfail(self, p):
        if not needSample:
            self.pFail = p
        else:
            print("The node needs to be sampled, can't set value directly")



    def eventTreeRelationship(self, ET):
        self.treeType = 'ET'
        self.ET = ET
        self.relateOrder = ET.top_events
        self.stateSpace = self.ET.simple_es
        self.hashStateLabel()

        # self.relateToP = ET.vectorTElogic['simpleES']
    def solveETrelatinoship(self, truncMethod, normET = None, truncateProb = None, ETcontESdic = None, contSet = None):  # truncMethod 'incAllES' or 'threshold'
        if self.treeType != 'ET':
            print('''The node can't be computed by an event tree, the relationship type is ''' + self.treeType + '.')
        else:
            for p in self.parent:
                teName = p.name
                self.ET.update_TE_fprob(teName, p.dist[p.stateStrToNum['f']])
            self.ET.solve()

            if truncMethod == 'threshold':  # if method is this, need to specify truncateProb as well
                if truncateProb != None:
                    self.ET.truncate(truncateProb, normET)
                    self.stateSpace = self.ET.truncES
                    self.leftSeqNum = self.ET.leftSeqNum

                # print('node, self.stateSpace in solveETrelatinoship', self.stateSpace)

                # self.hashStateLabel()
                # print('node self.stateStrToNum', self.stateStrToNum)
                dist = self.ET.truncatedESprob
                for key in dist:
                    self.dist[self.stateStrToNum[key]] = dist[key]

            elif truncMethod == 'incAllES': 
                self.ET.truncateIncAllES(normET, ETcontESdic , contSet)
                self.stateSpace = self.ET.truncES
                self.leftSeqNum = self.ET.leftSeqNum

                dist = self.ET.truncatedESprob
                for key in dist:
                    self.dist[self.stateStrToNum[key]] = dist[key]


            elif truncMethod == None:
                # self.stateSpace = self.ET.simple_es
                # self.hashStateLabel()
                dist = self.ET.endstate_prob
                self.dist = {self.stateStrToNum[x] : y for x, y in dist.items()}

        
        
    def addCondRelation(self, file = None, singleParentRelate = None, func = None):    ################## if add single parent relate, singleParentRelate = pNode
        if self.typ == 'sfn':  # self only has succ / fail states, simpler, update self.condDist, self.relateOrder
            with open(file) as f:
                seperator = None  # location to seperate parent and value
                rows = list(f)
                title = rows[0]
                titleList = title[:-1].split()

                if titleList[-1] == 'Median':
                    '''
                    This indicates dep. value is given by mean and median,
                    has uncertainty, need sample for value
                    '''
                    self.needSample = True
                    self.relateOrder = titleList[:-2]
                    seperator = -2

                else:
                    '''
                    This indicates dep. value is a const, no uncertainty.
                    '''
                    self.relateOrder = titleList[:-1]
                    seperator = -1

                for row in rows[1:]:
                    rowList = row[:-1].split()
                    if len(rowList) == 0:
                        break
                    if seperator == -2:
                        parentStats = rowList[:-2]  # This doesn't include parent names, only states like s, f, n..
                        median = float(rowList[-1])
                        mean = float(rowList[-2])
                        if len(parentStats) == 0:
                            self.mean = mean
                            self.median = median
                            self.needSample = True
                            return

                        for (i, stat) in enumerate(parentStats):
                            parentStats[i] = self.stateStrToNum[stat]

                        key = ''.join(parentStats)   ### key is num string combination here, not chars string. 

                        if mean == median:
                            self.condDist[key] = {}
                            self.condDist[key]['val'] = mean
                            self.needSample = False
                        else:
                            self.condDist[key] = {}
                            self.condDist[key]['mean'] = mean
                            self.condDist[key]['median'] = median

                    else:
                        parentStats = rowList[:-1]
                        val = float(rowList[-1])
                        if len(parentStats) == 0:
                            self.pFail = val
                            return
                        key = ''.join(parentStats)
                        self.condDist[key] = {}
                        self.condDist[key]['val'] = val
        elif self.typ == 'multi' and singleParent != None: #### only update self.parentSingleDist, don't update self.condDist or self.relateOrder. They are updated in compute graph
            self.addSinglePcondRel(file, singleParent)

        else:  # typ is not 'sfn' or 'multi'
            if func != None:
                self.condDist = func
            else:
                print('condDist function for node ' + self.uniqueName + ' is not defined, please input it.')



    def addSinglePcondRel(self, file, pNode):   ######### updates parentSingleDist
        df = pd.read_csv(file, delim_whitespace = True)
        self.stateSpace = list(df)
        self.hashStateLabel()
        dic = df.to_dict(orient = 'index')
        # print('node, dic', dic)
        # print(dic)
        singlePdic = {}
        for key in dic:
            # print('node key', key)
            # print('node pNode.stateStrToNum', pNode.stateStrToNum)
            # print('node pNode.uniqueName', pNode.uniqueName)

            if key in pNode.stateStrToNum:

                selfNumStr = pNode.stateStrToNum[key]
                # print('node selfNumStr is ', selfNumStr)
                singlePdic[selfNumStr] = {self.stateStrToNum[x]:y for x,y in dic[key].items() if x in self.stateStrToNum and y!= 0}
                # for x, y in dic[key].items():
                #     if 
                # {pNode.stateStrToNum[x]:y for x,y in dic[key].items() if y!= 0}
        self.parentSingleDist[pNode] = singlePdic


    def updateCondDist(self):   ####### Only after all single parent cond dists are defined. Update node.relateOrder and node.condDist here.
        if len(self.parent) != len(self.parentSingleDist):
            print('Single parent dist. is not complete for node ' + node.uniqueName)
            return
        else:
            for p in self.parent:
                self.relateOrder.append(p.name)
            if len(self.relateOrder) == 1:
                parentNode = self.parent[0]
                self.condDist = self.parentSingleDist[parentNode]
            else:
                print('More than one parent is not implemented yet.')
                return

    def updateTimeDist(self, func):   # update time dependent dist function. func is a func, this will make self.distT = func. 
        self.distT = func


    def getParentSetDist(self):
        "To be implemented!"
        return
























    # def faultTreeRelationship(self, FT):
    #     '''To be implemented'''

    # def solveFTrelationship(self, FT):
    #     '''To be implemented'''

# '''
# Node by state. Each node in this type is represented by single state, not by whole distribution. Good thing is it only takes short time to get from parent state to child state (basically only involves query in dictionary). Bad thing is we need to MonteCarlo to get the distribution. But this is hard because PRA involves many rare events that are very hard to simulate.
# '''

# class nodeByState(object)