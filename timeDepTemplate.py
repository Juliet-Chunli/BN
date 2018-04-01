# time depenent Peach Bottom power plant analysis. BWR, only consider RCIC as cooling method. 

# node.condDist should be defined for specific cases if it is not generic (can't be read from files). The point of node.condDist is that given parent states, define node distribution. condDist could be function and is not time dependent. Even though the value of node may depend on time, conditional distribution which defines the relationship to parent is not time dependent and should be defined in advance. node.condDist could be a function which can be updated by node.addCondRelation. The function should be defined as input node.parentSetDist (a dictionary as defined below), output dictionary = self.dist

# node.parentSetDist has key parent name, value dist of the parent. 

# node.updateTimeDist(func) will update node.distT to be func. So the definition of func is input t, output self.dist (a dictionary).

# nodes need to be added to G by G.addNods(nodes), list or single node
# prev nodes need to be added to G by G.addPrevNode(prevNodes), list or single node

# all node with prev needs to be initiated by G.initT0(initDic)

from node import *
from scipy.stats import lognorm
import numpy as np
import matplotlib.pyplot as plt
from computeGraph import diGraph
import matplotlib.ticker as mtick


############## U1 nodes ##############

batteryPrev_u1   = nodeByDist('batteryPrev_u1', 'batteryPrev', 1, typ = 'tIterPrev')
rcicManualPrev_u1 = nodeByDist('rcicManualPrev_u1', 'rcicManualPrev', 1, typ = 'tIterPrev')
dcPrev_u1        = nodeByDist('dcPrev_u1', 'dcPrev', 1, typ = 'tIterPrev')
rcicFuncPrev_u1  = nodeByDist('rcicFuncPrev_u1', 'rcicFuncPrev', 1, typ = 'tIterPrev' )
rcicFloodPrev_u1 = nodeByDist('rcicFloodPrev_u1', 'rcicFloodPrev', 1, typ = 'tIterPrev' )
coreStatePrev_u1 = nodeByDist('coreStatePrev_u1', 'coreStatePrev', 1, typ = 'tIterPrev')
portableDCPrev_u1    = nodeByDist('portableDCPrev_u1', 'portableDCPrev', 1, typ = 'tIterPrev')

battery_u1       = nodeByDist('battery_u1', 'battery', 1)
rcicManual_u1     = nodeByDist('rcicManual_u1', 'rcicManual', 1)
dc_u1            = nodeByDist('dc_u1', 'dc', 1)
rcicFunc_u1      = nodeByDist('rcicFunc_u1', 'rcicFunc', 1)
rcicFlood_u1     = nodeByDist('rcicFlood_u1', 'rcicFlood', 1)
coreState_u1     = nodeByDist('coreState_u1', 'coreState', 1)
portableDC_u1    = nodeByDist('portableDC_u1', 'portableDC', 1)

rcicControl_u1   = nodeByDist('rcicControT_u1', 'rcicControl', 1)
coreCool_u1      = nodeByDist('coreCool_u1', 'coreCool', 1)
cont_u1          = nodeByDist('cont_u1', 'cont', 1)

# U1 edges
battery_u1.addParent(batteryPrev_u1)
rcicManual_u1.addParent(rcicManualPrev_u1)
dc_u1.addParent(dcPrev_u1)
rcicFunc_u1.addParent(rcicFuncPrev_u1)
rcicFlood_u1.addParent(rcicFloodPrev_u1)
coreState_u1.addParent(coreStatePrev_u1)
portableDC_u1.addParent(portableDCPrev_u1)

dc_u1.addParent([portableDC_u1, battery_u1])
rcicControl_u1.addParent(dc_u1)
rcicManual_u1.addParent(battery_u1)
rcicFlood_u1.addParent([rcicControl_u1, rcicFuncPrev_u1])
rcicFunc_u1.addParent([rcicManual_u1, rcicFlood_u1])
coreCool_u1.addParent(rcicFunc_u1)
coreState_u1.addParent(coreCool_u1)
cont_u1.addParent(coreState_u1)

# U1  relationship
battery_u1.addPrev(batteryPrev_u1)
rcicManual_u1.addPrev(rcicManualPrev_u1)
dc_u1.addPrev(dcPrev_u1)
rcicFunc_u1.addPrev(rcicFuncPrev_u1)
rcicFlood_u1.addPrev(rcicFloodPrev_u1)
coreState_u1.addPrev(coreStatePrev_u1)
portableDC_u1.addPrev(portableDCPrev_u1)

# U1 all nodes
nodesU1 = [batteryPrev_u1, rcicManualPrev_u1, portableDCPrev_u1, dcPrev_u1, rcicFuncPrev_u1, rcicFloodPrev_u1, coreStatePrev_u1, battery_u1, rcicManual_u1, dc_u1, rcicFunc_u1, rcicFlood_u1, coreState_u1, portableDC_u1, rcicControl_u1, coreCool_u1, cont_u1]

# U1 prev nodes
prevNodesU1 = [batteryPrev_u1, rcicManualPrev_u1, portableDCPrev_u1, dcPrev_u1, rcicFuncPrev_u1, rcicFloodPrev_u1, coreStatePrev_u1]


############## U2 nodes ##############

batteryPrev_u2   = nodeByDist('batteryPrev_u2', 'batteryPrev', 2, typ = 'tIterPrev')
rcicManualPrev_u2 = nodeByDist('rcicManualPrev_u2', 'rcicManualPrev', 2, typ = 'tIterPrev')
dcPrev_u2        = nodeByDist('dcPrev_u2', 'dcPrev', 2, typ = 'tIterPrev')
rcicFuncPrev_u2  = nodeByDist('rcicFuncPrev_u2', 'rcicFuncPrev', 2, typ = 'tIterPrev' )
rcicFloodPrev_u2 = nodeByDist('rcicFloodPrev_u2', 'rcicFloodPrev', 2, typ = 'tIterPrev' )
coreStatePrev_u2 = nodeByDist('coreStatePrev_u2', 'coreStatePrev', 2, typ = 'tIterPrev')
portableDCPrev_u2    = nodeByDist('portableDCPrev_u2', 'portableDCPrev', 2, typ = 'tIterPrev')

battery_u2       = nodeByDist('battery_u2', 'battery', 2)
rcicManual_u2     = nodeByDist('rcicManual_u2', 'rcicManual', 2)
dc_u2            = nodeByDist('dc_u2', 'dc', 2)
rcicFunc_u2      = nodeByDist('rcicFunc_u2', 'rcicFunc', 2)
rcicFlood_u2     = nodeByDist('rcicFlood_u2', 'rcicFlood', 2)
coreState_u2     = nodeByDist('coreState_u2', 'coreState', 2)
portableDC_u2    = nodeByDist('portableDC_u2', 'portableDC', 2)

rcicControl_u2   = nodeByDist('rcicControT_u2', 'rcicControl', 2)
coreCool_u2      = nodeByDist('coreCool_u2', 'coreCool', 2)
cont_u2          = nodeByDist('cont_u2', 'cont', 2)

# U2 edges
battery_u2.addParent(batteryPrev_u2)
rcicManual_u2.addParent(rcicManualPrev_u2)
dc_u2.addParent(dcPrev_u2)
rcicFunc_u2.addParent(rcicFuncPrev_u2)
rcicFlood_u2.addParent(rcicFloodPrev_u2)
coreState_u2.addParent(coreStatePrev_u2)
portableDC_u2.addParent(portableDCPrev_u2)

dc_u2.addParent([portableDC_u2, battery_u2])
rcicControl_u2.addParent(dc_u2)
rcicManual_u2.addParent(battery_u2)
rcicFlood_u2.addParent([rcicControl_u2, rcicFuncPrev_u2])
rcicFunc_u2.addParent([rcicManual_u2, rcicFlood_u2])
coreCool_u2.addParent(rcicFunc_u2)
coreState_u2.addParent(coreCool_u2)
cont_u2.addParent(coreState_u2)

# U2  relationship
battery_u2.addPrev(batteryPrev_u2)
rcicManual_u2.addPrev(rcicManualPrev_u2)
dc_u2.addPrev(dcPrev_u2)
rcicFunc_u2.addPrev(rcicFuncPrev_u2)
rcicFlood_u2.addPrev(rcicFloodPrev_u2)
coreState_u2.addPrev(coreStatePrev_u2)
portableDC_u2.addPrev(portableDCPrev_u2)

# U2 all nodes
nodesU2 = [batteryPrev_u2, rcicManualPrev_u2, portableDCPrev_u2, dcPrev_u2, rcicFuncPrev_u2, rcicFloodPrev_u2, coreStatePrev_u2, battery_u2, rcicManual_u2, dc_u2, rcicFunc_u2, rcicFlood_u2, coreState_u2, portableDC_u2, rcicControl_u2, coreCool_u2, cont_u2]

# U2 prev nodes
prevNodesU2 = [batteryPrev_u2, rcicManualPrev_u2, portableDCPrev_u2, dcPrev_u2, rcicFuncPrev_u2, rcicFloodPrev_u2, coreStatePrev_u2]

# Interaction Edge
portableDC_u2.addParent(cont_u1)


############## Conditional Functions, use 1, 0 represent succ, fail ##############

# State battery availability at time t, battery state tuple (stateIndicator, time). 1 in stateIndicator means succ, 0 fail. time is the time into the state
def battery(t, deltaT, pdist, life = 4):  # t has to start from 0 in compute graph
	state, time = pdist['batteryPrev']
	newState = state
	newTime = time + deltaT
	if newState == 1 and newTime >= life:
		return (0, 0)
	else:
		return (newState, newTime)

# riciManual start states, tuple (state, needStart, time)
def rcicManual(t, deltaT, pdist, timeToStart = 1): 
	state, needStart, time = pdist['rcicManualPrev']
	if needStart == 0: # if don't need start, not start success, time increases
		return (0, 0, time + deltaT)
	else: # need start, check time
		if time + deltaT >= timeToStart: # reaches necessary working time, succ start, no need to start anymore
			return (1, 0, time + deltaT)
		else: # haven't reached the time, time accumulates, still need to start
			return (0, 1, time + deltaT)

# DC recover at time t. Triplet (state, time, recoveryStart). Need to define multiple recovery effort later. 
def portableDC(t, deltaT, pdist, mu = 0.3, sigma = 1.064): # unit hour, may add recStartDelay to include multiple restore efforts
	stateP, timeP, recStartTimeP =  pdist['portableDCPrev']
	if 'cont' in pdist and pdist['cont'] == 0: # if contS fails, damages dc recovery
		if stateP == 0: # didn't recover
			return (0, timeP + deltaT, t) # (fails, fail time accumulates, recovery start time set to now)
		else: # was ok
			return (0, 0, t) # (fails, fail time = 0, recovery start time set to now)

	else:
		if stateP == 1:  # succ at previous time, keep success
			return (stateP, timeP + deltaT, recStartTimeP) # if it's still ok, no need to recovery. Restore starting time is current. 
		elif t + deltaT < recStartTimeP: # fail at previous time, and t + deltaT < previous start time, recover fail, time passes, startT no change
			return (0, timeP + deltaT, recStartTimeP)
		else: # fail at previous time, sample to see whether recovery at this time step
			psucc = lognorm.cdf(t + deltaT - recStartTimeP, s = sigma, scale = np.exp(mu)) - lognorm.cdf(t - recStartTimeP, s = sigma, scale = np.exp(mu))  # correct time based on starting time of current resotre work
			state = np.random.choice(2, p = [1 - psucc, psucc])
			if state == 1: # succ by sample
				return (state, 0, recStartTimeP)
			else: # failure time accumulated
				return (state, timeP + deltaT, recStartTimeP)

# dc status at time t, tuple (state, time)
def dc(t, deltaT, pdist):
	stateP, timeP = pdist['dcPrev']
	batteryS = pdist['battery'][0]
	recoverS = pdist['portableDC'][0]
	if stateP == 1: # was succ
		if batteryS == 1 or recoverS == 1: # keep succ
			return (1, timeP + deltaT)
		else:
			return (0, 0)
	else: # was fail
		if batteryS == 1 or recoverS == 1: # restore succ
			return (1, 0)
		else:
			return (0, timeP + deltaT)

# RCIC control status, exactly the same as dc, tuple (state, time)
def rcicControl(t, deltaT, pdist):
	return (pdist['dc'])

# RCIC flood, state 1 or 0 for flood and not
def rcicFlood(t, deltaT, pdist, noCtrlWorkT = 1.5):
	sP = pdist['rcicFloodPrev']
	if sP == 1: # if ever flooded, always flooded
		return 1
	else: # never flooded
		rcicS, rcicT = pdist['rcicFuncPrev']
		rcicCtrlS, rcicCtrlT = pdist['rcicControl']
		if rcicS == 0: # rcic not working no flood
			return 0
		else: # rcic working
			if rcicCtrlS == 1: # rcic control ok no flood
				return 0
			else: # rcic control fail
				if min(rcicT, rcicCtrlT) >= noCtrlWorkT: # rcic work without control more than could, flood
					return 1
				else: # rcic work without control less than could, no flood
					return 0


# RCIC work status, (state, time)
def rcicFunc(t, deltaT, pdist):
	workS, workT = pdist['rcicFuncPrev']
	flood = pdist['rcicFlood']
	start, needStart, startTime = pdist['rcicManual']
	if flood == 1: # if flooded
		if workS == 1: # worked
			return (0, 0) # worked and flooded, (notwork, time = 0)
		else: # not work, time accumulates
			return (0, workT + deltaT)
	else: # never flooded
		if workS == 0 and start == 1: # didn't work but succ start
			return (1, 0)
		else: # if work succ, start fail, no change to workS; if work succ, no need to start, no change to workS
			return (workS, workT + deltaT) # nothing will change workS, (workS, t + dt)




# Core Cool State, (state, time)
def coreCool(t, deltaT, pdist):
	s = pdist['rcicFunc']  # same state as RCIC work
	return s


# Core state, (state, accumulateNoCoolTime), assume can stand if accumulated no cooling 10hr
def coreState(t, deltaT, pdist, noCoolThres = 10):
	state, noCoolTime = pdist['coreStatePrev']
	coolS, coolT = pdist['coreCool']
	if state == 1: # was ok
		if coolS == 1: # with cool
			return (1, noCoolTime) # was ok, with cool (ok, noCoolTime stays)
		elif noCoolTime + deltaT < noCoolThres: # accumulate no cool time < thres
			return (1, noCoolTime + deltaT) # was ok, no cool (ok, noCoolTime accumulates)
		else: # was ok, but accumulate noCoolTime > thres
			return (0, 0)
	else: # already fails
		if coolS == 1: # already fails with cooling now
			return (0, noCoolTime) # (fails, noCoolTime stays) 
		else: # already fails without cooling
			return (0, noCoolTime + deltaT) # (fails, noCoolTime accumlates) 


# Containment state, assume fail right after core fails for now
def cont(t, deltaT, pdist, coreFailThres = 0.5):
	coreState, noCoolTime = pdist['coreState']
	if coreState == 1:
		return 1
	else:
		return 0



def initGraph(prevNodes, batterySuccP = 0.8, timeToStartRCIC = 1): # make sure node order in list follows topological order

	def initBattery(batterySuccP):
		state = np.random.choice(2, p = [1 - batterySuccP, batterySuccP])
		return (state, 0)
	def initRCICmanual(batteryS):
		if batteryS == 0:
			return (0, 1, 0)  # (start not succ, needStart, time = 0)
		else: # with init battery, (start not succ, needStart, time = 0)
			return (0, 0, 0)
	def initportableDC(rcicNeedStart, timeToStartRCIC):
		if rcicNeedStart == 1: # rcic need start, start recover after starting rcic
			return (0, 0, timeToStartRCIC)  # (failRecover, inThisStateTime = 0, startTime = timeToStartRCIC)
		else:
			return (0, 0, 0) # (failRecover, inThisStateTime = 0, startTime = 0)

	batteryS = None
	rcicNeedStart = None
	d = {}
	for node in prevNodes:
		if node.name == 'batteryPrev':
			batteryS, batteryT = initBattery(batterySuccP)
			print(batteryS)
			if batteryS == 0: # no init battery, rcic need start
				rcicNeedStart = 1
			else: # with battery, rcic no need start
				rcicNeedStart = 0
			d[node] = (batteryS, batteryT)
		elif node.name == 'rcicManualPrev':
			d[node] = initRCICmanual(batteryS)
		elif node.name == 'portableDCPrev':
			d[node] = initportableDC(rcicNeedStart, timeToStartRCIC)
		elif node.name == 'dcPrev':
			d[node] = (batteryS, 0) # dc init same as battery
		elif node.name == 'rcicFuncPrev':
			d[node] = (batteryS, 0) # rcic work if battery works
		elif node.name == 'rcicFloodPrev': # init never flooded
			d[node] = 0
		elif node.name == 'coreStatePrev': # core init OK
			d[node] = (1, 0)

	return d


############## Add cond Dist for all non-prev nodes ##############


battery_u1.condDist = battery
rcicManual_u1.condDist = rcicManual
dc_u1.condDist = dc
rcicFunc_u1.condDist = rcicFunc
rcicFlood_u1.condDist = rcicFlood
coreState_u1.condDist = coreState
portableDC_u1.condDist = portableDC
rcicControl_u1.condDist = rcicControl
coreCool_u1.condDist = coreCool
cont_u1.condDist = cont


battery_u2.condDist = battery
rcicManual_u2.condDist = rcicManual
dc_u2.condDist = dc
rcicFunc_u2.condDist = rcicFunc
rcicFlood_u2.condDist = rcicFlood
coreState_u2.condDist = coreState
portableDC_u2.condDist = portableDC
rcicControl_u2.condDist = rcicControl
coreCool_u2.condDist = coreCool
cont_u2.condDist = cont


print('Unit 1 init battery state')
initDicU1 = initGraph(prevNodesU1)
print('Unit 2 init battery state')
initDicU2 = initGraph(prevNodesU2)
initDic = {**initDicU1, **initDicU2}




# if initDic['batteryPrev_u1'][0] == 1:
# 	print('u1 plant battery ok')
# else:
# 	print('u1 plant battery fails')


# if initDic['batteryPrev_u2'][0] == 1:
# 	print('u2 plant battery ok')
# else:
# 	print('u2 plant battery fails')



G = diGraph()
G.addNode(nodesU1)
G.addPrevNode(prevNodesU1)
G.addNode(nodesU2)
G.addPrevNode(prevNodesU2)
G.initT0(initDic)




nodeOfInterest = [portableDC_u1, dc_u1, coreState_u1, cont_u1, portableDC_u2, dc_u2, coreState_u2, cont_u2, rcicFunc_u1, rcicFunc_u2, rcicFlood_u1, rcicFlood_u2, rcicManual_u1, rcicManual_u2]

startTime = 0
endTime = 40
timeStep = 0.5



G.computeInTime(startTime, endTime, timeStep, nodeOfInterest)

# print('G.timeDepResult')

# print(G.timeDepResult)


caseNum = 6
T = G.timeSeq
dcRec_u1 = [tup[0] for tup in G.timeDepResult['portableDC_u1']]
dcStat_u1 = [tup[0] for tup in G.timeDepResult['dc_u1']]
coreS_u1 = [tup[0] for tup in G.timeDepResult['coreState_u1']]
contS_u1 = [s for s in G.timeDepResult['cont_u1']]

dcRec_u2 = [tup[0] for tup in G.timeDepResult['portableDC_u2']]
dcStat_u2 = [tup[0] for tup in G.timeDepResult['dc_u2']]
coreS_u2 = [tup[0] for tup in G.timeDepResult['coreState_u2']]
contS_u2 = [s for s in G.timeDepResult['cont_u2']]

Y = np.array([0, 1])
Ytick = ['f', 's']

plt.suptitle('Unit 1 situation Case ' + str(caseNum))
plt.subplot(2, 2, 1)
plt.plot(T, dcRec_u1)
plt.xlabel('t / hr')
plt.yticks(Y, Ytick)
plt.title('DC Restore Status')


plt.subplot(2, 2, 2)
plt.plot(T, dcStat_u1)
plt.xlabel('t / hr')
plt.yticks(Y, Ytick)
plt.title('DC Status')

plt.subplot(2, 2, 3)
plt.plot(T, coreS_u1)
plt.xlabel('t / hr')
plt.yticks(Y, Ytick)
plt.title('Core Status')


plt.subplot(2, 2, 4)
plt.plot(T, contS_u1)
plt.xlabel('t / hr')
plt.yticks(Y, Ytick)
plt.title('Containment Status')


plt.tight_layout()
plt.subplots_adjust(top=0.85)
plt.savefig('timeDepPeachBottomU1Case' + str(caseNum) + '.png', dpi = 300)

plt.show()


plt.suptitle('Unit 2 situation Case ' + str(caseNum))
plt.subplot(2, 2, 1)
plt.plot(T, dcRec_u2)
plt.xlabel('t / hr')
plt.yticks(Y, Ytick)
plt.title('DC Restore Status')


plt.subplot(2, 2, 2)
plt.plot(T, dcStat_u2)
plt.xlabel('t / hr')
plt.yticks(Y, Ytick)
plt.title('DC Status')

plt.subplot(2, 2, 3)
plt.plot(T, coreS_u2)
plt.xlabel('t / hr')
plt.yticks(Y, Ytick)
plt.title('Core Status')


plt.subplot(2, 2, 4)
plt.plot(T, contS_u2)
plt.xlabel('t / hr')
plt.yticks(Y, Ytick)
plt.title('Containment Status')


plt.tight_layout()
plt.subplots_adjust(top=0.85)
plt.savefig('timeDepPeachBottomU2Case' + str(caseNum) + '.png', dpi = 300)

plt.show()



result = G.timeDepResult
# print('rcicFunc_u1', result['rcicFunc_u1'])
# print('rcicFlood_u1', result['rcicFlood_u1'])
# print('rcicManual_u1', result['rcicManual_u1'])
# print('rcicFunc_u2', result['rcicFunc_u2'])
# print('rcicFlood_u2', result['rcicFlood_u2'])
# print('rcicManual_u2', result['rcicManual_u2'])
