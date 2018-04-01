#!/usr/bin/python3.6
""" Build event tree dictionary from input file. Input variable is the name of inputfile. eg: 'LOSS_OF_MAIN_FW.txt'
    return a list of list. # 0 element is Init. event, it is not a list, just string.
    #1 element TE, #2 and after tree logic and ES.
    dic['IE']: initiating event (string)
    dic['top_events']: list of top events
    dic['end_states']: list of end states, may include duplicates
    dic['tree_logic']: dictionary of tree logic, access a logic by self.tree_logic['05 SLFC']['top event']. 's' denotes the specific top event successes in path to end state; 'f' is fail; None is not considered. Whenever meets None, just times 1 in computing probability. dic['vectorTElogic']['fullES']: dic. of relationship of vector TE to full endstate(include seqNum). Key vectorTE = 'fsnssnf' is the combination of TE status in order of dic['top_events'] dic['vectorTElogic']['simpleES']: dic. of relationship of vector TE to simple endstate(no seqNum). Key vectorTE = 'fsnssnf' is the combination of TE status in order of dic['top_events']

"""
import numpy as np

def build_ET_dic(inputFileName):
    file = open(inputFileName, 'r')
    matrix = []

    for line in file:
        line = line.split()
        matrix.append(line)
    file.close()
    dic = {}
#     print(matrix[2:])
#     logic_matrix_withES = matrix[:][2:]
#     logic_matrix_withES = logic_matrix_withES[:]
#     print(logic_matrix_withES)
    dic['IE'] = matrix[0]
    dic['top_events'] = matrix[1]
    dic['end_states'] = []
    dic['tree_logic'] = {}
    dic['vectorTElogic'] = {}
    dic['vectorTElogic']['simpleES'] = {}
    dic['vectorTElogic']['fullES'] = {}
    for l in matrix[2:]:
        seqNum = l.pop()
        es = l.pop()
        dic['end_states'].append(' '.join([seqNum, es]))
    logic_matrix = matrix[2:]
    for i, es in enumerate(dic['end_states']):
        vectorTE = ''.join(logic_matrix[i])

        dic['vectorTElogic']['fullES'][vectorTE] = es
        dic['vectorTElogic']['simpleES'][vectorTE] = es.split(' ')[1]
        dic['tree_logic'][es] = {}
        for j, te in enumerate(dic['top_events']):
            dic['tree_logic'][es][te]=logic_matrix[i][j]

#     print(len(dic['vectorTElogic']['simpleES']))
#     print(len(dic['vectorTElogic']['fullES']))
#     print(len(dic['tree_logic']))
    return dic


'''
    This represents an event tree.
    self.IE: string indicating initiating event;
    self.top_events: list of top events;
    self.end_states: list of end states, eg: 01 LTFC;
    self.simple_es: end state list without duplicate
    self.tree_logic: dictionary of tree logic, access a logic by self.tree_logic[(i, 'end state')]['top event']. Note i in the tuple starts from 1.
    self.vectorTElogic: dictionary of corres. relat. of TE seq and endState. self.vectorTElogic['simpleES'] and
                     self.vectorTElogic['fullES'] are dic.s of simple ES and full ES, respectively. Key for those are strings
                     'fsnssnf' which is the combination of TE status in order of self.top_events

    The following atributes will only be nonempty after performing functions.
    self.te_fprob: dictionary of top event failure prob. Will only be nonempty after func update_TE_fprob or batchUpdateTEfprob.
    self.sequence_prob: dictionary of sequence probability. Only nonempty after func solve
    self.truncatedSeqProb: dic. of truncated sequence probability (unnormalized to make sure dist. shape doesn't change), non-empty after truncate. 
    self.endstate_prob: dictionary of simple ES probability. Only nonempty after func solve
    self.truncatedESprob: dictionary of trucated simple ES probability, only non-empty after self.truncate()
    self.TEdistDic: dictionary of TE distribution, assume lognormal, self.TEdistDic['TE Name'] is a list of [mu, std]. mu and std are mean and std for underling normal dist.
    self.MCresult: a dic of MC result with sequence name (num + ES name) as key. Include two parts: self.MCresult['bySeq'] and self.MCresult['bySample']. self.MCresult['bySeq'] is a dic that with seqName as key. Each element is a list of MC results for that specific sequence. MCresult['bySample] is a list of list. Each element of the list is the prob. of end sequencs corres. to one specifc sample. length of the element is the number of end states (seq), lenght of MCresult['bySample'] is number samples.

    The following are the functions.
    update_TE_fprob(topEventName, failprob)
        Input top event name as str and corres. failprob as float, will change self.te_fprob.
    batchUpdateTEfprob(failprobLlist):
        Input list of failprob in same order as TE, update all fprob at once.
    solve(self)
        Solve the ET. Will only work after updating all TE fail prob.
    truncate(self, riskRatio)
        Truncate the result based on riskRatio
    batchInputTEdist(distFileName)
        Input TE distribution file name, assume lognormal, assign value for self.TEdistDic to prepare for Monte Carlo.
    MonteCarloSeq(self, num)
        Monte Carlo based on pdf of ET top events, num is number of samples. Only applicable after batchInputTEdist. Assign value for self.MCresult['bySeq'] and self.MCresult['bySample'].
    '''

class Event_tree(object):
    def __init__(self, inputFileName): # end_state is list of end states
        dic = build_ET_dic(inputFileName)
        self.IE = dic['IE']
        self.top_events = dic['top_events']
        self.end_states = dic['end_states']
        endStateNoNum = []
        for es in self.end_states:
            num, esName = es.split()
            endStateNoNum.append(esName)



        self.simple_es = list(set(endStateNoNum))  # change to set first to remove redundencies
        self.tree_logic = dic['tree_logic']
        self.vectorTElogic = dic['vectorTElogic']
        self.te_fprob = {}
        self.sequence_prob = {}  # prob of each sequence
        self.truncatedSeqProb = {} 
        self.endstate_prob = {}
        self.truncatedESprob = {}  
        self.TEdistDic = {}
        self.MCresult = {}
        self.truncES = None #################### write !!!!!!
        self.leftSeqNum = None   ################### 


    def update_TE_fprob(self, topEventName, failprob):   # Input top event name as str and failprob as float
        self.te_fprob[topEventName] = failprob

    def batchUpdateTEfprob(self, failprobList):
        for i, fprob in enumerate(failprobList):
            self.update_TE_fprob(self.top_events[i],  fprob)

    def batchInputTEdist(self, distFileName):
        file = open(distFileName, 'r')
        next(file)
        for line in file:
            name, mean, median = line.split()
            mean = float(mean)
            median = float(median)
            mu = np.log(median)
            std = (2*(np.log(mean) - mu))**0.5
            self.TEdistDic[name] = [mu, std]
        file.close()
    def solve(self):
        for i, es in enumerate(self.tree_logic):
            self.sequence_prob[es] = 1
            for te in self.tree_logic[es]:
                if te in self.te_fprob and self.te_fprob[te] != None:
                    if self.tree_logic[es][te] == 's':
                        self.sequence_prob[es] *=  (1 - self.te_fprob[te])
                    elif self.tree_logic[es][te] == 'f':
                        self.sequence_prob[es] *=  self.te_fprob[te]
                elif te in self.te_fprob and self.te_fprob[te] == None:
                    print('Top event ' + te + 'has prob = None.')
                else:
                    print("Top event ' + te + 'doesn't have a prob. yet")
        for i, aES in enumerate(self.sequence_prob):
            num, es = aES.split()
            if es not in self.endstate_prob:
                self.endstate_prob[es] = self.sequence_prob[aES]
            else:
                self.endstate_prob[es] += self.sequence_prob[aES]

    def truncate(self, riskRatio, normalize):
        succP = 0
        normConst = 0
        for seq in self.sequence_prob:
            ind, es = seq.split(' ')
            if es == 'SUCCESS':
                succP += self.sequence_prob[seq]
                self.truncatedSeqProb[seq] = self.sequence_prob[seq]
                normConst += self.sequence_prob[seq]
        riskP = 1 - succP
        considered = 0
        riskThres = riskP * riskRatio

        for seq in sorted(self.sequence_prob, key = self.sequence_prob.get, reverse = True):
            ind, es = seq.split(' ')
            if es != 'SUCCESS':
                if considered + self.sequence_prob[seq] <= riskThres:
                    normConst += self.sequence_prob[seq]
                    considered += self.sequence_prob[seq]
                    self.truncatedSeqProb[seq] = self.sequence_prob[seq]
                else:
                    break
        self.leftSeqNum = len(self.truncatedSeqProb)
        if normalize:
            for seq in self.truncatedSeqProb:
                self.truncatedSeqProb[seq] /= normConst

        # return len(self.truncatedSeqProb)


        for seq in self.truncatedSeqProb:
            prob = self.truncatedSeqProb[seq]
            ind, es = seq.split(' ')
            if es not in self.truncatedESprob:
                self.truncatedESprob[es] = prob
            else:
                self.truncatedESprob[es] += prob
        # print('event_tree truncatedESprob len', len(self.truncatedESprob))
        self.truncES = list(self.truncatedESprob.keys())
        # print('event_tree, self.truncES', self.truncES)


    def truncateIncAllES(self, normalize, ETcontESdic, contSet):
        contES = contSet.copy() 
        succP = 0
        normConst = 0


        for seq in self.sequence_prob:
            ind, es = seq.split(' ')
            if es == 'SUCCESS':
                succP += self.sequence_prob[seq]
                self.truncatedSeqProb[seq] = self.sequence_prob[seq]
                normConst += self.sequence_prob[seq]

        riskP = 1 - succP
        considered = 0

        # for seq in sorted(self.sequence_prob, key = self.sequence_prob.get, reverse = True):

        #     if contES == set():
        #         break
        #     ind, es = seq.split(' ')
        #     if es != 'SUCCESS':
        #         # print('es in event_tree', es)
        #         if es in ETcontESdic:
        #             contES = contES.difference(ETcontESdic[es])
        #         normConst += self.sequence_prob[seq]
        #         considered += self.sequence_prob[seq]
        #         self.truncatedSeqProb[seq] = self.sequence_prob[seq]


        # Truncate but keep all ET end states

        ETes = set(self.simple_es)
        ETes.remove('SUCCESS')

        for seq in sorted(self.sequence_prob, key = self.sequence_prob.get, reverse = True):
            # print(ETes)

            if ETes == set():
                break
            ind, es = seq.split(' ')
            if es != 'SUCCESS':
                # print('es in event_tree', es)
                if es in ETes:
                    ETes.remove(es)
                normConst += self.sequence_prob[seq]
                considered += self.sequence_prob[seq]
                self.truncatedSeqProb[seq] = self.sequence_prob[seq]






        self.leftSeqNum = len(self.truncatedSeqProb)
        if normalize:
            for seq in self.truncatedSeqProb:
                self.truncatedSeqProb[seq] /= normConst

        # return len(self.truncatedSeqProb)


        for seq in self.truncatedSeqProb:
            prob = self.truncatedSeqProb[seq]
            ind, es = seq.split(' ')
            if es not in self.truncatedESprob:
                self.truncatedESprob[es] = prob
            else:
                self.truncatedESprob[es] += prob
        # print('event_tree truncatedESprob len', len(self.truncatedESprob))
        self.truncES = list(self.truncatedESprob.keys())
        # print('event_tree, self.truncES', self.truncES)
        




    def MonteCarloSeq(self, num):
        self.MCresult['bySeq'] = {}
        self.MCresult['bySample'] = []
        mu = []
        std = []
        for te in self.top_events:
            if te not in self.TEdistDic:
                return ('Please input pdf for top event ' + te )
            else:
                mu.append(self.TEdistDic[te][0])
                std.append(self.TEdistDic[te][1])
        mu = np.array(mu)
        std = np.array(std)
        for i in range(num):
            failprobList = np.random.lognormal(mu, std)
            self.batchUpdateTEfprob(failprobList)
            self.solve()
            seqProb = self.sequence_prob
            for key, value in seqProb.items():
                if key not in self.MCresult['bySeq']:
                    self.MCresult['bySeq'][key] = [value]
                else:
                    self.MCresult['bySeq'][key].append(value)
            aSample = []
            for es in self.end_states:
                aSample.append(seqProb[es])
            self.MCresult['bySample'].append(aSample)


    def __str__(self):
        return 'Initial event is: ' + str(self.IE) + '. Top events are ' + str(self.top_events) + '. End states are ' + str(self.simple_es) + '. Tree Logic is ' + str(self.tree_logic)



