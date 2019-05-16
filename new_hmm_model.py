from pomegranate import *
import pomegranate.utils
import math
import h5py
import numpy as np
#import matplotlib.pyplot as plt
#import nwalign
import Aligner
from ONTHelper import *
from scipy.signal import medfilt
import cPickle
import os
import time
import pickle

class HMMCaller:
    def __init__(self,currents_means_file = 'collected kmer mean stdvs.tsv', lconst = 25.0, sigmaconst = 2.5 , median_kernel_size = 1):
        self.lconstant = lconst
        self.sigmaconstant = sigmaconst
        self.currmeans = load_kmer_mean_spreads(currents_means_file)
        self.median_kernel = median_kernel_size

        print 'lconstant', self.lconstant
        print 'sigmaconstant', self.sigmaconstant
        print 'median_kernel', self.median_kernel
        # normalize the model shiet, and div spreads by 50
        kmertotalmean = sum([v[0] for v in self.currmeans.values()]) / float(len(self.currmeans))
        kmertotalspread = sum([abs(kmertotalmean - v[0]) for v in self.currmeans.values()]) / float(len(self.currmeans))
        if False: #TODO CHEAT
            for kmer in self.currmeans.keys():
                self.currmeans[kmer] = (
                (self.currmeans[kmer][0] - kmertotalmean) / kmertotalspread, self.currmeans[kmer][1] * self.sigmaconstant)
                # currmeans[kmer] = ((currmeans[kmer][0] - kmertotalmean) / kmertotalspread,sigmaconstant)
            print len(self.currmeans)

        self.kmers = generate_k_mers(5)
        for kmer in self.currmeans.keys():
            a=self.currmeans[kmer]
            self.currmeans[kmer] = (a[0]*5,a[1]*self.sigmaconstant)

        # print kmers
        def erlang_pdf(x, l, k):  # k is always 2, l is half of the traversal rate in num samples

            # desmos code: \frac{\left(l^k\cdot x^{k-1}\cdot \exp \left(-l\cdot x\right)\right)}{\left(k-1\right)!}
            return (math.pow(l, k) * math.pow(x, k - 1) * math.exp(-1.0 * l * x)) / 2  # factorial[k - 1]

        maxdwell = 2  # longest time to spend inside a state(after this, stay in the current state)
        states = {}
        invstates = {}

        self.modelfilename = '%s_%d_%f.pickle'%(currents_means_file,lconst,sigmaconst)
        if False and os.path.exists(self.modelfilename):
            print 'Loading model from', self.modelfilename
            self.model = cPickle.load(open(self.modelfilename))
        else:

            # 1024 states
            self.model = HiddenMarkovModel()
            states2 = False
            if states2:
                for kmer in self.kmers:
                    # for dwell in xrange(maxdwell):
                    state = State(NormalDistribution(self.currmeans[kmer][0], self.currmeans[kmer][1]), name=kmer + '0')
                    # print state
                    states[(kmer, 0)] = state
                    invstates[state] = (kmer, 0)
                    self.model.add_states(state)

                    state = State(NormalDistribution(self.currmeans[kmer][0], self.currmeans[kmer][1]), name=kmer + '1')
                    # print state
                    states[(kmer, 1)] = state
                    invstates[state] = (kmer, 1)
                    self.model.add_states(state)
                for kmer in self.kmers:
                    self.model.add_transition(self.model.start, states[(kmer, 0)], 1.0 / 1024)  #ACGTA0 -> CGTAC0
                    # model.add_transition(model.end,states[(kmer,1)],1.0/1024)
                    for dwell in xrange(maxdwell):
                        state = states[(kmer, dwell)]
                        if dwell == 1:
                            # add transitions to other kmers
                            for base in 'ACGT':
                                self.model.add_transition(state, states[(kmer[1:5] + base, 0)], (1 / self.lconstant) / 4.0)
                        else:  # only stay in this one
                            self.model.add_transition(state, states[(kmer, 1)], 1 / self.lconstant)
                        # we always have an exponential stay transition
                        self.model.add_transition(state, state, 1.0 - 1 / self.lconstant)
            if states2 == False:
                for base in 'ACGT': # ADD HP States AAAAA CCCCC GGGGG TTTTT
                    kmer = base*5
                    state = State(NormalDistribution(self.currmeans[kmer][0],self.currmeans[kmer][1]),name = kmer +'1')
                    states[(kmer,1)] = state
                    self.model.add_states(state)
                for kmer in self.kmers:
                    state = State(NormalDistribution(self.currmeans[kmer][0], self.currmeans[kmer][1]), name=kmer + '0')
                    # print state
                    states[(kmer, 0)] = state
                    self.model.add_states(state)
                for kmer in self.kmers:
                    self.model.add_transition(self.model.start, states[(kmer, 0)], 1.0 / 1024) # equal probs from start to here
                    # model.add_transition(model.end,states[(kmer,1)],1.0/1024)
                    #  add transitions to other kmers
                    state = states[(kmer, 0)]
                    if kmer in ['AAAAA','CCCCC','GGGGG','TTTTT']:
                        for base in 'ACGT':
                            if base in kmer:
                                self.model.add_transition(state,states[(kmer,1)], 1.0)
                                #print state.name,'->',states[(kmer,1)].name,1.0
                                self.model.add_transition(states[(kmer,1)],states[(kmer,0)],  1 / self.lconstant / 4.0)
                                #print states[(kmer,1)].name,'->',states[(kmer,0)].name, '1/lc/4'
                                self.model.add_transition(states[(kmer,1)], states[(kmer,1)], 1.0 - 1 / self.lconstant)
                                #print states[(kmer,1)].name,'->',states[(kmer,1)].name, '1-1/lc'
                            else:
                                self.model.add_transition(states[(kmer,1)],states[(kmer[0:4]+base,0)], 1/self.lconstant /4.0)
                                #print states[(kmer,1)].name,'->',states[(kmer[0:4]+base,0)].name, '1/lc/4'


                    else:
                        for base in 'ACGT': #lconstant/4 prob to each other state
                            self.model.add_transition(state, states[(kmer[1:5] + base, 0)], (1 / self.lconstant) / 4.0)
                            self.model.add_transition(state, states[(kmer[1:5] + base, 0)], (1 / self.lconstant) / 4.0)
                        self.model.add_transition(state, state, 1.0 - 1 / self.lconstant)   #stay probabilty of 1/lconstant
                    continue
                    for dwell in xrange(maxdwell):
                        if dwell == 1:
                            pass

                        #else:  # only stay in this one
                        #                            self.model.add_transition(state, states[(kmer, 1)], 1 / self.lconstant)
                        # we always have an exponential stay transition

            print 'Baking'
            #self.model.bake(verbose = True) #sloo
            self.model.bake(merge="None", verbose=True)  # quicker
            #cPickle.dump(self.model, open(self.modelfilename,'w'))
            self.model.freeze()
            self.model.freeze_distributions()
#            k = pomegranate.
            #print "Pomegranate gpu use is:", pomegranate.utils._is_gpu_enabled() #doesnt work
            print ' Model baked...'  # ,model

    def decode_states_to_seq(self, fwbackpath,raw_signal):
        prevnode = fwbackpath[0]
        prevnodename = self.model.states[fwbackpath[0]].name
        seq = prevnodename.strip('01')
        jumps = 0
        timings = {}
        lasttime = 0
        for i, node in enumerate(fwbackpath[1:]):
            #print i, node, self.model.states[node].name, self.currmeans[self.model.states[node].name[0:5]], raw_signal[i]

            nodename = self.model.states[node].name
            if nodename[-1] == '0':
                if nodename == prevnodename:
                    pass
                else:

                    if nodename[0:4] != seq[-4:]:
                        # print 'Warning, disallowed transition happened!', transmat[prevnode,node]
                        jumps += 1
                    seq += nodename[4]
                    if i - lasttime in timings:
                        timings[i - lasttime] += 1
                    else:
                        timings[i - lasttime] = 1
                    lasttime = i

            prevnodename = nodename
            prevnode = node
        return seq, jumps, timings

    def decode_states_to_seq_quick(self, fwbackpath,raw_signal, max_bases = 1000000):
        hpnodes = []
        for i,state in enumerate(self.model.states):
            if state.name.endswith('1'):
                hpnodes.append(i)
        prevnodeid = fwbackpath[0]
        prevnodename = self.model.states[fwbackpath[0]].name
        seq = prevnodename.strip('01')
        jumps = []
        homopolymer_time_index = 0
        dbgstates = []
        for i,nodeid in enumerate(fwbackpath[1:]):
            #print i, node, self.model.states[node].name, self.currmeans[self.model.states[node].name[0:5]], raw_signal[i]
            if nodeid != prevnodeid:
                if len(seq)>= max_bases:
                    break
                nodename = self.model.states[nodeid].name
                if prevnodeid in hpnodes: #meaning we transitioned out of a homopolymer node:
                    if len(seq) >0 and i-homopolymer_time_index> 0:
                        traversal_rate = i/len(seq)
                        if traversal_rate > 0:
                            seq += self.model.states[prevnodeid].name[0]*int(math.floor((i-homopolymer_time_index)/traversal_rate ))
                if nodename == 'GGGGG0' or nodename == 'GGGGG1':
                    prevnodename == self.model.states[prevnodeid].name
                    pass
                if nodeid not in hpnodes:
                    dbgstates.append(nodename)
                    if len(dbgstates)> 20:
                        dbgstates.pop(0)
                    if nodename[0:4] != seq[-4:]:
                        # print 'Warning, disallowed transition happened!', transmat[prevnode,node]
                        jumps.append((seq[-20:],nodename))
                    seq += nodename[4]
                else:
                    homopolymer_time_index = i
                    prevnodeid = nodeid
                    pass

                prevnodeid = nodeid
        return seq, jumps

    def decode_states_to_seq_map(self, fwbackpath,raw_signal, max_bases = 1000000):
        hpnodes = []
        for i,state in enumerate(self.model.states):
            if state.name.endswith('1'):
                hpnodes.append(i)
        prevnodeid = fwbackpath[0]
        prevnodename = self.model.states[fwbackpath[0]].name
        seq = prevnodename.strip('01')
        jumps = []
        homopolymer_time_index = 0
        dbgstates = []
        for i,nodeid in enumerate(fwbackpath[1:]):
            #print i, node, self.model.states[node].name, self.currmeans[self.model.states[node].name[0:5]], raw_signal[i]
            if nodeid != prevnodeid:
                if len(seq)>= max_bases:
                    break
                nodename = self.model.states[nodeid].name
                if prevnodeid in hpnodes: #meaning we transitioned out of a homopolymer node:
                    traversal_rate = i/len(seq)
                    seq += self.model.states[prevnodeid].name[0]*int(math.floor((i-homopolymer_time_index)/traversal_rate ))
                if nodename == 'GGGGG0' or nodename == 'GGGGG1':
                    prevnodename == self.model.states[prevnodeid].name
                    pass
                if nodeid not in hpnodes:
                    dbgstates.append(nodename)
                    if len(dbgstates)> 20:
                        dbgstates.pop(0)
                    if nodename[0:4] != seq[-4:]:
                        # print 'Warning, disallowed transition happened!', transmat[prevnode,node]
                        jumps.append((seq[-20:],nodename))
                    seq += nodename[4]
                else:
                    homopolymer_time_index = i
                    prevnodeid = nodeid
                    pass

                prevnodeid = nodeid
        return seq, jumps

    def callfast5(self, fast5path, truereference=None, max_samples = 1000000 ,algorithm = 'map',verbose = False, nwalign = False):
        fast5file = h5py.File(fast5path,'r+') #open in read and write!
        raw_signal = None
        for key in fast5file['Raw/Reads'].keys():
            if 'Read_' in key:
                raw_signal = np.copy(fast5file['Raw/Reads'][key]['Signal'])
        original_template_seq = None
        if 'SimuPore' in fast5file:
            original_template_seq = fast5file['SimuPore'].attrs['Original_template_seq']

        if raw_signal is None:
            print 'Failed to load raw signal from read'
            return None

        if self.median_kernel != 1:
            raw_signal = medfilt(raw_signal, self.median_kernel)

        sigmean = np.mean(raw_signal)
        sigdev = np.mean(np.abs(raw_signal - sigmean))
        #raw_signal = ((raw_signal - sigmean) / sigdev)*1.0 #TODO CHEAT
        #fwbackpath = self.model.predict(raw_signal, algorithm='viterbi')
        t_trace = time.clock()
        fwbackpath = self.model.predict(raw_signal[0:min(max_samples,raw_signal.shape[0])], algorithm='viterbi')[1:]
        t_trace = time.clock() - t_trace
        t_decode = time.clock()
        decoded_sequence, jumps = self.decode_states_to_seq_quick(fwbackpath, raw_signal, max_bases=100000)
        t_decode = time.clock() - t_decode
        t_align = time.clock()
        if nwalign:
            align_score_v = Aligner.parasail_nw_onlyscore(decoded_sequence, original_template_seq)
        else:
            align_score_v = -1
        if verbose:
            print 'Trace: %.3fs, Decode: %.3fs, Align: %.3fs' % (t_trace, t_decode, t_align)

        #add call result to fast5
        try:
            fast5file[u'Analyses/Basecall_1D_000/BaseCalled_template']
        except KeyError:
            Analyses = fast5file.create_group(u'Analyses/Basecall_1D_000/BaseCalled_template')
            #BaseCall_1D_000 = Analyses.create_group('BaseCall_1D_000')
            #BaseCalled_complement =  BaseCall_1D_000.create_group('BaseCalled_complement')
        fast5file[u'Analyses/Basecall_1D_000/BaseCalled_template'].attrs['Fastq'] = str('@%s\n%s\n'%(fast5path,decoded_sequence))
        fast5file.close()



        t_align = time.clock() - t_align
        return {'decoded_sequence':decoded_sequence,'align_score_v':align_score_v,'t_align':t_align,'t_trace':t_trace,'t_decode':t_decode}


        ##############################Below is dead code for MAP decoding######################
        t_trace = time.clock()
        logp, fwbackpath = self.model.maximum_a_posteriori(raw_signal)
        statesmatrix = self.model.predict_proba(raw_signal)
        fwbackpath = numpy.argmax(statesmatrix, 1)
        fwbackprobs = numpy.amax(statesmatrix, 1)
        t_trace = time.clock() - t_trace
        t_decode = time.clock()
        decoded_sequence, jumps = self.decode_states_to_seq_quick(fwbackpath, raw_signal, max_bases=max_bases)
        t_decode = time.clock() - t_decode
        t_align = time.clock()
        align_score_m = Aligner.parasail_nw_onlyscore(decoded_sequence, original_template_seq)
        t_align=time.clock() - t_align
        print 'Trace: %.3fs, Decode: %.3fs, Align: %.3fs' % (t_trace, t_decode, t_align)
        print 'v',align_score_v
        print 'm',align_score_m
        return decoded_sequence,align_score_v

        if algorithm == 'viterbi':
            fwbackpath = self.model.predict(raw_signal, algorithm='viterbi')[1:] # just returns a list of state IDs
            #fwbackpath = self.model.viterbi(raw_signal)[1:] # returns a full list of tuples of states, not just IDs (a bit heavy)
        elif algorithm == 'map':
            logp, fwbackpath =  self.model.maximum_a_posteriori(raw_signal)
            statesmatrix =  self.model.predict_proba(raw_signal)
            fwbackpath = numpy.argmax(statesmatrix,1)
            fwbackprobs = numpy.amax(statesmatrix,1)

            #logp: double The log p of the sequence under the viterbi path
            #path: list of Tuples of(state index, state object) of the posterior path.

            for i in range(100):
                print fwbackpath[i]
        t_trace = time.clock()-t_trace

        t_decode = time.clock()
        if algorithm == 'viterbi':
            decoded_sequence, jumps = self.decode_states_to_seq_quick(fwbackpath, raw_signal, max_bases = max_bases)
        elif algorithm == 'map':
            decoded_sequence, jumps = self.decode_states_to_seq_map(fwbackpath, raw_signal, max_bases = max_bases)
        align_score = 0
        t_decode = time.clock() - t_decode

        t_align = time.clock()
        if original_template_seq != None:
            align_score = Aligner.parasail_nw_onlyscore(decoded_sequence, original_template_seq)
        if truereference != None:
            align_score = Aligner.parasail_nw_onlyscore(decoded_sequence, truereference)
        t_align = time.clock() - t_align
        print 'Trace: %.3fs, Decode: %.3fs, Align: %.3fs'%(t_trace,t_decode,t_align)
        return decoded_sequence, align_score

def revcmp(seq):
    seq_dict = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    return "".join([seq_dict[base] for base in reversed(seq)])

if __name__ == "__main__":
    ref = revcmp(''.join(map(str.strip,open('read447_lomanlab.fasta').readlines()[1:])))
    resultsfile = open('new_hmm_map_model_results.tsv','w')

    # currmeans=(u, sigma for each kmer)
    currmeans = load_kmer_mean_spreads('collected kmer mean stdvs.tsv')
    #load kmer db:

    lconstant = 25.0
    sigmaconstant = 0.2
    ##iterate from here
    for median_kernel in [1,3,1,5]:
        for lconst in [100.0,5.0,7.0,10.0,15.0,25.0,30.0,35.0]:
            for sconst in [0.3,0.135,0.125,0.14,0.2,0.25,0.3,0.4,0.5]:
                lconstant= lconst
                sigmaconstant = sconst
                print 'lconstant',lconstant
                print 'sigmaconstant',sigmaconstant
                print 'median_kernel',median_kernel
                #normalize the model shiet, and div spreads by 50
                kmertotalmean = sum([v[0] for v in currmeans.values()])/float(len(currmeans))
                kmertotalspread = sum([abs(kmertotalmean-v[0]) for v in currmeans.values()]) /float(len(currmeans))
                for kmer in currmeans.keys():
                    currmeans[kmer] = ((currmeans[kmer][0] - kmertotalmean) / kmertotalspread,currmeans[kmer][1]*sigmaconstant)
                    #currmeans[kmer] = ((currmeans[kmer][0] - kmertotalmean) / kmertotalspread,sigmaconstant)
                print len(currmeans)

                kmers = generate_k_mers(5)

                #print kmers
                def erlang_pdf(x, l, k): #k is always 2, l is half of the traversal rate in num samples

                    # desmos code: \frac{\left(l^k\cdot x^{k-1}\cdot \exp \left(-l\cdot x\right)\right)}{\left(k-1\right)!}
                    return (math.pow(l, k) * math.pow(x, k - 1) * math.exp(-1.0 * l * x)) / 2# factorial[k - 1]

                maxdwell = 2 #longest time to spend inside a state(after this, stay in the current state)
                states = {}
                invstates = {}

                # 1024 states
                model = HiddenMarkovModel()
                states2 = False
                if states2:
                    for kmer in kmers:
                        #for dwell in xrange(maxdwell):
                            state = State(NormalDistribution(currmeans[kmer][0],currmeans[kmer][1]), name = kmer+'0')
                            #print state
                            states[(kmer, 0)] = state
                            invstates[state] = (kmer,0)
                            model.add_states(state)

                            state = State(NormalDistribution(currmeans[kmer][0], currmeans[kmer][1]),  name = kmer+'1')
                            # print state
                            states[(kmer, 1)] = state
                            invstates[state] = (kmer, 1)
                            model.add_states(state)
                    for kmer in kmers:
                        model.add_transition(model.start,states[(kmer,0)],1.0/1024)
                        #model.add_transition(model.end,states[(kmer,1)],1.0/1024)
                        for dwell in xrange(maxdwell):
                            state =   states[(kmer, dwell)]
                            if dwell == 1:
                                # add transitions to other kmers
                                for base in 'ACGT':
                                    model.add_transition(state, states[(kmer[1:5]+base, 0)],(1/lconstant)/4.0 )
                            else: #only stay in this one
                                model.add_transition(state,  states[(kmer, 1)], 1/lconstant)
                            #we always have an exponential stay transition
                            model.add_transition(state,  state, 1.0 - 1/lconstant)
                if states2 == False:
                    for kmer in kmers:
                        #for dwell in xrange(maxdwell):
                            state = State(NormalDistribution(currmeans[kmer][0],currmeans[kmer][1]), name = kmer+'0')
                            #print state
                            states[(kmer, 0)] = state
                            invstates[state] = (kmer,0)
                            model.add_states(state)

                            state = State(NormalDistribution(currmeans[kmer][0], currmeans[kmer][1]*2),  name = kmer+'1')
                            # print state
                            states[(kmer, 1)] = state
                            invstates[state] = (kmer, 1)
                            model.add_states(state)
                    for kmer in kmers:
                        model.add_transition(model.start,states[(kmer,0)],1.0/1024)
                        #model.add_transition(model.end,states[(kmer,1)],1.0/1024)
                        for dwell in xrange(maxdwell):
                            state =   states[(kmer, dwell)]
                            if dwell == 1:
                                # add transitions to other kmers
                                for base in 'ACGT':
                                    model.add_transition(state, states[(kmer[1:5]+base, 0)],(1/lconstant)/4.0 )
                            else: #only stay in this one
                                model.add_transition(state,  states[(kmer, 1)], 1/lconstant)
                            #we always have an exponential stay transition
                            model.add_transition(state,  state, 1.0 - 1/lconstant)

                print 'Baking'
                #model.bake(verbose = True) #sloo
                model.bake(merge = "None", verbose = True) #quicker
                print ' Model baked...'#,model











                #model.plot()
                transmat =  model.dense_transition_matrix()


                matrixf = open('transition_matrix.tsv','w')
                matrixf.write('P\t%s\n'%('\t'.join([node.name for node in model.states])))
                for row in range( transmat.shape[0]):
                    matrixf.write('%s\t%s\n'%(model.states[row].name,'\t'.join(map(str, transmat[row,:]))))
                fast5file = h5py.File('LomanLabz_PC_20160923_FNFAB39043_MN17250_sequencing_run_Human_1D_ligation_R9_4_61260_ch391_read447_strand.fast5')
                raw_signal = None
                for key in fast5file['Raw/Reads'].keys():
                    if 'Read_' in key:
                        raw_signal = np.copy(fast5file['Raw/Reads'][key]['Signal'])
                #normalize raw signal:
                old_raw = np.copy(raw_signal)

                for devtest in range(5):

                    sigstart = 10000
                    sigend = 12000
                    raw_signal = medfilt(raw_signal,median_kernel)
                    raw_signal = np.copy(old_raw[sigstart:sigend])
                    sigmean = np.mean(raw_signal)
                    sigdev = np.mean(np.abs(raw_signal-sigmean))
                    raw_signal = ((raw_signal - sigmean)/sigdev ) * (0.80+0.05*devtest)
                    unfiltered = np.copy(raw_signal)
                    #print raw_signal

                    fwbackpath =  model.predict(raw_signal,algorithm = 'viterbi')[1:]
                    #fwbackpath =  model.predict(raw_signal,algorithm = 'map')
                    #print fwbackpath
                    #for i,node in enumerate(fwbackpath):
                    #    print i,node, model.states[node].name, invstates[model.states[node]],currmeans[model.states[node].name[0:5]],raw_signal[i]

                    def decode_states_to_seq(path):
                        prevnode = path[0]
                        prevnodename = model.states[path[0]].name
                        seq = prevnodename.strip('01')
                        jumps = 0
                        timings = {}
                        lasttime = 0
                        for i, node in enumerate(path[1:]):
                            #print i, node, model.states[node].name, invstates[model.states[node]], currmeans[model.states[node].name[0:5]], \
                            raw_signal[i]

                            nodename = model.states[node].name
                            if nodename[-1] == '0':
                                if nodename == prevnodename:
                                    pass
                                else:

                                    if  nodename[0:4] != seq[-4:]:
                                        #print 'Warning, disallowed transition happened!', transmat[prevnode,node]
                                        jumps+=1
                                    seq+=nodename[4]
                                    if i-lasttime in timings:
                                        timings[i-lasttime] += 1
                                    else:
                                        timings[i-lasttime] = 1
                                    lasttime = i

                            prevnodename = nodename
                            prevnode = node
                        return seq,jumps,timings

                    def get_transition_times(m, path):
                        times = {}
                        prevkmer = m.states[path[0]].name[0:5]
                        prevstart = 0
                        for node in path:
                            return
                    decsec,jumps,timings = decode_states_to_seq(fwbackpath)
                    print [(key,timings[key]) for key in sorted(timings.keys())]
                    maxscore =  Aligner.align_local_maxscore(ref[:1200],decsec,match_score= 1)
                    resultstr = '%f_%f_%d_%d_%d_%d'%(lconstant,sigmaconstant,median_kernel,maxscore,len(decsec),jumps)
                    resultsfile.write('%s\t%s\n'%(resultstr.replace('_','\t'),decsec))
                    print maxscore,jumps,len(decsec),decsec
                    colormap = {'A':'b','C':'r','G':'g','T':'y'}
                    '''
                    if devtest == 2:

                        plt.clf()
                        fig = plt.figure(figsize=((sigend-sigstart)/300,4))
                        plt.plot(range(len(unfiltered)),raw_signal,'m-')
                        plt.plot(range(len(raw_signal)),raw_signal,'c-',linewidth = 1)
                        plt.scatter(range(sigend-sigstart),[model.states[node].distribution.parameters[0] for node in fwbackpath],s = 1, marker = '.',color = [colormap[model.states[node].name[4]] for node in fwbackpath])
                        for base,color in [('A','b'),('C','r'),('G','g'),('T','y')]:

                            spreadlow=np.zeros((sigend-sigstart))
                            spreadhigh=np.zeros((sigend-sigstart))

                            for i,node in enumerate(fwbackpath):
                                currkmer = model.states[node].name[0:5]
                                distparams = states[(currkmer[0:4]+base,0)].distribution.parameters
                                spreadlow[i] = distparams[0] - distparams[1] #uses stdev instead of squared!
                                spreadhigh[i] = distparams[0] + distparams[1]
                            plt.plot(range(sigend-sigstart),spreadhigh,color+'.',markersize = 1)
                            plt.plot(range(sigend-sigstart),spreadlow,color+'.',markersize = 1)
                        plt.xlim([0,sigend-sigstart])
                        plt.savefig('hmm_map_%s.png'%resultstr,dpi = 300)
                        # plt.show()
                        '''

    #print model.predict_proba(raw_signal[500:600])

    #this doesnt work, too many states, try to remodel:

    '''for kmer in kmers:
        for dwell in xrange(maxdwell):
            state = State(NormalDistribution(currmeans[kmer][0],currmeans[kmer][1]))
            #print state
            states[(kmer, dwell)] = state
            invstates[state] = (kmer,dwell)
            model.add_states(state)
    for kmer in kmers:
        model.add_transition(model.start,states[(kmer,0)],1.0/1024)
        model.add_transition(model.end,states[(kmer,0)],1.0/1024)
        for dwell in xrange(maxdwell):
            state =    states[(kmer, dwell)]
            if dwell < maxdwell -1:
            # add transitions to other kmers
                for base in 'ACGT':
                    model.add_transition(state,
                        states[(kmer[0:4]+base, 0)],
                        erlang_pdf(dwell, lconstant,2 ) / 4)
                # add stay transition:
                model.add_transition(state, states[(kmer, dwell + 1)], 1.0 - erlang_pdf(dwell,lconstant,2))
            else:
                for base in 'ACGT':
                    model.add_transition(state,
                                         states[(kmer[0:4] + base, 0)],
                                        0.25)'''