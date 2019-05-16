# TODO:
# 0. our example will be chr 6 complete, with MHC on-target
# 1. Create a reference to simulate - done
# 2. Simulate a read, but also pass it's first 100 nucleotides out, not just its fast5 file
# 3. Call the first 100 nucleotides with my HMM
# 4. Map the first 100 to reference with Mappy, output that to a SAM file
# 5. Decide on the alignment position wether to keep or reject
# 6. If read gets kept, get the rest of the read and align in to the SAM file
# 7. Each read should be quickstored as an alignstart-alignend, sorted based on alignstart
# 8. The in-ram SAM file should also be stored in a sorted order! should be a sorteddict :D  from sortedcontainers import SortedList
# 9. Periodically call GATK in a background thread (jeeeez....) to keep an up-to-date vcf of the variants
# 10. whenever background GATK finishes, call it anew
# 11. if the bam is kept in sorted order, and read group info is always kept, then a simple sam->bam conversion should not take very long...
# 12. repeat until, also use Albacore to perform the calling...
# 13. keep the decoding pomegranate hmm in memory (or pickle it or something...) to prevent having to reload and recreate the model each time we want to call stuff
# 14. look at GPU implementations of online pomegranate (probably requires some severe batching, so screw that :/)
# 15. out hmm has a nice O(T * (|S| + |E|)) complexity

from multiprocessing import Process
#import Mutator
import SimuPore
import mappy as mp #https://pypi.org/project/mappy/
from optparse import OptionParser
from time import clock
import time
import os
import pysam
import pyfaidx #https://pypi.org/project/pyfaidx/

from copy import deepcopy

import new_hmm_model #crappy pomegranate caller :/
import multiprocessing
import Queue
from subprocess import Popen


from sortedcontainers import SortedList

class Region:
	def __init__(self, chrom, start,end, name = ''):
		self.chrom = chrom
		self.start = start
		self.end = end
		self.name = name
		self.readcount = 0
		self.basecount = 0
		self.length = self.end - self.start
	def coveragedepth(self):
		return self.basecount / float(self.length)
	def __eq__(self, other):
		return self.start == other.start and self.chrom == other.chrom and self.end == other.end
	def __lt__(self, other):
		if self.chrom == other.chrom:
			if self.start == other.start:
				if self.end == other.end:
					if self.name == other.name:
						return False
					else:
						return self.name < other.name
				else:
					return self.end < other.end
			else:
				return self.start < other.start
		else:
			return self.chrom < other.chrom
	def __str__(self):
		return '%s:%d-%d/%s'%(self.chrom,self.start,self.end,self.name)

class SortableAlignment:
	def __init__(self,name,seq,qual, hit):
		self.name = name
		self.chrom = hit.ctg
		self.pos = hit.r_st
		cigar_str = hit.cigar_str
		if hit.q_st >0:  #MANUAL SOFT CLIPPING YAY
			cigar_str = '%dS'%hit.q_st + cigar_str
		if hit.q_en < len(seq):
			cigar_str = cigar_str+'%dS'%(len(seq)-hit.q_en)
		self.samstr = '%s\t%d\t%s\t%d\t%d\t%s\t*\t0\t%d\t%s\t%s\tRG:Z:Simupore\n'%(
			name,
			0 if hit.strand == '+' else 16,
			hit.ctg,
			hit.r_st+1,
			hit.mapq,
			cigar_str,
			abs(hit.q_st -hit.q_en),
			#seq[hit.q_st:hit.q_en], #only the mapped bits get stored, for some reason
			seq, #only the mapped bits get stored, for some reason
			'*' if qual == None else qual)

	def _cmp_key(self):
		return (self.chrom,self.pos)

	def __eq__(self, other):
		return self.chrom == other.chrom and self.pos == other.pos

	def __lt__(self,other):
		if self.chrom == other.chrom:
			return self.pos < other.pos
		else:
			return self.chrom < other.chrom


def Bed_to_SortedList(targetbed):
	target_regions = SortedList([])
	if os.path.exists(targetbed):
		for line in open(targetbed).readlines():
			line = line.strip().split('\t')
			try:
				if len(line)>3:
					target_regions.add(Region(line[0],int(line[1]),int(line[2]), ''.join(line[3:])))
				else:
					target_regions.add(Region(line[0],int(line[1]),int(line[2])))
			except ValueError:
				print 'Failed to parse line in bed file',targetbed,':',line
	return target_regions

def inRegion(regionlist,readregion):
	for region in regionlist:
		if region.chrom == readregion.chrom:
			if (readregion.start > region.start and readregion.start < region.end ) or (readregion.start > region.end and readregion.end < region.end ):
				return region
	return None

def samtobam(samfilename):
	samtoolspath = 'samtools'
	bamfilename = samfilename.replace('.sam','.bam')
	cmd = samtoolspath + ' view -bhS '+samfilename+' > '+bamfilename
	os.system(cmd)
	cmd = samtoolspath + ' index '+bamfilename
	os.system(cmd)
	return bamfilename

def callUG(bamfilename, reffilename,outvcfname,bedfilename):
	cmd = 'java -jar /opt/GATK_3.3-0/GenomeAnalysisTK.jar -T UnifiedGenotyper -R %s -I %s -o %s -glm BOTH -mbq 5 -deletions 0.4 -ploidy 2 -stand_call_conf 2 -L %s -dcov 1000 --defaultBaseQualities 10' % (reffilename,bamfilename, outvcfname, bedfilename)
	print cmd
	Popen(cmd.split(' '))
	##run(cmd)
	print 'Done',outvcfname
	return outvcfname

def CreateSamHeader(fastafile):
	header = '@HD\tVN:1.5\tGO:none\tSO:coordinate\n'
	print fastafile.keys()
	for region in sorted(fastafile.keys()):
		header += '@SQ\tSN:%s\tLN:%d\n' % (region, len(fastafile[region]))
	header += '@RG\tID:Simupore\tPU:Simupore\tLB:01\tSM:01\tPL:ILLUMINA\n'
	header += '@PG\tID:Simupore\tVN:Version0.01\tPN:simupore\n'
	return header

def writesamfile(sortedsamlist, samfilename, fastafile, idx=0):
	if idx == 0:
		outsamname = '%s.sam' % (samfilename)
	else:
		outsamname = '%s_%09d.sam' % (samfilename, idx)
	outsamfile = open(outsamname, 'w')
	outsamfile.write(CreateSamHeader(fastafile))
	# SortedList.update(sortedsamlist)
	for alignment in iter(sortedsamlist):
		outsamfile.write(alignment.samstr)
	outsamfile.close()
	return outsamname


def cigarlength(cigars):
	"""
	:param cigars: A list of (count,operator) tuples of a Cigar list
	:return: integer of the aligned length
	"""
	outl = 0
	for (count,op) in cigars:
		if op == 0: #M:
			outl += count
		if op == 1: #I
			outl += count
		if op == 7: #=
			outl += count
		if op == 8: #X
			outl += count
		if op == 4: #S
			outl += count
	return outl


def validate_mappy_hit (name,seq,qual,hit):
	if hit.q_st != 0 or hit.q_en != len(seq):
		#		if len(seq) != cigarlength(hit.cigar):
		print 'len seq = ', len(seq), 'cigarlength = ', cigarlength(hit.cigar)
		print alignment.samstr
		print hit
		print 'Ref start =', reffasta[hit.ctg][hit.r_st - 10:hit.r_st].seq.lower() + reffasta[hit.ctg][
																					 hit.r_st:hit.r_st + 10].seq.upper()
		print 'compare to=', '          ' + seq[0:10]
		print 'Ref end   =', reffasta[hit.ctg][hit.r_en - 10:hit.r_en].seq.upper() + reffasta[hit.ctg][
																					 hit.r_en:hit.r_en + 10].seq.lower()
		print 'compare to=', seq[-10:len(seq)]
		rs = reffasta[hit.ctg][hit.r_st:hit.r_en + 1].seq
		print rs, 'num Ns', rs.count('N')
		return -1
	else:
		return 0

def findregion(hit,regions):
	return 0

def SimProcess (Sim, Call, readindex, call_length = 1000000 ):
	t_total = clock()
	(id, template) = Sim.sample_reference()
	new_read = Sim.generatesequence(id, 'ACGTACGTACGTACGT' + 'G' * 50 + template,idstr = str(readindex))
	# call read
	basecalled, score = Call.callfast5(new_read.fast5path, max_bases= call_length)

	# align the first 500 nucs of it
	basecalled_early = basecalled[0:min(500, len(basecalled))]
	t_total = clock() - t_total
	#returns dict of readid
	result = {'simed_read':new_read,'basecalled':basecalled,'score':score,'totaltime':t_total}
	mpqueue.put(result)
	return result

class Pore:
	def __init__(self,Sim, HMM, inqueue,outqueue, id):
		self.Sim = Sim
		self.HMM = HMM
		self.inqueue = inqueue
		self.outqueue = outqueue
		self.id = id
		self.slept = 0
		self.ready = False
	def launch(self):
		print "Pore %i started"%self.id
		self.ready = True
		time.sleep(0.001)  # yield
		while (1):
			try:
				nextread = self.inqueue.get(block = False)
				time.sleep(0.001) #yield
				if type(nextread) == int:
					self.outqueue.put(self.SimProcess(nextread,max_samples=4000))
				else:
					self.outqueue.put(self.reCall(nextread))
				self.slept += 0.05
				time.sleep(0.001)
			except Queue.Empty:
				#print 'Queue empty'
				time.sleep(0.001)

	def SimProcess(self, readindex, max_samples=1000000):
		t_total = time.time()
		sampled_read  = self.Sim.sample_reference()

		time.sleep(0.001)  # yield
		new_read = self.Sim.generatesequence(sampled_read['id'], sampled_read['template'], idstr=str(readindex))

		time.sleep(0.001)  # yield
		# call read

		call_result = self.HMM.callfast5(new_read.fast5path, max_samples=max_samples)

		time.sleep(0.001)  # yield
		basecalled = call_result['decoded_sequence']
		# align the first 500 nucs of it
		t_total = time.time() - t_total
		# returns dict of readid
		result = {'basecalled': basecalled, 'score': call_result['align_score_v'], 't_total': t_total, 'poreid' : self.id, 't_sim':new_read.sim_time, 't_sample' : sampled_read['t_sample'],'readid':sampled_read['id'], 'reCalled':False, 'fast5path':new_read.fast5path, 'readlength':len(new_read.template)}
		for key in call_result.iterkeys():
			if key.startswith('t_'):
				result[key] = call_result[key]
		#print 'Pore %d finished in %fs'%(self.id, t_total)
		return result

	def reCall(self,fast5path):
		t_total = 0
		call_result = self.HMM.callfast5(fast5path,max_samples= 10000000)

		time.sleep(0.001)  # yield
		basecalled = call_result['decoded_sequence']
		# align the first 500 nucs of it
		t_total = time.time() - t_total
		# returns dict of readid
		result = {'basecalled': basecalled, 'score': call_result['align_score_v'], 't_total': t_total, 'poreid' : self.id, 't_sim':0, 't_sample' : 0, 'reCalled' :True, 'fast5path':fast5path, 'readid':fast5path}
		for key in call_result.iterkeys():
			if key.startswith('t_'):
				result[key] = call_result[key]
		#print 'Pore %d finished in %fs'%(self.id, t_total)
		return result

mpqueue = multiprocessing.Queue()
jobqueue = multiprocessing.Queue()
MT = 'v2'
if __name__ == '__main__':

	usage = "Just run it and hope for the best..."
	parser = OptionParser(usage=usage)
	parser.add_option("-f", "--fasta", action="store", type="string", dest="fasta",help="Fasta file to generate reads from", default='mutator_output_chr6_NA12878.fa')
	parser.add_option("-r", "--reference", action="store", type="string", dest="reference",help="Fasta file to use as a reference for alignment", default='chr6.fa')

	parser.add_option("-d", "--directory", action="store", type="string", dest="directory",help="Output directory where the .fast5 reads should be stored", default='simupore_output_lc19')
	parser.add_option("--outsam", action="store", type="string", dest="outsam",help="The name of the .sam file that should be created", default='lc19')
	parser.add_option("--targetbed", action="store", type="string", dest="targetbed",help="A .bed file of on target regions", default='MHC.bed')
	parser.add_option("--vcf", action="store", type="string", dest="vcf",help="Where to store the calls in .vcf format", default='lc19')
	parser.add_option("--targetcoverage", action="store", type="int", dest="targetcoverage",help="what coverage are we aiming for", default=50)
	parser.add_option("--numpores", action="store", type="int", dest="numpores",help="Number of parallel pores to simulate", default=6)



	(options,args ) = parser.parse_args()
	t0 = clock()
	maxprocess = options.numpores
	mpindexfilename = options.reference + '.mpidx'
	if os.path.exists(mpindexfilename):

		mprefseq = mp.Aligner(fn_idx_in = mpindexfilename, preset = 'map-ont')
	else:
		mprefseq = mp.Aligner(options.reference, fn_idx_out = mpindexfilename, preset = 'mp-ont')
	print "Alignment loaded in ", clock() - t0
	t0 = clock()

	UGJobs = []


	if os.path.exists(options.targetbed):
		target_regions = Bed_to_SortedList(options.targetbed)
	else:
		target_regions = []
		print "Target bed file does not exist:", options.targetbed


	outsam = open(options.outsam, 'w')

	reffasta = pyfaidx.Fasta(options.reference)

	alignments = SortedList()
	total_nucleotides = 0
	total_nucleotides_on_target = 0
	total_reads_on_target = 0
	sam_coverage_level = 0
	on_target_short = 0
	on_target_long = 0
	logtime = clock()

	Simulator = SimuPore.SimuPore(mean_read_length=4000,read_length_spread=0,noise= 1.0)
	Simulator.load_reference(options.fasta)
	Simulator.load_kmer_lookahead_time('timing_lookahead_table.tsv')
	Simulator.load_kmer_means_spreads('collected kmer mean stdvs.tsv')

	HMMCaller = new_hmm_model.HMMCaller()

	algo  = 'A'
	reads = 0
	processes = []
	while (1):
		try:
			if MT == 'v2':
				pores = []
				for i in range(maxprocess):
					#pores[-1].HMM.model = HMMCaller.model.copy()
					Simulatorprivate = SimuPore.SimuPore(mean_read_length=20000,read_length_spread=0,noise= 1.0) # add private simulators to each model, to reduce GIL overhead.
					Simulatorprivate.load_reference(options.fasta)
					Simulatorprivate.load_kmer_lookahead_time('timing_lookahead_table.tsv')
					Simulatorprivate.load_kmer_means_spreads('collected kmer mean stdvs.tsv')
					pores.append(Pore(Simulatorprivate,new_hmm_model.HMMCaller(),jobqueue,mpqueue,i))
					processes.append(multiprocessing.Process(target = pores[-1].launch , args = ()))
					jobqueue.put(reads)
					reads+=1
					processes[-1].start()
					time.sleep(0.1)
				for j in range(10):
					jobqueue.put(reads)
					reads += 1
				while (1):
					while not mpqueue.empty():
						#print "Job queue length = ",jobqueue.qsize(), "Output queue length = ", mpqueue.qsize()

						current_coverage = sum([region.coveragedepth() for region in target_regions]) /float(len(target_regions))
						t_get = clock()
						procresult = mpqueue.get()
						#print 'Get from queue time = %fs'%(clock()-t_get)
						t_align = clock()
						ontarget = False
						if 'readlength' in procresult:
							total_nucleotides += procresult['readlength']

						for hit in mprefseq.map(procresult['basecalled']):
							hit_region = Region(hit.ctg, hit.r_st, hit.r_en)
							#print 'Hit:' ,hit_region, procresult['readid'], 'dt=%.3fs' % (clock() - t_align)  # 2 milliseconds? Wow!
							readregion = inRegion(target_regions,hit_region)
							if readregion is not None:
								ontarget = True
								print 'READ ON TARGET',readregion
								readregion.basecount += len(procresult['basecalled'])
								total_nucleotides_on_target += len(procresult['basecalled'])
								print 'Hit:' ,hit_region, procresult['readid'], 'dt=%.3fs' % (clock() - t_align)  # 2 milliseconds? Wow!
								if procresult['reCalled']:
									on_target_long+=1
									alignments.add(SortableAlignment(procresult['readid'], procresult['basecalled'], None, hit))
									total_reads_on_target +=1
									print "RECALLED READ!","total_nucleotides_on_target=", total_nucleotides_on_target,'coverage = ',current_coverage
									jobqueue.put(reads)
									reads += 1
								else:
									on_target_short +=1
									jobqueue.put(procresult['fast5path'])

							else:
								alignments.add(SortableAlignment(procresult['readid'], procresult['basecalled'], None, hit))
							break
						if not ontarget:
							jobqueue.put(reads)
							reads += 1
						if jobqueue.qsize() < maxprocess:
							jobqueue.put(reads)
							reads += 1
						t_align = clock() - t_align
						timings = ''
						for k in procresult.keys():
							if k.startswith('t_'):
								timings+= ' '+ k + '='+str(procresult[k])

						#print "Rates:", procresult['t_total'], procresult['score'], procresult['poreid'], timings, 't_align_mappy=',t_align

						if int(current_coverage) > sam_coverage_level:
							sam_coverage_level = int(current_coverage)
							t0 = clock()
							smfn = writesamfile(alignments, options.outsam+'_%04d'%sam_coverage_level, reffasta)
							print "Wrote SAM in ", clock() - t0
							t0 = clock()
							bamfn = samtobam(smfn)
							print "converted Sam %s to BAM %s in %f"%(smfn,bamfn,clock() - t0)
							t0 = clock()
							callUG(bamfn, options.reference, 'test_out_%d.vcf'%sam_coverage_level,options.targetbed)

							print "Finished UG in ", clock() - t0
						if time.time()-logtime >5:
							logtime = time.time()
							print 'Stats: on_target_short=%d on_target_long=%d total_reads=%d queue_len=%d'%(on_target_short,on_target_long,reads,jobqueue.qsize())

					time.sleep(0.001)



			elif MT == 'dumb' :
				if len(processes) < maxprocess:
					reads +=1
					proc = multiprocessing.Process(target = SimProcess, args = (Simulator,HMMCaller,reads))
					processes.append(proc)
					proc.start()

				for process in processes:
					if process.exitcode is None:
						pass
					else:
						process.join()
						processes.remove(process)

				while not mpqueue.empty():
					procresult = mpqueue.get()
					t_align = clock()
					for hit in mprefseq.map(procresult['basecalled']):
						alignments.add(SortableAlignment(id,procresult['basecalled'],None,hit))
						hit_region = Region(hit.ctg,hit.r_st,hit.r_en)
						print hit_region, id, 'dt=%.3fs'%(clock()-t_align) #2 milliseconds? Wow!
						break
					print procresult['totaltime'],procresult['score']
				time.sleep(0.001)


			else:
				# simulate a read
				t_total = clock()
				(id,template) = Simulator.sample_reference()
				print 'Sim',id
				new_read = Simulator.generatesequence(id,'ACGTACGTACGTACGT'+'G'*50+template)
				#idregion = Region(id.partition(':')[0],id.partition(':')[])


				#call read
				print 'Call',
				basecalled, score = HMMCaller.callfast5(new_read.fast5path)
				print 'Simple called in %.3fs'%(clock()-t_total)
				#align the first 500 nucs of it
				basecalled_early = basecalled[0:min(500,len(basecalled))]
				print 'Align',
				for hit in mprefseq.map(basecalled_early):
					hit_region = Region(hit.ctg,hit.r_st,hit.r_en)
					print hit_region, id
				t_align = clock()
				for hit in mprefseq.map(basecalled):
					alignments.add(SortableAlignment(id,basecalled,None,hit))
					hit_region = Region(hit.ctg,hit.r_st,hit.r_en)
					print hit_region, id, 'dt=%.3fs'%(clock()-t_align) #2 milliseconds? Wow!
					break

				#ALGORITHM:
				#A: just align it
				#if algo == 'A':
				#	alignment = SortableAlignment(id, basecalled, )

				#B: keep if it hits
				#C: keep if it hits and target is not yet reached

				#determine if its a hit
				#if a hit, then add to SAM
				#when avg_coverage increases by 1, call UG!
				#slam out the bam
				#call UG as separate process!
				if reads%1000==0:
					writesamfile(alignments,options.outsam+'%d_.sam'%reads,reffasta)

		except KeyboardInterrupt:
			exit(1)


	for name,seq,qual in mp.fastx_read('simupore_references.fasta'):
		#print name
		for hit in mprefseq.map(seq):
			pass
			#print("{}\t{}\t{}\t{}".format(hit.ctg, hit.r_st, hit.r_en, hit.cigar_str))
			#print str(hit)
			alignment = SortableAlignment(name,seq,qual,hit) #10ms per 10k read :o
			alignments.add(alignment)


	print "Aligned ",len(alignments)," 3k reads in", clock() -t0
	t0 = clock()
	smfn = writesamfile(alignments, options.outsam, reffasta)
	print "Wrote SAM in ", clock() -t0
	t0 = clock()
	bamfn = samtobam(smfn)
	print "converted Sam to BAM in", clock()-t0
	t0 =  clock()
	#callUG(bamfn, options.fasta, 'test_out.vcf')

	print "Finished UG in ", clock()-t0
