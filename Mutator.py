# TODO: load reffasta, and generate phased, mutated chromosomes
# optionally index the chromosomes for mappy
import argparse
from pyfaidx import Fasta
import vcf
import textwrap



def mutate(vcfname, fastaname, outputname):
	refseq = Fasta(fastaname)
	outputfasta = {}
	vcf_reader = vcf.Reader(open(vcfname,'r'))
	for chromosome in refseq.keys():
		cm = chromosome + "_maternal"
		cp = chromosome + "_paternal"
		outputfasta[cm] = [] #maternal is first strand, paternal is second strand
		outputfasta[cp] = [] #join them later to avoid very large concatenations of strings...
		current_reference_position = 0
		for record in vcf_reader.fetch(chromosome): # this returns a 0 based system...
			call = record.genotype(args.sample)
			if call.gt_type > 0:
				refseq_refbases = refseq[chromosome][record.POS -1 :record.POS - 1 + len(record.REF )].seq
				if record.REF != refseq_refbases.upper():
					print record, 'Sanity check failed record.ref=', record.REF,' is not refseq=', refseq_refbases,call,call.gt_bases
					continue
				prefix_sequence = refseq[chromosome][current_reference_position:record.POS -1].seq
				bases = call.gt_bases.split('|') if call.phased else call.gt_bases.split('/')
				for base in bases:
					for char in base:
						if char not in 'ACGT':
							print 'Another sanity check failed, not ACGT',record,call,call.gt_bases
							continue
				outputfasta[cm].append(prefix_sequence+bases[0])
				outputfasta[cp].append(prefix_sequence+bases[1])
				current_reference_position = record.POS -1 + len(record.REF)
		suffix_sequence = refseq[chromosome][current_reference_position:0].seq
		outputfasta[cm].append(suffix_sequence)
		outputfasta[cp].append(suffix_sequence)

	outfile = open(args.output,'w')


	for chromosome in sorted(outputfasta.keys()):
		seq = ''.join(outputfasta[chromosome])
		print chromosome,len(seq)
		outfile.write('>'+chromosome+'\n')
		k = len(seq)
		for i in xrange(0,len(seq),70):
			outfile.write(seq[i:min(k,i+70)]+'\n')
			#if i %10000==0:
			#	print i

	#for chromosome in sorted(outputfasta.keys()):
	#	outfile.write('>'+chromosome+'\n'+textwrap.wrap(''.join(outputfasta[chromosome]),70)+'\n')
	outfile.close()






if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("-v", "--vcf", help="Input VCF file", default='NA12878.vcf.gz')
	parser.add_argument("-f", "--fasta", help="Input FASTA file", default='chr6.fa')
	parser.add_argument("-o", "--output", help="Output FASTA file", default='mutator_output_chr6_NA12878.fa')
	parser.add_argument("-s", "--sample", help="sample ID in vcf file", default='NA12878')
	# parser.add_argument("-i", help = "Index the results for Mappy", action = "store_true")

	args = parser.parse_args()
	print args
	mutate(args.vcf,args.fasta,args.output)