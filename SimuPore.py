# TODO:
# Add events to a group in the hdf
# Simulate raw signal
# Allow for an extra amount of 'drift', e.g. lf perlin noise
# Use lookahead when calculating the times
# simulation rate: 15kbases per sec per core
# if not storing the events, then simulation rate is aproxx 150kbps per core
from numpy.random import gamma
import math
import h5py
import numpy as np
from pyfaidx import Fasta
from optparse import OptionParser
from ONTHelper import *
#import matplotlib.pyplot as plt #BREAKS REMOTE DEBUGGING!
import random
import os
import time



usage = "Raw nanopore read analyzer"
parser = OptionParser(usage=usage, version="%prog 1.0")
def addSimuporeOptions(parser):
    parser.add_option("-f", "--fasta", action="store", type="string", dest="fasta",help="Fasta file to generate reads from", default='mutator_output_chr6_NA12878.fa')
    parser.add_option("-d", "--directory", action="store", type="string", dest="directory",help="Output directory where the .fast5 reads should be stored", default='simupore_output')
    parser.add_option("--fasta_out", action="store", type="string", dest="fasta_out",help="Output file where the created reference fasta reads should be stored", default='simupore_references_lc19.fasta')
    parser.add_option("-b", "--basename", action="store", type="string", dest="basename",help="Fast5 file basename, e.g. [BASENAME]_read_NN_ch_NN_strand.fast5", default='simulated')
    parser.add_option("-l", "--mean_read_length", action="store", type="int", dest="mean_read_length",help="Mean read length to be generated from input fasta", default=10000)
    parser.add_option("-s", "--read_length_spread", action="store", type="float", dest="read_length_spread",help="Spread of the read length to be generated", default=100.0)
    parser.add_option("-c", "--count", action="store", type="int", dest="count",help="Number of reads to simulate", default=10000)
    parser.add_option("-n", "--noise", action="store", type="float", dest="noise",help="Amount of gaussian noise to add to read signal, in pA", default=1.0)
    parser.add_option("-z", "--sampling_rate", action="store", type="int", dest="sampling_rate",help="Sampling rate of the signal, in Hz", default=4000)
    parser.add_option("-r", "--traversal_rate", action="store", type="float", dest="traversal_rate",help="Traversal rate of strand through pore, in bps", default=450.0)
    parser.add_option("-k", "--kmer_table", action = "store", type = "string", dest = "kmer_table", help="The configuration parameters (kmer table) for the simulator", default = "collected kmer mean stdvs.tsv")
    parser.add_option("-t", "--timing_table", action = "store", type = "string", dest = "timing_table", help="The configuration parameters (timing table) for the simulator lookahead times", default = "timing_lookahead_table.tsv")
    parser.add_option("--revcmp", action = "store_true", dest = "revcmp", help="Simulate reverse complement reads as well", default = True)
    parser.add_option("--events", action = "store_true", dest = "events", help ="Encode events into fast5 file")

addSimuporeOptions(parser)

def CreateHDFHeader(hdffile, simulator):
    # hdffile =  h5py.File("mytestfile.hdf5", "w") #testing only
    # Raw = hdffile.create_group('Raw')
    # Raw.create_group('Reads')
    hdffile.attrs['file_version'] = 1.0

    UniqueGlobalKey = hdffile.create_group('UniqueGlobalKey')
    # only channel_id attribute is string for this group
    channel_id = UniqueGlobalKey.create_group('channel_id')
    channel_id.attrs['channel_number'] = str(simulator.channel_number)
    channel_id.attrs['digitisation'] = 8192.0
    channel_id.attrs['offset'] = 4.0
    channel_id.attrs['range'] = 1337.3
    channel_id.attrs['sampling_rate'] = float(simulator.sampling_rate)

    context_tags = UniqueGlobalKey.create_group('context_tags')
    # all attributes are strings for this group
    context_tags.attrs['experiment_type'] = 'genomic_dna'
    context_tags.attrs['fast5_output_fastq_in_hdf'] = str(1)
    context_tags.attrs['fast5_reads_per_folder'] = str(4000)
    context_tags.attrs['fastq_enabled'] = str(1)
    context_tags.attrs['fastq_reads_per_file'] = str(400)
    context_tags.attrs['filename'] = 'XXX'
    context_tags.attrs['flowcell_type'] = 'flo-min106'
    context_tags.attrs['local_basecalling'] = str(1)
    context_tags.attrs['sample_frequency'] = str(simulator.sampling_rate)
    context_tags.attrs['sequencing_kit'] = 'sqk-lsk108'
    context_tags.attrs['user_filename_input'] = str(simulator.outdir)

    tracking_id = UniqueGlobalKey.create_group('tracking_id')
    # all attributes are strings for this group
    tracking_id.attrs['asic_id'] = '2970749871'
    tracking_id.attrs['asic_id_eeprom'] = '2063468'
    tracking_id.attrs['asic_temp'] = '32.255928'
    tracking_id.attrs['auto_update'] = '1'
    tracking_id.attrs['auto_update_source'] = 'https://mirror.oxfordnanoportal.com/software/MinKNOW/'
    tracking_id.attrs['bream_core_version'] = '1.7.10.1'
    tracking_id.attrs['bream_is_standard'] = '0'
    tracking_id.attrs['bream_ont_version'] = '1.7.10.1'
    tracking_id.attrs['device_id'] = 'MN18873'  # some sort of unique flow cell ID
    tracking_id.attrs['exp_script_name'] = '28'
    tracking_id.attrs['exp_script_purpose'] = 'mux_scan'
    tracking_id.attrs['exp_start_time'] = '2017-07-20T13:50:18Z'  # datetime format
    tracking_id.attrs['flow_cell_id'] = 'r9.4simulated'
    tracking_id.attrs['heatsink_temp'] = '34.394531'
    tracking_id.attrs['hostname'] = 'your_host'
    tracking_id.attrs['installation_type'] = 'map'
    tracking_id.attrs['local_firmware_file'] = '0'
    tracking_id.attrs['operating_system'] = 'Linux 3.13.0-43-generic'
    tracking_id.attrs['protocol_run_id'] = '747861b8-18ad-44df-9240-8618e7727920'
    tracking_id.attrs['run_id'] = '4ec9ef7dd8ad20651068ca189ac03ba8e6751b89'  # should this be unique?
    tracking_id.attrs['sample_id'] = 'No sample id given'  # fill this out later with refseq?
    tracking_id.attrs['usb_config'] = '1.1.1_ONT#MinION_fpga_1.1.0#ctrl#Auto'
    tracking_id.attrs['version'] = '1.7.10'

    return hdffile

class Simulated_Read:
    def __init__(self,refid, template, fast5path):
        self.refid = refid
        self.template = template
        self.fast5path = fast5path
        self.raw_currents = None
        self.sim_time = 0




class SimuPore:
    def __init__(self,seed = 1,mean_read_length = 20000, read_length_spread = 100.0, noise = 1.0, traversal_rate = 450, fast5_prefix = 'simulated', fasta_out = 'simupore_references.fasta', do_revcmp = True,sampling_rate = 4000, events = False,outdir = 'simupore_output'):
        self.seed = seed
        random.seed(seed)
        self.readcounter = 1
        self.channel_number = 1
        self.mean_read_length = mean_read_length
        self.read_length_spread = read_length_spread
        self.noise = noise
        self.traversal_rate = traversal_rate
        self.fast5_prefix = fast5_prefix
        self.fasta_out_path = fasta_out
        self.fasta_out = open(fasta_out,'w')
        self.do_revcmp = do_revcmp
        self.sampling_rate = sampling_rate
        self.outdir = outdir
        self.events = events
    def optionsstring(self):
        return ('--fasta ' + self.refseqpath +
                ' --directory ' + self.outdir +
                ' --fasta_out '+ self.fasta_out_path +
                ' --basename '+ self.fast5_prefix +
                ' --mean_read_length %d'%(self.mean_read_length) +
                ' --read_length_spread %f'%(self.read_length_spread) +
                ' --noise %f'%(self.noise) +
                ' --sampling_rate %d'%(self.sampling_rate) +
                ' --traversal_rate %d'%(self.traversal_rate) +
                ' --kmer_table %s'%(self.kmer_means_spreads_path) +
                ' --timing_table %s'%(self.kmer_lookahead_path) +
                (' --revcmp' if self.do_revcmp else '') +
                (' --events' if self.events else '')
                )
    def load_reference(self,refseqpath):
        self.refseqpath = refseqpath
        self.referenceSequence = Fasta(refseqpath)

    def load_kmer_means_spreads(self, kmer_mean_spreads_path):
        self.kmer_means_spreads_path = kmer_mean_spreads_path
        self.kmer_means_spreads = load_kmer_mean_spreads(kmer_mean_spreads_path)

    def load_kmer_lookahead_time(self, kmer_lookahead_path):
        """

        :param kmer_lookahead_path: a path to a file containing the lookahead times for each k-mer
        :return: None
        """
        self.kmer_lookahead_path = kmer_lookahead_path
        self.kmer_lookahead_times = load_kmer_lookahead_times(kmer_lookahead_path)

    def sample_reference(self, chrom = None, start = -1, end = -1, size = -1):
        """
        Sample the reference sequence. If no params specified, it will randomly sample the submitted ref.
        :param chrom: chromosome to sample from
        :param start: start position of the read
        :param end:  end position of the read
        :param size:  size of the desired read (not implemented)
        :return: returns a dict of: {'id':id,'template':template,'t_sample':t0}
        """
        t0 = time.time()
        if chrom is None:
            refseq_id = random.choice(self.referenceSequence.keys())
        else:
            refseq_id = chrom
        if len(self.referenceSequence[refseq_id]) > (self.mean_read_length + self.read_length_spread):
            if start == -1:
                startpos = random.randint(0, len(self.referenceSequence[refseq_id]) - self.mean_read_length)
            else:
                startpos = min(len(self.referenceSequence[refseq_id]), start)
        else:
            startpos = 0
        if end == -1:
            endpos = min(len(self.referenceSequence[refseq_id]), startpos + self.mean_read_length + random.randint(-1 * self.read_length_spread, self.read_length_spread))
        else:
            endpos = min(len(self.referenceSequence[refseq_id]),end)
        template = self.referenceSequence[refseq_id][startpos:endpos].seq.upper()
        reversed = False
        if 'N' in template:
            template = template.replace('N', random.choice('ACGT'))
        if self.do_revcmp and random.random < 0.5:
            template = revcmp(template)
            reversed = True

        id = refseq_id + ':' +str(startpos) + '-' + str(endpos) +  ('_rev' if reversed else '')
        t0 = time.time()-t0
        return {'id':id,'template':template,'t_sample':t0}

    def generatesequence(self,id,template, makehdf = True, idstr = '', verbose = False): #returns Simulated_Read class
        """

        :param id: A read identifier string
        :type id: str
        :param template: The template (true) sequence to simulate
        :param makehdf:
        :param idstr:
        :param verbose:
        :return: Simulated_Read object
        """


        t0 = time.time()
        hdffilename = "%s_%s_%06d.fast5" % (self.fast5_prefix, idstr, self.readcounter)
        if verbose:
            print 'Simulating read %s' % hdffilename

        self.readcounter += 1
        hdffile = h5py.File(self.outdir + '/' + hdffilename, "w")

        thisread = Simulated_Read(id,template, self.outdir + '/' +hdffilename )

        CreateHDFHeader(hdffile,self)
        simuporegroup = hdffile.create_group('SimuPore')
        simuporegroup.attrs['Original_template_id'] = id
        simuporegroup.attrs['Original_template_seq'] = template
        simuporegroup.attrs['SimuPore_args'] = str(self.optionsstring())
        rate = self.sampling_rate / self.traversal_rate

        # scale signals by 5 to get event -> raw transformation

        scale = 0.2
        currtime = 0
        events = []  # start,length,mean,stdv
        raw_currents = []
        move = 1
        for i in range(0, len(template) - 5):
            kmer = template[i:i + 5]
            if i < len(template) - 17:
                lookahead_kmer = template[i + 5 + 8:i + 5 + 12]
                lookahead_time = self.kmer_lookahead_times[lookahead_kmer]
            else:
                lookahead_time = 1.1
            (current_level, spread) = self.kmer_means_spreads[kmer]
            # lookahead_time = 1.1
            actual_time = gamma(2, (rate * lookahead_time / 2.0))
            #floor_time = int(math.ceil(actual_time)) #TODO MAJOR CHEATS
            floor_time = int(math.floor(actual_time))
            # print actual_time,floor_time
            if floor_time > 0:
                kmer_currents = np.random.normal(current_level, spread * self.noise, (floor_time))
                kmer_utf8 = kmer
                kmer_utf8.encode('utf8')
                if self.events:
                    events.append((currtime, floor_time, np.mean(current_level), np.std(kmer_currents), kmer_utf8, move))
                move = 1
                currtime += floor_time
                raw_currents += np.ndarray.tolist(kmer_currents / scale)
                # create base duration
            else:
                move += 1
        thisread.raw_currents = raw_currents
        # add raw currents to HDF
        Raw = hdffile.create_group('Raw')
        Reads = Raw.create_group('Reads')
        Read_1 = Reads.create_group('Read_1')
        Read_1.attrs['duration'] = len(raw_currents)
        Read_1.attrs['median_before'] = 512.0
        Read_1.attrs['read_number'] = 1
        Read_1.attrs['start_mux'] = 1
        Read_1.attrs['start_time'] = 0
        Read_1.attrs['read_id'] = generate_readid()
        # write refseqs to fasta file for future reference and ease of alignment
        self.fasta_out.write('>%s\n%s\n' % (hdffilename + ' ' + Read_1.attrs['read_id'], template))
        raw_current_array = np.array(raw_currents, dtype=np.uint16)
        dataset = Read_1.create_dataset(name='Signal', shape=raw_current_array.shape, dtype='int16',
                                        data=raw_current_array)

        # add events to HDF
        #print '#Template = %d, #Raw_signal = %d, Actual traversal rate = %f'%(len(template),len(raw_currents), len(template) / ((len(raw_currents) / float(self.sampling_rate))))
        # events_array = np.array(events, dtype = 'i4, i4, d, d, S, i4')
        if self.events:
            events_array = np.array(events, dtype=[('start', '<i8'), ('length', '<i8'), ('mean', '<f8'), ('stdv', '<f8'),
                                                   ('kmer', 'S5'), ('move', '<i4')])
            simuporegroup.create_dataset('Events', data=events_array)

        hdffile.close()
        thisread.sim_time = time.time()-t0
        return thisread

if __name__ == "__main__":
    (options, args) = parser.parse_args()
    print "The following options are used:", options

    readcounter = 1
    channel_number = 1
    ReferenceSequence = Fasta(options.fasta)
    random.seed(1)

    kmer_means_spreads = load_kmer_mean_spreads(options.kmer_table)
    kmer_lookahead_times = load_kmer_lookahead_times(options.timing_table)

    Simulator = SimuPore(seed = 1, mean_read_length= options.mean_read_length, read_length_spread= options.read_length_spread, noise = options.noise, traversal_rate = options.traversal_rate, fast5_prefix = options.basename, fasta_out = options.fasta_out, do_revcmp = options.revcmp, sampling_rate = options.sampling_rate, outdir = options.directory)
    Simulator.load_reference(options.fasta)
    Simulator.load_kmer_means_spreads(options.kmer_table)
    Simulator.load_kmer_lookahead_time(options.timing_table)
    if not os.path.exists(Simulator.outdir):
        os.mkdir(Simulator.outdir)
        print 'Created output directory',Simulator.outdir
    else:
        print 'Using existing directory',Simulator.outdir

    for c in range(options.count):
        t0 = time.clock()
        gensample = Simulator.sample_reference()
        Simulator.generatesequence(gensample['id'],gensample['template'])
        print 'Simulated %s in %fs'%(gensample['id'] ,time.clock()-t0)
        t0 = time.clock()

