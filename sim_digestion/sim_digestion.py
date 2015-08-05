def get_bs_seq(seq):
    ''' Return the bisulfite converted input sequence.

    Parameters
    ----------
    seq : str
        The input DNA sequence.

    Returns
    -------
    str
        The bisulfite converted sequence.
    '''
    seq = seq.upper()

    i = 0
    while i<len(seq):
        i = seq.find("C", i)
        # Cannot find any C.
        if i<0:
            break
        # The last C or non-CpG.
        if i+1>=len(seq) or seq[i+1]!="G":
            seq = seq[:i]+"T"+seq[i+1:]
        i += 1
    return seq

def get_cutting_sites(seq, recog_seq, cut_pos, offset=0):
    ''' Get the positions of the enzyme cutting sites.

    Parameters
    ----------
    seq : str
        The DNA sequence.
    recog_seq : str
        The recognition sequence of the enzyme.
    cut_pos : int
        The cutting position of the recognition sequence. recog_seq[cut_pos:] are in the right frangment.
    offset : int
        Add the number to cutting site. It is required if split the entire sequences
        into smaller sequences.

    Return
    -------
    Set of int
        the cutting site positions of forward strand.
    Set of int
        the cutting site positions of reverse strand.
    '''

    fw_sites = set()
    rev_sites = set()
    last_site = 0
    while True:
        try:
            site = seq.index(recog_seq, last_site) + cut_pos
            fw_sites.add(site + offset)
            rev_sites.add(site + len(recog_seq) - cut_pos + offset)
            last_site = site
        except ValueError:
            break
    return fw_sites, rev_sites

def get_cutting_sites_in_parallel(seq, recog_sites, cut_poss, process_num=0, max_nt_per_process=1000000):
    ''' Do in silico digestion of a sequence and return cutting sites.

    Parameters
    ----------
    seq : str
        The sequence to run simulation.
    recog_sites : List of str
        The enzyme recognition sites.
    cut_poss : int
        The cutting positions of the recognition sites. recog_seq[cut_pos:] are in the right frangment.
        For example, MspI is C|CGG so you want to use 1 for cut_pos.
    process_num : int, optional
        The number of mutliprocessing. 0 (default) means use all CPUs.
    max_nt_per_process : int, optional
        The maximum length of sequence per process run. Default is 1,000,000.
    Returns
    -------
    Set of int
        the cutting site positions of forward strand.
    Set of int
        the cutting site positions of reverse strand.
    '''
    from multiprocessing import Pool

    if not process_num:
        p = Pool()
    else:
        p = Pool(process_num)

    results = []
    for chrom_start in range(0, len(seq)+1, max_nt_per_process):
        for recog_site, cut_pos in zip(recog_sites, cut_poss):
            # Need some overlap with previous sub sequence.
            start = max(0, chrom_start - len(recog_site))
            end = chrom_start + max_nt_per_process
            sub_seq = seq[start:end]
            results.append(p.apply_async(
                get_cutting_sites, (sub_seq, recog_site, cut_pos, start)
            ))
    p.close()
    p.join()

    fw_sites = set()
    rev_sites = set()
    for r in results:
        fw_sub_sites, rev_sub_sites = r.get()
        fw_sites = fw_sites.union(fw_sub_sites)
        rev_sites = rev_sites.union(rev_sub_sites)
    return fw_sites, rev_sites

def write_fastq_files(r1_filename, r2_filename, chrom, fragments, read_len, mode='a'):
    ''' Write reads of fragments into Fastq files.

    Parameters
    ----------
    r1_filename : str
        The read 1 filename.
    r2_filename : str
        The read 2 filename. If single-end, set this to None.
    chrom : str
        The chromosome name.
    fragments : List
        A list of start, end, strand, and fragment sequence.
    mode : str, optional
        The output file handle mode. Default: 'a'.
    read_len : int
        The read length

    '''
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio import SeqIO

    fw1 = open(r1_filename, mode)
    if r2_filename:
        fw2 = open(r2_filename, mode)
    else:
        fw2 = None

    for start, end, strand, seq in fragments:
        read1 = Seq(seq[:read_len])
        record1 = SeqRecord(
            read1,
            id="%s:%s-%s_%s"%(chrom, start, end, strand),
            letter_annotations={"phred_quality": [40]*len(read1)},
            description='1',
        )
        SeqIO.write(record1, fw1, "fastq")

        if fw2:
            read2 = Seq(seq).reverse_complement()[:read_len]
            record2 = SeqRecord(
                read2,
                id="%s:%s-%s_%s"%(chrom, start, end, strand),
                letter_annotations={"phred_quality": [40]*len(read2)},
                description='2',
            )
            SeqIO.write(record2, fw2, "fastq")
    fw1.close()
    if fw2:
        fw2.close()

def get_fragments(seq, cutting_sites, size_selection_range, is_reverse):
    ''' Run the simulated digestion and yield valid fragments.

    Parameters
    ----------
    seq : str
        The sequence.
    cutting_sites : List of int
        A list of cutting site positions.
    size_selection_range : List of int
        The range of fragment size to filter valid fragments.
    is_reverse : bool
        Wether the sequence is in the reverse strand.

    Yield
    -------
    List
        start position, end position, strand, and fragment sequence.
    '''
    from Bio.Seq import Seq

    start = 0
    for site in sorted(cutting_sites):
        frag_len = site - start

        if size_selection_range[0]<=frag_len<=size_selection_range[1]:
            fragment = seq[start:site]
            if is_reverse:
                fragment = str(Seq(fragment).reverse_complement())
                strand = '-'
            else:
                strand = '+'
            # Need to do reverse complement first.
            if bisulfite_conversion:
                fragment = get_bs_seq(fragment)
            yield start, site, strand, fragment
        start = site


def main(fasta_filename, recong_sites_list, cut_poss_list, size_selection_range, out_dir, read_len, bisulfite_conversion):
    ''' The main function for simulated digestion.

    Parameters
    ----------
    fasta_filename : str
        The Fasta filename for simulation.
    recong_sites_list : List of List of str
        The enzyme recognition sites. The first dimention are enzyme(s), which will be
        used independently. The second dimention are enzyme(s), which will be used
        together.
    cut_poss_list : List of List of int
        The cut position. For example, MspI is C|CGG. the cut postition is 1.
        This 2-D list should be the same order as recong_sites_list.
    size_selection_range : List of int
        The range of fragment size to filter valid fragments.
    out_dir : str
        The output directory
    read_len : int
        The read length of simulated reads.
    bisulfite_conversion : bool
        Whether to simulate bisulfite conversion.

    '''
    import pysam
    import os

    fastafile = pysam.Fastafile(fasta_filename)
    r1_filename = os.path.join(
        out_dir,
        '{0}_R1.fastq'.format(os.path.splitext(os.path.basename(fasta_filename))[0])
    )
    r2_filename = os.path.join(
        out_dir,
        '{0}_R2.fastq'.format(os.path.splitext(os.path.basename(fasta_filename))[0])
    )

    for chrom in fastafile.references:
        seq = fastafile.fetch(chrom).upper()
        for recong_sites, cut_poss in zip(recong_sites_list, cut_poss_list):
            fw_sites, rev_sites = get_cutting_sites_in_parallel(
                seq,
                recong_sites,
                cut_poss
            )
            fragments = list(get_fragments(
                seq,
                fw_sites,
                size_selection_range,
                is_reverse=False
            ))
            fragments += list(get_fragments(
                seq,
                rev_sites,
                size_selection_range,
                is_reverse=True
            ))

            if bisulfite_conversion:
                fragments = [(s, e, strand, get_bs_seq(f))
                    for s, e, strand, f in fragments]

            write_fastq_files(r1_filename, r2_filename, chrom, fragments, read_len)

if __name__=="__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Simulated digestion of EpiQuest.")
    parser.add_argument('fasta_filename', help='The Fata filename for simulation.')
    parser.add_argument("-t", dest="seq_type", required=True, choices=["MiniSeq", "RRBS", "MidiSeq", "RRHP", "other"], help="Choose the sequence type.")
    parser.add_argument("-s", dest='recong_sites', default=None, help="If given other in -t, specify the recognition site here. For more than one site separate sequence by comma.")
    parser.add_argument("-c", dest='cut_poss', default=None, help="If given other in -t, specify the recognition site cutting position. For more than one site separate sequence by comma.")
    parser.add_argument("-b", dest='bisulfite_conversion', action='store_true', help="Whether to do bisulfite conversion.")
    parser.add_argument("-o", dest="out_dir", default="/mnt/", help="the output dir, which will store the results. (default: /mnt/)")
    parser.add_argument("-r", dest="range_str", default="40,350", help="The gel size selection range (<min>,<max>). default: 40,350")
    parser.add_argument("-l", dest="read_len", default=50, help="The simulated read length. default: 50bp")
    args = parser.parse_args()

    if args.seq_type in ("RRBS", "MiniSeq"):
        recong_sites_list = (("CCGG", "TCGA"),)
        cut_poss_list = ((1, 1),)
        bisulfite_conversion = True
    elif args.seq_type in ("MidiSeq"):
        recong_sites_list = (("CCGG", "TTAA"), ("CCGG", "CTAG"))
        cut_poss_list = ((1, 1),)
        bisulfite_conversion = True
    elif args.seq_type in ("RRHP"):
        recong_sites_list = (("CCGG",),)
        cut_poss_list = ((1,),)
        bisulfite_conversion = False
    elif args.seq_type == 'other':
        recong_sites_list = (args.recong_sites.upper().split(','),)
        cut_poss_list = (map(int, args.cut_poss.split(',')),)
        bisulfite_conversion = args.bisulfite_conversion
    else:
        raise Exception("Unknown seq_type (%s)"%(args.seq_type))

    size_selection_range = map(int, args.range_str.split(","))

    main(
        args.fasta_filename,
        recong_sites_list,
        cut_poss_list,
        size_selection_range,
        args.out_dir,
        args.read_len,
        bisulfite_conversion
    )
