import numpy as np
import multiprocessing
import pysam

#TODO: Add a flag that allows for binary counting (was there an insertion at this pos?) or add 1 for each inserted base.

def hash_readid(read_id, threads, thread_index):
    '''
    returns True if this read should be processed
    '''
    if int(read_id[:8], 16) % threads == thread_index:
        return True

    else:
        return False

def initialize_numpy_array(header):
    '''
    Create a dictionary using reference names as keys that holds a 2D numpy array initialized to 0 with length
    equal to the reference length and rows for match, mismatch, insertion, deletion.

    Pileup engine has a maximum per position coverage of 4294967295. Results will be skewed if that depth is
    exceeded.
    '''
    ref_arrays = {
    name: np.zeros((4, length), dtype=np.uint32)
    for name, length in zip(header.references, header.lengths)
    }

    return ref_arrays
    

def _pileup(bam,
            ref_path,
            threads,
            thread_index):

    ref_file = pysam.FastxFile(ref_path)
    ref_seqs = {}
    for seq in ref_file:
        ref_seqs[seq.name] = seq.sequence
    
    alignment_file = pysam.AlignmentFile(bam)
    ref_array = initialize_numpy_array(alignment_file.header)

    for read in alignment_file:
        if read.is_unmapped or read.is_secondary or read.is_supplementary or read.is_reverse:
            continue
        if not hash_readid(read.query_name, threads, thread_index):
            continue

        rpos = read.reference_start
        qpos = 0
        ref_name = read.reference_name
        ref_seq = ref_seqs[ref_name]
        q_seq = read.query_sequence
        arr = ref_array[ref_name]
        
        for op, length in read.cigartuples:
            if op == 7:  # =
                arr[0, rpos:rpos+length] += 1
                qpos += length
                rpos += length
            elif op == 8:  # X
                arr[1, rpos:rpos+length] += 1
                qpos += length
                rpos += length
            elif op == 0:  # M — ambiguous, must check
                for i in range(length):
                    if q_seq[qpos] == ref_seq[rpos]:
                        arr[0, rpos] += 1
                    else:
                        arr[1, rpos] += 1
                    qpos += 1
                    rpos += 1
            elif op == 1:  # I
                arr[2, rpos] += 1
                qpos += length
            elif op == 2:  # D
                arr[3, rpos:rpos+length] += 1
                rpos += length
            elif op == 4:  # S
                qpos += length
            elif op == 3:  # N
                rpos += length

    return ref_array
    
def _pileup_wrapper(args):
    bam, ref_path, threads, thread_index = args
    
    return _pileup(bam,
                   ref_path,
                   threads,
                   thread_index)

def _merge_arrays(arrays):
    template_array = arrays[0]
    for arr_dict in arrays[1:]:
        for key in arr_dict:
            template_array[key] += arr_dict[key]
    return template_array

def pileup(bam, ref_path, threads):

    args = [(bam, ref_path, threads, i) for i in range(threads)]

    with multiprocessing.Pool(threads) as p:
        arrays = p.map(_pileup_wrapper, args)

    return _merge_arrays(arrays)
        

        
    

    