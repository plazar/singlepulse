import os.path
import sys

import scipy.signal
import numpy as np

import presto
import infodata
import cands
import utils
import errors

USEBADBLOCKS = False

# Define some constants
FFTLEN = 8192     # Should be a power-of-two for best speed
CHUNKLEN = 8000   # FFTLEN - CHUNKLEN must be larger than largest downfactor
DETRENDLEN = 1000 # length of a linear piecewise chunk of data for detrending
BLOCKS_PER_CHUNK = CHUNKLEN / DETRENDLEN
OVERLAP = (FFTLEN - CHUNKLEN)/2
WORKLEN = CHUNKLEN + 2*OVERLAP  # currently it is FFTLEN...

DEFAULT_DOWNFACTS = [1,2,3,4,6,9,14,20,30,45,70,100,150]

# Define caches for kernels 
kern_cache = utils.KernelCache(utils.make_kern)
fftd_kern_cache = utils.KernelCache(lambda x: utils.make_fftd_kern(x, WORKLEN))


def find_single_pulses(datfn, downfacts, threshold, search_badblocks=False):
    """Given the path to a *.dat file and a list of downsampling
        factors, search for single pulses.

        Inputs:
            datfn: filename of a *.dat file. The corresponding
                    *.inf file must also exist.
            downfacts: A list of downsampling factors used for
                        smoothing timeseries data.
            threshold: sigma value above which events are
                        considered significant.
            search_badblocks: a boolean value, if True blocks with
                        excessively low or high sigmas are _not_ set
                        to 0 before searching. 
                        (Default: False, null out outlier blocks).

        Output:
            candlist: A CandidateList object containing the pulses 
                        found in the *.dat file searched.
    """
    timeseries, info = read_timeseries(datfn)
    timeseries, bad_blocks = prep_timeseries(timeseries)

    roundN = timeseries.size
    numchunks = roundN/CHUNKLEN
    
    # Mask off regions
    offregion_mask = np.zeros(timeseries.shape, dtype=np.bool)
    if info.breaks:
        offregions = zip([x[1] for x in info.onoff[:-1]],
                         [x[0] for x in info.onoff[1:]])
        for off, on in offregions:
            offregion_mask[off:on] = True

    blocknums = np.arange(len(timeseries))/DETRENDLEN
    bad_block_mask = np.in1d(blocknums, bad_blocks)
    
    # Mask bad blocks, or not
    if search_badblocks:
        tomask = np.bitwise_or(offregion_mask, bad_block_mask)
    else:
        tomask = offregion_mask
    masked_timeseries = np.ma.masked_array(timeseries, mask=tomask, fill_value=0.0)

    # Step through the data and search each chunk
    print "  Downfacts being used (bins): ", downfacts
    candlist = cands.CandidateList()
    for ii, chunk in enumerate(chunk_up_timeseries(masked_timeseries, \
                                                    numchunks, info)):
        for downfact in downfacts:
            tmpcandlist = chunk.search(threshold, downfact)
            tmpcandlist.add_info(info)
            tmpcandlist.prune_related1()
            candlist += tmpcandlist

        sys.stdout.write("\r  Now searching... (%6.2f%%)" % \
                            (float(ii+1)/(numchunks)*100))
        sys.stdout.flush()
    sys.stdout.write("\n")

    # A little pruning
    candlist.prune_related2(np.max(downfacts))
    candlist.prune_border_cases()
    if search_badblocks:
        print "  Found %d pulse candidates (%d are in bad blocks)" % \
                (len(candlist), len([c for c in candlist if c.from_badblock]))
    else:
        print "  Found %d pulse candidates" % len(candlist)
    
    # Finally return the candidates collected
    return candlist


def read_timeseries(datfn):
    """Read a dat file.
        Return the timeseries rounded down to the nearest DETRENLEN
        and the associated infodata object.
    """
    index = datfn.rfind(".dat")
    if index == -1:
        raise errors.SinglePulseError("File is not a *.dat file (%s)!" % datfn)
    filenmbase = datfn[:index]
    info = infodata.infodata(filenmbase+".inf")

    N = int(info.N)
    roundN = N/DETRENDLEN * DETRENDLEN
    numchunks = roundN/CHUNKLEN
    
    # Read in the file
    print 'Reading "%s"...' % os.path.split(datfn)[-1]
    timeseries = np.fromfile(datfn, dtype=np.float32, count=roundN)
    return timeseries, info 


def prep_timeseries(timeseries):
    """Prepare a timeseries for search.
        1) Detrend block-wise
        2) Compute statistics of each block
        3) Normalise blocks using stdev of block
        Return detrended, normalised timeries and list of bad blocks
    """
    # Split the timeseries into chunks for detrending
    numblocks = timeseries.size/DETRENDLEN
    timeseries.shape = (numblocks, DETRENDLEN)

    # Detrend
    print '  De-trending the data and computing statistics...'
    timeseries = np.ascontiguousarray(detrend_blocks(timeseries))

    # Calculate statistics for each block
    stds, std_stds, median_stds = calc_stats(timeseries)
    print "    pseudo-median block standard deviation = %.2f" % (median_stds)

    # Determine bad blocks
    lo_std = median_stds - 4.0 * std_stds
    hi_std = median_stds + 4.0 * std_stds
    bad_blocks = np.flatnonzero((stds < lo_std) | (stds > hi_std))
    print "    identified %d bad blocks out of %d (i.e. %.2f%%)" % \
          (len(bad_blocks), len(stds),
           100.0*float(len(bad_blocks))/float(len(stds)))
    
    # Normalise data and reshape
    stds[bad_blocks] = median_stds
    timeseries /= stds[:, np.newaxis]
    timeseries.shape = (numblocks*DETRENDLEN,)

    return (timeseries, bad_blocks)


def chunk_up_timeseries(timeseries, numchunks, info, skipchunks=0):
    """Generator function that yields Chunk objects
        made up of segments of the timeseries.

        Inputs:
            timeseries: a numpy array of timeseries data 
                        that will be chunked up.
            numchunks: the number of chunks to make.
            info: infodata object storing information from *.inf file.
            skipchunks: number of chunks to skip. (Default: 0)

            Notes: - Chunks are of size 'CHUNKLEN'
                   - 'skipchunks' > 0 means fewer than numchunks
                        will be yielded.
    """
    for chunknum in range(skipchunks, numchunks):
        yield get_chunk(timeseries, chunknum, info) 


def get_chunk(timeseries, chunknum, info):
    """Return a Chunk object corresponding to the
        chunknum'th chunk in timeseries.

        Inputs:
            timeseries: a numpy array of timeseries data 
                        that will be chunked up.
            chunknum: the chunk number to return.
            info: infodata object storing information from *.inf file.
        Output:
            chunk: a Chunk object.
    """
    masked_timeseries = np.ma.masked_array(timeseries)
    numchunks = timeseries.size/CHUNKLEN
    loind = chunknum*CHUNKLEN-OVERLAP
    hiind = (chunknum+1)*CHUNKLEN+OVERLAP
    # Take care of beginning and end of file overlap issues
    if (chunknum==0): # Beginning of file
        chunk = np.ma.zeros(WORKLEN, dtype=np.float32)
        chunk[OVERLAP:] = masked_timeseries[loind+OVERLAP:hiind]
    elif (chunknum==numchunks-1): # end of the timeseries
        chunk = np.ma.zeros(WORKLEN, dtype=np.float32)
        chunk[:-OVERLAP] = masked_timeseries[loind:hiind-OVERLAP]
    else:
        chunk = masked_timeseries[loind:hiind]
    return Chunk(chunk, loind+OVERLAP, info)
    


def detrend_blocks(blocks, usefast=False):
    """Detrend blocks of data using either a fast or slow method.
        
        Inputs:
            blocks: a 2D numpy array with blocks oriented along axis 1.
            usefast: a boolean value. If True use median removal.
                        If False remove linear trend.
    """
    if usefast:
        meds = np.median(blocks, axis=1)
        detrended = blocks - meds[:,np.newaxis]
    else:
        detrended = scipy.signal.detrend(blocks, type='linear', axis=1)
    return detrended


def calc_stats(blocks):
    """Calculate statistics for a series of blocks.
        
        The following gets rid of (hopefully) most of the 
        outlying values (i.e. power dropouts and single pulses)
        If you throw out 5% (2.5% at bottom and 2.5% at top)
        of random gaussian deviates, the measured stdev is ~0.871
        of the true stdev.  Thus the 1.0/0.871=1.148 correction below.
        The following is roughly .std() if the median is removed.
        
        Input:
            blocks: a 2D numpy array with blocks oriented along axis 1.

        Outputs:
            stds: Standard deviation for each block (after removing 
                    2.5% of smallest values and 2.5% of largest values
                    from each block).
            std_stds: Standard deviation of stds (after removing outliers)
            median_stds: Pseudo-median of stds (after removing outliers)
    """
    numblocks, blocklen = blocks.shape
    stds = np.sqrt((np.sort(blocks, axis=1)[:,blocklen/40:-blocklen/40]**2.0).sum(axis=1).astype(np.float64) / \
                    (0.95*blocklen))
    stds *= 1.148

    sort_stds = stds.copy()
    sort_stds.sort()

    # identify the differences with the largest values (this
    # will split off the chunks with very low and very high stds)
    locut = (sort_stds[1:numblocks/2+1] -
             sort_stds[:numblocks/2]).argmax() + 1
    hicut = (sort_stds[numblocks/2+1:] -
             sort_stds[numblocks/2:-1]).argmax() + numblocks/2 - 2
    std_stds = np.std(sort_stds[locut:hicut])
    median_stds = sort_stds[(locut+hicut)/2]

    return stds, std_stds, median_stds


class Chunk(object):
    """A Chunk object to represent part of a timeseries
        and efficiently search for single pulses.
    """
    def __init__(self, data, lobin, info):
        self.data = np.ma.masked_array(data, fill_value=0.0)
        self.lobin = lobin # Index of low bin of the Chunk with
                           # respect to the whole timeseries,
                           # excluding the overlap region.
        self.info = info

        self.fftd_data = None

    def smoothed(self, factor):
        """Return a smoothed numpy array of the chunk.
            Overlap is removed.

            Input:
                factor: number of bins to smooth by.
        """
        if factor == 1:
            # No smoothing
            smoothdata = self.data.filled()
        elif factor > 1:
            if self.fftd_data is None:
                self.fftd_data = presto.rfft(self.data.filled(), -1)
            kernel = fftd_kern_cache[factor]
            smoothdata = utils.fft_convolve(self.fftd_data, kernel)
        else:
            raise ValueError("Smoothing factor must be larger than 0.")
        return smoothdata[OVERLAP:-OVERLAP]

    def search(self, threshold, factor=1):
        """Search for significant events in the Chunk.
            
            Inputs:
                factor: number of bins to smooth by.
                threshold: sigma value above which events are
                        considered significant.

            Output:
                candlist: a CandidateList object containing
                            the candidates found.
        """
        candlist = cands.CandidateList()
        if self.data[OVERLAP:-OVERLAP].count(): # Search only if there are unmasked blocks
            smoothdata = self.smoothed(factor)
            hibins = np.flatnonzero(smoothdata>threshold)
            for localbin in hibins:
                bin = localbin+self.lobin
                from_badblock = self.data[OVERLAP:-OVERLAP].mask[localbin]
                c = cands.Candidate(self.info.DM, smoothdata[localbin], 
                                    self.info.dt*bin, bin, factor, from_badblock)
                candlist.append(c)
        return candlist
