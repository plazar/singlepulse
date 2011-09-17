#!/usr/bin/env python

"""
A script to search for single pulses in *.dat files.

Patrick Lazarus, Jan. 3, 2011
"""

import glob
import sys
import optparse

import numpy as np

import infodata
import singlepulse

def main():
    # Grab *.dat files
    datfns = args
    for g in options.globs:
        datfns += glob.glob(g)

    if not len(datfns):
        raise singlepulse.errors.SinglePulseError("No *.dat found! Exiting...\n")
    print "Number of *.dat files found: %d" % len(datfns)

    # search *.dat files
    for ii, datfn in enumerate(datfns):
        sys.stdout.write("Working on file %d of %d.\n" % (ii+1, len(datfns)))
        sys.stdout.flush()

        if options.maxwidth is not None:
            # Convert maxwidth to bins
            index = datfn.rfind(".dat")
            if index == -1:
                raise errors.SinglePulseError("File is not a *.dat file (%s)!" % datfn)
            filenmbase = datfn[:index]
            info = infodata.infodata(filenmbase+".inf")
            options.maxwidth_bins = options.maxwidth/info.dt
        
        downfacts = [d for d in options.downfacts if d <= options.maxwidth_bins]

        candlist = singlepulse.search.find_single_pulses(datfn, downfacts, \
                            options.threshold, \
                            search_badblocks=options.search_badblocks)
        # Write out info for singlepulse files found
        print "Writing *.singlepulse file..."
        candlist.write_singlepulses(datfn[:-4]+".singlepulse", \
                            bb_col=options.search_badblocks)
    

if __name__ == '__main__':
    parser = optparse.OptionParser()
    parser.add_option('-t', '--threshold', type='float', dest='threshold',
                        help="Threshold for detecting single pulse events. " \
                             "(Default: 5.0).", \
                        default=5.0)
    parser.add_option('-f', '--outfn', dest='outfn', default="tmp.singlepulse", \
                        help="File name to write out singlepulses to. " \
                             "(Default: tmp.singlepulse).")
    parser.add_option('-g', '--glob', action='append', dest='globs', \
                        help="Glob expression of files to read " \
                             "single pulse events from. Multiple " \
                             "-g/--glob expressions can be provided.", \
                        default=[])
    parser.add_option('--search-badblocks', dest='search_badblocks', \
                        action='store_true', default=False, \
                        help="If set blocks with excessively low or high sigmas " \
                                "are _not_ set to 0 before searching. " \
                                "(Default: null out outlier blocks).")
    parser.add_option('-m', '--maxwidth', dest='maxwidth', type='float', \
                        help="Max downsampling in seconds. Number of bins " \
                             "this corresponds to will be determined for " \
                             "each input *.dat file independently. " \
                             "(Default: Use downfacts up to 30 bins wide). " \
                             "NOTE: -m/--maxwidth takes precedence over " \
                             "-b/--maxwidth-bins.", \
                        default=None)
    parser.add_option('-b', '--maxwidth-bins', dest='maxwidth_bins', type='int', \
                        help="Max downsampling in bins. All input *.dat files " \
                             "will use the same downsampling widths regardless " \
                             "of their sample rates. (Default: Use downfacts up " \
                             "to 30 bins wide). " \
                             "NOTE: -m/--maxwidth takes precedence over " \
                             "-b/--maxwidth-bins.", \
                        default=30)
    parser.add_option('-d', '--downfacts', dest='downfacts', action='callback', \
                        callback=singlepulse.utils.parselist_callback, \
                        callback_kwargs={'cast':int}, type='string', \
                        help="A comma-separated list of downfacts to use. " \
                             "(Default: %s) " \
                             "NOTE: This option needs to be combined with " \
                             "-m/--maxwidth or -b/--maxwidth-bins to search " \
                             "using downsampling factors larger than 30 bins. " \
                             "Also, downfacts > 192 should not be used. " % \
                             singlepulse.search.DEFAULT_DOWNFACTS, \
                        default=singlepulse.search.DEFAULT_DOWNFACTS)
    options, args = parser.parse_args()
    main()

