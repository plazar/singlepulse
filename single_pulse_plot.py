#!/usr/bin/env python

"""
A script to read single pulses and plot them.

Patrick Lazarus, Dec. 30, 2010
"""

import glob
import sys
import optparse
import os.path

import numpy as np
import matplotlib.pyplot as plt

import infodata
import singlepulse

def main():
    # Grab *.singlepulse files
    spfns = args
    for g in options.globs:
        spfns += glob.glob(g)

    if not len(spfns):
        raise singlepulse.errors.SinglePulseError("No *.singlepulses found! " \
                                                    "Exiting...\n")
    print "Number of *.singlepulse files found: %d" % len(spfns)

    # Read *.single pulse files
    candlist = singlepulse.cands.CandidateList()
    allDMs = []
    for ii, spfn in enumerate(spfns):
        candlist += singlepulse.cands.read_singlepulses(spfn)
        inffn = os.path.splitext(spfn)[0]+".inf"
        info = infodata.infodata(inffn)
        allDMs.append(info.DM)
        sys.stdout.write("\rReading *.singlepulse files (%6.2f%%)" % \
                            (float(ii+1)/len(spfns)*100))
        sys.stdout.flush()
    
    sys.stdout.write("\n")
    sys.stdout.flush()
    
    allDMs = np.unique(allDMs)

    candlist.trim(dmlim=(options.lodm, options.hidm), \
                    timelim=(options.start, options.end), \
                    minsigma=options.threshold)

    if not len(candlist):
        raise singlepulse.errors.SinglePulseError("No candidates to plot!")
    print "Plotting..."
    fig = singlepulse.plot.plot(candlist, \
                            timelim=(options.start, options.end), \
                            minsigma=options.threshold, allDMs=allDMs)
    if options.savefn:
        print "Saving to file (%s)..." % options.savefn
        plt.savefig(options.savefn, papertype='letter', orientation='landscape')
    if options.interactive:
        plt.figtext(0.98, 0.02, "Press 'Q' to quit", ha="right", size='x-small')
        def quit(event):
            if event.key in ('q', 'Q'):
                print "Quitting..."
                plt.close(fig)
        fig.canvas.mpl_connect('key_press_event', quit)
        plt.show()
    

if __name__ == '__main__':
    parser = optparse.OptionParser()
    parser.add_option('-t', '--threshold', type='float', dest='threshold',
                        help="Only plot events with sigma above threshold. " \
                             "(Default: 5.0).", \
                        default=5.0)
    parser.add_option('-s', '--start', type='float', dest='start',
                        help="Only plot events occuring after this time. " \
                             "(in s; Default: 0 s).", \
                        default=0.0)
    parser.add_option('-e', '--end', type='float', dest='end',
                        help="Only plot events occuring before this time. " \
                             "(in s; Default: {end of observation}).", \
                        default=np.inf)
    parser.add_option('-l', '--lodm', type='float', dest='lodm',
                        help="Only plot events with DM larger than " \
                             "this DM. (Default: 0 pc cm-3).", \
                        default=0.0)
    parser.add_option('-d', '--hidm', type='float', dest='hidm',
                        help="Only plot events with DM smaller than " \
                             "this DM. (Default: 0 pc cm-3).", \
                        default=np.inf)
    parser.add_option('-n', '--noninteractive', action='store_false', \
                        dest='interactive', default=True, \
                        help="Do not interactively display the plot. " \
                             "(Default: show plot interactively).")
    parser.add_option('-f', '--savefn', dest='savefn', default=None, \
                        help="File name to save plot as. (Default: " \
                             "don't save plot).")
    parser.add_option('-g', '--glob', action='append', dest='globs', \
                        help="Glob expression of files to read " \
                             "single pulse events from. Multiple " \
                             "-g/--glob expressions can be provided.", \
                        default=[])
    options, args = parser.parse_args()
    main()

