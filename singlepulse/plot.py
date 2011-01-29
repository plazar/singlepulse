import os
import types

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

import cands

def plot(candlist, savefn=None, interactive=False, \
                timelim=None, minsigma=None):
    """Plot the candidates in 'candlist'.

        Inputs:
            candlist: a list of singlepulse Candidate objects.
                        (as defined in singlepulse/cands.py)
            savefn: file to save plot as. If 'savefn' is None
                    the plot is not saved. (Default: Do not
                    save the plot).
            interactive: show the plot interactively. (Default:
                    False, do not show plot).
            timelim: limits on time. (Default: None)
            minsigma: minimum sigma. (Default: None)
    """
    newplot(candlist, timelim, minsigma)

    if type(savefn) == types.StringType:
        print "Saving to file (%s)..." % savefn
        plt.savefig(savefn, papertype='letter', orientation='landscape')
    elif savefn is not None:
        raise ValueError("'savefn' must be a string.")
    if interactive:
        plt.figtext(0.98, 0.02, "Press 'Q' to quit", ha="right", size='x-small')
        def quit(event):
            if event.key in ('q', 'Q'):
                print "Quitting..."
                plt.close(event.canvas.figure)
        plt.gcf().canvas.mpl_connect('key_press_event', quit)
        plt.show()

def oldplot(candlist):
    """Plot the candidates in 'candlist' in the old style.

        Inputs:
            candlist: a list of singlepulse Candidate objects.
                        (as defined in singlepulse/cands.py)
    """
    raise NotImplementedError("Old-style plot is not fully implemented.")
    fig = plt.figure(figsize=(8.5, 11))
    sortedcands = cands.sorted_candlist(candlist, \
                            cmp=cands.cmp_cand('sigma'))
    sigmas = sortedcands.sigma
    times = sortedcands.time
    DMs = sortedcands.DM
    downfacts = sortedcands.downfact
    dmtime_ax = plt.axes((0.1, 0.1, 0.85, 0.5))
    plt.scatter(times, DMs, c=downfacts, \
                        s=np.clip((sigmas-5)**2+5, 5, 200), \
                        cmap=matplotlib.cm.GnBu)
    plt.xlabel("Time (s)", size="small")
    plt.ylabel(r"DM (pc cm$^{-3}$)", size="small")
    cb = plt.colorbar()
    cb.set_label("Downfactor (bins - CONVERT TO TIME!)", rotation=90, size="small")

    snrhist_ax = plt.axes((0.1, 0.65, 0.2, 0.25))
    plt.hist(sigmas, int(np.ptp(sigmas)+1), (sigmas.min(), sigmas.max()), \
                histtype='step', ec='k')
    plt.xlabel("Sigma")
    plt.ylabel("Number of pulses")

    dmhist_ax = plt.axes((0.35, 0.65, 0.2, 0.25))
    plt.hist(DMs, bins=sorted(candlist.infos.keys()), histtype='step', ec='k')

    snrdm_ax = plt.axes((0.6, 0.65, 0.2, 0.25))
    plt.plot(DMs, sigmas, 'k,')

def newplot(candlist, timelim=None, minsigma=None):
    """Plot the candidates in 'candlist' in the new style.

        Inputs:
            candlist: a list of singlepulse Candidate objects.
                        (as defined in singlepulse/cands.py)
            timelim: limits on time. (Default: None)
            minsigma: minimum sigma. (Default: None)
    """
    fig = plt.figure(figsize=(11, 8.5))
    sortedcands = cands.sorted_candlist(candlist, \
                            cmp=cands.cmp_cand('sigma'))
    sigmas = sortedcands.sigma
    times = sortedcands.time
    DMs = sortedcands.DM # DMs of Candidates
    allDMs = sorted(candlist.infos.keys()) # List of DMs searched (including 
                                           # those where no candidates were found).
    duration = sortedcands.duration
    dmtime_ax = plt.axes((0.1, 0.1, 0.5, 0.6))
    plt.scatter(times, DMs, c=duration*1000, \
                        s=(np.clip(sigmas, 0, 20)-np.min(sigmas)+1)**2, \
                        cmap=matplotlib.cm.GnBu, alpha=0.5)
    plt.xlabel("Time (s)", size="small")
    plt.ylabel(r"DM (pc cm$^{-3}$)", size="small")
    plt.setp(dmtime_ax.xaxis.get_ticklabels(), size='x-small')
    plt.setp(dmtime_ax.yaxis.get_ticklabels(), size='x-small')
    cb_ax = plt.axes((0.85, 0.705, 0.02, 0.24), frameon=False)
    cb = plt.colorbar(cax=cb_ax)
    cb.set_label("Pulse duration (ms)", rotation=90, size="small")
    plt.setp(cb_ax.yaxis.get_ticklabels(), size='x-small')
    
    snrhist_ax = plt.axes((0.6, 0.7, 0.2, 0.25))
    num_v_sigma, edges = np.histogram(sigmas, \
                            bins=np.arange(int(sigmas.min()), int(sigmas.max())+2))
    
    num_v_sigma[num_v_sigma<0.001] = 0.001
    plt.plot(edges.repeat(2), np.concatenate(([0.001], num_v_sigma.repeat(2), [0.001])), 'k-')
    plt.ylabel("Number of pulses", size='small')
    plt.yscale('log')
    plt.ylim(0.5, 2*num_v_sigma.max())
    plt.setp(snrhist_ax.xaxis.get_ticklabels(), visible=False)
    plt.setp(snrhist_ax.yaxis.get_ticklabels(), size='x-small')

    dmhist_ax = plt.axes((0.8, 0.1, 0.15, 0.6), sharey=dmtime_ax)
    vals, edges = np.histogram(DMs, bins=allDMs)
    plt.plot(np.concatenate(([0], vals.repeat(2), [0])), edges.repeat(2), 'k-')
    plt.xlabel("Number of pulses", size="small")
    plt.setp(dmhist_ax.yaxis.get_ticklabels(), visible=False)
    plt.setp(dmhist_ax.xaxis.get_ticklabels(), size='x-small')

    snrdm_ax = plt.axes((0.6, 0.1, 0.2, 0.6), sharex=snrhist_ax, \
                            sharey=dmtime_ax)
    plt.plot(sigmas, DMs, 'k,')
    plt.xlabel("Sigma", size="small")
    plt.setp(snrdm_ax.yaxis.get_ticklabels(), visible=False)
    plt.setp(snrdm_ax.xaxis.get_ticklabels(), size='x-small')
    
    dmtime_ax.set_ylim(allDMs[0], allDMs[-1]) # allDMs is sorted
    info0 = candlist.infos[allDMs[0]] # infodata object of lowest DM timeseries
    if timelim is not None:
        dmtime_ax.set_xlim(max(0, timelim[0]), min(info0.N*info0.dt, timelim[-1]))
    else:
        dmtime_ax.set_xlim(0, info0.N*info0.dt)
    if minsigma is not None:
        snrhist_ax.set_xlim(max(int(sigmas.min()), minsigma), int(sigmas.max())+2)
    else:
        snrhist_ax.set_xlim(int(sigmas.min()), int(sigmas.max())+2)
    
    dmtime_ax.format_coord = lambda x,y: "Time (s)=%.2f, DM (pc cm-3)=%.2f" % (x,y)
    snrhist_ax.format_coord = lambda x,y: "Sigma=%.2f, N=%d" % (x, np.round(y))
    dmhist_ax.format_coord = lambda x, y: "N=%d, DM (pc cm-3)=%.2f" % (np.round(x), y)
    snrdm_ax.format_coord = lambda x, y: "Sigma=%.2f, DM (pc cm-3)=%.2f" % (x,y)

    # Print text area
    plt.figtext(0.05, 0.96, "Single pulse results for '%s'" % \
                    os.path.split(info0.basenm)[-1].split("_DM")[0], size='large')
    plt.figtext(0.05, 0.90, "Source: %s" % info0.object)
    plt.figtext(0.05, 0.87, "Telescope: %s" % info0.telescope)
    plt.figtext(0.05, 0.84, "Instrument: %s" % info0.instrument)
    plt.figtext(0.05, 0.79, "RA (J2000): %s" % info0.RA)
    plt.figtext(0.05, 0.76, "Dec. (J2000): %s" % info0.DEC)
    if info0.bary:
        plt.figtext(0.05, 0.73, r"MJD$_{\mathrm{\mathsf{bary}}}$: %.12f" % info0.epoch)
    else:
        plt.figtext(0.05, 0.73, r"MJD$_{\mathrm{\mathsf{topo}}}$: %.12f" % info0.epoch)
    plt.figtext(0.3, 0.90, "N samples: %d" % info0.N)
    plt.figtext(0.3, 0.87, r"Sampling time: %.2f $\mu$s" % (info0.dt*1e6))
    plt.figtext(0.3, 0.84, r"Freq$_{\mathrm{\mathsf{ctr}}}$: %.1f MHz" % \
                    ((info0.numchan/2-0.5)*info0.chan_width+info0.lofreq))
    plt.figtext(0.3, 0.79, "Number of pulses: %d" % len(candlist))
    plt.figtext(0.3, 0.76, "Number of DMs: %d" % len(allDMs))
    if minsigma is not None:
        plt.figtext(0.3, 0.73, "Sigma threshold: %g" % minsigma)
    else:
        plt.figtext(0.3, 0.73, "Simga threshold: 0")
