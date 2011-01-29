import os.path
import types

import numpy as np

import infodata
import errors

class Candidate:
    """Single pulse candidate object.
    """
    def __init__(self, DM, sigma, time, bin, downfact):
        self.DM = DM
        self.sigma = sigma
        self.time = time
        self.bin = bin
        self.downfact = downfact
        self.duration = self.time/self.bin*self.downfact

    def __str__(self):
        return "%7.2f %7.2f %13.6f %10d     %3d\n"%\
               (self.DM, self.sigma, self.time, self.bin, self.downfact)

    def __cmp__(self, other):
	# Sort by bin by default
        return cmp(self.bin, other.bin)


class CandidateList(list):
    def __init__(self, *args):
        super(CandidateList, self).__init__(*args)
        self.infos = {} # Dictionary of infodata objects keyed by DM

    def get(self, key):
        return np.array([getattr(c, key) for c in self])

    def __getattr__(self, key):
        return self.get(key)

    def __add__(self, other):
        result = CandidateList(super(CandidateList, self).__add__(other))
        result.infos = self.infos
        for info in other.infos.values():
            result.add_info(info)
        return result

    def __iadd__(self, other):
        return self + other

    def add_info(self, info):
        self.infos[info.DM] = info

    def trim(self, dmlim=(0,np.inf), timelim=(0,np.inf), minsigma=0):
        """Remove candidates that are outside given DM or time
            ranges or below a sigma threshold.

            Note: trimming is non-reversible.
        """
        for jj in np.arange(len(self))[::-1]:
            c = self[jj]
            if c.sigma<minsigma or \
                    c.time<timelim[0] or c.time>timelim[1] or \
                    c.DM<dmlim[0] or c.DM>dmlim[1]:
                self.pop(jj)
        for dm in self.infos.keys():
            if dm<dmlim[0] or dm>dmlim[1]:
                self.infos.pop(dm)

    def prune_related1(self):
        """Remove candidates that are close to other candidates
            but less significant. Treats candidates from different
            downfacts independently.

            Note: Pruning is non-reversible.
        """
        toremove = set()
        downfacts = np.unique(self.downfact)
        for downfact in downfacts:
            indices = np.flatnonzero((self.downfact==downfact))
            hibins = self.bin[indices]
            hivals = self.sigma[indices]

            for ii in range(0, len(hibins)-1):
                if ii in toremove:  continue
                xbin, xsigma = hibins[ii], hivals[ii]
                for jj in range(ii+1, len(hibins)):
                    ybin, ysigma = hibins[jj], hivals[jj]
                    if (abs(ybin-xbin) > downfact/2):
                        break
                    else:
                        if jj in toremove:
                            continue
                        if (xsigma > ysigma):
                            toremove.add(indices[jj])
                        else:
                            toremove.add(indices[ii])
        # Now zap them starting from the end
        toremove = sorted(toremove, reverse=True)
        for ii in toremove:
            self.pop(ii)

    def prune_related2(self, maxdownfact):
        """Remove candidates that are close to other candidates
            but less significant.  This one works on the candidate 
            instances and looks at the different downfacts of the
            the different candidates.
 
            Note: Pruning is non-reversible
        """
        self.sort()
        toremove = set()
        for ii in range(0, len(self)-1):
            if ii in toremove:  continue
            xx = self[ii]
            xbin, xsigma = xx.bin, xx.sigma
            for jj in range(ii+1, len(self)):
                yy = self[jj]
                ybin, ysigma = yy.bin, yy.sigma
                if (abs(ybin-xbin) > maxdownfact/2):
                    break
                else:
                    if jj in toremove:
                        continue
                    prox = max([xx.downfact/2, yy.downfact/2, 1])
                    if (abs(ybin-xbin) <= prox):
                        if (xsigma > ysigma):
                            toremove.add(jj)
                        else:
                            toremove.add(ii)
        # Now zap them starting from the end
        toremove = sorted(toremove, reverse=True)
        for ii in toremove:
            self.pop(ii)
    
    def prune_border_cases(self):
        """Ignore those that are locate in a half-width
            of the boundary between data and padding.
        """
        info0 = self.infos[min(self.infos.keys())]
        if info0.breaks:
            offregions = zip([x[1] for x in info0.onoff[:-1]],
                             [x[0] for x in info0.onoff[1:]])
            toremove = set()
            for ii in range(len(self))[::-1]:
                cand = self[ii]
                loside = cand.bin-cand.downfact/2
                hiside = cand.bin+cand.downfact/2
                if hiside < offregions[0][0]: break
                for off, on in offregions:
                    if (hiside > off and loside < on):
                        toremove.add(ii)
            # Now zap them starting from the end
            toremove = sorted(toremove, reverse=True)
            for ii in toremove:
                self.pop(ii)

    def write_singlepulses(self, outfn):
        """Write single pulses in CandidateList to 'outfn'.
        
            Inputs:
                outfn: file to write to. Either an already opened 
                        file object or a filename.
        """
        if type(outfn) == types.StringType:
            outfile = open(outfn, 'w')
            toclose = True
        else:
            outfile = outfn
            toclose = False
 
        if len(self):
            outfile.write("# DM      Sigma      Time (s)     Sample    Downfact\n")
            for cand in self:
                outfile.write(str(cand))
        
        if toclose:
            # File was opened here, so close it.
            outfile.close()


def sorted_candlist(self, *args, **kwargs):
    return CandidateList(sorted(self, *args, **kwargs))


def cmp_cand(key):
    """Return a candidate comparer using the attribute key.
    """
    return lambda cand1, cand2: cmp(getattr(cand1, key), getattr(cand2, key))


def read_singlepulses(infile):
    """Read single pulses from *.singlepulse file.
    """
    index = infile.rfind(".singlepulse")
    if index == -1:
        raise errors.SinglePulseError("File is not a *.singlepulse file (%s)!" % \
                                        infile)
    filenmbase = infile[:index]
    
    info = infodata.infodata(filenmbase+".inf")
    if os.path.getsize(infile):
        canddata = np.loadtxt(infile, dtype=[('DM','f8'), ('sigma','f8'), ('time','f8'), ('bin','i8'), ('downfact','i8')])
        canddata = np.atleast_1d(canddata)
        candlist = CandidateList([Candidate(*data) for data in canddata])
    else:
        candlist = CandidateList()
    candlist.add_info(info)
    return candlist

