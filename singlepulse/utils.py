import numpy as np

import presto


class KernelCache(dict):
    """A cache for convolution kernels. If a requested kernel
        doesn't exist in the cache it will be calculated,
        stored and returned.
    """
    def __init__(self, kernel_factory, *args, **kwargs):
        """KernelCache object constructor.

            Inputs:
                kernel_factory: a function that constructs convolution
                                kernels. It should take 1 argument,
                                which is also used as the kernel's key.
                Other arguments are the usual arguments to dict(...)
        """
        super(KernelCache, self).__init__(*args, **kwargs)
        self.kernel_factory = kernel_factory

    def __missing__(self, key):
        """Call kernel_factory with key provided and its only argument.
            Insert the result into the dict and return the result.
        """
        return self.setdefault(key, self.kernel_factory(key))
       
    def set_kernel_factory(self, kf):
        self.kernel_factory = kf


def make_fftd_kern(downfact, fftlen):
    kern = np.zeros(fftlen, dtype=np.float32)
    # These offsets produce kernels that give results
    # equal to scipy.signal.convolve
    if downfact % 2:  # Odd number
        kern[:downfact/2+1] += 1.0
        kern[-(downfact/2):] += 1.0
    else:             # Even number
        kern[:downfact/2+1] += 1.0
        if (downfact > 2):
            kern[-(downfact/2-1):] += 1.0
    # The following normalization preserves the
    # RMS=1 characteristic of the data
    return presto.rfft(kern / np.sqrt(downfact), -1)


def make_kern(downfact):
    return np.ones(downfact, dtype=Num.float32) / np.sqrt(downfact)
    

def fft_convolve(fftd_data, fftd_kern):
    """Perform a convolution with the complex floating point vectors
            'fftd_data' and 'fftd_kern'.
    """
    # Note:  The initial FFTs should be done like:
    # fftd_kern = presto.rfft(kernel, -1)
    # fftd_data = presto.rfft(data, -1)
    prod = np.multiply(fftd_data, fftd_kern)
    prod.real[0] = fftd_kern.real[0] * fftd_data.real[0]
    prod.imag[0] = fftd_kern.imag[0] * fftd_data.imag[0]
    return presto.rfft(prod, 1).astype(np.float32)


def parselist_callback(option, opt, value, parser, cast=(lambda x:x)):
    """A function to parse comma-separated lists privded
        to the optparse.

        Inputs:
            cast: an optional function to cast elements of the list.
    """
    print "Parsing:", value
    parsedlist = [cast(el) for el in value.split(',')]
        
    setattr(parser.values, option.dest, parsedlist)
