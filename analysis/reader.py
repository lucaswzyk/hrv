# simple reader class for .txt files containing float RR intervals
#
# author: lucaswzyk

from random import shuffle
from numpy import array


# rrs         RR series
# rrs_cumul   accumulated RR values from rrs
# num_rr      length of RR series
class RRData:
    def __init__(self, settings):
        self.rrs, self.rrs_cumul, self.num_rr = self.prep(settings)

    # reads RR intervals from file and returns them as rrs
    @staticmethod
    def get_rrs(filename, min_beat, max_beat):
        rrs = []
        with open(filename) as f:
            i = min_beat
            for line in f:
                # ms to s
                rrs.append(float(line))
                i += 1
                if i > max_beat:
                    break
        f.close()

        # uncomment for verifying algorithm validity
        # shuffle(rrs)

        return rrs

    # calls get_rrs, set num_rr, accumulates rrs_cumul
    def prep(self, settings):
        rrs = self.get_rrs(settings.file, settings.min_beat, settings.max_beat)
        num_rr = len(rrs)

        s = 0
        rrs_cumul = []
        for j in rrs:
            rrs_cumul.append(s + j)
            s += j

        return array(rrs), array(rrs_cumul), num_rr
