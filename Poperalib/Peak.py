class Peak:
    """
        Peaks
    
    """    
   
    def __init__(self, start, end, chromosome, parent, pvalue, peakpoint, peakid, score):

        self.start = start

        self.end = end

        self.chromosome = chromosome

        self.parent = parent

        self.pvalue = float(pvalue)

        self.peakpoint = peakpoint

        self.peakid = peakid

        self.score = score


