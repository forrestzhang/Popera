class Hotspot:

    """
        Hotspots

    """

    def __init__(self, start, end, chromosome, hotspotid, peaks=list(), reads=0):

        region = chromosome+":"+str(start)+"-"+str(end)

        self.start = start

        self.end = end

        self.chromosome = chromosome

        self.hotspotid = hotspotid

        self.region = region

        self.reads = reads

        self.peaks = peaks


    def addpeaks(self, peak):

        self.peaks.append(peak)