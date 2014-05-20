Popera
======

DNase I hypersensitive sites identification


Options:

  --version             show program's version number and exit
  
  -h, --help            show this help message and exit.
  
  -d DATAFILES, --data=DATAFILES
                        data file, should be sorted bam format,
                        example=DH_sample1.bam,DH_sample2.bam
                        
  -n SAMPLENAMES, --name=SAMPLENAMES
                        NH sample name default=DH_sample1,DH_sample2
                        
  -b BW, --bandwidth=BW
                        kernel smooth band width, should >1, default=50
                        
  -t THRESHOLD, --threshold=THRESHOLD
                        Hot spots threshold, default=4.0
                        
  -l MINLENGTH, --minlength=MINLENGTH
                        minimum length of hot spots, default=5
                        
  -p PVALUE, --pavlue=PVALUE
                        p-value cutoff for peak identification, default=0.05
                        
  -i INITIAL, --initial=INITIAL
                        Peak's initial length, >1 and <minlength, default=5
                        
  --threads=NTHREADS    threads number or cpu number, default=4
  
  -w, --wig             whether out put wiggle file, default=False
  
  -f, --fdr             using FDR instead p-value
  
  -x EXCLUDECHR, --excludechr=EXCLUDECHR
                        Don't count those DHs, example='-x ChrM,ChrC'
                        
  -g, --gff             whether out put gff file, default=False
  
  --pe                  paired-end reads or single-end reads, default=False
                        (single end)
