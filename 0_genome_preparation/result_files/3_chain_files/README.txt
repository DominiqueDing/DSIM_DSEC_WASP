This directory contains the Drosophila sechellia 28 October 2005 assembly 
from the Broad Institute at MIT and Harvard.  The annotations are from 
UCSC and collaborators worldwide.

If you plan to download a large file or multiple files from 
this directory, we recommend that you use ftp rather than 
downloading the files via our website. To do so, ftp to 
hgdownload.cse.ucsc.edu, then go to the directory 
goldenPath/droSec1/bigZips. To download multiple files, use 
the "mget" command:

    mget <filename1> <filename2> ...
    - or -
    mget -a (to download all the files in the directory) 
 
------------------------------------------------------------
md5sum.txt - MD5 checksum of these files to verify correct transmission.

droSec1.2bit - contains the complete D. sechellia/droSec1 genome sequence
    in the 2bit file format.  Repeats from RepeatMasker and Tandem Repeats
    Finder (with period of 12 or less) are shown in lower case; non-repeating
    sequence is shown in upper case.  The utility program, twoBitToFa (available
    from the kent src tree), can be used to extract .fa file(s) from
    this file.  A pre-compiled version of the command line tool can be
    found at:
        http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/
    See also:
        http://genome.ucsc.edu/admin/git.html
	http://genome.ucsc.edu/admin/jk-install.html

scaffoldFa.gz - The working draft sequence in one FASTA record per scaffold.
    Repeats from RepeatMasker and Tandem Repeats Finder (with period
    of 12 or less) are in lower case while non-repeating sequence is
    in upper case.  

scaffoldFaMasked.gz - The working draft sequence in one FASTA record per
    scaffold. Repeats are masked by capital N's and non-repeating
    sequence is shown in upper case.

scaffoldOut.gz - RepeatMasker .out file for scaffolds.  These were 
    created with RepeatMasker at the -s sensitive setting.

scaffoldTrf.gz - Tandem Repeats Finder locations, filtered to keep repeats 
    with period less than or equal to 12, translated into one .bed file. 

xenoMrna.fa.gz - GenBank mRNAs from species other than that of 
    the genome. This sequence data is updated once a week via automatic 
    GenBank updates.
