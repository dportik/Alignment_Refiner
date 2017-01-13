# Alignment_Refiner

usage: python Alignment_Refinement.py [directory with phylip files] [Full Path to 'Master_Alignment_Assessment.txt']
example: python Alignment_Refinement.py dan/exons/phylip_files dan/exons/phylip_files/Alignment_Assessment/Master_Alignment_Assessment.txt


This script uses the 'Master_Alignment_Assessment.txt' generated from the 'alignment_assessment.py' script
and the same set of phylip alignment files.

Allows refinement of phylip alignment files via four options:

a. Reduce alignment set using a minimum base pair cutoff (ex. <100 bp).
- Will essentially copy and paste all alignment files passing threshold to 'Output_Refinement/' directory.

b. Remove sequences from all alignments using a missing data threshold (ex. only allow sequences with >50% bp present).
- Trims out all sequences above missing data threshold and creates new phylip files and 
   a log file to 'Output_Refinement/' directory. The log file records missing data levels for 
   the sequences that are removed per alignment. You should look at this, sometimes sequences
   have extremely high missing data levels (99%). 

c. For alignments above a missing data threshold (ex. >15% missing data), trim the alignments again.
- Uses trimal to do this with -gappyout -keepseqs options on. Creates new phylip files and
	a log file in 'Output_Refinement/' directory for alignments above missing data threshold. 
	Will ignore files below threshold and simply copy and paste them to the 'Output_Refinement/' directory.
	The log file records which alignments are trimmed and their missing data level. Check
	input alignment file and output alignment file for performance.

d. Perform all the above steps, (a, then c, then b).
- For alignments that are above the minimum alignment base pair cutoff, will trim alignments if above the
	overall missing data threshold (and ignore those below it), then perform sequence
	trimming across both re-trimmed and non-trimmed alignments. All output files and log files
	are created in 'Output_Refinement/' directory. You should look at the log files. 
	
I highly recommend re-running the alignment_assessment.py script and look through the new
Master_Alignment_Assessment.txt file to see how well this refinement performed. You can change the parameters
involved with option (d) to clean your data set, and I don't have default recommendations. 

written for Python 2.7.3

DEPENDENCIES:

numpy - Numerical Python

trimal v 1.4 (needs to be in path to call as 'trimal')


------------------------

Dan Portik

daniel.portik@berkeley.edu

August 2015

If you use this script, please cite:

Portik, D.M., Smith, L.L., and K. Bi. 2016. An evaluation of transcriptome-based exon capture for frog phylogenomics across multiple scales of divergence (Class: Amphibia, Order: Anura). Molecular Ecology Resources 16: 1069–1083.
