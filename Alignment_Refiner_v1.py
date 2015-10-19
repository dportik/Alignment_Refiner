import sys
import os
import subprocess as sp
import shutil
import numpy as np

'''
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


##############
DEPENDENCIES:
numpy - Numerical Python
trimal v 1.4 (needs to be in path to call as 'trimal')
##############

------------------------
written for Python 2.7.3
Dan Portik
daniel.portik@berkeley.edu
August 2015
------------------------
'''
#move to phylip directory
aln_in_directory = sys.argv[1]
os.chdir(aln_in_directory)

#read in assessment file
file_name = sys.argv[2]
assess_fh = open(file_name, 'r')

#intialize empty lists
contig_name_l = []
taxa_no_l = []
seq_length_l = []
perc_isites_l = []
perc_mdata_l = []
alnseq90_l = []
alnseq90 = int(0)
alnseq70_l = []
alnseq70 = int(0)

#read in lines from assessment file
lines = assess_fh.readlines()
aln_no = int(0)
for line in lines[1:]:
    cleanline = line.strip()
    sline = cleanline.split('\t')
    
    #name variables
    contig_name = sline[0]
    taxa_no = sline[1]
    seq_length = sline[2]
    perc_isites = sline[5]
    perc_mdata = sline[12]
    seqsabove90 = int(sline[13])
    seqsabove70 = int(sline[14])
    
    #add line variable data to lists
    contig_name_l.append(contig_name)
    taxa_no_l.append(taxa_no)
    seq_length_l.append(seq_length)
    perc_isites_l.append(perc_isites)
    perc_mdata_l.append(perc_mdata)
    alnseq90_l.append(seqsabove90)
    alnseq70_l.append(seqsabove70)

    if seqsabove90 >= int(1):
        alnseq90 += 1
    if seqsabove70 >= int(1):
        alnseq70 += 1

    aln_no+=1

phylip_files_here = int(0)
for filetype in os.listdir('.'):
    if filetype.endswith(".phylip"):
        phylip_files_here += 1

print '\n',"-----------------------------------------------------------------------------------------------------------",'\n',"There are {0} entries in 'Master_Alignment_Assessment.txt', and {1} phylip files in this directory.".format(aln_no, phylip_files_here),'\n',"-----------------------------------------------------------------------------------------------------------",'\n'

#create function to convert lists into arrays and do very basic stats
def quickstats(x):
    x_array = np.asarray(x,dtype=np.float64)
    x_avg = np.average(x_array)
    x_avg = np.around(x_avg, decimals = 1)
    x_avg = str(x_avg)
    x_min = np.amin(x_array)
    x_min = str(x_min)
    x_max = np.amax(x_array)
    x_max = str(x_max)
    x_median = np.median(x_array)
    x_median = str(x_median)
    return x_avg, x_min, x_max, x_median

#execute functions and print stats to screen
taxa_avg, taxa_min, taxa_max, taxa_med = quickstats(taxa_no_l)
print "TAXON stats:"
print "------------------------------------------------------------------"
print "Average number of taxa = {}.".format(taxa_avg)
print "Median number of taxa = {}.".format(taxa_med)
print "Minimum number of taxa = {}.".format(taxa_min)
print "Maximum number of taxa = {}.".format(taxa_max), '\n', '\n'

seq_avg, seq_min, seq_max, seq_med = quickstats(seq_length_l)
print "ALIGNMENT LENGTH stats:"
print "------------------------------------------------------------------"
print "Average alignment length (bp) = {}.".format(seq_avg)
print "Median alignment length (bp) = {}.".format(seq_med)
print "Minimum alignment length (bp) = {}.".format(seq_min)
print "Maximum alignment length (bp) = {}.".format(seq_max), '\n', '\n'

isites_avg, isites_min, isites_max, isites_med = quickstats(perc_isites_l)
print "INFORMATIVE SITES stats:"
print "------------------------------------------------------------------"
print "Average percentage of informative sites = {}%.".format(isites_avg)
print "Median percentage of informative sites = {}%.".format(isites_med)
print "Minimum percentage of informative sites = {}%.".format(isites_min)
print "Maximum percentage of informative sites = {}%.".format(isites_max), '\n', '\n'

mdata_avg, mdata_min, mdata_max, mdata_med = quickstats(perc_mdata_l)
print "MISSING DATA stats by ALIGNMENT:"
print "------------------------------------------------------------------"
print "Average percentage of missing data in alignments = {}%.".format(mdata_avg)
print "Median percentage of missing data in alignments = {}%.".format(mdata_med)
print "Minimum percentage of missing data in alignments = {}%.".format(mdata_min)
print "Maximum percentage of missing data in alignments = {}%.".format(mdata_max), '\n', '\n'

s90_avg, s90_min, s90_max, s90_med = quickstats(alnseq90_l)
s70_avg, s70_min, s70_max, s70_med = quickstats(alnseq70_l)
print "MISSING DATA stats by SEQUENCES:"
print "------------------------------------------------------------------"
print "There are {} alignments with at least one sequence with >90% missing data.".format(alnseq90)
print "Average number of sequences with >90% missing data per alignment = {}.".format(s90_avg)
print "Median number of sequences with >90% missing data per alignment = {}.".format(s90_med)
print "Minimum number of sequences with >90% missing data per alignment = {}.".format(s90_min)
print "Maximum number of sequences with >90% missing data per alignment = {}.".format(s90_max), '\n'
print "There are {} alignments with at least one sequence with >70% missing data.".format(alnseq70)
print "Average number of sequences with >70% missing data per alignment = {}.".format(s70_avg)
print "Median number of sequences with >70% missing data per alignment = {}.".format(s70_med)
print "Minimum number of sequences with >70% missing data per alignment = {}.".format(s70_min)
print "Maximum number of sequences with >70% missing data per alignment = {}.".format(s70_max), '\n'



assess_fh.close()
 
print '\n',"###############################################################################################################", '\n', '\n',"Refinement Options:", '\n', '-------------------','\n'
print "a. Reduce alignment set using a minimum base pair cutoff (ex. <100 bp)."
print "b. Remove sequences from all alignments using a missing data threshold (ex. only allow sequences with >50% bp present)."
print "c. For alignments above a missing data threshold (ex. >15% missing data), trim the alignments again."
print "d. Perform all the above steps, (a, then c, then b).", '\n'
print "####################################################################################################################",'\n'

out_dir = "Output_Refinement"
if not os.path.exists(out_dir):
    os.mkdir(out_dir)
        
decision_analysis = (raw_input("Please select one of the above options (a, b, c, d), or 'q' to quit: "))





###############################################################
#Begin definition of analysis 'a'

if decision_analysis == "a":
    decision_bp = None
    while decision_bp is None:
        try:
            decision_bp = int(raw_input("Please enter a minimum base pair value for alignments (ex. 100): "))
        except ValueError:
            print "That wasn't a number."
    contig_list = []
    a_fh = open(file_name, 'r')
    lines = a_fh.readlines()
    for line in lines[1:]:
        cleanline = line.strip()
        sline = cleanline.split('\t')
        contig_name = sline[0]
        taxa_no = sline[1]
        seq_length = sline[2]
        perc_isites = sline[5]
        perc_mdata = sline[12]
        if int(seq_length) >= decision_bp:
            contig_list.append(contig_name)
    fileshere = int(0)
    for filetype1 in os.listdir('.'):
        if filetype1.endswith(".phylip"):
            fileshere += 1
    filesmoved = int(0)
    for contig in contig_list:
        for filetype in os.listdir('.'):
            if filetype.endswith(".phylip"):
                name = filetype.split('.')
                if name[0] == contig:
                    filesmoved += 1
                    new_name = out_dir+'/'+contig+".phylip"
                    shutil.copyfile(filetype, new_name)
                    print "{} Passed!".format(contig)
    filesleft = (fileshere - filesmoved)
    print '\n', "Of the {0} alignments, {1} alignments were shorter than {2} bp, and {3} alignments that were greater than {2} bp.".format(fileshere, filesleft, decision_bp, filesmoved), '\n'

###############################################################
#Begin definition of analysis 'b'

elif decision_analysis == "b":
    taxa_removed_l = []
    undisturbed = int(0)
    fh_log = open('analysis_b_log.txt', 'a')
    print '\n',"Please select the missing data threshold to remove sequences across alignments.", '\n'
    print "For example, entering 60 will remove a sequence if it is missing more than 60% of the base pairs in the alignment,"
    print "that is, if less than 40% of the bp are present."
    print "In an alignment 1000 bp long, this requires > 400 bases to be present, < 600 bases missing.", '\n'
    decision_indmiss = None
    while decision_indmiss is None:
        try:
            decision_indmiss = int(raw_input("Please enter a sequence missing data threshold (ex. 60): "))
        except ValueError:
            print "That wasn't a number."
    fh_log.write("Missing data threshold for sequences = {}.".format(decision_indmiss)+'\n')
    for filetype in os.listdir('.'):
        if filetype.endswith(".phylip"):
            print "-------------------------------------------------------------------------------"
            print "Refining {}:".format(filetype)
            fh_log.write("Removed the following sequences from {}:".format(filetype)+'\n')
            name = filetype.split('.')
            temp_name = "temp."+name[0]+".phylip"
            fh_temp = open(temp_name, 'a')
            with open(filetype, 'r') as fh_phy:
                inds = int(0)
                ind_kept = int(0)
                ind_removed = int(0)
                lines = fh_phy.readlines()
                for line in lines[1:]:
                    inds += 1
                    line = line.strip()
                    splitline = line.split()
                    taxon = splitline[0]
                    sequence = splitline[1]
                    total_bps = int(0)
                    gaps = int(0)
                    bps = int(0)
                    for base in sequence:
                        base = base.upper()
                        if base == 'A' or base =='T' or base =='G' or base =='C' or base =='N' or base =='R' or base =='Y' or base =='S' or base =='W' or base =='K' or base =='M' or base =='B' or base =='D' or base =='H' or base =='V':
                            bps+=1
                            total_bps+=1
                        elif base == '-':
                            gaps+=1
                            total_bps+=1

                    cutoff = int( float(( float(gaps) / float(total_bps)) * float(100)))
                    present = int(100) - cutoff
                    if cutoff >= int(decision_indmiss):
                        fh_log.write("Removing {0}".format(taxon)+'\t'+'\t'+"[{0} bp in alignment, {1} missing, {2} present]".format(total_bps, gaps, bps)+'\t'+"[{}% bp present]".format(present)+'\t'+"[{}% missing data]".format(cutoff)+'\n')
                        #print "Removing {0} which is missing {1}% of the data: of {2} bp, {3} missing, {4} present ({5}% present).".format(taxon, cutoff, total_bps, gaps, bps, present)
                        ind_removed += 1
                    elif cutoff < int(decision_indmiss):
                        ind_kept += 1
                        fh_temp.write(line+'\n')
                    taxa_removed_l.append(ind_removed)
                if ind_removed == int(0):
                    undisturbed += 1
            fh_log.write('\n'+'\n')
            print '\n', "Of the {0} sequences, {1} were kept and {2} were removed".format(inds, ind_kept, ind_removed)
            fh_temp.close()
            ind_kept = str(ind_kept)
            total_bps = str(total_bps)
            temp_name2 = "temp2."+name[0]+".phylip"
            fh_temp2 = open(temp_name2, 'a')
            fh_temp2.write(ind_kept+" "+total_bps+'\n')
            with open(temp_name, 'r') as fh_temp:
                lines = fh_temp.readlines()
                for line in lines:
                    line = line.strip()
                    fh_temp2.write(line+'\n')
            fh_temp2.close()
            fh_temp.close()
            print "-------------------------------------------------------------------------------", '\n'
    for filetype in os.listdir('.'):
        if filetype.startswith('temp2.'):
            namesplit = filetype.split('.')
            name = namesplit[1]
            move_name = out_dir+'/'+name+".phylip"
            shutil.copyfile(filetype, move_name)
            os.remove(filetype)
        if filetype.startswith('temp.'):
            os.remove(filetype)
        if filetype.startswith('analysis_b_log'):
            shutil.move(filetype, out_dir)
    fh_log.close()
    print "Moving all output files to output directory."
    rem_avg, rem_min, rem_max, rem_med = quickstats(taxa_removed_l)
    print "REMOVED SEQUENCES stats:"
    print "----------------------------------------------------------------------------------"
    print "Average number of taxa removed across alignments = {}.".format(rem_avg)
    print "Median number of taxa removed across alignments = {}.".format(rem_med)
    print "Minimum number of taxa removed across alignments = {}.".format(rem_min)
    print "Maximum number of taxa removed across alignments = {}.".format(rem_max)
    print "There were {0} alignments with no sequences removed (of {1}).".format(undisturbed, phylip_files_here), '\n', '\n'

###############################################################
#Begin definition of analysis 'c'
   
elif decision_analysis == "c":
    log_name = "analysis_c_log.txt"
    fh_logc = open(log_name, 'a')
    
    print '\n',"Please select the missing data threshold of an entire alignment to select it for trimming."
    print "For example, entering 30 will cause all alignments with more than 30% overall missing data to be re-trimmed."
    print "If you want to re-trim all alignments with more the more stringent criteria here, just enter 1", '\n'
    decision_alnmiss = None
    while decision_alnmiss is None:
        try:
            decision_alnmiss = int(raw_input("Please enter an alignment overall missing data threshold (ex. 30): "))
        except ValueError:
            print "That wasn't a number."
    trim_list = []
    fine_list = []
    fh_logc.write("Missing data threshold for alignment = {}.".format(decision_alnmiss)+'\n')
    a_fh = open(file_name, 'r')
    lines = a_fh.readlines()
    for line in lines[1:]:
        cleanline = line.strip()
        sline = cleanline.split('\t')
        contig_name1 = sline[0]
        perc_mdata1 = sline[12]
        if float(perc_mdata1) >= float(decision_alnmiss):
            print "{0} has {1}% missing data.".format(contig_name1, perc_mdata1)
            fh_logc.write("{0} has {1}% missing data.".format(contig_name1, perc_mdata1)+'\n')
            trim_list.append(contig_name1)
        if float(perc_mdata1) < float(decision_alnmiss):
            fine_list.append(contig_name1)
    re_trim = len(trim_list)
    print '\n', '\n', "Number of alignments that require re-trimming = {}.".format(re_trim), '\n', '\n'
    
    trimmed_files = int(0)
    for aln in trim_list:
        for filetype in os.listdir('.'):
            if filetype.endswith(".phylip"):
                name = filetype.split('.')
                if name[0] == aln:
                    trimmed_files += 1
                    print "We've got a match, beginning trimming of {}.".format(filetype)
                    fh_logc.write("Trimmed {}.".format(filetype)+'\n')
                    input_f = filetype
                    output_f = aln+'.trimmed'+'.fasta'
                    program = 'trimal'
                    #trimal_call = "trimal -in {0} -out {1} -phylip -gappyout -keepseqs".format(input_f, output_f)
                    #proc_trimal = sp.call(trimal_call, shell=True)
                    proc_trimal = sp.Popen([program, '-in', '{}'.format(input_f), '-out', '{}'.format(output_f), '-fasta', '-gappyout', '-keepseqs'])
                    proc_trimal.wait()

                    out_name = aln+'.fasta2'
                    fh_fasta = open(out_name, 'a')
                    
                    fh_output = open(output_f,'r')
                    temp_lines = fh_output.readlines()
                    line_count = int(1)
                    for line in temp_lines:
                        clean = line.strip()
                        clean = clean.strip('\r')
                        if clean.startswith('>'):
                            if line_count < 2:
                                fh_fasta.write(clean+'\n')
                            else:
                                fh_fasta.write('\n'+clean+'\n')
                        else:
                            fh_fasta.write(clean)
                        line_count+=1
                    fh_output.close()
                    fh_fasta.close()

                    IDseq_dict = {}
                    with open(out_name, 'r') as fh_fasta2:
                        for line in fh_fasta2:
                            line = line.strip()
                            if line.startswith('>'):
                                name=line[1:]
                                IDseq_dict[name]=''
                            else:
                                IDseq_dict[name]+=line
                    fh_fasta2.close()
                    
                    temp_phylip = aln+".temp.phylip"  
                    fh_temp_phy = open(temp_phylip, 'a')

                    for seq in IDseq_dict:
                        fh_temp_phy.write(seq+' '+IDseq_dict[seq]+'\n')
                    fh_temp_phy.close()
                    
                    taxa_no = int(0)
                    with open(temp_phylip, 'r') as fh_read:
                        for line in fh_read:
                            seq_length = int(0)
                            line = line.strip()
                            split = line.split()
                            taxa = split[0]
                            seq = split[1]
                            for i in seq:
                                seq_length += 1
                            taxa_no += 1

                    taxa_no = str(taxa_no)
                    seq_length = str(seq_length)
                    final_phylip = aln+".final.phylip"
                    fh_final_phy = open(final_phylip, 'a')
                    fh_final_phy.write(taxa_no+' '+seq_length+'\n')
                    with open(temp_phylip, 'r') as fh_read:
                        for line in fh_read:
                            line = line.strip()
                            fh_final_phy.write(line+'\n')
                    fh_final_phy.close()

                    for filetype in os.listdir('.'):
                        if filetype.endswith('fasta2'):
                            os.remove(filetype)
                        if filetype.endswith('temp.phylip'):
                            os.remove(filetype)
                        if filetype.endswith('trimmed.fasta'):
                            os.remove(filetype)
                        if filetype.endswith('final.phylip'):
                            split = filetype.split('.')
                            prefix = split[0]
                            move_name =  out_dir+'/'+prefix+".phylip"
                            shutil.move(filetype, move_name)
                            
    fine_files = int(0)
    for aln in fine_list:
        for filetype in os.listdir('.'):
            if filetype.endswith(".phylip"):
                name = filetype.split('.')
                if name[0] == aln:
                    fine_files += 1
                    new_name = out_dir+'/'+aln+".phylip"
                    shutil.copyfile(filetype, new_name)

    move_name =  out_dir+'/'+log_name
    shutil.move(log_name, move_name)
        
    fh_logc.close()
    print '\n', '\n', "All {0} trimmed files and {1} non-trimmed alignments moved to output directory.".format(trimmed_files,fine_files)
    print "Check log file for the alignments that were trimmed.", '\n', '\n'

###############################################################
#Begin definition of analysis 'd'

elif decision_analysis == "d":
    #begin decisions of thresholds
    print '\n', '------------------------------------------------------------------------'
    decision_bp = None
    while decision_bp is None:
        try:
            decision_bp = int(raw_input("Please enter a minimum base pair value for alignments (ex. 100, or 1 to include all): "))
        except ValueError:
            print "That wasn't a number."
    #next decision
    print '\n','------------------------------------------------------------------------', '\n', "Please select the missing data threshold of an entire alignment to select it for trimming."
    print "For example, entering 30 will cause all alignments with more than 30% overall missing data to be re-trimmed."
    print "If you want to re-trim all alignments with more the more stringent criteria here, just enter 1", '\n'
    decision_alnmiss = None
    while decision_alnmiss is None:
        try:
            decision_alnmiss = int(raw_input("Please enter an alignment overall missing data threshold (ex. 30): "))
        except ValueError:
            print "That wasn't a number."
   
    #next decision
    print '\n','------------------------------------------------------------------------', '\n', "Please select the missing data threshold to remove sequences across alignments.", '\n'
    print "For example, entering 60 will remove a sequence if it is missing more than 60% of the base pairs in the alignment,"
    print "that is, if less than 40% of the bp are present."
    print "In an alignment 1000 bp long, this requires > 400 bases to be present, < 600 bases missing.", '\n'
    decision_indmiss = None
    while decision_indmiss is None:
        try:
            decision_indmiss = int(raw_input("Please enter a sequence missing data threshold (ex. 60): "))
        except ValueError:
            print "That wasn't a number."
    print '------------------------------------------------------------------------', '\n'

    contig_list = []
    trim_list = []
    fine_list = []
    
    a_fh = open(file_name, 'r')
    lines = a_fh.readlines()
    for line in lines[1:]:
        cleanline = line.strip()
        sline = cleanline.split('\t')
        contig_name = sline[0]
        taxa_no = sline[1]
        seq_length = sline[2]
        perc_isites = sline[5]
        perc_mdata = sline[12]
        if int(seq_length) >= decision_bp:
            contig_list.append(contig_name)
    length_pass = len(contig_list)
    length_pass = str(length_pass)
    print "There are {0} alignments that are longer than {1} bp.".format(length_pass, decision_bp), '\n'
    print "Now sorting through these alignments for those that have >{}% missing data to re-trim.".format(decision_alnmiss), '\n'

    log_name = "analysis_c_log.txt"
    fh_logc = open(log_name, 'a')            

    for aln in contig_list:
        for line in lines[1:]:
            cleanline = line.strip()
            sline = cleanline.split('\t')
            contig_name = sline[0]
            taxa_no = sline[1]
            seq_length = sline[2]
            perc_isites = sline[5]
            perc_mdata = sline[12]
            if contig_name == aln:
                if float(perc_mdata) >= float(decision_alnmiss):
                    print "{0} has {1}% missing data.".format(contig_name, perc_mdata)
                    fh_logc.write("{0} has {1}% missing data.".format(contig_name, perc_mdata)+'\n')
                    trim_list.append(contig_name)
                if float(perc_mdata) < float(decision_alnmiss):
                    fine_list.append(contig_name)
            
    re_trim = len(trim_list)
    print '\n', '\n', "Number of alignments that require re-trimming = {}.".format(re_trim), '\n', '\n'
    

    
    ########################################
    #start re-trimmmer
    trimmed_files = int(0)
    for aln in trim_list:
        for filetype in os.listdir('.'):
            if filetype.endswith(".phylip"):
                name = filetype.split('.')
                if name[0] == aln:
                    trimmed_files += 1
                    print "We've got a match, beginning trimming of {}.".format(filetype)
                    fh_logc.write("Trimmed {}.".format(filetype)+'\n')
                    input_f = filetype
                    output_f = aln+'.trimmed'+'.fasta'
                    program = 'trimal'
                    #trimal_call = "trimal -in {0} -out {1} -phylip -gappyout -keepseqs".format(input_f, output_f)
                    #proc_trimal = sp.call(trimal_call, shell=True)
                    proc_trimal = sp.Popen([program, '-in', '{}'.format(input_f), '-out', '{}'.format(output_f), '-fasta', '-gappyout', '-keepseqs'])
                    proc_trimal.wait()

                    out_name = aln+'.fasta2'
                    fh_fasta = open(out_name, 'a')
                    
                    fh_output = open(output_f,'r')
                    temp_lines = fh_output.readlines()
                    line_count = int(1)
                    for line in temp_lines:
                        clean = line.strip()
                        clean = clean.strip('\r')
                        if clean.startswith('>'):
                            if line_count < 2:
                                fh_fasta.write(clean+'\n')
                            else:
                                fh_fasta.write('\n'+clean+'\n')
                        else:
                            fh_fasta.write(clean)
                        line_count+=1
                    fh_output.close()
                    fh_fasta.close()

                    IDseq_dict = {}
                    with open(out_name, 'r') as fh_fasta2:
                        for line in fh_fasta2:
                            line = line.strip()
                            if line.startswith('>'):
                                name=line[1:]
                                IDseq_dict[name]=''
                            else:
                                IDseq_dict[name]+=line
                    fh_fasta2.close()
                    
                    temp_phylip = aln+".temp.phylip"  
                    fh_temp_phy = open(temp_phylip, 'a')

                    for seq in IDseq_dict:
                        fh_temp_phy.write(seq+' '+IDseq_dict[seq]+'\n')
                    fh_temp_phy.close()
                    
                    taxa_no = int(0)
                    with open(temp_phylip, 'r') as fh_read:
                        for line in fh_read:
                            seq_length = int(0)
                            line = line.strip()
                            split = line.split()
                            taxa = split[0]
                            seq = split[1]
                            for i in seq:
                                seq_length += 1
                            taxa_no += 1

                    taxa_no = str(taxa_no)
                    seq_length = str(seq_length)
                    final_phylip = aln+".final.phylip"
                    fh_final_phy = open(final_phylip, 'a')
                    fh_final_phy.write(taxa_no+' '+seq_length+'\n')
                    with open(temp_phylip, 'r') as fh_read:
                        for line in fh_read:
                            line = line.strip()
                            fh_final_phy.write(line+'\n')
                    fh_final_phy.close()

                    for filetype in os.listdir('.'):
                        if filetype.endswith('fasta2'):
                            os.remove(filetype)
                        if filetype.endswith('temp.phylip'):
                            os.remove(filetype)
                        if filetype.endswith('trimmed.fasta'):
                            os.remove(filetype)
                        if filetype.endswith('final.phylip'):
                            split = filetype.split('.')
                            prefix = split[0]
                            move_name =  out_dir+'/'+prefix+".temp.phylip"
                            shutil.move(filetype, move_name)
    fh_logc.close()
                            
    fine_files = int(0)
    for aln in fine_list:
        for filetype in os.listdir('.'):
            if filetype.endswith(".phylip"):
                name = filetype.split('.')
                if name[0] == aln:
                    fine_files += 1
                    new_name = out_dir+'/'+aln+".temp.phylip"
                    shutil.copyfile(filetype, new_name)

    move_name =  out_dir+'/'+log_name
    shutil.move(log_name, move_name)
        
    
    print '\n', '\n', "All {0} trimmed files and {1} non-trimmed alignments moved to output directory.".format(trimmed_files,fine_files)
    print "Check log file for the alignments that were trimmed.", '\n', '\n'


    ########################################
    #start taxon remover
    taxa_removed_l = []
    undisturbed = int(0)
    os.chdir(out_dir)
    fh_log = open('analysis_b_log.txt', 'a')
    for filetype in os.listdir('.'):
        if filetype.endswith("temp.phylip"):
            print "-------------------------------------------------------------------------------"
            print "Refining {}:".format(filetype)
            fh_log.write("Removed the following sequences from {}:".format(filetype)+'\n')
            name = filetype.split('.')
            temp_name = "temp2."+name[0]+".phylip"
            fh_temp = open(temp_name, 'a')
            with open(filetype, 'r') as fh_phy:
                inds = int(0)
                ind_kept = int(0)
                ind_removed = int(0)
                lines = fh_phy.readlines()
                for line in lines[1:]:
                    inds += 1
                    line = line.strip()
                    splitline = line.split()
                    taxon = splitline[0]
                    sequence = splitline[1]
                    total_bps = int(0)
                    gaps = int(0)
                    bps = int(0)
                    
                    for base in sequence:
                        base = base.upper()
                        if base == 'A' or base =='T' or base =='G' or base =='C' or base =='N' or base =='R' or base =='Y' or base =='S' or base =='W' or base =='K' or base =='M' or base =='B' or base =='D' or base =='H' or base =='V':
                            bps+=1
                            total_bps+=1
                        elif base == '-':
                            gaps+=1
                            total_bps+=1
                            
                    cutoff = int( float(( float(gaps) / float(total_bps)) * float(100)))
                    present = int(100) - cutoff
                    if cutoff >= int(decision_indmiss):
                        fh_log.write("Removing {0}".format(taxon)+'\t'+'\t'+"[{0} bp in alignment, {1} missing, {2} present]".format(total_bps, gaps, bps)+'\t'+"[{}% bp present]".format(present)+'\t'+"[{}% missing data]".format(cutoff)+'\n')
                        #print "Removing {0} which is missing {1}% of the data: of {2} bp, {3} missing, {4} present ({5}% present).".format(taxon, cutoff, total_bps, gaps, bps, present)
                        ind_removed += 1
                    elif cutoff < int(decision_indmiss):
                        ind_kept += 1
                        fh_temp.write(line+'\n')
                    taxa_removed_l.append(ind_removed)
                if ind_removed == int(0):
                    undisturbed += 1
            fh_log.write('\n'+'\n')
            print '\n', "Of the {0} sequences, {1} were kept and {2} were removed".format(inds, ind_kept, ind_removed)
            fh_temp.close()
            ind_kept = str(ind_kept)
            total_bps = str(total_bps)
            temp_name2 = name[0]+".phylip"
            fh_temp2 = open(temp_name2, 'a')
            fh_temp2.write(ind_kept+" "+total_bps+'\n')
            with open(temp_name, 'r') as fh_temp:
                lines = fh_temp.readlines()
                for line in lines:
                    line = line.strip()
                    fh_temp2.write(line+'\n')
            fh_temp2.close()
            fh_temp.close()
            print "-------------------------------------------------------------------------------", '\n'
    for filetype in os.listdir('.'):
        if filetype.endswith('temp2.phylip'):
            os.remove(filetype)
        if filetype.endswith('temp.phylip'):
            os.remove(filetype)
        if filetype.startswith('temp2'):
            os.remove(filetype)
    fh_log.close()

    
    fh_decisions = open("analysis_decisions_and_screen_log.txt", 'a')
    fh_decisions.write("USER DECISIONS:"+'\n'+"Minimum length for alignments = {}bp.".format(decision_bp)+'\n'+"Missing data threshold for sequences = {}%.".format(decision_indmiss)+'\n'+"Missing data threshold for alignments = {}%.".format(decision_alnmiss)+'\n'+'\n')
    fh_decisions.write("There were {0} alignments that are longer than {1} bp.".format(length_pass, decision_bp)+'\n'+"Number of alignments that required re-trimming = {}.".format(re_trim)+'\n'+"All {0} trimmed files and {1} non-trimmed alignments moved to output directory.".format(trimmed_files,fine_files)+'\n')
    print '\n', "############################################################################################"
    print "There were {0} alignments that are longer than {1} bp.".format(length_pass, decision_bp), '\n'
    print '\n', "Number of alignments that required re-trimming = {}.".format(re_trim), '\n'
    print '\n', "All {0} trimmed files and {1} non-trimmed alignments moved to output directory.".format(trimmed_files,fine_files)
    print "Check log file for the alignments that were trimmed.", '\n', '\n'

    rem_avg, rem_min, rem_max, rem_med = quickstats(taxa_removed_l)
    print "REMOVED SEQUENCES stats:"
    print "----------------------------------------------------------------------------------"
    print "Average number of taxa removed across alignments = {}.".format(rem_avg)
    print "Median number of taxa removed across alignments = {}.".format(rem_med)
    print "Minimum number of taxa removed across alignments = {}.".format(rem_min)
    print "Maximum number of taxa removed across alignments = {}.".format(rem_max)
    print "There were {} alignments with no sequences removed.".format(undisturbed), '\n'
    print "----------------------------------------------------------------------------------",'\n'
    print '\n', "############################################################################################" '\n', '\n'

    fh_decisions.write("Average number of taxa removed across alignments = {}.".format(rem_avg)+'\n'+"Median number of taxa removed across alignments = {}.".format(rem_med)+'\n'+"Minimum number of taxa removed across alignments = {}.".format(rem_min)+'\n'+"Maximum number of taxa removed across alignments = {}.".format(rem_max)+'\n'+"There were {} alignments with no sequences removed.".format(undisturbed))
    fh_decisions.close()

    print "All finished, time to re-run the 'Alignment Assessment' again to check the alignment refinement.", '\n', '\n'
    
elif decision_analysis == "q":
    print "Quitting program now", '\n'
else:
    print "Not a valid decision, now quitting", '\n'
