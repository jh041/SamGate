###       #! usr/bin/python

print(	"\nWelcome to SamGate (an anagram of MegaSat).\n\n"
		"This script takes the output files\n"
		"from msatcommander and formats a \n"
		"primer file for MegaSat.\n\n"
		"Before we begin, I have to tell you\n"
		"that you'll need the biopython module\n"
		"installed on your machine.\n\n"
		"If you are a windows user,\n"
		"watch this YouTube video on\n" 
		"how to install -  "
		"https://www.youtube.com/watch?v=ddpYVA-7wq4\n\n"
		"If you are a unix user simply type:\n"
		"\"easy_install --user biopython\"\n"
		"Or, if you've downloaded the source code,\n"
		"use the setup.py script like this:\n"
		"python setup.py build\n"
		"python setup.py install --user\n\n"
		"Ok, once biopython is ready I'm gonna need\n"
		"three things from you...\n")

fasta_file = raw_input("First, enter the name of,\n"
						"or pathway to, a fasta file,\n"
						"the one used as input for msatcommander.\n\n"
						"The file should have an extentsion like\n"
						"\".fasta\", \".fna\", or \".fa\"\n\n"
						"If there was more than one fasta file,\n"
						"you should concatenate them into\n"
						"a single file first.\n\n"
						"I'll take that filename from you now:  ")

msats_file = raw_input("\n\nNext I need the name of the\n"
						"microsatellites output file from msatcommander,\n"
						"the one that looks like \"msatcommander.microsatellites\"\n"
						"Tab separated format only please.\n\n"
						"Enter the filename here:   ")

primer_file = raw_input("\n\nLastly, I need the primers file from msatcommander,\n"
						"It should be something like \"msatcommander.primers\"\n\n"
						"Enter the file name here:  ")

###################### compare msat and primer files #####################################

#some msats do not have corresponding primers, these need to be removed first.
#also, I want to make sure that each line in each file corresponds to the other.

new_msats_file = []
new_primer_file = []   

for msat in open(msats_file):
	msat = msat.strip().split()
	for primer in open(primer_file):
		primer = primer.strip().split()
		if msat[0]==primer[0] and msat[2]==primer[2]:  ### name column, change col numbers
			new_msats_file.append(msat)                     ### if your output is different
			new_primer_file.append(primer)
			
#############################Parse fasta file and extract relevant scaffolds##############

scaffold = []

for line in new_msats_file:
	scaffold.append(line[0])
				
def fasta_reader(filename):
	from Bio.SeqIO.FastaIO import FastaIterator
	with open(filename) as handle:
		for record in FastaIterator(handle):
			yield record

#str(entry.id) is the header of fasta entry
#str(entry.seq) is the sequence of specific fasta entry
	
useful_scaffolds = []
		
for entry in fasta_reader(fasta_file):
	if str(entry.id) in scaffold:
		useful_scaffolds.append(entry)

#### Define class type MegaSat ###########################################################

class MegaSat:
	def __init__(
			self, scaffold_name, locus_name, left_primer, rev_comp_right_primer,
			three_flank, five_flank, repeat):
		self.name = scaffold_name
		self.id = locus_name
		self.left = left_primer
		self.right = rev_comp_right_primer
		self.three = three_flank
		self.five = five_flank
		self.repeat = repeat
		
######### extract variables from the fasta, microsatellites and primer files #############

### indices correspond to the column numbers of the mac osx msatcommander output
### if you have a different columns in your output change the indices.   

gates = []

for lineA,lineB in zip(new_msats_file,new_primer_file):
	name = lineA[0]      ### record name (DNA scaffold or sequence), as per fasta file
	id = lineA[2]        ### locus id, not the same thing as a record id
	repeat = lineA[3]    ### repeat motif
	left = lineB[5]      ### forward primer
	right = lineB[14].replace("A","B").replace("T","A").replace("B","T")  ### reverse primer
	right = right.replace("C","D").replace("G","C").replace("D","G")
	right = right[::-1]
	entL = lineB[4].split(',')                   ###forward primer start position,length
	left_start = int(entL[0]) + int(entL[1])
	entR = lineB[13].split(',')                  ###reverse primer start position, length 
	right_start = int(entR[0])
	for fasta in useful_scaffolds:
		if str(fasta.id)==lineB[0]:
			five = str(fasta.seq[left_start:int(lineA[4])])   ###repeat start
			three = str(fasta.seq[int(lineA[5]):right_start])     ###repeat end
			gate = MegaSat(name, id, left, right, three, five, repeat)
			gates.append(gate)
    
#######################print gates into a tab delimited text file#########################

with open('MegaSat_primer.txt', mode = 'w') as file:
	file.write("Locus-Name\tno_tail_left\tRev_Comp_of_no_tail_Right\tThreeFlank\tFiveFlank\tRepeat\n")
	for gate in gates:
		file.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (gate.id, gate.left, gate.right, gate.three, gate.five, gate.repeat))
		
print("\n\nAll done here.\n\n"
		"Move or rename the output file\n"
		"if you don't want it overwritten\n"
		"the next time you use SamGate.\n\n")

##########################################The End#########################################		
