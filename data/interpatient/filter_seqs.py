from Bio import SeqIO
import glob
import subprocess
import os

os.system('mkdir fasta')

## filter by length
total = 0
length  = 0
for record in SeqIO.parse('ncbi_mar20.fna', 'fasta'):
	total += 1

	if len(record.seq) >29000:
		length += 1
		f = open('fasta/' + str(record.id) + ".fasta", 'w+')
		f.write(">" + record.description + "\n")
		f.write(str(record.seq) + "\n")
		f.close()
print(str(total) + " Total sequences")
print(str(length) + " Were greater than 29 Kb")


### now run fastANI

ani = 0
os.system('rm ./final_filtered_seqs.fna')
for fn in glob.glob('./fasta/*.fasta'):
	out = './fasta/' + fn.split("/")[-1].split(".fasta")[0] + ".out"
	with open(os.devnull, 'w') as devnull:
		subprocess.call(['./fastANI', '-q', fn, '-r', 'reference.fna', '-o', out], stdout=devnull, stderr= devnull)

	f2 = open(out)
	result = f2.readline()
	f2.close()
	if float(result.split()[2]) > 99:
		os.system('cat ' + fn + " >> ./final_filtered_seqs.fna")
		ani += 1

print(str(ani) + " Were >99% nucleotide identity to the reference.")
