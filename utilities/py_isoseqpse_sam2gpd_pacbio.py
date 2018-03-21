#!/usr/bin/env python
import sys,re,time,argparse
from multiprocessing import cpu_count,Pool

def main(args):
	sys.stdout.write("Start analysis: " + time.strftime("%a,%d %b %Y %H:%M:%S") + "\n")
	sys.stdout.flush()
	f1p1t1,f1p0t1,f1p1t0,f1p0t0,f0p1t1,f0p0t1,f0p1t0,f0p0t0 = parse_primer_info(args.csv)
	dic_chr_seq = parse_reference_genome(args.genome)
	gpd_file = args.output
	p = Pool(processes=args.cpu)
	csize = 1000
	results = p.imap(func=convert,iterable=generate_tx(args.input),chunksize=csize)
	for res in results:
		if not res: continue
		gpd_file.write(res+"\n")
	gpd_file.close()
	sys.stdout.write("Finish analysis: " + time.strftime("%a,%d %b %Y %H:%M:%S") + "\n")
	sys.stdout.flush()

def parse_primer_info(primer_csv_fl):
	global f1p1t1,f1p0t1,f1p1t0,f1p0t0,f0p1t1,f0p0t1,f0p1t0,f0p0t0
	f1p1t1 = []
	f1p0t1 = []
	f1p1t0 = []
	f1p0t0 = []
	f0p1t1 = []
	f0p0t1 = []
	f0p1t0 = []
	f0p0t0 = []
	head = 1
	for line in primer_csv_fl:
		if head:
			head -= 1
			continue
		id,strand,fiveseen,polyAseen,threeseen,fiveend,polyAend,threeend,primer,chimera = line.rstrip("\n").split(",")
		if chimera == "0":
			if fiveseen == "1" and polyAseen == "1" and threeseen == "1":
				f1p1t1.append(id)
			elif fiveseen == "1" and polyAseen == "0" and threeseen == "1":
				f1p0t1.append(id)
			elif fiveseen == "1" and polyAseen == "1" and threeseen == "0":
				f1p1t0.append(id)
			elif fiveseen == "1" and polyAseen == "0" and threeseen == "0":
				f1p0t0.append(id)
			elif fiveseen == "0" and polyAseen == "1" and threeseen == "1":
				f0p1t1.append(id)
			elif fiveseen == "0" and polyAseen == "0" and threeseen == "1":
				f0p0t1.append(id)
			elif fiveseen == "0" and polyAseen == "1" and threeseen == "0":
				f0p1t0.append(id)
			else:
				f0p0t0.append(id)
	primer_csv_fl.close()
	return f1p1t1,f1p0t1,f1p1t0,f1p0t0,f0p1t1,f0p0t1,f0p1t0,f0p0t0

def parse_reference_genome(genome_fasta):
	global dic_chr_seq
	dic_chr_seq = {}
	chr = ""
	seq_list = [""]
	for line in genome_fasta:
		if line.startswith(">"):
			dic_chr_seq[chr] = ("".join(seq_list)).upper()
			seq_list = []
			chr = re.search("^>(\w+)",line.strip()).group(1)
		else:
			seq_list.append(line.strip())
	dic_chr_seq[chr] = ("".join(seq_list)).upper()
	del dic_chr_seq[""]
	genome_fasta.close()
	return dic_chr_seq

def extract_exon_length_from_cigar(cigar):
	cigar_m = ["0"] + re.findall(r"(\d+)M",cigar)
	cigar_d = ["0"] + re.findall(r"(\d+)D",cigar)
	cigar_m_s,cigar_d_s = [0,0]
	for m in cigar_m:
		cigar_m_s += int(m)
	for d in cigar_d:
		cigar_d_s += int(d)
	exon_length = cigar_m_s+cigar_d_s
	return exon_length

def extract_soft_clip_from_cigar(cigar):
	cigar_5 = ["0"] + re.findall(r"^(\d+)S",cigar)
	cigar_3 = ["0"] + re.findall(r"(\d+)S$",cigar)
	cigar_5_s,cigar_3_s = [0,0]
	for s5 in cigar_5:
		cigar_5_s += int(s5)
	for s3 in cigar_3:
		cigar_3_s += int(s3)
	return cigar_5_s,cigar_3_s

def extract_indel_from_cigar(cigar):
	cigar_ins = ["0"] + re.findall(r"(\d+)I",cigar)
	cigar_del = ["0"] + re.findall(r"(\d+)D",cigar)
	cigar_ins_s,cigar_del_s = 0,0
	for inso in cigar_ins:
		cigar_ins_s += int(inso)
	for delo in cigar_del:
		cigar_del_s += int(delo)
	return cigar_ins_s,cigar_del_s

def get_aligned_read_sequence(cigar,seq):
	cigar_list = re.split("(\d+)",cigar)[1:]
	seq_list = []
	for i in range(1,len(cigar_list),2):
		if cigar_list[i] == "S" or cigar_list[i] == "I":
			seq = seq[int(cigar_list[i-1]):]
		elif cigar_list[i] == "M":
			seq_new = seq[:int(cigar_list[i-1])]
			seq_list.append(seq_new)
			seq = seq[int(cigar_list[i-1]):]
		elif cigar_list[i] == "D":
			seq_new = "-"*int(cigar_list[i-1])
			seq_list.append(seq_new)
		else:
			pass
	seq_output = "".join(seq_list)
	return seq_output

def get_aligned_reference_sequence(dic_chr_seq,chr,exon_number,exon_start_list,exon_end_list):
	seq_list = []
	for i in range(0,int(exon_number)):
		start = int(exon_start_list[i])
		end = int(exon_end_list[i])
		seq_list.append(dic_chr_seq[chr][start:end])
	seq = "".join(seq_list)
	return seq


def generate_tx(inf_list):
	z = 0
	for inf in inf_list:
		for line in inf:
			if line[0] == "@": continue
			z += 1
			yield (line,z)
		inf.close()

def convert(inputs):
	(line,z) = inputs
	qname,flag,rname,pos,mapq,cigar,rnext,pnext,tlen,seq = line.strip().split("\t")[:10]
	tag = "\t".join(line.strip().split("\t")[11:])
	if (rname != "*") and ("M" in cigar) and ("H" not in cigar): # aligned and non-hard-clipped reads
		# soft-clip insertion and deletion info
		s5,s3 = extract_soft_clip_from_cigar(cigar)
		ins_s,del_s = extract_indel_from_cigar(cigar)
		sf = str(s5)+"_"+str(s3)
		# exon info
		cigar_n_l = 0
		exon_length = 0
		exon_start = int(pos)-1
		exon_end = 0
		exon_start_list = []
		exon_end_list = []
		if "N" in cigar:
			for exon in cigar.split("N"):
				exon = exon + "N"
				exon_start = exon_start + exon_length + cigar_n_l
				exon_length = extract_exon_length_from_cigar(exon)
				exon_end = exon_start + exon_length
				if re.search(r"(\d+)N",exon):
					cigar_n_l = int((re.search(r"(\d+)N",exon)).group(1))
				exon_start_list.append(str(exon_start))
				exon_end_list.append(str(exon_end))
		else:
			exon_start = exon_start
			exon_length = extract_exon_length_from_cigar(cigar)
			exon_end = exon_start + exon_length
			exon_start_list.append(str(exon_start))
			exon_end_list.append(str(exon_end))
		exon_start_list.append("")
		exon_end_list.append("")

		# sequence identity between read and reference info
		exon_number = str(len(exon_start_list)-1)
		read_seq = get_aligned_read_sequence(cigar,seq)
		ref_seq = get_aligned_reference_sequence(dic_chr_seq,rname,exon_number,exon_start_list,exon_end_list)
		diff_count = sum(1 for a, b in zip(read_seq,ref_seq) if a!=b)
		err_in_read = str(s5 + s3) + "_" + str(ins_s) + "_" + str(del_s) + "_" + str(diff_count-del_s) + "_" + str(len(seq))
		err_in_ref = str(del_s) + "_" + str(diff_count-del_s) + "_" + str(len(ref_seq))
		err_rate_in_read = str(float(sum([s5,s3,ins_s,diff_count]))/len(seq))
		err_rate_in_ref = str(float(diff_count)/len(ref_seq))
		# strand info
		if qname not in (f0p1t0 + f0p0t0):
			if flag == "0":
				strand = "+"
			else:
				strand = "-"
		else:
			if "N" in cigar:
				if re.search(r"XS:A:(\S)",tag):
					strand = (re.search(r"XS:A:(\S)",tag)).group(1)
				else:
					strand = "*"
			else:
				strand = "*"
		# rename read ID
		if qname in f1p1t1:
			qname = qname + "_F1P1T1"
		elif qname in f1p1t0:
			qname = qname + "_F1P1T0"
		elif qname in f1p0t1:
			qname = qname + "_F1P0T1"
		elif qname in f1p0t0:
			qname = qname + "_F1P0T0"
		elif qname in f0p1t1:
			qname = qname + "_F0P1T1"
		elif qname in f0p0t1:
			qname = qname + "_F0P0T1"
		elif qname in f0p1t0:
			qname = qname + "_F0P1T0"
		else:
			qname = qname + "_F0P0T0"

		gpd_output = "\t".join([qname,qname,rname,strand,str(int(pos)-1),str(exon_end),mapq,sf,str(len(exon_start_list)-1),",".join(exon_start_list),",".join(exon_end_list),flag,err_in_read,err_in_ref,err_rate_in_read,err_rate_in_ref])
	else:
		gpd_output = None
	return gpd_output

def do_inputs():
	output_gpd_format = '''
1. read id
2. read id
3. chromosome id
4. strand
5. start site of alignment
6. end site of alignment
7. MAPQ 
8. number of nucleotides that are softly-clipped by aligner; left_right
9. exon count
10. exon start set
11. exon end set
12. sam flag
13. error in read sequence (soft-clip,insertion,deletion,mismatch)
14. error in reference genome sequence (deletion,mismatch)
15. error rate in read sequence (soft-clip,insertion,deletion,mismatch)
16. error rate in reference genome sequence (deletion,mismatch)'''
	parser = argparse.ArgumentParser(description="Function: convert sam to gpd, specifically for PacBio transcriptome long-read data after the processes: 1) extract ROI by PacBio SMRTanalysis2.3 'ConsensusTools.sh CircularConsensus'; 2) classifiy ROI by PacBio SMRTanalysis2.3 'pbtranscript.py classify' with the parameter '--detect_chimera_nfl' to get full-length (FL) and non-full-length (nFL), non-chimera (NC) ROI; 3) separate the nFLNC by the script './utilities/py_isoseqcon_separate_nflnc_fasta.py' to get nFLNC_sense.fasta and nFLNC_unknown.fasta; and 4) align ROI by GMAP aligner with the parameter '-z sense_force -f samse -n 5 --split-output' for FLNC and nFLNC_sense and the parameter '-f samse -n 5 --split-output' for nFLNC_unknown.",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-i','--input',type=argparse.FileType('r'),required=True,nargs="+",help="Input: sam file(s) (Note: multiple sam file can be as inputs split by 'space'. Currently, only the sam output of GMAP aligner is tested.")
	parser.add_argument('-c','--csv',type=argparse.FileType('r'),required=True,help="Input: csv file (i.e., *.primer_info.csv file generated by 'pbtranscript.py classify' process)")
	parser.add_argument('-g','--genome',type=argparse.FileType('r'),required=True,help="Input: reference genome file (fasta format) for calculating the identity between aligned read and reference genome")
	parser.add_argument('-o','--output',type=argparse.FileType('w'),required=True,help="Output: gpd file")
	parser.add_argument('-p','--cpu',type=int,default=cpu_count(),help="Number of cpus")
	args = parser.parse_args()
	return args

if __name__=="__main__":
	args = do_inputs()
	main(args)
