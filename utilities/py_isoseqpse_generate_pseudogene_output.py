#!/usr/bin/env python
import sys,re,time,argparse

def main(args):
	sys.stdout.write("Start analysis: " + time.strftime("%a,%d %b %Y %H:%M:%S") + "\n")
	sys.stdout.flush()
	dic_read_label = parse_alignment_label(args.input_lr)
	pgene_list = parse_pseudogene_list(args.input_pg)
	generate_pseudogene_output(args.input_iso,args.output)
	sys.stdout.write("Finish analysis: " + time.strftime("%a,%d %b %Y %H:%M:%S") + "\n")
	sys.stdout.flush()

def parse_alignment_label(input_lr_gpd):
	global dic_read_label
	dic_read_label = {}
	for line in input_lr_gpd:
		read_id = line.strip().split("\t")[0]
		alignment_label = line.strip().split("\t")[-1]
		if read_id.endswith("F1P1T1"):
			if alignment_label == "U":
				dic_read_label[read_id] = "FL_U"
			elif alignment_label == "M":
				dic_read_label[read_id] = "FL_M"
			else:
				sys.stderr.write("Error, no alignment label for the read below\n"+read_id+"\n")
		else:
			if alignment_label == "U":
				dic_read_label[read_id] = "nFL_U"
			elif alignment_label == "M":
				dic_read_label[read_id] = "nFL_M"
			else:
				sys.stderr.write("Error, no alignment label for the read below\n"+read_id+"\n")
	input_lr_gpd.close()
	return dic_read_label

def parse_pseudogene_list(input_pseu_list):
	global pgene_list
	pgene_list = []
	for inf	in input_pseu_list:
		for line in inf:
			pgene = line.strip().split("\t")[0]
			pgene_list.append(pgene)
		inf.close()
	return pgene_list

def generate_pseudogene_output(input_iso_gpd,output_gpd):
	for line in input_iso_gpd:
		gene_id_set,iso_id_set = line.strip().split("\t")[:2]
		lr_id_set = line.strip().split("\t")[-1]
		fl_u,fl_m,nfl_u,nfl_m = 0,0,0,0
		for gene_iso in (gene_id_set.split(",") + iso_id_set.split(",")):
			if gene_iso in pgene_list:
				for read in lr_id_set.split(","):
					if dic_read_label[read] == "FL_U":
						fl_u += 1
					elif dic_read_label[read] == "FL_M":
						fl_m += 1
					elif dic_read_label[read] == "nFL_U":
						nfl_u += 1
					else:
						nfl_m += 1
				print >>output_gpd, "\t".join(["\t".join(line.strip().split("\t")[:-1]),str(fl_u),str(fl_m),str(fl_u+nfl_u),str(fl_m+nfl_m)])
				break
	input_iso_gpd.close()
	output_gpd.close()

def do_inputs():
	output_gpd_format = '''
1. gene id
2. isoform id
3. chromosome id
4. strand
5. TSS (+)
6. TTS (+)
7. number of support full-length long reads
8. number of support total long reads
9. exon count
10. exon start set
11. exon end set
12. For novel isoform, derived genic locus
13. For novel isoform, overlap percentage with derived genic locus
14. For novel singleton isoform, if it is located at the last exon of any known isoform. If yes, isoform ID otherwise '-'
15. For novel singleton isoform, the overlap percentage with the the last exon
16. For novel multi-exon isoform, number of splice sites are detected by anno and/or short reads; and the total number of splice sites
17. For novel multi-exon isoform, if the multi-exon isoform is the subset (based on splice junction combination) of known multi-exon isoform, isoform ID if yes otherwise '-'
18. For novel isoform, maximal length of polyA track in defined region
19. For novel isoform, maximal percentage of nucleotide A in defined region
20. Number of supporting full-length long read with unique alignment
21. Number of supporting full-length long read with multiple alignment
22. Number of supporting total long read with unique alignment
23. Number of supporting total long read with multiple alignment'''

	parser = argparse.ArgumentParser(description="Function: get best one alignment for each read",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-i','--input_iso',type=argparse.FileType('r'),required=True,help="Input: constructed isoform (GPD file) generated by 'py_isoseqpse_generate_output.py'")
	parser.add_argument('-l','--input_lr',type=argparse.FileType('r'),required=True,help="Input: long read (GPD file) with best alignment generated by 'py_isoseqpse_best_alignment.py'")
	parser.add_argument('-g','--input_pg',type=argparse.FileType('r'),required=True,nargs="+",help="Input: pseudogene list (TXT files, tab-split), the first column showing the pseudogene gene id and/or transcript id")
	parser.add_argument('-o','--output',type=argparse.FileType('w'),required=True,help="Output: constructed pseudogenes (GPD file)")
	args = parser.parse_args()
	return args

if __name__=="__main__":
	args = do_inputs()
	main(args)
