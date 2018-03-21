#!/usr/bin/env python
import sys,time,argparse

def main(args):
	sys.stdout.write("Start analysis: " + time.strftime("%a,%d %b %Y %H:%M:%S") + "\n")
	sys.stdout.flush()
	extract_uniq_read_id(args.input,args.output)
	sys.stdout.write("Finish analysis: " + time.strftime("%a,%d %b %Y %H:%M:%S") + "\n")
	sys.stdout.flush()

def extract_uniq_read_id(input_uniq_gpd_list,output_txt):
	for inf in inf_list:
		for line in inf:
			if line[0] == "@": continue
			qname = line.strip().split("\t")[0]
			print >>output_txt, qname
		inf.close()
	output_txt.close()

def do_inputs():
	output_gpd_format = '''
1. read id'''

	parser = argparse.ArgumentParser(description="Function: extract uniq aligned read id",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-i','--input',type=argparse.FileType('r'),required=True,nargs="+",help="Input: sam file(s) (Note: multiple unique aligned sam file can be as inputs split by 'space'. Currently, only the sam output of GMAP aligner is tested.")
	parser.add_argument('-o','--output',type=argparse.FileType('w'),required=True,help="Output: txt file")
	args = parser.parse_args()
	return args

if __name__=="__main__":
	args = do_inputs()
	main(args)
