#!/usr/bin/python

import re
import sys
import optparse
import csv
from rpy import *

def main():
	parser = optparse.OptionParser()
	parser.add_option('-s', '--snp_vcf', dest = 'snp_vcf', action = 'store', type = 'string', default = None, help = "VCF of SNPs")
	parser.add_option('-c', '--hist_color', dest = 'hist_color', action = 'store', type = 'string', default = "darkgray", help = "Color for 1Mb histograms") 
	parser.add_option('-y', '--ylim', dest = 'ylim', action = 'store', type = 'int', default= 100, help = "Upper limit of Y axis")
	parser.add_option('-z', '--standardize', dest = 'standardize', default= 'false', help = "Standardize X-axis")
	parser.add_option('-e', '--ems', dest = 'ems', default= 'false', help = "Whether EMS variants should be filtered for")
	parser.add_option('-o', '--output', dest = 'plot_output', action = 'store', type = 'string', default = 'EMS_Variant_Density_Plot.pdf', help = "Output file name of plot")
	(options, args) = parser.parse_args()


	i, ii, iii, iv, v, x = parse_snp_vcf(snp_vcf = options.snp_vcf, ems=options.ems)

	create_histograms(plot_output = options.plot_output, hist_color=options.hist_color, ylim=options.ylim, ems=options.ems, standardize=options.standardize, i = i, ii = ii, iii = iii, iv = iv, v = v, x = x)

def create_histograms(plot_output = None, hist_color=None, ylim=None, ems=None, standardize=None , i = None, ii = None, iii = None, iv = None, v = None, x = None):
	breaks = { 'I' : 16 , 'II' : 16,  'III' : 14, 'IV' : 18, 'V' : 21, 'X' : 18 }

	try:
		r.pdf(plot_output, 8, 8)
		if len(i) > 0:
		        plot_data(position_list = i, chr = "I", breaks = breaks["I"], hist_color=hist_color, ylim=ylim, ems=ems, standardize=standardize)
        	if len(ii) > 0:
			plot_data(position_list = ii, chr = "II", breaks = breaks["II"], hist_color=hist_color, ylim=ylim, ems=ems, standardize=standardize)
		if len(iii) > 0:
		        plot_data(position_list = iii, chr = "III", breaks = breaks["III"], hist_color=hist_color, ylim=ylim, ems=ems, standardize=standardize)
        	if len(iv) > 0:
			plot_data(position_list = iv, chr = "IV", breaks = breaks["IV"], hist_color=hist_color, ylim=ylim, ems=ems, standardize=standardize)
		if len(v) > 0:
		        plot_data(position_list = v, chr = "V", breaks = breaks["V"], hist_color=hist_color, ylim=ylim, ems=ems, standardize=standardize)
		if len(x) > 0:
	        	plot_data(position_list = x, chr = "X", breaks = breaks["X"], hist_color=hist_color, ylim=ylim, ems=ems, standardize=standardize)
	        r.dev_off()
    	except Exception as inst:
        	print inst
        	print "There was an error creating the plot pdf... Please try again"

def parse_snp_vcf(snp_vcf = None, ems=None):
	i_file = open(snp_vcf, 'rU')
	reader = csv.reader(i_file, delimiter = '\t', quoting = csv.QUOTE_NONE)
	skip_headers(reader = reader, i_file = i_file)

	i_position_list = []
	ii_position_list = []
	iii_position_list = []
	iv_position_list = []
	v_position_list = []
	x_position_list = []

	for row in reader:
		chromosome = row[0].upper()
		chromosome = re.sub("CHROMOSOME_", "", chromosome, flags = re.IGNORECASE)
		chromosome = re.sub("chr", "", chromosome, flags = re.IGNORECASE)


		position = row[1]
		ref_allele = row[3]
		alt_allele = row[4]

		if  (ems=='true'):
			if (ref_allele =="G" or ref_allele =="C") and (alt_allele =="A" or alt_allele =="T"):
				if chromosome == "I":
					i_position_list.append(position)
				elif chromosome == "II":
					ii_position_list.append(position)
				elif chromosome == "III":
					iii_position_list.append(position)
				elif chromosome == "IV":
					iv_position_list.append(position)
				elif chromosome == "V":
					v_position_list.append(position)
				elif chromosome == "X":
					x_position_list.append(position)
		elif (ems=='false'):
			if chromosome == "I":
				i_position_list.append(position)
			elif chromosome == "II":
				ii_position_list.append(position)
			elif chromosome == "III":
				iii_position_list.append(position)
			elif chromosome == "IV":
				iv_position_list.append(position)
			elif chromosome == "V":
				v_position_list.append(position)
			elif chromosome == "X":
				x_position_list.append(position)

	return i_position_list, ii_position_list, iii_position_list, iv_position_list, v_position_list, x_position_list

def skip_headers(reader = None, i_file = None):
	# count headers
	comment = 0
	while reader.next()[0].startswith('#'):
		comment = comment + 1
	
	# skip headers
	i_file.seek(0)
	for i in range(0, comment):
		reader.next()

def plot_data(position_list = None, chr = None, breaks = None, hist_color=None, ylim = None, ems=None, standardize=None):
	positions = ",".join(map(str, map(lambda x: float(x) / 1000000, position_list)))
	positions = "c(" + positions + ")"
	
	if (standardize=='true'):
		r("hist(" + positions + ", xlim=c(0,21), ylim=c(0, %d "%ylim +"),col='"+ hist_color + "', breaks = seq(0, as.integer( ' " + str(breaks) + " '), by=1), main = 'LG " + chr + "', ylab = 'Frequency Of SNPs', xlab = 'Location (Mb)')")
		r("hist(" + positions + ", xlim=c(0,21), add=TRUE,  ylim=c(0, %d "%ylim +"), col=rgb(1, 0, 0, 1), breaks = seq(0, as.integer( ' " + str(breaks) + " '), by=.5), main = 'Chr " + chr + "', ylab = 'Number Of SNPs', xlab = 'Location (Mb)')")
		r("axis(1, at=seq(0, 21, by=1), labels=FALSE, tcl=-0.5)")
		r("axis(1, at=seq(0, 21, by=0.5), labels=FALSE, tcl=-0.25)")
	elif (standardize=='false'):
		r("hist(" + positions + ", xlim=c(0,as.integer( ' " + str(breaks) + " ')), ylim=c(0, %d "%ylim +"),col='"+ hist_color + "', breaks = seq(0, as.integer( ' " + str(breaks) + " '), by=1), main = 'LG " + chr + "', ylab = 'Frequency Of SNPs', xlab = 'Location (Mb)')")
		r("hist(" + positions + ", xlim=c(0,as.integer( ' " + str(breaks) + " ')), add=TRUE,  ylim=c(0, %d "%ylim +"), col=rgb(1, 0, 0, 1), breaks = seq(0, as.integer( ' " + str(breaks) + " '), by=.5), main = 'Chr " + chr + "', ylab = 'Number Of SNPs', xlab = 'Location (Mb)')")
		r("axis(1, at=seq(0, as.integer( ' " + str(breaks) + " '), by=1), labels=FALSE, tcl=-0.5)")
		r("axis(1, at=seq(0, as.integer( ' " + str(breaks) + " '), by=0.5), labels=FALSE, tcl=-0.25)")



if __name__ == "__main__":
	main()
