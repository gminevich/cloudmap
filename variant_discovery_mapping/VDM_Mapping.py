#!/usr/bin/python

import re
import sys
import optparse
import csv
import re
import pprint
from decimal import *
from rpy import *

def main():
	csv.field_size_limit(1000000000)

	parser = optparse.OptionParser()
	parser.add_option('-v', '--sample_vcf', dest = 'sample_vcf', action = 'store', type = 'string', default = None, help = "Sample VCF from GATK Unified Genotyper")
	parser.add_option('-l', '--loess_span', dest = 'loess_span', action = 'store', type = 'float', default = .1, help = "Loess span")
	parser.add_option('-d', '--d_yaxis', dest = 'd_yaxis', action = 'store', type = 'float', default = 1, help = "y-axis upper limit for dot plot")  
	parser.add_option('-y', '--h_yaxis', dest = 'h_yaxis', action = 'store', type = 'int', default = 0, help = "y-axis upper limit for histogram plot")   
	parser.add_option('-c', '--points_color', dest = 'points_color', action = 'store', type = 'string', default = "gray27", help = "Color for data points") 
	parser.add_option('-k', '--loess_color', dest = 'loess_color', action = 'store', type = 'string', default = "green2", help = "Color for loess regression line")        
	parser.add_option('-z', '--standardize', dest = 'standardize', default= 'true', help = "Standardize X-axis")
	parser.add_option('-b', '--break_file', dest = 'break_file', action = 'store', type = 'string', default = 'C.elegans', help = "File defining the breaks per chromosome")
	parser.add_option('-x', '--bin_size', dest = 'bin_size', action = 'store', type = 'int', default = 1000000, help = "Size of histogram bins, default is 1mb")
	parser.add_option('-n', '--normalize_bins', dest = 'normalize_bins', default= 'true', help = "Normalize histograms")

	parser.add_option('-o', '--output', dest = 'output', action = 'store', type = 'string', default = None, help = "Output file name")
	parser.add_option('-s', '--location_plot_output', dest = 'location_plot_output', action = 'store', type = 'string', default = "SNP_Mapping_Plot.pdf", help = "Output file name of SNP plots by chromosomal location")

	(options, args) = parser.parse_args()

	vcf_info = parse_vcf(sample_vcf = options.sample_vcf)

	output_vcf_info(output = options.output, vcf_info = vcf_info)
	
	#output plot with all ratios
	rounded_bin_size = int(round((float(options.bin_size) / 1000000), 1) * 1000000)
	
	normalized_histogram_bins_per_mb = calculate_normalized_histogram_bins_per_xbase(vcf_info = vcf_info, xbase = rounded_bin_size, normalize_bins = options.normalize_bins)
	max_y_hist_mb = normalized_histogram_bins_per_mb[max(normalized_histogram_bins_per_mb, key = lambda x: normalized_histogram_bins_per_mb.get(x) )]

	normalized_histogram_bins_per_5kb = calculate_normalized_histogram_bins_per_xbase(vcf_info = vcf_info, xbase = (rounded_bin_size / 2), normalize_bins = options.normalize_bins)
	max_y_hist_5kb = normalized_histogram_bins_per_5kb[max(normalized_histogram_bins_per_5kb, key = lambda x: normalized_histogram_bins_per_5kb.get(x) )]

	max_y_hist_overall = myround(max(max_y_hist_mb, max_y_hist_5kb) + int(round(round(max(max_y_hist_mb, max_y_hist_5kb)) * .1)))


	break_dict = parse_breaks(break_file = options.break_file)

	output_scatter_plots_by_location(location_plot_output = options.location_plot_output, vcf_info = vcf_info, loess_span=options.loess_span, d_yaxis=options.d_yaxis, h_yaxis=options.h_yaxis, points_color=options.points_color, loess_color=options.loess_color, standardize =options.standardize, normalized_hist_per_mb = normalized_histogram_bins_per_mb, normalized_hist_per_5kb = normalized_histogram_bins_per_5kb, breaks = break_dict, rounded_bin_size = rounded_bin_size, max_y_hist_overall = max_y_hist_overall)

def myround(x, base=10):
    return int(base * round(float(x)/base))

def skip_headers(reader = None, i_file = None):
	# count headers
	comment = 0
	while reader.next()[0].startswith('#'):
		comment = comment + 1
	
	# skip headers
	i_file.seek(0)
	for i in range(0, comment):
		reader.next()

def parse_breaks(break_file = None):
	if break_file == 'C.elegans':
		break_dict = { 'I' : 16 , 'II' : 16,  'III' : 14, 'IV' : 18, 'V' : 21, 'X' : 18 }
		return break_dict
	elif break_file == 'Brachypodium':
		break_dict = { '1' : 75 , '2' : 60,  '3' : 60, '4' : 50, '5' : 30 }
		return break_dict
	elif break_file == 'Arabadopsis':
		break_dict = { '1' : 31 , '2' : 20,  '3' : 24, '4' : 19, '5' : 27 }
		return break_dict
	else:
		i_file = open(break_file, 'rU')
		break_dict = {}
		reader = csv.reader(i_file, delimiter = '\t')
		for row in reader:
			chromosome = row[0].upper()
			chromosome = re.sub("CHROMOSOME_", "", chromosome, flags = re.IGNORECASE)
			chromosome = re.sub("chr", "", chromosome, flags = re.IGNORECASE)
			break_count = row[1]
			break_dict[chromosome] = int(break_count)
		return break_dict


def location_comparer(location_1, location_2):
	chr_loc_1 = location_1.split(':')[0]
	pos_loc_1 = int(location_1.split(':')[1])

	chr_loc_2 = location_2.split(':')[0]
	pos_loc_2 = int(location_2.split(':')[1])

	if chr_loc_1 == chr_loc_2:
		if pos_loc_1 < pos_loc_2:
			return -1
		elif pos_loc_1 == pos_loc_1:
			return 0
		elif pos_loc_1 > pos_loc_2:
			return 1
	elif chr_loc_1 < chr_loc_2:
		return -1
	elif chr_loc_1 > chr_loc_2:
		return 1

def output_vcf_info(output = None, vcf_info = None):
	o_file = open(output, 'wb')
	writer = csv.writer(o_file, delimiter = '\t')

	writer.writerow(["#Chr\t", "Pos\t", "Alt Count\t", "Ref Count\t", "Read Depth\t", "Ratio\t"])

	location_sorted_vcf_info_keys = sorted(vcf_info.keys(), cmp=location_comparer)

	for location in location_sorted_vcf_info_keys:
		alt_allele_count, ref_allele_count, read_depth, ratio = vcf_info[location]
		
		location_info = location.split(':')
		chromosome = location_info[0]
		position = location_info[1]

		writer.writerow([chromosome, position, alt_allele_count, ref_allele_count, read_depth, ratio])

	o_file.close()

def output_scatter_plots_by_location(location_plot_output = None, vcf_info = None, loess_span="", d_yaxis="", h_yaxis="", points_color="", loess_color="", standardize=None, normalized_hist_per_mb = None, normalized_hist_per_5kb = None, breaks = None, rounded_bin_size = 1000000, max_y_hist_overall = ""):
	positions = {}
	current_chr = ""
	prev_chr = ""

	x_label = "Location (Mb)"
	filtered_label = ''

	location_sorted_vcf_info_keys = sorted(vcf_info.keys(), cmp=location_comparer)
	
	break_unit = Decimal(rounded_bin_size) / Decimal(1000000)
	max_breaks = max(breaks.values())

	try:
		r.pdf(location_plot_output, 8, 8)
	
		for location in location_sorted_vcf_info_keys:
			current_chr = location.split(':')[0]
			position = location.split(':')[1]

			alt_allele_count, ref_allele_count, read_depth, ratio = vcf_info[location]
		
			if prev_chr != current_chr:
				if prev_chr != "":
					hist_dict_mb = get_hist_dict_by_chr(normalized_hist_per_xbase = normalized_hist_per_mb, chr = prev_chr)
					hist_dict_5kb = get_hist_dict_by_chr(normalized_hist_per_xbase = normalized_hist_per_5kb, chr = prev_chr)
					
					if h_yaxis == 0:
						plot_data(chr_dict = positions, hist_dict_mb = hist_dict_mb, hist_dict_5kb = hist_dict_5kb, chr = prev_chr + filtered_label, x_label = "Location (Mb)", divide_position = True, draw_secondary_grid_lines = True, loess_span=loess_span, d_yaxis=d_yaxis, h_yaxis=max_y_hist_overall, points_color=points_color, loess_color=loess_color, breaks = breaks[prev_chr], standardize=standardize, max_breaks = max_breaks, break_unit = break_unit)
					else:
						plot_data(chr_dict = positions, hist_dict_mb = hist_dict_mb, hist_dict_5kb = hist_dict_5kb, chr = prev_chr + filtered_label, x_label = "Location (Mb)", divide_position = True, draw_secondary_grid_lines = True, loess_span=loess_span, d_yaxis=d_yaxis, h_yaxis=h_yaxis, points_color=points_color, loess_color=loess_color, breaks = breaks[prev_chr], standardize=standardize, max_breaks = max_breaks, break_unit = break_unit)

				prev_chr = current_chr
				positions = {}
		
			positions[position] = ratio

		hist_dict_mb = get_hist_dict_by_chr(normalized_hist_per_xbase = normalized_hist_per_mb, chr = current_chr)
		hist_dict_5kb = get_hist_dict_by_chr(normalized_hist_per_xbase = normalized_hist_per_5kb, chr = current_chr)
			
		if h_yaxis == 0:					
			plot_data(chr_dict = positions, hist_dict_mb = hist_dict_mb, hist_dict_5kb = hist_dict_5kb, chr = current_chr + filtered_label, x_label = "Location (Mb)", divide_position = True, draw_secondary_grid_lines = True, loess_span=loess_span, d_yaxis=d_yaxis, h_yaxis=max_y_hist_overall, points_color=points_color, loess_color=loess_color, breaks = breaks[current_chr], standardize=standardize, max_breaks = max_breaks, break_unit = break_unit)
		else:
			plot_data(chr_dict = positions, hist_dict_mb = hist_dict_mb, hist_dict_5kb = hist_dict_5kb, chr = current_chr + filtered_label, x_label = "Location (Mb)", divide_position = True, draw_secondary_grid_lines = True, loess_span=loess_span, d_yaxis=d_yaxis, h_yaxis=h_yaxis, points_color=points_color, loess_color=loess_color, breaks = breaks[current_chr], standardize=standardize, max_breaks = max_breaks, break_unit = break_unit)

		r.dev_off()
		
	except Exception as inst:
        	print inst
        	print "There was an error creating the location plot pdf... Please try again"

def get_hist_dict_by_chr(normalized_hist_per_xbase = None, chr = ''):
	hist_dict = {}	

	for location in normalized_hist_per_xbase:
		chromosome = location.split(':')[0]		
		if chromosome == chr:
			position = int(location.split(':')[1])
			hist_dict[position] = normalized_hist_per_xbase[location]
	
	max_location = max(hist_dict.keys(), key=int)
	for i in range(1, max_location):
		if i not in hist_dict:
			hist_dict[i] = 0	
	
	return hist_dict	

def plot_data(chr_dict =  None, hist_dict_mb = None, hist_dict_5kb = None, chr = "", x_label = "", divide_position = False, draw_secondary_grid_lines = False, loess_span=None, d_yaxis=None, h_yaxis=None, points_color="", loess_color="", breaks = None, standardize= None, max_breaks = 1, break_unit = 1):
	ratios = "c("
	positions = "c("
	
	for position in chr_dict:
		ratio = chr_dict[position]
		if divide_position:
		       	position = float(position) / 1000000.0
	        positions = positions + str(position) + ", "
		ratios = ratios + str(ratio) + ", "

	if len(ratios) == 2:
		ratios = ratios + ")"
	else:
		ratios = ratios[0:len(ratios) - 2] + ")"

	if len(positions) == 2:
		positions = positions + ")"
	else:
		positions = positions[0:len(positions) - 2] + ")"

	r("x <- " + positions)
	r("y <- " + ratios)

	hist_mb_values = "c("
    	for position in sorted(hist_dict_mb):
		hist_mb_values = hist_mb_values + str(hist_dict_mb[position]) + ", "
	
	if len(hist_mb_values) == 2:
		hist_mb_values = hist_mb_values + ")"
	else:
		hist_mb_values = hist_mb_values[0:len(hist_mb_values) - 2] + ")"

	hist_5kb_values = "c("
	for position in sorted(hist_dict_5kb):
		hist_5kb_values = hist_5kb_values + str(hist_dict_5kb[position]) + ", "	

	if len(hist_5kb_values) == 2:
		hist_5kb_values = hist_5kb_values + ")"
	else:
		hist_5kb_values = hist_5kb_values[0:len(hist_5kb_values) - 2] + ")"

	r("xz <- " + hist_mb_values)
	r("yz <- " + hist_5kb_values)


	max_break_str = str(max_breaks)
	break_unit_str = str(Decimal(break_unit)) 	
	half_break_unit_str = str(Decimal(break_unit) / Decimal(2))
	break_penta_unit_str = str(Decimal(break_unit) * Decimal(5))

	if (standardize=='true'):  
		r("plot(x, y, ,cex=0.60, xlim=c(0," + max_break_str + "), main='LG " + chr + " (Variant Discovery Mapping)', xlab= '" + x_label + "', ylim = c(0, %f " %d_yaxis + "), ylab='Ratios of variant reads/total reads (at variant positions)', pch=10, col='"+ points_color +"')")
		r("lines(loess.smooth(x, y, span = %f "%loess_span + "), lwd=5, col='"+ loess_color +"')")
		r("axis(1, at=seq(0, " + max_break_str + ", by=" + break_unit_str + "), labels=FALSE, tcl=-0.5)")
		r("axis(1, at=seq(0, " + max_break_str + ", by=" + half_break_unit_str + "), labels=FALSE, tcl=-0.25)")
		r("axis(2, at=seq(floor(min(y)), 1, by=0.1), labels=FALSE, tcl=-0.2)")
	elif (standardize=='false'):
		r("plot(x, y, cex=0.60, main='LG " + chr + " (Variant Discovery Mapping)', xlab= '" + x_label + "', ylim = c(0, %f " %d_yaxis + "), ylab='Ratios of variant reads/total reads (at variant positions)', pch=10, col='"+ points_color +"')")
		r("lines(loess.smooth(x, y, span = %f "%loess_span + "), lwd=5, col='"+ loess_color +"')")    
		r("axis(1, at=seq(0, as.integer( ' " + str(breaks) + " '), by= " + break_unit_str + "), labels=FALSE, tcl=-0.5)")
		r("axis(1, at=seq(0, as.integer( ' " + str(breaks) + " '), by= " + half_break_unit_str + "), labels=FALSE, tcl=-0.25)")	
		r("axis(2, at=seq(floor(min(y)), 1, by=0.1), labels=FALSE, tcl=-0.2)")

	if draw_secondary_grid_lines:
		r("abline(h = seq(floor(min(y)), 1, by=0.1), v = seq(floor(min(x)), length(x), by= 1), col='gray')")
	else:
		r("grid(lty = 1, col = 'gray')")

	if (standardize=='true'):
		r("barplot(xz, xlim=c(0, " + max_break_str + "), ylim = c(0, " + str(h_yaxis) + "), yaxp=c(0, " + str(h_yaxis) + ", 1), space = 0, col='darkgray', width = " + break_unit_str + ", xlab='Location (Mb)', ylab='Normalized frequency of pure parental alleles ', main='LG " + chr + " (Variant Discovery Mapping)')")
		r("barplot(yz, space = 0, add=TRUE, width = " + half_break_unit_str + ", col='green2')")	
		r("axis(1, hadj = 1, at=seq(0, " + max_break_str + ", by= " + break_unit_str + "), labels=FALSE, tcl=-0.5)")
		r("axis(1, at=seq(0, " + max_break_str + ", by= " + break_penta_unit_str + "), labels=TRUE, tcl=-0.5)")
		r("axis(1, at=seq(0, " + max_break_str + ", by= " + half_break_unit_str + "), labels=FALSE, tcl=-0.25)")
	elif (standardize=='false'):
		r("barplot(xz, ylim = c(0, " + str(h_yaxis) + "), yaxp=c(0, " + str(h_yaxis) + ", 1), space = 0, col='darkgray', width = 1, xlab='Location (Mb)', ylab='Normalized frequency of pure parental alleles ', main='LG " + chr + " (Variant Discovery Mapping)')")	
		r("barplot(yz, space = 0, add=TRUE, width = 0.5, col='green2')")	
		r("axis(1, at=seq(0, as.integer( ' " + str(breaks) + " '), by= " + break_unit_str + "), labels=FALSE, tcl=-0.5)")
		r("axis(1, at=seq(0, as.integer( ' " + str(breaks) + " '), by= " + break_penta_unit_str + "), labels=TRUE, tcl=-0.5)")
		r("axis(1, at=seq(0, as.integer( ' " + str(breaks) + " '), by= " + half_break_unit_str + "), labels=FALSE, tcl=-0.25)")


		
def calculate_normalized_histogram_bins_per_xbase(vcf_info = None, xbase = 1000000, normalize_bins = None):
	normalized_histogram_bins_per_xbase = {}

	ref_snp_count_per_xbase = get_ref_snp_count_per_xbase(vcf_info = vcf_info, xbase = xbase)
	mean_one_ratio_snp_count_per_chromosome = get_mean_one_ratio_snp_count_per_chromosome(vcf_info = vcf_info, xbase = xbase)
	one_ratio_snp_count_per_xbase = get_one_ratio_snp_count_per_xbase(vcf_info = vcf_info, xbase = xbase)

	for location in ref_snp_count_per_xbase:
		chromosome = location.split(':')[0]
		mean_one_ratio_snp_count = mean_one_ratio_snp_count_per_chromosome[chromosome]
		ref_snp_count = ref_snp_count_per_xbase[location]

		one_ratio_snp_count = 0
		if location in one_ratio_snp_count_per_xbase:
			one_ratio_snp_count = one_ratio_snp_count_per_xbase[location]	

		if normalize_bins == 'true':							
			if one_ratio_snp_count == 0:
				normalized_histogram_bins_per_xbase[location] = 0
			elif one_ratio_snp_count == ref_snp_count:
				normalized_histogram_bins_per_xbase[location] = (Decimal(one_ratio_snp_count**2) * Decimal(mean_one_ratio_snp_count))					
			else:
				normalized_histogram_bins_per_xbase[location] = (Decimal(one_ratio_snp_count**2) / (Decimal(ref_snp_count)-Decimal(one_ratio_snp_count))) * Decimal(mean_one_ratio_snp_count)					
		else:
			normalized_histogram_bins_per_xbase[location] = one_ratio_snp_count

	return normalized_histogram_bins_per_xbase


def get_ref_snp_count_per_xbase(vcf_info = None, xbase = 1000000):
	ref_snps_per_xbase = {}

	for location in vcf_info:
		location_info = location.split(':')

		chromosome = location_info[0].upper()
		chromosome = re.sub("CHROMOSOME_", "", chromosome, flags = re.IGNORECASE)
		chromosome = re.sub("chr", "", chromosome, flags = re.IGNORECASE)

		position = location_info[1]
		xbase_position = (int(position) / xbase) + 1

		location = chromosome + ":" + str(xbase_position)
		if location in ref_snps_per_xbase:
			ref_snps_per_xbase[location] = ref_snps_per_xbase[location] + 1
		else:
			ref_snps_per_xbase[location] = 1

	return ref_snps_per_xbase



def get_mean_one_ratio_snp_count_per_chromosome(vcf_info, xbase = 1000000):
	sample_snp_count_per_xbase = {}
	
	for location in vcf_info:
		alt_allele_count, ref_allele_count, read_depth, ratio = vcf_info[location]
			
		location_info = location.split(':')
		chromosome = location_info[0]
		position = location_info[1]
		xbase_position = (int(position) / xbase) + 1
		xbase_location = chromosome + ":" + str(xbase_position)
		
		if ratio == 1:
			if xbase_location in sample_snp_count_per_xbase:
				sample_snp_count_per_xbase[xbase_location] = sample_snp_count_per_xbase[xbase_location] + 1
			else:
				sample_snp_count_per_xbase[xbase_location] = 1
		
		elif ratio != 1 and xbase_location not in sample_snp_count_per_xbase:
			sample_snp_count_per_xbase[xbase_location] = 0

	mean_one_ratio_snp_count_per_chromosome = {}
	for location in sample_snp_count_per_xbase:
		chromosome = location.split(':')[0]
		sample_count = sample_snp_count_per_xbase[location]
		if chromosome in mean_one_ratio_snp_count_per_chromosome:
			mean_one_ratio_snp_count_per_chromosome[chromosome].append(sample_count)
		else:
			mean_one_ratio_snp_count_per_chromosome[chromosome] = [sample_count]

	for chromosome in mean_one_ratio_snp_count_per_chromosome:
		summa = sum(mean_one_ratio_snp_count_per_chromosome[chromosome])
		count = len(mean_one_ratio_snp_count_per_chromosome[chromosome])

		mean_one_ratio_snp_count_per_chromosome[chromosome] = Decimal(summa) / Decimal(count)

	return mean_one_ratio_snp_count_per_chromosome


def get_one_ratio_snp_count_per_xbase(vcf_info = None, xbase = 1000000):
	one_ratio_snp_count_per_xbase = {}

	for location in vcf_info:
		alt_allele_count, ref_allele_count, read_depth, ratio = vcf_info[location]
			
		location_info = location.split(':')
		chromosome = location_info[0]
		position = location_info[1]
		xbase_position = (int(position) / xbase) + 1
		xbase_location = chromosome + ":" + str(xbase_position)
		
		if ratio == 1:
			if xbase_location in one_ratio_snp_count_per_xbase:
				one_ratio_snp_count_per_xbase[xbase_location] = one_ratio_snp_count_per_xbase[xbase_location] + 1
			else:
				one_ratio_snp_count_per_xbase[xbase_location] = 1

		elif ratio != 1 and xbase_location not in one_ratio_snp_count_per_xbase:
			one_ratio_snp_count_per_xbase[xbase_location] = 0

	return one_ratio_snp_count_per_xbase


def parse_vcf(sample_vcf = None):
	i_file = open(sample_vcf, 'rU')
	reader = csv.reader(i_file, delimiter = '\t', quoting = csv.QUOTE_NONE)	

	skip_headers(reader = reader, i_file = i_file)
	vcf_info = {}

	for row in reader:
		chromosome = row[0].upper()
		chromosome = re.sub("CHROMOSOME_", "", chromosome, flags = re.IGNORECASE)
		chromosome = re.sub("chr", "", chromosome, flags = re.IGNORECASE)


		if chromosome != 'MTDNA':
			position = row[1]
			#ref_allele = row[2]
			#read_depth = row[3]
			#read_bases = row[4]

			vcf_format_info = row[8].split(":")
			vcf_allele_freq_data = row[9] 
			
			read_depth_data_index = vcf_format_info.index("DP")
			read_depth = vcf_allele_freq_data.split(":")[read_depth_data_index]

			ref_and_alt_counts_data_index = vcf_format_info.index("AD")
			ref_and_alt_counts = vcf_allele_freq_data.split(":")[ref_and_alt_counts_data_index]	
			ref_allele_count = ref_and_alt_counts.split(",")[0]
			alt_allele_count = ref_and_alt_counts.split(",")[1]

			location = chromosome + ":" + position
	
			if Decimal(read_depth!=0):		
				getcontext().prec = 6	
				ratio = Decimal(alt_allele_count) / Decimal(read_depth)
	
				vcf_info[location] = (alt_allele_count, ref_allele_count, read_depth, ratio)
		
				#debug line
				#print chromosome, position, read_depth, ref_allele_count, alt_allele_count, ratio, id

	i_file.close()

	return vcf_info

def parse_read_bases(read_bases = None, alt_allele = None):
	read_bases = re.sub('\$', '', read_bases)
	read_bases = re.sub('\^[^\s]', '', read_bases)

	ref_allele_matches = re.findall("\.|\,", read_bases)
	ref_allele_count = len(ref_allele_matches)

	alt_allele_matches = re.findall(alt_allele, read_bases, flags = re.IGNORECASE)
	alt_allele_count = len(alt_allele_matches)

	#debug line
	#print read_bases, alt_allele, alt_allele_count, ref_allele_count

	return ref_allele_count, alt_allele_count

if __name__ == "__main__":
	main()
