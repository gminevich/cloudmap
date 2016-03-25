#!/usr/bin/python

import sys
import optparse
import csv
import operator
import pprint

def main():
	parser = optparse.OptionParser()
	parser.add_option("-i", dest = "input_files", action = "callback", callback = split_input_files)
	parser.add_option("-n", dest = "sample_names", action = "callback", callback = split_sample_names)
	parser.add_option("-s", dest = "summary_output_file", action = "store", type = "string")
	parser.add_option("-o", dest = "data_output_file", action = "store", type = "string")
	(options, args) = parser.parse_args()

	#get the counts
	gene_count = {}
	for input_file in options.input_files:
		gene_count = summarize_counts(input_file = input_file, gene_count = gene_count)
	
	#delete the single counts
	for k, v in gene_count.items():
		if v == 1:
			del gene_count[k]

	# we do two passes because we don't know the order yet and don't want to hold onto everything in memory

	file_sample_dictionary = create_file_sample_dictionary(input_files = options.input_files, sample_names = options.sample_names)

	#get snpeff lines
	gene_data = {}
	all_headers = []
	for input_file in options.input_files:
		headers = []
		headers, gene_data = grab_snpeff_data(input_file = input_file, gene_data = gene_data, gene_count = gene_count, file_sample_dictionary = file_sample_dictionary)
		
		if len(all_headers) == 0 :
			all_headers.extend(headers)

	sorted_gene_count_list = sorted(gene_count.iteritems(), key = operator.itemgetter(1))
	sorted_gene_count_list.reverse()

	print_summary(gene_count_list = sorted_gene_count_list, summary_output_file = options.summary_output_file)
	print_sorted_snpeff_data(gene_data = gene_data, gene_count_list = sorted_gene_count_list, data_output_file = options.data_output_file, all_headers = all_headers)

def create_file_sample_dictionary(input_files = None, sample_names = None):
	file_sample_dictionary = {}

	for i in range(0, len(input_files)):
		if i < len(sample_names):
			file_sample_dictionary[input_files[i]] = sample_names[i]
		else:
			file_sample_dictionary[input_files[i]] = ""

	return file_sample_dictionary
		
def split_input_files(option, opt_str, value, parser):
	input_files = [] 
	collect = True  
	for i in range(0, len(parser.rargs)):
		if "-s" == parser.rargs[i] or "-o" == parser.rargs[i] or "-n" == parser.rargs[i]:
			collect = False
		elif collect == True:
			input_files.append(parser.rargs[i])

	setattr(parser.values, option.dest, input_files)

def split_sample_names(option, opt_str, value, parser):
	input_files = [] 
	collect = True  
	for i in range(0, len(parser.rargs)):
		if "-s" == parser.rargs[i] or "-o" == parser.rargs[i] or "-i" == parser.rargs[i]:
			collect = False
		elif collect == True:
			input_files.append(parser.rargs[i])

	setattr(parser.values, option.dest, input_files)


def print_summary(gene_count_list = None, summary_output_file = ""):
	f = open(summary_output_file, 'w')
	f.write('\n'.join('%s\t%s' % g for g in gene_count_list))
	f.close()

def print_sorted_snpeff_data(gene_data = {}, gene_count_list = None, data_output_file = "", file_sample_dictionary = None, all_headers = None):
	f = open(data_output_file, 'w')

	for header in all_headers:
		f.write('Sample Name')
		f.write('\t')
		f.write(('\t').join(header))
		f.write('\n')

	for i in range(0, len(gene_count_list)):
		gene_id = gene_count_list[i][0]
		f.write(gene_data[gene_id] + "\n")

	f.close()

def grab_snpeff_data(input_file = "", gene_data = {}, gene_count = {}, file_sample_dictionary = None):
	i_file = open(input_file, 'rU')
	reader = csv.reader(i_file, delimiter = '\t')
	headers = skip_and_collect_headers(reader = reader, i_file = i_file)

	sample_name = file_sample_dictionary[input_file]

	for row in reader:
		gene_id = row[9].upper()
		if gene_id in gene_count:
			if gene_id in gene_data:
				gene_data[gene_id] = gene_data[gene_id] + '\n' + sample_name + '\t' + '\t'.join(row).rstrip()
			else:
				gene_data[gene_id] = sample_name + '\t' + '\t'.join(row).rstrip()

	return headers, gene_data

def summarize_counts(input_file= "", gene_count = {}):
	i_file = open(input_file, 'rU')
	reader = csv.reader(i_file, delimiter = '\t')
	skip_and_collect_headers(reader = reader, i_file = i_file)

	counted_genes = []

	for row in reader:
		gene_id = row[9].upper()
		if gene_id not in counted_genes:
			if gene_id in gene_count and gene_id not in counted_genes:
				gene_count[gene_id] += 1
			else:
				gene_count[gene_id] = 1
			
			counted_genes.append(gene_id)

	return gene_count
	
def skip_and_collect_headers(reader = None, i_file = None):
	headers = []	

	# count headers
	comment = 0
	while True:
		row = reader.next()
		if len(row) == 0:
			comment = comment + 1
		elif row[0].startswith('#'):
			comment = comment + 1
			if comment == 3:
				headers.append(row)
		else:
			break
	
	# skip headers
	i_file.seek(0)
	for i in range(0, comment):
		reader.next()

	return headers

if __name__ == "__main__":
	main()
