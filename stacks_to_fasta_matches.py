#!/usr/bin/env python
# Python 2.7
# stacks_to_fasta_matches.py v1.1
# 5/21/13
# Ryan Waples



# Command line outline
# > python   stacks_to_fasta_matches.py   stacks_path   catalog_name individual_name  subset_file   FASTA_file

# Usage example
# > python ./stacks_to_fasta_matches.py /stacks/data/ batch_3.catalog   PHOOD_2131  ./subset.txt  ./output.fa


# Inputs
	# A set of stacks catalog and individual files (*.tags.tsv, *.snp.tsv, *.alleles.tsv)
		# specified by a [path] + [catalog_name] + [individual_name]
	# A 'whitelist' of tags.
		# specefied by a [filename]
		# tags_subset_file  format - one tag_id per line.
			# 1
			# 43
			# 12395
			# 6549
			# etc..

# Outputs	
	# An FASTA file	
		# For each white listed tag:
			# Each haplotype present in the catalog will generate one FASTA sequence entry.
		# FASTA header: 
			# >[TAG#]_[ALLELE]_[SNP1_pos]-[SNP2_pos]-[SNP3_pos]-[...]

import sys

def stacks_to_fasta_matches(stacks_path, catalog_name, individual_name, subset_file, out_file):
	"""
	
	Takes a whitelist of catalog tags and prints each observed haplotype in the supplied indiviudal
	into a sequence to a FASTA file.
	
	Function also returns a dictionary of FASTA entries.
	
	"""
	print "\n"
	if not stacks_path.endswith('/'):
		stacks_path = stacks_path + "/"
	print ("\n".join(['stacks_path ->\t' + stacks_path, 
						'catalog_name ->\t' + catalog_name, 
						"individual_name ->\t" + individual_name,
						'subset_file ->\t' + subset_file, 
						'out_file ->\t'+ out_file]))
	print "\n"	
	
	#used variables
	tags_to_keep 	= set()		# strings
	seq_of_tag		= dict()	# strings
	snp_pos_of_tag	= dict()	# lists of ints
	alleles_of_tag	= dict()	# lists of strings
	ind_ids_of_cat_id = dict()	# set of strings
	FASTA_dict		= dict()	# strings
	


	with open(subset_file, "r") as SUBSET:
		for line in SUBSET:
			tag = line.strip()
			tags_to_keep.add(tag)

	with open("{0}{1}{2}".format(stacks_path, catalog_name, ".tags.tsv" ), "r") as TAGS:
		for line in TAGS:
			columns = line.strip().split("\t")
			id = columns[2] 
			seq = columns[9]
			seq_type = columns[6]
			if (id in tags_to_keep and seq_type == 'consensus'):
				seq_of_tag[id] = seq

	with open("{0}{1}{2}".format(stacks_path, catalog_name, ".snps.tsv" ), "r") as SNPS:	
		for line in SNPS:
			columns = line.strip().split("\t")
			id = columns[2]
			pos = columns[3]
			if id in tags_to_keep and id in seq_of_tag:
				snp_pos_of_tag.setdefault(id, [])
				snp_pos_of_tag[id].append(int(pos))
			

	with open("{0}{1}{2}".format(stacks_path, individual_name, ".matches.tsv" ), "r") as MATCHES:
		for line in MATCHES:
			columns = line.strip().split("\t")
			cat_id = columns[2]
			ind_id = columns[4]
			allele = columns[5]
			if (cat_id in tags_to_keep and cat_id in seq_of_tag):
				alleles_of_tag.setdefault(cat_id, [])
				alleles_of_tag[cat_id].append(allele)
				ind_ids_of_cat_id.setdefault(cat_id, set())
				ind_ids_of_cat_id[cat_id].add(ind_id)
							
	with open(out_file, "w") as FASTA:
		for cat_id in sorted(alleles_of_tag.keys(), key=int):
			base_seq = seq_of_tag[cat_id]
			haplotype = base_seq[:]
			for allele in alleles_of_tag[cat_id]:
				snp_index = 0
				try:
					for snp in sorted(snp_pos_of_tag[cat_id]):
						haplotype = haplotype[:snp] + allele[snp_index] + haplotype[snp+1:]
						snp_index += 1
				except KeyError: # happens on 'consensus' genotypes
					haplotype = haplotype
				try:
					snp_positions = "_".join([str(x) for x in snp_pos_of_tag[cat_id]] )
				except KeyError:
					snp_positions = "none"
				ind_id = "_".join(ind_ids_of_cat_id[cat_id])
				header = ">cat:{}|ind:{}|allele:{}|pos:{}".format(cat_id, ind_id, allele, snp_positions)
				FASTA_dict[header] = haplotype
				header_line = header + "\n"
				seq_line = haplotype + "\n"
				FASTA.write(header_line)
				FASTA.write(seq_line)

	return FASTA_dict	


def main(argv=None):
	if argv is None:
		argv = sys.argv[1:]

	stacks_to_fasta_matches(*argv)
	
if __name__ == "__main__":
	main(sys.argv[1:])
			
