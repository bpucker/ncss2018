### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.1 ###

__usage__ = """
					python compare_variation_rates.py
					--vcf <FULL_PATH_TO_VCF_FILE>
					--ncss <FULL_PATH_TO_SPLICE_SITE_FILE>
					--gff <FULL_PATH_TO_ANNOTATION_FILE>
					--fasta <FULL_PATH_TO_GENOME_FILE>
					"""

import sys, re

# --- end of imports --- #

def load_chr_per_gene( gff_file ):
	"""! @brief load chromosomes of all genes """
	
	chr_per_gene = {}
	with open( gff_file, "r" ) as f:
		line = f.readline()
		while line:
			if line[0] != '#':
				parts = line.strip().split('\t')
				if parts[2] == "mRNA":
					try:
						ID = re.findall( "rna\d+", parts[-1] )[0]
						try:
							chr_per_gene[ ID ]
						except KeyError:
							chr_per_gene.update( { ID: parts[0] } )
					except IndexError:
						print line
			line = f.readline()
	return chr_per_gene


def load_splice_site_pos( ncss_file, chr_per_gene ):
	"""! @brief load all splice site positions """
	
	css_splice_site_pos = {}
	major_ncss_splice_site_pos = {}
	minor_ncss_splice_site_pos = {}
	
	css_total_len = 0
	major_ncss_total_len = 0
	minor_ncss_total_len = 0
	
	with open( ncss_file, "r" ) as f:
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			if len( parts ) > 3:
				
				subparts = parts[4].split(',')
				chromosome = chr_per_gene[ parts[0] ]
				if parts[1:3] == [ "GT", "AG" ]:
					css_total_len += 4
					for each in subparts:
						start = int( each.split('_')[0] )
						css_splice_site_pos.update( { chromosome + "_%_" + str( start+1 ): None } )
						css_splice_site_pos.update( { chromosome + "_%_" + str( start+2 ): None } )
				elif parts[1:3] == [ "GC", "AG" ] or parts[1:3] == [ "AT", "AC" ]:
					major_ncss_total_len += 4
					for each in subparts:
						start = int( each.split('_')[0] )
						major_ncss_splice_site_pos.update( { chromosome + "_%_" + str( start+1 ): None } )
						major_ncss_splice_site_pos.update( { chromosome + "_%_" + str( start+2 ): None } )
				else:
					minor_ncss_total_len += 4
					for each in subparts:
						start = int( each.split('_')[0] )
						minor_ncss_splice_site_pos.update( { chromosome + "_%_" + str( start+1 ): None } )
						minor_ncss_splice_site_pos.update( { chromosome + "_%_" + str( start+2 ): None } )
			line = f.readline()
	return css_splice_site_pos, major_ncss_splice_site_pos, minor_ncss_splice_site_pos, css_total_len, major_ncss_total_len, minor_ncss_total_len


def load_variants( vcf_file ):
	"""! @brief load all variants from given VCF """
	
	variants = {}
	
	with open( vcf_file, "r" ) as f:
		line = f.readline()
		while line:
			if line[0] != '#':
				parts = line.strip().split('\t')
				variants.update( { parts[0] + "_%_" + parts[1]: None } )
			line = f.readline()
	return variants


def get_genome_size( fasta_file ):
	"""! @brief get genome size """
	
	genome_size = 0
	with open( fasta_file, "r" ) as f:
		line = f.readline()
		while line:
			if line[0] != '>':
				genome_size += len( line.strip() )
			line = f.readline()
	return genome_size


def get_matches( variants, splice_site_pos, total_len ):
	"""! @brief calcualte matches """
	
	match = 0
	mismatch = 0
	for variant in variants.keys():
		try:
			splice_site_pos[ variant ]
			match += 1
		except KeyError:
			mismatch += 1
	return ( 100.0 * match ) / total_len


def main( arguments ):
	"""! @brief compare general variation rates against variation rates in splice sites """
	
	vcf_file = arguments[ arguments.index('--vcf')+1 ]		#VCF file containing variants of the same species
	ncss_file = arguments[ arguments.index('--ncss')+1 ]	#splice sites file
	gff_file = arguments[ arguments.index('--gff')+1 ]	#GFF3 file
	fasta_file = arguments[ arguments.index('--fasta')+1 ]	#FASTA file of genome
	
	
	# --- check the mutational load of splice sites --- #
	chr_per_gene = load_chr_per_gene( gff_file )
	css_splice_site_pos, major_ncss_splice_site_pos, minor_ncss_splice_site_pos, css_total_len, major_ncss_total_len, minor_ncss_total_len = load_splice_site_pos( ncss_file, chr_per_gene )
	variants = load_variants( vcf_file )
	genome_size = get_genome_size( fasta_file )
	
	print len( css_splice_site_pos.keys() )
	print len( major_ncss_splice_site_pos.keys() )
	print len( minor_ncss_splice_site_pos.keys() )
		
	print "percentage of canonical splice site positions with variants: " + str( get_matches( variants, css_splice_site_pos, css_total_len ) ) + "%"
	print "percentage of major ncss positions with variants: " + str( get_matches( variants, major_ncss_splice_site_pos, major_ncss_total_len ) ) + "%"
	print "percentage of minor ncss positions with variants: " + str( get_matches( variants, minor_ncss_splice_site_pos, minor_ncss_total_len ) ) + "%"


if __name__ == '__main__':
	if '--vcf' in sys.argv and '--ncss' in sys.argv and '--gff' in sys.argv and '--fasta' in sys.argv:
		main( sys.argv )
	else:
		sys.exit( __usage__ )
