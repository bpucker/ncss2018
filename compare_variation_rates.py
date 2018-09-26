### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.1 ###

__usage__ = """
					python compare_variation_rates.py
					--vcf <FULL_PATH_TO_VCF_FILE>
					--ncss <FULL_PATH_TO_SPLICE_SITE_FILE>
					"""

from scipy import stats
import sys

# --- end of imports --- #


def load_ncss( ncss_file ):
	"""! @brief load all minor ncss """
	
	minor_ncss = []
	total_splice_site_length = 0
	with open( ncss_file, "r" ) as f:
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			if len( parts ) > 3:
				total_splice_site_length += 4
				if parts[3] == "ncss":
					if parts[1] == "GC" and parts[2] == "AG":
						pass
					elif parts[1] == "AT" and parts[2] == "AC":
						pass
					else:
						minor_ncss.append( parts[1]+parts[2] )
			line = f.readline()
	return minor_ncss, total_splice_site_length


def load_variation_rate_from_vcf( vcf_file ):
	"""! @brief load variation rate from given VCF file
	@warning only homozygous variants against the reference are considered
	"""
	
	# --- generate all keys --- #
	variation_counts = {}
	for nt1 in [ "A", "C", "G", "T" ]:
		for nt2 in [ "A", "C", "G", "T" ]:
			if nt1 != nt2:
				variation_counts.update( { nt1+nt2: 0 } )
	
	# --- process file --- #
	with open( vcf_file, "r" ) as f:
		line = f.readline()
		while line:
			if line[0] != '#':
				parts = line.strip().split('\t')
				if len( parts[3] ) == 1:
					if len( parts[4] ) == 1:
						if parts[3] != parts[4]:
							if parts[3] != "N" and parts[4] != "N":
								variation_counts[ parts[3]+parts[4] ] += 1
			line = f.readline()
	
	return variation_counts


def get_variants( query, ref ):
	"""! @brief get variants from ref to query """
	
	variants = []
	for idx, nt in enumerate( query ):
		if nt != ref[ idx ]:
			variants.append( ref[ idx ]+nt )
	return variants


def calculate_variation_rate_for_ncss( minor_ncss ):
	"""! @brief calcualte variation ratios for minor ncss """
	
	# --- generate all keys --- #
	variation_counts = {}
	for nt1 in [ "A", "C", "G", "T" ]:
		for nt2 in [ "A", "C", "G", "T" ]:
			if nt1 != nt2:
				variation_counts.update( { nt1+nt2: 0 } )
	
	# --- investigate data --- #
	for ncss in minor_ncss:
		gtag = get_variants( ncss, "GTAG" )
		gcag = get_variants( ncss, "GCAG" )
		atac = get_variants( ncss, "ATAC" )
		
		if gtag >= gcag and gtag >= atac:
			for each in gtag:
				if not "N" in each:
					variation_counts[ each ] += 1
		elif gcag >= gtag and gcag >= atac:
			for each in gcag:
				if not "N" in each:
					variation_counts[ each ] += 1
		else:
			for each in atac:
				if not "N" in each:
					variation_counts[ each ] += 1
	
	return variation_counts


def main( arguments ):
	"""! @brief compare general variation rates against variation rates in splice sites """
	
	vcf_file = arguments[ arguments.index('--vcf')+1 ]		#VCF file containing variants of the same species
	ncss_file = arguments[ arguments.index('--ncss')+1 ]	#splice sites file
	
	minor_ncss, total_splice_site_length = load_ncss( ncss_file )
	variation_rate = load_variation_rate_from_vcf( vcf_file )
	minor_ncss_variation_rate = calculate_variation_rate_for_ncss( minor_ncss )
	
	ncss_values = []
	global_values = []
	for key in sorted( minor_ncss_variation_rate.keys() ):
		ncss_values.append( minor_ncss_variation_rate[ key ] )
		global_values.append( variation_rate[ key ] )
	
	print sorted( minor_ncss_variation_rate.keys() )
	print ncss_values
	print global_values
	
	print "comparison of splice site diversity to overall variation pattern:"
	print stats.chisquare( ncss_values,  f_exp=global_values )


if __name__ == '__main__':
	if '--vcf' in sys.argv and '--ncss' in sys.argv:
		main( sys.argv )
	else:
		sys.exit( __usage__ )
