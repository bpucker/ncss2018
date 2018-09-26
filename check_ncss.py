### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.1 ###

### MEMORY SAVING VERSION ###

__usage__ = """
	python check_ncss.py
	--cov <FULL_PATH_TO_INPUT_COVERAGE_FILE>
	--doc <FULL_PATH_TO_INPUT_DOCUMENTATION_FILE>
	--out <FULL_PATH_TO_OUTPUT_DIRECTORY>
	
	bug reports and feature requests: bpucker@cebitec.uni-bielefeld.de
					"""

import glob, re, sys, os, shelve
import matplotlib.pyplot as plt
import numpy as np
from operator import itemgetter

# --- end of imports --- #

def load_pos_from_gff( gff_file ):
	"""! @brief load positions from GFF3 file """
	
	pos = {}
	with open( gff_file, "r" ) as f:
		line = f.readline()
		while line:
			if line[0] != '#':
				parts = line.strip().split('\t')
				if parts[2] in [ "gene", "mRNA" ]:
					try:
						ID = re.findall( "rna\d+", parts[-1] )[0]
						parent = re.findall( "gene\d+", parts[-1] )[0]
					except IndexError:
						ID = re.findall( "gene\d+", parts[-1] )[0]
						parent = re.findall( "gene\d+", parts[-1] )[0]
					if parts[6]:
						orientation = True
					else:
						orientation = False
					try:
						pos.update( { ID: { 'id': ID, 'parent': parent, 'chr': parts[0], 'start': parts[3], 'end': parts[4], 'orientation': orientation } } )
					except:
						pass
			line = f.readline()
	return pos


def load_all_ncss_regions( docf, transcript_gene_pos ):
	"""! @brief load all ncss candidate regions from given documentation file """
	
	ncss_regions = []
	
	with open( docf, "r" ) as f:
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			if len( parts ) > 3 and "css" in line:	#could be edited back to 'ncss' to focus on those only
				info = transcript_gene_pos[ parts[0] ]
				ex3_pos, in5_pos, in3_pos, ex5_pos = parts[4].split(',')
				if info['orientation']:
					ncss_regions.append( { 	'id': info['parent'],
															'chr': info['chr'],
															'e3': int( ex3_pos.split('_')[1] )-1,
															'i5': int( in5_pos.split('_')[0] ),
															'i3': int( in3_pos.split('_')[1] )-1,
															'e5': int( ex5_pos.split('_')[0] ),
															'5prime': parts[1],
															'3prime': parts[2]
														} )
				else:
					ncss_regions.append( { 	'id': info['parent'],
															'chr': info['chr'],
															'e3': int( ex3_pos.split('_')[0] ),
															'i5': int( in5_pos.split('_')[1] )-1,
															'i3': int( in3_pos.split('_')[0] ),
															'e5': int( ex5_pos.split('_')[1] )-1,
															'5prime': parts[1],
															'3prime': parts[2]
														} )
			line = f.readline()
	return sorted( ncss_regions, key=itemgetter( 'chr' ) )


def load_coverage( cov_file, chromosome ):
	"""! @brief load coverage from given file """
	
	with open( cov_file, "r" ) as f:
		line = f.readline()
		cov = []
		while line:
			parts = line.strip().split('\t')
			if parts[0] == chromosome:
				cov.append( float( parts[2] ) )
			line = f.readline()
	return cov


def analyze_cov_at_ncss( cov_file, ncss_regions, result_file ):
	"""! @brief analyze coverage around ncss """
	
	with open( result_file, "w", 0 ) as out:
		out.write( "GeneID\tExon3prime\tIntron5prime\tIntron3prime\tExon5prime\t5prime_splice_site\t3pime_splice_site\n" )
		current_chr = ""
		for ncss in ncss_regions:
			if current_chr != ncss['chr']:
				chr_cov = load_coverage( cov_file, ncss['chr'] )
				current_chr = "" + ncss['chr']
			e3 = chr_cov[ ncss['e3'] ]
			i5 = chr_cov[ ncss['i5'] ]
			i3 = chr_cov[ ncss['i3'] ]
			e5 = chr_cov[ ncss['e5'] ]
			out.write( "\t".join( map( str, [ ncss['id'], e3, i5, i3, e5, ncss['5prime'], ncss['3prime'] ] ) ) + '\n' )


def get_cov_at_ncss( doc_file, cov_file, output_dir ):
	"""! @brief get coverage at all putative ncss """
	
	gff_file = doc_file.replace( ".txt", ".gff" )
	transcript_gene_pos = load_pos_from_gff( gff_file )
	ncss_regions = load_all_ncss_regions( doc_file, transcript_gene_pos )
	print "number of ncss regions: " + str( len( ncss_regions ) )
	result_file = output_dir + "splice_site_coverage_check.txt"
	analyze_cov_at_ncss( cov_file, ncss_regions, result_file )
	return result_file


def load_all_splice_sites( input_file ):
	"""! @brief load all splice sits from given input file """
	
	splice_sites = []
	with open( input_file, "r" ) as f:
		f.readline()	#header
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			splice_sites.append( { 'id': parts[0], 'e3': float( parts[1] ), 'i5': float( parts[2] ), 'i3': float( parts[3] ), 'e5': float( parts[4] ), '5prime': parts[5], '3prime': parts[6] } )
			line = f.readline()
	return splice_sites


def analyze_splice_sites( splice_sites, output_dir, min_cov_cutoff ):
	"""! @brief analyze splice sites by group """
	
	GTAG_5prime_ratio = []
	GTAG_3prime_ratio = []
	GTAG_not_expressed = 0
	
	GCAG_5prime_ratio = []
	GCAG_3prime_ratio = []
	GCAG_not_expressed = 0
	
	ATAC_5prime_ratio = []
	ATAC_3prime_ratio = []
	ATAC_not_expressed = 0
	
	others_5prime_ratio = []
	others_3prime_ratio = []
	others_not_expressed = 0
	
	for site in splice_sites:
		if site['5prime'] == "GT" and site['3prime'] == "AG":
			if site['e3'] > 0 and site['e5'] > 0:
				GTAG_5prime_ratio.append( ( site['e3']-site['i5'] ) / site['e3'] )
				GTAG_3prime_ratio.append( ( site['e5']-site['i3'] ) / site['e5']  )
			else:
				GTAG_not_expressed += 1
		elif site['5prime'] == "GC" and site['3prime'] == "AG":
			if site['e3'] > 0 and site['e5'] > 0:
				GCAG_5prime_ratio.append( ( site['e3']-site['i5'] ) / site['e3'] )
				GCAG_3prime_ratio.append( ( site['e5']-site['i3'] ) / site['e5']  )
			else:
				GCAG_not_expressed += 1
		elif site['5prime'] == "AT" and site['3prime'] == "AC":
			if site['e3'] > 0 and site['e5'] > 0:
				ATAC_5prime_ratio.append( ( site['e3']-site['i5'] ) / site['e3'] )
				ATAC_3prime_ratio.append( ( site['e5']-site['i3'] ) / site['e5']  )
			else:
				ATAC_not_expressed += 1
		else:
			if site['e3'] > 0 and site['e5'] > 0:
				others_5prime_ratio.append( ( site['e3']-site['i5'] ) / site['e3'] )
				others_3prime_ratio.append( ( site['e5']-site['i3'] ) / site['e5']  )
			else:
				others_not_expressed += 1
	
	# --- construct ratio distribution plots --- #
	fig_5prime_out = output_dir + "5prime_ratios.png"
	fig, ax = plt.subplots()
	#ax.hist( GTAG_5prime_ratio, bins=10000, color="green", label="GT-AG" )
	ax.hist( GCAG_5prime_ratio, bins=1000, color="blue", label="GC-AG" )
	ax.hist( ATAC_5prime_ratio, bins=1000, color="purple", label="AT-AC" )
	ax.hist( others_5prime_ratio, bins=1000, color="red", label="others" )
	
	ax.set_ylim( 0, 100 )
	ax.set_xlim( 0, 1 )
	
	ax.set_xlabel("coverage ratio over 5prime intron border")
	ax.set_ylabel("counts")
	
	fig.savefig( fig_5prime_out, dpi=300 )


def get_supported_ncss( splice_sites, supported_ncss_file, cutoff ):
	"""! @brief get supported ncss """
	
	supported_ncss = []
	
	with open( supported_ncss_file, "w" ) as out:
		for site in splice_sites:
			if site['5prime'] != "GT": 
				if site['e3'] > 0 and site['e5'] > 0:
					if ( ( site['e3']-site['i5'] ) / site['e3'] ) > cutoff and ( ( site['e5']-site['i3'] ) / site['e5'] ) > cutoff:
						out.write( "\t".join( [ site['id'], site['5prime'], site['3prime'] ] ) + '\n' )
						supported_ncss.append( site['5prime']+".."+site['3prime'] )
			elif site['3prime'] != "AG":
				if site['e3'] > 0 and site['e5'] > 0:
					if ( ( site['e3']-site['i5'] ) / site['e3'] ) > cutoff and ( ( site['e5']-site['i3'] ) / site['e5'] ) > cutoff:
						out.write( "\t".join( [ site['id'], site['5prime'], site['3prime'] ] ) + '\n' )
						supported_ncss.append( site['5prime']+".."+site['3prime'] )
			else:
				if site['e3'] > 0 and site['e5'] > 0:
					if ( ( site['e3']-site['i5'] ) / site['e3'] ) > cutoff and ( ( site['e5']-site['i3'] ) / site['e5'] ) > cutoff:
						out.write( "\t".join( [ site['id'], site['5prime'], site['3prime'] ] ) + '\n' )
						supported_ncss.append( site['5prime']+".."+site['3prime'] )
	print "number of supported ncss: " + str( len( supported_ncss ) )
	overview_file = supported_ncss_file + "_overview.txt"
	with open( overview_file, "w" ) as out:
		for splice_site in list( set( supported_ncss ) ):
			out.write( splice_site + "\t" + str( supported_ncss.count( splice_site ) ) + '\n' )


def investigate_ncss( input_file, output_dir, cutoff, min_cov_cutoff ):
	"""! @brief check coverage drop at predicted ncss """
	
	supported_ncss_file = output_dir + "supported_ncss.txt"
	
	splice_sites = load_all_splice_sites( input_file )
	#analyze_splice_sites( splice_sites, output_dir, min_cov_cutoff )
	
	get_supported_ncss( splice_sites, supported_ncss_file, cutoff )


def main( arguments ):
	"""! @brief get coverage at ncss and check all putative ncss if they are real """
	
	output_dir = arguments[ arguments.index('--out')+1 ]
	cov_file = arguments[ arguments.index('--cov')+1 ]
	doc_file = arguments[ arguments.index('--doc')+1 ]
	
	if output_dir[-1] != '/':
		output_dir += "/"
	if not os.path.exists( output_dir ):
		os.makedirs( output_dir )
	
	splice_site_support_cutoff = 0.2	#difference in coverage between exon and intron to call valid splice site
	min_cov_cutoff = 3
	
	result_file = get_cov_at_ncss( doc_file, cov_file, output_dir )
	investigate_ncss( result_file, output_dir, splice_site_support_cutoff, min_cov_cutoff )


if __name__ == "__main__":
	if '--out' in sys.argv and '--cov' in sys.argv and '--doc' in sys.argv:
		main( sys.argv )
	else:
		sys.exit( __usage__ )
