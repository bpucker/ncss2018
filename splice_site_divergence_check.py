### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.1 ###

__usage__ = """
	python splice_site_divergence_check.py
	--data_dir <FULL_PATH_TO_DATA_FROM_NCBI>
	--output_dir <FULL_PATH_TO_OUTPUT_DIRECTORY>
	--species_file <FULL_PATH_TO_FILE_WITH_SPECIES_ORDER>
					"""

import glob, re, os, sys
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
from operator import itemgetter
from scipy import stats
from math import log as ln
from scipy import stats

# --- end of imports --- #


def load_sequences( multiple_fasta_file ):
	"""! @brief load candidate gene IDs from file """
	
	sequences = {}
	
	with open( multiple_fasta_file ) as f:
		header = f.readline()[1:].strip().split(' ')[0]
		seq = ""
		line = f.readline()
		while line:
			if line[0] == '>':
					sequences.update( { header: seq } )
					header = line.strip()[1:].split(' ')[0]
					seq = ""
			else:
				seq += line.strip()
			line = f.readline()
		sequences.update( { header: seq } )
	return sequences


def load_results( filename ):
	"""! @brief load ncss genes and ncss distribution from doc file """
	
	ncss_genes = []
	ncss_distribution = {}
	rnas = []
	
	status = False
	with open( filename, "r" ) as f:
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			if len( parts ) >= 3:
				try:
					rnas.append( re.findall( "rna\d+", parts[0] )[0] )
				except IndexError:
					print line
			elif '..' in line:
				ncss_distribution.update( { parts[0]: int( parts[1] ) } )
			elif '#' in line:
				status = True
			elif status:
				ncss_genes.append( line.strip() )
			line = f.readline()
	return ncss_genes, ncss_distribution, rnas


def combine_all_results( results, output_file, species ):
	"""! @brief combine all results in one file """
	
	# --- get all possible splice site combinations --- #
	order = []
	all_splice_sites = {}
	for one in [ "A", "C", "G", "T" ]:
		for two in [ "A", "C", "G", "T" ]:
			for three in [ "A", "C", "G", "T" ]:
				for four in [ "A", "C", "G", "T" ]:			
					all_splice_sites.update( { one+two+"..."+three+four: [] } )
					order.append( one+two+"..."+three+four )
	
	# --- combine all data --- #
	#species = sorted( results.keys() )
	for spec in species:
		for splice_site in all_splice_sites.keys():
			try:
				all_splice_sites[ splice_site ].append( results[ spec ][ splice_site ] )
			except KeyError:
				all_splice_sites[ splice_site ].append( 0 )
	
	# --- writing output --- #
	with open( output_file, "w" ) as out:
		out.write( "splice_site\t" + "\t".join( species ) + '\n' )
		for splice_site in sorted( all_splice_sites.keys() ):
			out.write( "\t".join( [ splice_site ] + map( str, all_splice_sites[ splice_site ] ) ) + '\n' )
	return all_splice_sites, order


def find_outliers( all_splice_sites, species, outlier_file ):
	"""! @brief identify species with uncommon splice site distributions """
	
	outliers = []
	for splice_site in all_splice_sites.keys():
		values = all_splice_sites[ splice_site ]
		median = np.median( values )
		for idx, val in enumerate( values ):
			if val > 3*median and median > 0:
				outliers.append( { 'splice_site': splice_site, 'spec': species[ idx ], 'val': val } )
			elif val < 0.3*median and splice_site in [ "GT...AG", "GC...AG" ]:
				outliers.append( { 'splice_site': splice_site, 'spec': species[ idx ], 'val': val } )
	with open( outlier_file, "w" ) as out:
		for outlier in sorted( outliers, key=itemgetter('splice_site', 'val', 'spec') ):
			out.write( "\t".join( map( str, [ outlier['splice_site'], outlier['spec'], outlier['val'] ] ) ) + '\n' )


def construct_boxplot( all_splice_sites, order, fig_file ):
	"""! @brief construct box plot for all splice site combinations """
	
	values = []
	labels = []
	for idx, key in enumerate( order ):
		if key not in [ "GT...AG", "GC...AG" ]:		#"GT...AG", "GC...AG"
			values.append( all_splice_sites[ key ] )
			labels.append( key )
	fig, ax = plt.subplots( figsize=(30,5) )
	ax.boxplot( values  )
	ax.set_ylabel( "counts" )
	ax.set_xticklabels( labels, rotation=90 )
	ax.spines['right'].set_color('none')
	ax.spines['top'].set_color('none')
	plt.subplots_adjust( left=0.02, right=0.99, top=1.0, bottom=0.2 )
	fig.savefig( fig_file, dpi=600 )
	plt.close("all")


def calc_cor( all_splice_sites, a, b ):
	"""! @brief calculate correlation of ncss pattern between species A and species B """
	
	spec_a_values = []
	spec_b_values = []
	for key in all_splice_sites.keys():
		spec_a_values.append( all_splice_sites[ key ][ a ] )
		spec_b_values.append( all_splice_sites[ key ][ b ] )
	r, p = stats.spearmanr( spec_a_values, spec_b_values )
	return r,p


def correlate_patterns_between_species( all_splice_sites, species, spec_cor_file ):
	"""! @brief compare ncss distribution between species """
	
	# --- load and calculate data --- #
	sorted_spec = sorted( species )
	cor_results = []
	p_value_results = []
	for spec1 in species:
		tmp_cor = []
		tmp_p = []
		a = sorted_spec.index( spec1 )
		for spec2 in species:
			b = sorted_spec.index( spec2 )
			cor, p_value = calc_cor( all_splice_sites, a, b )
			tmp_cor.append( cor )
			tmp_p.append( p_value )
		cor_results.append( tmp_cor )
		p_value_results.append( tmp_p )
	
	# --- write data into output file --- #
	species_labels = []
	with open( spec_cor_file, "w" ) as out:
		out.write( "\t".join( [ "x" ] + species ) + '\n' )
		for idx, spec in enumerate( species ):
			out.write( "\t".join( [ spec ] + map( str, cor_results[ idx ] ) ) + '\n' )
			species_labels.append( spec.replace( "_", " " ) )
	
	# --- construct correlation plot --- #
	fig_file = spec_cor_file.replace(".txt", ".png")
	fig, ax = plt.subplots( )
	
	x_values = []
	y_values = []
	correlation = []
	
	for y, values in enumerate( cor_results ):
		for x, value in enumerate( values ):
			if value == "-":
				pass
			else:
				x_values.append( x )
				y_values.append( y )
				correlation.append( value )
	
	ax.scatter( x_values, y_values, c=correlation, s=5, marker="s", cmap="bwr", zorder=1 )
	
	for tick in ax.xaxis.get_major_ticks():
		tick.label.set_fontsize(3) 
	for tick in ax.yaxis.get_major_ticks():
		tick.label.set_fontsize(3) 
	
	ax.set_xticklabels( [ "" ] + species_labels, rotation=90, fontsize=2 )
	ax.set_yticklabels( [ "" ] + species_labels, fontsize=2 )
	
	ax.set_xlim( -1, len( cor_results ) )
	ax.set_ylim( -1, len( cor_results ) )
	
	start, end = ax.get_xlim()
	ax.xaxis.set_ticks( np.arange( start, end, 1 ) )
	
	start, end = ax.get_ylim()
	ax.yaxis.set_ticks( np.arange( start, end, 1 ) )
	
	ax.spines["top"].set_visible(False)
	ax.spines["right"].set_visible(False)
	ax.set_frame_on(False)
	
	plt.subplots_adjust( left=0.1, right=0.995, top=0.995, bottom=0.15 )
	fig.savefig( fig_file, dpi=900 )
	plt.close("all")


def splice_sites_per_species( all_splice_sites, species, total_per_spec_file, genes_per_spec ):
	"""! @brief visualize total number of splice sites per species """
	
	sorted_specs = sorted( species )
	total_ncss = []
	total_gene_numbers = []
	total_css = []
	ratios = []
	css_per_spec = {}
	ncss_per_spec = {}
	for spec in species:
		total = []
		idx = sorted_specs.index( spec )
		for key in all_splice_sites.keys():
			if key not in [ "GT...AG", "GC...AG", "AT...AC" ]:
				total.append( all_splice_sites[ key ][ idx ] )
		total_ncss.append( sum( total ) )
		total_gene_numbers.append( genes_per_spec[ spec ] )
		css = sum( total ) + all_splice_sites[ "GT...AG" ][ idx ] + all_splice_sites[ "GC...AG" ][ idx ] + all_splice_sites[ "AT...AC" ][ idx ]
		total_css.append( css )
		try:
			ratios.append( sum( total ) / float( css ) )
		except ZeroDivisionError:
			ratios.append( 0 )
		print spec + ": " + str( sum( total ) ) + " (ncss) - " + str( css ) + "(all)"
		css_per_spec.update( { spec: css } )
		ncss_per_spec.update( { spec: sum( total ) } )
	
	fig, ax = plt.subplots( figsize=(20,5) )
	ax2 = ax.twinx()
	ax3 = ax.twinx()
	ax4 = ax.twinx()
	ax.plot( np.arange( 0, len( sorted_specs ), 1 ), total_ncss, color="red", label="total_ncss", linestyle=":" )
	#ax2.plot( np.arange( 0, len( sorted_specs ), 1 ), total_css, color="blue", label="total_css", linestyle=":" )
	ax3.plot( np.arange( 0, len( sorted_specs ), 1 ), ratios, color="green", label="splice_site_ratio", linestyle=":" )
	ax4.plot( np.arange( 0, len( sorted_specs ), 1 ), total_gene_numbers, color="black", label="total_genes", linestyle=":" )
	
	print "ncss per species ranges from " + str( min( total_ncss ) ) + " to " + str( max( total_ncss ) )
	print "total number of splice sites ranges from " + str( min( total_css ) ) + " to " + str( max( total_css ) )
	
	for tick in ax.xaxis.get_major_ticks():
		tick.label.set_fontsize(10) 
	for tick in ax.yaxis.get_major_ticks():
		tick.label.set_fontsize(5) 
	
	ax.set_xlim( 0, len( species ) )
	start, end = ax.get_xlim()
	ax.xaxis.set_ticks( np.arange( start, end, 1 ) )
	ax.set_xticklabels( species, rotation=90 )
	
	ax.set_ylabel( "number of non-canonical splice sites" )
	
	ax.spines["top"].set_visible(False)
	ax.spines["right"].set_visible(False)
	ax.set_frame_on(False)
	
	my_handle = [ 	mpatches.Patch(color='red', label='total_ncss'),
								#mpatches.Patch(color='blue', label='total_css'),
								mpatches.Patch(color='green', label='ncss / css'),
								mpatches.Patch(color='black', label='total_genes')
							]
	
	ax.legend( handles=my_handle, bbox_to_anchor=(0.5, 0.9), fontsize=10 )
	
	plt.subplots_adjust( left=0.03, right=0.98, top=0.99, bottom=0.38 )
	
	fig.savefig( total_per_spec_file, dpi=300 )
	plt.close('all')
	return css_per_spec, ncss_per_spec


def get_genome_sizes( input_dir ):
	"""! @brief calculates the genome size per species based on provided assembly file """
	
	genome_sizes = {}
	assembly_files = glob.glob( input_dir + "*.fna" )
	for filename in assembly_files:
		ID = filename.split('/')[-1].split('.')[0]
		seqs = "".join( load_sequences( filename ).values() )
		genome_sizes.update( { ID: len( seqs ) } )
	return genome_sizes


def genome_size_splice_site_cor( genome_sizes, css, ncss, fig_file ):
	"""! @brief plots the correlation of css and ncss with the genome size """
	
	fig, ax = plt.subplots()
	x_values = []	#genome size
	y_values_css = []
	y_values_ncss = []
	labels = []
	for spec in genome_sizes.keys():
		x_values.append( genome_sizes[ spec ] / 1000000.0 )
		y_values_css.append( css[ spec ] )
		y_values_ncss.append( ncss[ spec ] )
		labels.append( spec )
	
	ax2 = ax.twinx()
	ax.scatter( x_values, y_values_css, color="green", label="css" )
	ax2.scatter( x_values, y_values_ncss, color="red", label="ncss" )
	
	ax.set_xlabel( "genome size [Mbp]" )
	ax.set_ylabel( "css" )
	ax2.set_ylabel( "ncss" )
	
	c1, p1 = stats.spearmanr( x_values, y_values_css )
	c2, p2 = stats.spearmanr( x_values, y_values_ncss )
	
	print "Spearman correlation css: " + str( c1 ) + " (p-value=" + str( p1 ) + ")"
	print "Spearman correlation ncss: " + str( c2 ) + " (p-value=" + str( p2 ) + ")"
	
	ax.set_title( "css: r="+str( round( c1, 4) )+", p="+str( p1 ) + "; ncss: r="+str( round( c2, 4 ) )[:4] +", p="+str( p2 ), fontsize=5 )
	ax.legend(  handles=[ mpatches.Patch(color='green', label='css'), mpatches.Patch(color='red', label='ncss') ], bbox_to_anchor=( 0.9, 0.9 ), fontsize=10 )
	
	fig.savefig( fig_file, dpi=300 )


def correlate_with_divergence( all_splice_sites, cor_fig_file ):
	"""! @brief correlate splice sites with distance to GT-AG canonical splice site """
	
	divergence = []
	counts = []
	labels = []
	for key in all_splice_sites.keys():
		if key not in [ "GT...AG" ]:	#,  "GC...AG"
			dist = 0
			for idx, nt in enumerate( "GT...AG" ):
				if nt != key[ idx ]:
					dist += 1
			divergence.append( dist )
			counts.append( np.mean( map( float, all_splice_sites[ key ] ) ) )
			labels.append( key )
	
	r, p = stats.spearmanr( divergence, counts )
	
	fig, ax = plt.subplots()
	ax.scatter( divergence, counts, color="green", s=10 )
	ax.set_yscale('log')
	ax.set_xlabel( "divergence from canonical splice site GT...AG" )
	ax.set_ylabel( "average number of observed splice sites across species" )
	
	ax.set_title( "r="+str( round( r, 4 ) )+", p-value="+str( p ) )
	
	plt.subplots_adjust( left=0.1, right=0.95, top=0.8, bottom=0.1 )
	fig.savefig( cor_fig_file, dpi=600 )
	plt.close("all")


def ncss_total_splice_site_correlation(  css_per_spec, ncss_per_spec, total_vs_ncss_fig_file ):
	"""! @brief check correlation between total number of splice sites and ncss """
	
	fig, ax = plt.subplots()
	x_values = []	#canonical splice sites
	y_values = []	#non-canonical splice sites
	for spec in css_per_spec.keys():
		x_values.append( css_per_spec[ spec ])
		y_values.append( ncss_per_spec[ spec ] )
	
	ax.scatter( x_values, y_values, color="green" )
	
	ax.set_xlabel( "number of total splice sites" )
	ax.set_ylabel( "number of non-canonical splice sites" )
	
	c1, p1 = stats.spearmanr( x_values, y_values )
	
	print "Spearman correlation between total splice sites and ncss: " + str( c1 ) + " (p-value=" + str( p1 ) + ")"
	
	ax.set_title( "total splice sites and ncss correlation: r="+str(c1)+", p="+str( p1 ), fontsize=5 )
	ax.legend( bbox_to_anchor=( 0.9, 0.9 ), fontsize=5 )
	
	fig.savefig( total_vs_ncss_fig_file, dpi=600 )


def main( arguments ):
	"""! @brief runs all analyses """
	
	input_dir = arguments[ arguments.index('--data_dir')+1 ]
	output_dir = arguments[ arguments.index('--output_dir')+1 ]
	spec_order_file = arguments[ arguments.index('--species_file')+1 ]
	
	if input_dir[-1] != '/':
		input_dir += "/"
	
	if output_dir[-1] != "/":
		output_dir += "/"
	
	if not os.path.exists( output_dir ):
		os.makedirs( output_dir )
	
	output_file = output_dir + "all_splite_sites.txt"
	combined_result_file = output_dir + "combined_div.txt"
	box_fig_file = output_dir + "boxplot.png"
	spec_cor_file = output_dir + "spec_cor.txt"
	total_per_spec_file = output_dir +  "total_per_spec.png"
	genome_size_cor_file = output_dir + "genome_size_cor.png"
	outlier_file = output_dir + "outliers2.txt"
	divergence_cor_fig_file = output_dir + "divergence_cor.png"
	total_vs_ncss_fig_file = output_dir + "total_splice_sites_vs_ncss.png"
	
	# ---- get species order --- #
	with open( spec_order_file, "r" ) as f:
		species = f.read().strip().replace(' ', '_').split('\n')
	
	# --- combine data --- #
	result_files = glob.glob( input_dir + "*.txt" )
	results = {}
	genes_per_spec = {}
	for spec in species:
		filename = input_dir + spec + ".txt"
		ncss_genes, ncss_distribution, rnas = load_results( filename )
		results.update( { spec:  ncss_distribution} )
		genes_per_spec.update( { spec: len( list( set( rnas ) ) ) } )
	all_splice_sites, order = combine_all_results( results, combined_result_file, species )
	
	# --- identify interesting outliers --- #
	find_outliers( all_splice_sites, species, outlier_file )
	
	# --- check correlation with divergence of splice sites --- #
	correlate_with_divergence( all_splice_sites, divergence_cor_fig_file )
	
	# --- check correlation with genome size --- #
	css_per_spec, ncss_per_spec = splice_sites_per_species( all_splice_sites, species, total_per_spec_file, genes_per_spec )
	genome_sizes = get_genome_sizes( input_dir )
	genome_size_splice_site_cor( genome_sizes, css_per_spec, ncss_per_spec, genome_size_cor_file )
	
	# --- correlate ncss with total number of splice sites --- #
	ncss_total_splice_site_correlation(  css_per_spec, ncss_per_spec, total_vs_ncss_fig_file )
	
	# --- compare patterns between species --- #
	correlate_patterns_between_species( all_splice_sites, species, spec_cor_file )
	
	# --- construct nice figures --- #
	construct_boxplot( all_splice_sites, order, box_fig_file )


if __name__ == '__main__':
	
	if '--data_dir' in sys.argv and '--output_dir' in sys.argv and '--species_file' in sys.argv:
		main( sys.argv )
	else:
		sys.exit( __usage__ )
