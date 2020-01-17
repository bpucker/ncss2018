### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.1 ###

__usage__ = """
					python analyze_supported_splice_sites.py
					--data <FULL_PATH_TO_DATA_FOLDER>
					--sssd <FULL_PATH_TO_RNA_SEQ_SUPPORT_FILES> #path to all species (not one particular)
					--cov_rep <FULL_PATH_TO_RNA_SEQ_COVERAGE_REPORT_FILE> #names have to be the same as in NCBI folder; RNA-Seq amount
					--out <FULL_PATH_TO_OUTPUT_FOLDER>
					"""

import glob, re, sys, os
import matplotlib.pyplot as plt
from scipy import stats
import numpy as np
import warnings

# --- end of imports --- #

def load_splice_site_counts( filename ):
	"""! @brief load counts for all splice sites """
	
	gtag, gcag, atac, others = 0, 0, 0, 0
	
	with open( filename, "r" ) as f:
		line = f.readline()
		while line:
			if 'GT..AG' in line:
				gtag = int( re.findall( "\d+", line )[0] )
			elif 'GC..AG' in line:
				gcag = int( re.findall( "\d+", line )[0] )
			elif 'AT..AC' in line:
				atac = int( re.findall( "\d+", line )[0] )
			elif not 'N' in line:
				others += int( re.findall( "\d+", line )[0] )
			line = f.readline()
	return gtag, gcag, atac, others


def generate_brief_overview( input_dir, output_file ):
	"""! @brief generate overview based on top (non-canonical) splice sites """
	
	overview_files = glob.glob( input_dir + "*/*_overview.txt" )
	print(overview_files)
	specs = []
	gtags = []
	gcags = []
	atacs = []
	all_others = []
	
	for filename in sorted( overview_files ):
		spec = filename.split('/')[-2]
		gtag, gcag, atac, others = load_splice_site_counts( filename )
		specs.append( spec )
		gtags.append( gtag )
		gcags.append( gcag )
		atacs.append( atac )
		all_others.append( others )
	
	with open( output_file, "w" ) as out:
		out.write( "splice_site\t" + "\t".join( specs ) + '\n' )
		out.write( "GT-AG\t" + "\t".join( map( str, gtags ) ) + '\n' )
		out.write( "GC-AG\t" + "\t".join( map( str, gcags ) ) + '\n' )
		out.write( "AT-AC\t" + "\t".join( map( str, atacs ) ) + '\n' )
		out.write( "others\t" + "\t".join( map( str, all_others ) ) + '\n' )



def correlate_support_and_coverage( data_dir, splice_site_support_dir, cov_report_file, output_dir ):
	"""! @brief correlate support of splice sites with available RNA-Seq read coverage """
	
	# --- load available RNA-Seq read coverage from file --- #
	cov_per_spec = {}
	with open( cov_report_file, "r" ) as f:
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			cov_per_spec.update( { parts[0].lower(): float( parts[-1] ) } )
			line = f.readline()
	
	# --- get total number of annotated splice sites per species --- #
	total_splice_sites_per_spec = {}
	total_ncss_per_spec = {}
	total_splice_site_files = glob.glob( data_dir + "*.txt" )
	for filename in total_splice_site_files:
		counter = 0
		ncss_counter = 0
		ID = filename.split('/')[-1].split('.')[0].lower()
		with open( filename, "r" ) as f:
			line = f.readline()
			while line:
				parts = line.strip().split('\t')
				if len( parts ) > 3:
					counter += 1
					if parts[3] == "ncss":
						ncss_counter += 1
				line = f.readline()
		total_splice_sites_per_spec.update( { ID: counter } )
		total_ncss_per_spec.update( { ID: ncss_counter } )
		
	# --- get total number of supported splice sites per species --- #
	supported_splice_sites_per_spec = {}
	supported_ncss_per_spec = {}
	supported_splice_site_files = glob.glob( splice_site_support_dir + "*/supported_ncss.txt" )
	#name should be changed to "supported_splice_sites.txt"
	for filename in supported_splice_site_files:
		counter = 0
		ncss_counter = 0
		ID = filename.split('/')[-2].lower()
		with open( filename, "r" ) as f:
			line = f.readline()
			while line:
				parts = line.strip().split('\t')
				if len( parts ) == 3:
					if ( ( parts[1] == "GT" ) + ( parts[2] == "AG" ) ) < 2:
						ncss_counter += 1
				counter += 1
				line = f.readline()
		supported_splice_sites_per_spec.update( { ID: counter } )
		supported_ncss_per_spec.update( { ID: ncss_counter } )
	
	# --- constructing general correlation figure --- #
	fig_file = output_dir + "RNA_seq_cov_splice_site_support_correlation.png"
	fig, ax = plt.subplots()
	
	x_values = []	#canonical splice sites
	y_values = []	#non-canonical splice sites
	y_values2 = []	#percentage
	for spec in cov_per_spec.keys():
		try:
			y_values.append( supported_splice_sites_per_spec[ spec ] )
			percent = (100.0*supported_splice_sites_per_spec[ spec ] ) / total_splice_sites_per_spec[ spec ] 
			print spec + "\t" + str( percent ) + "%"
			y_values2.append( percent )
			x_values.append( cov_per_spec[ spec ]/ 1000000000.0 )
		except KeyError:
			print spec
	
	ax.scatter( x_values, y_values, color="green", label="counts" )
	ax.scatter( [], [], color="blue", label="percent" )
	ax2 = ax.twinx()
	ax2.scatter( x_values, y_values2, color="blue", label="percent" )
	
	ax.set_xlabel( "total number of sequenced nt in RNA-Seq data [billion]" )
	ax.set_ylabel( "number of supported splice sites" )
	ax2.set_ylabel( "percent of supported splice sites" )
	
	c1, p1 = stats.spearmanr( x_values, y_values )
	
	print "Spearman correlation between number supported splice sites and nucleotides sequenced: " + str( c1 ) + " (p-value=" + str( p1 ) + ")"
	
	ax.set_title( "RNA-Seq coverage and splice site support: r="+str(c1)+", p="+str( p1 ), fontsize=5 )
	ax.legend( bbox_to_anchor=( 0.9, 0.9 ), fontsize=5 )
	
	plt.subplots_adjust( left=0.15, top=0.95, right=0.9, bottom=0.12 )
	
	fig.savefig( fig_file, dpi=300 )
	
	
	# --- construct ncss correlation figure --- #
	fig_file = output_dir + "RNA_seq_cov_ncss_support_correlation.png"
	fig, ax = plt.subplots()
	
	x_values = []	#canonical splice sites
	y_values = []	#non-canonical splice sites
	y_values2 = []	#percentage
	for spec in cov_per_spec.keys():
		try:
			try:
				percent = (100.0*supported_ncss_per_spec[ spec ] ) / total_ncss_per_spec[ spec ] 
				y_values.append( supported_ncss_per_spec[ spec ] )
				print spec + "\t" + str( percent ) + "%"
				y_values2.append( percent )
				x_values.append( cov_per_spec[ spec ] / 1000000000.0 )
			except ZeroDivisionError:
				print spec + " ZeroDivisionError"	
		except KeyError:
			print spec + " KeyError"
	
	ax.scatter( x_values, y_values, color="green", label="counts" )
	ax.scatter( [], [], color="blue", label="percent" )
	ax2 = ax.twinx()
	ax2.scatter( x_values, y_values2, color="blue", label="percent" )
	
	ax.set_xlabel( "total number of sequenced nt in RNA-Seq data [billion]" )
	ax.set_ylabel( "number of supported non-canonical splice sites" )
	ax2.set_ylabel( "percent of supported non-canonical splice sites" )
	
	c1, p1 = stats.spearmanr( x_values, y_values )
	
	print "Spearman correlation between number supported ncss and nucleotides sequenced: " + str( c1 ) + " (p-value=" + str( p1 ) + ")"
	
	ax.set_title( "RNA-Seq coverage and splice site support: r="+str(c1)+", p="+str( p1 ), fontsize=5 )
	ax.legend( bbox_to_anchor=( 0.9, 0.9 ), fontsize=5 )
	
	plt.subplots_adjust( left=0.15, top=0.95, right=0.9, bottom=0.12 )
	
	fig.savefig( fig_file, dpi=300 )
	

def construct_overview_figure( overview_file, overview_figure ):
	"""! @brief construct a figure to illustrate the ratio between differen splice site combinations """
	
	with open( overview_file, "r" ) as f:
		specs = f.readline().strip().split('\t')[1:]
		gt_ag = map( int, f.readline().strip().split('\t')[1:] )
		gc_ag = map( int, f.readline().strip().split('\t')[1:] )
		at_ac = map( int, f.readline().strip().split('\t')[1:] )
		all_others = map( int, f.readline().strip().split('\t')[1:] )
	
	gtag = []
	gcag = []
	atac = []
	others = []
	for idx, each in enumerate( specs ):
		gtag.append( float( gt_ag[ idx ] ) / ( gt_ag[ idx ] + gc_ag[ idx ] + at_ac[ idx ] +all_others[ idx ] ) )
		gcag.append( float( gc_ag[ idx ] ) / ( gt_ag[ idx ] + gc_ag[ idx ] + at_ac[ idx ] +all_others[ idx ] ) )
		atac.append( float( at_ac[ idx ] ) / ( gt_ag[ idx ] + gc_ag[ idx ] + at_ac[ idx ] +all_others[ idx ] ) )
		others.append( float( all_others[ idx ] ) / ( gt_ag[ idx ] + gc_ag[ idx ] + at_ac[ idx ] +all_others[ idx ] ) )
	
	fig, ax = plt.subplots()
	
	ax.boxplot( [ gtag, gcag, atac, others ] )
	ax.set_yscale('log')
	
	ax.set_xticklabels( [ "GT-AG", "GC-AG", "AT-AC", "others" ] )
	ax.set_ylabel( "proportion of total splice sites" )
	
	fig.savefig( overview_figure, dpi=300 )
	plt.close('all')
	
	
	print "GT-AG (median): " + str( 100*np.median( gtag ) )[:5] + "%"
	print "GC-AG (median): " + str( 100*np.median( gcag ) )[:5] + "%"
	print "AT-AC (median): " + str( 100*np.median( atac ) )[:5] + "%"
	print "others (median): " + str( 100*np.median( others ) )[:5] + "%"
	


def construct_combined_file_with_splice_site_support( splice_site_support_dir, output_dir ):
	"""! @brief combine all files with splice site support values """
	
	output_file = output_dir + "all_supported_splice_sites.txt"
	
	input_files = sorted( glob.glob( splice_site_support_dir + "*/splice_site_coverage_check.txt" ) )
	with open( output_file, "w" ) as out:
		out.write( "Species\tGeneID\tExon3prime\tIntron5prime\tIntron3prime\tExon5prime\t5prime_splice_site\t3prime_splice_site\n" )
		for filename in input_files:
			ID = filename.split('/')[-2]
			with open( filename, "r" ) as f:
				f.readline()	#header
				line = f.readline()
				while line:
					out.write( ID + '\t' + line )
					line = f.readline()


def main( arguments ):
	"""! @brief assess correlation between RNA-Seq coverage and supported splice sites """
	
	data_dir = arguments[ arguments.index( '--data' )+1 ]	#data folder with all NCBI files after processing
	splice_site_support_dir = arguments[ arguments.index( '--sssd' )+1 ]	#RNA-Seq ncss support folder
	cov_report_file = arguments[ arguments.index( '--cov_rep' )+1 ]	#RNA-seq coverage overview file
	output_dir = arguments[ arguments.index( '--out' )+1 ]	#output folder
	
	if output_dir[-1] != '/':
		output_dir += "/"
	if not os.path.exists( output_dir ):
		os.makedirs( output_dir )
	
	overview_file = output_dir + "overview.txt"
	generate_brief_overview( splice_site_support_dir, overview_file )
	overview_figure = output_dir + "overview.png"
	construct_overview_figure( overview_file, overview_figure )
	
	
	#analyze percentage of supported splice sites and correlate it with coverage
	correlate_support_and_coverage( data_dir, splice_site_support_dir, cov_report_file, output_dir )
	
	#construct combined splice site coverage file
	construct_combined_file_with_splice_site_support( splice_site_support_dir, output_dir )
	print "all done!"


if __name__ == "__main__":
	if '--data' in sys.argv and '--sssd' in sys.argv and '--cov_rep' in sys.argv and '--out' in sys.argv:
		main( sys.argv )
	else:
		sys.exit( __usage__ )
