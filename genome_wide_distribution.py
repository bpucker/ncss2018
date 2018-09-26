### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.1 ###

#note: some code parts are derived from other scripts e.g. in https://github.com/bpucker/script_collection

__usage__ = """
					python genome_wide_distribution.py
					--in <FULL_PATH_TO_DATA_FOLDER>
					--out <FULL_PATH_TO_OUTPUT_FOLDER>
					info: needs results of splice site diversity analysis
					(identified non-canonical splice sites)
					"""

import os, glob, re
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from operator import itemgetter

# --- end of imports --- #

def load_annotation( annotation_file ):
	"""! @brief load all gene positions from gff3 file """
	
	annotation = {}
	
	with open( annotation_file, "r" ) as f:
		line = f.readline()
		while line:
			if line[0] != '#':
				parts = line.strip().split('\t')
				if parts[2] in [ "mRNA" ]:
					ID = re.findall( "rna\d+", parts[-1] )[0]
					annotation.update( { ID: { 'id': ID, 'chr': parts[0], 'start': int( parts[3] ), 'end': int( parts[4] ), 'orientation': parts[6] } } )
			line = f.readline()
	return annotation


def load_seq_lengths( multiple_fasta_file, len_cutoff=10000000 ):
	"""! @brief load candidate gene IDs from file """
	
	seq_lens = {}
	
	with open( multiple_fasta_file ) as f:
		header = f.readline()[1:].strip().split(' ')[0]
		seq = ""
		line = f.readline()
		while line:
			if line[0] == '>':
				if len( seq ) > len_cutoff:
					seq_lens.update( { header: len( seq ) } )
				header = line.strip()[1:].split(' ')[0]
				seq = ""
			else:
				seq += line.strip()
			line = f.readline()
		if len( seq ) > len_cutoff:
			seq_lens.update( { header: len( seq ) } )
	return seq_lens


def plot_genome_wide_distribution( genes_to_plot, chr_lengths, fig_output_file ):
	"""! @brief show genome wide distribution of genes """
	
	# --- construct plot --- #	
	fig, ax = plt.subplots( figsize=( 10, int( len( chr_lengths.keys() ) / 2 ) ) )
	chr_names = sorted( chr_lengths.keys() )
	y_offset = len( chr_names )
	
	
	# --- adding chromosomes --- #
	for idx, each in enumerate( chr_names ):
		ax.plot( [ 0,  chr_lengths[ each ]/1000000.0 ], [ y_offset-idx, y_offset-idx ] , color="black", linewidth=.5 )
		ax.text( 0, y_offset-idx+0.2, each, fontsize=5 )	#chr_lengths[ each ]/1000000.0
	
	
	# --- adding gene positions --- #
	redx = []
	redy = []
	blackx = []
	blacky = []
	for gene in genes_to_plot:
		try:
			y = y_offset-chr_names.index( gene['chr'] )
			x = ( gene['start']+gene['end'] ) / 2000000.0
			if gene['color'] == "red":
				redx.append( x )
				redy.append( y )
			else:
				blackx.append( x )
				blacky.append( y )
		except ValueError:
			pass	#print gene
			
	ax.scatter( redx, redy, s=1, color="red", zorder=3 )
	ax.scatter( blackx, blacky, s=1, color="black", zorder=2 )
	
	# --- improving overall layout --- #
	ax.set_xlabel( "chromosome position [Mbp]" )
	
	ax.spines["top"].set_visible(False)
	ax.spines["left"].set_visible(False)
	ax.spines["right"].set_visible(False)
	ax.set_frame_on(False)
	ax.axes.get_yaxis().set_visible(False)
	
	ax.set_ylim( 0, len( chr_names )+1 )
	
	if int( len( chr_lengths.keys() ) / 2 ) <= 3:
		bottom_value = 0.25
	elif int( len( chr_lengths.keys() ) / 2 ) <= 5:
		bottom_value = 0.15
	elif int( len( chr_lengths.keys() ) / 2 ) <= 8:
		bottom_value = 0.1
	elif int( len( chr_lengths.keys() ) / 2 ) <= 10:
		bottom_value = 0.075
	else:
		bottom_value = 0.05
	
	ax.legend(  handles=[ mpatches.Patch(color='black', label='css genes'), mpatches.Patch(color='red', label='ncss genes') ], bbox_to_anchor=( 0.9, 0.9 ), fontsize=5, alpha=0.5 )
	ax.xaxis.set_tick_params(labelsize=5)
	plt.subplots_adjust( left=0.0, right=0.98, top=1.0, bottom=bottom_value )
	
	fig.savefig( fig_output_file, dpi=300 )
	plt.close("all")


def load_splice_site_genes( input_file ):
	"""! @brief load and classify all transcript IDs """
	
	css_IDs = []
	ncss_IDs = []
	data = {}
	with open( input_file, "r" ) as f:
		f.readline()	#header
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			if len( parts ) > 3:
				if parts[3] == 'ncss':
					try:
						data[ parts[0] ] = "ncss"
					except KeyError:
						data.update( { parts[0]: "ncss" } )
				else:
					try:
						data[ parts[0] ]
					except KeyError:
						data.update( { parts[0]: "css" } )
			line = f.readline()
	for key in data.keys():
		if data[ key ] == "ncss":
			ncss_IDs.append( key )
		else:
			css_IDs.append( key )
	return css_IDs, ncss_IDs


def main( arguments ):
	"""! @brief construct figures to illustrate genome-wide distribution of non-canonical splice sites """
	
	input_dir = arguments[ arguments.index( '--in' )+1 ]	#NCBI genome ceck folder
	output_dir = arguments[ arguments.index( '--out' )+1 ]	#output folder
	
	input_files = glob.glob( input_dir + "*.txt" )
	for filename in input_files:
		ID = filename.split('/')[-1].split('.')[0]
		try:
			fig_output_file = output_dir + ID + "_genome_wide_distribution.png"
			
			annotation_file = filename.replace( ".txt", ".gff" )
			genome_seq_file = filename.replace( ".txt", ".fna" )
			
			css_IDs, ncss_IDs = load_splice_site_genes( filename )
			
			annotation = load_annotation( annotation_file )
			
			# --- construct figure of genome wide distribution --- #
			genes_to_plot = []
			for gene in css_IDs:
				try:
					gene = annotation[ gene ]
					gene.update( { 'color': "black" } )
					genes_to_plot.append( gene )
				except KeyError:
					print gene
			for gene in ncss_IDs:
				gene = annotation[ gene ]
				gene.update( { 'color': "red" } )
				genes_to_plot.append( gene )
			
			chr_lengths = load_seq_lengths( genome_seq_file )
			plot_genome_wide_distribution( genes_to_plot, chr_lengths, fig_output_file )
		except:
			print ID
	print "all done!"


if __name__ == '__main__':
	if '--in' in sys.argv and '--out' in sys.argv:
		main( sys.argv )
	else:
		sys.exit( __usage__ )
	
