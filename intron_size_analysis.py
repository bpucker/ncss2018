### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.1 ###

__usage__ = """
					python intron_size_analysis.py
					--in <FULL_PATH_TO_INPUT_FOLDER>
					--out <FULL_PATH_TO_OUTPUT_FOLDER>
					info: runs for a long time!
					"""

import glob, os, sys
import matplotlib.pyplot as plt
from scipy import stats

# --- end of imports --- #

def calculate_bins( intron_sizes, normalization=True ):
	"""! @brief calculate size and positions of bins """
	
	bin_number = 1000.0
	upper_cutoff = 5000
	
	# --- load data from file --- #
	min_size = 20	#min( intron_sizes )
	max_size = upper_cutoff	#min( [ upper_cutoff, max( intron_sizes ) ] )
	
	step = (max_size-min_size) / bin_number
	start = min_size + 0
	end  = min_size + step
	
	x_values = []
	y_values = []
	
	while end <= max_size:
		counter = 0
		for size in intron_sizes:
			if size < end and size >= start:
				counter += 1
		x_values.append( ( start+end )/2.0 )
		y_values.append( counter )
		start += step
		end += step
	counter = 0
	for size in intron_sizes:
		if size >= end:
			counter += 1
	x_values.append( end )
	y_values.append( counter )
	
	# --- normalization --- #
	if normalization:
		y_norm = []
		total = 100.0 * sum( y_values )
		for val in y_values:
			y_norm.append( val / total )
		return x_values, y_norm
	else:
		return x_values, y_values


def run_intron_size_analysis( filename, output_file ):
	"""! @brief fun intron size analysis """
	
	# --- loading data --- #
	css = []
	ncss = []
	with open( filename, "r" ) as f:
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			if len( parts ) > 3:
				if parts[ 3 ] == "css":
					positions = sorted( map( int, parts[4].replace('_', ",").split(',') ) )
					css.append( positions[ 4 ]-positions[ 3 ] )
					if positions[ 4 ]-positions[ 3 ] == 1:
						print line
				elif parts[3] == "ncss":
					positions = sorted( map( int, parts[4].replace('_', ",").split(',') ) )
					ncss.append( positions[ 4 ]-positions[ 3 ] )
			line = f.readline()
	
	return css, ncss


def construct_combined_plot( data, fig_file, fig_file2 ):
	"""! @brief construct consensus plot """
	
	css = []
	ncss = []
	for entry in data:
		css += entry['css']
		ncss += entry['ncss']
	
	# --- plot distribution --- #
	x_css = []
	y_css = []
	x_ncss = []
	y_ncss  = []
	
	x_css, y_css = calculate_bins( css )
	x_ncss, y_ncss = calculate_bins( ncss )
	
	fig, ax = plt.subplots()
	ax.plot( x_css, y_css, color="green", label="css", alpha=0.5 )
	ax.plot( x_ncss, y_ncss, color="red", label="ncss", alpha=0.5 )
	
	ax.legend( bbox_to_anchor=( 0.9, 0.9 ) )
	
	ax.set_xlabel("intron length [bp]")
	ax.set_ylabel( "frequency [%]" )
	
	fig.savefig( fig_file, dpi=300 )
	plt.close("all")
	
	data_output_file = fig_file + "_values.txt"
	with open( data_output_file, "w" ) as out:
		out.write( "\t".join( map( str, x_css ) ) + '\n' )
		out.write( "\t".join( map( str, y_css ) ) + '\n' )
		out.write( "\t".join( map( str, y_ncss ) ) + '\n' )
	
	print stats.wilcoxon( y_css, y_ncss, zero_method="pratt" )	
	
	# --- analyze how many introns lengths are multiples of 3 --- #
	css = []
	ncss = []
	for each in data:
		css += each['css']
		ncss += each['ncss']
	
	ncss_no3 = 0
	ncss_3 = 0
	
	css_no3 = 0
	css_3 = 0
	
	for each in css:
		if each % 3 == 0:
			css_3 += 1
		else:
			css_no3 += 1
	for each in ncss:
		if each % 3 == 0:
			ncss_3 += 1
		else:
			ncss_no3 += 1
	
	print ncss_no3
	print ncss_3
	
	print css_no3
	print css_3
	
	print stats.chisquare( [ ncss_3, ncss_no3 ], [ css_3, css_no3 ] )


def main( arguments ):
	"""! @brief runs all parts """
	
	data_dir = arguments[ arguments.index('--in')+1 ]	#NCBI splice site check result folder
	output_dir = arguments[ arguments.index('--out')+1 ]	#output folder
	
	if not os.path.exists( output_dir ):
		os.makedirs( output_dir )
	
	filenames = glob.glob( data_dir + "*.txt" )
	data = []
	fig_file = output_dir + "consensus_plot.png"
	fig_file2 = output_dir + "frame_check.png"
	for filename in filenames:
		ID = filename.split('/')[-1].split('.')[0]
		output_file = output_dir + ID + ".png"
		css, ncss = run_intron_size_analysis( filename, output_file )
		data.append( { 'css': css, 'ncss': ncss } )
	construct_combined_plot( data, fig_file, fig_file2 )


if __name__ == '__main__':
	if '--in' in sys.argv and '--out' in sys.argv:
		main( sys.argv )
	else:
		sys.exit( __usage__ )
