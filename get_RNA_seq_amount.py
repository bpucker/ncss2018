### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.1 ###

__usage__ = """
					python get_RNA_seq_amount.py
					--cov_dir <FULL_PATH_TO_FOLDER_WITH_COVERAGE_FILES>
					--out <FULL_PATH_TO_OUTPUT_FILE>
					"""

import glob, os, sys

# --- end of imports --- #

def main( arguments ):
	"""! @brief get amount of RNA-Seq reads """
	
	cov_file_dir = arguments[ arguments.index( '--cov_dir' ) +1 ]
	output_file = arguments[ arguments.index( '--out' ) +1 ]
	
	cov_files = glob.glob( cov_file_dir + "*.cov" )
	print "number of detected cov files: " + str( len( cov_files ) )
	
	with open( output_file, "a", 0 ) as out:
		for filename in sorted( cov_files ):
			total_cov = []
			with open( filename, "r" ) as f:
				line = f.readline()
				while line:
					parts = line.strip().split('\t')
					total_cov.append( float( parts[-1] ) )
					line = f.readline()
					if len( total_cov ) == 1000000:
						total_cov = [ sum( total_cov ) ]
			out.write( filename.split('/')[-1].split('.')[0] + '\t' + str( sum( total_cov ) ) + '\n' )


if __name__ == '__main__':
	if '--cov_dir' in sys.argv and '--out' in sys.argv:
		main( sys.argv )
	else:
		sys.exit( __usage__ )
