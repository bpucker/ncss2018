### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.1 ###

__usage__ = """
					python construct_NCBI_sample_overview_file.py
					--in <FULL_PATH_TO_DATA_INPUT_FOLDER>
					--out <FULL_PATH_TO_OUTPUT_FOLDER>
					"""

import glob, sys, os

# --- end of imports --- #


def load_sequences( multiple_fasta_file ):
	"""! @brief load candidate gene IDs from file """
	
	sequences = {}
	
	with open( multiple_fasta_file ) as f:
		header = f.readline()[1:].strip().split(' ')[0]
		seq = []
		line = f.readline()
		while line:
			if line[0] == '>':
					sequences.update( { header: "".join( seq ) } )
					header = line.strip()[1:].split(' ')[0]
					seq = []
			else:
				seq.append( line.strip() )
			line = f.readline()
		sequences.update( { header: "".join( seq ) } )
	return sequences


def get_intron_number( intron_file ):
	"""! @brief get intron number per species """
	
	counter = 0
	with open( intron_file, "r" ) as f:
		line = f.readline()
		while line:
			if len( line.split('\t') ) > 3:
				counter += 1
			line = f.readline()
	return counter


def main( arguments ):
	"""! @brief calculates statistics about data sets """
	
	input_dir = arguments[ arguments.index('--in')+1 ]		#NCBI data folder
	output_file = arguments[ arguments.index('--out')+1 ]		#statistics output folder
	
	genomes = glob.glob( input_dir + "*.fna" )
	
	with open( output_file, "w", 0 ) as out:
		out.write( 'Species\tGenomeSize\tNumberOfSequences\tNumberOfGenes\tNumberOfIntrons\tIntronsPerGene\tGC[%]\tCodingProportion\n' )
		for filename in sorted( genomes ):
			try:
				seqs = load_sequences( filename )
				peps = load_sequences( filename.replace( '.fna', '.pep.fa' ) )
				intron_number = get_intron_number( filename.replace( '.fna', '.txt' ) )
				genome_seq = "".join( seqs.values() ).upper()
				gc = ( genome_seq.count( 'G' ) + genome_seq.count( 'C' ) ) / ( float( len( genome_seq ) ) - genome_seq.count('N') )
				total_cds_len = 3.0*len( "".join( peps.values() ) )
				out.write( "\t".join( map( str, [ 	filename.split('/')[-1].split('.')[0],
																	len( genome_seq ),
																	len( seqs.keys() ),
																	len( peps.keys() ),
																	intron_number,
																	round( float( intron_number ) / len( peps.keys() ), 3),
																	round( gc, 5)*100,
																	round( total_cds_len / len( genome_seq ), 5 )*100
																] ) ) + '\n' )
			except:
				print filename


if __name__ == '__main__':
	if '--in' in sys.argv and '--out' in sys.argv:
		main( sys.argv )
	else:
		sys.exit( __usage__ )
