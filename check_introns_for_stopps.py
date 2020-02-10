### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.2 ###

__usage__ = """
					python check_introns_for_stopps.py
					--in <FULL_PATH_TO_INPUT_FOLDER>
					--out <FULL_PATH_TO_OUTPUT_FOLDER>
					
					bug reports and feature requests: bpucker@cebitec.uni-bielefeld.de
					"""

import glob, os, sys, re
from operator import itemgetter

# --- end of imports --- #

def load_sequences( fasta_file ):
	"""! @brief load candidate gene IDs from file """
	
	sequences = {}
	
	with open( fasta_file ) as f:
		header = f.readline()[1:].strip().split(" ")[0]
		seq = ""
		line = f.readline()
		while line:
			if line[0] == '>':
					sequences.update( { header: seq } )
					header = line.strip()[1:].split(" ")[0]
					seq = ""
			else:
				seq += line.strip()
			line = f.readline()
		sequences.update( { header: seq } )
	return sequences


def run_intron_size_analysis( filename ):
	"""! @brief fun intron size analysis """
	
	# --- loading data --- #
	css = []
	gcag = []
	atac = []
	other = []
	
	total_css_counter = 0
	total_gcag_counter = 0
	total_atac_counter = 0
	total_other_counter = 0
	
	with open( filename, "r" ) as f:
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			if len( parts ) > 3:
				status = True
				for each in parts[1] + parts[2]:
					if each not in "ACGT":
						status = False
				if status:
					control_length = int( parts[-1] )
					intron_blocks = parts[4].split(',')
					left_values = map( int, intron_blocks[1].split('_') )
					right_values = map( int, intron_blocks[2].split('_') )
					if left_values[0] < right_values[0]:	#gene on forward strand
						start = left_values[0] - 1
						end = right_values[1] - 1
						orientation = True
					else:	#gene on reverse strand
						start = right_values[0] - 1
						end = left_values[1] -1
						orientation = False
					
					if end-start != control_length:
						print "ERROR: " + str( control_length ) + "\t" + str( start ) + "\t" + str( end )
					
					if parts[1] + parts[2] == "GTAG":
						total_css_counter += 1
						if control_length % 3 == 0:
							css.append( { 'ID': parts[0], 'start': start, 'end': end, 'orientation': orientation } )
					elif parts[1] + parts[2] == "GCAG":
						total_gcag_counter += 1
						if control_length % 3 == 0:
							gcag.append( { 'ID': parts[0], 'start': start, 'end': end, 'orientation': orientation } )
					elif parts[1] + parts[2] == "ATAC":
						total_atac_counter += 1
						if control_length % 3 == 0:
							atac.append( { 'ID': parts[0], 'start': start, 'end': end, 'orientation': orientation } )
					else:
						total_other_counter += 1
						if control_length % 3 == 0:
							other.append( { 'ID': parts[0], 'start': start, 'end': end, 'orientation': orientation } )
					
			line = f.readline()
	
	return css, gcag, atac, other, total_css_counter, total_gcag_counter, total_atac_counter, total_other_counter 


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
					if parts[6] == "+":
						orientation = True
					else:
						orientation = False
					try:
						pos.update( { ID: { 'id': ID, 'parent': parent, 'chr': parts[0], 'orientation': orientation } } )
					except:
						pass
			line = f.readline()
	return pos


def analyze_introns( data, rna_info, seqs ):
	"""! @brief screen introns of length 3x for in frame stopp codons """
	
	clean_counter = 0
	stopp_counter = 0
	
	for each in data:
		chromosome = rna_info[ each['ID'] ]['chr']
		seq = seqs[ chromosome ][ each['start']: each['end'] ].upper()
		if seq.count( "A" ) + seq.count( "C" ) +  seq.count( "G" ) + seq.count( "T" ) == len( seq ):
			# ---- count stopp codons --- #
			codons = [ seq[ i : i+3 ] for i in range( 0, len( seq ), 3 ) ]
			if each['orientation']:
				stopps = codons.count( "TAA" ) + codons.count( "TGA" ) + codons.count( "TAG" )
			else:
				stopps = codons.count( "TTA" ) + codons.count( "TCA" ) + codons.count( "CTA" )
			
			# --- count introns with stopp / without stopp codons --- #
			if stopps == 0:
				clean_counter += 1
			else:
				stopp_counter += 1
		else:
			pass	#print "ERROR: ambiguity character detected - skipping sequence"
	return clean_counter, stopp_counter


def main( arguments ):
	"""! @brief runs all parts """
	
	data_dir = arguments[ arguments.index('--in')+1 ]	#NCBI splice site check result folder
	output_dir = arguments[ arguments.index('--out')+1 ]	#output folder
	
	if not os.path.exists( output_dir ):
		os.makedirs( output_dir )
	
	cutoff = 5000
	small_intron_ncss = {}
	large_intron_ncss = {}
	
	filenames = sorted( glob.glob( data_dir + "*.txt" ) )
	data = []
	report_file = output_dir + "report.txt"
	
	
	global_css_counter = 0
	global_gcag_counter = 0
	global_atac_counter = 0
	global_other_counter = 0
	
	
	d3_global_css_counter = 0
	d3_global_gcag_counter = 0
	d3_global_atac_counter = 0
	d3_global_other_counter = 0
	
	d3_global_css_clean_counter = 0
	d3_global_gcag_clean_counter = 0
	d3_global_atac_clean_counter = 0
	d3_global_other_clean_counter = 0
	
	
	with open( report_file, "w" ) as out:
		out.write( "\t".join( [ "SpeciesID",
										"Number_GTAG_3D", "Number_GTAG", "GTAG_CLEAN", "GTAG_STOPP",
										"Number_GCAG_3D", "Number_GCAG", "GCAG_CLEAN", "GCAG_STOPP",
										"Number_ATAC_3D", "Number_ATAC", "ATAC_CLEAN", "ATAC_STOPP",
										"Number_OTHER_3D", "Number_OTHER", "OTHER_CLEAN", "OTHER_STOPP"
										] ) + '\n'  )
		for filename in filenames:
			ID = filename.split('/')[-1].split('.')[0]
			print ID
			gff_file = data_dir + ID + ".gff"
			fasta_file = data_dir + ID + ".fna"
			css, gcag, atac, other, total_css_counter, total_gcag_counter, total_atac_counter, total_other_counter = run_intron_size_analysis( filename )
			
			global_css_counter += total_css_counter
			global_gcag_counter += total_gcag_counter
			global_atac_counter += total_atac_counter
			global_other_counter += total_other_counter			
			
			d3_global_css_counter += len( css )
			d3_global_gcag_counter += len( gcag )
			d3_global_atac_counter += len( atac )
			d3_global_other_counter += len( other )
			
			rna_info = load_pos_from_gff( gff_file )
			seqs = load_sequences( fasta_file )
			css_info = analyze_introns( css, rna_info, seqs )
			gcag_info = analyze_introns( gcag, rna_info, seqs )
			atac_info = analyze_introns( atac, rna_info, seqs )
			other_info = analyze_introns( other, rna_info, seqs )
			
			d3_global_css_clean_counter += css_info[0]
			d3_global_gcag_clean_counter += gcag_info[0]
			d3_global_atac_clean_counter += atac_info[0]
			d3_global_other_clean_counter += other_info[0]
			
			out.write( "\t".join( map( str, [ 	ID,
																len( css ), total_css_counter,	css_info[0], css_info[1], 
																len( gcag ),  total_gcag_counter, gcag_info[0], gcag_info[1],
																len( atac ),  total_atac_counter, atac_info[0], atac_info[1],
																len( other ),  total_other_counter, other_info[0], other_info[1]																
															] ) ) + "\n" )
		
	
	print "GT-AG introns divisible by 3: " + str(  d3_global_css_counter / float( global_css_counter ) )
	print "GC-AG introns divisible by 3: " + str(  d3_global_gcag_counter / float( global_gcag_counter ) )
	print "AT-AC introns divisible by 3: " + str(  d3_global_atac_counter / float( global_atac_counter ) )
	print "OTHER introns divisible by 3: " + str(  d3_global_other_counter / float( global_other_counter ) )

	
	print "GT-AG introns divisible by 3: " + str(  d3_global_css_clean_counter / float( global_css_counter ) ) + "\tn=" + str( d3_global_css_clean_counter )
	print "GC-AG introns divisible by 3: " + str(  d3_global_gcag_clean_counter / float( global_gcag_counter ) ) + "\tn=" + str( d3_global_gcag_clean_counter )
	print "AT-AC introns divisible by 3: " + str(  d3_global_atac_clean_counter / float( global_atac_counter ) ) + "\tn=" + str( d3_global_atac_clean_counter )
	print "OTHER introns divisible by 3: " + str(  d3_global_other_clean_counter / float( global_other_counter ) ) + "\tn=" + str( d3_global_other_clean_counter )


if __name__ == '__main__':
	if '--in' in sys.argv and '--out' in sys.argv:
		main( sys.argv )
	else:
		sys.exit( __usage__ )
