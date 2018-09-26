### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.5 ###

__usage__ = """
	python NCBI_splice_site_check.py
	--in <FULL_PATH_TO_INPUT_DIR>
	info: all .fna files will be processed; gff3 files with same name are expected
	some functions are loosly based on code used in Pucker et al., 2017; doi:10.1186/s13104-017-2985-y
	bug reports and feature requests: bpucker@cebitec.uni-bielefeld.de
					"""


import re, sys, glob
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


def revcomp( seq ):
	"""! @brief construct reverse complement of sequence """
	
	new_seq = []
	
	bases = { 'a':'t', 't':'a', 'c':'g', 'g':'c', 'n':'n', '-':'-' }
	for nt in seq.lower():
		try:
			new_seq.append( bases[nt] )
		except:
			new_seq.append( 'n' )
	return ''.join( new_seq[::-1] ).upper()


def remove_non_repr_transcripts( transcripts, rna_to_gene ):
	"""! @brief remove non representative transcripts from collection """
	#each transcript has list of exon/CDS features
	#start, end of all CDS features are used for pep length calculation
	#take transcript with longest encoded peptide as representative
	
	# --- assign transcripts to genes --- #
	transcripts_per_genes = {}
	for key in transcripts.keys():
		try:
			transcripts_per_genes[ rna_to_gene[ key ] ].append( key )
		except KeyError:
			try:
				transcripts_per_genes.update( { rna_to_gene[ key ]: [ key ] } )
			except KeyError:
				transcripts_per_genes.update( { key: [ key ] } )
	
	# --- get representative transcripts --- #
	repr_transcripts = {}
	for gene in transcripts_per_genes.keys():
		if len( transcripts_per_genes[ gene ] ) == 1:
			repr_transcripts.update( { transcripts_per_genes[gene][0]: transcripts[ transcripts_per_genes[gene][0] ] } )
		else:
			CDS_len_per_transcript = []
			for ID in transcripts_per_genes[ gene ]:
				counter = 0
				for element in transcripts[ ID ]:
					if element['type'] == "CDS":
						counter += element['end']-element['start']
				CDS_len_per_transcript.append( { 'id': ID, 'len': counter } )
			repr_trans = sorted( CDS_len_per_transcript, key=itemgetter('len') )[-1]
			repr_transcripts.update( { repr_trans['id']: transcripts[ repr_trans['id'] ] } )
	return repr_transcripts


def construct_codingseqs( transcripts, transcript_file, genome_seq ):
	"""! @brief construct all transcripts based on given CDS/exon positions """
	
	with open( transcript_file, "w" ) as t_out:
		for ID in transcripts.keys():
			features = transcripts[ ID ]
			CDS_parts = []
			for feature in features:
				if feature['type'] == "CDS":
					CDS_parts.append( feature )
			if len( CDS_parts ) > 0:
				seq = []
				if CDS_parts[0]['orientation'] == "+":
					for CDS in CDS_parts:
						seq.append( genome_seq[ CDS['chr'] ][ CDS['start']-1:CDS['end'] ] )
					t_out.write( '>' + ID + '\n' + "".join( seq ).upper() + '\n' )
				else:
					for CDS in CDS_parts[::-1]:
						seq.append( genome_seq[ CDS['chr'] ][ CDS['start']-1:CDS['end'] ] )
					t_out.write( '>' + ID + '\n' + revcomp( "".join( seq )).upper() + '\n' )
			else:
				print "ERROR: no CDS features detected - " + ID


def get_all_intron_borders_per_gene( gff_file, ref_seq, out, transcript_file ):
	"""! @brief check all transcript for non-canonical splice sites """
	
	# --- identify all protein coding genes --- #
	relevant_genes = {}
	rna_to_gene = {}
	with open( gff_file, "r" ) as f:
		line = f.readline()
		while line:
			if line[0] != '#':
				parts = line.strip().split('\t')
				if parts[2] == "mRNA":
					ID = re.findall( "rna\d+", parts[-1] )[0]
					parent = re.findall( "gene\d+", parts[-1] )[0]
					relevant_genes.update( { parent: None } )
					rna_to_gene.update( { ID: parent } )
			line = f.readline()
	
	# --- loading all data --- #
	transcripts = {}
	with open( gff_file, "r" ) as f:
		line = f.readline()
		while line:
			if line[0] != '#':
				parts = line.strip().split('\t')
				if parts[2] in [ "CDS" ]:	#"exon",
					try:
						status = False
						try:
							mRNA_id = re.findall( "rna\d+", line )[0]
							relevant_genes[ rna_to_gene[ mRNA_id ] ]
							status = True
						except IndexError:
							try:
								mRNA_id = re.findall( "gene\d+", line )[0]
								relevant_genes[ mRNA_id ]
								status = True
							except IndexError:
								print parts[-1][:50]
						if status:
							try:
								transcripts[ mRNA_id ].append( { 'type': parts[2], 'start': int( parts[3] ), 'end': int( parts[4] ), 'orientation': parts[6], 'chr': parts[0] } )
							except KeyError:
								transcripts.update( { mRNA_id: [ { 'type': parts[2], 'start': int( parts[3] ), 'end': int( parts[4] ), 'orientation': parts[6], 'chr': parts[0] } ] } )
					except KeyError:	#not protein coding gene
						pass
			line = f.readline()
	print "number of identified transcripts: " + str( len( transcripts.keys() ) )
	
	# --- remove non-representative transcripts --- #
	transcripts = remove_non_repr_transcripts( transcripts, rna_to_gene )
	print "number of representative transcripts: " + str( len( transcripts.keys() ) )
	
	# --- constructing file with corresponding mRNAs --- #
	construct_codingseqs( transcripts, transcript_file, ref_seq )
	
	# --- screen all transcripts --- #
	normal_splice_events = 0
	u12_splice_events = 0
	freq_deviation_event = 0
	other_splice_events = 0	
	too_short_intron_events = 0
	
	for mRNA in transcripts.keys():
		raw_exons = sorted( transcripts[ mRNA ], key=itemgetter( 'start' ) )
		chromosome = raw_exons[0]['chr']
		exons = [ raw_exons[0] ]
		if len( raw_exons ) > 1:
			for each in raw_exons[1:]:
				if each['start'] < exons[-1]['end']:
					if each['end'] > exons[-1]['end']:
						del exons[-1]
						exons.append( each )
				else:
					exons.append( each )
		
		if exons[0]['orientation'] == '+':	#transcript on forward strand
			chr_seq = ref_seq[ chromosome ]
			for idx, exon in enumerate( exons[:-1] ):
				splice_donor = chr_seq[ exon['end']:exon['end']+2 ].upper()
				splice_acceptor = chr_seq[ exons[idx+1]['start']-3:exons[idx+1]['start']-1 ].upper()
				normal_counter, u12_counter, freq_deviation, other_counter , short_counter = check_splice_sites( splice_donor, splice_acceptor, mRNA, idx, exon['end'], exons[idx+1]['start']-3, out, exon['orientation'] )
				normal_splice_events += normal_counter
				u12_splice_events += u12_counter
				freq_deviation_event += freq_deviation
				other_splice_events += other_counter
				too_short_intron_events += short_counter
		elif exons[0]['orientation'] == "-":	#transcript on reverse strand
			exons = exons[::-1]
			chr_seq = ref_seq[ chromosome ]
			for idx, exon in enumerate( exons[:-1] ):
				splice_donor = revcomp( chr_seq[ exon['start']-3:exon['start']-1 ] ).upper()
				splice_acceptor = revcomp( chr_seq[ exons[idx+1]['end']:exons[idx+1]['end']+2 ] ).upper()
				normal_counter, u12_counter, freq_deviation, other_counter, short_counter = check_splice_sites( splice_donor, splice_acceptor, mRNA, idx, exon['start']-3, exons[idx+1]['end'], out, exon['orientation'] )
				normal_splice_events += normal_counter
				u12_splice_events += u12_counter
				freq_deviation_event += freq_deviation
				other_splice_events += other_counter
				too_short_intron_events += short_counter
		else:
			print "ERROR: orientation! "
	
	out.write( "number of normal splice events: " + str( normal_splice_events ) + '\n' )
	out.write( "number of U12 splice events: " + str( u12_splice_events ) + '\n' )
	out.write( "number of frequent deviation type events: " + str( freq_deviation_event ) + '\n' )
	out.write( "number of other splice events: " + str( other_splice_events ) + '\n' )
	return rna_to_gene


def check_splice_sites( donor_site, acceptor_site, gene_ID, exon_idx, d_pos, a_pos, out, exon_orientation ):
	"""! @brief check the donor and acceptor splice site """
	
	min_intron_size = 20
	
	# --- these positions are indices to extract coverage values from a list of coverage values per chromosome --- #
	if exon_orientation == "+":
		ex3_pos = str( d_pos-2 ) + "_" + str( d_pos )
		in5_pos = str( d_pos ) + "_" + str( d_pos+2 )
		in3_pos = str( a_pos ) + "_" + str( a_pos+2 )
		ex5_pos = str(a_pos+2 ) + "_" + str( a_pos+4 )
		intron_size = a_pos+2-d_pos
	else:
		ex3_pos = str( d_pos+2 ) + "_" + str( d_pos+4 )
		in5_pos = str( d_pos ) + "_" + str( d_pos+2 )
		in3_pos = str( a_pos ) + "_" + str( a_pos+2 )
		ex5_pos = str( a_pos-2 ) + "_" + str( a_pos )
		intron_size = d_pos+2-a_pos
	
	if intron_size >= min_intron_size:
		# --- classify splice site as canonical or non-canonical --- #
		if ( donor_site == "GT" ) + ( acceptor_site == "AG" ) == 2:
			out.write( gene_ID + "\tGT\tAG\tcss\t" + ",".join( [ ex3_pos, in5_pos, in3_pos, ex5_pos ] ) + '\t' + str( intron_size ) + "\n" )
			return ( 1, 0, 0, 0, 0 )
		elif ( donor_site == "AT" ) + ( acceptor_site == "AC" ) == 2:
			out.write( gene_ID + "\tAT\tAC\tncss\t" + ",".join( [ ex3_pos, in5_pos, in3_pos, ex5_pos ] ) + '\t' + str( intron_size ) + "\n")
			return ( 0, 1, 0, 0, 0 )
		elif ( donor_site == "GC" ) + ( acceptor_site == "AG" ) == 2:
			out.write( gene_ID + "\tGC\tAG\tncss\t" + ",".join( [ ex3_pos, in5_pos, in3_pos, ex5_pos ] ) + '\t' + str( intron_size ) + "\n" )
			return ( 0, 1, 1, 0, 0 )
		else:
			out.write( gene_ID + '\t' + donor_site + '\t' + acceptor_site + "\tncss\t" + ",".join( [ ex3_pos, in5_pos, in3_pos, ex5_pos ] ) + '\t' + str( intron_size ) + "\n" )
			return ( 0, 0, 0, 1, 0 )
	else:
		return ( 0, 0, 0, 0, 1 )


def final_analysis( doc_file, rna_to_gene ):
	"""! @brief run statistic analysis on output of first function """
	
	# --- analyze results --- #	
	data = []
	ncss_genes = []
	with open( doc_file, "r" ) as f:
		line = f.readline()
		while line:
			if not "number" in line:
				parts = line.strip().split('\t')
				data.append( { 'gene_id': rna_to_gene[ parts[0] ], 'splice_site': parts[1]+"..."+parts[2] } )
				if "ncss" in line and not "N" in parts[1]+parts[2]:
					ncss_genes.append( parts[0] )
			line = f.readline()
	splice_sites = []
	N_splice_site_counter = 0
	for each in data:
		if not "N" in each['splice_site']:
			splice_sites.append( each['splice_site'] )
		else:
			N_splice_site_counter += 1
	
	unique_splice_sites = list( set( splice_sites ) )
	occurences = [ 0 ] * len( unique_splice_sites )
	for each in data:
		try:
			occurences[ unique_splice_sites.index( each['splice_site'] ) ] += 1
		except ValueError:	#splice sites containing Ns are filtered out
			pass
	
	with open( doc_file, "a" ) as out:
		out.write("\n\n\nnumber of different splice site types: " + str( len( unique_splice_sites ) ) + '\n' )
		out.write( "total number of spliced introns: " + str( len( data ) ) + '\n' )
		out.write( "number of splice sites with N: " + str( N_splice_site_counter ) + '\n' )
		for idx, each in enumerate( unique_splice_sites ):
			#print each + "\t" + str( occurences[ idx ] )
			out.write( each + '\t' + str( occurences[ idx ] ) + '\n' )
		out.write( "\n\n\n# --- genes with non-canonical splice sites --- #\n" )
		out.write( "\n".join( sorted( list( set( ncss_genes ) ) ) ) + '\n' )


def load_genetic_code():
	"""! @brief return standard genetic code """
	
	genetic_code = {	'CTT': 'L',
									'ATG': 'M',
									'AAG': 'K',
									'AAA': 'K',
									'ATC': 'I',
									'AAC': 'N',
									'ATA': 'I',
									'AGG': 'R',
									'CCT': 'P',
									'ACT': 'T',
									'AGC': 'S',
									'ACA': 'T',
									'AGA': 'R',
									'CAT': 'H',
									'AAT': 'N',
									'ATT': 'I',
									'CTG': 'L',
									'CTA': 'L',
									'CTC': 'L',
									'CAC': 'H',
									'ACG': 'T',
									'CCG': 'P',
									'AGT': 'S',
									'CAG': 'Q',
									'CAA': 'Q',
									'CCC': 'P',
									'TAG': '*',
									'TAT': 'Y',
									'GGT': 'G',
									'TGT': 'C',
									'CGA': 'R',
									'CCA': 'P',
									'TCT': 'S',
									'GAT': 'D',
									'CGG': 'R',
									'TTT': 'F',
									'TGC': 'C',
									'GGG': 'G',
									'TGA': '*',
									'GGA': 'G',
									'TGG': 'W',
									'GGC': 'G',
									'TAC': 'Y',
									'GAG': 'E',
									'TCG': 'S',
									'TTA': 'L',
									'GAC': 'D',
									'TCC': 'S',
									'GAA': 'E',
									'TCA': 'S',
									'GCA': 'A',
									'GTA': 'V',
									'GCC': 'A',
									'GTC': 'V',
									'GCG': 'A',
									'GTG': 'V',
									'TTC': 'F',
									'GTT': 'V',
									'GCT': 'A',
									'ACC': 'T',
									'TTG': 'L',
									'CGT': 'R',
									'TAA': '*',
									'CGC': 'R'
								}
	return genetic_code


def translate( seq, genetic_code ):
	"""! @brief translates the given nucleotide sequence into peptide and splits at each star (stop codon) """
	
	seq = seq.upper()
	
	peptide = []
	
	for i in range( int( len( seq ) / 3.0 ) ):
		codon = seq[i*3:i*3+3]
		try:
			peptide.append( genetic_code[ codon ] )
		except:
			peptide.append( "*" )
	return "".join( peptide )


def construct_pep_file( genetic_code, CDS_file, pep_file ):
	"""! @brief construct peptide sequences for all constructed CDS entries """
	
	CDS = load_sequences( CDS_file )
	with open( pep_file, "w" ) as out:
		for key in CDS.keys():
			out.write( '>' + key + '\n' + translate( CDS[ key ], genetic_code ) + '\n' )


def main( arguments ):
	"""! @brief controls processing of multiple files """
	
	input_dir = arguments[ arguments.index('--in')+1 ]
	if input_dir[-1] != "/":
		input_dir += "/"
	
	fasta_files = glob.glob( input_dir + "*.fna" )
	print "number of : " + str( len( fasta_files ) )
	for ref_seq_file in sorted( fasta_files ):
		ID = ref_seq_file.split('/')[-1].split('.')[0]
		try:
			print "processing .... " + ID
			doc_file = ref_seq_file.replace( ".fna", ".txt" )
			gff_file = ref_seq_file.replace( ".fna", ".gff" )
			transcript_file = ref_seq_file.replace( ".fna", ".mRNA.fasta" )
			pep_file = ref_seq_file.replace( ".fna", ".pep.fa" )
			
			# -- loading sequences --- #
			ref_seq = load_sequences( ref_seq_file )
			
			# --- checking for intron borders --- #
			with open( doc_file, "w" ) as out:
				rna_to_gene = get_all_intron_borders_per_gene( gff_file, ref_seq, out, transcript_file )
			
			g2t_file = ref_seq_file.replace( ".fna", ".g2t" )
			with open( g2t_file, "w" ) as g2t:
				for key in rna_to_gene.keys():
					g2t.write( key + '\t' + rna_to_gene[ key ] + '\n' )
			
			# --- calculating stats of different intron borders --- #
			final_analysis( doc_file, rna_to_gene )
			
			# ---- translate extracted CDS --- #
			genetic_code = load_genetic_code()
			construct_pep_file( genetic_code, transcript_file, pep_file )
			
		except:
			print "ERROR: " + ID

if __name__ == '__main__':
	
	if '--in' in sys.argv:
		main( sys.argv )
	else:
		sys.exit( __usage__ )

