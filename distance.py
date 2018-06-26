#!/usr/bin/env python3
import sys
import optparse

import protein_oligo_library as oligo

def main():

    usage = "usage: %prog [options]"
    option_parser = optparse.OptionParser( usage )
    add_program_options( option_parser )
    options, arguments = option_parser.parse_args()


    if options.query is None:
        print( "ERROR: Fasta query file must be provided." )
        sys.exit( 0 )

    names, sequences = oligo.read_fasta_lists( options.query )
    names, sequences = oligo.sort_sequences_by_length( names, sequences )

    # Get the sequences sorted in increasing order
    names.reverse()
    sequences.reverse()

    out_file = open( options.output, 'w' )
    
    for index in range( len( sequences ) ):
        current_distances = oligo.get_distance_from_other_sequences( sequences[ index ],
                                                               sequences,
                                                               options.XmerWindowSize, 1
                                                             )
        for item in range( len( current_distances ) ):
            first_seq_name = ''.join( names[ index ].split()[ 0 ] ) 
            second_seq_name = ''.join( names[ item ].split()[ 0 ] ) 
            out_file.write( "%s\t%s\t%f\n" % ( first_seq_name, second_seq_name, 100 - current_distances[ item ] ) )
            

    out_file.close()        

def add_program_options( options ):
    options.add_option( '-q', '--query', help = "Fasta query file to perform calculations on. [None, required]" )

    # TODO: Add description of format of output file  
    options.add_option( '-o', '--output', help = "File to write program output to. [out.txt]", default = "out.txt" )

    options.add_option( '-x', '--XmerWindowSize', help = "Size of xmers to grab from each sequence to do the comparisons [19]", type = int,
                        default = 19 )


    
if __name__ == '__main__':
    main()

