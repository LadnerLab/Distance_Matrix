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
    names, sequences = sort_sequences_by_length( names, sequences )

    names.reverse()
    sequences.reverse()


def add_program_options( options ):
    options.add_option( '-q', '--query', help = "Fasta query file to perform calculations on. [None, required]" )

    # TODO: Add description of format of output file  
    options.add_option( '-o', '--output', help = "File to write program output to. [out.txt]", default = "out.txt" )

    options.add_option( '-x', '--XmerWindowSize', help = "Size of xmers to grab from each sequence to do the comparisons [19]", type = int,
                        default = 19 )


def sort_sequences_by_length( names_list, sequence_list ):
    sequence_dict = {}

    out_names = list()
    out_seqs = list()

    for index in range( len( names_list ) ):
        current_seq = sequence_list[ index ]
        current_name = names_list[ index ]

        sequence_dict[ current_seq ] = current_name

    out_seqs = sorted( list( sequence_dict.keys() ), key = len ) 

    for item in range( len( out_seqs ) ):
        current_item = out_seqs[ item ]
        out_names.append( sequence_dict[ current_item ] )

    return out_names, out_seqs





    
if __name__ == '__main__':
    main()

