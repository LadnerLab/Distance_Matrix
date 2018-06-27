#!/usr/bin/env python3

from matplotlib import pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
import numpy as np 
import optparse
import sys
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

    cluster_dict = {}

    # Get the sequences sorted in decreasing order
    names.reverse()
    sequences.reverse()

    num_seqs = len( sequences )
    out_list = np.empty( 1 )
    for current_seq in range( num_seqs ):
        for inner_index in range( current_seq + 1, num_seqs ):
            out_list = np.append( out_list, oligo.get_single_sequence_dist( sequences[ current_seq ], sequences[ inner_index ], options.XmerWindowSize, 1 ) )
            
    out_list = np.delete( out_list, 0 )

    Z = linkage( out_list, 'single' )

    cluster = fcluster( Z, options.clusters, criterion ='maxclust')

    out_file = open( options.output, 'w' )
    
    for sequence in range( len( names ) ):
        if cluster[ sequence ] not in cluster_dict:
            cluster_dict[ cluster[ sequence ] ] = list()
        cluster_dict[ cluster[ sequence ] ].append( sequences[ sequence ] )
            
        out_file.write( "%d %s\n" % ( cluster[ sequence ], names[ sequence ] ) )

    display_cluster_information( cluster_dict )

    out_file.close()

                                        

def display_cluster_information( cluster_dict ):
    num_clusters = len( cluster_dict.keys() )
    min_cluster_size = min( [ len( item ) for item in cluster_dict.values() ] )
    max_cluster_size = max( [ len( item ) for item in cluster_dict.values() ] )
    avg_cluster_size = sum( [ len( item ) for item in cluster_dict.values() ] ) 

    print( "Number of clusters: %d." % num_clusters )
    print( "Minimum Cluster Size: %d." % min_cluster_size )
    print( "Maximum Cluster Size: %d." % max_cluster_size )
    print( "Average Cluster Size: %d." % avg_cluster_size )
    


def add_program_options( options ):
    options.add_option( '-q', '--query', help = "Fasta query file to perform calculations on. [None, required]" )

    # TODO: Add description of format of output file  
    options.add_option( '-o', '--output', help = "File to write program output to. [out.txt]", default = "out.txt" )
    options.add_option( '-c', '--clusters', help = "Maximum number of clusters to produce in output. [4]", default = 4, type = int )

    options.add_option( '-x', '--XmerWindowSize', help = "Size of xmers to grab from each sequence to do the comparisons [19]", type = int,
                        default = 19 )


    
if __name__ == '__main__':
    main()
