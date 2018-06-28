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
            out_list = np.append( out_list,
                                  oligo.get_single_sequence_dist( sequences[ current_seq ],
                                                                  sequences[ inner_index ], options.XmerWindowSize, 1
                                                                )
                                )
            
    out_list = np.delete( out_list, 0 )

    Z = linkage( out_list, 'single' )

    cluster = fcluster( Z, options.clusters, criterion ='maxclust')

    out_file = open( options.output, 'w' )
    
    for sequence in range( len( names ) ):
        if cluster[ sequence ] not in cluster_dict:
            cluster_dict[ cluster[ sequence ] ] = list()
        cluster_dict[ cluster[ sequence ] ].append( ( names[ sequence ], sequences[ sequence ] ) )
            
        out_file.write( "%d %s\n" % ( cluster[ sequence ], names[ sequence ] ) )

    display_cluster_information( cluster_dict, out_list, options.XmerWindowSize, 1 )

    out_file.close()

                                        

def display_cluster_information( cluster_dict, list_of_distances, window_size, step_size ):

    dict_values = cluster_dict.values()

    clusters_total = 0
    cluster_seqs = 0

    avg_distance = 0 
    max_distance = 0
    min_distance = len( dict_values )

    # Cluster stats
    num_clusters = len( cluster_dict.keys() )
    min_cluster_size = min( [ len( item ) for item in dict_values ] )
    max_cluster_size = max( [ len( item ) for item in dict_values ] )
    avg_cluster_size = sum( [ len( item ) for item in dict_values ] ) / num_clusters

    # Distance stats
    for item in dict_values:
        if len( item ) > 1:
            current_matrix = oligo.create_distance_matrix_of_sequences( [ seq[ 1 ] for seq in item ], window_size,
                                                                        step_size

                                                                      )

            matrix_array = list()
            for current_distance in range( len( current_matrix ) ):
                for local_distance in range( current_distance + 1, len( current_matrix ) ):
                    matrix_array.append( current_matrix[ current_distance ][ local_distance ] )
                    
            clusters_total += sum( matrix_array )
            cluster_seqs += len( matrix_array )

            local_max = max( matrix_array )
            local_min = min( matrix_array )

            max_distance = max( local_max, max_distance )
            min_distance = min( local_min, min_distance )

    avg_distance = clusters_total / cluster_seqs
    
    print( "Number of clusters: %d." % num_clusters )
    print( "Minimum Cluster Size: %.2f." % min_cluster_size )
    print( "Maximum Cluster Size: %.2f." % max_cluster_size )
    print( "Average Cluster Size: %.2f." % avg_cluster_size )

    print( "Minimum distance between any two sequences within the clusters: %.2f" % min_distance )
    print( "Average distance between any two sequences within the clusters: %.2f" % avg_distance )
    print( "Maximum distance between any two sequences within the clusters: %.2f" % max_distance )

    get_species_from_file( 'names.dmp' )
    

def get_taxid_from_name( in_name ):
    split_name = in_name.split( 'TaxID=' )

    return_name = None
    if len( in_name ) > 1:
        split_name = split_name[ 1 ].split()
        return_name = split_name[ 0 ]
    return return_name

def get_species_from_file( in_file ):
    open_file = open( in_file, 'r' )
    taxid_dict = {}

    for line in open_file:
        line = line.split( '|' )
        taxID = line[ 0 ].strip()
        species = line[ 1 ].strip()

        if taxID not in taxid_dict:
            taxid_dict[ int( taxID ) ] = list()
        taxid_dict[ taxID ].append( species )
    open_file.close()

    return taxid_dict

def add_program_options( options ):
    options.add_option( '-q', '--query', help = "Fasta query file to perform calculations on. [None, required]" )

    options.add_option( '-o', '--output',
                        help = "File to write program output to. Output is a tab-delimited file containing cluster number, and then sequence name [out.txt]",
                        default = "out.txt"
                      )
    options.add_option( '-c', '--clusters', help = "Maximum number of clusters to produce in output. [4]", default = 4, type = int )

    options.add_option( '-x', '--XmerWindowSize', help = "Size of xmers to grab from each sequence to do the comparisons [19]", type = int,
                        default = 19 )


    
if __name__ == '__main__':
    main()
