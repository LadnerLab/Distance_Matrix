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

    ymer_dict = {}
    for index in range( num_seqs ):
        ymer_dict[ sequences[ index ] ] = oligo.subset_lists_iter( sequences[ index ], options.XmerWindowSize, 1 )

    num_seqs = len( sequences )
    out_list = np.empty( 1 )

    for current_seq in range( num_seqs ):
        for inner_index in range( current_seq + 1, num_seqs ):
            out_list = np.append( out_list,
                                  oligo.get_single_sequence_dist( ymer_dict[ sequences[ current_seq ] ],
                                                                  ymer_dict[ sequences[ inner_index ] ], options.XmerWindowSize, 1
                                                                )
                                )
            
    out_list = np.delete( out_list, 0 )

    if options.verbose:
        print( "Distance Matrix complete" )

    Z = linkage( out_list, 'single' )

    cluster = fcluster( Z, options.clusters, criterion ='maxclust' )

    out_file = open( options.output, 'w' )

    if options.verbose:
        print( "Clustering Complete" )
    
    for sequence in range( len( names ) ):
        if cluster[ sequence ] not in cluster_dict:
            cluster_dict[ cluster[ sequence ] ] = list()
        cluster_dict[ cluster[ sequence ] ].append( ( names[ sequence ], sequences[ sequence ] ) )
            
        out_file.write( "%d %s\n" % ( cluster[ sequence ], names[ sequence ] ) )

    if options.verbose:
        display_cluster_information( cluster_dict, out_list, options.XmerWindowSize, 1, ymer_dict )

    out_file.close()

                                        

def display_cluster_information( cluster_dict, list_of_distances, window_size, step_size, ymer_dict = None ):

    dict_values = cluster_dict.values()

    clusters_total = 0
    cluster_seqs = 0

    avg_distance = 0 
    max_distance = 0
    min_distance = len( dict_values )

    # Viral stats
    avg_species_per_cluster = 0
    num_species = 0
    species_from_sequences = set()
    species_per_cluster = {}

    

    # Cluster stats
    num_clusters = len( cluster_dict.keys() )
    min_cluster_size = min( [ len( item ) for item in dict_values ] )
    max_cluster_size = max( [ len( item ) for item in dict_values ] )
    avg_cluster_size = sum( [ len( item ) for item in dict_values ] ) / num_clusters

    # Distance stats
    for key, item in cluster_dict.items():
        names = [ seq[ 0 ] for seq in item ]
        for current_name in names:
            id = get_taxid_from_name( current_name )
            species_from_sequences.add( id )

            if key not in species_per_cluster:
                species_per_cluster[ key ] = set()
            species_per_cluster[ key ].add( id )

        current_matrix = oligo.create_distance_matrix_of_sequences( [ seq[ 1 ] for seq in item ], window_size,
                                                                    step_size, ymer_dict 

                                                                  )
        if len( item ) > 1:
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

    num_species = len( species_from_sequences )
    avg_distance = clusters_total / cluster_seqs
    avg_species_per_cluster = ( sum( [ len( item ) for item in species_per_cluster.values() ] ) / num_clusters )
    avg_cluster_per_species = num_clusters / num_species
    
    print( "\nNumber of clusters: %d." % num_clusters )
    print( "Minimum Cluster Size: %.2f." % min_cluster_size )
    print( "Maximum Cluster Size: %.2f." % max_cluster_size )
    print( "Average Cluster Size: %.2f.\n" % avg_cluster_size )

    print( "Minimum distance between any two sequences within the clusters: %.2f" % min_distance )
    print( "Average distance between any two sequences within the clusters: %.2f" % avg_distance )
    print( "Maximum distance between any two sequences within the clusters: %.2f\n" % max_distance )

    print( "Number of species found in file: %d" % num_species )
    print( "Average species per cluster: %.2f" % avg_species_per_cluster )
    print( "Average clusters per species: %.2f" % avg_cluster_per_species )

    

def get_taxid_from_name( in_name ):
    split_name = in_name.split( 'TaxID=' )

    return_name = None
    if len( split_name ) > 1:
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

    options.add_option( '-v', help = "Display statistical output of clusters, disabled by default because this is very slow. [False]",
                        action = 'store_true', dest = 'verbose'
                      )  
    
if __name__ == '__main__':
    main()
