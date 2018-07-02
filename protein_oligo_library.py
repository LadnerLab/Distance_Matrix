#!/usr/bin/env python3
   
def read_fasta_lists( file_to_read ):
    """
       Reads a list of fastas from file_to_read
    
       Returns:
        names- a list of names of the sequences found in the fasta file
        sequences- a list of the sequences found in the fasta file
    """

    file_in = open( file_to_read, 'r' )
    count = 0

    names = []
    sequences = []
    current_sequence = ''

    for line in file_in:
        line = line.strip()
        if line and line[ 0 ] == '>':
            count += 1
            names.append( line[ 1: ] )
            if count > 1:
                sequences.append( current_sequence )
            current_sequence = ''

        else:
            current_sequence += line

    sequences.append( current_sequence )
    file_in.close()

    return names, sequences
 
def write_fastas( names_list, sequence_list, output_name="out.txt" ):
    """
        Writes a fasta file from a list of names and sequences to output file provided

    """
    out_file = open( output_name, 'w+' )
    for index in range( len( names_list ) ):
        out_file.write( '>' + names_list[ index ] + '\n' + 
                        sequence_list[ index ] + '\n'
                      )
    out_file.close()
        

def char_in_string( test_string, character ):
    """
        Checks if a character is found within a given string
    """
    for index in range( len( test_string ) ):
        if test_string[ index ] == character:
            return True
    return False

def percentage_of_char_in_string( test_string, character ):
    """
        Calculates what percent of a given test_string is given character
    
        Params:
           test_string- string to test for character
           character- character to test for in string
       Returns:
           floating point value percent of the string that
           is character
    """
    length = len( test_string )
    char_count = 0.0

    for index in range( length ):
        if test_string[ index ] == character:
            char_count += 1
    return ( char_count / length ) * 100

def count_char_in_string( test_string, character ):
    """
        Counts how much of a certain character are in test_string

        Params:
            test_string- string to count
            character- character to count in string
        Returns:
            integer value representing the number of character 
                were found in test_string
    """
    length = len( test_string )
    count = 0

    for index in range( length ):
        if test_string[ index ] == character:
            count += 1
    return count

def min_concurrent_chars( test_string, delimeter_char ):
    """
        Finds the minimum number of concurrent non-delimiter char in 
        test_string
    
        Params:
          test_string- string to search
          delimeter_char- character to reset the count
        Returns:
          integer value, the smallest amount of concurrent characters
               between delimeter character
    """

    split_string = test_string.split( delimeter_char )
    min_length = len( split_string[ 0 ] )

    for substring in split_string[ 1: ]:
        current_length = len( substring )
        if current_length > 0 and current_length < min_length:
            min_length = current_length

    return min_length


def remove_char_from_string( test_string, to_remove ):
    """
        Removes character to_remove from string test_string
        Note: Case sensitive method, 'a' != 'A'
        Returns:
           test_string, minus any instance of to_remove character
    """
    output = ""
    for index in range( len( test_string ) ):
        if test_string[ index ] != to_remove:
            output += test_string[ index ]
    return output


def create_list_of_uniques( names, sequences ):
    """
       Removes duplicates from a list
       Params:
          names- names of the sequences
       sequences:
          a list of sequences who may or may not be unique
    """
    return_names = []
    return_sequences = []

    unique_values = set()

    for index in range( len( sequences ) ):
        starting_length = len( unique_values )
        unique_values.add( sequences[ index ] )

        if len( unique_values ) > starting_length:
            return_names.append( names[ index ] )
            return_sequences.append( sequences[ index ] )
    return return_names, return_sequences

def create_valid_sequence_list( names_list, sequence_list, min_length, percent_valid, start, end ):
   """
       Creates a sequence list of valid sequences.
       A valid sequence is defined by not having any 'X' characters,
       and not violating the parameters of either min_length or percent_valid 
       
       Returns:
           a list of names of those sequences that were valid, with the new bounds appended to the name
           a list of the sequences that were found valid
   """
   valid_names = []
   valid_sequences = []

   for sequence in range( len( sequence_list ) ):
      current_sequence = sequence_list[ sequence ][ start : end ]

      if is_valid_sequence( current_sequence, min_length, percent_valid ):
           valid_names.append( names_list[ sequence ] )
           current_sequence = remove_char_from_string( current_sequence, '-' )
           valid_sequences.append( current_sequence[ start : end ] )

   names_list = append_suffix( valid_names, start, end )

   return names_list, valid_sequences

def is_valid_sequence( sequence, min_length, percent_valid ):
   """
       Determines whether a given sequence is valid 
       A valid sequence is defined by not having any 'X' characters,
           and not violating the parameters of either min_length or percent_valid 
   """
   if not char_in_string( sequence, 'X' ):
       if min_length is None:
           return percentage_of_char_in_string( sequence, '-' ) < ( 100 - percent_valid )
       else:
           return ( min_concurrent_chars( sequence, '-' ) >= min_length )
   return False

         
def append_suffix( string, start, end ):
   """
       Appends _start_end to a string
   """
   return "%s_%s_%s" % ( string, str( start ), str( end ) ) 


def subset_lists_iter( sequence, window_size, step_size ):
    xmer_set= set()

    start = 0
    end = window_size

    while end <= len( sequence ):
        xmer = sequence[ start:end ]

        if not 'X' in xmer:
            xmer_set.add( xmer )

        start += step_size
        end = start + step_size + window_size - 1

    return xmer_set

def subset_lists( name, sequence, window_size, step_size ):
   """
       Creates a list of subsets of windowSize size in intervals of stepSize
       Note: Uses recursive subset_lists_helper for operations
   
       Params:
            name: String name of sequence to be split up
            sequence: String sequence to be split up into a list
       Returns:
            a list of the split up names, with a suffix applied, and a list of the segments 
            of the list, as specified
   """
   new_names = []
   new_seqs = []
   return subset_lists_helper( name, sequence, new_names, new_seqs, window_size, step_size, 0, window_size )

def subset_lists_helper( name, sequence, name_arr, seq_arr, window_size, step_size, start, end ):
    """
        Recursive helper method called by subset_lists
    """
    if start + window_size < len( sequence ):
       if len( sequence[ start: end ] ) == 1:
           return
       seq_arr.append( sequence[ start : end ] ) 
       name_arr.append( append_suffix( name, start + 1, end ) )

       subset_lists_helper( name, sequence, name_arr, seq_arr, window_size, step_size, start + step_size, start + step_size + window_size )
    return name_arr, seq_arr
   


def component_xmer_locs(ymer, xdict, xmer_size, step_size):
    xmer_locs = []
    ymer_xdict = subset_lists_iter( '', ymer, xmer_size, step_size)
    for each in ymer_xdict.keys():
        xmer_locs += xdict[each]
    return set(xmer_locs)
    
def get_single_sequence_dist( first_seq_ymers, second_seq_ymers, window_size, step_size, return_subset_ymers = False ):
    """
        Computes the 'distance' between two sequences, where distance is the 
        number of shared xmers / the average of the xmers of both sequences

        :param first_seq_ymers: list of ymers to compare to second_seq_ymers
        :param second_seq_ymers: list of ymers to compare to first_seq_ymers
        :param return_subset_ymers: boolean flag to return the set of a sequence's ymers
       
        :returns: percentage of average ymers shared between between sequence first_seq and sequence second_seq,
                  and also first_seq's subset ymers, if return_subset_ymers is set to true
    """
    intersection = len( first_seq_ymers & second_seq_ymers )
    average_length = min( len( first_seq_ymers ), len( second_seq_ymers ) )

    if average_length == 0:
        average_length = max( len( first_seq_ymers ), len( second_seq_ymers ) )
        

    if return_subset_ymers:
        return ( ( intersection / average_length ) * 100 ), first_seq_ymers
    return 100 - ( ( intersection / average_length ) * 100 )


def get_distance_from_other_sequences( in_seq, sequence_list, ymer_dict, window_size, step_size, return_ymer_list = False ):
    """
        Computes the distance from one sequence to all of the sequence in sequence_list,
        where distance is defined as in get_single_sequence_dist
    
        :param in_seq: string sequence to be compared to all sequences in sequence_list
        :param sequence_list: list of sequences to compare to in_seq
        :param window_size: integer window_size to capture in subset_lists
        :param step_size: integer step_size to use in subset_lists
        :param ymer_dict: dictionary of sequence: list of ymers
    
        :returns: list of distances from in_seq to all of the sequences in sequence_list

    """
    return_difference_list = list()

    # in_seq_ymers = subset_lists_iter( in_seq, window_size, step_size )
    in_seq_ymers = ymer_dict[ in_seq ]
    for current_seq in sequence_list:
        # current_seq_ymers = subset_lists_iter( current_seq, window_size, step_size )
        current_seq_ymers = ymer_dict[ current_seq ]
        if return_ymer_list:
            difference, ymers = get_single_sequence_dist( in_seq_ymers, current_seq_ymers, window_size, step_size, return_ymer_list )
        else:
            difference = get_single_sequence_dist( in_seq_ymers, current_seq_ymers, window_size, step_size, return_ymer_list )


        return_difference_list.append( difference )

    return return_difference_list

def create_distance_matrix_of_sequences( sequence_list, window_size, step_size, ymer_dict = None ):
    """
        Create a distance matrix from a list of sequences, where distance is defined as in 
        get_single_sequence_dist

        :param sequence_list: list of string sequences from which to construct the distance matrix
        :param window_size: integer window_size to capture at a time from each sequence
        :param step_size: integer step_size to use in subset_lists
        :param ymer_dict: dictionary of sequence: list of ymers, will be created from sequence_list
                          if none is provided
    

        :returns: matrix of sequence distances, adjacency matrix of a graph
    """
    distance_matrix = list()

    if ymer_dict is None:
       ymer_dict = {}
       for current_sequence in range( len( sequence_list) ):
           ymer_dict[ sequence_list[ current_sequence ] ] = \
                   subset_lists_iter( sequence_list[ current_sequence ], window_size, step_size )

    
    for current_sequence in range( len( sequence_list) ):
        distance_matrix.append(
                               get_distance_from_other_sequences(
                                                                 sequence_list[ current_sequence ],
                                                                 sequence_list, ymer_dict, window_size, step_size
                                                                 )
                              )
    return distance_matrix

        


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

def get_taxid_from_name( in_name ):
    """
        Gets the taxonomic id from a sequence name
        Note: Capable of parsing uniref and uniprot types

        :param in_name: string name of genetic sequence
        :returns: string taxonomic id of sequence
    """

    # uniref naming convention
    split_name = in_name.split( 'TaxID=' )

    # uniprot naming convention
    if len( split_name == 1 ):
        split_name = in_name.split( 'OX=' )

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
