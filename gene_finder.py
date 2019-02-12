# -*- coding: utf-8 -*-
"""
gene_finder_mini_project

@author: Anne.j
"""

import random
from amino_acids import aa, codons, aa_table   # you may find these useful
from load import load_seq


def shuffle_string(s):
    """Shuffles the characters in the input string
        NOTE: this is a helper function, you do not
        have to modify this in any way
    """
    return ''.join(random.sample(s, len(s)))

# YOU WILL START YOUR IMPLEMENTATION FROM HERE DOWN ###


def get_complement(nucleotide):
    """ Returns the complementary nucleotide

        nucleotide: a nucleotide (A, C, G, or T) represented as a string
        returns: the complementary nucleotide
    >>> get_complement('A')
    'T'
    >>> get_complement('C')
    'G'
    """

    if nucleotide == 'A':
        return 'T'
    elif nucleotide == 'T':
        return 'A'
    elif nucleotide == 'G':
        return 'C'
    elif nucleotide == 'C':
        return 'G'


def get_reverse_complement(dna):
    """ Computes the reverse complementary sequence of DNA for the specfied DNA
        sequence

        dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as a string
        tl;dr first flip, then find complement
    >>> get_reverse_complement("ATGCCCGCTTT")
    'AAAGCGGGCAT'
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'
    """

    dna_reversed = dna[::-1]
    reversed_complement=''

    for nucleotide in dna_reversed:
        reversed_complement += get_complement(nucleotide)
    return reversed_complement


def rest_of_ORF(dna):
    """ Takes a DNA sequence that is assumed to begin with a start
        codon and returns the sequence up to but not including the
        first in frame stop codon.  If there is no in frame stop codon,
        returns the whole string.

        dna: a DNA sequence
        returns: the open reading frame represented as a string

    >>> rest_of_ORF("ATGTGAA")
    'ATG'
    >>> rest_of_ORF("ATGAGATAGG")
    'ATGAGA'
    >>> rest_of_ORF("TAGTAGTAG")
    ''
    >>> rest_of_ORF("ATGCCCGGGAAA")
    'ATGCCCGGGAAA'
    >>> rest_of_ORF("")
    """

    stop_codon= ['TAG', 'TAA', 'TGA']
    for i in range(0, len(dna), 3): #the for loop
        if dna[i:i+3] in stop_codon: #not the for loop, and is a conditional statement for cutting sequence UP TO stop_codon
            return(dna[:i])
    return dna


def find_all_ORFs_oneframe(dna):
    """ Finds all non-nested open reading frames in the given DNA
        sequence and returns them as a list. This function should
        only find ORFs that are in the default frame of the sequence
        (i.e. they start on indices that are multiples of 3).
        By non-nested we mean that if an ORF occurs entirely within
        another ORF, it should not be included in the returned list of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")
    ['ATGCATGAATGTAGA', 'ATGTGCCC']
    >>> find_all_ORFs_oneframe("TAGTAGTAGTAG")
    ['', '', '', '']
    >>> find_all_ORFs_oneframe("")
    none
    """

    list_of_finalized_protien_strings=[]
    i=0
    start='ATG'

    while i< len(dna):
        a = dna[i:i+3]
        if a == start:
            placestartcut=dna[i:]
            placeendcut=rest_of_ORF(placestartcut)
            list_of_finalized_protien_strings.append(placeendcut)
            i=i+len(placeendcut)
        else:
            i=i+3
    return list_of_finalized_protien_strings


def find_all_ORFs(dna):
    """
        Finds all non-nested open reading frames in the given DNA sequence in
        all 3 possible frames and returns them as a list.  By non-nested we mean
        that if an ORF occurs entirely within another ORF and they are both in
        the same frame, it *should not* be included in the returned list of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

    >>> find_all_ORFs("ATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']
    """

    newlist=[]

    zerothh= dna
    firstt= dna[1:]
    secondd= dna[2:]

    a= find_all_ORFs_oneframe(zerothh)
    b= find_all_ORFs_oneframe(firstt)
    c= find_all_ORFs_oneframe(secondd)

    newlist=a+ b+ c
    return newlist


def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']

    """
    new_list=[]

    dna_strand_one= dna
    dna_strand_two=get_reverse_complement(dna)


    all_orf_of_strand_one=find_all_ORFs(dna_strand_one)
    all_orf_of_strand_two=find_all_ORFs(dna_strand_two)

    new_list=all_orf_of_strand_one+ all_orf_of_strand_two

    return new_list


def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    """
    a=find_all_ORFs_both_strands(dna)
    i=0
    max_length=0
    max_index=0 #stores location
    while i< len(a):
        if len(a[i])>max_length:
            max_length= len(a[i])
            max_index = i
        i += 1
    return a[max_index]


def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF
    """

    max_length=0

    for i in range(num_trials):
        k=shuffle_string(dna)
        longest_of_shuffling=longest_ORF(k)
        if len(longest_of_shuffling)>max_length:
            max_length=len(longest_of_shuffling)

        print(i)
    return max_length


def coding_strand_to_AA(dna):
    """ Computes the Protein encoded by a sequence of DNA.  This function
        does not check for start and stop codons (it assumes that the input
        DNA sequence represents an protein coding region).

        dna: a DNA sequence represented as a string
        returns: a string containing the sequence of amino acids encoded by the
                 the input DNA fragment

        >>> coding_strand_to_AA("ATGCGA")
        'MR'
        >>> coding_strand_to_AA("ATGCCCGCTTT")
        'MPA'
    """

    amino_acids=''
    for i in range(0, len(dna),3):
        if dna[i:i+3] in aa_table.keys():
            amino_acids += aa_table[dna[i:i+3]]
    return amino_acids


def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.

        tl;dr find longest thing in dna, convert it into amino acids
    """
    orf_orf =[]
    threshold = longest_ORF_noncoding(dna, 1500)

    print(threshold)
    #length_of_threshold=len(threshold)

    listoforfs=find_all_ORFs(dna)
    len_listoforfs=len(listoforfs)
    # print(len(listoforfs))
    i=0

    while i < len_listoforfs:
        if len(listoforfs[i])>threshold:
            orf_orf.append(coding_strand_to_AA(listoforfs[i]))
            # print("hi")
        i+=1
        #     convertedacids=coding_strand_to_AA(listoforfs[i])
        #     orf_orf.append(convertedacids)
        # orf_orf+=orf_orf.append(convertedacids)
    return orf_orf

    # print(threshold)

if __name__ == "__main__":
    import doctest
    # doctest.testmod()
    dna=load_seq("./data/X73525.fa")
    ok= gene_finder(dna)
    print(ok)
