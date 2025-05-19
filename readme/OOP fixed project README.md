## RENAME this file Sarsam_OOP_FinalProject_2023.py

##Assignment: Add to the constructor and methods of a parent class and child classes
##            which inherit the base class properties. NOTE: You are not allowed
##            to import any specialized libraries for this project (e.g., no Biopython)
##            The idea is for you to write these methods from scratch.

## Begin with the parent Seq class and the child DNA class we created in lecture below.
## 

### Seq Class
#
#  Constructor:
#  (1) Use the string functions upper and strip to clean up self.sequence.
#  (2) Add a variable self.kmers to the constructor and make it equal to an empty list.

#  Methods:
#  (1) Add a method called make_kmers that makes overlapping kmers of a given length from self.sequence
#      appends these to self.kmers. Default kmer parameter=3.
#  (2) Add a method called fasta that returns a fasta formatted string like this:
#      >species gene
#      AGATTGATAGATAGATAT


### DNA Class: INHERITS Seq class
#   
#  Constructor:
#  Use re.sub to change any non nucleotide characters in self.sequence into an 'N'.
#      re.sub('[^ATGCU]','N',sequence) will change any character that is not a
#      capital A, T, G, C or U into an N. (Seq already uppercases and strips.)

#  Methods:
#  (1) Add a method called print_info that is like print_record, but adds geneid and an
#      empty space to the beginning of the string.
#  (2) Add a method called reverse_complement that returns the reverse complement of
#      self.sequence
#  (3) Add a method called six_frames that returns all 6 frames of self.sequence
#      This include the 3 forward frames, and the 3 reverse complement frames

### RNA Class:  INHERITS DNA class
#  
#  Construtor:
#  Use the super() function (see DNA Class example).
#  (1) Automatically change all Ts to Us in self.sequence. 
#  (2) Add self.codons equals to an empty list

#  Methods:
#  (1) Add make_codons which breaks the self.sequence into 3 letter codons
#      and appends these codons to self.codons unless they are less than 3 letters long.
#  (2) Add translate which uses the Global Variable standard_code below to
#      translate the codons in self.codons and returns a protein sequence.

### Protein Class: INHERITS Seq class
#
#  Construtor:
#  Use the super() function (see DNA Class example).
#  Use re.sub to change any non LETTER characters in self.sequence into an 'X'.

#  Methods:
#  The next 2 methods use a kyte_doolittle and the aa_mol_weights dictionaries.
#  (2) Add total_hydro, which return the sum of the total hydrophobicity of a self.sequence
#  (3) Add mol_weight, which returns the total molecular weight of the protein
#      sequence assigned to the protein object.print("Sarsam_Fixed virsion with doctest")

import re  #first step, i import re to use later to clean the DNA and Protein sequence 
import doctest                #python doctest use to test the code by going over embedded examples in the docstring. 

standard_code = {
    "UUU": "F", "UUC": "F", "UUA": "L", "UUG": "L", "UCU": "S", "UCC": "S",
    "UCA": "S", "UCG": "S", "UAU": "Y", "UAC": "Y", "UAA": "*", "UAG": "*",
    "UGA": "*", "UGU": "C", "UGC": "C", "UGG": "W", "CUU": "L", "CUC": "L",
    "CUA": "L", "CUG": "L", "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAU": "H", "CAC": "H", "CAA": "Q", "CAG": "Q", "CGU": "R", "CGC": "R",
    "CGA": "R", "CGG": "R", "AUU": "I", "AUC": "I", "AUA": "I", "AUG": "M",
    "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T", "AAU": "N", "AAC": "N",
    "AAA": "K", "AAG": "K", "AGU": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GUU": "V", "GUC": "V", "GUA": "V", "GUG": "V", "GCU": "A", "GCC": "A",
    "GCA": "A", "GCG": "A", "GAU": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G"
}

kyte_doolittle = {
    'A': 1.8, 'C': 2.5, 'D': -3.5, 'E': -3.5, 'F': 2.8, 'G': -0.4, 'H': -3.2,
    'I': 4.5, 'K': -3.9, 'L': 3.8, 'M': 1.9, 'N': -3.5, 'P': -1.6, 'Q': -3.5,
    'R': -4.5, 'S': -0.8, 'T': -0.7, 'V': 4.2, 'W': -0.9, 'X': 0, 'Y': -1.3
}

aa_mol_weights = {
    'A': 89.09, 'C': 121.15, 'D': 133.1, 'E': 147.13, 'F': 165.19, 'G': 75.07,
    'H': 155.16, 'I': 131.17, 'K': 146.19, 'L': 131.17, 'M': 149.21, 'N': 132.12,
    'P': 115.13, 'Q': 146.15, 'R': 174.2, 'S': 105.09, 'T': 119.12, 'V': 117.15,
    'W': 204.23, 'X': 0, 'Y': 181.19
}


class Seq:    #making class called seq
    def __init__(self, sequence, gene, species, **kwargs):  #initializes the sequence their gene and species 
        """
        >>> s = Seq("acgtac", "my_gene", "Mouse")
        >>> s.sequence
        'ACGTAC'
        """
        self.sequence = sequence.upper().strip()   ##sequence their gene and species and added **kwargs for an empty list. 
        self.gene = gene
        self.species = species
        self.kmers = []     ##making an empty list of kmers

    def make_kmers(self, k=3):     ##this to generate a list of k-mers 
        """
        >>> s = Seq("ACGTAC", "gene", "sp")
        >>> s.make_kmers(3)
        >>> s.kmers
        ['ACG', 'CGT', 'GTA', 'TAC']
        """
        self.kmers = [self.sequence[i:i + k] for i in range(len(self.sequence) - k + 1)]

    def fasta(self):                                 
        """
        >>> s = Seq("acgt", "gene1", "human")
        >>> print(s.fasta())
        >human gene1
        ACGT
        """
        return f">{self.species} {self.gene}\n{self.sequence}"            #this is to return the sequence in a fasta format 

    def __str__(self):                                                #this one to return the sequence as a string
        return self.sequence

    def print_record(self):                                            #this one to print out the results. 
        print(f"{self.species} {self.gene}: {self.sequence}")


class DNA(Seq):                                                             #this one here for DNA class
    def __init__(self, sequence, gene, species, geneid, **kwargs):              ## initializes the sequence, gene, species and genid.
        """
        >>> d = DNA("attgggxxx", "geneA", "Human", "G1")
        >>> d.sequence
        'ATTGGGNNN'
        """
        cleaned_seq = re.sub('[^ATGCU]', 'N', sequence.upper().strip())
        super().__init__(cleaned_seq, gene, species)                              #this one to clean the sequence
        self.geneid = geneid

    def reverse_complement(self):                                             #this code to returns the reverse complement of the DNA sequence 
        """   
        >>> d = DNA("ATGC", "g", "s", "id")
        >>> d.reverse_complement()
        'GCAT'
        """
        complement_seq = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'U': 'A', 'N': 'N'}           #these are the complement so each A=T, G=C, with N as unknown base
        return ''.join(complement_seq.get(base, 'N') for base in reversed(self.sequence))

    def six_frames(self):                                                               #this one to generate the siz reading frames for the DNA
        """
        >>> d = DNA("ATGCGT", "g", "s", "id")
        >>> len(d.six_frames())
        6
        """
        frames = [self.sequence[i:] for i in range(3)]
        rev_comp = self.reverse_complement()
        frames += [rev_comp[i:] for i in range(3)]
        return frames

    def analysis(self):                                                          #this one to get the sequence data
        """
        >>> d = DNA("GGCCAATT", "g", "s", "id")
        >>> d.analysis()
        4
        """
        return self.sequence.count('G') + self.sequence.count('C')

    def print_info(self):
        print(f"{self.geneid} {self.species}\n{self.gene}:\n{self.sequence}")


class RNA(DNA):                                                                  #this is RNA class that is going to inhirit what ever done for the the DNA class.
    def __init__(self, sequence, gene, species, geneid, **kwargs):                
        """
        >>> r = RNA("atgctt", "r", "s", "id")
        >>> r.sequence
        'AUGCUU'
        """
        super().__init__(sequence, gene, species, geneid, **kwargs)                   #Because it is RNA, it will inherits the DNA and convert each T to U and make the RNA seq
        self.sequence = self.sequence.replace('T', 'U')
        self.codons = []                                                         #Add self.codons equals to an empty list

    def make_codons(self):
        """
        >>> r = RNA("AUGCGU", "r", "s", "id")
        >>> r.make_codons()
        >>> r.codons
        ['AUG', 'CGU']
        """
        self.codons = [self.sequence[i:i + 3] for i in range(0, len(self.sequence), 3)           # to have the sequence codons as three 
                       if len(self.sequence[i:i + 3]) == 3]

    def translate(self):                                                                 
        """
        >>> r = RNA("AUGGAA", "r", "s", "id")                                     
        >>> r.make_codons()
        >>> r.translate()
        'ME'
        """
        protein = ''
        for codon in self.codons:
            protein += standard_code.get(codon, 'X')                     ### do the translation based on the given dictionary 
        return protein


class Protein(Seq):                                                     ### this is the protein class that inherits fro the seq class
    def __init__(self, sequence, gene, species, geneid, **kwargs):                      
        """
        >>> p = Protein("M@EEK", "prot", "sp", "ID1")                                          
        >>> p.sequence
        'MXEEK'
        """
        clean_seq = re.sub(r'[^A-Za-z]', 'X', sequence.upper().strip())                  ##to clean the seq and replace the non alphabet with X for the unknowen acid
        super().__init__(clean_seq, gene, species)
        self.geneid = geneid

    def total_hydro(self):                                                      ##this one to calculate the total hydrophobicity of the protein to sum the values of the AA from the dictionary 
        """
        >>> p = Protein("MEK", "prot", "sp", "id")
        >>> round(p.total_hydro(), 2)
        -5.5
        """
        return sum(kyte_doolittle.get(aa, 0) for aa in self.sequence)

    def mol_weight(self):                                                      ##this one does the same as the previus one but this time to sum the mol_weight by summing the molecular weighs of the AA from the dictionary
        """
        >>> p = Protein("MEK", "prot", "sp", "id")
        >>> round(p.mol_weight(), 2)
        442.53
        """
        return sum(aa_mol_weights.get(aa, 0) for aa in self.sequence)


if __name__ == "__main__":                                                        ##this is the code for doctes: the doctest code is embded with each code as a green color called docstrings. this is very helpful with knowing if any code is broken or working 
  
    doctest.testmod(verbose=True)
