#!/usr/bin/python
import argparse
import os
import sys
import subprocess
try:
    import libbowtie
    libbowtie_loaded=True
except ImportError:
    libbowtie_loaded=False
try:
    import pysam
    pysam_loaded=True
except ImportError:
    pysam_loaded=False
    
from Bio import SeqIO

def parseFilename(filename):
    """
    Returns the base filename and extension as a tuple.
    """
    splitfilename=filename.split('.')
    base="".join(splitfilename[:-1])
    extension=splitfilename[-1]
    return base,extension

def callBowtie2(reference_genome,seq_files,output_filename,input_type_arg,bowtie2_arguments=[]):
    """
    Uses libbowtie to call bowtie2 if available, otherwise calls it from the command line (assumes on path)
    """
    if bowtie2_arguments==[]:
        bowtie2_arguments=['-t','--threads 4','--local','--very-sensitive-local']
    if libbowtie_loaded:
        print
        print 'Calling libbowtie'
        print
        libbowtie.bowtie(reference_genome,output_filename,unpaired=",".join(seq_files),extra_arguments=bowtie2_arguments+[input_type_arg])
    else:
        bw_commandline='bowtie2 {0} {1} -x {2} -U {3} -S {4}.sam'.format(" ".join(bowtie2_arguments),input_type_arg,reference_genome,",".join(seq_files),output_filename).split(" ")
        #make the call to bowtie2 (assumes the executable is on the path and the index is reachable)   
        print
        print 'Calling bowtie2 with command line:'
        print " ".join(bw_commandline)
        print
        try:
            print subprocess.check_output(bw_commandline,stderr=subprocess.STDOUT)
            return True
        except Exception as ex:
            print "Error({0}): {1}".format(ex.errno, ex.strerror)
            return False

def makeFasta(input_filename, output_filename=None, sequence_name=None):
    """
    Converts raw sequence files into FASTA files by prepending a sequence name (defaults to filename) on the first line.
    """
    base_filename= parseFilename(input_filename)[0]
    if sequence_name==None:
        sequence_name = base_filename
    if output_filename==None:
        output_filename=base_filename+".fasta"
    
    try:
        with open(input_filename,'rb') as raw_file:
            file_contents=raw_file.read()
            raw_file.close()
        with open(output_filename,'w') as new_file:
            new_file.write('>'+sequence_name+'\n')
            new_file.write(file_contents)
            new_file.close()
        return 1
    except Exception as ex:
        print ''.join(ex.args())
        return 

def convertABItoFASTQ(input_filename,output_filename):
    """
    Converts an ABI sequence file containing called bases and quality scores into a FASTQ file.
    """
    with open(input_filename,'rb') as infile:
        sequences = SeqIO.parse(infile, "abi")
        with open(output_filename,'w') as outfile:
            count=SeqIO.write(sequences,outfile,'fastq')
    return count

def SamtoBam(sam_filename, bam_filename=None):
    """
    Converts a sam file to a bam file using pysam.
    """
    base_filename=parseFilename(sam_filename)[0]
    if bam_filename==None:
        bam_filename=base_filename+".bam"
    try:
        infile = pysam.Samfile( sam_filename, "r" )
        outfile = pysam.Samfile( bam_filename, "wb", template = infile ) 
        for s in infile: outfile.write(s)
    except Exception as ex:
        return "Error converting {0} to {1}: {2}".format(sam_filename,bam_filename," ".join(ex.args()))
        raise ex

def convertSortAlign(output_filename):
        #Pregenerate file names for all the intermediate steps (output_filename is the output of the Bowtie2 alignment)
    #Note that the file extension is not always given depending on the input conventions of the tool
    sam_filename=output_filename+'.sam'
    bam_filename=output_filename+'.bam'
    sorted_filename_input=output_filename+'_sorted'
    sorted_filename_output=output_filename+'_sorted.bam'
    
    #convert sam to bam
    print 'Converting {0} to {1} . . .'.format(sam_filename,bam_filename)
    try:
        SamtoBam(sam_filename,bam_filename)
    except Exception as ex:
        print "Error converting sam to bam ({0}): {1}".format(ex.errno, ex.strerror)
        return False   
    
    #sort
    print 'Sorting {0} -> {1}'.format(bam_filename,sorted_filename_output)
    try:
        pysam.sort(bam_filename,sorted_filename_input)
    except Exception as ex:
        print "Error sorting bam file ({0}): {1}".format(ex.errno, ex.strerror)
        return False   
    
    #index
    print 'Indexing {0} . . .'.format(sorted_filename_output)
    try:
        pysam.index(sorted_filename_output)
    except Exception as ex:
        print "Error indexing bam file ({0}): {1}".format(ex.errno, ex.strerror)
        return False   
    
    print
    print 'Done'
    return True

def alignSangerFromRaw(sequence_directory,reference_genome,output_filename,bowtie2_args=None):
    """
    Generates a sorted, indexed alignment against reference_genome from all the Sanger 
        sequencing reads (in ABI format with quality scores) in sequence_directory.
        Sequence files must be have .AB1 extensions.
        
        First it uses BioPython to generate a FASTQ file for each .AB1 file in sequence_directory.
        Next it calls Bowtie2 to generate an alignment against the given reference_genome 
            (assumes the executable is on the path and the index is reachable)
            bowtie2: any additional arguments to pass to bowtie2. If none given, defaults to:
            -t --threads 4 --local --very-sensitive-local
                
        Lastly, it uses Pysam to convert the .sam output of bowtie2 to .bam, sort and index the .bam file
        
        Output files will have prefix output_filename and be generated in sequence_directory.
    """
    input_extension='seq'
    output_extension='fasta'
    seq_files=[]
        
    try:
        os.chdir(sequence_directory)
        dir_listing=os.listdir(sequence_directory)
    except Exception as ex:
        print "Error changing to {0} ({1}): {2}".format(sequence_directory,ex.errno, ex.strerror)
        return False   

    try:
        for filename in dir_listing:
            base_filename,extension=parseFilename(filename)
            if extension==input_extension:
                if makeFasta(base_filename+"."+input_extension,base_filename+"."+output_extension):
                    seq_files.append(base_filename+"."+output_extension)
                    print 'Converted {0} to {1}'.format(base_filename+"."+input_extension,base_filename+"."+output_extension)
    except Exception as ex:
        print "Error converting sequence files ({0}): {1}".format(ex.errno, ex.strerror)
        return False   
    
    if callBowtie2(reference_genome,seq_files,output_filename,'-f',bowtie2_args):
        return convertSortAlign(output_filename)
    else:
        return False

def alignSangerFromABI(sequence_directory,reference_genome,output_filename,bowtie2_args=None):
    """
    Generates a sorted, indexed alignment against reference_genome from all the Sanger 
        sequencing reads (in ABI format with quality scores) in sequence_directory.
        Sequence files must be have .AB1 extensions.
        
        First it uses BioPython to generate a FASTQ file for each .AB1 file in sequence_directory.
        Next it calls Bowtie2 to generate an alignment against the given reference_genome 
            (assumes the executable is on the path and the index is reachable)
            bowtie2: any additional arguments to pass to bowtie2. If none given, defaults to:
            -t --threads 4 --local --very-sensitive-local
                
        Lastly, it uses Pysam to convert the .sam output of bowtie2 to .bam, sort and index the .bam file
        
        Output files will have prefix output_filename and be generated in sequence_directory.
    """
    input_extension='ab1'
    output_extension='fastq'
    seq_files=[]

    try:
        os.chdir(sequence_directory)
        dir_listing=os.listdir(sequence_directory)
    except Exception as ex:
        print "Error changing to {0} ({1}): {2}".format(sequence_directory,ex.errno, ex.strerror)
        return False   
    try:
        for filename in dir_listing:
            base_filename,extension=parseFilename(filename)
            if extension==input_extension:
                if convertABItoFASTQ(base_filename+"."+input_extension,base_filename+"."+output_extension):
                    print 'Converted {0} to {1}'.format(base_filename+"."+input_extension,base_filename+"."+output_extension)
                    seq_files.append(base_filename+"."+output_extension)
    except Exception as ex:
        print "Error converting sequence files ({0}): {1}".format(ex.errno, ex.strerror)
        return False   

    if callBowtie2(reference_genome,seq_files,output_filename,'-q',bowtie2_args):
        return convertSortAlign(output_filename)
    else:
        return False
   
def insA_test():
    try:
        alignSangerFromABI('/home/phageghost/sequencing/Sanger/insA_uspC10seq','MG1655','insA_uspC10_Sanger_bowtie2')
    except subprocess.CalledProcessError as ex:
        print ex.returncode
        print ex.output

def ABItest():
    os.chdir('/home/phageghost/sequencing/Sanger/insA_uspC10seq')
    print convertABItoFASTQ('GK-1F-F_B09.ab1','GK-1F-F_B09.fsq')

def main():
    arg_parser=argparse.ArgumentParser(description='Generates a sequence alignment from Sanger sequencing reads. Currently supports raw sequence data and ABI format.')
    
    arg_parser.add_argument('sequence_directory', help='Directory containing the Sanger sequencing reads')
    arg_parser.add_argument('bowtie2_index', help='Name of the genome index to be used by Bowtie2 (must be accessible on the bowtie2 path)')
    arg_parser.add_argument('output_filename', help='Base filename for the aligned output')
    arg_parser.add_argument('-r', '--raw', action='store_true', help='specify that input sequences are raw reads')
    arg_parser.add_argument('-a','--abi', action='store_true', help='specify that input sequences are abi chromatograms')
    arg_parser.add_argument('-b','--bowtie2args',help='a string containing any additional arguments to pass to bowtie2')

    args=arg_parser.parse_args()
    
    if args.raw and args.abi:
        #these are mutually exclusive!
        print "Either --raw or --abi may be chosen, but not both."
    if pysam_loaded:
        print "pysam found"
    if libbowtie_loaded:
        print "libbowtie found"
    
    print 'args.bowtie2args={}'.format(args.bowtie2args)
    
    if args.raw:
        alignSangerFromRaw(args.sequence_directory, args.bowtie2_index, args.output_filename, args.bowtie2args)
    else:
        alignSangerFromABI(args.sequence_directory, args.bowtie2_index, args.output_filename, args.bowtie2args)
if __name__ == "__main__":
    sys.exit(main())
