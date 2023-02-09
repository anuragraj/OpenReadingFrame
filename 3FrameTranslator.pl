#! /usr/bin/perl

# --------------------------------------------------------------------------------------
#Program to translate DNA sequence fasta file to three frame translations and writing all the ORFs 
#Author: Anurag Raj [anurag.igib@gmail.com]
#Date: 15/06/2022
#usgae: perl 3FrameTranslator.pl <input.fasta> <output.fasta>
# --------------------------------------------------------------------------------------

use strict;
use warnings;

use Bio::SeqIO;
use Bio::Seq;
use Bio::Tools::CodonTable;

my $SCRIPT_NAME = $0;
my $DESCRIPTION = "Translate sequences in three frame seqeunces from fasta file.";
my $AUTHOR = "Anurag Raj [anurag.igib\@gmail.com]";

#input files
my $infile = $ARGV[0];
my $outfile = $ARGV[1];

# Usage information
die "Usage: $0 <fasta-file> <output-name>\n", if (@ARGV != 2);

#Write the output file
open OUT , ">$outfile" or die "Cannot open $outfile for writing results:$!\n";

#Read the input file
my $seqio = Bio::SeqIO -> new(-file => "$infile", -format => 'Fasta');

my $x = 0; #counter for number of sequences

#translate the sequences
print "Translating the DNA transcripts to three frame...\n";
print "Translated sequences: ";

while(my $seq = $seqio -> next_seq){
	my $header = $seq -> id;

	#skip mitochondrial sequences
	if($header =~ /Mt_tRNA/ or $header =~ /Mt_rRNA/){
		next;
	}
	
	#translating to standard codon table
	my $prot_obj1 = $seq-> translate(-frame => 0, -terminator => '', -codontable_id => 1);
	my $prot_seq1 = $prot_obj1 -> seq;
	ORFS($prot_seq1,$header,1);
	
	my $prot_obj2 = $seq-> translate(-frame => 1, -terminator => '', -codontable_id => 1);
	my $prot_seq2 = $prot_obj2 -> seq;
	ORFS($prot_seq2,$header,2);
	
	my $prot_obj3 = $seq-> translate(-frame => 2, -terminator => '', -codontable_id => 1);
	my $prot_seq3 = $prot_obj3 -> seq;
	ORFS($prot_seq3,$header,3);
	$x++;

	print "$x \r" if ($x % 1000 == 0);
}
print "\nTotal $x sequences translated.\n";

close OUT;

##################################################
# Subroutine to write the ORFs to the output file
##################################################
sub ORFS{
    my ($string,$header,$frame) = @_;
    my @data = split(/\*/,$string);
    my $i = 1;
    foreach my $seq(@data){

		#write ORFs only when length more than 50aa
        if(length($seq)>=50){
        	my @header_data = split(/\|/,$header);
        	my $header_new = $header_data[0].".Frame".$frame.".".$i."|".join( "|",@header_data[1..$#header_data]); 
            print OUT ">$header_new\n$seq\n";
            $i++;
        }
    }
}