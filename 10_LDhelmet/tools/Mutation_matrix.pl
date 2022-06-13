#! /usr/bin/env perl
# Prepare mutation matrixinput for LDHELMET
# input: vcf file for target population. It is required to be phased.
# Author: Fan Han
use strict;
use warnings;
use Data::Dumper;
use v5.10;
use List::Util qw(max);

if(!@ARGV){
	die "Mutation_matrix.pl population.vcf > population_mutation.matrix\n";
}

my $ref = "outgroup_ancestral.allele";
my $data = shift @ARGV;

# read ancestral allele
my %ancestor;
open IN, '<', $ref or die "No ancestral allele is found!\n";
while (<IN>){
	chomp;
	my @row = split(/\t/);
	if($row[2] == 0.91 ){
		$ancestor{$row[0]}{$row[1]} = "A";
	}
	if($row[3] == 0.91 ){
		$ancestor{$row[0]}{$row[1]} = "C";
	}
	if($row[4] == 0.91 ){
		$ancestor{$row[0]}{$row[1]} = "G";
	}
	if($row[5] == 0.91 ){
		$ancestor{$row[0]}{$row[1]} = "T";
	}
}
close IN;


my %matrix;
my %Fi;
open IN, '<', $data or die "No vcf is found!\n";
while(<IN>){
	chomp;
	next if (/^#/);
	my @row = split(/\t/, $_, 10);
	if(exists $ancestor{$row[0]}{$row[1]}){
		$Fi{$ancestor{$row[0]}{$row[1]}}++;

		if($row[3] eq $ancestor{$row[0]}{$row[1]}){
			$matrix{$row[3]}{$row[4]}++;

		}elsif($row[4] eq $ancestor{$row[0]}{$row[1]}){
			$matrix{$row[4]}{$row[3]}++;
		}
	}
	


}
close IN;


foreach my $k1 (sort keys %matrix){
	print $k1;
	for my $k2 (sort keys %{$matrix{$k1}}){
		
		print "\t", "$k2:", $matrix{$k1}{$k2};
	}
	print "\t$Fi{$k1}\n";
}
