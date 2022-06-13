#! /usr/bin/env perl
# Prepare ancestral allele input for LDHELMET
# input: tped file for all SNPs; phased data
# Author: Fan Han
use strict;
use warnings;
use Data::Dumper;
use v5.10;

my $data = "outgroup.tped";

open IN, '<', $data or die "No tped is found!\n";
while(<IN>){
	chomp;
	next if (/^#/);
	my @row = split(/\t| /, $_, 5);
	my $chr = $row[0];
	my $pos = $row[3];
	my $geno = $row[4];


	my $A = () = $geno =~ /A/g;
	my $C = () = $geno =~ /C/g;
	my $G = () = $geno =~ /G/g;
	my $T = () = $geno =~ /T/g;
	
	if($A+$C+$G+$T == 0){
		say $chr, "\t", $pos, "\t", "0\t0\t0\t0";
		}else{

	if($A/($A+$C+$G+$T) >= 0.6){
		say $chr, "\t", $pos, "\t", "0.91\t0.03\t0.03\t0.03";
	}elsif($C/($A+$C+$G+$T) >= 0.6){
		say $chr, "\t", $pos, "\t", "0.03\t0.91\t0.03\t0.03";
	}elsif($G/($A+$C+$G+$T) >= 0.6){
		say $chr, "\t", $pos, "\t", "0.03\t0.03\t0.91\t0.03";
	}elsif($T/($A+$C+$G+$T) >= 0.6){
		say $chr, "\t", $pos, "\t", "0.03\t0.03\t0.03\t0.91";
	}else{
		say $chr, "\t", $pos, "\t", "0\t0\t0\t0";
	}
}

}
close IN;
