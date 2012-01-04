#!/usr/bin/perl
#####################################
# Program: clean.pl  -  Date: Thu Sep 15 22:53:28 BRT 2011
# Autor: Fabio C. P. Navarro - Ludwig
# Goal:
#
# Input:
#
# Output:
#
#####################################


use Getopt::Long;
use strict;
use warnings;

my ( $file, $file2, $usage ) = ( "", "", 0 );
GetOptions( "f|file=s"      => \$file,
#           "f2|file2=s"      => \$file2,
           "h|help|usage"  => \$usage )
    || die "Error while parsing comand line arguments";

usage() if( $usage );

my %reads = ();

if( $file eq "" ) {
    usage();
}
#if( $file2 eq "" ) {
#    usage();
#}



#1749_1812_1463_F3_16    r       chr10   91805   86      50M     *       0       0       GGACAAGTTCTACAAGAGCCAGGGAGAACACTGTCCGTGCATGCGTCCCA      @G@=:IIIIIIDH7EIIIDIIBI
#2306_1099_1892_F3_16            chr10   92584   21      50M     *       0       0       ATAAGCCAAGGCTTGCCTTCAGCAACAGGTTTTCCACCACGTCACTGCCC      IIIIIIIIIIIIIB>CCIIIIIII
open ( IN0, "<$file" );
while ( <IN0> ) {
	chomp $_;
	$_ =~ s/r|BC_//;
	my @tokens = split ( /[ \t\cI]+/,$_ );
	my @tmp = split ( /[\_]/,$tokens[0] );
	$reads{"$tmp[0]_$tmp[1]_$tmp[2]"} .= "$tokens[0] $tokens[1] $tokens[2] $tokens[3] $tokens[4] $tokens[8]\n";
}
#open ( IN1, "<$file2" );
#while ( <IN1> ) {
#	chomp $_;
#	my @tokens = split ( /[ \t\cI]+/,$_ );
#}

foreach my $val ( values %reads ) {
	my @tmp = split ( /[\n]/, $val);
	if ( scalar (@tmp) >= 2 ) {
		print "$val";
	}
}

sub usage {
    die "Usage: $0 [options]
Available Options :
   -f  | --file         : Filename in.
   -h  | --help --usage : Print this message.
";

}
