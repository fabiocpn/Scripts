#!/usr/bin/perl 

use Getopt::Long;
use strict;
use Switch;

my $file = "";
my $orig = "";
my $usage = "";
my $ref = "";
my $limite = 0;
my $chrref = "";

my $nlinha = 4;

GetOptions("f|file=s"      => \$file,
           "l|limit=i"     => \$limite,
		   "r|ref=s"       => \$ref,
		   "o|orig=s"      => \$orig,
           "h|help|usage"  => \$usage)
    || die "Error while parsing comand line arguments";

usage() if($usage);

if($file eq "" || $limite == 0 || $orig eq "" || $ref eq "") {
    usage();
}

if(defined($ARGV[0])) {
    usage();
}

    open(ORIG,"<$orig");
	my ( $first, $index ) = ( 1, 0 );
	my ( $read, $seq ) = ( "", "" );
	my %horig = ();
	
	while ( <ORIG> ){
		chomp $_;
		my @tokens = split (/[ \t\cI]+/,$_);
		if ( $index == 0 ) {
		#if ( $_ =~ /^\>/ ) {
			$read = $tokens[0];
			$index = 1;
		}
		else {
			$seq = $_;
			$horig{$read} = $seq;
			$index = 0;
		}
	}
	

	

&main($file, $limite); 

sub main {
    my ($filein, $limite) = @_;

    my $begin = 0;
    my $nexon = 0;
    my @arraylinha = ();
    my @arraystrand = ();
    my $strand = "";
	my $lastchr = "";
	my $readseq = "";
	my $lastreadseq = "";

	open(DEBUG,">DEBUG");
    open(IN,"<$filein");
	my ( $last_read, $read ) = ( "", "" );
    while(<IN>) {
	print DEBUG "$_";
	chomp $_;

	my $line = "$_\n";
	#Verifica se eh comeco de sequencia
	if($_ =~ /^>|chr/) {
		
		#my $origline = "";
		if ( $_ =~ /[>]/ ) {
			my @tokens = split (/[ \t\cI]+/,$_);
			$last_read = $read;
			$read = $tokens[0];
		#	while ( $origline ne $read ) {
		#		$origline = <ORIG>;
		#		if ( ! $origline ) {
		#			seek( ORIG, 0, 0 );
		#			$origline = <ORIG>;
		#		}
		#		chomp($origline);
		#	}
		#	$lastreadseq = $readseq;
		#	$readseq = <ORIG>;
		#	chomp($readseq);
		}	

		if ( $_ =~ /chr/ && $_ ne $lastchr ) {
				$lastchr = $_;
				&loadchr($ref, $_, \$chrref);
		}
	    if(scalar(@arraylinha) && $chrref ) {
			&processa($strand, \@arraylinha, $limite, \$chrref, $horig{$last_read});
	    }
		else {
	    	for(my $y=0;$y<scalar(@arraylinha);$y++) {
				$arraylinha[$y][0] =~ s/^\t//g;
				print "\t".$arraylinha[$y][0]."\n";
				print "\t".$arraylinha[$y][1]."\n";
				print "\t".$arraylinha[$y][2]."\n";
				print "\t".$arraylinha[$y][3]."\n";
			}
		}
		
	   	print "$line";
	    $begin = 0;
	    $nexon = 0;
	    @arraylinha = ();
	    @arraystrand = ();
	    $strand = "";
	} else {
	    ###Se for a primeira linha, faca
	    if($_ ne "") {
		if($begin == 0 && $_ ne "") {
		    @arraystrand = split(/[ \t\cI]+/, $_);
		    $strand = $arraystrand[3];
		    $arraylinha[$nexon][$begin] = $_;
		} else {
		    $arraylinha[$nexon][$begin] = $_;
		}
		$begin++;
		###Se for a ultima linha do exon faca
		if($begin == 4 ) {
		    $begin = 0;
		    $nexon++;
		}
	    }
	}
    }
    if(scalar(@arraylinha)) {
	&processa($strand, \@arraylinha, $limite, \$chrref, $horig{$read});
    close (IN);
    }
}

sub processa {
	my ( $strand, $arraylinha, $limite, $ref, $readseq )  = @_;
    my $strand = shift;
    my $arraylinha = shift;
    my $limite = shift;
	my $ref = shift;
	my $chrref = ${$ref};
	my $readseq = shift;

    my @current = ();
    my @next = ();
    my @arrayret = ();
    my @pos = ();

    my $seq1 = "";
    my $seq2 = "";
    my $seq3 = "";
    my $cab = "";

    my $distgbdna = 0;
    my $distgenom = 0;

    ###Verifica se a seq tem 1 exon, se tiver mais do que 1, concatena, caso contrario, somente imprime###
    if(scalar(@{$arraylinha}) > 1) {
	###Aqui roda ate o penultimo exon###
	for(my $i = 0;$i<(scalar(@{$arraylinha})-1);$i++) {
	    ###Pegando os dados do atual e o proximo exon###
	    @current = split(/[ \t\cI]+/, $arraylinha->[$i][0]);
	    @next = split(/[ \t\cI]/, $arraylinha->[$i+1][0]);
	    
	    ###Tirando os tabels das sequencias###
	    $arraylinha->[$i][1] =~ s/^\t+//g;
	    $arraylinha->[$i][2] =~ s/^\t+//g;
	    $arraylinha->[$i][3] =~ s/^\t+//g;
	    
	    $arraylinha->[$i+1][1] =~ s/^\t+//g;
	    $arraylinha->[$i+1][2] =~ s/^\t+//g;
	    $arraylinha->[$i+1][3] =~ s/^\t+//g;
	    
	    ###Pegando a distancia dos exons. Se for plus eh uma conta, se for minus eh outra###
	    if($strand eq "+") {
			$distgbdna = ($next[1]-$current[2])-1;
			$distgenom = ($next[4]-$current[5])-1;
	    } else {
			$distgbdna = ($current[1]-$next[2])-1;
			$distgenom = ($next[4]-$current[5])-1;
	    }
	    
	    ###Verificando se a distancia eh maior do que o limite passado como argumento. Se for, faca###
	    if( ($distgbdna <= $limite && $distgbdna > 0 && $distgenom == 0) || ($distgenom <= $limite && $distgenom > 0 && $distgbdna == 0) ) {
		###Pegando o cabecalho
		my @basesc = split(/\//, $current[6]);
		my @basesn = split(/\//, $next[6]);
		
		###Pegando os dados para mais tarde montar o cabecalho do exon###
		###Se for plus, pega de um jeito, se for minus, pega de outro###
		if($strand eq "+") {
		    if(!defined($pos[0]) && !defined($pos[2]) && !defined($pos[4]) && !defined($pos[5])) {
			$pos[0] = $current[1];
			$pos[2] = $current[4];
			$pos[4] = $basesc[0];
			$pos[5] = $basesc[1];
		    }
		    $pos[1] = $next[2];
		    $pos[3] = $next[5];
		    
		    $pos[4] += $basesn[0];
		    $pos[5] += $basesn[1];
		} else {
		    if(!defined($pos[0]) && !defined($pos[2]) && !defined($pos[4]) && !defined($pos[5])) {
			$pos[1] = $current[2];
			$pos[2] = $current[4];
			$pos[4] = $basesc[0];
			$pos[5] = $basesc[1];
		    }
		    $pos[0] = $next[1];
		    $pos[3] = $next[5];
		    
		    $pos[4] += $basesn[0];
		    $pos[5] += $basesn[1];
		}
		
		###Definindo as variaveis aonde vao ficar os gaps###
		my $gapgdna = "";
		my $gapcent = "";
		my $gapgeno = "";
		my $gapgdna1 = "";
		my $gapcent1 = "";
		my $gapgeno1 = "";
	
		###Nesta parte concatena as sequencias. Se for o primeiro exon, faca###
		if(defined($seq1) && $seq1 eq "") {
		    $seq1 .= $arraylinha->[$i][1];
		}
		if(defined($seq2) && $seq2 eq "") {
		    $seq2 .= $arraylinha->[$i][2];
		}
		if(defined($seq3) && $seq3 eq "") {
		    $seq3 .= $arraylinha->[$i][3];
		}
		

		my $genoindex = "";		
		my $gdnaindex = "";
	

		if($strand eq "+") {	
			my $tempgeno = $seq3;
			my $tempgdna = $seq1;
			$tempgeno =~ s/-//g;
			$tempgdna =~ s/-//g;

			$genoindex = $pos[2]+length($tempgeno)-1;		
			$gdnaindex = $pos[0]+length($tempgdna)-1;
		}
		else {
            #$genoindex = $pos[3]+length($arraylinha->[$i+1][3])-1;
			$genoindex = $current[5];
            $gdnaindex = $pos[0]+length($arraylinha->[$i+1][1])-1;
		}
		###Nesta parte esta montando os gaps, tanto do gdna, quanto do genoma###
		for(my $x=0;$x<$distgbdna;$x++) {
		    #$gapgdna .= "";
		    $gapcent .= "Y";
		    $gapgeno .= "-";
		}
		for(my $x=0;$x<$distgenom;$x++) {
		    $gapgdna1 .= "-";
		    $gapcent1 .= "X";
		    #$gapgeno1 .= "x";
			#my $temp = $genoindex+$x;
			#print "$temp\n";
		}

		if($strand eq "+") {
			$gapgeno1 = substr(${$ref},$genoindex,$distgenom);
			$gapgdna = lc(substr($readseq,$gdnaindex,$distgbdna));
		} else {
			$gapgeno1 = &invcomp ( substr(${$ref},$genoindex,$distgenom) );
			#$gapgdna = &invcomp ( substr($readseq,$gdnaindex,$distgbdna) );
			$gapgdna = lc(substr($readseq,$gdnaindex,$distgbdna));
			#$gapgdna = substr(&invcomp($readseq),$gdnaindex,$distgbdna);
		}
#		my $debug = substr($readseq,$gdnaindex-5,$distgbdna+10);

		###Concatenando os exons###
		$seq1 .= $gapgdna.$gapgdna1;
		$seq2 .= $gapcent.$gapcent1;
		$seq3 .= $gapgeno.$gapgeno1;
		
		###Concatenando o segundo exon###
		$seq1 .= $arraylinha->[$i+1][1];
		$seq2 .= $arraylinha->[$i+1][2];
		$seq3 .= $arraylinha->[$i+1][3];
	    } else {
		###Neste else, significa que a distancia entre os exons eh maior do que o limite###
		###Este if signifca que eh o penultimo exon####
		if($i == scalar(@{$arraylinha})-2) {
		    ###Se existir sequencias concatenadas, insira em um array####
		    if(defined($seq1) && $seq1 ne "" && defined($seq3) && $seq3 ne "") {
			if(scalar(@pos)) {
			    $cab = "$pos[0] $pos[1] $strand $pos[2] $pos[3] $pos[4]/$pos[5]";
			}
			push(@arrayret, [$cab, $seq1, $seq2, $seq3]);
			push(@arrayret, [$arraylinha->[$i+1][0], $arraylinha->[$i+1][1], $arraylinha->[$i+1][2], $arraylinha->[$i+1][3]]);
		    } else {
			###Senao, insira as sequencias normais###
			push(@arrayret, [$arraylinha->[$i][0], $arraylinha->[$i][1], $arraylinha->[$i][2], $arraylinha->[$i][3]]);
			push(@arrayret, [$arraylinha->[$i+1][0], $arraylinha->[$i+1][1], $arraylinha->[$i+1][2], $arraylinha->[$i+1][3]]);
		    }
		} else {
		    ###Se nao for o penultimo exon, faca####
		    ###Se existir sequencias concatenadas, faca###
		    if(defined($seq1) && $seq1 ne "" && defined($seq3) && $seq3 ne "") {
			if(scalar(@pos)) {
			    $cab = "$pos[0] $pos[1] $strand $pos[2] $pos[3] $pos[4]/$pos[5]";
			}
			push(@arrayret, [$cab, $seq1, $seq2, $seq3]);
		    } else {
			###Senao, insira as sequencias normais###
			push(@arrayret, [$arraylinha->[$i][0], $arraylinha->[$i][1], $arraylinha->[$i][2], $arraylinha->[$i][3]]);
		    }
		}
		
		###Limpando as variaveis###
		$cab = "";
		$seq1 = "";
		$seq2 = "";
		$seq3 = "";
		@pos = ();
	    }
	}
	###Nesta parte, verifica se existe sequencias concatenadas do comeco ao fim da sequencia###
	if(defined($seq1) && $seq1 ne "" && defined($seq3) && $seq3 ne "") {
	    if(scalar(@pos)) {
		$cab = "$pos[0] $pos[1] $strand $pos[2] $pos[3] $pos[4]/$pos[5]";
	    }
	    push(@arrayret, [$cab, $seq1, $seq2, $seq3]);
	}
	
	###Se existir algo no array que foi criado com as concatenacoes, faca###
	if(scalar(@arrayret)) {
	    ###Rode o array, exon a exon imprimindo o resultado###
	    for(my $y=0;$y<scalar(@arrayret);$y++) {
		$arrayret[$y][0] =~ s/^\t//g;
		print "\t".$arrayret[$y][0]."\n";
		print "\t".$arrayret[$y][1]."\n";
		print "\t".$arrayret[$y][2]."\n";
		print "\t".$arrayret[$y][3]."\n";
	    }
	}
    } else {
	###Aqui significa que somente existe 1 exon. Ae eh soh para imprimir###
	print $arraylinha->[0][0]."\n";
	print $arraylinha->[0][1]."\n";
	print $arraylinha->[0][2]."\n";
	print $arraylinha->[0][3]."\n";
    }
}

sub loadchr {
	my ( $ref, $chr, $schrref ) = @_;
	
	$chr =~ s/[\cI\s ]+//g;
	${$schrref} = "";
	open(CHR,"<$ref/$chr.fa");
    while(<CHR>) {
		chomp $_;
		my $seq = lc( $_ );
		if ( $_ !~ /[>]/ ){
			${$schrref} .= $seq;
			#push( @{$chrref}, split(//,$seq) );
		}
	}
}

sub invcomp {
	my $seq = shift;
	$seq = lc ( $seq );
	my $seq_aux = reverse ( $seq );
	$seq_aux =~ s/a/T/g;
	$seq_aux =~ s/t/A/g;
	$seq_aux =~ s/c/G/g;
	$seq_aux =~ s/g/C/g;
	$seq = lc ( $seq_aux );
	return $seq;
}

sub usage {
    die "Usage: $0 [options]
Available Options :
   -f  | --file         : Filename in.
   -l  | --limit        : Limit of the gaps.
   -r  | --ref       	: Reference Genome.
   -o  | --orig        : Original Reads.
   -h  | --help --usage : Print this message.
";
}
