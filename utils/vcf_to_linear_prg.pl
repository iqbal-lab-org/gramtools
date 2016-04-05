#!/usr/bin/perl -w
use strict;

use Getopt::Long; 


my %vars = ( "vcf" => "zam",
	     "ref" => "",
             "min_freq" =>0,
             "help"=>'');


## For 1000 genomes VCF, can use
## /Net/banyan/data0/users/zam/results/20150429_build_1000g_for_gramtools/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5a.20130502.sites.vcf
## For human ref genome my  /Net/birch/data/zam/ref/hs/fasta/grc37/Homo_sapiens.GRCh37.60.dna.WHOLE_GENOME.fa";
&GetOptions(
    ##mandatory args
    'vcf:s' =>\$vars{"vcf"},
    'ref:s' =>\$vars{"ref"},
    'min_freq:s' => \$vars{"min_freq"},
    'help'  =>\$vars{"help"},
    );


check_args(\%vars);


## load the reference into memory:
my %refseq = ();#chr--> long string
my @chroms=();#collect list of chromosomes
get_ref_seq($vars{"ref"}, \%refseq, \@chroms); 

## parse the VCF and print a linearised PRG in gramtools format
print_linearised_poa(\%refseq, $vars{"ref"}, 
		     $vars{"vcf"}, $vars{"min_freq"},
		     \@chroms);


sub print_linearised_poa_for_one_chr
{
    my ($href_refsequence, $reff, $vcf_file, $chrom, $nextvar, $min_freq)= @_;


    if (!exists $href_refsequence->{$chrom})
    {
	die("Cannot find sequence for chromosome $chrom");
    }

    my $seq = $href_refsequence->{$chrom};
    open(VCF, $vcf_file)||die("Cannot open VCF file $vcf_file");
    my $curr_pos=1; ## 1-based

    while (<VCF>)
    {
	my $lyne  = $_;
	chomp $lyne;

	if ($lyne !~ /^\#/)
	{
	    ## I will work entirely in 1-based coordinates, except at the point of extracting substrings.


	    my @sp = split(/\t/, $lyne);

	    if ($sp[4] !~ /^[ACGTacgt]+$/)
	    {
		##excluding lines which do not properly specify the alternate allele.
		next;
	    }


	    my $info = $sp[7];
	    if ($min_freq>0)
	    {

		if ($info =~ /\;AF=([0123456789\.]+)/)
		{
		    my $freq = $1;

		    if ($freq<$min_freq)
		    {
			next; #ignore this variant if too rare
		    }
		}
		else
		{
		    #if no allele frequency annotation, do not filter by frequency
		    
		}
	    }

	    if ($sp[0] eq $chrom)
	    {

		if ($curr_pos < $sp[1] )
		{
		    my $len = $sp[1]-$curr_pos;
		    print substr($seq, $curr_pos-1, $len);
		    #$curr_pos=$sp[1];
		}

		#replace N with C
		$sp[3]=~ s/[^ACGTacgt]/C/g;



		print $nextvar;#left marker before the site starts
		print $sp[3];		##print the ref allele first
		print $nextvar+1;#even numbers between alleles

		##Now work our way through the alternate alleles
		if ($sp[4]=~ /,/)
		{
		    my @sp2 = split(/,/, $sp[4]);
		    my $i;
		    for ($i=0; $i<scalar(@sp2); $i++)
		    {
			my $allele = $sp2[$i];
			$allele =~ s/[^ACGTacgt]/C/g;
			print $allele;
			if ($i<scalar(@sp2)-1)
			{
			    print $nextvar+1;#even number between alleles
			}
			else
			{
			    print $nextvar;#last one goes back to nextvar (odd)
			}
		    }
		}
		else #we have just one alternate allele
		{
		    $sp[4]=~ s/[^ACGTacgt]/C/g;
		    print $sp[4];
		    print $nextvar;
		}
		$nextvar+=2;
		$curr_pos=$sp[1]+length($sp[3]);
	    }
	    else
	    {
		#ignore
	    }
	}
    }
    close(VCF);

    if ($curr_pos<length($seq)+1)
    {
	print substr($seq, $curr_pos, length($seq)-$curr_pos-1);
    }




    return $nextvar;
}
sub print_linearised_poa
{
    my ($href_refseq, $reference, $vcf, $min_f, $aref_chroms) = @_;

    my @chrs = @$aref_chroms; #(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,"X","Y");

    my $next_var_number = 5;
    foreach my $chr (@chrs)
    {
	## so this is not ideal - parse the VCF once for each chromosome
	$next_var_number = 
	    print_linearised_poa_for_one_chr($href_refseq, $reference, 
					     $vcf, $chr, 
					     $next_var_number, $min_f);
    }

}
sub get_ref_seq
{
    my ($fasta, $href, $aref_chroms) = @_;

    my $chr = "";
    my $seq = "";
    open(FILE, $fasta)||die();
    my $first=1;

    while (<FILE>)
    {
	my $line = $_;
	chomp $line;

	if ($line =~ /^>(\S+)/)
	{
	    if ($first !=1)
	    {
		$seq =~ s/[^ACGTacgt]/C/g;
		$href->{$chr}=$seq;
	    }
	    $first=0;
	    $chr = $1;
	    push @$aref_chroms, $chr;
	    $seq="";
	}
	else
	{
	    $seq .= $line;
	    #$href->{$chr}=($href->{$chr}).$line;
	}
    }
    close(FILE);

    ##now do the final chromosome in the file
    ##replacing N with C
    $seq =~ s/[^ACGTacgt]/C/g;
    $href->{$chr}=$seq;
}



sub check_args
{
    my ($href) = @_;

    if ($href->{"help"})
    {
	print "Usage: perl vcf_to_linear_prg.pl --vcf <VCF> --ref species.fasta --min_freq 0.01\n";
	print "\n";
	print "This script is not super-sophisticated - it builds a lin-PRG\n";
	print "as it sweeps once through the VCF\n";
	print "If it meets a new VCF record that overlaps an old one,\nit will ignore it.\n";
	print "The most important consequence is that it won't encode SNPs \"underneath\" a long deletion\n";
	exit(0);
    }

    if ($href->{"vcf"} eq "")
    {
	die("You must specify a VCF file with --vcf \n");
    }

    if ($href->{"ref"} eq "")
    {
	die("You must specify a reference fasta file with --ref \n");
    }
	
    if (!(-e $href->{"vcf"}))
    {
	print "Specified VCF file ";
	print $href->{"vcf"};
	die(" does not exist");
    }


    if (!(-e $href->{"ref"}))
    {
	print "Specified reference fasta file ";
	print $href->{"ref"};
	die(" does not exist");
    }

    ##let's just avoid any mess with tiny numbers
    if ($href->{"min_freq"}<0.0001)
    {
	$href->{"min_freq"}=0;
    }

}

