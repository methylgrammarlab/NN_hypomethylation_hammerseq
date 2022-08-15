#!/usr/bin/env perl

use strict;
use List::Util qw/sum/;
my $WHITELIST = "/Users/bermanb/genomic-data-misc/blacklists/hg19-excluding-blacklist.v2.bed.gz";

my @fns = @ARGV;
$::USAGE = "eventsToBedgraph.pl f1.events.txt f2.events.txt ..";

die "$::USAGE\n" unless (@fns);

foreach my $fn (@fns)
{
    die "$fn is not correct filename format\n" unless ($fn =~ /^(.*)\.txt$/);
    my $base = $1;
    my $bg_fn = "${base}.hg19.bedgraph";

    die "Can't write output file $bg_fn\n" unless (open(OUTF,">$bg_fn"));
    die "Can't open file $fn\n" unless (open(F,$fn));

    my $linenum = 0;
    LINE: while (my $line = <F>)
    {
	$linenum++;
	next LINE if ($line =~ /^\#chrom/); # Header

	my $outline = processLine($line);
	print OUTF $outline."\n" if ($outline);
    }
    close(F);
    close(OUTF);

    my $sorted_fn = "${base}.hg19.sorted.bedgraph";

    runCmd("bedtools sort -i $bg_fn | bedtools intersect -a stdin -b ${WHITELIST} > $sorted_fn");
    runCmd("rm -f $bg_fn");
}

sub processLine
{
    my ($inline) = @_;

    chomp $inline;
    my @f = split("\t",$inline);
    my $nf = scalar(@f);
    print STDERR "Wrong number of columns ($nf): $inline\n" unless ($nf == 11);

    # MM, MU, UU, UM
    my @counts = (List::Util::sum(@f[3,7]), List::Util::sum(@f[4,8]), List::Util::sum(@f[5,9]), List::Util::sum(@f[6,10]));

    my $mm = $counts[0];
    my $mu = $counts[1];
    my $total = List::Util::sum($mm,$mu);

    return 0 if ($total<=0);
    
    my @outf = (@f[0..2], sprintf("%s","${mm}/${total}")); #, $total);

    my $outline = join("\t",@outf);
    return $outline;
}

sub runCmd
{
    my ($cmd) = @_;

    print STDERR "$cmd\n";
    print STDERR `$cmd`;
}



