#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;

my %opts;
GetOptions(
    "i|input=s"     => \$opts{input},
    "p|pool-out=s"  => \$opts{pool_out},
    "o|pairs-out=s" => \$opts{pairs_out},
) or die "Invalid arguments\n";

if (!defined $opts{input} || !defined $opts{pool_out} || !defined $opts{pairs_out}) {
    die "Usage: perl $0 --input <speci_genes.txt> --pool-out <gene_pool.txt> --pairs-out <drug_gene_pairs.txt>\n";
}

open(my $in_fh, "<", $opts{input}) or die "Cannot open input $opts{input}: $!\n";
open(my $pool_fh, ">", $opts{pool_out}) or die "Cannot write $opts{pool_out}: $!\n";
open(my $pair_fh, ">", $opts{pairs_out}) or die "Cannot write $opts{pairs_out}: $!\n";

my %gene_pool;
while (my $line = <$in_fh>) {
    chomp $line;
    next if $line =~ /^\s*$/;
    my ($drug, $genes_str) = split(/\t/, $line, 2);
    next unless defined $drug && defined $genes_str;
    $drug =~ s/^\s+|\s+$//g;
    my @genes = split(/,\s*/, $genes_str);
    for my $g (@genes) {
        $g =~ s/^\s+|\s+$//g;
        next unless length $g;
        $gene_pool{$g} = 1;
        print $pair_fh "$drug\t$g\n";
    }
}

for my $gene (sort keys %gene_pool) {
    print $pool_fh "$gene\n";
}

close($in_fh);
close($pool_fh);
close($pair_fh);
