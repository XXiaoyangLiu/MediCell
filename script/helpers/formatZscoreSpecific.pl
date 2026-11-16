#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;

my %opts;
GetOptions(
    "p|pool=s"     => \$opts{pool},
    "w|twas=s"     => \$opts{twas},
    "o|out=s"      => \$opts{out},
    "gene-col=s"   => \$opts{gene_col},
    "ct-col=s"     => \$opts{ct_col},
    "z-col=s"      => \$opts{z_col},
) or die "Invalid arguments\n";

if (!defined $opts{pool} || !defined $opts{twas} || !defined $opts{out}) {
    die "Usage: perl $0 --pool <gene_pool.txt> --twas <twas_table> --out <zscore_out.txt> [--gene-col COL --ct-col COL --z-col COL]\n";
}

open(my $pool_fh, "<", $opts{pool}) or die "Cannot open gene pool $opts{pool}: $!\n";
my %keep_genes;
while (my $line = <$pool_fh>) {
    chomp $line;
    $line =~ s/^\s+|\s+$//g;
    next unless length $line;
    my ($g) = split(/\s+/, $line);
    $keep_genes{$g} = 1 if defined $g && length $g;
}
close($pool_fh);

sub parse_csv_line {
    my ($line) = @_;
    my @fields;
    my $cur = "";
    my $in_quotes = 0;
    for my $i (0 .. length($line) - 1) {
        my $ch = substr($line, $i, 1);
        if ($ch eq q{"}) {
            $in_quotes = !$in_quotes;
        } elsif ($ch eq "," && !$in_quotes) {
            push @fields, $cur;
            $cur = "";
        } else {
            $cur .= $ch;
        }
    }
    push @fields, $cur;
    for my $f (@fields) {
        $f =~ s/^"|"$//g;
        $f =~ s/^\s+|\s+$//g;
    }
    return @fields;
}

open(my $twas_fh, "<", $opts{twas}) or die "Cannot open TWAS file $opts{twas}: $!\n";
my $header = <$twas_fh>;
die "Empty TWAS file: $opts{twas}\n" unless defined $header;
$header =~ s/[\r\n]+$//;

my $is_csv = ($header =~ /,/ && $header !~ /\t/) ? 1 : 0;
my @cols = $is_csv ? parse_csv_line($header) : ($header =~ /\t/ ? split(/\t/, $header) : split(/\s+/, $header));

my ($gene_idx, $ct_idx, $z_idx) = (-1, -1, -1);
my $req_gene = defined $opts{gene_col} && length($opts{gene_col}) ? lc($opts{gene_col}) : "";
my $req_ct   = defined $opts{ct_col}   && length($opts{ct_col})   ? lc($opts{ct_col})   : "";
my $req_z    = defined $opts{z_col}    && length($opts{z_col})    ? lc($opts{z_col})    : "";

for my $idx (0 .. $#cols) {
    my $name = lc($cols[$idx]);
    $name =~ s/^"|"$//g;
    $name =~ s/^\s+|\s+$//g;
    if ($req_gene ne "") {
        $gene_idx = $idx if $gene_idx < 0 && $name eq $req_gene;
    } elsif ($gene_idx < 0 && ($name eq "gene" || $name eq "gene_name" || $name eq "symbol")) {
        $gene_idx = $idx;
    }
    if ($req_ct ne "") {
        $ct_idx = $idx if $ct_idx < 0 && $name eq $req_ct;
    } elsif ($ct_idx < 0 && ($name eq "celltype" || $name eq "cell_type")) {
        $ct_idx = $idx;
    }
    if ($req_z ne "") {
        $z_idx = $idx if $z_idx < 0 && $name eq $req_z;
    } elsif ($z_idx < 0 && ($name eq "zscore" || $name eq "z")) {
        $z_idx = $idx;
    }
}

if (($req_gene ne "" && $gene_idx < 0) || ($req_ct ne "" && $ct_idx < 0) || ($req_z ne "" && $z_idx < 0)) {
    die "Specified column name not found in TWAS header: $opts{twas}\n";
}

# Fallback to legacy 13+ column whitespace format if header names were not matched.
if ($gene_idx < 0 || $ct_idx < 0 || $z_idx < 0) {
    $gene_idx = 1;
    $z_idx    = 2;
    $ct_idx   = 12;
}

my %best_z;
while (my $line = <$twas_fh>) {
    $line =~ s/[\r\n]+$//;
    next if $line =~ /^\s*$/;
    my @a = $is_csv ? parse_csv_line($line) : ($line =~ /\t/ ? split(/\t/, $line) : split(/\s+/, $line));
    next if @a <= $gene_idx || @a <= $ct_idx || @a <= $z_idx;
    my $gene = $a[$gene_idx];
    my $ct   = $a[$ct_idx];
    my $z    = $a[$z_idx];
    next unless defined $gene && defined $ct && defined $z;
    next unless exists $keep_genes{$gene};
    next unless $z =~ /^[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?$/;
    my $key = "$gene\t$ct";
    if (!exists $best_z{$key} || abs($z) > abs($best_z{$key})) {
        $best_z{$key} = $z;
    }
}
close($twas_fh);

open(my $out_fh, ">", $opts{out}) or die "Cannot write $opts{out}: $!\n";
print $out_fh "Gene\tCellType\tZscore\n";
for my $key (sort keys %best_z) {
    print $out_fh "$key\t$best_z{$key}\n";
}
close($out_fh);
