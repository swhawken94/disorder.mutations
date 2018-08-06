#!/usr/bin/perl -w                                                                                                                                                      \
                                                                                                                                                                         

    use strict;
use warnings;

my $genetable=$ARGV[0];  # gene_table.csv                                                                                                                                
my $outPrefix=$ARGV[1];
my $geneName;
my $line;
my @lineArr;
my @domains_list;
my $domain_site;
my $start;
my $end;
my $protein_ID;
my @site;
my $pfam;
my $sequence;

open(GENETABLE_CLEAN, ">$outPrefix.genetable_clean.bed");

open(GENETABLE, '<', $genetable);

while(<GENETABLE>) {
    chomp($_);
    $line=$_;
    @lineArr = split("\t", $line);
    $geneName = $lineArr[0];
    $protein_ID = $lineArr[4];
    $pfam = $lineArr[21];
    $sequence = $lineArr[6];
    @domains_list = split("; ", $pfam);
    $domain_site = $domains_list[-1];
    @site = split("-", $domain_site);
    $start = $site[0];
    $end = $site[1];
    if($start ne "May" || "Jan" || "Aug"){
        print GENETABLE_CLEAN $protein_ID . "\t" . $start . "\t" . $end . "\t" . $geneName . "\t" . $sequence . "\n";
    }
}
