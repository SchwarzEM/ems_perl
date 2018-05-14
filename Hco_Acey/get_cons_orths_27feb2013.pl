#!/usr/bin/env perl

use strict;
use warnings;

my $header        = "Gene\tStrict_orthologs [consensus]\n";
my @infiles       = @ARGV;
my $threshold     = @infiles;
my $gene2orth_ref = ();

foreach my $infile (@infiles) {
    open my $INFILE, '<', $infile or die "Can't open orthology table $infile: $!";
    while (my $input = <$INFILE>) {
        chomp $input;
        if ( $input =~ /\A (\S+) \t (\S+) \b/xms ) { 
            my $gene = $1;
            my $orth = $2;
            if ( $gene ne 'Gene' ) { 
                # Set things up so I can count whether an ortholog was seen in X files:
                $gene2orth_ref->{$gene}->{'orth'}->{$orth}->{'file'}->{$infile} = 1;
            }
        }
    }
    close $INFILE or die "Can't close filehandle to orthology table $infile: $!";
}

my @observed_genes = sort keys %{ $gene2orth_ref };
my @eligible_genes = ();

foreach my $obs_gene (@observed_genes) {
    if ( exists $gene2orth_ref->{$obs_gene}->{'orth'} ) { 
        push @eligible_genes, $obs_gene;
    }
}

foreach my $elig_gene (@eligible_genes) { 
    my @elig_orths = sort keys %{ $gene2orth_ref->{$elig_gene}->{'orth'} };
    foreach my $elig_orth (@elig_orths) { 
        my @obs_files = sort keys %{ $gene2orth_ref->{$elig_gene}->{'orth'}->{$elig_orth}->{'file'} };
        my $file_count = @obs_files;
        if ( $file_count == $threshold ) {
            print $header if $header;
            $header = q{};
            print "$elig_gene\t$elig_orth\n";
        }
    }
}

