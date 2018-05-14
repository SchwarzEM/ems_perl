#!/usr/bin/env perl

# clumsy_blast2wbgenes.pl -- Erich Schwarz <emsch@its.caltech.edu>, 11/26/2008.
# Purpose: get clean line summary of gene hits in TBlastN vs. PS1010 tiny contigs.

my %gene_lines = ();
my $stored     = q{};

while (my $input = <> ) { 
    chomp $input;
    if ( $stored =~ /\S/xms ) { 
        $input =~ s/\A\s*//;
        $input =~ s/\s*\z//;
        $stored .= $input;
        if ( $stored =~ /\A ( .+ WBGene\d+ ( \s+ locus:\S+ )* ) /xms ) {
            $stored = $1;
        }
        $stored =~ s/\sCE\d+\s/ /;
        $stored =~ s/\s+/\t/g;
        $gene_lines{$stored} = 1;
        $stored = q{};
    }
    if ( ( $stored !~ /\S/xms ) and ( $input =~ / \A \s* Query= \s* (\S+) ( .+ WBGene\d+ .* ) \z /xms ) ) {     
        my $cds = $1;
        my $line = $2;
        $cds =~ s/[a-z]\z//;
        $line =~ s/\A\s*//;
        $line =~ s/\s*\z//;
        $stored = $cds . q{ } . $line;
    }
}

foreach my $line ( sort keys %gene_lines ) { 
    print "$line\n";
}


