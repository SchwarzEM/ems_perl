#!/usr/bin/env perl

# prot_feats2tbl.pl -- Erich Schwarz <emsch@its.caltech.edu>, 11/23/2010.
# Purpose: get simple ordered table of protein features for genes.

my $wbgene           = q{};
my $feature          = q{};

my %OK_features = ( SignalP => 1,
                    Signalp => 1,
                    Tmhmm   => 1, 
                    Seg     => 1, 
                    Ncoils  => 1, 

                    # Silently allow these, but don't use them in later output:
                    'N-Glycosylation site' => 1,
                    'Oxidation site'       => 1,
                    'Phosphorylation site' => 1,
);

my $gene2feature_ref = ();

while (my $input = <>) { 
    chomp $input;
    if ( $input =~ / \A \" (WBGene\d+) \" \t \" ([^\t\"]+) \" \s* \z /xms ) { 
        $wbgene  = $1;
        $feature = $2;
        if (! $OK_features{$feature} ) {
            warn "Unauthorized feature in: $input\n";
        }
        # Standardize spelling:
        if ( $feature eq 'Signalp' ) { 
            $feature = SignalP;
        }
        $gene2feature_ref->{$wbgene}->{$feature} = 1;
    }
    elsif ( $input !~ / \A \" (WBGene\d+) \" \t \s* \z/xms ) { 
        warn "Can't parse: $input\n";
    }
}

my @printlist = qw( SignalP Tmhmm Seg Ncoils );
foreach my $gene (sort keys %{ $gene2feature_ref } ) { 
    my @features = ();
    foreach my $feature (@printlist) { 
        if ( $gene2feature_ref->{$gene}->{$feature} ) { 
            push @features, $feature;
        }
    }
    my $text = join "; ", @features;
    print "$gene\t\"$text\"\n";
}

