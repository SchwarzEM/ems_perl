#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;

my $seg      = q{};
my $proteome = q{};

my $gene = q{};
my $tx   = q{};

my $data_ref;
my $verbose;
my $help;

GetOptions ( 'seg=s'      => \$seg,
             'proteome=s' => \$proteome,
             'verbose'    => \$verbose,
             'help'       => \$help, );

if ( $help or (! $seg) or (! $proteome) ) { 
    die "Format: summarize_psegs_14nov2012.pl\n",
        "    --seg|-s       [FASTA of low-complexity pseg domains]\n",
        "    --proteome|-p  [proteome; both -s and -p assume AUGUSTUS gene/trx names]\n",
        "    --verbose|-b   [print *all* isoforms values for SEG, rather than default of the highest %]\n",
        "    --help|-h      [print this message]\n",
        ;
}

open my $PROTEOME, '<', $proteome or die "Can't open proteome: $proteome\n";
while (my $input = <$PROTEOME>) { 
    chomp $input;
    if ( $input =~ /\A > /xms ) { 
        if ( $input =~ /\A > ((\S+\.g\d+)\.t\d+) /xms ) { 
            $tx   = $1;
            $gene = $2;
        }
        else { 
            die "Can't parse: $input\n";
        }
    }
    else {
        $input =~ s/\s//g;
        my $len = length($input);
        $data_ref->{'gene'}->{$gene}->{'tx'}->{$tx}->{'len'} += $len;
    }
}
close $PROTEOME or die "Can't close filehandle to proteome: $proteome\n";

open my $PSEG, '<', $seg or die "Can't open filehandle to PSEG FASTA of low-complexity domains: $seg\n";
while (my $input = <$PSEG>) { 
    chomp $input;
    if ( $input =~ /\A > /xms ) { 
        if ( $input =~ /\A > ((\S+\.g\d+)\.t\d+) \( (\d+) \- (\d+) \) /xms ) { 
            $tx       = $1;
            $gene     = $2;
            my $start = $3;
            my $stop  = $4;
            my $seg_span = $stop - $start + 1;
            $data_ref->{'gene'}->{$gene}->{'tx'}->{$tx}->{'seg_span'} += $seg_span;
        }
        else { 
            die "Can't parse: $input\n";
        }
    }
}
close $PSEG or die "Can't close filehandle to PSEG FASTA of low-complexity domains: $seg\n";

my @genes = grep { /\S/ } sort keys %{ $data_ref->{'gene'} };
my $header = "Gene\tPsegs\n";
foreach my $gene1 (@genes) { 
    my $annots_ref;
    my @annots = ();
    my @txs    = grep { /\S/ } sort keys %{ $data_ref->{'gene'}->{$gene1}->{'tx'} };
    foreach my $tx1 (@txs) { 
        if ( exists $data_ref->{'gene'}->{$gene1}->{'tx'}->{$tx1}->{'len'}) { 
            my $length   = $data_ref->{'gene'}->{$gene1}->{'tx'}->{$tx1}->{'len'};
            my $seg_span = 0;
            if ( exists $data_ref->{'gene'}->{$gene1}->{'tx'}->{$tx1}->{'seg_span'}) { 
                $seg_span = $data_ref->{'gene'}->{$gene1}->{'tx'}->{$tx1}->{'seg_span'};
            }
            if ( $seg_span > 0 ) { 
                my $orig_fraction = ($seg_span/$length);
                my $fraction = $orig_fraction;
                $fraction = sprintf "%.2f", $fraction;
                my $annot = "$fraction ($seg_span/$length)";
                $annots_ref->{'annot'}->{$annot} = $orig_fraction;
            }
        }
        else { 
            die "Failed to get length for protein product of $tx1\n";
        }
    }
    my $annot_text = q{};
    if (exists $annots_ref->{'annot'}) {
        @annots = grep { /\S/ } sort keys %{ $annots_ref->{'annot'} };
        if ($verbose) {  
            $annot_text = join '; ', @annots;
        }
        else { 
            $annot_text = $annots[-1];        
        }

        # Note that I earlier used:
        #     @annots = sort { $annots_ref->{'annot'}->{$b} <=> $annots_ref->{'annot'}->{$a} } grep { /\S/ } sort keys %{ $annots_ref->{'annot'} };
        #     $annot_text = $annots[0];
        # for the non-verbose option.  But it turns out that a straight ASCII-betical sort of the texts puts my preferred summary at the end anyway!
        # I.e., ASCII-betical picks the largest SEG%, and then the largest AA span if SEG% values tie.  Perfect.
    }
    print $header if $header;
    $header = q{};
    if ( $annot_text =~ /\S/xms ) { 
        $annot_text = 'SEG: ' . $annot_text;
    }
    print "$gene1\t$annot_text\n";
}

