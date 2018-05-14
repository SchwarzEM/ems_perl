#!/usr/bin/env perl

# sequin_pep2tbl_v01.pl -- Erich Schwarz <emsch@its.caltech.edu>, 10/6/2008.
# Purpose: for each two CONTIG.pep and CONTIG.tbl.part files, make a CONTIG.tbl file.

use strict;
use warnings;
use Cwd;

my $usable_contigs_ref;
my $annots_ref;

# These four variables all cleared by zero_out_stuff():
my $gene      = q{};
my $protein   = q{};
my $prot_desc = q{};
my $waiting   = 0;

my $dir = getcwd;
opendir my $DIR, $dir 
    or die "Can't open filehandle to current directory $dir: $!";

foreach my $file ( readdir($DIR) ) { 
    if ($file =~ / \A ([^\.]+) \. (pep|tbl\.part) \z /xms ) {
        my $contig = $1;
        my $suffix = $2;
        $usable_contigs_ref->{$contig}->{$suffix} = 1;
    }
}
closedir $DIR
    or die "Can't close filehandle to current directory $dir: $!";

foreach my $contig (sort keys %{ $usable_contigs_ref } ) { 
    zero_out_stuff();
    if ( ( exists $usable_contigs_ref->{$contig}->{'pep'} ) 
         and ( exists $usable_contigs_ref->{$contig}->{'tbl.part'} ) ) { 
        open my $PEP, '<', "$contig.pep" 
            or die "Can't open Sequin-pep file $contig.pep: $!";
        while ( my $input = <$PEP> ) { 
            chomp $input;

            # E.g.:
            # [protein=Cbre_JD01.001]  [gene=lcl|Cbre_JD01.001]  [prot_desc=Orthologous to C. elegans WBGene00016652|C44E4.3.]
            if ($input =~ /\A > \S+ .* \[protein=(.+)\] .* \[gene=(.+)\] .* \[prot_desc=(.+)\] \z/xms ) {
                $protein   = $1;
                $gene      = $2;
                $prot_desc = $3;
                $annots_ref->{$contig}->{$gene}->{'protein'}   = $protein;
                $annots_ref->{$contig}->{$gene}->{'prot_desc'} = $prot_desc;
            }
            elsif ( $input =~ /\A > /xms ) { 
                die "Can't parse header line: $input!\n";
            }
        }
        close $PEP or die "Can't close filehandle to Sequin-pep file $contig.pep: $!";
        zero_out_stuff();

        open my $PART, '<', "$contig.tbl.part"
            or die "Can't open partial .tbl input file $contig.tbl.part: $!";
        open my $TBL, '>', "$contig.tbl" 
            or die "Can't open complete .tbl output file $contig.tbl: $!";

        while ( my $input = <$PART> ) { 
            chomp $input;
            if ( $waiting and ( $input =~ / \A (?:<|>)? \d+ \t (?:<|>)? \d+ \t gene \z /xms ) ) { 
                $protein   = $annots_ref->{$contig}->{$gene}->{'protein'};
                $prot_desc = $annots_ref->{$contig}->{$gene}->{'prot_desc'};
                print {$TBL} "\t\t\t", 'product',   "\t", $protein,   "\n";
                print {$TBL} "\t\t\t", 'prot_desc', "\t", $prot_desc, "\n";
            }
            if ( $input =~ /\A \t \t \t gene \t (\S+) \z/xms ) { 
                $gene = $1;
            }
            if ( $input =~ / \A \d+ \t \d+ \t CDS \z/xms ) { 
                $waiting = 1;
            }
            print {$TBL} "$input\n";
        }
        $protein   = $annots_ref->{$contig}->{$gene}->{'protein'};
        $prot_desc = $annots_ref->{$contig}->{$gene}->{'prot_desc'};
        print {$TBL} "\t\t\t", 'product',   "\t", $protein,   "\n";
        print {$TBL} "\t\t\t", 'prot_desc', "\t", $prot_desc, "\n";
        close $PART or die "Can't close filehandle to Sequin-pep file $contig.tbl.part: $!";
    }
}

sub zero_out_stuff { 
    $gene      = q{};
    $protein   = q{};
    $prot_desc = q{};
    $waiting   = 0;
}

