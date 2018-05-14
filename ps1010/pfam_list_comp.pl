#!/usr/bin/env perl 

# pfam_list_comp.pl -- Erich Schwarz <emsch@its.caltech.edu>, 11/11/2009.
# Purpose: filter all-present versus optional all-absent PFAM domains in a roster, given outputs of pfam_tbl2roster.pl

use strict;
use warnings;
use Getopt::Long;

my @neg_files = ();
my @pos_files = ();
my $key_file  = q{};

my %neg_names = ();

my %pos_names = ();
my $pos_level = 0;

GetOptions ( 'neg:s{,}' => \@neg_files,
             'pos=s{,}' => \@pos_files, 
             'key=s'    => \$key_file,  );

if ( (! @pos_files) or (! $key_file) ) { 
    die "Format: pfam_list_comp.pl",
        " --pos|-p [positive PFAM/gene rosters]",
        " --neg|-n [optional negative PFAM/gene rosters]",
        " --key|-k [key roster to be filtered]\n",
        ;
}

$pos_level = scalar(@pos_files);

if (@neg_files) {
    foreach my $neg_file (@neg_files) { 
        open my $NEG, '<', $neg_file or die "Can't open negative file $neg_file: $!";
        while (my $input = <$NEG>) { 
            chomp $input;
            if ( $input =~ /\A \d+ \t (PF\d+\S*) \t /xms ) { 
                my $name = $1;
                # This assumes a well-formed report, with no more than one line per PFAM-A domain:
                $neg_names{$name}++;
            }
        }
        close $NEG or die "Can't close filehandle to negative file $neg_file: $!";
    }
}

foreach my $pos_file (@pos_files) { 
    open my $POS, '<', $pos_file or die "Can't open positive file $pos_file: $!";
    while (my $input = <$POS>) {
        chomp $input;
        if ( $input =~ /\A \d+ \t (PF\d+\S*) \t /xms ) {
            my $name = $1;
            # This assumes a well-formed report, with no more than one line per PFAM-A domain:
            $pos_names{$name}++;
        }
    }
    close $POS or die "Can't close filehandle to positive file $pos_file: $!";
}

# Censor any names that weren't found in all positive files:
foreach my $pos_name (sort keys %pos_names) {
    if ( $pos_names{$pos_name} < $pos_level ) {
        delete $pos_names{$pos_name};
    }
}

open my $KEY, '<', $key_file or die "Can't open key roster $key_file: $!";
while (my $input = <$KEY>) { 
    chomp $input;
    if ( $input =~ /\A \d+ \t (PF\d+\S*) \t /xms ) { 
        my $name = $1;
        if ( ( exists $pos_names{$name} ) and (! exists $neg_names{$name} ) ) { 
            print "$input\n";
        }
    }
}
close $KEY or die "Can't close filehandle to key roster $key_file: $!";

