#!/usr/bin/env perl

# fimohits2beds.pl -- Erich Schwarz <emsch@its.caltech.edu>, 4/11/2010.
# Purpose: given a *.bed file, and a fimo hitlist of motifs against that *.bed's FASTA, generate bed subsets.

use strict;
use warnings;
use Getopt::Long;
use File::Basename;

my $subset_data;

my $bed_file       = q{};
my $fimo_file      = q{};
my $motif_set_name = q{};

my $mot_id      = q{};
my $ele_id      = q{};
my @motif_list  = ();
my %element_ids = ();

my $help;

GetOptions ( 'bed=s'    => \$bed_file,
             'fimo=s'   => \$fimo_file,
             'prefix=s' => \$motif_set_name,
             'help'     => \$help,           );

if ( $help 
     or (! (-r $bed_file  ) ) 
     or (! (-r $fimo_file ) ) ) { 
    die "Format: fimotxts2bed_subsets.pl",
        " --bed|-b [BED file]",
        " --fimo|-f [fimo.txt] file",
        " --prefix|-p [motif set name]",
        "\n",
        ;
}

# Read BED for the first time, just to get IDs for elements.
open my $BED1, '<', $bed_file or die "Can't open BED file $bed_file: $!";
while (my $input = <$BED1>) { 
    chomp $input;

    # Silently ignore lines starting with 'track name="', but require parsing of all others or death.
    if (      ( $input !~ /\A track \s+ name = \" /xms   )
          and ( $input =~ / \A (?: \S+ \t){3} (\S+) /xms ) ) { 
        $ele_id = $1;
        $element_ids{$ele_id} = 1;
    }
    else { 
        if ( $input !~ /\A track \s+ name = \" /xms ) { 
            die "In BED file $bed_file, can't parse input line: $input\n";
        }
    }
}
close $BED1 or die "Can't close handle to BED file $bed_file: $!";

# Read fimo.txt file to get each subset of elements:
open my $FIMO, '<', $fimo_file or die "Can't open FIMO text file $fimo_file: $!";
while (my $input = <$FIMO>) { 
    chomp $input;

    # Silently ignore lines starting with '#', but require parsing of all others or death.
    if (     ( $input !~ /\A \# /xms                            ) 
         and ( $input =~ /\A (\S+) \t [^\s]+_([^\s_]+) \t /xms ) ) { 
        $mot_id = $1;
        $ele_id = $2;
        if (! ( exists $element_ids{$ele_id} ) ) { 
            die "In FIMO file $fimo_file, can't identify element \"$ele_id\" in input line: $input\n";
        }

        # Record mapping of elements to motifs (can be one-to-many) for later BED reading and subset-producing.
        $subset_data->{'element'}->{$ele_id}->{'motif'}->{$mot_id} = 1;
    }
    elsif ( $input !~ /\A \# /xms ) { 
        die "In FIMO file $fimo_file, can't parse input line: $input\n";
    }
}
close $FIMO or die "Can't close handle to FIMO text file $fimo_file: $!";

# Read BED for second time, to build individual subset line lists.
open my $BED2, '<', $bed_file or die "Can't open BED file $bed_file: $!";
while (my $input = <$BED2>) {
    chomp $input;

    # This time, store any line starting with 'track name="' -- once, or death.
    if ( $input =~ /\A track \s+ name= /xms ) { 
        if ( exists $subset_data->{'bed_header'} ) { 
            die "Multiple header lines:\n",
                "Already stored: $subset_data->{'bed_header'}\n",
                "Cannot also store: $input\n",
                ;
        }
        $subset_data->{'bed_header'} = $input;
    }

    # Otherwise, require parsing of lines or death.
    if (      ( $input !~ /\A track \s+ name = \" /xms   )
          and ( $input =~ / \A (?: \S+ \t){3} (\S+) /xms ) ) {
        $ele_id = $1;
        if (! ( exists $element_ids{$ele_id} ) ) {
            die "In BED file $bed_file, can't identify element \"$ele_id\" in input line: $input\n";
        } 
        if ( exists $subset_data->{'element'}->{$ele_id} ) { 
            @motif_list = sort keys %{ $subset_data->{'element'}->{$ele_id}->{'motif'} };
            foreach my $motif1 (@motif_list) { 
                push @{ $subset_data->{'subset_bed_lines'}->{$motif1} }, $input;
            }
        }
    }
    elsif ( $input !~ /\A track \s+ name = \" /xms ) {
        die "In BED file $bed_file, can't parse input line: $input\n";
    }
}
close $BED2 or die "Can't close handle to BED file $bed_file: $!";

foreach my $motif2 (sort keys %{ $subset_data->{'subset_bed_lines'} } ) { 
    my $precursor_name = q{};
    if (! $motif_set_name) {
        $precursor_name = basename($bed_file);
        $precursor_name =~ s/\.bed\z//;
    }
    if ($motif_set_name) { 
        $precursor_name = $motif_set_name;
    }
    $precursor_name = $precursor_name . q{_motif.} . $motif2;
    my $bed_sub_name = $precursor_name . q{_subset.bed};
    $bed_sub_name = failsafe_name($bed_sub_name);
    open my $BED_SUBSET, '>', $bed_sub_name or die "Can't open BED subset file $bed_sub_name: $!";
    my $bedsub_header = $subset_data->{'bed_header'};
    if ( $bedsub_header =~ / \A track \s+ name = \" [^\"]+ \" \s+ description = \" [^\"]+ \" (.*) \z /xms ) { 
        $bedsub_header = $1;
        $bedsub_header =    qq{track name=\"$precursor_name} 
                          .  q{_elmts}
                          . qq{\" description=\"$precursor_name}
                          .  q{_elmts"}
                          . $bedsub_header 
                          ;
    }
    print $BED_SUBSET "$bedsub_header\n";
    foreach my $subset_bed_line ( @{ $subset_data->{'subset_bed_lines'}->{$motif2} } ) { 
        print $BED_SUBSET "$subset_bed_line\n";
    }
    close $BED_SUBSET or die "Can't close filehandle to BED subset file $bed_sub_name: $!";
}

sub failsafe_name {
    my $filename = $_[0];
    if (-e $filename) {
        my $suffix = 0;
        while (-e $filename) {
            $suffix++;
            $filename =~ s/\.\d+\z//xms;
            $filename = $filename . ".$suffix";
        }
    }
    return $filename;
}

