#!/usr/bin/env perl

# patch_defective_rnaiace.pl -- Erich Schwarz <emsch@its.caltech.edu>, 11/23/2010.
# Purpose: patch .ace which has a handful of errors that nevertheless should be fixed.

use strict;
use warnings;
use Getopt::Long;

my $object      = q{};
my $text        = q{};
my $ok_to_print = 1;  
my $lines_added = 0;  

my $acefile          = q{};
my @patch_files      = ();
my %objects_to_patch = ();

my $delete_ref;
my $add_ref;   
my $help;

GetOptions ( 'acefile=s'       => \$acefile,
             'patchfiles=s{,}' => \@patch_files,
             'help'            => \$help, );

if ($help or (! $acefile) or (! @patch_files) ) { 
    die "Format: patch_defective_rnaiace.pl --patchfiles|-p [1+ patch .ace files] --acefile|-a  <input .ace file or stream>\n";
}

foreach my $patch_file (@patch_files) { 
    open my $PATCH_FILE, '<', $patch_file or die "Can't open patch file $patch_file: $!";
    while (my $input = <$PATCH_FILE>) { 
        chomp $input;
        if ( $input =~ / \A [^:\s\"]+ \s+ : \s \" ([^\"]+) \" /xms ) {
            $object = $1;
        }
        elsif ( $input =~ /\A \-D \s+ (.+) \s* \z/xms ) { 
            $text = $1;
            # To enforce possibility of later text-match, trim dangling empty spaces.
            $text =~ s/\s*\z//;
            $delete_ref->{$object}->{$text} = 1;
            $objects_to_patch{$object} = 1;
        }
        elsif ( $input =~ /\A(\S.*)\s*\z/xms ) {
            $text = $1;
            $add_ref->{$object}->{$text} = 1;
        }
    }
    close $PATCH_FILE or die "Can't close filehandle to patch file $patch_file: $!";
}

open my $ACEFILE, '<', $acefile or die "Can't open .ace file $acefile: $!";
while (my $input = <$ACEFILE>) { 
    chomp $input;
    if ( $input =~ / \A [^:\s\"]+ \s+ : \s \" ([^\"]+) \" /xms ) { 
        $object = $1;
        $lines_added = 0;
    }
    if (! $objects_to_patch{$object} ) { 
        print "$input\n";
    }
    if ( $objects_to_patch{$object} ) { 
        $ok_to_print = 1;
        # Again, to make text-match at least possible, trim dangling empty spaces.
        $input =~ s/\s*\z//;
        foreach my $reject ( sort keys %{ $delete_ref->{$object} } ) { 
            if ( $input eq $reject ) { 
                $ok_to_print = 0;
                # Print wanted text exactly once per object:
                if (! $lines_added ) { 
                    foreach my $patchline ( sort keys %{ $add_ref->{$object} } ) { 
                        print "$patchline\n";
                    }
                    $lines_added = 1;
                }
            } 
        }
        if ($ok_to_print) { 
            print "$input\n";
        }
    }
}
close $ACEFILE or die "Can't close filehand to .ace file $acefile: $!";

