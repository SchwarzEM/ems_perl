#!/usr/bin/perl

# list_all_mot_files.pl -- Erich Schwarz <emsch@its.caltech.edu>, 11/29/2007.
# Purpose: extract out a very long list of Ali's motif files.

use strict;
use warnings;
use File::Spec::Functions;

my %files = ();

# This can mine out several directories, listed in @ARGV:
foreach my $start_dir (@ARGV) { 

    # Require a real, readable directory with dir_check:
    if ( dir_check($start_dir) ) { 
        opendir(my $DIR_START, $start_dir) 
            or die "Couldn't open directory $start_dir: $!";
        my $next1 = q{};

        # Assume [1] that every named file is a directory.
        while ( defined ( $next1 = readdir($DIR_START) ) ) { 
            # Exclude '.', '..'
            if ( $next1 =~ /\w/xms ) { 

                # Need to build a full path from $next1:
                my $next2 = catdir($start_dir, $next1, "tabs");

                # But test assumption [1] here with dir_check:
                if ( dir_check($next2) ) { 
                    opendir(my $TABS_DIR, $next2) 
                        or die "Couldn't open directory $next2: $!";
                    my $next3 = q{};

                    # Now get only the '-ce.1.tab' files:
                    while ( defined ( $next3 = readdir($TABS_DIR) ) ) {
                        if ( $next3 =~ /\-ce\.1\.tab\z/xms ) {
                            $next3 = catfile($next2, $next3);
                            $files{$next3} = 1;
                        }
                    }
                    closedir($TABS_DIR);
                }
            }
        }
        closedir($DIR_START);
    }
}

foreach my $list_file (sort keys %files) {
    print "$list_file\n";
}

sub dir_check { 
    my $_dir_query = $_[0];
    if ( (-d $_dir_query) and (-r $_dir_query) ) { 
        return 1;
    }
    else {
        return 0;
    }
}

