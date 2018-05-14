#!/usr/bin/perl

# move_many_dirconts.pl -- Erich Schwarz <emsch@its.caltech.edu>, 11/29/2007.
# Purpose: move files from 1+ dirs. to a target dir.; useful for MANY files (hi, Ali!).

use strict;
use warnings;
use File::Spec::Functions;

# Initialize directory variables.
my @directories = @ARGV;
my @target_dirs = ();
my $target_dir = q{};

# Initialize nonfatal error tracking:
my $no_sources = 1;
my $no_moves   = 1;

# Define source(s) and target directories, or die/complain.
# N.B.: grep returns an *array*, so:
@target_dirs 
    = grep { $_ =~ /\A \-\-target=\S+ \z/xms } @directories 
        or die_loudly();

# Then get scalar.  This is a kludge; make more elegant!
$target_dir = $target_dirs[0];

# Correct variables.
$target_dir =~ s/\A\-\-target=//;
@directories 
    = grep { $_ !~ /\A \-\-target=\S+ \z/xms } @directories;

# Check for fatal errors:
if (  (! $target_dir            ) 
   or (! dir_check($target_dir) ) 
   or (! @directories           ) ) {
    die_loudly();
}
if (! -w $target_dir ) { 
    die "Can't write to target directory $target_dir: $!";
}

# Initialize relocatable contents of source directories.
my %relocatables = ();

# This can mine out several members of @directories:
foreach my $start_dir (@directories) { 

    # Require a real, readable directory with dir_check:
    if ( dir_check($start_dir) ) { 
        opendir(my $DIR_START, $start_dir) 
            or die "Couldn't open directory $start_dir: $!";

        # Only if we get here, do we zero this flag:
        $no_sources = 0;
        my $next1 = q{};

        # Assume [1] that every named file is worth moving.
        while ( defined ( $next1 = readdir($DIR_START) ) ) { 
            # Exclude '.', '..'
            if ( $next1 =~ /\w/xms ) { 

                # Need to build a full path from $next1:
                my $next2 = catdir($start_dir, $next1);
                my $unwanted = catdir($target_dir, $next1);

                # Move, if everything's OK:
                if (  (  -w $next2      ) 
                  and (  -w $target_dir ) 
                  and (! -e $unwanted   ) ) { 
                    system "mv -u $next2 $target_dir";
                    # Only now, zero this flag:
                    $no_moves = 0;
                }
            }
        }
        closedir($DIR_START);
    }
}

# Announce any silent failures:
if ($no_sources) { 
    print "None of the putative source directories were usable!\n";
}
if ( (! $no_sources) and ($no_moves) ) {
    print "Nice directories ... but nothing in them was movable!\n";
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

sub die_loudly {
    die "Format: ",
        "  ./move_many_motdirs.pl",
        "  --target=[target_dir]",
        "  [source_dir1] ... [source_dirN]\n",
    ;
}

