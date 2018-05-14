#!/usr/bin/env perl

# recurs_dirs_777_to_775.pl -- Erich Schwarz <emsch@its.caltech.edu>, 4/12/2009.
# Purpose: given old and new permission states, recurs. change old->new in dir.: default is 777->755.

use strict;
use warnings;
use Getopt::Long;
use File::Find;
use Cwd;

my $directory = q{};
my $old_perms = 0777;
my $new_perms = 0755;

# N.B.: commenting the permission flags out, for now, until I figure out inputs.

GetOptions ( "dir=s" => \$directory,
#            "old=s" => \$old_perms,
#            "new=s" => \$new_perms,
           );

$old_perms = $old_perms & 07777;
$new_perms = $new_perms & 07777;

if (! $directory) { 
    die "Format: recurs_dirs_777_to_775.pl --dir|-d [directory] \n";
}

switch_perms($directory);

find(\&switch_perms_of_dirs, ($directory));

sub get_mode { 
    my $_file = $_[0];
    # Bitwise, octal-aware wackiness from Programming Perl, ch. 29, p. 746:
    my $_mode = (stat $_file)[2] & 07777;
    if (! $_mode ) { 
        warn "get_mode failed to get a mode for file $_file!\n";
        return;
    }
    return $_mode;
}

sub switch_perms { 
    my $_directory = $_[0];
    if (-d $_directory) {
        my $mode = get_mode($_directory);
        if ($mode == $old_perms) {
            chmod($new_perms, $_directory)
                or warn "Can't change mode of",
                        " directory $_directory",
                        " from $old_perms to $new_perms!\n",
                        ;
        }
    }
    return;
}

sub switch_perms_of_dirs { 
    my $dir = getcwd;
    opendir my $DIR, $dir or die "Can't open handle for my own directory $dir: $!";
    my @files = readdir $DIR;
    closedir $DIR or die "Can't close handle for my own directory $dir: $!";
    foreach my $file (@files) { 
        switch_perms($file);
    }
    return;
}
