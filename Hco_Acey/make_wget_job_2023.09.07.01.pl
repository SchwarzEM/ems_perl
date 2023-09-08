#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use Cwd;

my $work_dir   = getcwd;

my $url_file   = q{};
my $login_data = q{};
my $file_list  = q{};

my $user       = q{};
my $password   = q{};
my $target_url = q{};

$url_file   = $ARGV[0] if $ARGV[0];
$login_data = $ARGV[1] if $ARGV[1];
$file_list  = $ARGV[2] if $ARGV[2];

if ( (! $url_file ) or (! $login_data ) or (! $file_list ) ) {
    die "Format: make_wget_job_2023.09.07.01.pl [file with target URL] [file with name and password] [file with file list]\n";
}

open my $LOGIN_DATA, '<', $login_data;
while ( my $input = <$LOGIN_DATA> ) {
    chomp $input;
    if ( $input =~ /\A user \t (.*) \z/xms ) {
        $user = $1;
    }
    elsif ( $input =~ /\A password \t (.*) \z/xms   ) {
        $password = $1;
    }
    else {
        die "From login data file $login_data, cannot parse: $input\n";
    }
}
close $LOGIN_DATA;

open my $URL_FILE, '<', $url_file;
while ( my $target = <$URL_FILE> ) {
    chomp $target;
    if ( (! $target_url ) and ( $target =~ /\A http \S+ \/ \z/xms ) ) {
        $target_url = $target;
    }
    elsif ( ( $target =~ /\S/xms ) and ( $target !~ /\A http \S+ \/ \z/xms ) ) {
        die "From URL file $url_file, cannot parse: $target\n";
    }
}
close $URL_FILE;

open my $FILE_LIST, '<', $file_list;
while ( my $file = <$FILE_LIST> ) {
    chomp $file;

    print "cd $work_dir ;\n" if $work_dir;
    $work_dir = q{};

    print "wget --no-verbose --no-clobber";
    if ( $user and $password ) {
        print " --user=$user --password=$password";
    }
    print " $target_url$file ;\n";
}

close $FILE_LIST;
