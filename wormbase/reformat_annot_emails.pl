#!/usr/bin/env perl

# reformat_annot_emails.pl -- Erich Schwarz <emsch@its.caltech.edu>, 8/25/2009.
# Purpose: take a blob of e-mailed annots, and reformat it into a useful textfile.

use strict;
use warnings;
use Text::Wrap qw($columns &wrap);

my $print_stream    = 1;
my $print_this_line = 1;
my $columns         = 78;

my %seen_lines  = ();
my %seen_wbrefs = ();

my $text_breaker = "\n"
                   . '************************************************************'
                   . "\n"
                   ;

if (! $ARGV[0]) { 
    die "Format: ./reformat_annot_emails.pl [ugly e-mail blobfile] > [less-ugly text]\n";
}

while (my $input = <>) { 
    chomp $input;

    # Two separate issues.
    #     Do I want to just ignore a long stream of boilerplate?
    #     Assuming it's not boilerplate, do I still want to skip *this* line?
    # The first decision can be made in long blocks.
    # The second decision has to be made line-by-line.
    #     So, *default* the second decision to 'yes' in each loop:
    $print_this_line = 1;

    # Skip over streams of longish e-mail boilerplate:
    if ( ( $input =~ / \A From \s+ postgres\@tazendra.caltech.edu /xms ) 
         or ( $input =~ / \A From \s+ genenames-bounces\@brie4.cshl.edu /xms ) ) { 
        $print_stream = 0;
    }

    # Resume printing streams at subject line:
    if  ( $input =~ / \A Subject: \s /xms ) { 
        $print_stream = 1;
    }

    # But do not ever print exactly the same input line twice:
    if ($seen_lines{$input}) {
        $print_this_line = 0;
    }
    if (! $seen_lines{$input}) { 
        $seen_lines{$input} = 1;
    }

    # Put a text-breaker if the WBPaper is new:
    if ( ( $input =~ / \A Subject: \s+ (WBPaper\d+) /xms ) 
         and (! $seen_wbrefs{$1} ) ) { 
        my $papername = $1;
        print "$text_breaker\n";
        $seen_wbrefs{$papername} = 1;
    }
    
    # Enforce maximum line lengths of 78 characters:
    if ( length($input) > 78 ) { 
        $input = wrap(q{}, q{}, $input) . "\n";
    }

    # Still OK to print?  Go for it!
    if ( ($print_stream) and ($print_this_line) ) { 
        print "$input\n";
    }
}

