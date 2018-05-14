#!/usr/bin/env perl
    
use strict;
use warnings;

# Not 'foreach'! which only lets the first script get fixed.
while (my $main_script = <>) {
    chomp $main_script;
    my $backup_script = $main_script . '.backup';
    my $inline_perl_text_patch = q{ perl -ne ' s/08:00:00/12:00:00/; s/mem=7gb/mem=20gb/; print ' };
    system "    mv $main_script $backup_script ;\n";
    system "    cat $backup_script | $inline_perl_text_patch > $main_script ;\n";
}   

