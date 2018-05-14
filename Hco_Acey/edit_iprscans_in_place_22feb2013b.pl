#!/usr/bin/env perl

use strict;
use warnings;

while (my $main_script = <>) { 
    chomp $main_script;
    my $backup_script = $main_script . '.backup';
    my $inline_perl_text_patch = q{ perl -ne ' s/ppn\=4/ppn=8/; s/mem\=10gb\n/mem=10gb\n\#PBS -l feature=intel10\n/; print ' };
    system "    mv -i $main_script $backup_script ;\n";
    system "    cat $backup_script | $inline_perl_text_patch > $main_script ;\n";
}

