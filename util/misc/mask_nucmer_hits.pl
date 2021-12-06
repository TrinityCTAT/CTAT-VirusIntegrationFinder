#!/usr/bin/env perl

use strict;
use warnings;
use lib ($ENV{EUK_MODULES});
use Fasta_reader;
use Overlap_piler;


=notes

mummer-4.0.0rc1/nucmer ref_genome.fa Virus_db_Dec062021.nonUnq1kbMsk.filt90.fasta 2>&1 | tee nucmer.mine.log

mummer-4.0.0rc1/show-coords -T out.delta >  nucmer.coords


=cut



my $coords_file = "nucmer.coords";
my $fasta_file = "Virus_db_Dec062021.nonUnq1kbMsk.filt90.fasta";


my %acc_to_coords;

open(my $fh, $coords_file);
while(my $line = <$fh>) {
    chomp $line;
    unless ($line =~ /^\d/) {
        next;
    }
    my @x = split(/\t/, $line);
    
    
    my $viral_acc = $x[8];
    my $begin = int($x[2]);
    my $end = int($x[3]);
    ($begin, $end) = sort {$a<=>$b} ($begin, $end);
    
    push (@{$acc_to_coords{$viral_acc}}, [$begin, $end]);

}


my $fasta_reader = new Fasta_reader($fasta_file);
while (my $seq_obj = $fasta_reader->next()) {
    my $accession = $seq_obj->get_accession();
    
    if (exists $acc_to_coords{$accession}) {
        my @coords = @{$acc_to_coords{$accession}};
        @coords = &Overlap_piler::simple_coordsets_collapser(@coords);
        
        # mask them:
        my $sequence = $seq_obj->get_sequence();
        my @chars = split(//, $sequence);
        foreach my $coordset (@coords) {
            my ($lend, $rend) = @$coordset;
            for (my $i = $lend - 1; $i < $rend; $i++) {
                $chars[$i] = 'N';
            }
        }
        $seq_obj->{sequence} = join("", @chars);
        
        
    }

    my $fasta_record = $seq_obj->get_FASTA_format();
    print($fasta_record);
    
}
