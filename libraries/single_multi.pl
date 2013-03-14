#!/usr/bin/env perl
use strict;
use warnings;
use Bio::Seq;
use Bio::SearchIO::blast;
use Bio::Tools::Run::StandAloneBlast;
use autodie;
## Sorted library provided on the command line by copy number according to BLAST searches
## The criteria for a hit are not very strong: (E: 0.00001  Score: 100)
## This assumes an argument which is the name of one of the BLAST databases
## esmer, non, unas, esmer-non, esmer-unas, non-unas, all
## It will then compare the provided library against all possible libraries

## This script assumes a name 'esmer' 'non' 'unas' 'esmer-non' ... 'all'
my @possible_libraries = ('esmer','non','unas','esmer-non','esmer-unas','non-unas','all');
my $query = $ARGV[0];
$query = qq"${query}.fasta";
## Search libraries are the various combinations of esmer, nonesmer, and unassigned which are of potential interest
## thus all 7 of the possible states
## These libraries are generated in make_blast.sh
my @search_libraries = ('esmer',
			'non', 
			'unas',
			'esmer-non',
			'esmer-unas',
			'non-unas',
			'all',);
## States to save: singles, doubles, few (3-10), many (10+)
my $blast_output = new Bio::SearchIO(-format => 'blast',);
for my $library (@search_libraries) {
    my $query_library = new Bio::SeqIO(-file => $query, -format => 'Fasta');
    my $output_directory = qq"${query}_vs_${library}";
    if (!-d $output_directory) {
	system("mkdir $output_directory");
    }
    open(COUNTS, ">${output_directory}/counts.txt");
    open(ZEROS, ">${output_directory}/zeros.fasta");
    open(SINGLES, ">${output_directory}/singles.fasta");
    open(DOUBLES, ">${output_directory}/doubles.fasta");
    open(FEW, ">${output_directory}/few.fasta");
    open(MANY, ">${output_directory}/many.fasta");
    
    my $seq_count = 0;
    my @params = (-program => 'blastn', -I => 't', -database => $library);
    while (my $query_seq = $query_library->next_seq()) {
	$seq_count++;
	my $id = $query_seq->id;
	my $desc = $query_seq->desc;
	my $seq = $query_seq->seq;
	my $search = new Bio::Tools::Run::StandAloneBlast(@params);
	my $blast_output;
	eval {
	    $blast_output = $search->blastall($query_seq);
	};
	if ($@) {
	    print "Error? $@\n";
	}
	my $result_count = 0;
	while (my $result = $blast_output->next_result()) {
	    $result_count++;
	    my $query_name = $result->query_name();
	    my $query_length = $result->query_length();
	    my $query_descr = $result->query_description();
	    my $stats = $result->available_statistics();
	    my $hits = $result->num_hits();
	    
	    my $hit_count = 0;
	    my $score_cutoff = 100;
	    my $sig_cutoff = 0.00001;
	  HITLOOP: while (my $hits = $result->next_hit()) {
	      my $hit_name = $hits->name();
	      my $hit_length = $hits->length();
	      my $hit_acc = $hits->accession();
	      my $hit_descr = $hits->description();
	      my $hit_score = $hits->score();
	      my $hit_sig = $hits->significance();
	      my $hit_bits = $hits->bits();
	      next HITLOOP unless ($hit_sig < $sig_cutoff and $hit_score > $score_cutoff);
	      $hit_count++;
	  }
	    my $fasta_entry = ">$id $desc
$seq
";
	    if ($hit_count == 1) {
		print SINGLES $fasta_entry;
	    } elsif ($hit_count == 2) {
		print DOUBLES $fasta_entry;
	    } elsif ($hit_count >= 3 and $hit_count <= 10) {
		print FEW $fasta_entry;
	    } elsif ($hit_count > 10) {
		print MANY $fasta_entry;
	    } else {
		print ZEROS $fasta_entry;
	    }
	    
	    print COUNTS "$id\t$hit_count\n";
	} ## End of the individual blast search
    } ## End of this individual sequence
} ## End of this search library  (the 7 states to search against)
close(COUNTS);
close(ZEROS);
close(SINGLES);
close(DOUBLES);
close(FEW);
close(MANY);

