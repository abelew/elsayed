#!/usr/bin/env perl
use strict;
use warnings;
use Bio::Seq;
use Bio::SearchIO::fasta;
use autodie;


## Need to read in the sequences from the search library
## Because ggsearch is dropping a portion of the description
my %descriptions_by_accession = ();
use Bio::SeqIO;
my $library = new Bio::SeqIO(-file => 'non.fasta', -format => 'Fasta');
while (my $query_seq = $library->next_seq()) {
    my $id = $query_seq->id;
    my $desc = $query_seq->desc;

    my @desc_array = split(/\;/, $desc);

    foreach my $element (@desc_array) {
    my $description = '';
	if ($element =~ /^description/) {
	    my ($d, $value) = split(/=/, $element);
	    $description = $value;
	    $description =~ s/%2C /, /g;
	    $descriptions_by_accession{$id} = $description;
	    last;
	}
    }
}



# format can be 'fasta', 'blast', 'exonerate', ...
my $searchio = new Bio::SearchIO(-format => 'fasta', -file => 'ggout2.txt', -best => 1, -signif => 0.0001);
print "Query Name\tQuery length\tHit ID\tHit Description\tHit Length\tScore\tE\tIdentity\tHit length\tHit Matches\n";
while (my $result = $searchio->next_result()) {
    while(my $hit = $result->next_hit) {
	my $query_name = $result->query_name();
	my $query_length = $result->query_length();
	my $accession = $hit->accession();
	my $acc2 = $hit->name();
	my $length = $hit->length();
	my $score = $hit->raw_score();
	my $sig = $hit->significance();
	my $ident = $hit->frac_identical();


	my $hit_len;
	my $hit_matches;
	while (my $hsp = $hit->next_hsp) {
	    $hit_len = $hsp->length('total');
	    my @matches = $hsp->matches(-seq => 'hit');
	    $hit_matches = $matches[1];
	}

	print "$query_name\t$query_length\t$acc2\t$descriptions_by_accession{$acc2}\t$length\t$score\t$sig\t$ident\t$hit_len\t$hit_matches\n";
	# process the Bio::Search::Hit::HitI object
#	while(my $hsp = $hit->next_hsp) {
	    # process the Bio::Search::HSP::HSPI object
#	}
    }
}
