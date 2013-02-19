#!/usr/bin/perl
use strict;
use warnings;
use autodie qw":all";

use Bio::Seq;
use Bio::SeqIO;
use Bio::Tools::GFF;

my @hypo = Get_Hypothetical_IDs();
Extract_Hypothetical("high_copy.fasta", \@hypo);
Extract_Hypothetical("low_copy.fasta", \@hypo);


sub Get_Hypothetical_IDs {
    my @hypothetical_ids = ();
    open(HYP, "<hypothetical.gff");
    my $input_hypo = new Bio::Tools::GFF(-fh => \*HYP, -gff_version => 3);
  LOOP: while(my $feature = $input_hypo->next_feature) {
      push(@hypothetical_ids, $feature->each_tag_value("ID"));
  }
    close(HYP);
    return(@hypothetical_ids);
}


sub Extract_Hypothetical {
    my $filename = shift;
    my $ids = shift;
    my $input = new Bio::SeqIO(-file => "<$filename", -format => 'Fasta');
    my $out_hypo = new Bio::SeqIO(-file => ">${filename}-hypo", -format => 'Fasta');
    my $out_nothypo = new Bio::SeqIO(-file => ">${filename}-nohypo", -format => 'Fasta');

#    my $hypo_p = 0;
  LOOP: while (my $in_seq = $input->next_seq()) {
      my $input_id = $in_seq->id;
      foreach my $id (@{$ids}) {
	  if ($id eq $input_id) {
#	      $hypo_p++;
	      $out_hypo->write_seq($in_seq);
	      next LOOP;
	  }
      } ## End foreach
      ## If we get here, then it should not be hypothetical
      $out_nothypo->write_seq($in_seq);
  } ## End each sequence of the library
}
