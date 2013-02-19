#!/bin/env perl
use strict;
use warnings;
use autodie qw":all";
use PerlIO;
use Bio::Seq;
use Bio::SeqIO;
use Bio::Tools::GFF;
use Data::Dumper;
use Bio::Perl;

## Purpose:  Read in the Esmeraldo, Non-Esmeraldo, and Unassigned genome files/annotations
##   Create a library from them containing just the annotated begin .. end

my $start_padding = 0;
my $end_padding = 0;
my @genomes = ('esmer', 'non', 'unas');

foreach my $genome (@genomes) {
    print "Starting $genome\n";
    my $chromosomes = Read_Chromosomes($genome);
    Make_Library($chromosomes, $genome, $start_padding, $end_padding);
}

sub Read_Chromosomes {
    my $type = shift;
    my $chromosomes = {};
    my $input_fasta = qq"genome/${type}.fasta";
    my $input_genome = new Bio::SeqIO(-file => $input_fasta, -format => 'Fasta');
    while (my $genome_seq = $input_genome->next_seq()) {
	next unless(defined($genome_seq->id));
	my $id = $genome_seq->id;
# Get rid of the spurious 'GeneDB' crap
	$id =~ s/^GeneDB\|//g;
	$chromosomes->{$id} = $genome_seq;
    }
    return($chromosomes);
}


sub Make_Library {
    my $chromosomes = shift;
    my $type = shift;
    my $start_padding = shift;
    my $end_padding = shift;
    
    my $input_gff = qq"gff/${type}.gff";
    my $output_fasta = qq"${type}_annotated_genes.fasta";
    
    ## I don't think Bio::Tools::GFF has a -file accessor...
    open(FASTA, ">$output_fasta");
    open(GFF, "<$input_gff");
    my $input_annotation = new Bio::Tools::GFF(-fh => \*GFF, -gff_version => 3);
    # loop over the input stream
    #    print Dumper($input_annotation);
    #   my $test_num = 10;
  LOOP: while(my $feature = $input_annotation->next_feature()) {
      # do something with feature
      next LOOP unless ($feature->{_primary_tag} eq 'gene');
      ## The fields of $feature include:
      ## _gsf_frame .
      ## _gsf_seq_id TcChr17-S
      ## _gsf_tag_hash HASH(0x4c56958)
      ## _location Bio::Location::Simple=HASH(0x4c566a0)
      ##   Bio::Location::Simple has: ->strand ->start ->end ->length
      ## _parse_h HASH(0x4c56988)
      ## _primary_tag mRNA
      ## _root_cleanup_methods ARRAY(0x4c56bb0)
      ## _source_tag EuPathDB
      ## $input_annotation->gff_string($feature) gets access to the raw gff line
      ## $feature->all_tags(); gets all the annotation tags
      ## $feature->each_tag_value("ID"); gets the ID tag as an array
      
      my $location = $feature->{_location};
      
      my $start = $location->start();
      my $end = $location->end();
      my $strand = $location->strand();
      my @ids = $feature->each_tag_value("ID");
      my $id = "";
      my $gff_chr = $feature->{_gsf_seq_id};
      my $gff_string = $input_annotation->gff_string($feature);
      if (!defined($chromosomes->{$gff_chr})) {
	  print "Something is wrong with $gff_chr\n";
	  next LOOP;
      }
      my $gene_seq = $chromosomes->{$gff_chr}->subseq($start, $end);

      foreach my $i (@ids) {
	  $i =~ s/^cds_//g;
	  $i =~ s/\-\d+$//g;
	  $id .= "$i ";
      }

      my @gff_information = split(/\t+/, $gff_string);
      my $description_string = $gff_information[8];
#      print "$id";
#      if ($strand == -1) {
#	  $gene_seq = reverse($gene_seq);
#	  $gene_seq =~ tr/ATGCatgc/TACGtacg/;
#      }
      print FASTA ">$id $description_string
$gene_seq
";

  } ## End looking at every feature
    
    close(FASTA);
}  ## End Make_Library
