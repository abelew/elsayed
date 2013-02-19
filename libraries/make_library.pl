#!/usr/bin/perl
use strict;
use warnings;
use autodie qw":all";
use PerlIO;
use Bio::Seq;
use Bio::SeqIO;
use Bio::Tools::GFF;
use Data::Dumper;


## Purpose:  Read in the Esmeraldo, Non-Esmeraldo, and Unassigned genome files/annotations
##   Create a library from them containing all open reading frames with some small amount of UTR sequence
##   Fix a subset of errors noticed in the GFF files.
##   Perhaps generate a genbank flatfile version of these for future work so that I can extend the annotation information.

my $padding = 60;
my @genomes = ('esmer', 'non', 'unas');



foreach my $genome (@genomes) {
    print "Starting $genome\n";
    my $chromosomes = Read_Chromosomes($genome);
    Make_Library($chromosomes, $genome, $padding);
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
    my $padding = shift;

    my $input_gff = qq"gff/${type}.gff";
    my $output_fasta = qq"${type}_plus_${padding}.fasta";
    my $output_genbank = qq"genome/{$type}.genbank";

    ## I don't think Bio::Tools::GFF has a -file accessor...
    open(GFF, "<$input_gff");
#    my $input_annotation = new Bio::Tools::GFF(-file => $input_gff, -gff_version => 2);
    my $input_annotation = new Bio::Tools::GFF(-fh => \*GFF, -gff_version => 3);
    # loop over the input stream
#    print Dumper($input_annotation);

    my %counters = (
	plus_total => 0.0,
	minus_total => 0.0,
	total => 0.0,
	plus_aug => 0.0,
	minus_aug => 0.0,
	plus_stop => 0.0,
	minus_stop => 0.0,
	plus_perfect => 0.0,
	minus_perfect => 0.0,
	);
    LOOP: while(my $feature = $input_annotation->next_feature()) {
	# do something with feature
	next LOOP unless ($feature->{_primary_tag} eq 'CDS');
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

	foreach my $i (@ids) {
	    $i =~ s/^cds_//g;
	    $i =~ s/\-\d+$//g;
	    $id .= "$i ";
	}

	$counters{total}++;

	if (!defined($chromosomes->{$gff_chr})) {
	    print "Something is wrong with $gff_chr\n";
	    next LOOP;
	}
	my $initial_CDS = $chromosomes->{$gff_chr}->subseq($start, $end);
	if ($strand == -1) {
	    $counters{plus_total}++;
	    $initial_CDS = reverse($initial_CDS);
	    $initial_CDS =~ tr/ATGC/TACG/;
	} else {
	    $counters{minus_total}++;
	}
	my $start_codon = substr($initial_CDS, 0, 3);
	my $stop_codon = substr($initial_CDS, -3);

	if ($start_codon eq 'ATG') {
	    if ($strand == 1) {
		$counters{plus_aug}++;
		if ($stop_codon eq 'TAA' or $stop_codon eq 'TAG' or $stop_codon eq 'TGA') {
		    $counters{plus_perfect}++;
		}
	    } else {  ## strand -1
		$counters{minus_aug}++;
		if ($stop_codon eq 'TAA' or $stop_codon eq 'TAG' or $stop_codon eq 'TGA') {
		    $counters{minus_perfect}++;
		}
	    } ## End strand test
	} ## End testing for AUG

	if ($stop_codon eq 'TAA' or $stop_codon eq 'TAG' or $stop_codon eq 'TGA') {
	    if ($strand == 1) {
		$counters{plus_stop}++;
	    } else {
		$counters{minus_stop}++;
	    }
	} ## End checking stop codon

    } ## End looking at every feature

    ## Calculate numbers of ORFS
    my $plus_perfect_pct = ($counters{plus_perfect} / $counters{plus_total}) * 100.0;
    my $minus_perfect_pct = ($counters{minus_perfect} / $counters{minus_total}) * 100.0;
    my $plus_pct = ($counters{plus_total} / $counters{total}) * 100.0;
    my $plus_aug = ($counters{plus_aug} / $counters{plus_total}) * 100.0;
    my $minus_aug = ($counters{minus_aug} / $counters{minus_total})  * 100.0;
    print "Total ORFs on plus strand: $counters{plus_total}
Total ORFs on minus strand: $counters{minus_total}
Percent ORFs on plus strand: $plus_pct
Percent ORFs which start AUG on plus strand: $plus_aug
Percent ORFs which start AUG on minus strand: $minus_aug
Percent perfect ORFs on plus strand: $plus_perfect_pct
Percent perfect ORFs on minus strand: $minus_perfect_pct
";

}


# 	my $chr_length = $chromosomes->{$gff_chr}->length();
# 	if (($end + 243) > $chr_length) {
# 	    next LOOP;
# 	}
# 	if ($start < 240) {
# 	    next LOOP;
# 	}
# 	my $start_library_sequence = $chromosomes->{$gff_chr}->subseq($start - 240, $end + 240);
# 	if ($strand == -1) {
# 	    $start_library_sequence = reverse($start_library_sequence);
# 	    $start_library_sequence =~ tr/ATGC/TACG/;
# 	}
# 	my @tests = (0, 1, 2, -1, -2);
# 	my $library_sequence;
# 	TST: foreach my $test (@tests) {
# 	    my $test_library_sequence = Check_Frame($start_library_sequence, $test);
# 	    if (defined($test_library_sequence)) {
# 		$library_sequence = $test_library_sequence;
# 		last TST;
# 	    }
# 	}
# 	if (!defined($library_sequence)) {
# 	    $library_sequence = substr($start_library_sequence, 60, -60);
# 	}
# #	print "TESTME: <$start_codon> <$stop_codon>\n";
# 	print "TESTME: $strand $id 
# $library_sequence\n\n";
#     }
#     $input_annotation->close();
#     close(GFF);
#}

sub Check_Frame {
    my $seq = shift;
    my $frame = shift;
#    print "$seq\n";
    my $test_start = substr($seq, (240 + $frame), 3);
    my $test_stop = substr($seq, (-243 + $frame), 3);
    print "CheckTST: <$test_start> <$test_stop>\n";
    if ($test_start eq 'ATG' and ($test_stop eq 'TAA' or $test_stop eq 'TAG' or $test_stop eq 'TGA')) {
	print "Found ORF in $frame frame\n";
	return(substr($seq, (60 + $frame), (-60 + $frame)));
    } else {
	return(undef);
    }
}
