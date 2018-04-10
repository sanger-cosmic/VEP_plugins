=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=head1 NAME

  DIAS

=head1 SYNOPSIS

  mv DIAS.pm ~/.vep/Plugins
  perl variant_effect_predictor.pl -i variations.vcf --cache --plugin DIAS

=head1 DESCRIPTION

  This is a plugin for the Ensembl Variant Effect Predictor (VEP) that
  adds extra fields to the output like AA_START, AA_STOP.
  
  This is written with the intention to be incorporated to the
  COSMIC annotation pipeline.

=cut

package DIAS;

use strict;
use warnings;
use Sanger::Cosmic::Dias::VEPAnnotationFormatter;
use Sanger::Cosmic::Dias::GenomeFormatter;
use Sanger::Cosmic::Dias::GenomicVariantDIAS;
use Data::Dumper;

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);

#--------------------------------------------------------------------------------#
sub new {
	my $class = shift;
    my $self = $class->SUPER::new(@_);
#	
#	# check config is OK
#	
#	# FASTA file defined, optimal
#	if (!defined($self->{config}->{fasta})) {
#	  
#	  # offline mode won't work without FASTA
#	  die("ERROR: Cannot generate HGVS without either a FASTA file (--fasta) or a database connection (--cache or --database)\n") if defined($self->{config}->{offline}) and !defined($self->{config}->{quiet});
#	  
#	  # cache mode will work, but DB will be accessed
#	  warn("WARNING: Database will be accessed using this plugin; use a FASTA file (--fasta) for optimal performance\n") if defined($self->{config}->{cache}) and !defined($self->{config}->{quiet});
#	}
#	
	if (!defined($self->{config}->{hgvs})) {
	  warn("WARNING: Plugin is enabling --hgvs\n") unless defined($self->{config}->{quiet});
	  $self->{config}->{hgvs} = 1;
	}
	return $self;
}
#--------------------------------------------------------------------------------#
sub version {
	return '90.0';
}
#--------------------------------------------------------------------------------#
sub feature_types {
	return ['Transcript', 'Intergenic'];
}
#--------------------------------------------------------------------------------#
sub variant_feature_types {
	return ['VariationFeature'];
	#return ['BaseVariationFeature'];
}
#--------------------------------------------------------------------------------#
sub get_header_info {
	my $self = shift;

	return {
		#RECORD 					=> 'input variant',
		CDS_START 				=> 'start position of transcript change in CDS coordinates',
		CDS_STARTOFFSET 		=> 'intron offset of CDS_START',
		CDS_STOP 				=> 'stop position of transcript change in CDS coordinates',
		CDS_STOPOFFSET 			=> 'intron offset of CDS_STOP',
		UTR_START 				=> 'UTR offset of CDS_START',
		UTR_STOP 				=> 'UTR offset of CDS_STOP',
		CDS_WT 					=> 'reference CDS sequence',
		CDS_MT 					=> 'mutant CDS sequence',
		CDS_SYNTAX 				=> 'CDS mutation syntax in HGVS format',
		AA_START 				=> 'start position of protein change in AA coordinates',
		AA_STOP 				=> 'stop position of protein change in AA coordinates',
		AA_WT					=> 'reference AA sequence',
		AA_MT 					=> 'mutant AA sequence',
		AA_SYNTAX 				=> 'AA mutation in HGVS syntax (1 letter code)',
		GENOME_START			=> 'start position of genomic change in genomic coordinates',
		GENOME_STOP	 			=> 'stop position of genomic change in genomic coordinates',
		GENOME_WT 				=> 'reference genomic sequence',
		GENOME_MT				=> 'mutant genomic sequence',
		GENOME_SYNTAX 			=> 'genomic mutation syntax in HGVS format',
		GENOME_VER 				=> 'assembly build',
		PERCENT_MUT_ALLELE 		=> 'Percentage of the mutant allele',
		#MUT_BURDEN 				=> '',
		CHR 					=> 'chromosome',
		STRAND 					=> 'strand of the transcript annotation',
		VARIANT_ONTOLOGY		=> 'SO accession for the variant type',
		CONSEQUENCES_ONTOLOGY 	=> 'SO accessions for consequence terms',
		ID_SAMPLE 				=> 'Cosmic sample',
		ID_STUDY 				=> 'Cosmic study',
		ID_PAPER 				=> 'Cosmic paper',
		USERNAME 				=> 'user inserting the mutation',
		GENE_NAME 				=> 'gene name',
		ACCESSION 				=> 'Ensembl stable_id of the transcript',
		DB 						=> 'source database of annotations',
		DBVERSION 				=> 'version of source database',
		#ID_VARIANT 				=> '',#TODO
		#ID_VARIANT_TYPE 		=> '',#TODO
		#ID_QUALITY 				=> '',#TODO
        ID_FEATURE_TYPE_COSMIC  => 'type of variant (coding=1/non-coding=4)',
		ID_MUT_SOMATIC_STATUS 	=> 'Cosmic somatic status',
		ID_MUT_VERIF_STATUS 	=> 'Cosmic verification status',
		ID_MUTATION_CURRENT 	=> 'current Cosmic mutation ID',
		#INTRA_INTERGENIC 		=> 'location of mutation on the genome relative to genes (intergenic/intragenic)',
		#WITHIN_GENE_FOOTPRINT 	=> 'located within a gene boundary (y/n)',
		ID_MUT_TYPE 			=> 'SO accession for the variant type',
		PARENT_MUT_LENGTH 		=> 'length of genomic mutant sequence',
		CDS_MUT_LENGTH 			=> 'length of CDS mutant sequence',
		AA_MUT_LENGTH 			=> 'length of AA mutant sequence',
		#NCV_REMARK 				=> 'remark specifying the cDNA location of mutation',
	};
}
#--------------------------------------------------------------------------------#
sub run {
	my ($self, $vfoa, $line_hash) = @_;
	#my ($self, $tva, $line_hash) = @_;
	my $input_var = get_variant_cosmic_data($line_hash->{Uploaded_variation});
	
	if ($vfoa->isa('Bio::EnsEMBL::Variation::IntergenicVariationAllele')) {
		my $genomic = get_genomic_data($vfoa);
		return {
			GENOME_START			=> $genomic->{START},
			GENOME_STOP	 			=> $genomic->{STOP},
			GENOME_WT 				=> $genomic->{WT},
			GENOME_MT				=> $genomic->{MT},
			GENOME_SYNTAX 			=> $line_hash->{HGVSg},
			GENOME_VER 				=> $self->{config}->{assembly},
			PERCENT_MUT_ALLELE 		=> $input_var->{percent_mut_allele} || '',
			CHR 					=> $genomic->{CHR},
			STRAND 					=> $genomic->{STRAND},
			ID_SAMPLE 				=> $input_var->{id_sample},
			ID_STUDY 				=> $input_var->{id_study} || '',
			ID_PAPER 				=> $input_var->{id_paper} || '',
			USERNAME 				=> $ENV{USER},
			DB 						=> $genomic->{DB},
			DBVERSION 				=> $genomic->{DBVERSION},
			ID_FEATURE_TYPE_COSMIC  => 4,
			ID_MUT_SOMATIC_STATUS 	=> $input_var->{confirmed},
			ID_MUT_VERIF_STATUS 	=> $input_var->{verified},
			ID_MUTATION_CURRENT 	=> $input_var->{id_mutation_current},
			ID_MUT_TYPE 			=> $genomic->{VARIANT_ONTOLOGY},
			PARENT_MUT_LENGTH 		=> $genomic->{MT} ne '-' ? length($genomic->{MT}) : '',
			VARIANT_ONTOLOGY		=> $genomic->{VARIANT_ONTOLOGY},		#TODO
			CONSEQUENCES_ONTOLOGY 	=> $genomic->{CONSEQUENCES_ONTOLOGY},	#TODO
		};
		
	} elsif ($vfoa->isa('Bio::EnsEMBL::Variation::TranscriptVariationAllele')) {
		my $tva = $vfoa;
		my $annotation_formatter = Sanger::Cosmic::Dias::VEPAnnotationFormatter->new(transcript_variation_allele => $tva,
																					 assembly_version => $self->{config}->{assembly});
		
		#TODO
		if (!defined $tva->hgvs_transcript) {	# Skip annotations where the variant lies outside the transcript cDNA
			return undef;				# these will not have a valid CDS syntax (i.e. no 'upstream_gene_variant' or 'downstream_gene_variant')
		}
		
		my $cds = $annotation_formatter->get_cds($tva);
		my $mrna = $annotation_formatter->get_mrna($tva);
		my $protein;
		
		if ($tva->transcript->translation) {	# If protein-coding transcript
			$protein = $annotation_formatter->get_protein($tva);
		}
		else {									# All other transcript-types (pseudogenes, lincRNAs etc)
			$protein = undef;
		}
		
		my $aa_mut = get_aa_mut_allele($protein, $line_hash);
		
		return {
			CDS_START 				=> $cds->{START},
			CDS_STARTOFFSET 		=> $cds->{STARTOFFSET} || '',
			CDS_STOP 				=> $cds->{STOP},
			CDS_STOPOFFSET 			=> $cds->{STOPOFFSET} || '',
			UTR_START 		 		=> $cds->{UTR_START} || '',
			UTR_STOP 	 			=> $cds->{UTR_STOP} || '',
			CDS_WT 					=> $cds->{WT} || '',
			CDS_MT 					=> $cds->{MT} || '',
			CDS_SYNTAX 				=> $cds->{SYNTAX},
			AA_START 				=> $protein->{START} || '',
			AA_STOP 				=> $protein->{STOP} || '',
			AA_WT					=> $protein->{WT} || '',
			AA_MT 					=> $aa_mut,
			AA_SYNTAX 				=> $protein->{SYNTAX} || '',
			GENOME_START			=> $tva->variation_feature->start,
			GENOME_STOP	 			=> $tva->variation_feature->end,
			GENOME_WT 				=> $cds->{WT} || '',
			GENOME_MT				=> $cds->{MT} || '',
			GENOME_SYNTAX 			=> $line_hash->{HGVSg},
			GENOME_VER 				=> $self->{config}->{assembly},
			PERCENT_MUT_ALLELE 		=> $input_var->{percent_mut_allele} || '',
			#MUT_BURDEN 				=> '',#TODO
			CHR 					=> $tva->variation_feature->seq_region_name,
			STRAND 					=> $cds->{STRAND},
			VARIANT_ONTOLOGY		=> $cds->{VARIANT_ONTOLOGY},
			CONSEQUENCES_ONTOLOGY 	=> $cds->{CONSEQUENCES_ONTOLOGY},
			ID_SAMPLE 				=> $input_var->{id_sample},
			ID_STUDY 				=> $input_var->{id_study} || '',
			ID_PAPER 				=> $input_var->{id_paper} || '',
			USERNAME 				=> $ENV{USER},
			GENE_NAME 				=> $cds->{GENE},
			ACCESSION 				=> $cds->{ACCESSION},
			CCDS 					=> $cds->{CCDS},
			DB 						=> $cds->{DB},
			DBVERSION 				=> $cds->{DBVERSION},
			#ID_VARIANT 				=> '',#TODO
			#ID_VARIANT_TYPE 		=> '',#TODO
			#ID_QUALITY 				=> '',#TODO
			ID_FEATURE_TYPE_COSMIC  => defined $cds->{SYNTAX} && $cds->{SYNTAX} =~ /^c\./ ? 1 : 4,
			ID_MUT_SOMATIC_STATUS 	=> $input_var->{confirmed},
			ID_MUT_VERIF_STATUS 	=> $input_var->{verified},
			ID_MUTATION_CURRENT 	=> $input_var->{id_mutation_current},
			ID_MUT_TYPE 			=> $cds->{VARIANT_ONTOLOGY},
			PARENT_MUT_LENGTH 		=> $cds->{MT} ne '-' ? length($cds->{MT}) : '',	#TODO - Do we need all 3 mut length fields?
			CDS_MUT_LENGTH 			=> $cds->{MT} ne '-' ? length($cds->{MT}) : '',
			AA_MUT_LENGTH 			=> defined $aa_mut ? length($aa_mut) : '',
		};
	}
}
#--------------------------------------------------------------------------------#
#TODO - Use Text::CSV_XS to parse this
sub get_variant_cosmic_data {
	my $input_csv = shift;
	my @cols = split(',', $input_csv);
	my $var = Sanger::Cosmic::Dias::GenomicVariantDIAS->new(chr 				=> $cols[0],
															genome_version		=> $cols[1],
															start				=> $cols[2],
															stop				=> $cols[3],
															wt 					=> $cols[4],
															mt 					=> $cols[5],
															percent_mut_allele 	=> $cols[6] || undef,
															confirmed 			=> $cols[7],
															verified			=> $cols[8] || undef,
															id_sample 			=> $cols[9],
															id_mutation_current	=> $cols[10],
															#id_study 			=> $cols[11] ne 'cosu' ? $cols[11] =~ s/cosu//gi : undef,	# empty values = 'cosu'
															#id_paper 			=> $cols[12] ne 'cosp' ? $cols[12] =~ s/cosp//gi : undef,	# empty values = 'cosp'
															id_study 			=> $cols[11] || undef,	# empty values = 'cosu'
															id_paper 			=> $cols[12] || undef,	# empty values = 'cosp'
															);
	return $var;
}
#--------------------------------------------------------------------------------#
sub get_genomic_data {
	my $vfoa = shift;
	my $genome_formatter = Sanger::Cosmic::Dias::GenomeFormatter->new(variation_allele => $vfoa);
	return {WT => $genome_formatter->wt_allele,
			MT => $genome_formatter->mut_allele,
			START => $genome_formatter->start,
			STOP => $genome_formatter->stop,
			CHR => $genome_formatter->chromosome,
			STRAND => $genome_formatter->strand,
			DB => $genome_formatter->db,
			DBVERSION => $genome_formatter->db_version,
			VARIANT_ONTOLOGY => $genome_formatter->variant_ontology,
			CONSEQUENCES_ONTOLOGY => $genome_formatter->consequences_ontology,
			};
}
#--------------------------------------------------------------------------------#
#VEP does not provide mutant aa for frameshift variants, so use Downstream column returned from the 'Downstream' plugin
sub get_aa_mut_allele {
	my ($protein, $line_hash) = @_;
	my $allele = $protein->{MT};
	if (/stop_gained/ ~~ $protein->{CONSEQUENCES_ONTOLOGY} || /frameshift_variant/ ~~ $protein->{CONSEQUENCES_ONTOLOGY}) {
		$allele = $line_hash->{DownstreamProtein}."*";	#adding the stop codon is cosmic notation
	}
	return $allele || undef;
}
#--------------------------------------------------------------------------------#
1;
