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
  vep -i variations.tsv --cache --plugin DIAS

=head1 DESCRIPTION

  This is a plugin for the Ensembl Variant Effect Predictor (VEP) that
  adds extra fields to the output like AA_START, AA_STOP.
  
  This is written with the intention of being incorporated into the
  COSMIC annotation pipeline.

=cut

package DIAS;

use strict;
use warnings;
use Sanger::Cosmic::Dias::VEPAnnotationFormatter;
use Sanger::Cosmic::Dias::GenomeFormatter;
use Sanger::Cosmic::Dias::GenomicVariantDIAS;
use Bio::EnsEMBL::VEP::Utils qw(get_version_data);
use Sanger::Cosmic::Dias::Constants qw(%HEADER_DESCRIPTIONS %DEFAULT_COLUMN_VALUES);
use Try::Tiny;
use Data::Dumper;
use Hash::Util qw(lock_keys);
use Hash::Merge qw(merge);

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);

#use constant ID_FEATURE_TYPE_CODING => 1;		# variants within the boundaries (i.e. UTR + CDS + introns) of a CODING transcript
#use constant ID_FEATURE_TYPE_NONCODING => 4;	# variants within the boundaries (i.e. UTR + CDS + introns) of a NONCODING transcript
#use constant ID_FEATURE_TYPE_INTERGENIC => 5;	# variants within promoter regions of genes or intergenic regions

#--------------------------------------------------------------------------------#
#TODO - No need to check for hgvs option?
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
	return '92.0';
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
	return \%Sanger::Cosmic::Dias::Constants::HEADER_DESCRIPTIONS;
}
#--------------------------------------------------------------------------------#
sub run {
	my ($self, $vfoa, $line_hash) = @_;
	#my ($self, $tva, $line_hash) = @_;
	
	my $input_var;
	my $success = try {
		$input_var = get_input_variant_data($line_hash->{Uploaded_variation});
		1;
	} catch {
		warn "ERROR : $_\n";
		return 0;
	};
	return undef unless $success; 		# an undef returned from plugin->run() will filter that line from the output
	
	my $genomic = $self->get_genomic_data($vfoa);
	
	my %default_data = %Sanger::Cosmic::Dias::Constants::DEFAULT_COLUMN_VALUES;
	#TODO - lock_keys() provides no advantage because we merge
	lock_keys(%default_data);	
	my $merge_hash = Hash::Merge->new('RIGHT_PRECEDENT');
	
	if ($vfoa->isa('Bio::EnsEMBL::Variation::IntergenicVariationAllele')) {
		my %intergenic_data = (
			GENOME_START			=> $genomic->{START},
			GENOME_STOP	 			=> $genomic->{STOP},
			GENOME_WT 				=> $genomic->{WT},
			GENOME_MT				=> $genomic->{MT},
			GENOME_SYNTAX 			=> $line_hash->{HGVSg},
			HGVSg_OFFSET			=> $genomic->{HGVS_OFFSET},
			GENOME_VER 				=> $self->{config}->{assembly},
			PERCENT_MUT_ALLELE 		=> $input_var->percent_mut_allele,
			CHR 					=> $genomic->{CHR},
			STRAND 					=> $genomic->{STRAND},	#VEP produces a 'STRAND' column, so no need for this
			ID_SAMPLE 				=> $input_var->id_sample,
			ID_STUDY 				=> $input_var->id_study,
			ID_PAPER 				=> $input_var->id_paper,
			USERNAME 				=> $ENV{USER},
			DB 						=> $genomic->{DB},
			DBVERSION 				=> $genomic->{DBVERSION},
			ID_FEATURE_TYPE_COSMIC  => $Sanger::Cosmic::Dias::Constants::ID_FEATURE_TYPE_INTERGENIC,
			ID_MUT_SOMATIC_STATUS 	=> $input_var->confirmed,
			ID_MUT_VERIF_STATUS 	=> $input_var->verified,
			#ID_MUTATION_CURRENT 	=> join(';', @{$input_var->id_mutation_current}),
			ID_MUTATION_CURRENT 	=> $input_var->id_mutation_current,
			ID_MUT_TYPE 			=> $genomic->{VARIANT_ONTOLOGY},
			PARENT_MUT_LENGTH 		=> defined $genomic->{MT} ? length($genomic->{MT}) : '',
			VARIANT_ONTOLOGY		=> $genomic->{VARIANT_ONTOLOGY},
			CONSEQUENCES_ONTOLOGY 	=> $genomic->{CONSEQUENCES_ONTOLOGY},
			ANNOTATOR_VERSION 		=> $self->get_annotator_version,
		);
		
		my $data = $merge_hash->merge(\%default_data, \%intergenic_data);
		return $data;
		
	} elsif ($vfoa->isa('Bio::EnsEMBL::Variation::TranscriptVariationAllele')) {
		my $tva = $vfoa;
		my $annotation_formatter = Sanger::Cosmic::Dias::VEPAnnotationFormatter->new(transcript_variation_allele => $tva,
																					 assembly_version => $self->{config}->{assembly});
		
		my $cds = $annotation_formatter->get_cds($tva);
		my $mrna = $annotation_formatter->get_mrna($tva);
		my $protein = undef;
		
		if ($tva->transcript->translation) {	# If protein-coding transcript then populate. Leave empty for all other transcript-types (pseudogenes, lincRNAs etc)
			$protein = $annotation_formatter->get_protein($tva);
		}
		
		my $aa_mut = get_aa_mut_allele($protein, $line_hash);
		my $id_feature_type_cosmic = get_id_feature_type($cds);
			
		#TODO - $genomic->{WT} eq $cds->{WT}?
	
		my %transcript_data = (
			CDS_START 				=> $id_feature_type_cosmic == $Sanger::Cosmic::Dias::Constants::ID_FEATURE_TYPE_CODING ? $cds->{START} : undef,		# populate only for protein-coding transcripts
			CDNA_START 				=> $id_feature_type_cosmic == $Sanger::Cosmic::Dias::Constants::ID_FEATURE_TYPE_NONCODING ? $cds->{START} : undef,	# populate only for non-coding transcripts
			CDS_STARTOFFSET 		=> $cds->{STARTOFFSET},
			CDS_STOP 				=> $id_feature_type_cosmic == $Sanger::Cosmic::Dias::Constants::ID_FEATURE_TYPE_CODING ? $cds->{STOP} : undef,		# populate only for protein-coding transcripts
			CDNA_STOP 				=> $id_feature_type_cosmic == $Sanger::Cosmic::Dias::Constants::ID_FEATURE_TYPE_NONCODING ? $cds->{STOP} : undef,	# populate only for non-coding transcripts
			CDS_STOPOFFSET 			=> $cds->{STOPOFFSET},
			UTR_START 		 		=> $cds->{UTR_START},
			UTR_STOP 	 			=> $cds->{UTR_STOP},
			CDS_WT 					=> $cds->{WT},
			CDS_MT 					=> $cds->{MT},
			CDS_SYNTAX 				=> $cds->{SYNTAX},
			AA_START 				=> $protein->{START},
			AA_STOP 				=> $protein->{STOP},
			AA_WT					=> $protein->{WT},
			AA_MT 					=> $aa_mut,
			AA_SYNTAX 				=> $protein->{SYNTAX},
			GENOME_START			=> $genomic->{START},
			GENOME_STOP	 			=> $genomic->{STOP},
			GENOME_WT 				=> $genomic->{WT},
			GENOME_MT				=> $genomic->{MT},
			GENOME_SYNTAX 			=> $line_hash->{HGVSg},
			HGVSg_OFFSET 			=> $genomic->{HGVS_OFFSET},
			GENOME_VER 				=> $self->{config}->{assembly},
			PERCENT_MUT_ALLELE 		=> $input_var->percent_mut_allele,
			CHR 					=> $tva->variation_feature->seq_region_name,
			STRAND 					=> $cds->{STRAND},
			VARIANT_ONTOLOGY		=> $cds->{VARIANT_ONTOLOGY},
			CONSEQUENCES_ONTOLOGY 	=> $cds->{CONSEQUENCES_ONTOLOGY},
			ID_SAMPLE 				=> $input_var->id_sample,
			ID_STUDY 				=> $input_var->id_study,
			ID_PAPER 				=> $input_var->id_paper,
			GENE_NAME 				=> $cds->{GENE},
			ACCESSION 				=> $cds->{ACCESSION},
			CCDS 					=> $cds->{CCDS},
			DB 						=> $cds->{DB},
			DBVERSION 				=> $cds->{DBVERSION},
			ID_FEATURE_TYPE_COSMIC  => $id_feature_type_cosmic,
			ID_MUT_SOMATIC_STATUS 	=> $input_var->confirmed,
			ID_MUT_VERIF_STATUS 	=> $input_var->verified,
			#ID_MUTATION_CURRENT 	=> join(';', @{$input_var->id_mutation_current}),
			ID_MUTATION_CURRENT 	=> $input_var->id_mutation_current,
			ID_MUT_TYPE 			=> $cds->{VARIANT_ONTOLOGY},
			PARENT_MUT_LENGTH 		=> defined $cds->{MT} ? length($cds->{MT}) : '',
			CDS_MUT_LENGTH 			=> defined $cds->{MT} ? length($cds->{MT}) : '',
			AA_MUT_LENGTH 			=> defined $aa_mut ? length($aa_mut) : '',
			ANNOTATOR_VERSION 		=> $self->get_annotator_version,
		);
		my $data = $merge_hash->merge(\%default_data, \%transcript_data);
		return $data;
	}
}
#--------------------------------------------------------------------------------#
sub get_id_feature_type {
	my $cds = shift;
	if (is_coding($cds->{SYNTAX})) {
		return $Sanger::Cosmic::Dias::Constants::ID_FEATURE_TYPE_CODING;
	}
	else {
		if (is_within_gene_boundary($cds)) {
			return $Sanger::Cosmic::Dias::Constants::ID_FEATURE_TYPE_INTERGENIC;
		} else {
			return $Sanger::Cosmic::Dias::Constants::ID_FEATURE_TYPE_NONCODING;
		}
	}
}
#--------------------------------------------------------------------------------#
sub is_coding {
	my $syntax = shift;
	if (defined $syntax) {
		if ($syntax =~ /^c\./) {
			return 1;
		}
	}
	return 0;
}
#--------------------------------------------------------------------------------#
#These variants lie outside the cDNA (i.e. gene footprint) but are still within the gene boundaries. e.g. promoter regions
sub is_within_gene_boundary {
	my $cds = shift;
	if (/upstream_gene_variant/ ~~ $cds->{CONSEQUENCES_ONTOLOGY} || /downstream_gene_variant/ ~~ $cds->{CONSEQUENCES_ONTOLOGY}) {
		return 1;
	} else {
		return 0;
	}
}
#--------------------------------------------------------------------------------#
#TODO - Use Text::CSV_XS to parse this
sub get_input_variant_data {
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
															id_mutation_current    => $cols[10],
															#id_mutation_current    => [split(';', $cols[10])],
															id_study 			=> $cols[11] || undef,
															id_paper 			=> $cols[12] || undef,
															);
	return $var;
}
#--------------------------------------------------------------------------------#
sub get_genomic_data {
	my $self = shift;
	my $vfoa = shift;
	#my $genome_formatter = Sanger::Cosmic::Dias::GenomeFormatter->new(variation_allele => $vfoa, _slice_adaptor => $self->{config}->{sa});
	my $genome_formatter = Sanger::Cosmic::Dias::GenomeFormatter->new(variation_allele => $vfoa);
	return {WT 						=> $genome_formatter->wt_allele,
			MT 						=> $genome_formatter->mut_allele,
			START 					=> $genome_formatter->start,
			STOP 					=> $genome_formatter->stop,
			CHR 					=> $genome_formatter->chromosome,
			STRAND 					=> $genome_formatter->strand,
			DB 						=> $genome_formatter->db,
			DBVERSION 				=> $genome_formatter->db_version,
			VARIANT_ONTOLOGY 		=> $genome_formatter->variant_ontology,
			CONSEQUENCES_ONTOLOGY 	=> $genome_formatter->consequences_ontology,
			HGVS_OFFSET				=> $genome_formatter->hgvs_offset,
			};
}
#--------------------------------------------------------------------------------#
#VEP does not provide mutant aa for frameshift variants, so use Downstream column returned from the 'Downstream' plugin
sub get_aa_mut_allele {
	my ($protein, $line_hash) = @_;
	my $allele = $protein->{MT};
	if (/frameshift_variant/ ~~ $protein->{CONSEQUENCES_ONTOLOGY}) {
		$allele = $line_hash->{DownstreamProtein}."*";
	}
	return $allele || undef;
}
#--------------------------------------------------------------------------------#
sub get_annotator_version {
	my $self = shift;
	my $version = get_version_data()->{'ensembl-vep'};
	return $version->{release}.".".$version->{sub};
}
#--------------------------------------------------------------------------------#
1;
