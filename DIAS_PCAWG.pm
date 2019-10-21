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

package DIAS_PCAWG;

use strict;
use warnings;
use Sanger::Cosmic::Dias::VEPAnnotationFormatter;
use Sanger::Cosmic::Dias::GenomeFormatter;
#use Sanger::Cosmic::Dias::GenomicVariantDIAS;
use Bio::EnsEMBL::VEP::Utils qw(get_version_data);
use Sanger::Cosmic::Dias::Constants qw(%HEADER_DESCRIPTIONS_PCAWG %DEFAULT_COLUMN_VALUES_PCAWG);
use Try::Tiny;
use File::Basename;
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
	return '93.0';
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
	return \%Sanger::Cosmic::Dias::Constants::HEADER_DESCRIPTIONS_PCAWG;
}
#--------------------------------------------------------------------------------#
sub run {
	my ($self, $vfoa, $line_hash) = @_;
	#my ($self, $tva, $line_hash) = @_;
	my $assembly_ver = $self->config->{assembly};
	my $input_file = $self->config->{input_file};
	my $input_format = $self->config->{format};
	
	my $input_var;
	my $success = try {
		$input_var = parse_input_variant($line_hash->{Uploaded_variation}, $input_file, $input_format, $vfoa);
		1;
	} catch {
		warn "WARNING : $_\n";
		warn $line_hash->{Uploaded_variation}."\n";
		return 0;
	};
	return undef unless $success; 		# an undef returned from plugin->run() will filter that line from the output
	
	my $genomic = $self->get_genomic_data($vfoa);
	
	#my %default_data = %Sanger::Cosmic::Dias::Constants::DEFAULT_COLUMN_VALUES_PCAWG;
	#lock_keys(%default_data);
	#my $merge_hash = Hash::Merge->new('RIGHT_PRECEDENT');
	
	my %data;
	
	
	if ($vfoa->isa('Bio::EnsEMBL::Variation::IntergenicVariationAllele')) {
		%data = (
			GENOME_START			=> $genomic->{START},
			GENOME_STOP	 			=> $genomic->{STOP},
			GENOME_START_ORIGINAL	=> $genomic->{START_ORIGINAL},
			GENOME_STOP_ORIGINAL	=> $genomic->{STOP_ORIGINAL},
			GENOME_WT 				=> $genomic->{WT},
			GENOME_WT_ORIGINAL		=> $genomic->{WT_ORIGINAL},
			GENOME_MT				=> $genomic->{MT},
			GENOME_MT_ORIGINAL		=> $genomic->{MT_ORIGINAL},
			GENOME_STRAND_ORIGINAL	=> $genomic->{STRAND},
			GENOME_SYNTAX 			=> $line_hash->{HGVSg},
			HGVSg_OFFSET			=> $genomic->{HGVS_OFFSET},
			GENOME_VER 				=> $assembly_ver,
			CHR 					=> $genomic->{CHR},
			STRAND 					=> $genomic->{STRAND},	#VEP produces a 'STRAND' column, so no need for this
			ID_SAMPLE				=> $input_var->{id_sample},
			USERNAME 				=> $ENV{USER},
			DB 						=> $genomic->{DB},
			DBVERSION 				=> $genomic->{DBVERSION},
			ID_FEATURE_TYPE_COSMIC  => $Sanger::Cosmic::Dias::Constants::ID_FEATURE_TYPE_INTERGENIC,
			ID_MUT_TYPE 			=> $genomic->{VARIANT_ONTOLOGY},
			PARENT_MUT_LENGTH 		=> defined $genomic->{MT} ? length($genomic->{MT}) : '',
			VARIANT_ONTOLOGY		=> $genomic->{VARIANT_ONTOLOGY},
			CONSEQUENCES_ONTOLOGY 	=> $genomic->{CONSEQUENCES_ONTOLOGY},
			ANNOTATOR_VERSION 		=> $self->get_annotator_version,
			THOUSAND_GENOMES_AF		=> $input_var->{'1000GENOMES_AF'},
			THOUSAND_GENOMES_ID		=> $input_var->{'1000GENOMES_ID'},
			CALLERS					=> $input_var->{CALLERS},
			PCAWG_COSMIC					=> $input_var->{COSMIC},
			DBSNP					=> $input_var->{DBSNP},
			DBSNP_SOMATIC			=> $input_var->{DBSNP_SOMATIC},
			MODEL_SCORE				=> $input_var->{MODEL_SCORE},
			N_ALT_COUNT				=> $input_var->{N_ALT_COUNT},
			N_REF_COUNT				=> $input_var->{N_REF_COUNT},
			NUMCALLERS				=> $input_var->{NUMCALLERS},
			N_VAF					=> $input_var->{N_VAF},
			REPEAT_MASKER			=> $input_var->{REPEAT_MASKER},
			SIGNATURE_N3			=> $input_var->{SIGNATURE_N3},
			SIGNATURE_R1			=> $input_var->{SIGNATURE_R1},
			SIGNATURE_R2			=> $input_var->{SIGNATURE_R2},
			SNV_NEAR_INDEL			=> $input_var->{SNV_NEAR_INDEL},
			T_ALT_COUNT				=> $input_var->{T_ALT_COUNT},
			T_REF_COUNT				=> $input_var->{T_REF_COUNT},
			VAF						=> $input_var->{VAF},
			VARIANT_CLASSIFICATION	=> $input_var->{VARIANT_CLASSIFICATION}
		);
		
		#my $data = $merge_hash->merge(\%default_data, \%intergenic_data);
		#return $data;
		
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
	
		%data = (
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
			GENOME_START_ORIGINAL	=> $genomic->{START_ORIGINAL},
			GENOME_STOP_ORIGINAL	=> $genomic->{STOP_ORIGINAL},
			GENOME_WT 				=> $genomic->{WT},
			GENOME_WT_ORIGINAL		=> $genomic->{WT_ORIGINAL},
			GENOME_MT				=> $genomic->{MT},
			GENOME_MT_ORIGINAL		=> $genomic->{MT_ORIGINAL},
			GENOME_STRAND_ORIGINAL 	=> $genomic->{STRAND},
			GENOME_SYNTAX 			=> $line_hash->{HGVSg},
			HGVSg_OFFSET 			=> $genomic->{HGVS_OFFSET},
			GENOME_VER 				=> $assembly_ver,
			CHR 					=> $tva->variation_feature->seq_region_name,
			STRAND 					=> $cds->{STRAND},
			VARIANT_ONTOLOGY		=> $cds->{VARIANT_ONTOLOGY},
			CONSEQUENCES_ONTOLOGY 	=> $cds->{CONSEQUENCES_ONTOLOGY},
			ID_SAMPLE				=> $input_var->{id_sample},
			GENE_NAME 				=> $cds->{GENE},
			ACCESSION 				=> $cds->{ACCESSION},
			CCDS 					=> $cds->{CCDS},
			USERNAME 				=> $ENV{USER},
			DB 						=> $cds->{DB},
			DBVERSION 				=> $cds->{DBVERSION},
			ID_FEATURE_TYPE_COSMIC  => $id_feature_type_cosmic,
			ID_MUT_TYPE 			=> $cds->{VARIANT_ONTOLOGY},
			PARENT_MUT_LENGTH 		=> defined $cds->{MT} ? length($cds->{MT}) : '',
			CDS_MUT_LENGTH 			=> defined $cds->{MT} ? length($cds->{MT}) : '',
			AA_MUT_LENGTH 			=> defined $aa_mut ? length($aa_mut) : '',
			ANNOTATOR_VERSION 		=> $self->get_annotator_version,
			THOUSAND_GENOMES_AF		=> $input_var->{'1000GENOMES_AF'},
			THOUSAND_GENOMES_ID		=> $input_var->{'1000GENOMES_ID'},
			CALLERS					=> $input_var->{CALLERS},
			PCAWG_COSMIC					=> $input_var->{COSMIC},
			DBSNP					=> $input_var->{DBSNP},
			DBSNP_SOMATIC			=> $input_var->{DBSNP_SOMATIC},
			MODEL_SCORE				=> $input_var->{MODEL_SCORE},
			N_ALT_COUNT				=> $input_var->{N_ALT_COUNT},
			N_REF_COUNT				=> $input_var->{N_REF_COUNT},
			NUMCALLERS				=> $input_var->{NUMCALLERS},
			N_VAF					=> $input_var->{N_VAF},
			REPEAT_MASKER			=> $input_var->{REPEAT_MASKER},
			SIGNATURE_N3			=> $input_var->{SIGNATURE_N3},
			SIGNATURE_R1			=> $input_var->{SIGNATURE_R1},
			SIGNATURE_R2			=> $input_var->{SIGNATURE_R2},
			SNV_NEAR_INDEL			=> $input_var->{SNV_NEAR_INDEL},
			T_ALT_COUNT				=> $input_var->{T_ALT_COUNT},
			T_REF_COUNT				=> $input_var->{T_REF_COUNT},
			VAF						=> $input_var->{VAF},
			VARIANT_CLASSIFICATION	=> $input_var->{VARIANT_CLASSIFICATION},
		);
		#my $data = $merge_hash->merge(\%default_data, \%transcript_data);
		#return $data;
	}
	return \%data;
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
sub parse_input_variant {
	my ($input, $input_file, $input_format, $vfoa) = @_;	
	
	die "ERROR: Cannot parse $input_format. Expecting VCF\n" unless $input_format eq 'vcf';
	
	my ($in_filename, $in_dir, $suffix) = fileparse($input_file);
	my $icgc_UUID = (split(/\./, $in_filename))[0];		# This is the tumour aliquot ID that will be mapped to the ICGC sample name using COSMIC-ICGC mapping tables
	
	my $vcf_input = $vfoa->variation_feature->{_line};	# Get the input VCF line
	my $vcf_info = $vcf_input->[7];						# We want to keep the key-value data in the INFO column of the VCF
	
	my %vcf_info_hash = ();
	for my $key_value (split(';', $vcf_info)) {
		my ($key, $value) = split('=', $key_value);
		$vcf_info_hash{uc($key)} = $value;
	}
	
	my @vals = split('_', $input);
	my ($wt, $mut) = split('/', $vals[2]);
	my %var = (chr => $vals[0],
			   start => $vals[1],
			   stop => $vals[1],
			   wt => $wt,
			   mut => $mut,
			   id_sample => $icgc_UUID,
			   );
	
	@var{keys %vcf_info_hash} = values %vcf_info_hash;	#Add %vcf_info_hash into %var
	return \%var;
}
#--------------------------------------------------------------------------------#
sub get_genomic_data {
	my $self = shift;
	my $vfoa = shift;
	#my $genome_formatter = Sanger::Cosmic::Dias::GenomeFormatter->new(variation_allele => $vfoa, _slice_adaptor => $self->{config}->{sa});
	my $genome_formatter = Sanger::Cosmic::Dias::GenomeFormatter->new(variation_allele => $vfoa);
	return {WT 						=> $genome_formatter->wt_allele,
			WT_ORIGINAL				=> $genome_formatter->wt_allele_original,
			MT 						=> $genome_formatter->mut_allele,
			MT_ORIGINAL				=> $genome_formatter->mut_allele_original,
			START 					=> $genome_formatter->start,
			START_ORIGINAL			=> $genome_formatter->start_original,
			STOP 					=> $genome_formatter->stop,
			STOP_ORIGINAL			=> $genome_formatter->stop_original,
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
	if (/frameshift_variant/ ~~ $protein->{CONSEQUENCES_ONTOLOGY}) {	# Only want to annotate frameshifts which have a known protein mutation
		if ($protein->{SYNTAX} ne 'p.?' && $line_hash->{HGVSp} !~ /p.Met1.*?\?/) {	# Skip frameshifts with unknown protein mutation and start_lost consequences that look like p.Met1? (the latter don't have a known syntax, therefore we cannot provide an AA_MUT allele)
			$allele = $line_hash->{DownstreamProtein};
			warn "mut_allele is empty for frameshift variant - ".$line_hash->{Uploaded_variation}."\n" unless $allele;
			$allele .= "*";	# Add the stop codon
		}
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
