package InsertionSeq::ContigBlastHit;

#############################################################################
#
#    ISseeker    - Finds portions of contigs flanking IS sequences
#                  and blasts them against annotated references 
#                  to infer IS insertion points
#
#    Written by Brian Bishop
# 
#
#    Copyright (C) 2015  JCVI
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#    See example in doc directory
#
#
#############################################################################

use parent q(InsertionSeq::BlastHit);

use strict;

use Log::Log4perl;

use InsertionSeq::Flank;


##
## Whole contig is IS
##
our $TYPE_ENTIRE		= "ENTIRE";

##
## IS in the middle with at least ~50 bases extra on each side
##
our $TYPE_EMBED			= "EMBED";

##
## Embedded, but doesn't go completely to IS ends
##
our $TYPE_EMBED_TRUNC	= "EMBED_TRUNC";

##
## Beginning Edge
##
our $TYPE_EDGE_BEGIN	= "EDGE_BEGIN";
##
## Ending Edge
##
our $TYPE_EDGE_END		= "EDGE_END";

##
## An Edge, but IS seq too long (gaps?)
##
our $TYPE_EDGE_LONG		= "EDGE_LONG";

##
## Edge but doesn't go completely to IS end
##
our $TYPE_EDGE_BEGIN_TRUNC	= "EDGE_BEGIN_TRUNC";
our $TYPE_EDGE_END_TRUNC	= "EDGE_END_TRUNC";

##
## Flank type codes:
## Whole match is entire flank +/- $ANNOT_BASE_WIGGLE_ROOM
##
our $TYPE_FLANK_WHOLE         = "W";
our $TYPE_FLANK_PARTIAL       = "P";
our $TYPE_FLANK_NONE          = "X";

##
## Don't know what this is
our $TYPE_UNKN			= "UNKN";


##
## What direction is the IS with respect to the contig
##
our $DIRECTION_FWD		= "F";
our $DIRECTION_REV		= "R";

##
## Distance between flank match positions must be within this percentage of IS coordinates on the reference
##
our $WIGGLE_OVERLAP_PERCENT     = 90;

our $CONTIG_WIGGLE_ROOM		= 10;
##
## Match must be this close to one (edge) or both (embedded) IS ends or is rejected
##
our $IS_WIGGLE_ROOM			= 10;
##
## Match must be this far away from Contig edge to be considered embedded
##
our $BASE_EDGE_MARGIN		= 50;

##
## Bases to extract on each side, where possible
##
our $FLANK_EXTRACT_BASES			= 500;
##
## Less than this, we don't have an interesting flank
##
our $MINIMUM_FLANK_EXTRACT_BASES	= 20;

##
## Tokens to indicate Beginning/Ending Flanks
##
our $BEGINNING_FLANK			= "B";
our $ENDING_FLANK				= "E";




my $log = Log::Log4perl->get_logger();


sub new 
{
    my $class = shift;
    my $self = $class->SUPER::new(@_);
    
    
    #die ("ContigBlastHit missing IS name.\n") if ( !defined($self->{is_name}) );
    die ("ContigBlastHit missing IS name.\n") unless ( $self->get_is_name() );
 
 	$self->evaluate();
    
    return $self;
}





sub type {
	my $self = shift;
	if (@_) {
		$self->{type} = shift;
	}
	return $self->{type};
}

sub evaluate {
	my $self = shift;
	
	#
	#   identify scaffolds comprised entirely of IS element
	#
	if ( $self->{slength} - $self->{shitlength} <= $CONTIG_WIGGLE_ROOM )
	{
		$self->{type} = $TYPE_ENTIRE;
	}
	## s_start & s_end are >50 bases from the scaffold edges
    ## q_start and q_end are 1 and IS_length, respectively (adding IS wiggle room)
	elsif (    $self->{sstart}   >  $BASE_EDGE_MARGIN 
			&& $self->{send}     <  $self->{slength} - $BASE_EDGE_MARGIN)
	{
			if ($self->is_full_length_is())
			{
				$self->{type} = $TYPE_EMBED;
			}
			else
			{
				$self->{type} = $TYPE_EMBED_TRUNC;
			}
	}
	## 2) identify scaffolds with IS at edge
    ##    s_start or s_end = 1 or scaffold_length (from info file created above)
    ##    match length usually 17-50 bases but < IS_length
    elsif ($self->{sstart} <= $BASE_EDGE_MARGIN 
        || $self->{send} >=  $self->{slength} - $BASE_EDGE_MARGIN)   
	{
		if ($self->{qhitlength} <= $self->{qlength})
		{
				if ($self->is_partial_is() || $self->is_full_length_is())
				{
					$self->{type} = $TYPE_EDGE_BEGIN if ($self->{sstart} <= $BASE_EDGE_MARGIN);
					$self->{type} = $TYPE_EDGE_END if ($self->{send} >=  $self->{slength} - $BASE_EDGE_MARGIN);
				}
				else
				{
					$self->{type} = $TYPE_EDGE_BEGIN_TRUNC if ($self->{sstart} <= $BASE_EDGE_MARGIN);
					$self->{type} = $TYPE_EDGE_END_TRUNC if ($self->{send} >=  $self->{slength} - $BASE_EDGE_MARGIN);
				}

		}
		else
		{
			## Query Hit length is greater than Query Sequence Length. 
			## Caused by too many gaps?
			## Probably shouldn't happen with our criteria:
			$self->{type} = $TYPE_EDGE_LONG;
		}
	}
	else
	{
		$self->{type} = $TYPE_UNKN;
	}

}

sub get_flank
{
	my $self = shift;
	my $location = shift;
	my ($left,$right) = $self->get_flank_coords($location);
	
	my $flank;
	
	if (defined($left))
	{
#		$log->info("Flank returned.\n");
		$flank =  InsertionSeq::Flank->new(
			is_name => $self->{is_name},
			location => $location,
			contig_direction => $self->{sdir},
			contig_lower_coord => $left,
			contig_upper_coord => $right,
			contig_name => $self->{name},
			is_pct_id => $self->{pct_id},
			contig_file_name => $self->{filename});
	}
#	else
#	{
#		$log->info("No flank returned.\n");
#	}
	
	return $flank;
}

sub get_flank_coords
{
	my $self = shift;
	my $direction = shift;

	if ($direction eq $BEGINNING_FLANK)
	{
    	return $self->get_upstream_flank_coords();
	}
	elsif ($direction eq $ENDING_FLANK)
	{
    	return $self->get_downstream_flank_coords();
	}
	else
	{
		die "Unknown flank coord extraction direction: $direction\n";	
	}
}

sub get_upstream_flank_coords
{
	my $self = shift;

	my ($start,$end);

	##
	## Rev match at start edge
	## or rev middle match 
	if ($self->{sdir} eq $DIRECTION_REV
	    && ($self->{type} eq $TYPE_EDGE_BEGIN || $self->{type} eq $TYPE_EMBED))
	{
      	$start = $self->{send} + 1;
		$end = ($start + $FLANK_EXTRACT_BASES) - 1;
		## Don't go past end of subject seq
		$end = $self->{slength} if ($end > $self->{slength});
	} 
	##
	## Fwd match at end edge
	## or fwd middle match
	elsif ($self->{sdir} eq $DIRECTION_FWD
	    && ($self->{type} eq $TYPE_EDGE_END || $self->{type} eq $TYPE_EMBED))
	{
      	$end = $self->{sstart} - 1;
    	$start = ($end - $FLANK_EXTRACT_BASES) + 1; 
		$start = 1 if ($start < 1);
	}

	##
	## If too small, forget about it
	##
	if ( ($end - $start) + 1 < $MINIMUM_FLANK_EXTRACT_BASES)
	{
    	undef $start;
		undef $end;
	}

	($start,$end);
}

sub get_downstream_flank_coords
{
	my $self = shift;

	my ($start,$end);

	##
	## Fwd match at start edge
	## or middle match 
	if ($self->{sdir} eq $DIRECTION_FWD 
		&& ($self->{type} eq $TYPE_EDGE_BEGIN || $self->{type} eq $TYPE_EMBED))
	{
      	$start = $self->{send} + 1;
		$end = ($start + $FLANK_EXTRACT_BASES) - 1;
		## Don't go past end of subject seq
		$end = $self->{slength} if ($end > $self->{slength});
	} 
	##
	## Rev match at end edge
	## or rev middle match
	elsif ($self->{sdir} eq $DIRECTION_REV
	        && ($self->{type} eq $TYPE_EDGE_END || $self->{type} eq $TYPE_EMBED))
	{
      	$end = $self->{sstart} - 1;
    	$start = ($end - $FLANK_EXTRACT_BASES) + 1; 
		$start = 1 if ($start < 1);
	}

	##
	## If too small, forget about it
	##
	if ( ($end - $start) + 1 < $MINIMUM_FLANK_EXTRACT_BASES)
	{
    	undef $start;
		undef $end;
	}


	($start,$end);
}

sub is_annotatable
{
	my $self = shift;
	
	return 1 if ( ( $self->{type} eq  $TYPE_EDGE_BEGIN 
				|| $self->{type} eq  $TYPE_EDGE_END 
				|| $self->{type} eq $TYPE_EMBED) );
	
	return 0;
}


sub is_full_length_is
{
	my $self = shift;
	return 1 if ($self->{qstart}  <=  $IS_WIGGLE_ROOM + 1
				&& $self->{qend}    >=  $self->{qlength} - $IS_WIGGLE_ROOM);
		
	return 0;
}

sub is_partial_is
{
	my $self = shift;
	return 1 if ($self->{qstart}  <=  $IS_WIGGLE_ROOM + 1
				||  $self->{qend}    >=  $self->{qlength} - $IS_WIGGLE_ROOM);
		
	return 0;
}



sub to_feat_mysql
{
    my $self = shift;
    my $q_genome = shift;
    my $s_genome = shift;
    my $reference = shift;


	my $sql = "";
	if ($self->{type} eq $TYPE_ENTIRE || $self->{type} eq $TYPE_EMBED || $self->{type} eq $TYPE_EDGE_BEGIN || $self->{type} eq $TYPE_EDGE_END || $self->{type} eq $TYPE_EMBED_TRUNC)
	#if ($self->{type} eq $TYPE_ENTIRE  || $self->{type} eq $TYPE_EMBED)
	{
    	$sql = "INSERT INTO is_query_feature (is_run_id,is_element, reference, q_genome, s_genome, is_pct_id, flank_pct_id, contig_name, feat_type, flank_begin_end, orientation, begin_base, end_base, is_annotated ) VALUES (";
        $sql .= '@is_run_id,';
        $sql .= "'$self->{is_name}',";
        $sql .= "'$reference',";
        $sql .= "'$q_genome',";
        $sql .= "'$s_genome',";
        $sql .= "$self->{pct_id},";
        $sql .= "NULL,";
		$sql .= "'$self->{name}',";
        $sql .= "'E',";
        $sql .= "'X',";
        $sql .= "'$self->{sdir}',";
        $sql .= "'$self->{sstart}',";
        $sql .= "'$self->{send}',";
        $sql .= "0";
        $sql .= ");\n" ;
        #$sql .= "SET $SQL_LAST_FEAT_ID_VAR = LAST_INSERT_ID();\n";
	}
	else
	{
    	$log->error("Illegal Contig Blast Hit type sent to routine to_feat_mysql() for contig $self->{name}: $self->{type}\n");
	}

    $sql;


}


##our $CSV_HEADER = "id, mate_id, offset_from_previous, is_element, genome, contig_name, %_id, match_len, feat_type, flank_begin_end, orientation, contig_flank_begin_base, contig_flank_end_base, is_annotated, reference, nearest_base, after_gene, in_gene, before_gene, match_quality, contig_count\n";
sub to_csv
{
    my $self = shift;
    my $genome = shift;

	#$TYPE_EDGE_BEGIN $TYPE_EDGE_END
	my $csv;
	if ($self->{type} eq $TYPE_ENTIRE || $self->{type} eq $TYPE_EMBED || $self->{type} eq $TYPE_EDGE_BEGIN || $self->{type} eq $TYPE_EDGE_END || $self->{type} eq $TYPE_EMBED_TRUNC)
	{
		$csv .= "$self->{id},";
		$csv .= ","; 			# No mate
		$csv .= ","; 			# No offset
        $csv .= "$self->{is_name},";
        $csv .= "$genome,";
        $csv .= ",";
		$csv .= "$self->{pct_id},";
		$csv .= "$self->{matchlen},";
        $csv .= "$self->{type},";
        $csv .= "X,";
        $csv .= "$self->{sdir},";
        $csv .= "$self->{sstart},";
        $csv .= "$self->{send},";
        $csv .= "0";
        $csv .= "\n" ;

	}
	else
	{
    	$log->logdie("Illegal Contig Blast Hit type sent to routine to_csv() for contig $self->{name}: $self->{type}\n");
	}

    $csv;


}


## Check to see if IS flanking sequences correspond to a known IS element in the reference genome(s)
sub check_reference_for_IS
{
    my $type = shift;
    my $name = shift;
    my $start = shift;
    $start++; # need to correct for flank coordinates versus IS coordinates
    my $end = shift;
    $end--; # need to correct for flank coordinates versus IS coordinates
    my $length = ($end - $start) + 1;
    my $referenceISHits = shift;

    #print STDERR "check_reference_for_IS: $name $start $end $length\n";
    for my $ref_hit (@{$referenceISHits}) {
	#print STDERR "refIS: $ref_hit->{name} $ref_hit->{sstart} $ref_hit->{send} $ref_hit->{shitlength}\n";
	if ($name ne $ref_hit->{name}) {
	    #skip hits to different contigs
	    next;
	}
	#check for overlap
	if (($end < $ref_hit->{sstart}) || ($ref_hit->{send} < $start)) {
	    #no overlap so skip
	    next;
	}
	my $start_overlap = ($start > $ref_hit->{sstart}) ? $start : $ref_hit->{sstart};
	my $end_overlap = ($end < $ref_hit->{send}) ? $end : $ref_hit->{send};
	my $overlap_length_100 = (($end_overlap - $start_overlap) + 1) * 100;
	if ((($overlap_length_100 / $length) >= $WIGGLE_OVERLAP_PERCENT) && (($overlap_length_100 / $ref_hit->{shitlength}) >= $WIGGLE_OVERLAP_PERCENT)) {
	    #print STDERR "return 1\n";
	    return 1;
	}
    }

    return 0;
}


1;
