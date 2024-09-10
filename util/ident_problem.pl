#!/usr/bin/perl
use warnings;
use strict;
use File::Basename;
use Class::Struct;

struct (PRIMERINFO => {id => '$', count => '$'});
struct (CHRBED => {id => '$', fwd_regions => '@', rev_regions => '@'});
struct (SEQREGION => {pos => '$', primer_idx => '$'});
struct (BLASTINFO => {primers => '@', chrs => '@'});

###Declarations

my $MAX_REGION_LEN = 1000;
my $MAX_HIT_PERCENTILE = 80;

sub ParseBlastRes($);		#my $blastInfoObj = ParseBlastRes($file);
sub FindPrimerPairs($); 	#my %pairs = FindPrimerPairs($blastInfo);
sub BlacklistByPair($$); 	#my %blacklistSet = BlacklistByPair(\%pairDict,$blastInfo);
sub BlacklistByPercentile($);	#my %blacklistSet = BlacklistByPercentile($blastInfo);

#Test routines:
sub LogBlastInfoObj();		#LogBlastInfoObj();
sub LogPrimerPairs();		#LogPrimerPairs();

###Handling Input

if(@ARGV < 1){
    die "Usage: ".basename($0). " blastres.tab > primers.list\n";
}

##MAIN Script

my $BlastInfoObj = ParseBlastRes(shift(@ARGV));
#LogBlastInfoObj();
my %PrimerPairs = FindPrimerPairs($BlastInfoObj);
#LogPrimerPairs();
my %BlacklistSet = BlacklistByPair(\%PrimerPairs,$BlastInfoObj);
%BlacklistSet = (%BlacklistSet,BlacklistByPercentile($BlastInfoObj));
foreach my $idx (sort keys %BlacklistSet){
    print $BlastInfoObj->primers->[$idx]->id, "\n";
}

##Subroutien Definitions

#Given a path to a blast formated file in (default outfmt 6) parses a datastructure to be
#   used to identify problem primers
#Assumes that the blast results are sorted by the primer id
#Input - path to blastres file
#Output - BLASTINFO Object
sub ParseBlastRes($){
    my $file = shift;
    open(my $fh, $file) or die "Could not open $file: $!\n";
    my @primerList;
    my @chrBEDList;
    my %chrIdx;
    while(my $line = <$fh>){
        chomp($line);
        my ($pID,$sID,undef,undef,undef,undef,undef,undef,$s,$e,undef,undef) = split(/\t/,$line);
        my $strand = ($s > $e) ? '-' : '+';
        #Get the information on the primer
        if(@primerList and $pID eq $primerList[$#primerList]->id){
            $primerList[$#primerList]->count($primerList[$#primerList]->count + 1);
        } else {
            push(@primerList,PRIMERINFO->new(id => $pID, count => 1));
        }
        #Assemble the information on the sequence
        my $sReg = SEQREGION->new(pos => $s, primer_idx => $#primerList);
        #Find the CHRBED object, or create it if necessary
        my $idx = $#chrBEDList + 1;
        if(exists $chrIdx{$sID}){
            $idx = $chrIdx{$sID};
        } else {
            push(@chrBEDList,CHRBED->new(id => $sID));
            $chrIdx{$sID} = $idx;
        }
        #Add the sequence region to the appropriate set of regions
        if($strand eq '+'){
            push(@{$chrBEDList[$idx]->fwd_regions},$sReg);
        } else{
            push(@{$chrBEDList[$idx]->rev_regions},$sReg);
        }
    }
    close($fh);
    #Ensure each set of regions is sorted in order by position to save comparisons later
    foreach my $chrBed (@chrBEDList){
        @{$chrBed->fwd_regions} = sort {$a->pos <=> $b->pos} @{$chrBed->fwd_regions};
        @{$chrBed->rev_regions} = sort {$a->pos <=> $b->pos} @{$chrBed->rev_regions};
    }
    my $obj = BLASTINFO->new(primers => \@primerList,chrs => \@chrBEDList);
    return $obj;
}

#Given the blast information structure finds all pairs of primers with small enough
#amplicons
#Input - a BLASTINFO obj
#Output - a primer_idx keyed hash of sets of primer_idx
sub FindPrimerPairs($){
    my $blastInfo = shift;
    my %pairDict;
    foreach my $chrBed (@{$blastInfo->chrs}){
        #compare all fwd regions to rev regions creating pairs whenever they overlap
        #Note: fwd and revRegs are sorted by ascending position so each subsequent fwd region
        #needn't be compared to every rev region, only thosw which could overlap
        my $first = 0;
        foreach my $fwdReg (@{$chrBed->fwd_regions}){
            for(my $i = $first; $i < @{$chrBed->rev_regions}; $i++){
                my $revReg = $chrBed->rev_regions($i);
                #Skip if to far apart, and note for future comparisons
                if($revReg->pos - $MAX_REGION_LEN  > $fwdReg->pos){
                    $first++;
                    next;
                }
                #Stop looking if fwd is after reverse
                last if($fwdReg->pos > $revReg->pos);
                #If this point is reached the regions overlap
                $pairDict{$fwdReg->primer_idx} = {} unless(exists $pairDict{$fwdReg->primer_idx});
                $pairDict{$fwdReg->primer_idx}->{$revReg->primer_idx} = 1;
                $pairDict{$revReg->primer_idx} = {} unless(exists $pairDict{$revReg->primer_idx} );
                $pairDict{$revReg->primer_idx}->{$fwdReg->primer_idx} = 1;
            }
        }
    }
    return %pairDict;
}

#Given a hash of sets of primers each primer forms pairs with, return a minimal set of
#   primers which eliminate all pairs
#Proceeds by a greedy algorithm, at each iteration the primer which participates in the
#   most pairs is eliminated; ties broken by the number of hits for the primer; if there
#   are still ties, both are removed
#Input - reference to a primer_idx keyed hash of sets of primer_idx
#      - a BLASTINFO object
#Output - a set of primer_idx
sub BlacklistByPair($$){
    my $pairDictRef = shift;
    my $blastInfo = shift;
    my %blacklist;
    while(%{$pairDictRef}){
        my %worst = (val => -1, idx => [], count => -1);
        while(my ($idx,$setRef) = each %{$pairDictRef}){
            my $val = scalar(keys(%{$setRef}));
            my $count = $blastInfo->primers->[$idx]->count;
            if($val > $worst{val} or ($val == $worst{val} and $count > $worst{count})){
                $worst{val} = $val;
                $worst{idx} = [$idx];
                $worst{count} = $count;
            } elsif($val == $worst{val} and $count == $worst{count}){
                push(@{$worst{idx}},$idx);
            }
        }
        foreach my $idx (@{$worst{idx}}){
            foreach my $sub_idx (keys %{$pairDictRef->{$idx}}){
                delete $pairDictRef->{$sub_idx}->{$idx};
                if(!%{$pairDictRef->{$sub_idx}}){
                    delete $pairDictRef->{$sub_idx};
                }
            }
            delete $pairDictRef->{$idx};
            $blacklist{$idx} = 1;
        }
    }
    return %blacklist;
}

#Given a BLASTINFO object blacklists primers whith hit counts above the MAX_HIT_PERCENTILE
#Input - a BLASTINFO object
#Output - a set of primer_idx
sub BlacklistByPercentile($){
    my $blastInfo = shift;
    my %uniqCounts = map {($_->count,1)} @{$blastInfo->primers};
    my @counts = sort {$b <=> $a} (keys %uniqCounts);
    my $threshold = $counts[int(@counts * (100 - $MAX_HIT_PERCENTILE) / 100)];
    my %blacklist;
    for(my $i = 0; $i < @{$blastInfo->primers}; $i++){
        $blacklist{$i} = 1 if($blastInfo->primers->[$i]->count > $threshold);
    }
    return %blacklist;
}

sub LogBlastInfoObj(){
    foreach my $primer (@{$BlastInfoObj->primers}){
        print STDERR $primer->id,"\t",$primer->count,"\n";
    }
    foreach my $chr (@{$BlastInfoObj->chrs}){
        print STDERR $chr->id,"\n+:";
        foreach my $reg (@{$chr->fwd_regions}){
            print STDERR $BlastInfoObj->primers->[$reg->primer_idx]->id,'@',$reg->pos,"\t";
        }
        print STDERR "\n-:";
        foreach my $reg (@{$chr->rev_regions}){
            print STDERR $BlastInfoObj->primers->[$reg->primer_idx]->id,'@',$reg->pos,"\t";
        }
        print STDERR "\n";
    }
}

sub LogPrimerPairs(){
    foreach my $idx (sort keys %PrimerPairs){
        print STDERR $BlastInfoObj->primers->[$idx]->id,"\n\t";
        foreach my $sub_idx (sort keys %{$PrimerPairs{$idx}}){
            print STDERR $BlastInfoObj->primers->[$sub_idx]->id,"\t";
        }
        print STDERR "\n";
    }
}
