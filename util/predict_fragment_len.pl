#!/usr/bin/perl
use warnings;
use strict;
use Bio::SeqIO;
use File::Basename;
use Zach::Util::File qw(OpenFileHandle);
use Zach::Util::Array qw(Unique);

sub LoadPrimers($);
sub ProcessContext($$\%);
sub LocatePrimer($$$);
sub ProcessRxn($\%);

if(@ARGV < 3){
    die "Usage: ".basename($0). " SuggestedRxns PrimerSeq.fna TargetInfo ContextSeq.fna > AmpliconSizes.tab\n".
        "\tOutput has rows with the first element of RXNID and the subsequent elements are each an amplicon len\n";
}

my %FileDict;
@FileDict{qw(rxn primer targetinfo context)} = @ARGV[0 .. 3];
while (my ($type,$file) = each %FileDict){
    die "$type file ($file) does not exists\n" unless(-f $file);
}

sub main {
    my %primerDict = LoadPrimers($FileDict{primer});
    my %ampLenDict = ProcessContext($FileDict{targetinfo},$FileDict{context},%primerDict);
    ProcessRxn($FileDict{rxn},%ampLenDict);
} main();


sub LoadPrimers($){
    my $fh = OpenFileHandle(shift,"primer","ERROR");
    my $in = Bio::SeqIO->new(-fh => $fh, -format => 'fasta');
    my %primerDict; #primerID keyed hash of arrays with two elements (fwd and rev);
    while(my $seqObj = $in->next_seq){
        $primerDict{$seqObj->id} = [$seqObj->seq,$seqObj->revcom()->seq];
    }
    return %primerDict;
}

sub ProcessContext($$\%){
    my ($tgtInfoFile,$contextFile,$rPrimerDict) = @_;
    my $fh = OpenFileHandle($tgtInfoFile,"tgtInfo","ERROR");
    <$fh>;
    my %tgtInfo; #targetID keyed hash of arrays of two element arrays (fwd pID, rev pID)
    while(my $line = <$fh>){
        chomp($line);
        my ($tgt,undef,$fwdStr,$revStr) = split(/\t/,$line);
        my @fwdList = split(/,/,$fwdStr);
        my @revList = split(/,/,$revStr);
        $tgtInfo{$tgt} = [];
        foreach my $fPID (@fwdList){
            foreach my $rPID (@revList){
                push(@{$tgtInfo{$tgt}}, [$fPID,$rPID]);
            }
        }
    }
    close($fh);
    my %ampliconLenDict; #tgt/fwd/rev keyed hash of amplicon lengths
    $fh = OpenFileHandle($contextFile,"context","ERROR");
    my $in = Bio::SeqIO->new(-fh => $fh, -format => 'fasta');
    while(my $seqObj = $in->next_seq){
        next unless(exists $tgtInfo{$seqObj->id});
        foreach my $rPair (@{$tgtInfo{$seqObj->id}}){
            #Skip unless primer sequences are loaded for both the fwd and reverse primer
            #   Implemented as skip if there are any primer ids which don't exists
            next unless(!grep {!exists $rPrimerDict->{$_}} @{$rPair});
            #   rPrimerDict->{rPair->[0]}->[0] => using the fwd primer id, use the fwd sequence of that primer
            my @primerPos = map {LocatePrimer($rPrimerDict->{$rPair->[$_]}->[$_],$seqObj->seq,$_)} (0,1);
            #Skip if any positions are not found
            next if(grep {!defined $_} @primerPos);
            #Skip if the fwd pos is too close to revPos
            my $primerSum = eval(join("+",(map {length($rPrimerDict->{$rPair->[$_]}->[0])} (0,1))));
            my $ampLen = -1*eval(join("-",@primerPos))+1;
            next if($ampLen < $primerSum);
            my $key = join("/",($seqObj->id,@{$rPair}));
            $ampliconLenDict{$key} = $ampLen;
        }
    }
    return %ampliconLenDict;
}

#Given a primer sequence and its target, finds the outermost position of an exact match
#   in the appropriate half of the sequences
#If fwd, starts at 0 and slides right until the halfway point
#If rev, starts at the end and slides left until the halfway points
#If no match undef is returned
sub LocatePrimer($$$){
    my($primer,$context,$bRev) = @_;
    if(!$bRev){
        for(my $i = 0; $i < length($context)/2 + 1; $i++){
            return $i if($primer eq substr($context,$i,length($primer)));
        }
    } else {
        for(my $i = length($context)-1; $i > length($context)/2 - 1; $i--){
            return $i if($primer eq substr($context,$i-length($primer)+1,length($primer)));
        }
    }
    return undef;
}

sub ProcessRxn($\%){
    my $fh = OpenFileHandle(shift,"rxn","ERROR");
    my $rAmpDict = shift;
    <$fh>;
    print "Reaction\tAmplicon\tLen\n";
    while(my $line = <$fh>){
        chomp($line);
        my($rID,undef,$pIDStr,$tIDStr) = split(/\t/,$line);
        my @pIDList = Unique(split(/,/,$pIDStr));
        my @tIDList = Unique(split(/,/,$tIDStr));
        foreach my $tID (@tIDList){
            foreach my $fPID (@pIDList){
                foreach my $rPID (@pIDList){
                    my $key = join("/",($tID,$fPID,$rPID));
                    next unless(exists $rAmpDict->{$key});
                    print join("\t",($rID,$key,$rAmpDict->{$key})),"\n";
                }
            }
        }

    }
}
