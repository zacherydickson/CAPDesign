#!/bin/bash

if [ "$#" -lt 2 ]; then
    >&2 echo "Usage: $(basename $0) Targets.fna BlasklistBlastDb [Threads=1]> Blacklist.fna";
    exit 1;
fi

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
CAPDESIGN="$SCRIPT_DIR/../CAPDesign"

TargetFile=$1;
BlastDB=$2;
Threads=$3;
[ -z $Threads ] && Threads=1;
CandFile="/tmp/CAP_Black_$(rndstr 8).fna";
BlastFile="/tmp/CAP_Black_$(rndstr 8).tab";
BlackFile="/tmp/CAP_Black_$(rndstr 8).list";

"$CAPDESIGN" -C "$TargetFile" >| $CandFile;
blastn -task blastn -query "$CandFile" -db "$BlastDB" 
blastn -task blastn -query "$CandFile" -db "$BlastDB" -num_threads "$Threads" -outfmt 6 -out "$BlastFile"
perl "$SCRIPT_DIR/ident_problem" "$BlastFile" >| "$BlackFile";
awk '
    (ARGIND == 1){InSet[$1]=1;next}
    /^>/{
        p=0;
        if(InSet[substr($0,2)]){
            p = 1;
        }
    }
    (p==1)
' "$BlackFile" "$CandFile"

rm -rf "$CandFile" "$BlastFile";
