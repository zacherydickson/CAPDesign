#CAPDesign

Program to design a minimal set of multiplex primers to amplify a set of targets.

#Installation

```
    git clone --recurse-submodules https://github.com/seqan/seqan3.git
    cmake src
    make 
```

As a note The current version was developed with seqan3.1, and ranges-v3 0.11

#BasicUsage

```
    ./CAPDesign targets.fasta
```

The fasta file is assumed to have a set of sequences with flanks. The targets are assumed to be approximately centered in the provided sequences with approximately the same amount of padding.

For detailed options use `./CAPDesign -h`

