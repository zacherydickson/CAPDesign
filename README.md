#CAPDesign

Program to design a minimal set of multiplex primers to amplify a set of targets.

#Installation

```
    cmake src
    make 
```

#BasicUsage

```
    ./CAPDesign targets.fasta
```

The fasta file is assumed to have a set of sequences with flanks. The targets are assumed to be approximately centered in the provided sequences with approximately the same amount of padding.

