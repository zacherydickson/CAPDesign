double initGibbs = 2.09;

double deltaGibbs[] = {
    -0.29, -0.73, -0.59, -0.19,
    -0.73, -1.16, -1.33, -0.59,
    -0.59, -1.46, -1.16, -0.73,
    0.11, -0.59, -0.73, -0.29
};
//Indexes are based on decimal interpretations of binary encodings of dinucleotides
//  A = 00, C = 01, G = 10, T = 00
//  AA = 0000 -> 0, AC = 0001 -> 1, AG = 0010 -> 2, AT = 0011 -> 3,
//  ...,  GC -> 1001 -> 9, ...
//  TA = 1100 -> 12, ...
//At a temperature of 60 degrees Celsius and a Salinity of 0.18

//Source: Supplement 9, and Table S22 from Xie et al (2022) doi: s41467-022-29500-4
