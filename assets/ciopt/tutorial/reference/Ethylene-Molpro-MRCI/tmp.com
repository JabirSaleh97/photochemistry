***ethylene SA-2-CAS(2/2)
! One step of S0/S1 conical intersection search in Cartesian coordinates
memory,50,m
FILE,2,jkl.wfn
punch,all.pun
geometry={
nosymm,noorient,angstrom;
C,,         0.04100073,      0.09369438,     -0.44185365
C,,        -0.29480700,     -0.08739467,      0.95135555
H,,         1.10937842,     -0.14303472,     -0.54485365
H,,        -0.09069704,      1.11745263,     -0.80372136
H,,        -0.49973689,     -0.57081966,     -1.12302139
H,,        -0.10163859,     -0.32631263,      1.96904732
}
basis=6-31g**
 
INT
HF
MULTI
config,csf;
OCC,9;
CLOSED,7;
WF,16,1,0;
STATE,2;
WEIGHT,1,1;
MAXITER,100;
 
MRCI
OCC,9;
CLOSED,7;
WF,16,1,0;
STATE,2;
---
