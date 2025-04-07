***ethylene SA-2-CAS(2/2)
! One step of S0/S1 conical intersection search in Cartesian coordinates
memory,10,m
FILE,2,but5.wfn
punch,all.pun
geometry={
nosymm,noorient;
C,,        -3.47939045,     -0.23510132,      0.39418477
C,,        -0.98623629,     -0.06566735,      0.69519811
C,,         1.14072498,      0.47728046,     -0.90614600
C,,         3.55350561,      0.49915922,     -0.30754951
H,,        -4.92545923,      1.31018850,      0.69214101
H,,        -4.59523497,     -2.00099479,     -0.06430543
H,,        -0.50491568,     -0.50893756,      2.71799450
H,,         0.62841234,      0.97573875,     -2.82297692
H,,         4.20731343,      0.01635085,      1.55827083
H,,         4.96785745,      0.97701411,     -1.68532983
}
basis=6-31g**
 
INT
 
UHF
 
MULTI
config,csf;
OCC,17;
CLOSED,13;
wf,30,0,0
STATE,4;
MAXITER,100;
CPMCSCF,GRAD,  4.1,record=5100.1
CPMCSCF,GRAD,  3.1,record=5101.1
CPMCSCF,GRAD,  2.1,record=5102.1
FORCES
SAMC,5100.1
FORCES
SAMC,5101.1
FORCES
SAMC,5102.1
---
