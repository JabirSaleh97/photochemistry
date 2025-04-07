***ethylene SA-2-CAS(2/2)
! One step of S0/S1 conical intersection search in Cartesian coordinates
memory,50,m
FILE,2,uio.wfn
punch,all.pun
geometry={
nosymm,noorient,angstrom;
C,,        -0.00420004,      0.08511218,     -0.43764416
C,,        -0.33680536,     -0.09715874,      0.94681648
H,,         1.06913897,     -0.16610410,     -0.52387809
H,,        -0.11059882,      1.11072186,     -0.80564341
H,,        -0.53205718,     -0.57947709,     -1.12921215
H,,        -0.08547758,     -0.35209412,      1.95056131
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
 
RS2,shift=.3,mix=2,root=  2
OCC,9;
CLOSED,7;
WF,16,1,0;
STATE,2;
 
forces
 
RS2,shift=.3,mix=2,root=  1
OCC,9;
CLOSED,7;
WF,16,1,0;
STATE,2;
 
forces
---
