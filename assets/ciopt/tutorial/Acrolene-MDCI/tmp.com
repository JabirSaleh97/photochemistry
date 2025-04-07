***ethylene SA-2-CAS(2/2)
! One step of S0/S1 conical intersection search in Cartesian coordinates
memory,25,m
FILE,2,acro.wfn
punch,all.pun
geometry={
nosymm,noorient;
C,,        -1.23004629,     -0.47542035,     -0.00040534
C,,         1.66073827,     -0.15983205,     -0.00046867
C,,         2.79274382,      1.68405615,      0.00013939
O,,         5.42699582,      2.08667337,      0.00040726
H,,        -2.04529679,     -2.21772428,     -0.00066036
H,,        -2.28742858,      1.19201892,      0.00037311
H,,         2.60392223,     -1.98442716,     -0.00092085
H,,         2.86130937,      4.03176216,      0.00084171
}
basis=6-31g*
 
INT
KS
 
MULTI
config,csf;
OCC,17;
CLOSED,12;
STATE,3;
MAXITER,100;
CPMCSCF,GRAD,  2.1,record=5100.1
CPMCSCF,GRAD,  1.1,record=5101.1
FORCES
SAMC,5100.1
FORCES
SAMC,5101.1
---
