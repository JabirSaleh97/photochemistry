***ethylene SA-2-CAS(2/2)
! One step of S0/S1 conical intersection search in Cartesian coordinates
memory,25,m
FILE,2,acro.wfn
punch,all.pun
geometry={
nosymm,noorient;
C,,        -1.23007789,     -0.47542260,     -0.00040535
C,,         1.66074825,     -0.15978637,     -0.00046894
C,,         2.79274466,      1.68394950,      0.00013967
O,,         5.42693612,      2.08671374,      0.00040704
H,,        -2.04530608,     -2.21773564,     -0.00066035
H,,        -2.28743227,      1.19202982,      0.00037352
H,,         2.60400056,     -1.98432869,     -0.00092076
H,,         2.86077572,      4.03161234,      0.00084190
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
