# ORCA Tutorial: Geometry Optimization & UV-vis Spectra of Uracil, Cytosine and Thymine

This tutorial outlines the steps to optimize the geometry of U, T and C nucleobases and compute their vertical excitation spectra using ORCA 6.0.1 package.

## Prerequisites

- **ORCA Installation**: Ensure ORCA is installed on your system. [Installation Guide](https://www.faccts.de/docs/orca/6.0/tutorials/installing-orca.html)
- **Input Files**: Prepare the molecular structure of uracil in an appropriate format (e.g., `.xyz` or `.pdb`).

## Step 1: Geometry Optimization

Create an input file named `uracil_s0.inp` with the following content:

```plaintext
! pbe0 d4 def2-SVP Opt RIJCOSX Def2/J defgrid3 TightSCF

%geom
    MaxIter 200
    #Convergence tight
end

* xyz 0 1
C   0.7161260117   0.9937054230   0.5446110515  
N   1.5295535430   0.0964202102  -0.1547034571  
C   1.1486238151  -1.1902113468  -0.5167547652  
N  -0.1416009934  -1.5743016712  -0.1391606333  
C  -1.0041440815  -0.7553883162   0.5491894544  
C  -0.6222418826   0.5095480268   0.9031977090  
O   1.8970629335  -1.9704453009  -1.1408395844  
H  -1.2681643439   1.1883874522   1.4453994906  
O   1.1473685821   2.1235378362   0.8224609717  
H   2.4633281615   0.3946765571  -0.4201853029  
H  -0.4305344904  -2.5115764067  -0.3988396958  
H  -2.0011471316  -1.1451087139   0.8039473321  
*
```

**Explanation**:

- `! pbe0 d4 def2-SVP Opt`: Specifies the use of the PBE0 hybrid functional with the def2-SVP basis set and D4 empirical dispersion for geometry optimization.
- The `* xyz 0 1` block defines the molecular geometry in Cartesian coordinates with a charge of 0 and a multiplicity of 1.

Run the optimization. We could either start an interactive session, load appropriate module and run interactive job:

```bash
sub-interactive
module load ORCA/6.0.1_avx2
orca uracil_s0.inp > uracil_s0.out &
tail -f uracil_s0.out
```
or submit orca job to the SLURM queuing system. We shall normally use the latter approach. In this case the geometry we'll be given in an external file and we'll use 16 cores to speed the calculations and use up to 2000mb per core:

```plaintext
%pal nprocs 16 end
%maxcore 2000
! pbe0 d4 def2-SVP Opt RIJCOSX Def2/J defgrid3 TightSCF

%geom
    MaxIter 200
    #Convergence tight
end

* xyzfile 0 1 uracil_s0.xyz 

```

where `uracil_s0.xyz` is a regular XYZ file. 

We may submit the job usng `sub-orca-6.0.1-avx2` command, which issued without parameters returns the syntax:

```bash
sub-orca-6.0.1-avx2

orca-6.0.1-avx

Usage: /usr/local/bin/bem2/sub-orca-6.0.1-avx2 FILE PARAMETERS
Parameters:
	-p PARTITION 		 Set partition (queue). Default = normal
	-n NODES 		 Set number of nodes. Default = 1
	-c CORES 		 Up to 48. Default = 6
	-m MEMORY 		 In GB, up to 180 (must be integer value). Default = 30
	-t TIME_LIMIT 		 In hours. Default = 12
```

So to submit our job we should issue the following command (16 cpus times roughly 2gb gives 32gb and we extend the timelimit to 72h of walltime):

```bash
sub-orca-6.0.1-avx2 uracil_s0.inp -c 16 -m 32 -t 72
```

This produces the following output:

```bash
orca-6.0.1-avx

Submitted batch job 3469149
The job is being submitted with the following parameters:

	file            uracil_s0.inp
	partition    normal
	nodes        1
	cores per node 16
	memory       32 GB (per node)
	time limit   72 hours
```

We may monitor our job using `squeue` command:

```bash
squeue --user $USER
```
or even fancier;
```bash
squeue --user $USER --Format jobid:10,name,state:10,timeleft:12,minmemory:10,stdout:100'
```

After completion, the optimized geometry can be found in the output file `uracil_s0.out` (actually the sub-orca script is producing *.orc files but please change that to *.out).

## Step 2: Vertical Excitation Spectrum Calculation

Using the optimized geometry, create a new input file named `uracil_spectrum.inp`:

```plaintext
! B3LYP def2-SVP TDDFT

* xyzfile 0 1 uracil_opt.xyz
```

**Explanation**:

- `! B3LYP def2-SVP TDDFT`: Specifies the use of the B3LYP functional with the def2-SVP basis set to perform a Time-Dependent Density Functional Theory (TDDFT) calculation for excited states.
- `* xyzfile 0 1 uracil_opt.xyz`: Imports the optimized geometry from the `uracil_opt.xyz` file.

Run the TDDFT calculation:

```bash
orca uracil_spectrum.inp > uracil_spectrum.out
```

The vertical excitation energies and oscillator strengths will be reported in the `uracil_spectrum.out` file.

## References

For more detailed information, refer to the [ORCA 6.0 Tutorials](https://www.faccts.de/docs/orca/6.0/tutorials/).
