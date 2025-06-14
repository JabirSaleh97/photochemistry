## 🚀 Turbomole Quick Start Tutorial

<details>
<summary><strong>⚙️ Geometry optimization at MP2 level and vertical excitation spectra at ADC(2) level</strong></summary>

### 📦 Prerequisites

Before you begin, make sure you have:

- Access to a Unix/Linux shell
- Turbomole installed and the environment set (e.g., `source $TURBODIR/Config_turbo_env`)
- Basic knowledge of terminal commands

### 🧪 1. Create a Working Directory

```bash
mkdir u0 && cd u0
```

### 📄 2. Create a Coord File

Create a file named `coord` with the geometry in Turbomole format or use `x2c` program to convert xyz file.
We start by preparing a starting geometry for S0 minimum from Lan's paper. Once we have the `uracil_s0.xyz` file prepared:

```
12

C   0.716126   0.993705   0.544611
N   1.529554   0.096420  -0.154703
C   1.148624  -1.190211  -0.516755
N  -0.141601  -1.574302  -0.139161
C  -1.004144  -0.755388   0.549189
C  -0.622242   0.509548   0.903198
O   1.897063  -1.970445  -1.140840
H  -1.268164   1.188387   1.445399
O   1.147369   2.123538   0.822461
H   2.463328   0.394677  -0.420185
H  -0.430534  -2.511576  -0.398840
H  -2.001147  -1.145109   0.803947
```
we may convert it to TM format:
```bash
sub-interactive
module load turbomole/7.6
x2t uracil_s0.xyz > coord
cat coord
```
```
$coord
    1.35328201072375      1.87783029867125      1.02916563445857       c
    2.89043815282585      0.18220739293642     -0.29234630065798       n
    2.17058478017215     -2.24917282051696     -0.97652542353099       c
   -0.26758710897313     -2.97499961745060     -0.26297617722905       n
   -1.89755714968622     -1.42747643782881      1.03781680065710       c
   -1.17586696323939      0.96290616735081      1.70679685630975       c
    3.58492951116094     -3.72360139363823     -2.15587515201806       o
   -2.39648264110992      2.24572596006565      2.73140825080796       h
    2.16821317388575      4.01290523523556      1.55422603818584       o
    4.65501527512215      0.74583143768893     -0.79403457167588       h
   -0.81359134733963     -4.74619078118309     -0.75369836754574       h
   -3.78161976511649     -2.16394239284409      1.51923964871451       h
$end
```

### ⚙️ 2. Setup calculations with `define`

Using the `define` interactively is rather tedious, thus it is convenient to prepare a plain text file containing the required commands and redirect it to STDIN. For S0 ground state optimization, we could use the following `scs-mp2.def` file (please note that the formatting including blank lines is crucial). In this example, we use the SCS variant of MP2 with a def2-SV(P) basis set and EHT guess. We will use `ricc2` program thus we need also an auxiliary basis set (cbas).

```


a coord
sy c1
ired
*
b all SV(P)
*
eht



scf
iter
200

cc
freeze
*
cbas
*
memory
4000
ricc2
model mp2
scs  cos=1.20000   css=0.33333
maxiter 200
geoopt mp2
*
*
*
*
```
We are ready to run `define`

```bash
define < scs-mp2.def
```

### ⚛️ 3. Run a Single Point HF Calculation (optional)

We could now run preliminary HF calculations interactively:

```bash
dscf > dscf.out
```
or submit this job to a queue using [Turbomole](./assets/scripts/Turbomole) (I will be using my scripts available in this repo in [assets/scripts](./assets/scripts) folder).

```bash
To setup turbomole job use:

Turbomole [options] program

To prepare turbomole job (x2t or define)

Turbomole -s x2t coordinates.xyz
Turbomole -s define define_script

Where options are:
  -s toggle sequential/batch execution
  -m memory [32, total: 32]
  -n # nodes [1]
  -p # cores [16]
  -q queue [normal]
  -t wall time [72]
  -x test execution scripts
  -h print this help
  -e print turbomole environment
  -v version: 7.6, 7.8 [7.6]
  -a architecture [smp]
```
Let us submit `dscf` calculations first:
```bash
Turbomole dscf
```
The output is redirected to `$cwd_program.log` file (in this case `u0_dscf.log`).

### 🔁 5. Optimize the Geometry

To optimize the geometry we could submit an interactive `jobex` run 
```bash
jobex -c 200 -level cc2 | tee jobex.out
```
or submit it to a queue
```bash
Turbomole jobex -c 200 -level cc2
```

### 📊 6. Analyze Results

- The successful run produces `GEO_OPT_CONVERGED` file.
- Check `energy` for total energy
- Look into `control`, `job.last` and output files for further analysis
- Use `tm2molden` to visualize orbitals
- Use `t2x > trj.xyz` to check the optimization trajectory

### 🧼 7. Clean Up

```bash
rm -f *.tmp *.bak CC* *cao
```

### ⚙️ 8. Calculate the vertical excitation energies at ADC(2) level

Let us now prepare calculations of vertical excitation energies at the ADC(2) level. We can simply copy the coordinates file with optimized geometry: 
```bash
mkdir vee && cd vee
cp ../coord .
```
Now we need to run define using the following script which we put in `scs-adc2.def` file:

```


a coord
sy c1
ired
*
b all SV(P)
*
eht



scf
iter
200

cc
freeze
*
cbas
*
memory
4000
ricc2
model adc(2)
scs  cos=1.20000   css=0.33333
maxiter 200
*
exci
irrep=a nexc=3 multiplicity=1
irrep=a nexc=3 multiplicity=3
spectrum states=all
exprop states=all
*
*
*
```
We are ready to run `define` and the calculations. First `dscf` and then `ricc2`

```bash
sub-interactive 
module load turbomole/7.6
define < scs-adc2.def
dscf | tee dscf.log
ricc2 |tee ricc2.log
```

Alternatively we can exit the interactive job once the calculations are set up and submit the job to the queue.

```bash
sub-interactive 
module load turbomole/7.6
define < scs-adc2.def
exit # to quit the interactive job and return to access node
Turbomole -v 7.6 ricc2 # The Turbomole script will issue the `dscf` calculations automatically.
```

---

### 📚 Resources

- [Turbomole Documentation](https://www.turbomole.org)

</details>

<details>
<summary><strong>✨ Location of uracil MECPs using ADC(2) method</strong></summary>

### Penalty-constrained MECP optimization of Levine & Martínez

[Levine and coworkers](https://pubs.acs.org/doi/10.1021/jp0761618) developed a penalty-constrained MECP optimization algorithm based on minimizing the objective function:

$$F_{IJ}(\mathbf{R}) = \frac{1}{2}\left[ E_I(\mathbf{R}) + E_J(\mathbf{R}) \right] + \sigma \left( \frac{[E_I(\mathbf{R}) - E_J(\mathbf{R})]^2}{E_I(\mathbf{R}) - E_J(\mathbf{R}) + \alpha} \right)$$

where $\alpha$ is a small fixed smoothing parameter to avoid singularities and $\sigma$ is a Lagrange multiplier for a penalty function that drives the energy gap to zero. Minimization of $F_{IJ}$ is performed iteratively for increasingly large values of $\sigma$. This approach is implemented in [CIOpt](./assets/ciopt) program by Levine and Martínez. The precompiled executable should work on most Linux workstations.


### 🧪 1. Create a Working Directory

```bash
mkdir u01a && cd u01a
```

### 📄 2. Create a coord File

We start by preparing a starting geometry for CI denoted u01a in Lan's paper. Once we have the `uracil_ci_u01a.xyz` file prepared we may proceed further.

```
12

C   0.672439   1.019467   0.655871
N   1.517725   0.066593  -0.034713
C   1.104961  -1.154804  -0.508462
N  -0.220285  -1.589127  -0.124376
C  -0.867678  -0.832307   0.775829
C  -0.748312   0.590670   0.596781
O   1.774326  -1.919924  -1.222101
H  -1.130910   0.918647  -0.403385
O   1.196255   2.047388   1.066268
H   2.493484   0.298659  -0.177008
H  -0.419681  -2.573868  -0.293920
H  -1.496061  -1.261484   1.570978
```

### ⚙️ 2. Setup calculations with `xyz2ciopt.py`

To locate the MECP we will use CIOpt program of Martinez et al. and Turbomole to calculate energies and gradients for the ground state at MP2 level and the S1 state at ADC(2) level. CIOpt extracts these data from specified output files using template files. The preparation of these templates is rather tedious and we shall use my [xyz2ciopt.py](./assets/scripts/xyz2ciopt.py) script (needs numpy which can be installed via `pip install numpy` use it at your own risk):

```
This script generates CIOpt templates from xyz file provided based on the
available template file. If there's no template the default template shall be
generated. If -g option is given the data_file should be the CIOpt log file.

Usage: xyz2ciopt.py [options] data_file

Options:
  -h, --help       show this help
  -r, --rs2        prepare molpro rs2  templates
  -c, --rs2c       prepare molpro rs2c templates
  -x, --xms        prepare molpro xms-rs2 templates
  -e, --eomcc      prepare molpro eom-ccsd templates
  -a, --adc2       prepare turbomole adc2 templates (requires adcmp2.sh in ~/bin dir)
                   (in case of CC2 change last line to:&%08(f20.10)00230)
  -l, --log        in case of adc2 read data from ricc2.log instead of gradient files
  -g, --grep       extract trajectory in xyz format from CIOpt.log
```
Usually the following should do the trick:
```bash
xyz2ciopt.py -a uracil_ci_u01a.xyz
```
You should see the following files in the directory:
```bash
ls
adcmp2.sh*   template.read   template.readg2  template.writeg
Control.dat  template.readg  template.write   uracil_ci_u01a.xyz
```
The `adcmp2.sh` is a bash shell script that runs the ground and excited state calculations. The `xyz2ciopt.py` looks for it in `~/bin` directory; it contains `define` scripts that may need to be tweaked. The `upper_def` is for the ADC(2) energy and gradient calculations for S1 state and `lower_def` is for ground-state MP2 calculations. In this example, we should add the `scs  cos=1.20000   css=0.33333` line to both scripts to use the SCS variants. 

```bash
upper_def() {
cat <<EOF | $Turbomole define


a coord
sy c1
*
no
b all SV(P)
*
eht



scf
iter
500

cc
freeze
*
cbas
*
memory
4000
ricc2
model adc(2)
scs  cos=1.20000   css=0.33333
maxiter 100
geoopt adc(2) (a 1)
*
exci
irrep=a nexc=1
*
*
*
EOF
}

lower_def() {
cat <<EOF | $Turbomole define


a coord
sy c1
*
no
b all SV(P)
*
eht



scf
iter
500

cc
freeze
*
cbas
*
memory
2000
ricc2
model mp2
scs  cos=1.20000   css=0.33333
maxiter 100
geoopt mp2
*
*
*
*
EOF
}

```
The `Control.dat` is an input for CIOpt program which may need some customization but usually the defaults are fine.
```
 &control
 nopt=3
 natoms=12
 nstates=2
 istate=2
 nefunc=7
 dlambdagap=3.5
 alpha=0.02
 tol=1.0d-06
 gtol=5.0d-03
 cigap=0.001
 znoncart=.false.
 zangrad=.true.
 zmultigrad=.true.
 cinpdeck='last.xyz'
 coutfile='gradients'
 crunstr='./adcmp2.sh -v smp -n 1 -p 4 -m 800mb -r last.xyz '
/
C   0.672439   1.019467   0.655871
N   1.517725   0.066593  -0.034713
C   1.104961  -1.154804  -0.508462
N  -0.220285  -1.589127  -0.124376
C  -0.867678  -0.832307   0.775829
C  -0.748312   0.590670   0.596781
O   1.774326  -1.919924  -1.222101
H  -1.130910   0.918647  -0.403385
O   1.196255   2.047388   1.066268
H   2.493484   0.298659  -0.177008
H  -0.419681  -2.573868  -0.293920
H  -1.496061  -1.261484   1.570978
```

### 💥 3. Run the CIOpt

We may now submit the `CIOpt` job to a queue:

```bash
usage: CIOpt [options]

Where options are:
  -c code: molpro, adcmp2 [adcmp2]
  -v code version [7.6] molpro: 2012.1.25, turbo: 7.5
  -a code arch [smp]
  -s toggle sequential/batch execution
  -m memory [32]
  -n # nodes [1]
  -p # cores [4]
  -q queue [normal]
  -x test execution scripts
  -t wall time [72]
  -h print this help"
```

In this case, I'd recommend the following (please note that this script expects `CIOpt.e` executable in particular directory, for instance the SLURM function uses `~/appl/ciopt/bin/CIOpt.e` - you may need to change that to reflect your settings):

```bash
CIOpt -c adcmp2 -p 16
```

## 🔬 4. Analyze Results

- The output of `CIOpt` lands in `CIOpt.log` you may use `xyz2ciopt.py -g CIOpt.log` command to extract the trajectory to `CIOpt.xyz` file.
- The `iter.log` file contains information on the model function energy, step and gradients
- `full.log` contains all the output generated by Turbomole. 
- `last.xyz` contains the last geometry i.e. final MECP geometry upon completion

---

## 🧽 Cleanup

```bash
# these were scratch directories for adcmp2.sh script
rm -rf lower upper 
```

---

## 📚 References & Tips

- [Turbomole User Manual](https://www.turbomole.org)

</details>


<details>
<summary><strong>⚙️ Calculation of TPA spectra of PNA using CC2 method</strong></summary>

### 📦 Prerequisites

Before you begin, make sure you have:

- Access to a Unix/Linux shell
- Turbomole installed and the environment set (e.g., `source $TURBODIR/Config_turbo_env`)
- Basic knowledge of terminal commands

### 🧪 1. Create a Working Directory

```bash
mkdir pna && cd pna
```

### 📄 2. Create a Coord File

To make thing simple we will use molget to prepare initial structure
but this may be available in the structure library of TURBOMOLE. 

```bash
sub-interactive
module load OpenBabel
molget.x pna
```

We should now have the `pna.xyz` file in the current directory:

```
16
C6H6N2O2
C         -0.66049        1.19850       -0.00025
C         -0.66049       -1.19850        0.00025
C          0.72051        1.20210       -0.00025
C          0.72051       -1.20210       -0.00015
H          3.29841        0.84000       -0.00095
H          3.29841       -0.84000       -0.00055
C         -1.35039        0.00000       -0.00025
C          1.41661       -0.00000       -0.00015
O         -3.43939        1.05480       -0.00045
N          2.81341       -0.00000       -0.00055
O         -3.43939       -1.05480       -0.00045
N         -2.83039        0.00000       -0.00035
H         -1.20239        2.13270        0.00425
H         -1.20239       -2.13270        0.00015
H          1.25871        2.13850       -0.00055
H          1.25871       -2.13850        0.00025
```

No we need to create a `coord` file with the geometry in Turbomole format:

```bash
module purge
module load turbomole/7.9
x2t pna.xyz > coord
cat coord
```
```
$coord
   -1.24814520805407      2.26483676036399     -0.00047243153116       c
   -1.24814520805407     -2.26483676036399      0.00047243153116       c
    1.36156657005411      2.27163977441264     -0.00047243153116       c
    1.36156657005411     -2.27163977441264     -0.00028345891869       c
    6.23309154672689      1.58736994468565     -0.00179523981839       h
    6.23309154672689     -1.58736994468565     -0.00103934936854       h
   -2.55186726143339      0.00000000000000     -0.00047243153116       c
    2.67700492540611     -0.00000000000000     -0.00028345891869       c
   -6.49950513577663      1.99328311625526     -0.00085037675608       o
    5.31657437628339     -0.00000000000000     -0.00103934936854       n
   -6.49950513577663     -1.99328311625526     -0.00085037675608       o
   -5.34866192587953      0.00000000000000     -0.00066140414362       n
   -2.27218779498878      4.03021890598938      0.00803133602966       h
   -2.27218779498878     -4.03021890598938      0.00028345891869       h
    2.37861717032770      4.04117931751221     -0.00103934936854       h
    2.37861717032770     -4.04117931751221      0.00047243153116       h
$end
```

### ⚙️ 2. Setup calculations with `define`

Using the `define` interactively is rather tedious, thus it is
convenient to prepare a plain text file containing the required
commands and redirect it to STDIN. For S0 ground state optimization,
we could use the following `cc2.def` file (please note that the
formatting including blank lines is crucial). In this example, we
use the regular CC2 method with a def2-SV(P) basis set and EHT
guess.

```


a coord
sy c1
ired
*
b all SV(P)
*
eht



scf
iter
200

cc
freeze
*
cbas
*
memory
4000
ricc2
model cc2
maxiter 200
geoopt cc2
*
*
*
*
```
We are ready to run `define`

```bash
define < cc2.def
```

### 🔁 4. Optimize the Geometry

In this example we will optimize the geometry assuming PTE variant
of COSMO solvation model.  For that we need to edit the control file
(there is another method using `cosmoprep` but this one is way
simpler) and add the following lines at the top:

```
$cosmo
  gauss
  nleb = 3
  solvent=chloroform
$reaction_field
  post-SCF
  ccs-like
```

To optimize the geometry we could submit an interactive `jobex` run 
```bash
jobex -c 200 -level cc2 | tee jobex.out
```
or submit it to a queue
```bash
Turbomole jobex -c 200 -level cc2
```

### 📊 6. Analyze Results

- The successful run produces `GEO_OPT_CONVERGED` file.
- Check `energy` for total energy
- Look into `control`, `job.last` and output files for further analysis
- Use `tm2molden` to visualize orbitals
- Use `t2x > trj.xyz` to check the optimization trajectory

### 🧼 7. Clean Up

This time we will keep `cosmo_transfer.tmp`
```bash
rm -f *.bak CC* *cao
```

### ⚙️ 8. Calculate the TPA spectrum 

Let us now prepare calculations of two-photon spectra 
at the CC2 level. We can simply copy the `coord`
file with optimized geometry and other relevant files: 

```bash
mkdir tpa 
cp coord auxbasis basis control cosmo_transfer.tmp mos out_cc2.ccf out.ccf tpa
cd tpa
```

Again, we need to modify the control file and add these lines:

```
$excitations
 irrep=a nexc=1
 twophoton states=all operators=(diplen,diplen)
$laplace
 conv = 8
```

Now we are ready to submit the job (in principle `dscf` and then `ricc2` but `Turbomole` script does this by default.

```bash
Turbomole ricc2
```

The two-photon absorption cross section is given by

$$\sigma_{2PA}(\omega) = \frac{4\pi^3 \alpha a_0^5 \omega^2}{c} g(2\omega) \delta_{2PA}$$

where $\omega$ is the angular frequency of the photon ($2\hbar\omega$ corresponds to the transition energy to the excited state (denoted hereafter as $\hbar\omega_f$). One may use either a Lorentzian or Gaussian profile for $g(2\omega)$:
- Lorentzian profile:

    $$g(2\omega) = \frac{1}{\pi} \frac{1}{2\Gamma_f} \frac{1}{(\omega_f - 2\omega)^2 + \left(\frac{1}{2} \Gamma_f\right)^2}$$

    For $\omega_f = 2\omega$:

    $$\sigma_{2PA}\left(\omega = \frac{1}{2} \omega_f \right) = \frac{4\pi^3 \alpha a_0^5 \omega^2}{c\pi(\Gamma_f/2)} \delta_{2PA}$$

    $$\sigma_{2PA}\left(\omega = \frac{1}{2} \omega_f \right) = \frac{8\pi^2 \alpha a_0^5 \omega^2}{c \Gamma_f} \delta_{2PA}$$

- Gaussian profile:

    $$g(2\omega) = \frac{1}{\sigma \sqrt{2\pi}} \exp \left[ -\frac{1}{2} \left( \frac{2\omega - \omega_f}{\sigma} \right)^2 \right]$$

    For $\omega_f = 2\omega$:

    $$g(2\omega) = \frac{1}{\sigma \sqrt{2\pi}} = \frac{2}{\Gamma_f \sqrt{2\ln 2} \sqrt{2\pi}}$$

    $$\sigma_{2PA}\left(\omega = \frac{1}{2} \omega_f \right) = \frac{8\pi^3 \sqrt{2\ln 2} \alpha a_0^5 \omega^2}{c \Gamma_f \sqrt{2\pi}} \delta_{2PA}$$

Comparison of the two line profiles at $\omega_f = 2\omega$ leads to:

$$\frac{2}{\pi \Gamma_f} : \frac{2}{\Gamma_f \sqrt{2\ln 2} \sqrt{2\pi}} \rightarrow 0.637 \cdot \frac{1}{\Gamma_f} : 0.939 \cdot \frac{1}{\Gamma_f}$$


---

### 📚 Resources

- [Turbomole Documentation](https://www.turbomole.org)

</details>

