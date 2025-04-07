##🚀 Turbomole Quick Start Tutorial
<details>
<summary><strong>⚙️ Geometry optimization at MP2 level</strong></summary>

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

### ⚛️ 3. Run a Single Point HF Calculation

We could now run preliminary HF calculations interactively:

```bash
dscf > dscf.out
```
or submit this job to a queue:



### 🔁 5. Optimize the Geometry (Optional)

```bash
jobex > jobex.out
```

Use `jobex -ri` for RI-DFT geometry optimizations (faster).

### 📊 6. Analyze Results

- Check `energy` for total energy
- Look into `control` and output files for further analysis
- Use tools like `aoforce`, `eiger`, or `escf` for vibrational and electronic properties

### 🧼 7. Clean Up

```bash
rm -f *.tmp *.bak
```

---

### 📚 Resources

- [Turbomole Documentation](https://www.turbomole.org)
- `man define`, `man ridft`, etc., for module-specific help

</details>







<details>
<summary><strong>✨ Location of uracil MECPs using ADC(2) method</strong></summary>

### 📄 1. Create a coord File

We start by preparing a starting geometry for CI denoted u01a in Lan's paper. Once we have the `uracil_ci_u01a.xyz` file prepared:

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
we may convert it to TM format:
```bash
sub-interactive
module load turbomole/7.6
x2t uracil_ci_u01a.xyz > coord
```

### ⚙️ 2. Setup calculations with `define`

Using the `define` interactively is rather tedious, thus it is convenient to prepare a plain text file containing the required commands and redirect it to STDIN. For S1 excited state optimization we could use the following `adc2.def` file:

```


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
maxiter 100
geoopt adc(2) (a 1)
*
exci
irrep=a nexc=1
*
*
*
```

```bash
define < adc2.def
```

Inside `define`, do the following:

1. Load coordinates (`a coord`)
2. Define atom types (`*`)
3. Choose basis set (e.g., `basis -> def2-TZVP`)
4. Enable resolution of identity (`ri`)
5. Set method as MP2 (`mp`)
6. Save and exit (`*`, then `yes`)

Then edit `control` manually to switch from MP2 to ADC(2):

```bash
$scfmo    mos
$closed shells
a b
$adc
  adc(2)
$soes
  5
$excitations
  singlet
```

- `$soes 5` → number of excited states to compute
- `singlet` → restricts to singlet states

You can also use `esprep` to help prepare this.

---

## 💥 2. Excited State Spectrum Calculation

Use the `escf` module to compute excitation energies and oscillator strengths:

```bash
escf > escf.out
```

Key outputs:

- **Excitation energies** (in eV)
- **Oscillator strengths** for spectra
- **Transition dipole moments**

📌 _Plotting tip_: Extract excitation energies + oscillator strengths and plot a stick spectrum using Python, Origin, etc.

---

## 📉 3. Excited-State Geometry Optimization (ADC(2)-S1)

Turbomole can optimize the first excited state with:

```bash
jobex -level cc2 -es x -statpt > jobex_s1.out
```

Where `-es x` means optimize the **first** excited state (S1). To optimize higher states:

```bash
jobex -level cc2 -es 2 -statpt
```

Also make sure the following are in your `control` file:

```bash
$statpt
  update
$exopt
  adc(2)
```

🔁 Repeat until convergence. Monitor energy and gradients.

---

## 🔬 4. Analyze Results

- Optimized excited-state geometry in `geoopt.xyz`
- Excited-state energies in `escf.out`
- Use `tm2x` or external tools for visualization (e.g., orbitals, transitions)

---

## 🧽 Cleanup

```bash
rm -f *.tmp *.bak
```

---

## 📚 References & Tips

- [Turbomole User Manual](https://www.turbomole.org)
- Use `ricc2` instead of `escf` for large systems with ADC(2)
- Set parallel execution via `export PARA_ARCH=SMP` and `export TURBOTMPDIR=/tmp`

</details>
