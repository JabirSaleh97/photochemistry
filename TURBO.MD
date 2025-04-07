<details>
<summary><strong>🚀 Turbomole Quick Start Tutorial</strong></summary>

A short guide to help you get started with **Turbomole**, a powerful quantum chemistry software package.

### 📦 Prerequisites

Before you begin, make sure you have:

- Access to a Unix/Linux shell
- Turbomole installed and the environment set (e.g., `source $TURBODIR/Config_turbo_env`)
- Basic knowledge of terminal commands

### 🧪 1. Create a Working Directory

```bash
mkdir h2o_test && cd h2o_test
```

### 📄 2. Create a Coord File

Create a file named `coord` with the geometry in Turbomole format or use `x2c` program to convert xyz file"

```bash
x2c 
```
```
$coord
  0.000000000  0.000000000  0.000000000  o
  0.758602000  0.000000000  0.504284000  h
 -0.758602000  0.000000000  0.504284000  h
$end
```

### ⚙️ 3. Run `define` to Set Up the Calculation

```bash
define
```

Follow the menu to:

- Read coordinates (`a coord`)
- Specify the atomic types (`*`)
- Choose basis set (e.g., `basis -> def2-SVP`)
- Select method (e.g., DFT -> `d`)
- Save and exit (`*`, `yes` to save)

### ⚛️ 4. Run a DFT Single Point Calculation

```bash
ridft > ridft.out
```

Or for other modules:

```bash
dscf > dscf.out   # For HF or post-HF methods
```

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

Using the `define` interactively is rather tedious, thus it is convenient to prepare a plain text file containing the required commands and redirect it to STDIN:

for S1 excited state optimization we could use the following `adc2.def` file:

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
