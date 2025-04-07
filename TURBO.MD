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

Create a file named `coord` with the geometry in Turbomole format:

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
