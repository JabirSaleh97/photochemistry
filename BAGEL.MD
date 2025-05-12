## 🚀 BAGEL Quick Start Tutorial

<details>
<summary><strong>⚙️ Geometry optimization at XMS-CAPT2 level of ethylene in the ground (N) and excited (V) state and location of conical intersetion geometries</strong></summary>

### 📦 Prerequisites

Before you begin, make sure you have:

- Access to a Unix/Linux shell
- Basic knowledge of terminal commands
- BAGEL package installed. It's available as a package in Debian and Ubuntu (`apt install BAGEL`) alternatively it may be downloaded and compiled - check out the official [installation guide](https://nubakery.org/quickstart/installation_guide.html).
- At WCSS you may copy a precompiled package from my directory. Using the same directory structure you won't have to edit the scripts.
```bash
mkdir -p ~/appl/bagel
cp -r ~rgora/appl/bagel/1.2.0-patch ~/appl/bagel
cp ~rgora/bin/Bagel ~/bin
cp ~rgora/bin/xyz2bagel.py ~/bin
```
- MOLDEN package installed

### 🧪 1. Create a Working Directory

```bash
mkdir -p ethylene/s0 && cd ethylene/s0
```

### 📄 2. Create an Input File

We need a starting geometry for ground state optimization. In such a case I often use [molget](https://github.com/jensengroup/molget) package from Jan Jensen's group (requires access to Open Babel package).
```bash
sub-interactive
module load openbabel
molget ethylene
exit
```

