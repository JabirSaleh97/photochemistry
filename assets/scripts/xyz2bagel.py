#!/usr/bin/env python3
"""
This script generates bagel input file from XYZ coordinates.
"""
import sys
import os
from string import Template

# Define the BAGEL input JSON template
bagel_template = Template('''
{
  "bagel": [
    {
      "title" : "molecule",
      "basis" : "svp",
      "df_basis" : "svp-jkfit",
      "angstrom" : true,
      "geometry" : [
$xyz
      ]
    },
    {
      "title" : "hf"
    },
    {
      "title" : "print",
      "file" : "hf.molden",
      "orbitals" : true
    },
    {
      "title" : "casscf",
      "nstate" : 2,
      "nact" : 2,
      "nclosed" : 7,
      "natocc" : true,
      "maxiter": 200,
      "maxiter_micro": 200,
      "active" : [ 8, 9 ]
    },
    {
      "title" : "print",
      "file" : "casscf.molden",
      "orbitals" : true
    },
    {
      "title" : "optimize",
      "target" : 1,
      "method" : [
        {
          "title" : "caspt2",
          "smith" : {
            "method" : "caspt2",
            "ms" : "true",
            "xms" : "true",
            "sssr" : "true",
            "shift" : 0.2,
            "frozen" : true,
            "maxiter" : 200
          },
          "nstate" : 2,
          "nact" : 2,
          "nclosed" : 7,
          "natocc" : true,
          "maxiter" : 400,
          "maxiter_micro" : 200,
          "active" : [ 8, 9 ]
        }
      ]
    },
    {
      "title" : "print",
      "file" : "final.molden",
      "orbitals" : true
    }
  ]
}
''')

def xyz_to_atoms_bad(xyz_file):
    atoms = []
    try:
        with open(xyz_file, 'r') as f:
            lines = (line.strip() for line in f if line.strip())
            next(lines)  # Skip number of atoms
            next(lines)  # Skip comment line

            for line in lines:
                print(line)
                parts = line.split()
                if len(parts) >= 4:
                    atom = {
                        'element': parts[0],
                        'x': float(parts[1]),
                        'y': float(parts[2]),
                        'z': float(parts[3])
                    }
                    atoms.append(atom)
    except (IOError, ValueError) as e:
        print(f"Error reading {xyz_file}: {e}", file=sys.stderr)
    return atoms

def xyz_to_atoms(xyz_file):
    atoms = []
    try:
        with open(xyz_file, 'r') as f:
            lines = f.readlines()  # Read all lines exactly, no filtering

            if len(lines) < 2:
                raise ValueError("XYZ file is too short.")

            atom_lines = lines[2:]  # Skip first two lines (count and comment)

            for line in atom_lines:
                line = line.strip()
                if not line:
                    continue  # Skip completely empty lines safely
                parts = line.split()
                if len(parts) >= 4:
                    atom = {
                        'element': parts[0],
                        'x': float(parts[1]),
                        'y': float(parts[2]),
                        'z': float(parts[3])
                    }
                    atoms.append(atom)
    except (IOError, ValueError) as e:
        print(f"Error reading {xyz_file}: {e}", file=sys.stderr)
    return atoms


def format_geometry(atoms):
    """
    Format atoms into JSON-style string entries.
    """
    geometry_entries = []
    for atom in atoms:
        entry = f'        {{ "atom": "{atom["element"]}", "xyz": [{atom["x"]}, {atom["y"]}, {atom["z"]}] }}'
        geometry_entries.append(entry)
    return ',\n'.join(geometry_entries)

def save_bagel_input(content, output_file):
    try:
        with open(output_file, 'w') as f:
            f.write(content)
    except IOError as e:
        print(f"Error writing {output_file}: {e}", file=sys.stderr)

def main():
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} input.xyz output.json")
        sys.exit(1)

    input_file, output_file = sys.argv[1], sys.argv[2]

    if not os.path.isfile(input_file):
        print(f"Input file {input_file} does not exist.", file=sys.stderr)
        sys.exit(1)

    atoms = xyz_to_atoms(input_file)
    if atoms:
        geometry = format_geometry(atoms)
        bagel_input = bagel_template.substitute(xyz=geometry)
        save_bagel_input(bagel_input, output_file)

if __name__ == '__main__':
    main()

