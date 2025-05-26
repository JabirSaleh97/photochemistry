#!/usr/bin/env python3
import sys
import re
from collections import defaultdict

if len(sys.argv) < 2:
    print(f"Usage: {sys.argv[0]} file1 [file2 ...]")
    sys.exit(1)

state_groups = defaultdict(list)

for filename in sys.argv[1:]:
    with open(filename, 'r') as f:
        lines = f.readlines()

    in_fci_block = False
    i = 0
    while i < len(lines):
        line = lines[i].rstrip('\n')

        if '=== FCI iteration ===' in line:
            in_fci_block = True
            states = {}
            i += 1
            continue

        if in_fci_block:
            m = re.match(r'^\s*\d+\s+(\d+)\s+\*\s+(-?\d+\.\d+)', line)
            if m:
                state_label = m.group(1)
                energy = m.group(2)
                states[state_label] = {'energy': energy, 'vectors': []}
                i += 1
                continue

            m = re.match(r'^\s*\* ci vector, state\s+(\d+)', line)
            if m:
                state = m.group(1)
                i += 1
                while i < len(lines):
                    vec_line = lines[i].strip()
                    if not vec_line or vec_line.startswith('*'):
                        i += 1
                        break
                    parts = vec_line.split()
                    if len(parts) >= 2:
                        label, value = parts[0], parts[1]
                        states[state]['vectors'].append((label, value))
                    i += 1

                if state in states:
                    entry_parts = [state, states[state]['energy']]
                    for label, value in states[state]['vectors']:
                        entry_parts.extend([label, value])
                    state_groups[state].append(' '.join(entry_parts))
                continue

        i += 1

# Final output: Global grouping by state label
for idx, state_label in enumerate(sorted(state_groups.keys(), key=int)):
    for entry in state_groups[state_label]:
        print(entry)
    if idx < len(state_groups) - 1:
        print()  # Blank line between groups

