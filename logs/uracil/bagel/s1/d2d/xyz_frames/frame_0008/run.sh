#!/bin/bash

# Columns x Rows
tile="5x"


# Start jmol without display
jmol -o -g 1500x1500 -n casscf.molden.spt
# Crop the generated MO pictures
mogrify -verbose -trim +repage mo_*.png

# Create the montages
# PNG with transparent background

montage -verbose -label "MO 8 (1.49)" mo_8.png -label "MO 9 (0.51)" mo_9.png \
-geometry "600x600>" \
-shadow \
-pointsize 60 \
-background none \
montage.0.PNG

# JPG with white background
mogrify -verbose -background white -alpha remove -format jpg montage.0.PNG




