import os

# Path to your input .xyz file
input_file = "liic.xyz"

# Create a directory to store the output frames
output_dir = "xyz_frames"
os.makedirs(output_dir, exist_ok=True)

# Read the entire file
with open(input_file, "r") as file:
    lines = file.readlines()

# Determine frame size from the first line
frame_size = int(lines[0].strip()) + 2  # number of atoms + 2 header lines
num_frames = len(lines) // frame_size

# Split and save each frame
for i in range(num_frames):
    frame_lines = lines[i * frame_size: (i + 1) * frame_size]
    output_filename = f"frame_{i + 1:04d}.xyz"
    output_path = os.path.join(output_dir, output_filename)
    with open(output_path, "w") as out_file:
        out_file.writelines(frame_lines)

print(f"Extracted {num_frames} frames into '{output_dir}' directory.")

