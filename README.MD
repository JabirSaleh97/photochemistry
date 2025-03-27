<details>
<summary><strong>📘 Getting Started with MobaXterm</strong></summary>

## 📥 Download and Install
1. Go to the [MobaXterm official website](https://mobaxterm.mobatek.net/).
2. Download the **Home Edition** (Installer or Portable version).
3. Install or unzip the package and launch `MobaXterm`.

## 🔌 Connect to a Remote Server (via SSH)
1. Open MobaXterm.
2. Click on **Session** (top left corner).
3. Select **SSH**.
4. Fill in the remote host details:
   - **Remote host**: `ui.wcss.pl`
   - **Specify username**: your-username
5. (Optional) Check **Use private key** if you're using an SSH key.
6. Click **OK** to connect.

## 📁 Upload/Download Files
- Once connected via SSH, a file browser will appear on the left side.
- You can drag and drop files between your local machine and the server.

## 🛠 Tips
- Use the built-in terminal for Linux commands.
- Save sessions for quick reconnects.
- You can also use tools like `scp`, `rsync`, or Git directly in the terminal.

---

> 🧠 **Note**: MobaXterm is available only on Windows. For macOS or Linux, consider using Terminal or other SSH clients.

</details>

<details>
<summary><strong>🐧 Basic Linux / UNIX Commands</strong></summary>

## 📄 File Permissions in UNIX

Each file in UNIX has:
- a type (directory, file, or link),
- access rights for:
  - the owner,
  - the group,
  - others.

Example output of `ls -l`:

```bash
drwxr-xr-x 2 rgora kdm 2048 Mar 26 16:39 n2h2/
-rw-r--r-- 1 rgora kdm 4552 Mar 26 12:12 scan.dat
lrw-r--r-- 1 rgora kdm   75 Mar 26 12:14 scan.gpl
```

Legend:
- `d` – directory, `-` – regular file, `l` – symbolic link
- `r` – read, `w` – write, `x` – execute

File permissions are displayed as:
```bash
-rwx r-x r--
|---|---|---|
  U   G   O
```
U = User, G = Group, O = Others

## 🔧 Essential Commands

```bash
# Navigation
cd                # Go to home directory
cd /path/to/dir   # Change to specific directory

# File listing
ls                # List files
ls -la            # Long list incl. hidden files
ls -1 /dir        # One file per line

# File removal
rm file           # Delete a file
rm -r folder      # Recursively delete a folder

# Copy/move files
cp file1 file2    # Copy file1 to file2
cp -r dir1 dir2    # Copy directory
mv file newname   # Rename or move file

# Create/delete directories
mkdir new_folder
rmdir empty_folder

# File viewing
cat file.txt      # Print content
more file.txt     # Page through file
less file.txt     # Like more, but scrollable

# Changing permissions
chmod +x file     # Make file executable
chmod a-w file    # Remove write permission for all
chmod ugo=rx file # Read/execute access to all
```

## 💡 Tips

- Use `man <command>` or `<command> --help` to read about commands.
- Use wildcards:
  - `*` matches any string.
  - `?` matches a single character.
- Use `|` (pipe) to connect commands:
  ```bash
  ls -l | less
  ```
- Use `>` and `<` to redirect output:
  ```bash
  ls -l > output.txt
  sort < output.txt
  ```

</details>

<details>
<summary><strong>📝 Modern <code>vim</code> Basics</strong></summary>

## ✍️ Modes in `vim`

- `Esc` – Command mode
- `i` – Insert before cursor
- `a` – Append after cursor
- `o` / `O` – Open new line below / above
- `v` – Visual selection mode (character-wise)
- `V` – Visual line selection
- `Ctrl + v` – Visual block selection (column mode)

## 🧭 Movement (you can also use arrow keys)

- `↑ ↓ ← →` – Move cursor
- `0` / `$` – Start / end of line
- `gg` / `G` – Top / bottom of file
- `:n` – Go to line `n`
- `Ctrl + G` – Show file and line info

## 🧾 Editing & Copy-Paste

- `x` / `X` – Delete character under / before cursor
- `dd` – Delete current line
- `yy` – Copy line
- `p` – Paste after cursor
- `u` – Undo
- `Ctrl + r` – Redo

## 🟦 Block Editing

- `Ctrl + v` – Start block selection
- Use arrow keys to expand selection
- `d`, `y`, `p`, etc. work with blocks

Example:  
- Select a vertical block with `Ctrl + v`, move with arrows, press `I` to insert text at the beginning of every selected line (then press `Esc` twice to apply).

## 🔍 Search and Replace

- `/pattern` – Search forward
- `?pattern` – Search backward
- `n` / `N` – Next / previous match
- `:s/old/new/g` – Replace in current line
- `:%s/old/new/g` – Replace in entire file

## 🖫 Save and Exit

- `:w` – Save
- `:q` – Quit
- `:wq` or `ZZ` – Save and quit
- `:q!` – Quit without saving

## 🪟 Working with Splits

- `:split` or `:sp filename` – Horizontal split
- `:vsplit` or `:vsp filename` – Vertical split
- `Ctrl + w + w` – Switch between windows
- `Ctrl + w + q` – Close current split
- `Ctrl + w + =` – Make splits equal size

> 💡 Tip: To open files in splits from the command line:  
> `vim -O file1 file2` (vertical)  
> `vim -o file1 file2` (horizontal)

### 📝 Vim Cheat Sheet

📄 [Download the Cheat Sheet PDF](./assets/mysheet.pdf)

</details>

<details>
<summary><strong>🔐 Remote Access with SSH</strong></summary>

## 🔌 Connect to Remote Server

```bash
ssh -Y your_username@ui.wcss.pl
```

- `-Y` allows GUI applications over SSH

## ⚙️ Execute Remote Commands

```bash
ssh ui.wcss.pl "qstat -u $USER"
```

## 📂 File Transfers with SCP

```bash
scp file.txt ui.wcss.pl:             # upload to home directory
scp -r folder/ ui.wcss.pl:target/    # upload directory
scp ui.wcss.pl:file.txt ./           # download
```

## 🔑 SSH Key Setup

```bash
ssh-keygen -t rsa
```

Then copy your public key:

```bash
cat ~/.ssh/id_rsa.pub | ssh your_username@ui.wcss.pl "cat >> ~/.ssh/authorized_keys"
```

- Next time you connect: no password needed.
- Use `ssh-add` if your key is passphrase protected.

</details>

<details>
<summary><strong>📋 Job Scheduling with PBS</strong></summary>

## 🔎 Check Queue Status

```bash
qstat
qstat -u $USER           # your jobs only
qstat -f Job_ID          # detailed info
pbsnodes -a              # check compute nodes
```

## 🚀 Submit Interactive Job

```bash
qsub -I -X -l select=1:ncpus=4:mem=4gb -l walltime=6:00:00
```

## 🧾 Submit Batch Job with Script

```bash
cat << EOF | qsub -r n
#!/bin/bash
#PBS -N my_job
#PBS -l ncpus=4
#PBS -l mem=2mb
#PBS -q main
#PBS -m n
#PBS -l walltime=1:00:00

cd ~
ls -l >& output.txt
EOF
```

- Use `qstat -u $USER` to check job status: `Q` = queued, `R` = running

</details>

<details>
<summary><strong>📅 Job Scheduling with SLURM</strong></summary>

## 📊 Check Job and Node Status

```bash
squeue                     # list current jobs
squeue -u $USER            # your jobs
sinfo                      # info about nodes
```

## 🚀 Submit a Job

```bash
sbatch job_script.sh
```

Example `job_script.sh`:

```bash
#!/bin/bash
#SBATCH --job-name=my_job
#SBATCH --output=result.out
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --time=01:00:00

cd $HOME
ls -l > result.txt
```

## 🧪 Interactive Job Session

```bash
srun --ntasks=1 --cpus-per-task=4 --mem=4G --time=01:00:00 --pty bash
```

## ❌ Cancel Job

```bash
scancel <job_id>
```

## 📊 WCSS has a comprehensive users' manual at [man.e-science.pl](https://man.e-science.pl/en/home)

</details>


