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
-rwxr-xr--
|---|---|--|
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

