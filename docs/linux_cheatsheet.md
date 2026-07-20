# Linux / vi cheatsheet

Quick reference for common vi and bash commands. Not exhaustive -- covers enough to view, edit, and navigate files on an HPC.

## vi commands

**Insert mode:**
- `i` / `I` -- insert at cursor / line start
- `a` / `A` -- append after cursor / line end
- `o` / `O` -- open new line below / above

**Save and quit:**
- `:w` -- write (save)
- `:q` -- quit
- `:wq` or `:x` -- write + quit
- `:q!` -- force quit, discard changes

**Navigation:**
- `h` / `j` / `k` / `l` -- left / down / up / right
- `w` / `b` -- word forward / backward
- `0` / `$` -- line start / line end
- `gg` / `G` -- file start / file end

**Editing:**
- `x` -- delete character
- `dd` -- delete line
- `yy` -- yank (copy) line
- `p` / `P` -- paste below / above cursor
- `u` -- undo
- `Ctrl+r` -- redo

**Search:**
- `/pattern` -- search forward
- `?pattern` -- search backward
- `n` / `N` -- next / previous match

**Visual mode:**
- `v` -- character select
- `V` -- line select
- then `d` (delete), `y` (yank)

**Command mode:**
- `:set nu` -- show line numbers
- `:set nonu` -- hide line numbers
- `:s/old/new/g` -- find and replace in current line
- `:%s/old/new/g` -- find and replace in file

**Shell escape:**
- `:!command` -- run shell command without leaving vi

## Bash terminal

**File listing:**
- `ls` -- list files/dirs
- `ls -la` -- detailed + hidden files
- `ls -lh` -- human-readable sizes
- `ls -lt` -- sort by modification time

**Directory navigation:**
- `cd path` -- change directory
- `cd ..` -- up one level
- `cd ~` or `cd` -- go to home directory
- `cd -` -- go to previous directory

**Viewing files:**
- `cat file` -- print file contents
- `less file` / `more file` -- page through file (space=next, b=back, q=quit)
- `head -n N file` -- first N lines
- `tail -n N file` -- last N lines

**File management:**
- `mkdir -p path` -- create directory (with parents)
- `cp src dst` -- copy file (`-r`: recursive for directories)
- `mv src dst` -- move or rename
- `rm file` -- remove file (`-r`: recursive, `-f`: force)
- `git clone <url>` -- clone a repository (e.g., `git clone https://github.com/UMN-GDC/GDCGenomicsQC.git`)

**Searching:**
- `grep pattern file(s)` -- search file contents (`-i`: case-insensitive, `-r`: recursive, `-n`: show line numbers)
- `find path -name "pattern"` -- find files by name

**Utilities:**
- `which command` -- show path to executable
- `wc -l file` -- count lines (`-w`: words, `-c`: characters)

**Pipes and redirection:**
- `|` -- pipe output to next command (e.g., `ls -la | grep foo`)
- `>` -- redirect output to file (overwrite)
- `>>` -- redirect output to file (append)

**Keyboard shortcuts:**
- `Ctrl+c` -- interrupt running command
- `Ctrl+d` -- exit shell / end input
- `Tab` -- auto-complete path or command
- `Up` / `Down` arrows -- scroll through command history
