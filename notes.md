# Shell commands
## Hotkeys

- `ctrl + l`
	- clears the terminal window (= shell output window)
- `ctrl + d`
	- Ends the input
- `ctrl + c`
	- Interrupts the program

## Operators

- Pipe symbol `|`
	- Sends output from one program to another
- `<` symbol
	- sends the content of the file after the `<` into the standard input of the next program. E.g. `wc < file.py`.
- `>` 
	- write (and potentially overwrite) output in a file
- `>>`
	- append output to a file

## Programs
- `man "command"`
	- Shows the manual page of a shell command
	- Use space to move down a page and `q` to exit the manual
- `ls`
	- = List
	- Lists files and folders in the current directory
	- options
		- `-l` prints a long list
- `pwd`
	- =print working directory
	- Prints the current directory
- `cd "directory"`
	-  = change directory
	- options 
		- `..` to move one folder up
		- `-` navigate to the last directory
- `ln -s "existing-thing" "new-name"`
	- = link
	- generates a link
	- options
		- `-s` generates a symbolic link. As opposed to a hard link, which we don't need.
		- "existing-thing" is the thing (directory or file) the link should point to
		- "new-name" is the name of the link (optional, will be the same if not given)
- `rm "file or folder"`
	- =remove
- `cat "file"`
	- =catenate
	- show the contents of the file in the terminal
- `less`
	- shows content of a file in a paginated way (page by page)
- `head`
	- shows you the first `n` (default 10) lines of a file
	- options
		- `-n` defines the number of lines to be printed (e.g., `head -n 20 file`)
- `grep "pattern" "file"`
	- maybe stands for global regular expression and print
	- outputs lines that match the pattern in the file
	- Example
		```shell
		# find GGG in FASTA
		grep GGGG data1.fasta
		
		# sort our FASTA data ids
		grep id data1.fasta | sort
		```
	- options
		- `-v` print non-matches instead of matches
- `awk`
	- awkward little programming language.
- `cut -f "field number" -d "seperator"`
	- seperates an input according to a "seperator" and prints out the "field number"
- `gzcat`
	- Uncompresses file and writes the content to stdoutput
- `time`
	- times the command that follows
- `cp old_name new_name`
	- = copy
- `ssh-add`
	- Adds a ssh-key to the key manager. You can confirm that the agent runs via 
- `tr "from" "to"`
	- `tr a-z A-Z` translates all small letters to large letters
	- options
		- -d deletes the following characters from the input
		- -dc delete everything but
- `comm file1 file2`
	- compares 2 (sorted) list. The output is 3 columns, the first is lines unique to file1, the second to file2 and the third lines that occur in both 

### More shell commands

Metacharacters - *, [], etc.

* `*` is for matching multiple things
* `{word1,word2,word3}xxx` expands to the three words, each followed by `xxx`.
* `[xyz]` matches any of `x` `y` or `z`.

`echo $PATH` - see the directories in your shell's path (where it looks for
commands). I.e., `PATH` is a list of directories that your shell will look
in when you type a command.

Variables are not expanded inside hard quotes (`'`) but they are in soft
quotes (`"`).

Add two directories to the start of your PATH

```sh
PATH=.:$HOME/bin:$PATH
```

float are floating point numbers.

* chmod - change mode of a file (or directory)
* history
* mkdir - make a new directory
* cd `.` and `..`
* cp, rm, mv, mkdir
* grep
* cut
* echo
* cat
* less
* tr
* awk
* sed

* perl
* test
* make
* tmux
* paste
* comm


### Shell command line editing keystrokes

* Control-a takes you to the beginning of the line
* Control-e takes you to the end of the line
* Control-r incremental reverse search of your command history
* Control-y yank back killed (deleted text)
* Alt-. insert the last word from the previous command
* Control-k kill (delete) the characters to the end of the line

* Control-f move forward one character
* Alt-f move forward one word

* Control-b move backward one character
* Alt-b move backward one word

* Control-d delete one character
* Alt-d delete one word (after the cursor)
* Alt-backspace delete the word before the cursor



# Git
## Cloning a repository
The first time you download a local copy of a repositiory you have to use the command
```shell
git clone git@github.com:path-to-repo
```

## Pull changes
To get the newest version of the repo use
```shell
git pull
```

## add files to version control
```shell
git add "file or folder"
```

## status
```shell
git status
```

## commit changes
```shell
git commit
```
Options
	- -a commits all files that have ever been added via git add 

## push changes to the github repo
Can only be done after a commit
```shell
git push
```

## add new files

# Python

## Core functionality
### Error handling / Exceptions
Exceptions are unexpected errors that make the program stop. You can catch errors by
```python
try:
	number1, number2 = command.split()[1:]
except ValueError:
	print("There was a value error")
```


### Data structures and types
#### sets

Sets are grups of items that cannot contain duplicates. It is very fast to check if an item is already in a set.
Example

```python
 names = set()
 names.add("Terry")
 names.add("Eva")
 names.add("Eva") # will not be added again
```

#### dictionaries

Dictionaries are a bunch of keys paired to values.
The dictionary has a hash table associated to it or is a hash table or whatever, but this makes the lookup really fast
Example:
```python 
import sys

# To create a dictionary write
wordCoundts = {}

# the program reads the console input and loops though its lines
for line in sys.stdin:
	# it then loops through lines 
	for word in line.split():
		if word in wordCounts:
			# To access a dictionary by its key, write the dictionary name and in [] the key
			wordCounts[word] = wordCounts[word] + 1 
		else:
			wordCounts[word] = 1
```


#### list

```python
example_list = []
```

Lists can be accesses by writing the variable name and in `[]` the number of the item, for example `example_list[1]`. You can also access sublists by proving a range, like `example_list[0:3]` or to access the last elements starting from element 4 `example_list[3:]`  

Methods
	- append() adds an element of the list 

#### strings
Methods
	- split(seperator) 
	- ", ".join(list)   merges the elemtns of the list with the string between items of the list
	- .lower() returns lowercase of the string
	- .upper() returns uppercase of the string

#### integers
Full numbers.
Strings can be converted into integers by `int()`

#### float
Natural number.
Strings can be converted into integers by `float()`

### Custom functions

```python
#define the function
def calculateToFarenheit(c):
	return (9/5)*c+32

# call the function
print("100??C is equal to: ", calculateToFarenheit(100), "?? F")
```

### import

load a package. You can either load a complete package or a module from a package 

```python
# import complete package
import Bio
Bio.SeqIO.parse(sys.stdin, "fasta")

# Import one module only
from Bio import SeqIO
SeqIO.parse(sys.stdin, "fasta")
```

### range function

`range()` returns a so-called iterator that provides the next number in the
range everytime its called.  Python does not include end_value so for
example `range(10, 20)` is the number range \[10, 19\].

```python
for number in range(start_value,end_value):
	print(number)

#or
for number in range(end_value):
	print(number)

```

### len function
`len()` outputs the length of a list, dictionary, set or string


### input

```python
variable = input('Input prompt')
```
Waits for input from the console. Once the input is given, the input is saved in variable

### f-strings

f-strings can contain variables in curled brackets. They are written with
an f before the quotes surrounding the string

```python
print(f"Variable 1 has the value: {variable_1}")
```

To print a raw representation of a variable put `!r` after a variable

```python
print(f"Variable 1 has the value: {variable_1!r}")
```

## Flow of control

### for loop

Example:

```python
for line in sys.stdin: 
	print(line)
```


You can also nest loops like this

```python
for number1 in 1,2,3,4,5:
	print("Number 1: ", number1)
	for number2 in 1,2,3,4:
		print("   ","Number 2: ", number2)
```
### while loop
```python
while count < 20 :
	print(count)
	count = count + 1
```
A loop can be exited by using the `break` statement and one iteration can be skipped by using the `continue` statement.

### if statement

Example:

```python
if (line.startswith(">")):
    print(line)
```

### in statement

with the in statement it can be checked if an element is in a data
structure

```python
ids = set()
if line in ids:
	print("Found duplicate: ", line, end="")
else
	ids.add(line)
```

## sys package

The sys package contains some basic functions to interact with the system

```python
import sys
```

### sys.stdin

Stands for standard input and is used to receive input from the terminal

## Argparse package

A package to accept input from the console

Define argument parser
```python
parser = argparse.ArgumentParser(description = "Filter FASTA files")
```

add arguments to the parser
```python
parser.add_argument(
	"--translate", action="store_true", help="Translate the sequence to AA"
)
```

read an argument from the parser
```python
print("translate: ", args.translate)
```

Provide arguments from the shell
```shell
python3 script.py --translate
```
This will print "translate True"


## Biopython package

### SeqIO module

#### SeqIO.parse(input, format)
