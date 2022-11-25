## Here are some high-level git notes

### Seeing status

`git status` lists known files with changes, unknown files, deletions,
renames, etc.

### Adding a file

Use `git add` to tell git that you want it to look after a file for you.

### Where are the remote copies?

```sh
$ git remote -v
origin	git@github.com:VirologyCharite/2022-berlin-python-course.git (fetch)
origin	git@github.com:VirologyCharite/2022-berlin-python-course.git (push)
```

### Save your changes

`git commit`.  Use `-a` to commit changes to all files that have previously
been added via `git add`. Use `-m` to give a log message.

### See what has changed

`git diff`

### Throw away local changes to a file

`git checkout -- FILENAME`

Warning: do this at your own risk!

`git reset --hard` - do this with all files.

### Make a new branch

Here called "add-documentation". The `-b` tells `git checkout` to create
the branch.

`git checkout -b add-documentation`

### Switch to a branch

`git checkout BRANCH-NAME`

### Switch back to previous branch

`git checkout -`

### Delete a branch locally 

`git branch -D BRANCH-NAME`
