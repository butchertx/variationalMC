# spin-1-vmc
C++ implementation of VMC for spin-1 and SU(3) lattice models

# Installation instructions

## Cloning the repository

Click the green "Code" button at the top right of the homepage for the repository.  Make sure the HTTPS tab is selected, copy the https address, then on your system navigate to the directory where you want to keep the code.  Type the command (in Git Bash for Windows, or terminal for Linux/Unix/Mac):

``` bash
git clone https://github.com/butchertx/spin-1-vmc.git
```

This will copy the spin-1-vmc directory into your current working directory.

## Install Intel MKL

## Link Intel libraries

## Build the spin-1-vmc library


# Git Instructions and Workflow

  If you are new to git or want to refresh yourself on the Git Workflow, see [this online tutorial](https://www.atlassian.com/git/tutorials/setting-up-a-repository).

  ## Create a new branch

  To get started working on a feature, use the command

  ``` bash
  git checkout -b <descriptive-branch-name>
  ```

  This command creates a new branch <descriptive-branch-name> (you don't have to include the angular brackets), and switches to that new branch as the current branch.

  ## Edit the code

  Make edits to the code to work on your new feature.
  
  ## Commit your changes to your development branch
  
  Use ``git add <file>`` to stage changes in a given file.  This marks them to be committed.  Then run ``git commit`` to commit your changes.  This is like pressing ``ctrl-S`` when working on a file in a text editor; it saves the changes you've made to your git branch.  Git will open a terminal editor (vim is the default) where you can type a commit message describing which changes you've made in this commit.  To exit vim, type ``:wq`` and press Enter.  You can also use the shorthand ``git commit -m "<message>"`` to create a commit message in one step.
  
  ## Useful commands
  
  Use
  
  ``` bash
  git status
  ```
  
  to view which branch you're on, and any changed files that are staged or unstaged:
  
  ![image](https://user-images.githubusercontent.com/16173143/202441913-6b84ee3f-2be5-4a42-9c47-04f218fcfe5e.png)

  The above image demonstrates a series of actions you would take to verify you're on your development branch, and add and commit any files you've changed.
  
  Another useful command is 
  
  ``` bash
  git log
  ```
  
  This will make a list of the previous commits from your current branch:
  
  ![image](https://user-images.githubusercontent.com/16173143/202442235-9a63538f-520f-4948-b4c3-7e19dbcbe1a4.png)


  ## Keep your branch up-to-date with ``main``
  
  ### Sync your branch with ``main`` with a rebase
  
  Run 

  ``` bash
  git fetch
  ``` 

  to fetch any changes to other branches/code on Github.  Then, run
  
  ``` bash
  git rebase origin/main
  ```
  
  This will combine any commits from ``main`` into your commit history.  There may be merge conflicts that you will need to sort out.  I find this is easiest to do in VS Code with the visual interface (add the Git extension).  It will highlight changes and allow you to point-and-click to keep or discard changes.
  
  When ``git rebase`` is finished, you will have a branch that can be merged without any conflicts.  Run
  
  ``` bash
  git push origin <branch-name>
  ```
  
  This will push your branch to Github.  (Note: you can and should do this regularly, even if you don't want to create a pull request yet.  This will allow other users to look at your code by pulling from Github).
  
    
  ## Create a pull request
  
  To get your branch merged back into the main branch, you will need to create a "pull request" on Github.  When you push your branch, and log in to the Github repo, there will be a banner at the top asking you to "Compare & pull request".  Click this banner.  You will come to a screen that looks like this:
  
  ![image](https://user-images.githubusercontent.com/16173143/202440722-1e4832e8-3529-4148-a5a6-e8feedb3a901.png)
  
  At the top you can select which branches you want to merge.  You can also add a title and a description, using Markdown for formatting if desired.  Then click "Create pull request".  This pull request can be viewed, edited, commented on, and eventually merged (by me, @butchertx) by clicking the "Pull requests" tab at the top of the repository site page.

