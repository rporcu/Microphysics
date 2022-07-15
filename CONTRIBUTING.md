## Development Model

Development generally follows the following ideas:

  * New features are merged into to the `develop` branch using
    Merge Requests (MRs).

    Nightly regression testing is used to ensure that no answers
    change (or if they do, that the changes were expected).

  * Bug fixes, questions and contributions of new features are welcome!

       * Bugs should be reported through GitLab Issues.
       * All contributions should be done via merge requests.
         A merge request should be generated from your fork of
         MFIX-Exa and target the `develop` branch. See below for
         details on how this process works.

         In general we squash commits upon merge to have a clean history.
         *Please ensure that your MR title and first post are descriptive,
         since these will be used for a squashed commit message.*

         Please note the following:
            If you choose to make contributions to the code
            then you hereby grant a non-exclusive, royalty-free perpetual license
            to install, use, modify, prepare derivative works,
            incorporate into other computer software,
            distribute, and sublicense such enhancements or derivative works
            thereof, in binary and source code form.

  * On the first workday of each month, we make a tagged release.

## Git workflow

MFIX-Exa uses [git](https://git-scm.com) for version control. If you
are new to git, you can follow one of these tutorials:
- [Learn git with bitbucket](https://www.atlassian.com/git/tutorials/learn-git-with-bitbucket-cloud)
- [git - the simple guide](http://rogerdudler.github.io/git-guide/)

### Make your own fork and create a branch on it

The basic workflow is:
- Fork the main repo (or update it if you already created it).
- Implement your changes and push them on a new branch `<branch_name>` on
your fork.
- Create a Merge Request from branch `<branch_name>` on your fork to branch
`develop` on the main MFIX-Exa repository.

First, let us setup your local git repo. To make your own fork of the main
(`upstream`) repository, press the fork button on the [MFIX GitLab page](https://mfix.netl.doe.gov/gitlab/exa/mfix).


```
git clone https://mfix.netl.doe.gov/gitlab/<myGitLabUsername>/mfix.git

# Then, navigate into your repo, add a new remote for the main MFIX-Exa repo, and fetch it:
cd mfix
git remote add upstream https://mfix.netl.doe.gov/gitlab/exa/mfix
git remote set-url --push upstream https://mfix.netl.doe.gov/gitlab/<myGitLabUsername>/mfix.git
git fetch upstream

# We recommend setting your development branch to track the upstream one instead of your fork:
git branch -u upstream/develop
```

> Note: you do not have to re-do the setup above every time.
> Instead, in the future, you need to update the `develop` branch
> on your fork with
> ```
> git checkout develop
> git pull
> ```

Make sure you are on the `develop` branch with
```
git checkout develop
```
in the MFIX-Exa directory.

Create a branch `<branch_name>` (the branch name should reflect the piece
of code you want to add, like `high_order_interpolation`) with
```
git checkout -b <branch_name>
```
and do the coding you want.
Add the files you work on to the git staging area with
```
git add <file_I_created> <and_file_I_modified>
```
### Commit & push your changes

Periodically commit your changes with
```
git commit -m "This is a 50-char description to explain my work"
```

The commit message (between quotation marks) is super important in order to
follow the developments and identify bugs.

For the moment, commits are on your local repo only. You can push them to
your fork with
```
git push --set-upstream origin <branch_name>
```

If you want to synchronize your branch with the `develop` branch (this is useful
when `develop` is being modified while you are working on
`<branch_name>`), you can use
```
git pull upstream develop
```
and fix any conflicts that may occur.

Do not merge your branch for MR into your local `develop` branch,
because it will make your local `develop` branch diverge from the
matching branch in the main repository after your MR is merged.

### Submit a Merge Request

A Merge Request is the way to efficiently visualize the changes you made
and to propose your new feature/improvement/fix to MFIX-Exa.
Right after you push changes, a banner should appear on the GitLab page of
your fork, with your `<branch_name>`.
- Click on the `Create merge request` button to prepare your MR.
- It is time to communicate your changes: write a title and a description for
your MR. People who review your MR are happy to know
  * what feature/fix you propose, and why
  * how you made it (created a new class than inherits from...)
  * and anything relevant to your MR (performance tests, images, *etc.*)
- Press `Create merge request`. Now you can navigate through your MR, which
highlights the changes you made.

Please DO NOT write large Merge Requests, as they are very difficult and
time-consuming to review. As much as possible, split them into small,
targeted MRs.
For example, if find typos in the documentation open a merge request that only fixes typos.
If you want to fix a bug, make a small merge request that only fixes a bug.
If you want to implement a large feature, write helper functionality first, test it and submit those as a first merge request.
If you want to implement a feature and are not too sure how to split it,
just open a discussion about your plans and ping other MFIX-Exa developers on it to chime in.

Even before your work is ready to merge, it can be convenient to create a MR
(so you can use GitLab tools to visualize your changes). In this case, select
the option "Mark as draft" after creating the MR.

Once your merge request is made, we will review and potentially merge it.
We recommend always creating a new branch for each merge request, as per the above instructions.
Once your merge request is merged, you can delete your local MR branch with
```
git branch -D <branch_name>
```

and you can delete the remote one on your fork with
```
git push origin --delete <branch_name>
```

Generally speaking, you want to follow the following rules.

  * Do not merge your branch into your local `develop` branch that tracks MFIX-Exa
    `develop` branch.  Otherwise your local `develop` branch will diverge from MFIX-Exa
    `develop` branch.

  * Do not commit in your `develop` branch that tracks MFIX-Exa `develop` branch.

  * Always create a new branch based off `develop` branch for each merge request, unless you are
    going to use git to fix it later.

If you have accidentally committed in `develop` branch, you can fix it as follows,
```
git checkout -b new_branch
git checkout develop
git reset HEAD~2  # Here 2 is the number of commits you have accidentally committed in development
git checkout .
```
After this, the local `develop` should be in sync with MFIX-Exa `develop` and your recent
commits have been saved in `new_branch` branch.

If for some reason your MR branch has diverged from MFIX-Exa, you can try to fix it as follows.  Before
you try it, you should back up your code in case things might go wrong.
```
git fetch upstream   # assuming upstream is the remote name for the official MFIX-Exa repo
git checkout -b xxx upstream/develop # replace xxx with whatever name you like
git branch -D develop
git checkout -b develop upstream/develop
git checkout xxx
git merge yyy  # here yyy is your MR branch with unclean history
git rebase -i upstream/develop
```
You will see something like below in your editor,
```
pick 7451d9d commit message a
pick c4c2459 commit message b
pick 6fj3g90 commit message c
```
This now requires a bit of knowledge on what those commits are, which commits have been merged,
which commits are actually new.  However, you should only see your only commits.  So it should be
easy to figure out which commits have already been merged.  Assuming the first two commits have been
merged, you can drop them by replace `pick` with `drop`,
```
drop 7451d9d commit message a
drop c4c2459 commit message b
pick 6fj3g90 commit message c
```
After saving and then exiting the editor, `git log` should show a clean history based on top of
`develop` branch.  You can also do `git diff yyy..xxx` to make sure nothing new was dropped.  If
all goes well, you can submit a MR using `xxx` branch.
Don't worry, if something goes wrong during the rebase, you an always `git rebase --abort` and start over.

## MFIX-Exa Coding Style Guide

### Code Guidelines

MFIX-Exa style guidelines follow the AMReX coding style.

Developers should adhere to the following coding guidelines:
  * Indentations should use 4 spaces, not tabs.
  * Use curly braces for single statement blocks. For example:
```cpp
       for (int n=0; n<10; ++n) {
           Print() << "Like this!";
       }
```
  or
```cpp
       for (int n=0; n<10; ++n) { Print() << "Like this!"; }
```
  but not
```cpp

       for (int n=0; n<10; ++n) Print() << "Not like this.";
```
  or
```cpp
       for (int n=0; n<10; ++n)
          Print() << "Not like this.";
```
  * Add a space after the function name and before the
parenthesis of the parameter list (but
not when simply calling the function). For example:
```cpp
        void CorrectFunctionDec (int input)
```
  Not
```cpp
        void IncorrectFunctionDec(int input)
```
  This makes it easy to find where functions are defined with grep.
  * Member variables should be prefixed with `m_`. For example:
```cpp
       amrex::Real m_variable;
```
These guidelines should be adhered to in new contributions, but
please refrain from making stylistic changes to unrelated sections of code in your MRs.



## Core Developers

People who make a number of substantive contributions will be named
"core developers" of MFIX-Exa.  The criteria for becoming a core
developer are flexible, but generally involve one of the following:

  * 100 non-trivial commits to `mfix/src/`  *and/or*

  * addition of a new algorithm / module  *and/or*

  * substantial input into the code design process or testing

If a core developer is inactive for multiple years, we may reassess their
status as a core developer.

The current list of core developers is:
Ann Almgren (LBNL),
John Bell (LBNL),
William Fullmer (NETL),
Jordan Musser (NETL),
Andrew Myers (LBNL),
Roberto Porcu (NETL),
Deepak Rangarajan (NETL),
Michele Rosso (LBNL),
Weiqun Zhang (LBNL)
