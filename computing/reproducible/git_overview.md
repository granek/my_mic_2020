# Why Version Control
* [R and version control for the solo data analyst](https://stackoverflow.com/questions/2712421/r-and-version-control-for-the-solo-data-analyst)

## Why Git
* Popular
* Powerful
* Distributed
* Free
* Open Source

# Using Git in RStudio
RStudio provides a nice GUI for using git, which can reduce the
learning curve for git.
	
## Basic Setup in RStudio
1. Go to Global Options (from the Tools menu)
2. Click Git/SVN
3. Click Enable version control interface for RStudio projects

## Intro to Git in RStudio

### Create a Project with Version Control

1. Select the "New Project" command from the Project menu
2. Select the "New Directory" option in the "Create Project" dialog
box
3. Select "Empty Project" for "Project Type"
4. In the "Create New Project" dialog do the following:
   * Directory name: "planets"
   * Create project as subdirectory of: "~"
   * Be sure there is a check in the box for "Create a git repository"
   * Click "Create Project"
5. Switch to the "Git" tab in the upper right pane.
6. Check the "Staged" boxes next to .gitignore and planets.Rproj, then
click "Commit"
7. In the **Commit message** box, type "First commit!"
8. Oops, we have gotten our first git message - it wants to know your
name and email because this is incorporated in the commit history!
9. Copy the following text from the error message and paste it into a
new text file:
```
git config --global user.email "you@example.com"
git config --global user.name "Your Name"
```
10. Replace "you@example.com" and "Your Name" with your information,
then copy the revised text.
11. Open the RStudio Terminal (under the **Tools** menu)
12. Paste the revised text at the command prompt, then hit the
**return** key on your keyboard
13. You can confirm that you successfully configured your name and email
with this command: `git config --list`
14. Now go back to the **Review Changes** dialog and close the **Git
Commit** message.
15. Click **Commit** again, then close the **Git Commit** and **Review
    Changes** dialog boxes.

### Tracking Changes

#### Adding a New File
1. Make a new text file in the **planets** project.
2. Add some text: "Cold and dry, but everything is my favorite color"
3. Save as `mars.txt`.
4. Notice that `mars.txt` is now in the Git pane.  The two question
   marks under **Status** mean that git doesn't know anything about
   this file (and is not responsible for it)
5. Click the checkbox under **Staged** to add the file to the repo,
notice that **Status** changes to "A" (for "added").
6. Now click **Commit** to open the **Review Changes** dialog, put in the
   commit message "Start notes on Mars as a base", commit the
   file, and close the dialog boxes.
7. To see our project history, click on the clock icon in the Git pane.

#### Making changes
1. Add a second line to `mars.txt`: "The two moons may be a problem
for Wolfman" and save.
2. Notice that `mars.txt` pops up again in the **Git** pane with an
**M** for "modified".
3. Click on the **Diff** button to open the **Review Changes**
   dialog. It shows us the line we added highlighted in green, and
   shows us the unchanged previous line in gray.
4. To commit these changes to git, add the commit message "Add concerns about effects of Mars' moons
on Wolfman", click **Commit**. Notice the error message from git:

```
Changes not staged for commit:
	modified:   mars.txt

no changes added to commit
```

Git messages are usually very informative! Git is telling us that we
can't commit because we didn't *add* anything to the commit!  We
forgot to click the **Staged** button.

So, now click on the **Staged** checkbox, be sure our commit message
is still there,  click **Commit**, check the message from git, then close the dialogs.

#### Staging
Git insists that we add files to the set we want to commit before actually committing anything. This allows us to commit our changes in stages and capture changes in logical portions rather than only large batches. For example, suppose we’re adding a few citations to our supervisor’s work to our thesis. We might want to commit those additions, and the corresponding addition to the bibliography, but not commit the work we’re doing on the conclusion (which we haven’t finished yet).

To allow for this, Git has a special staging area where it keeps track of things that have been added to the current changeset but not yet committed.

> ##### Staging Area
> If you think of Git as taking snapshots of changes over the life of a project, `git add` specifies what will go in a snapshot (putting things in the staging area), and `git commit` then actually takes the snapshot, and makes a permanent record of it (as a commit). If you don’t have anything staged when you type git commit, Git will prompt you to use `git commit -a` or `git commit --all`, which is kind of like gathering everyone for the picture! However, it’s almost always better to explicitly add things to the staging area, because you might commit changes you forgot you made. (Going back to snapshots, you might get the extra with incomplete makeup walking on the stage for the snapshot because you used -a!) Try to stage things manually, or you might find yourself searching for “git undo commit” more than you would like!

Let's examine this process more carefully . . .

1. Add a third line to `mars.txt`: "But the Mummy will appreciate the lack of humidity" and save.
2. Open the **Review Changes** dialog.
3. Notice that the **Show: Unstaged** button is checked. We can toggle
   between viewing staged changes and unstaged changes by clicking on
   the buttons next to **Show: Staged** and **Show: Unstaged**.  Right now no
   changes are staged.
4. Click on the **Staged** button next to `mars.txt` to stage the
   changes.  Notice that **Show: Staged** is automatically selected,
   but we can flip back to **Show: Unstaged** to see that there are no
   longer any unstaged changes.
5. Add the commit message "Discuss concerns about Mars' climate for Mummy", click **Commit**, check the message from git, then close the dialogs.
6. Now we can check the history and see the changes made at each
commit and the commit message.

#### Challenge: Choosing a Commit Message
Which of the following commit messages would be most appropriate for the last commit made to mars.txt?

1. “Changes”
2. “Added line ‘But the Mummy will appreciate the lack of humidity’ to
mars.txt”
3. “Discuss effects of Mars’ climate on the Mummy”

### Reverting Changes

"If I could turn back time . . ." -- Cher

One of the key benefits of using version control is the ability to do
something not usually possible in life - going back in time.  Lot's of
software has **Undo**, but when version control is used well, it is
like an infinite super-duper undo with a safety net.

#### Basic Reversion
1. Add a forth line to `mars.txt`: "An ill-considered change" and save.
2. Open the **Review Changes** dialog.
3. If we want to get rid of this change, we can click the **Revert**
button, then click **Yes** . . . voila!

#### Revert Some Changes, but not All
1. Add a these two lines `mars.txt`:
```
No wind
A little bit dusty
```
2. Open the **Review Changes** dialog.
3. Oops, it turns out that there **is** wind on mars, so we need to
revert line 4, but commit line 5!
4. Click on "No wind" to select it, then click **Discard line**.
5. Now click on **Stage All** to add the remaining line, add a commit
   message and commit

#### Going Further Back
Git is super powerful!  It will let you do all sorts of things.  For example, you can revert some or all files back to any previous commit . . . but you cannot access all of git's features with the RStudio interface.  To do fancier things you have to use the command line interface.  Within RStudio you can do this using the **Terminal** under the **Tools** menu.  The Software Carpentry module [Version Control with Git](http://swcarpentry.github.io/git-novice/) covers this, and google will help you figure out how to do just about anything with git.

### Ignoring Things

What if we have files that we do not want Git to track for us, like backup files created by our editor or intermediate files created during data analysis. Let’s create a few dummy files! Make new text files with the following names (you can leave them empty): 
```
a.dat
b.dat
c.dat
results/a.out
results/b.out
```

Putting these files under version control would be a waste of disk space. What’s worse, having them all listed could distract us from changes that actually matter, so let’s tell Git to ignore them.

RStudio makes this pretty easy to do:

1. Click on the `a.dat` file in the git pane to select it.
2. Click on the gear in the git pane (it might say "More" next to it), and select **Ignore**.  This will open a **Git ignore** dialog box with a list of files that git is ignoring, with `a.dat` added to the list (the others were automatically added by RStudio when we made this project).
3. Click on **save** to accept the addition of `a.dat` to the list. 
> Notice that `a.dat` has now disappeared from the git pane, but `.gitignore` has just appeared as modified.  `.gitignore` is the file where we specify what git should ignore - that **Git ignore** dialog box was actually a mini text editor, and by adding `a.dat` we changed the file (you may remember that the first thing we did after creating the **planets** project was to commit `.gitignore` to the repo).

4. What about the other `.dat` files? Select `b.dat` in the git pane then open the **Git ignore** dialog.  Now replace "a.dat" and "b.dat" with "*.dat", this will tell git to ignore *any* file that ends in ".dat"
5. Let's also add `results` directory, but it is a good idea to add it as "results/", the forward slash at the end means that we are specifically refering to a directory (and its contents)
6. Now lets stage and commit our changes to `.gitignore`
> Sometimes you will want to include `.gitignore` in the repo, and sometimes you won't.  In cases where you don't want to include it in the repo, you will want to add ".gitignore" to the list of files to ignore in `.gitignore`, which is very meta!

Now let's work through the challenges at the end of [Version Control with Git: 6. Ignoring Things](https://swcarpentry.github.io/git-novice/06-ignore/index.html)


# References
The bulk of this set of lessons is a translation from Unix command line to RStudio GUI of the The Software Carpentry module [Version Control with Git](http://swcarpentry.github.io/git-novice/), specifically:

- [Create a Project with Version Control](#create-a-project-with-version-control) is based on [Version Control with Git: 3. Creating a Repository](http://swcarpentry.github.io/git-novice/03-create/)
- [Tracking Changes](#tracking-changes) is based on: [Version Control with Git: 4. Tracking Changes](http://swcarpentry.github.io/git-novice/04-changes/)
- [Reverting Changes](#reverting-changes) is loosely based on: [Version Control with Git: 5. Exploring History](http://swcarpentry.github.io/git-novice/05-history/)
- [Ignoring Things](#ignoring-things) is based on [Version Control with Git: 6. Ignoring Things](http://swcarpentry.github.io/git-novice/06-ignore/)
- [Remotes in GitHub](#remotes-in-github) is based on [Version Control with Git: 7. Remotes in GitHub](http://swcarpentry.github.io/git-novice/07-github/)
- [Collaborating](#collaborating) is based on [Version Control with Git: 8. Collaborating](http://swcarpentry.github.io/git-novice/08-collab/) and [Version Control with Git and SVN](https://support.rstudio.com/hc/en-us/articles/200532077-Version-Control-with-Git-and-SVN)
- [Conflicts](#conflicts) is based on [Version Control with Git: 8. Conflicts](http://swcarpentry.github.io/git-novice/09-conflict/)

[Collaborating](#collaborating) is also based on  [Version Control with Git and SVN](https://support.rstudio.com/hc/en-us/articles/200532077-Version-Control-with-Git-and-SVN).

The [Team Conflicts](#team-conflicts) section is based on [We'll git there, slowly but surely](https://github.com/mine-cetinkaya-rundel/github-class-sigcse-2018).
