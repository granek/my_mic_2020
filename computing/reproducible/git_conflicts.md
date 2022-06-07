# Conflicts in Git

As soon as people can work in parallel, it’s likely someone’s going to step on someone else’s toes. This will even happen with a single person: if we are working on a piece of software on both our laptop and a server in the lab, we could make different changes to each copy. Version control helps us manage these conflicts by giving us tools to resolve overlapping changes.

To see how we can resolve conflicts, we must first create one.

## Creating a Conflict
Continue working with the same partner you collaborated with above.

1. In the pairs, person **A** should add a line to the `mars.txt` file in **B's** repo: "This line added to Wolfman's copy".
2. Save, commit, and push.
3. Now **B** should add the following line to the `mars.txt` file *in their own repo*: "We added a different line in the other copy", then save, commit, and push. **Note: B** should not pull before making the edit.

When **B** tries to push, git should give an error!
```
 ! [rejected]        master -> master (fetch first)
error: failed to push some refs to 'git@gitlab.oit.duke.edu:mic-course/singularity-example.git'
hint: Updates were rejected because the remote contains work that you do
hint: not have locally. This is usually caused by another repository pushing
hint: to the same ref. You may want to first integrate the remote changes
hint: (e.g., 'git pull ...') before pushing again.
hint: See the 'Note about fast-forwards' in 'git push --help' for details.
```

Git detects that the changes made in one copy overlap with those made in the other and stops us from trampling on our previous work. What we have to do is pull the changes from GitHub, merge them into the copy we’re currently working in, and then push that. 

## Merging
1. **B** should do a pull. Git will give us a message that there’s a conflict, and mark that conflict in the affected file.
```
Auto-merging README.md
CONFLICT (content): Merge conflict in README.md
Automatic merge failed; fix conflicts and then commit the result.
```
2. **B** should close and reopen **mars.txt** in the Jupyter Text Editor so that these changes are loaded.
3. When **B** looks at `mars.txt` in the editor, there should be markings from git indicating where the conflict is:
**B's** change is preceded by `<<<<<<< HEAD`. Git has then inserted `=======` as a separator between the conflicting changes and marked the end of the content downloaded from GitHub with `>>>>>>>`. (The string of letters and digits after that marker is a unique identifier for the commit we’ve just pulled.)
It is now up to us to edit this file to remove these markers and reconcile the changes. We can do anything we want: keep the change made in the local repository, keep the change made in the remote repository, write something new to replace both, or get rid of the change entirely. 
4. Let’s replace both of the added lines, and the git markings, and add this as the last line: "We removed the conflict on this line". so now `mars.txt` should look something like:
```
Cold and dry, but everything is my favorite color
The two moons may be a problem for Wolfman
But the Mummy will appreciate the lack of humidity
A little bit dusty
Rhymes with cars
We removed the conflict on this line
```
5. Once that is saved `git status` should report:
```
Your branch and 'origin/maing' have diverged,
and have 1 and 1 different commits each, respectively.
  (use "git pull" to merge the remote branch into yours)

All conflicts fixed but you are still merging.
  (use "git commit" to conclude merge)

Changes to be committed:
        modified:   mars.txt
```

Note the important guidance "All conflicts fixed but you are still merging."

6. Stage and commit as normal. In your commit message you might want to indicate that you did a merge, something like: "mars.txt: replaced conflicting lines between A and B".

7. Push your commit

Git keeps track of what we’ve merged with what, so we don’t have to fix things by hand again when the collaborator who made the first change pulls again


We don’t need to merge again because Git knows someone has already done that.

Git’s ability to resolve conflicts is very useful, but conflict resolution costs time and effort, and can introduce errors if conflicts are not resolved correctly. If you find yourself resolving a lot of conflicts in a project, consider these technical approaches to reducing them:

* Pull from upstream more frequently, especially before starting new work
* Use topic branches to segregate work, merging to master when complete
* Make smaller more atomic commits
* Where logically appropriate, break large files into smaller ones so that it is less likely that two authors will alter the same file simultaneously

Conflicts can also be minimized with project management strategies:

* Clarify who is responsible for what areas with your collaborators
* Discuss what order tasks should be carried out in with your collaborators so that tasks expected to change the same lines won’t be worked on simultaneously
* If the conflicts are stylistic churn (e.g. tabs vs. spaces), establish a project convention that is governing and use code style tools (e.g. htmltidy, perltidy, rubocop, etc.) to enforce, if necessary

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
