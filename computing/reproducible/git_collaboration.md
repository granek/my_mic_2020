# Collaborating

## Cloning someone else's repo
For the next step, get into pairs. One person will be the “Owner” and the other will be the “Collaborator”. The goal is that the Collaborator add changes into the Owner’s repository. We will switch roles at the end, so both persons will play Owner and Collaborator.

### Github Setup
1. The Owner needs to give the Collaborator access. The owner needs to find their repository on GitHub, click the settings button on the right, then select Collaborators, and enter your partner’s username.
![Github Collaborator Settings](http://swcarpentry.github.io/git-novice/fig/github-add-collaborators.png)
2. To accept access to the Owner’s repo, the Collaborator needs to go to https://github.com/notifications. Once there she can accept access to the Owner’s repo. If nothing shows up there, the collaborator should check their email to see if they received an email notification.

### Cloning the Repository
 Next, the Collaborator needs to download a copy of the Owner’s repository to her machine. This is called “cloning a repo”.

1. To get the URL for the Owner's repository, the Collaborator should go to the Owner’s repository on Github and click the green *Clone or download* button.  In the dialog box that pops up, be sure it says "Clone with SSH", then copy the Repository URL (it should be **git@github.com:OWNER_GITHUB_USERNAME/planets.git**, but with the Owner's Github username instead of "OWNER_GITHUB_USERNAME")
2. Do `cd ~` to move to your home directory.
3. To clone run `git clone URL ~/planets_shared_1`, replacing **URL** with the Repository URL you just copied. The last part, **~/planet_shared**, tells git to clone the repo into a directory with this name. If you don't supply a name git  will just use the name of the repository as it is on GitHub (planets), but this would create problems because you already have a directory named planets.

## Modifying someone else's repo
Now the Collaborator (i.e. not the owner) can now make a change in the clone of the Owner’s repository, exactly the same way as we’ve been doing before, but make sure you do this in the **planets_shared_1** directory.

1. Make a new text file and put in the text "It IS a planet!" and save as `pluto.txt` 
2. Stage `pluto.txt`, and commit it with the message "Add notes about Pluto"
3. Push the change to the Owner’s repository on GitHub by running `git push`
> Note that we didn’t have to create a remote called origin: Git uses this name by default when we clone a repository. (This is why origin was a sensible choice earlier when we were setting up remotes by hand.)

4. Take a look to the Owner’s repository on its GitHub website now (maybe you need to refresh your browser.) You should be able to see the new commit made by the Collaborator.
5. To download the Collaborator’s changes from GitHub, the Owner should do `git pull` in their own terminal. Now the three repositories (Owner’s local, Collaborator’s local, and Owner’s on GitHub) are back in sync.
6. Switch roles, whoever was Owner first time around should repeat this as Collaborator, and vice versa. Just to reduce confusion the person who is now collaborator can call the cloned repository **~/planets_shared_2**

### > A Basic Collaborative Workflow
> In practice, it is good to be sure that you have an updated version of the repository you are collaborating on, so you should git pull before making our changes. The basic collaborative workflow would be:
> 
> 1. update your local repo with git pull origin master
> 2. make your changes and stage them
> 3. commit your changes
> 4. upload the changes to GitHub with git push

> It is better to make many commits with smaller changes rather than of one commit with massive changes: small commits are easier to read and review.


# References
The bulk of this set of lessons is a translation from Unix command line to RStudio GUI of the The Software Carpentry module [Version Control with Git](http://swcarpentry.github.io/git-novice/), specifically:

- [Collaborating](#collaborating) is based on [Version Control with Git: 8. Collaborating](http://swcarpentry.github.io/git-novice/08-collab/) and [Version Control with Git and SVN](https://support.rstudio.com/hc/en-us/articles/200532077-Version-Control-with-Git-and-SVN)

[Collaborating](#collaborating) is also based on  [Version Control with Git and SVN](https://support.rstudio.com/hc/en-us/articles/200532077-Version-Control-with-Git-and-SVN).
