
# Remote Git Repository

## Why a Remote Repository
* Collaborate
* Remote backup
* Share/Publish

## Setup SSH Key

### Create Key in RStudio
1. Go to Global Options (from the Tools menu)
2. Click Git/SVN
3. Click **Create RSA Key**
4. Enter a passphrase (Optional)
5. Click **Create**
6. Click **View public key**
7. Select and copy all the text in the box
8. Follow these instructions to
	[Add an SSH key to your GitHub account](https://help.github.com/articles/adding-a-new-ssh-key-to-your-github-account/), starting with step #2.

## Remotes in GitHub

1. Follow  [Version Control with Git: 7. Remotes in GitHub](http://swcarpentry.github.io/git-novice/07-github/) through where it says "The next step is to connect the two repositories. We do this by making the GitHub repository a remote for the local repository. The home page of the repository on GitHub includes the string we need to identify it:"
2. In Github, under "Quick setup — if you’ve done this kind of thing before", click on "SSH"
3. Below that, copy the text in the section **. . . or push an existing repository from the command line** by clicking on the "copy to clipboard" button on the right side.
3. In RStudio, open the Terminal (in the **Tools** menu)
4. Paste the command you copied from Github.  It should look like this (except your github username should be there instead of "GITHUB_USER_NAME"):
```
git remote add origin git@github.com:GITHUB_USER_NAME/planets.git
git push -u origin master
```
If you entered a passphrase when you generated your SSH key, you will be prompted to enter it.
> The "git remote . . ." command associates your local repo with the repo you just made on Github.  The "git push . . ." command pushes everything from your local repo to the Github repo. Now the Github repo should be a perfect copy of your local repo.

5. Let's check that these two commands worked . . . go to your Github account in your webbrower and click on "planets" near the top.  You should see the three files that we commited to our local repo: `.gitignore`, `mars.txt`, and `planets.Rproj`
6. Let's add another line to `mars.txt` in RStudio: "Rhymes with cars". Then stage and commit it.
7. Notice that near the top of the **Review Changes** dialog it now says "Your branch is ahead of 'origin/master' by 1 commit.".  This is because we have commited local changes, but haven't updated our remote repo (on Github).  Let's sync the remote repo by clicking the **Push** button.  Check on Github for the changes!

Do the following challenges from the bottom of [Version Control with Git: 7. Remotes in GitHub](http://swcarpentry.github.io/git-novice/07-github/):

  * GitHub GUI
  * Push vs. Commit



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
