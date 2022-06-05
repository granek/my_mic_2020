# Cloning
Cloning a repository creates a local copy of a repository.  Usually one clones from a remote repository (e.g. GitHub, Gitlab, Bitbucket).  You can clone a remote repository that you own, a remote repository that you are contributing to, and any other remote repository that is public.  You generally clone remote repositories that you own or are contributing to so that you you can work on them - adding material, fixing bugs, or otherwise updating.  A common reason to clone public remote repositories is to install software that is published as a repository.

## How to Clone the Intern Repository:
1. In a web browser, go to the [repository on GitLab](https://gitlab.oit.duke.edu/mic-course/2022-mic).
2. Click on the blue **Clone** button in the top right part of the repository page.
3. A dialog box will pop up that says "Clone with SSH", "Clone with HTTPS", etc. 
4. In this case you will want "Clone with HTTPS", so click on the clipboard icon to copy the URL.
5. In RStudio click on the **File** menu, then select **New Project**.
4. Choose **Version Control** in the *Create Project* dialog box.
5. Choose **Git** in the *Create Project from Version Control* dialog box.
6. Paste the repository URL you just copied into the *Repository URL* box.
7. For *Project directory name* you will generally want to use the repo name.
8. For *Create project as subdirectory of:* click on the **Browse** button
9. In the *Choose Directory* dialog that opens, click on **Home** then **New Folder**
10. In the *New Folder* dialog box type "workshop", then click **OK**.
11. Now in the *Choose Directory* dialog click **Choose**
12. Now that you are back in the *Clone Git Repository* dialog, click **Create Project**.

## Avoiding Merge Conflicts
To avoid merge conflicts, if you want to edit a file in your fork I recommend that you copy the file, and either rename the copy so that you can distinguish if from the original, or keep a separate directory for the copies that you edit.

### Rename the copy
If you wanted to edit this file, you could rename the original `git_cloning_and_forking.Rmd` to `git_cloning_and_forking_MINE.Rmd`.

### Directory for editted copies
You could make a new directory in your fork called `my_edits` and put a copy of this file in `my_edits`.


# References
This lesson borrows heavily from [Github Help](https://help.github.com/en) and [GitHub Community Forum](https://github.community),  specifically:
- [Cloning a repository](https://help.github.com/en/articles/cloning-a-repository)

## Additional References
- [Pro Git: GitHub - Contributing to a Project](https://git-scm.com/book/en/v2/GitHub-Contributing-to-a-Project)
- [Pro Git: Cloning an Existing Repository](https://git-scm.com/book/en/v2/Git-Basics-Getting-a-Git-Repository#_git_cloning)
