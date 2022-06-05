Git will often get upset if you try to `git pull` after modifying files in the git repository. The first option below prevents this from happening. The second option below will allow you to fix things if you forgot to do the second option before taking notes.



# Before you start taking notes
When a session starts, if you want to take notes in an Rmd notebook, do the following before you do anything else.

1. Find the Rmarkdown file in the RStudio **Files** pane.
2. Click the checkbox next to the file.
3. Click the **More** button near the top of the **Files** pane and select **Copy**
4. By default, RStudio will put "CopyOf" at the beginning of the duplicate file's name. Click **OK** in the dialog box that pops up to accept this.
5. Open the "CopyOf" file and take notes there.

# If you took notes in the original Rmarkdown file
If you have started to modify the original Rmarkdown file, you can do the following, then you should be able to `git pull`
## Rename and Restore
1. Rename the Rmarkdown file: In the **Files** pane, click the checkbox next to the file, click **Rename**, then add "CopyOf" at the beginnning of the file name (so "myfile.Rmd" becomes "CopyOfmyfile.Rmd")
2. Try `git pull` again

## Git Stash
1. In the terminal, run `git stash` to stash changes
2. Try `git pull` again
