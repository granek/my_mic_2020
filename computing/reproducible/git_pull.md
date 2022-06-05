# Updating Workshop Content: Git Pull
To download changes that have been made to the wokshop content you can do a `git pull`

1. Click on the **Pull** button in the Git pane.
2. Git will put some cryptic text, but as long as it doesn't say "error", you should be all set.
2. Once the git pull is done, click on **Close** in the dialog box.

If you get an error message, try the [When Git Pull Fails](#when-git-pull-fails) section below.

# When Git Pull Fails
If you get the following error when you do a `git pull`:
*error: Your local changes to the following files would be overwritten by merge . . .*

There are a few options for solving the problem, which are given below. After doing one of the following try doing a **git pull** again.

## Copy and Revert

1. Look in the Git pane (top right) for the files that you have changed.
2. Do the following for each changed file, one at a time:
    1. Find the changed files in the Files pane (bottom right).
    2. Click the checkbox next to the file
    3. Click on **Rename** in the Files pane
    4. In the Rename dialog box, add "_MINE" to the end of filename before the ".Rmd" extension. So "fastqc.Rmd" would become "fastqc_MINE.Rmd"
    5. Click **OK**
    6. In the Git pane, the blue **M** next to the changed file should now be a red **D**. Click on the filename to highlight it.
    7. Click on the **More** menu in the Git pane and select **Revert...**
    8. In the **Revert Changes** dialog box click **Yes**
    9. The changed file should now disappear from the Git pane
3. Now you should be able to do a "git pull" by clicking the **Pull** button in the Git pane

## Stashing

1. Click on the **Tools** menu, select **Terminal** then **New Terminal**. This should open a Terminal pane in the bottom left.
2. In the Terminal type `git stash` then run the command with the "Return/Enter" key.
3. The changed file should now disappear from the Git pane. If it does not, try clicking the refresh button (circle arrow) in the top right corner of the git pane.
4. Now you should be able to do a "git pull" by clicking the **Pull** button in the Git pane

There is more information about `git stash` in the [Pro Git book](https://git-scm.com/book/en/v2/Git-Tools-Stashing-and-Cleaning)

