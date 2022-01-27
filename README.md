
# lapalombella_whipp

This repository holds the analysis scripts for the collaborative project with the Lapalombella Lab looking at Richter's transformation from CLL.

The best approach is to source the file "R/dependencies.R" and then "R/configs.R" when you start the session.  This will load the data using the command above and also load several libraries we always are using.

The main data object you want to work with is "cds_main".  This stores aligned UMAP coordinates internally, but prealignment coordinates are stored in metadata columns if you wish to use them.

R scripts should be in the "R" directory.  If you want to output graphics or stats, you can create "plots" or "stats" folders in the root directory of the project.  This is good for scrolling through preliminary outputs.  Once you start putting figures together for the paper, you can write to an external directory. 

### Modifed workflow for reducing GIT cognitive burden

The repository includes a new **scratch** directory.  This directory is not tracked by git so anything in there is kept local.  Don't edit the .gitignore file.  This is where you can do work etc without worrying about causing git conflicts.  *Since this is not tracked by git, it won't be archived and if you delete it or save over it, it will be hard or impossible to retrieve*.  

We will keep the files in the R directory as the current "best available work".  If you have something you would like to add to that, then you can email the R script to me and I will review and incorporate.  Or you can make it a new, unique file (something like ethan_work_mmddyy.R) and push it to github.  That way I can pull it down, review and add into the best available work.  Once it makes it out of the scratch directory it will be tracked by git.  Making the files unique prevents either of us from having to deal with merge conflicts.

If I generate code for you like during our meetings, I will label it uniquely and put it in the R directory so you can pull it from github.

The ultra low-drag way for you to pull and push updates is like so:

* go to the terminal panel
* Check to make sure you are in the project directory.  If not, use the dropdown menu to start a new terminal.
* run ```git add .``` to have git index everything you have changed outside of the ignored files
* run ```git stash``` to return the repository to it's last committed state
* run ```git pull``` to pull any changes I may have made
* if you want to commit and push something then move it out of your scratch directory into the R directory
* run ```git add .``` to add everything new or ```git add <insert filename here>``` to add a single file
* run ```git commit -m "<add a concise message here explaining what you did>"```
* run ```git push```
* let me know by email that you have done so

That should do it.

There are more complex ways of using git that may be more "appropriate", but I think this should keep us on the same page without unneccessary duplication of files etc.  Let me know if you have any questions.





