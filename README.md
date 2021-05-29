
# lapalombella_whipp

This repository holds the analysis scripts for the collaborative project with the Lapalombella Lab looking at Richter's transformation from CLL.

The companion data repository is at (https://github.com/blaserlab/lapalombella.whipp.datapkg).  Due to file size limitations, the .rda objects are stored using git-lfs.  This means that the data repo cannot be installed directly.  The data are available as precompiled binaries "<package name>.tar.gz", which you can install from R studio.

The best approach is to install interactively using the "packages" tab.  Navigate to "~/network/X/Labs/Blaser/single_cell/lapalombella_whipp/datapkg" and install the latest version.

Check "NEWS.md" for information about new versions of the data.

Installing the datapackage should install the package "lazyData" as a dependency.  This provides functionality to "lazy load" the data objects.  This means they will reside in a hidden environment in your R session and be available when called in your script.  

In order to activate the hidden environment you must run the following once and only once:

```{r}
requireData("lapalombella.whipp.datapkg")
```
If you run it again it will unload the data.  Do this only if you want to free up memory.

The best approach is to source the file "R/dependencies.R" when you start the session.  This will load the data using the command above and also load several libraries we always are using.

The main data object you want to work with is "cds_main".  This stores aligned UMAP coordinates internally, but prealignment coordinates are stored in metadata columns if you wish to use them.

If you make modifications to data objects, they will appear in your global environment and will supercede the stored version when you subsequently use them.  For lightweight modifications like adding convenience metadata columns, no changes to the stored data are required.  If we want to incorporate them in the future we can do that.  If you develop something that takes longer to compute, we can talk about adding it as a new data object to the data package.  Generally speaking, there should be no reason to save data objects (.rda, .RData, .rds) in this repo unless we are going to immediately move them over to the data package.

R scripts should be in the "R" directory.  If you want to output graphics or stats, you can create "plots" or "stats" folders in the root directory of the project.  This is good for scrolling through preliminary outputs.  Once you start putting figures together for the paper, you can write to an external directory. 




