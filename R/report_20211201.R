#' ---
#' title: "Followup on Discussion"
#' author: "Brad Blaser"
#' date: "Dec 12 2021"
#' output: pdf_document
#' ---
#'
#+ include=FALSE
knitr::opts_chunk$set(warning = FALSE, message = FALSE)
source("R/dependencies.R")
source("R/configs.R")
#' ## Update
#' I did some "modernization" to the repository.  You should pull changes or re-clone and reapply your changes if that is easier.
#' 
#' We are now using a package manager called renv.  This gives us control and reproducibility regarding the software packages we use.  Otherwise if everything was operating from a central package library on our computer, if I update something, it may break something else in your analysis.  We want to try to avoid that.  So what we have is a central library on the server to pull packages from.  Renv copies those packages to a directory within your copy of the analysis project.  This is at "renv/library/R-4.1/x86_64-pc-linux-gnu/".  
#' 
#' When you first clone this analysis that directory will not exist or will be empty.  Renv also generates a lockfile which is just a text-based document listing all of the packages in use. **When you first set up the analysis project, uncomment ```renv::hydrate()``` in dependencies.R and run it**.  This will go to the appropriate locations (local server or remote repositories) and install the versions you need.  It may take a while.
#' 
#' Once hydrate is complete, you should uncomment and run all of the renv::install() commands.  This will install a bunch of things that usually get missed by hydrate.  Then you should be able to run the library commands.
#' 
#' The last two commands on dependencies.R will install and if necessary update the prepackaged data and then load it into the session.
#' 
#' ## Looking at Size Factors
#' 
#' You can find the size factors in the cell metadata for the cds_main object.  This you can get by running ```SingleCellExperiment::colData(cds_main)``` or ```bb_cellmeta(cds_main)``` which is a simple wrapper around colData that returns a tibble which is often easier to work with.  Like all "bb_" functions it should be loaded with blaseRtools.
#' 
#+ echo = FALSE
pander(head(bb_cellmeta(cds_main))[1:3], caption = "Cell Metadata")
#'
#' We can then plot the distribution of size factors by sample
#+ echo = TRUE, dev = "png", fig.align = "center", fig.width = 6.5, fig.height = 3
ggplot(data = bb_cellmeta(cds_main), 
       mapping = aes(x = Size_Factor, color = sample)) +
  geom_density()

#' It looks like there are definitely differences in the size factor distribution but I'm not sure if this is unexpectedly bad.
#' 
#' We can also look at the qc results from the individual samples:
#' 
#+ echo = TRUE, dev = "png", fig.align = "center", fig.height = 3, fig.width = 3.5
map(.x = bb_cellmeta(cds_main) %>% pull(sample) %>% unique(),
     .f = ~ind_qc_res[[.x]][[2]])

#' Here you can see that L37 has a bimodal distribution which also seems like it is reflected in the size factors above.  
#' 
#' As we discussed, these are things we normalize for before making dimension reductions or doing expression analysis.  
#' 
#' Let me know if you have any questions.
#' 
#' Brad

