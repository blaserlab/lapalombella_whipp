source("renv/activate.R")
renv::settings$snapshot.type("all")
options(renv.config.auto.snapshot = TRUE)
options(renv.config.install.verbose = TRUE)

# set default git protocol to https
options(usethis.protocol  = "https")

# the following line sets the default editor to nano
# if you want to use something else like vim, put it in your ~/.Rprofile
options(editor="nano")

# remove the automatic biostrings coloring which is really slow
options(Biostrings.coloring = FALSE)

# make the R prompt show the active git branch
if (prompt::is_git_dir()) prompt::set_prompt(paste0("[ ", gert::git_branch(), " ] > "))

if (file.exists("~/.Rprofile")) source("~/.Rprofile")

