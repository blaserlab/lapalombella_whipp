library(blaseRtools)

renv::install("blaserlab/blaseRtemplates")

blaseRtemplates::setup_git_collab()

# blasertemplates::fork_github_project(repo = <owner>/<repo>)


#restore package library from lock file
renv::restore()

# create a working branch for your day's work
blaseRtemplates::git_easy_branch(branch = "user_working")

#enter ssh
# #gitcreds::gitcreds_set()

####10.221.152.46?
####git@github.com:blaserlab/lapalombella_whipp.git
####2f42mt73GKxAxb00kn8BSMz7InNJZ+Nw2lnsczdUzTQ
####ssh-ed25519 AAAAC3NzaC1lZDI1NTE5AAAAIHS/WzMt/jZIgsc+irXi6hZXuEL4BbDGPSEv+61GzWjQ ethan.whipp@osumc.edu

usethis::edit_r_environ()

gert::git_pull()

bb_renv_datapkg("~/network/X/Labs/Blaser/collaborators/lapalombella_whipp_network/datapkg")
requireData("lapalombella.whipp.datapkg")


#use the monocle3 instructions for pseudotime