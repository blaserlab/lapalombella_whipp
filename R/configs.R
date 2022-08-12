theme_set(theme_cowplot(font_size = 10))

# conflicts ---------------------------------------------------------------
# resolve conflicting function names here

conflict_prefer("filter", "dplyr")
conflict_prefer("lag", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("rename", "dplyr")
conflict_prefer("count", "dplyr")


# output -------------------------------------------------
figures_out <- "~/network/X/Labs/Blaser/share/collaborators/lapalombella_whipp_network/figures"

# source local configs ----------------------------------------------------
# these are sourced after main configs and will overwrite duplicate entries if
# present. The file local_configs.R is ignored by git and so is useful for user-
# specific configurations such as output directories or formatting.

source("R/local_configs.R")
