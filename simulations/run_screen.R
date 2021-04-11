## in taupo
R310
.libPaths(c("/home/Shared/Rlib/release-3.10-lib", 
            "/usr/local/R/R-3.6.1/library", 
            "/home/fiona/tmp/R310_lib/"))
setwd("correlationtree_analysis/")
source("simulations/non_parametric/simu_np.R")
source("simulations/parametric/simu_p.R")
source("simulations/non_parametric/simu_np_hacking.R")
