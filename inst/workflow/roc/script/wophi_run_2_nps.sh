#!/bin/sh

echo "wophi_run_2_nps.sh"

### Check if configuration file exists.
if [ -e 00-set_env.r ]; then
  echo "00-set_env.r found."
else
  echo "00-set_env.r does not exist."
  exit
fi

### Root of all outputs.
ALL_OUT=`Rscript -e 'source("00-set_env.r");cat(prefix$all.out)'`

### Codes.
CODE_PATH=`Rscript -e 'source("00-set_env.r");cat(prefix$code)'`
CODE_PLOT_PATH=`Rscript -e 'source("00-set_env.r");cat(prefix$code.plot.nps)'`

### Plotting fitted results.
Rscript ${CODE_PLOT_PATH}/06-plotsingle_model.r > \
          ${ALL_OUT}/log_nps/06-plotsingle_model 2>&1 &
Rscript ${CODE_PLOT_PATH}/06-plotsingle_prxy.r > \
          ${ALL_OUT}/log_nps/06-plotsingle_prxy 2>&1 &
Rscript ${CODE_PLOT_PATH}/06-plotsingle_prxy_wci.r > \
          ${ALL_OUT}/log_nps/06-plotsingle_prxy_wci 2>&1 &

### Plotting diagnoses.
Rscript ${CODE_PLOT_PATH}/07-plotdiag_acceptvsEPhi.r > \
          ${ALL_OUT}/log_nps/07-plotdiag_acceptvsEPhi 2>&1 &
Rscript ${CODE_PLOT_PATH}/07-plotdiag_bin.r > \
          ${ALL_OUT}/log_nps/07-plotdiag_bin 2>&1 &
Rscript ${CODE_PLOT_PATH}/07-plotdiag_EPhi_hist.r > \
          ${ALL_OUT}/log_nps/07-plotdiag_EPhi_hist 2>&1 &

### Plotting traces.
Rscript ${CODE_PLOT_PATH}/07-plottrace_quantile_Phi.r > \
          ${ALL_OUT}/log_nps/07-plottrace_quantile_Phi 2>&1 &
