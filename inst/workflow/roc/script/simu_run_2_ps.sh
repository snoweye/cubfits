#!/bin/sh

echo "simu_run_2_ps.sh"

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
CODE_PLOT_PATH=`Rscript -e 'source("00-set_env.r");cat(prefix$code.plot.ps)'`

### Dump tsv files.
Rscript ${CODE_PLOT_PATH}/05-subset_tsv.r > \
          ${ALL_OUT}/log_ps/05-subset_tsv 2>&1 &

### Plot fitted results.
Rscript ${CODE_PLOT_PATH}/s6-plotsingle_model_true.r > \
          ${ALL_OUT}/log_ps/s6-plotsingle_model_true 2>&1 &
Rscript ${CODE_PLOT_PATH}/06-plotsingle_model.r > \
          ${ALL_OUT}/log_ps/06-plotsingle_model 2>&1 &
Rscript ${CODE_PLOT_PATH}/06-plotsingle_prxy.r > \
          ${ALL_OUT}/log_ps/06-plotsingle_prxy 2>&1 &
Rscript ${CODE_PLOT_PATH}/06-plotsingle_prxy_wci.r > \
          ${ALL_OUT}/log_ps/06-plotsingle_prxy_wci 2>&1 &

### Plot diagnoses.
Rscript ${CODE_PLOT_PATH}/07-plotdiag_acceptvsEPhi.r > \
          ${ALL_OUT}/log_ps/07-plotdiag_acceptvsEPhi 2>&1 &
Rscript ${CODE_PLOT_PATH}/07-plotdiag_bin.r > \
          ${ALL_OUT}/log_ps/07-plotdiag_bin 2>&1 &
Rscript ${CODE_PLOT_PATH}/07-plotdiag_EPhi_hist.r > \
          ${ALL_OUT}/log_ps/07-plotdiag_EPhi_hist 2>&1 &

### Plot traces.
Rscript ${CODE_PLOT_PATH}/07-plottrace_quantile_Phi.r > \
          ${ALL_OUT}/log_ps/07-plottrace_quantile_Phi 2>&1 &

### Plot for fake data only.
Rscript ${CODE_PLOT_PATH}/s6-plotsingle_b_corr_true.r > \
          ${ALL_OUT}/log_ps/s6-plotsingle_b_corr_true 2>&1 &
Rscript ${CODE_PLOT_PATH}/s6-plotsingle_prxy_true.r > \
          ${ALL_OUT}/log_ps/s6-plotsingle_prxy_true 2>&1 &
Rscript ${CODE_PLOT_PATH}/s6-plotsingle_prxy_true_wci.r > \
          ${ALL_OUT}/log_ps/s6-plotsingle_prxy_true_wci 2>&1 &
Rscript ${CODE_PLOT_PATH}/s6-plotsingle_scu_mscu.r > \
          ${ALL_OUT}/log_ps/s6-plotsingle_scu_mscu 2>&1 &

### Plot for matched cases only.
Rscript ${CODE_PLOT_PATH}/m6-plot_b_corr.r > \
          ${ALL_OUT}/log_ps/m6-plot_b_corr 2>&1 &
Rscript ${CODE_PLOT_PATH}/m6-plot_b_corr_negsel.r > \
          ${ALL_OUT}/log_ps/m6-plot_b_corr_negsel 2>&1 &
Rscript ${CODE_PLOT_PATH}/m6-plot_prxy.r > \
          ${ALL_OUT}/log_ps/m6-plot_prxy 2>&1 &
Rscript ${CODE_PLOT_PATH}/m6-plot_prxy_wci.r > \
          ${ALL_OUT}/log_ps/m6-plot_prxy_wci 2>&1 &
Rscript ${CODE_PLOT_PATH}/m7-plot_bin.r > \
          ${ALL_OUT}/log_ps/m7-plot_bin 2>&1 &

### Plot for multiple figures.
Rscript ${CODE_PLOT_PATH}/s7-plotaa_deltat_true.r > \
          ${ALL_OUT}/log_ps/s7-plotaa_deltat_true 2>&1 &
Rscript ${CODE_PLOT_PATH}/s8-plotmulti_true.r > \
          ${ALL_OUT}/log_ps/s8-plotmulti_true 2>&1 &
