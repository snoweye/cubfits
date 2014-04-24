#!/bin/sh

echo "simu_run_2.sh"

# Check if configuration file exists.
if [ -e 00-set_env.r ]; then
  echo "00-set_env.r found."
else
  echo "00-set_env.r does not exist."
  exit
fi

# Root of all outputs.
ALL_OUT=`Rscript -e 'source("00-set_env.r");cat(prefix$all.out)'`

# Codes.
CODE_PATH=`Rscript -e 'source("00-set_env.r");cat(prefix$code)'`
CODE_PLOT_PATH=`Rscript -e 'source("00-set_env.r");cat(prefix$code.plot)'`

# Plotting for fake data only.
Rscript ${CODE_PLOT_PATH}/s1-plotdiag_simu_phi.r > \
          ${ALL_OUT}/log/s1-plotdiag_simu_phi 2>&1 &
Rscript ${CODE_PLOT_PATH}/s2-plotdiag_bin_true.r > \
          ${ALL_OUT}/log/s2-plotdiag_bin_true 2>&1 &
Rscript ${CODE_PLOT_PATH}/s3-plotdiag_bin_est.r > \
          ${ALL_OUT}/log/s3-plotdiag_bin_est 2>&1 &

# Plotting.
Rscript ${CODE_PLOT_PATH}/03-plotdiag_bin_est.r > \
          ${ALL_OUT}/log/03-plotdiag_bin_est 2>&1 &
Rscript ${CODE_PLOT_PATH}/03-plotdiag_init.r > \
          ${ALL_OUT}/log/03-plotdiag_init 2>&1 &


# Subset MCMC results.
NP=7
mpiexec -np ${NP} Rscript ${CODE_PATH}/05-subset-tp.r > \
                            ${ALL_OUT}/log/05-subset-tp 2>&1
Rscript ${CODE_PATH}/05-subset_tsv.r > \
          ${ALL_OUT}/log/05-subset_tsv 2>&1 &

# Plotting fitted results.
Rscript ${CODE_PLOT_PATH}/s6-plotsingle_model_true.r > \
          ${ALL_OUT}/log/s6-plotsingle_model_true 2>&1 &
Rscript ${CODE_PLOT_PATH}/06-plotsingle_model.r > \
          ${ALL_OUT}/log/06-plotsingle_model 2>&1 &
Rscript ${CODE_PLOT_PATH}/06-plotsingle_prxy.r > \
          ${ALL_OUT}/log/06-plotsingle_prxy 2>&1 &
Rscript ${CODE_PLOT_PATH}/06-plotsingle_prxy_wci.r > \
          ${ALL_OUT}/log/06-plotsingle_prxy_wci 2>&1 &

# Plotting diagnoses.
Rscript ${CODE_PLOT_PATH}/07-plotdiag_scuo_cai.r > \
          ${ALL_OUT}/log/07-plotdiag_scuo_cai 2>&1 &
Rscript ${CODE_PLOT_PATH}/07-plotdiag_accept_hist.r > \
          ${ALL_OUT}/log/07-plotdiag_accept_hist 2>&1 &
Rscript ${CODE_PLOT_PATH}/07-plotdiag_acceptvsEPhi.r > \
          ${ALL_OUT}/log/07-plotdiag_acceptvsEPhi 2>&1 &
Rscript ${CODE_PLOT_PATH}/07-plotdiag_bin.r > \
          ${ALL_OUT}/log/07-plotdiag_bin 2>&1 &
Rscript ${CODE_PLOT_PATH}/07-plotdiag_sigmaW_hist.r > \
          ${ALL_OUT}/log/07-plotdiag_sigmaW_hist 2>&1 &
Rscript ${CODE_PLOT_PATH}/07-plotdiag_EPhi_hist.r > \
          ${ALL_OUT}/log/07-plotdiag_EPhi_hist 2>&1 &
Rscript ${CODE_PLOT_PATH}/07-plotdiag_medPhi_EPhi.r > \
          ${ALL_OUT}/log/07-plotdiag_medPhi_EPhi 2>&1 &

# Plotting traces.
Rscript ${CODE_PLOT_PATH}/07-plottrace_param_meanEPhi.r > \
          ${ALL_OUT}/log/07-plottrace_param_meanEPhi 2>&1 &
Rscript ${CODE_PLOT_PATH}/07-plottrace_prior.r > \
          ${ALL_OUT}/log/07-plottrace_prior 2>&1 &
Rscript ${CODE_PLOT_PATH}/07-plottrace_quantile_Phi.r > \
          ${ALL_OUT}/log/07-plottrace_quantile_Phi 2>&1 &

# Plotting for fake data only.
Rscript ${CODE_PLOT_PATH}/s6-plotsingle_b_corr_true.r > \
          ${ALL_OUT}/log/s6-plotsingle_b_corr_true 2>&1 &
Rscript ${CODE_PLOT_PATH}/s6-plotsingle_prxy_true.r > \
          ${ALL_OUT}/log/s6-plotsingle_prxy_true 2>&1 &
Rscript ${CODE_PLOT_PATH}/s6-plotsingle_prxy_true_wci.r > \
          ${ALL_OUT}/log/s6-plotsingle_prxy_true_wci 2>&1 &
Rscript ${CODE_PLOT_PATH}/s6-plotsingle_scu_mscu.r > \
          ${ALL_OUT}/log/s6-plotsingle_scu_mscu 2>&1 &

# Plotting for matched cases only.
Rscript ${CODE_PLOT_PATH}/m6-plot_b_corr.r > \
          ${ALL_OUT}/log/m6-plot_b_corr 2>&1 &
Rscript ${CODE_PLOT_PATH}/m6-plot_b_corr_negsel.r > \
          ${ALL_OUT}/log/m6-plot_b_corr_negsel 2>&1 &
Rscript ${CODE_PLOT_PATH}/m6-plot_prxy.r > \
          ${ALL_OUT}/log/m6-plot_prxy 2>&1 &
Rscript ${CODE_PLOT_PATH}/m6-plot_prxy_wci.r > \
          ${ALL_OUT}/log/m6-plot_prxy_wci 2>&1 &

# Plotting for multiple figures.
Rscript ${CODE_PLOT_PATH}/s8-plotmulti_true.r > \
          ${ALL_OUT}/log/s8-plotmulti_true 2>&1 &
