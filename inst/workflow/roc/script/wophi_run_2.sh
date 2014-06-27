#!/bin/sh

echo "wophi_run_2.sh"

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
CODE_PLOT_PATH=`Rscript -e 'source("00-set_env.r");cat(prefix$code.plot)'`

### Check parallel.
NP=5
PARALLEL=`Rscript -e 'source("00-set_env.r");cat(run.info$parallel)'`
if [ "$PARALLEL" = "task.pull" ] || [ "$PARALLEL" = "pbdLapply" ]; then
  MPI_EXEC="mpiexec -np ${NP}"
else
  MPI_EXEC=
fi

### Plot.
Rscript ${CODE_PLOT_PATH}/03-plotdiag_bin_est.r > \
          ${ALL_OUT}/log/03-plotdiag_bin_est 2>&1 &
Rscript ${CODE_PLOT_PATH}/03-plotdiag_init.r > \
          ${ALL_OUT}/log/03-plotdiag_init 2>&1 &


### Subset MCMC results.
${MPI_EXEC} Rscript ${CODE_PATH}/05-subset-tp.r > \
                      ${ALL_OUT}/log/05-subset-tp 2>&1

### Dump tsv files.
Rscript ${CODE_PLOT_PATH}/05-subset_tsv.r > \
          ${ALL_OUT}/log/05-subset_tsv 2>&1 &

### Plot fitted results.
Rscript ${CODE_PLOT_PATH}/06-plotsingle_model.r > \
          ${ALL_OUT}/log/06-plotsingle_model 2>&1 &
Rscript ${CODE_PLOT_PATH}/06-plotsingle_prxy.r > \
          ${ALL_OUT}/log/06-plotsingle_prxy 2>&1 &
Rscript ${CODE_PLOT_PATH}/06-plotsingle_prxy_wci.r > \
          ${ALL_OUT}/log/06-plotsingle_prxy_wci 2>&1 &

### Plot diagnoses.
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

### Plot traces.
Rscript ${CODE_PLOT_PATH}/07-plottrace_param_meanEPhi.r > \
          ${ALL_OUT}/log/07-plottrace_param_meanEPhi 2>&1 &
Rscript ${CODE_PLOT_PATH}/07-plottrace_prior.r > \
          ${ALL_OUT}/log/07-plottrace_prior 2>&1 &
Rscript ${CODE_PLOT_PATH}/07-plottrace_quantile_Phi.r > \
          ${ALL_OUT}/log/07-plottrace_quantile_Phi 2>&1 &
