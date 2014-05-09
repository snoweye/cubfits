#!/bin/sh

echo "wphi_wophi_run_1.sh"

### Check if configuration file exists.
if [ -e 00-set_env.r ]; then
  echo "00-set_env.r found."
else
  echo "00-set_env.r does not exist."
  exit
fi

### Root of all outputs.
ALL_OUT=`Rscript -e 'source("00-set_env.r");cat(prefix$all.out)'`
if [ "X${ALL_OUT}" = "X" ]; then
  echo "ALL_OUT does not exit."
  exit
else
  echo "export ALL_OUT=${ALL_OUT}"
fi

### Codes.
CODE_PATH=`Rscript -e 'source("00-set_env.r");cat(prefix$code)'`
CODE_PLOT_PATH=`Rscript -e 'source("00-set_env.r");cat(prefix$code.plot)'`

### Read in fake data or real data.
Rscript ${CODE_PATH}/02-get_data.r > \
          ${ALL_OUT}/log/02-get_data 2>&1

### Get initial values.
NP=7
mpiexec -np ${NP} Rscript ${CODE_PATH}/03-get_init-tp.r > \
                            ${ALL_OUT}/log/03-get_init-tp 2>&1

### Plotting.
Rscript ${CODE_PLOT_PATH}/03-plotdiag_bin_est.r > \
          ${ALL_OUT}/log/03-plotdiag_bin_est 2>&1 &
Rscript ${CODE_PLOT_PATH}/03-plotdiag_init.r > \
          ${ALL_OUT}/log/03-plotdiag_init 2>&1 &

### Run MCMC for wophi.
NP=7
nohup mpiexec -np ${NP} Rscript ${CODE_PATH}/04-ad_wophi_pm-tp.r > \
                                  ${ALL_OUT}/log/04-ad_wophi_pm-tp 2>&1 &
RUN_1=$!
nohup mpiexec -np ${NP} Rscript ${CODE_PATH}/04-ad_wophi_scuo-tp.r > \
                                  ${ALL_OUT}/log/04-ad_wophi_scuo-tp 2>&1 &
RUN_2=$!

### Wait to finish all jobs.
wait $RUN_1
wait $RUN_2

### Run MCMC for wphi.
nohup mpiexec -np ${NP} Rscript ${CODE_PATH}/04-ad_wphi_pm-tp.r > \
                                  ${ALL_OUT}/log/04-ad_wphi_pm-tp 2>&1 &
RUN_3=$!
nohup mpiexec -np ${NP} Rscript ${CODE_PATH}/04-ad_wphi_scuo-tp.r > \
                                  ${ALL_OUT}/log/04-ad_wphi_scuo-tp 2>&1 &
RUN_4=$!

### Wait to finish all jobs.
wait $RUN_3
wait $RUN_4

### Subset MCMC wphi results.
NP=5
mpiexec -np ${NP} Rscript ${CODE_PATH}/05-subset-tp.r > \
                            ${ALL_OUT}/log/05-subset-tp 2>&1

### Run MCMC for wophi.
NP=7
nohup mpiexec -np ${NP} Rscript ${CODE_PATH}/04-ad_wphi_wophi_pm-tp.r > \
                                  ${ALL_OUT}/log/04-ad_wphi_wophi_pm-tp 2>&1 &
RUN_5=$!
nohup mpiexec -np ${NP} Rscript ${CODE_PATH}/04-ad_wphi_wophi_scuo-tp.r > \
                                  ${ALL_OUT}/log/04-ad_wphi_wophi_scuo-tp 2>&1 &
RUN_6=$!

### Wait to finish all jobs.
wait $RUN_5
wait $RUN_6

