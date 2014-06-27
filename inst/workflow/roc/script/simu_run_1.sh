#!/bin/sh

echo "simu_run_1.sh"

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

### Check parallel.
NP=5
PARALLEL=`Rscript -e 'source("00-set_env.r");cat(run.info$parallel)'`
if [ "$PARALLEL" = "task.pull" ] || [ "$PARALLEL" = "pbdLapply" ]; then
  MPI_EXEC="mpiexec -np ${NP}"
else
  MPI_EXEC=
fi

### Generate fake data for simulation only.
Rscript ${CODE_PATH}/s1-generate_data.r > \
          ${ALL_OUT}/log/s1-generate_data 2>&1

### Read in fake data or real data.
Rscript ${CODE_PATH}/02-get_data.r > \
          ${ALL_OUT}/log/02-get_data 2>&1

### Get initial values.
${MPI_EXEC} Rscript ${CODE_PATH}/03-get_init-tp.r > \
                      ${ALL_OUT}/log/03-get_init-tp 2>&1

### Plot.
Rscript ${CODE_PLOT_PATH}/03-plotdiag_bin_est.r > \
          ${ALL_OUT}/log/03-plotdiag_bin_est 2>&1 &
Rscript ${CODE_PLOT_PATH}/03-plotdiag_init.r > \
          ${ALL_OUT}/log/03-plotdiag_init 2>&1 &

### Run MCMC.
nohup ${MPI_EXEC} Rscript ${CODE_PATH}/04-wophi_pm-tp.r > \
                            ${ALL_OUT}/log/04-wophi_pm-tp 2>&1 &
RUN_1=$!
nohup ${MPI_EXEC} Rscript ${CODE_PATH}/04-wophi_scuo-tp.r > \
                            ${ALL_OUT}/log/04-wophi_scuo-tp 2>&1 &
RUN_2=$!
wait $RUN_1
wait $RUN_2
nohup ${MPI_EXEC} Rscript ${CODE_PATH}/04-wphi_pm-tp.r > \
                            ${ALL_OUT}/log/04-wphi_pm-tp 2>&1 &
RUN_3=$!
nohup ${MPI_EXEC} Rscript ${CODE_PATH}/04-wphi_scuo-tp.r > \
                            ${ALL_OUT}/log/04-wphi_scuo-tp 2>&1 &
RUN_4=$!
wait $RUN_3
wait $RUN_4
# nohup ${MPI_EXEC} Rscript ${CODE_PATH}/s4-wophi_true-tp.r > \
#                             ${ALL_OUT}/log/s4-wophi_true-tp 2>&1 &
# RUN_5=$!
# nohup ${MPI_EXEC} Rscript ${CODE_PATH}/s4-wphi_true-tp.r > \
#                             ${ALL_OUT}/log/s4-wphi_true-tp 2>&1 &
# RUN_6=$!

### Wait to finish all jobs.
# wait $RUN_5
# wait $RUN_6
