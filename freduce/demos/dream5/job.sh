#!/bin/sh

#PBS -l nodes=1,walltime=00:30:00,mem=4G

source ../../test/set_environment_yeti.sh
bash fitModel.sh | tee logfile
