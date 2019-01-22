#!/usr/bin/env bash

if [ -d /opt/conda/envs/morphct ] && [ $BRANCH != 'master' ]; then
	echo "Using Cache";
	source activate morphct
else
	echo "Rebuilding Conda Env";
    rm -rf /opt/conda/envs/morphct
	conda env create -f environment.yml;
        source activate morphct
fi
