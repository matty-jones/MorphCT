#!/usr/bin/env bash

if [ -d /opt/conda/envs/morphct ]; then
	echo "Using Cache";
	source activate morphct
else
	echo "Rebuilding Conda Env";
	conda env create -f environment.yml;
        source activate morphct
fi
