#!/usr/bin/env bash

# TODO: export environment variable to then have gams read the path from there.
# TODO: also call gams from here and then set the log file path here too


cd ${snakemake_params[modelpath]}
pwd -P

${snakemake_params[gamspath]}gams \
${snakemake_params[sharedcodepath]}/highres.gms \
gdxCompress=1 \
--weather_yr "${snakemake_wildcards[year]}" \
--dem_yr "${snakemake_wildcards[year]}" \
--codefolderpath "${snakemake_params[sharedcodepath]}" \
--co2intensity "${snakemake_params[co2intensity]}"
