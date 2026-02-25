# MGA rules start here

localrules:
    set_cost_optimal,
    create_mga_parameters,
    convert_MGAresults,
    mga_finished,
    run_mga, # only temporarily

mgapath= modelpath / "mga/{slack}_{objective}_{tech_group}_{zone_group}/"

rule set_cost_optimal:
    input:
        resultspath=rules.convert_results.output.resultsdb,
    output:
        cost_opt=modelpath / "mga/cost_optimal.tsv",
    script:
        "../scripts/set_cost_optimal.py"

rule create_mga_parameters:
    input: 
        cost_optimal=rules.set_cost_optimal.output.cost_opt,
    output:
        mga_params=mgapath / "mga_parameters.dd",
        mgaMode=mgapath / "mgaMode.gms"
    params:
        techs = lambda wc: tech_groups[wc.tech_group],
        zones = lambda wc: zone_groups[wc.zone_group],
        all_zones = aggregated_regions,
    script:
        "../scripts/create_mga_params.py"

rule run_mga:
    input:
        shared_code_path + "highres_data_input.gms",
        shared_code_path + "highres_hydro.gms",
        shared_code_path + "highres_results.gms",
        shared_code_path + "highres_storage_setup.gms",
        shared_code_path + "highres_storage_uc_setup.gms",
        shared_code_path + "highres_uc_setup.gms",
        shared_code_path + "highres_mga.gms",
        modelpath / "inputs.finished",
        modelpath / "vre_{year}_.gdx",
        rules.create_mga_parameters.output.mga_params,
        rules.create_mga_parameters.output.mgaMode,
    params:
        gamspath=gamspath,
        modelpath=str(modelpath),
        mgapath=str(mgapath),
        co2_target_type=config["parameters"]["co2_target_type"],
        co2_target_extent=config["parameters"]["co2_target_extent"],
        sharedcodepath=abs_gams_code_path,
        outname=str(mgapath) + "/results",
        store_initial_level=config["parameters"]["store_initial_level"],
        store_final_level=config["parameters"]["store_final_level"],
        hydro_res_initial_fill=config["parameters"]["hydro_res_initial_fill"],
        hydro_res_min=config["parameters"]["hydro_res_min"],
        pgen=config["parameters"]["pgen"],
        emis_price=config["parameters"]["emis_price"],
        trans_inv=config["parameters"]["trans_inv"],
        trans_cap_lim=config["parameters"]["trans_cap_lim"],
        unit_commitment=config["parameters"]["unit_commitment"],
        frequency_response=config["parameters"]["frequency_response"],
        store_uc=config["parameters"]["store_uc"],
    resources:
        mem_mb=50000,
    log:
        str(mgapath) + "/highres_mga.lst",
        str(mgapath) + "/highres_mga.log",
    output:
        modelresults=protected(mgapath / "results.gdx"),
        # modelresultsdd=protected(modelpath / "results.db"),
    script:
        "../scripts/run_gams_mga.py"

rule convert_MGAresults:
    input:
        rules.run_mga.output.modelresults,
    output:
        resultsdb=ensure(protected(mgapath / "results.db"), non_empty=True),
    shell:
        gamspath + "gamstool sqlitewrite gdxIn={input} o={output}"

rule mga_finished:
    input:
        rules.set_cost_optimal.output.cost_opt,
        rules.create_mga_parameters.output.mga_params,
        rules.run_mga.output.modelresults,
        rules.convert_MGAresults.output.resultsdb,
    output:
        touch(mgapath / "mga.finished"),
