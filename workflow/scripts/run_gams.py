import subprocess

# TODO - test on linux

args = [
    snakemake.params.gamspath + "gams",
    snakemake.params.sharedcodepath / "highres.gms",
    "gdxCompress=1",
    "--codefolderpath=" + str(snakemake.params.sharedcodepath),
    # May bring this back at some stage but commented for now
    #"--co2intensity=" + str(snakemake.params.co2intensity),
    "--co2_target_type=" + str(snakemake.params.co2_target_type),
    "--co2_target_extent=" + str(snakemake.params.co2_target_extent),
    "--weather_yr=" + str(snakemake.wildcards.year),
    "--outname=" + snakemake.params.outname,
    "--store_initial_level=" + str(snakemake.params.store_initial_level),
    "--store_final_level=" + str(snakemake.params.store_final_level),
    "--hydro_res_initial_fill=" + str(snakemake.params.hydro_res_initial_fill),
    "--hydro_res_min=" + str(snakemake.params.hydro_res_min),
    "--pgen=" + str(snakemake.params.pgen),
    "--emis_price=" + str(snakemake.params.emis_price),
    "--trans_inv=" + str(snakemake.params.trans_inv),
    "--trans_cap_lim=" + str(snakemake.params.trans_cap_lim),
]

process = subprocess.Popen(
    args,
    stdout=subprocess.PIPE,
    stderr=subprocess.STDOUT,
    cwd=snakemake.params.modelpath,
)

# Below writes gams output to terminal in realtime and log

with open(snakemake.log[0], "w") as f:
    while True:
        line = process.stdout.readline()

        if not line and process.poll() is not None:
            break
        print(line.decode(), end="")
        f.write(line.decode().rstrip("\n"))
