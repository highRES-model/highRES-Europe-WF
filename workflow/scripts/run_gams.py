import subprocess

# TODO - test on linux

args = [
    snakemake.params.gamspath + "gams",
    snakemake.params.sharedcodepath + "/highres.gms",
    "gdxCompress=1",
    "--co2intensity=" + str(snakemake.params.co2intensity),
    "--weather_yr=" + str(snakemake.wildcards.year),
    "--vre_restrict="+str(snakemake.wildcards.quantile),
    "--outname="+snakemake.params.outname
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
