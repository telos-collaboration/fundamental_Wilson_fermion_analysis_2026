from glob import glob


rule package_smeared:
    input:
        files=glob("raw_data/spectrum/*"),
        script="spectrum/src/parse_mesons.jl",
        julia_instantiated="intermediary_data/julia_ready",
    output:
        h5file="data_assets/spectrum/corr_sp4_FUN.h5",
    conda:
        "../envs/scattering.yml"
    # Start packaging early,
    # since it is time consuming and many other processes depend on it
    priority: 10
    shell:
        "julia {input.script} --h5file {output.h5file} {input.files}"


rule package_gflow:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
    input:
        files=glob("raw_data/gradient_flow/*/topology/out/out_flow"),
        script="spectrum/src/package_flows.py",
    output:
        h5="data_assets/spectrum/nf2_gflow.h5",
    conda:
        "../envs/flow_analysis.yml"
    shell:
        "python -m {params.module} {input.files} --h5_filename {output.h5}"
