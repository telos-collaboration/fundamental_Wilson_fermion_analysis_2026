from functools import partial


metadata = pd.read_csv("metadata/spectrum/ensemble_metadata.csv")


def all_samples(wildcards, observables):
    return [
        f"intermediary_data/{dir_template}/{observable}_samples.json".format(**row)
        for observable in observables
        for row in metadata.to_dict(orient="records")
        if row["use_in_main_plots"]
    ]


rule plot_w0_vs_mpcac:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
    input:
        data=partial(all_samples, observables=["w0", "mpcac"]),
        script="spectrum/src/plots/w0_vs_mpcac.py",
    output:
        plot="assets/spectrum/plots/w0_vs_mpcac.{plot_filetype}",
    conda:
        "../envs/flow_analysis.yml"
    shell:
        "python -m {params.module} {input.data} --plot_styles {plot_styles} --plot_file {output.plot}"


def plaquette_data(wildcards):
    return [
        f"intermediary_data/{dir_template}_{{start_type}}start/plaquette_mean.csv".format(
            **row
        )
        for row in metadata.to_dict(orient="records")
        if row["use_in_plaquette_phase_diagram"]
    ]


rule plaquette_phase_diagram:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
    input:
        data=plaquette_data,
        script="spectrum/src/plots/plaquette_phasediagram.py",
    output:
        plot=f"assets/spectrum/plots/plaquette_phasediagram.{plot_filetype}",
        definitions="assets/spectrum/definitions/plaquette_phasediagram_betas.tex",
    conda:
        "../envs/flow_analysis.yml"
    shell:
        "python -m {params.module} {input.data} --plot_styles {plot_styles} --plot_file {output.plot} --definitions_file {output.definitions}"
