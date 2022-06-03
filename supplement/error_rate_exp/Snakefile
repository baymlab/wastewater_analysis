import json
import pprint
import numpy as np
from math import log10, floor


configfile: "snek_config.yaml"


pangolin = [x + "_EPI_ISL_" + y for x, y in zip(config["vocs"], config["lineages"])]


def round_sig(x, sig=2):
    if x == 0:
        return x
    return round(x, sig - int(floor(log10(abs(x)))) - 1)


abus = [
    round(float(x), 2)
    for x in (
        config["abundances"]
        if not config["abundances_scale_log"]
        else np.geomspace(
            config["abundances_log"][0],
            config["abundances_log"][1],
            config["abundances_log"][2],
        )
    )
]
errors = [
    round_sig(float(x), 3)
    for x in (
        config["errors"]
        if not config["errors_scale_log"]
        else np.geomspace(
            config["errors_log"][1],
            config["errors_log"][0],
            config["errors_log"][2],
        )
    )
]


wildcard_constraints:
    ref="[^_]+",
    dataset="(?!des)[^/]+",
    format="(?!contamination)[^_]+",


# format="[^_]+",


rule create_benchmark_error_compare:
    input:
        fasta="genome_data/sequences.fasta",
        metadata="genome_data/metadata.tsv",
        voc=expand("genome_data/{voc}.fasta", voc=pangolin),
    output:
        expand(
            "benchmarks/{{dataset}}_{{format}}/wwsim_{voc}_ab{ab}_er{er}_1.fastq",
            voc=pangolin,
            ab=abus,
            er=errors,
        ),
        expand(
            "benchmarks/{{dataset}}_{{format}}/wwsim_{voc}_ab{ab}_er{er}_2.fastq",
            voc=pangolin,
            ab=abus,
            er=errors,
        ),
        snek=touch("benchmarks/{dataset}_{format}/snek"),
    threads: 2
    params:
        vocs=lambda wildcards, input: ",".join(input.voc),
        errs=lambda wildcards: ",".join([str(er) for er in errors]),
        percs=lambda wildcards: ",".join([str(ab) for ab in abus]),
        spike=lambda wildcards: "--spike_only" if config["spike_only"] else "",
        json=lambda wildcards: json.dumps(config),
    run:
        if "chimeric" in wildcards.format:
            shell(
                "python benchmarking/create_chimeric_benchmarks.py "
                "--voc_perc {params.percs} "
                "--chim_perc {params.errs} "
                "-m {input.metadata} -fr {input.fasta} "
                "-fv {params.vocs} "
                "-o benchmarks/{wildcards.dataset}_{wildcards.format} "
                "--total_cov {config[tot_cov]} "
                "-d '-' "
                "{params.spike} "
                "-d '-' "
            )
        else:
            sid = "--sub_error " if "s" in wildcards.format else ""
            sid += "--ins_error " if "i" in wildcards.format else ""
            sid += "--del_error " if "d" in wildcards.format else ""
            shell(
                "python benchmarking/create_error_benchmarks.py "
                "--voc_perc {params.percs} "
                "--err_perc {params.errs} "
                "-m {input.metadata} -fr {input.fasta} "
                "-fv {params.vocs} "
                "-o benchmarks/{wildcards.dataset}_{wildcards.format} "
                "--total_cov {config[tot_cov]} "
                "-d '-' "
                "{params.spike} "
                "{sid} "
                "-d '-' "
            )
        shell("echo {params.json} > {output.snek}")


# rule create_benchmark_designer:
#     input:
#         fasta="reference_sets/designer/{ref}/sequences.fasta",
#         metadata="reference_sets/designer/{ref}/metadata.tsv",
#         voc=expand(
#             "reference_sets/designer/{{ref}}/{voc}.fa", voc=config["designer_vocs"]
#         ),
#     output:
#         expand(
#             "benchmarks/des_{{ref}}_{{dataset}}_{{format}}/wwsim_{voc}_er{er}_1.fastq",
#             voc=config["designer_vocs"],
#             er=errors,
#         ),
#         expand(
#             "benchmarks/des_{{ref}}_{{dataset}}_{{format}}/wwsim_{voc}_er{er}_2.fastq",
#             voc=config["designer_vocs"],
#             er=errors,
#         ),
#         snek=touch("benchmarks/des_{ref}_{dataset}_{format}/snek"),
#     threads: 2
#     params:
#         vocs=lambda wildcards, input: ",".join(input.voc),
#         errs=lambda wildcards: ",".join([str(er) for er in errors]),
#         spike=lambda wildcards: "--spike_only" if config["spike_only"] else "",
#         json=lambda wildcards: json.dumps(config),
#     run:
#         sid = "--sub_error " if "s" in wildcards.format else ""
#         sid += "--ins_error " if "i" in wildcards.format else ""
#         sid += "--del_error " if "d" in wildcards.format else ""
#         shell(
#             "python benchmarking/create_designer_benchmarks.py "
#             "--voc_perc {config[designer_abu]} "
#             "--err_perc {params.errs} "
#             "-m {input.metadata} -fr {input.fasta} "
#             "-fv {params.vocs} "
#             "-o benchmarks/des_{wildcards.ref}_{wildcards.dataset}_{wildcards.format} "
#             "--total_cov {config[tot_cov]} "
#             "-d '-' "
#             "{params.spike} "
#             "{sid} "
#         )
#         shell("echo {params.json} > {output.snek}")


rule create_benchmark_contamination:
    input:
        fasta="genome_data/sequences.fasta",
        metadata="genome_data/metadata.tsv",
        voc=expand("genome_data/{voc}.fasta", voc=pangolin),
        cont=expand(
            "genome_data/contamination/{contaminant}.fasta",
            contaminant=config["contaminants"],
        ),
    output:
        expand(
            "benchmarks/{{dataset}}_contamination/wwsim_{voc}_ab{ab}_1.fastq",
            voc=pangolin,
            ab=abus,
        ),
        expand(
            "benchmarks/{{dataset}}_contamination/wwsim_{voc}_ab{ab}_2.fastq",
            voc=pangolin,
            ab=abus,
        ),
        snek=touch("benchmarks/{dataset}_contamination/snek"),
    params:
        vocs=lambda wildcards, input: ",".join(input.voc),
        contaminants=lambda wildcards, input: ",".join(input.cont),
        percs=lambda wildcards: ",".join([str(ab) for ab in abus]),
        spike=lambda wildcards: "--spike_only" if config["spike_only"] else "",
        json=lambda wildcards: json.dumps(config),
    run:
        shell(
            "python benchmarking/create_contamination_benchmarks.py "
            "-m {input.metadata} -fr {input.fasta} "
            "-fv {params.vocs} "
            "-fc {params.contaminants} "
            "-o benchmarks/{wildcards.dataset}_contamination "
            "--sars2_perc {params.percs} "
            "--total_sars2_cov {config[tot_cov]} "
            "--conts_amount {config[conts_amount]} "
            "{params.spike} "
            "-d '-' "
            "--no_errors "
        )
        shell("echo {params.json} > {output.snek}")


rule run_kallisto_batch_jobs:
    input:
        idx=expand("{ref}/sequences.kallisto_idx", ref=config["ref"]),
        wwsim1=expand(
            "benchmarks/{{dataset}}_{{format}}/wwsim_{voc}_ab{ab}_er{{er}}_1.fastq",
            voc=pangolin,
            ab=abus,
        ),
        wwsim2=expand(
            "benchmarks/{{dataset}}_{{format}}/wwsim_{voc}_ab{ab}_er{{er}}_2.fastq",
            voc=pangolin,
            ab=abus,
        ),
    output:
        preds=expand(
            "benchmarks/{{dataset}}_{{format}}/out/{voc}_ab{ab}_er{{er}}/predictions_m{min_ab}.tsv",
            voc=pangolin,
            ab=abus,
            min_ab=config["min_ab"],
        ),
    params:
        vocs=lambda wildcards, input: ",".join(config["vocs"]),
    threads: 2
    run:
        outdir = "/".join(str(output.preds[0]).split("/")[0:-2])
        shell("mkdir -p {outdir}")
        for voc in pangolin:
            for ab in abus:
                shell(
                    # "srun --ntasks=1 --cpus-per-task={threads} "
                    "kallisto quant -t {threads} -b {config[bootstraps]} "
                    "-i {input.idx} -o {outdir}/{voc}_ab{ab}_er{wildcards.er} "
                    "benchmarks/{wildcards.dataset}_{wildcards.format}/wwsim_{voc}_ab{ab}_er{wildcards.er}_1.fastq "
                    "benchmarks/{wildcards.dataset}_{wildcards.format}/wwsim_{voc}_ab{ab}_er{wildcards.er}_2.fastq "
                )
                shell(
                    # "srun --ntasks=1 --cpus-per-task=1 "
                    "python pipeline/output_abundances.py "
                    "-m {config[min_ab]} "
                    "-o {outdir}/{voc}_ab{ab}_er{wildcards.er}/predictions_m{config[min_ab]}.tsv "
                    "--metadata {config[ref]}/metadata.tsv "
                    "--voc {params.vocs} "
                    "{outdir}/{voc}_ab{ab}_er{wildcards.er}/abundance.tsv "
                )


# rule run_kallisto_batch_jobs_designer:
#     input:
#         idx="reference_sets/designer/{ref}/sequences.kallisto_idx",
#         ref=config["ref"],
#         wwsim1=expand(
#             "benchmarks/des_{{ref}}_{{dataset}}_{{format}}/wwsim_{{voc}}_er{er}_1.fastq",
#             er=errors,
#         ),
#         wwsim2=expand(
#             "benchmarks/des_{{ref}}_{{dataset}}_{{format}}/wwsim_{{voc}}_er{er}_2.fastq",
#             er=errors,
#         ),
#     output:
#         preds=expand(
#             "benchmarks/des_{{ref}}_{{dataset}}_{{format}}/out/{{voc}}_er{er}/predictions_m{min_ab}.tsv",
#             er=errors,
#             min_ab=config["min_ab"],
#         ),
#     params:
#         vocs=lambda wildcards, input: ",".join(config["designer_vocs"]),
#     threads: 12
#     run:
#         outdir = "/".join(str(output.preds[0]).split("/")[0:-2])
#         shell("mkdir -p {outdir}")
#         for er in errors:
#             shell(
#                 # "srun --ntasks=1 --cpus-per-task={threads} "
#                 "kallisto quant -t {threads} -b {config[bootstraps]} "
#                 "-i {input.idx} -o {outdir}/{wildcards.voc}_er{er} "
#                 "benchmarks/des_{wildcards.ref}_{wildcards.dataset}_{wildcards.format}/wwsim_{wildcards.voc}_er{er}_1.fastq "
#                 "benchmarks/des_{wildcards.ref}_{wildcards.dataset}_{wildcards.format}/wwsim_{wildcards.voc}_er{er}_2.fastq "
#             )
#             shell(
#                 # "srun --ntasks=1 --cpus-per-task=1 "
#                 "python pipeline/output_abundances.py "
#                 "-m {config[min_ab]} "
#                 "-o {outdir}/{wildcards.voc}_er{er}/predictions_m{config[min_ab]}.tsv "
#                 "--metadata reference_sets/designer/{wildcards.ref}/metadata.tsv "
#                 "--voc {params.vocs} "
#                 "{outdir}/{wildcards.voc}_er{er}/abundance.tsv "
#             )


rule run_kallisto_batch_jobs_contamination:
    input:
        idx=expand("{ref}/sequences.kallisto_idx", ref=config["ref"]),
        wwsim1=expand(
            "benchmarks/{{dataset}}_contamination/wwsim_{{voc}}_ab{ab}_1.fastq",
            ab=abus,
        ),
        wwsim2=expand(
            "benchmarks/{{dataset}}_contamination/wwsim_{{voc}}_ab{ab}_2.fastq",
            ab=abus,
        ),
    output:
        preds=expand(
            "benchmarks/{{dataset}}_contamination/out/{{voc}}_ab{ab}/predictions_m{min_ab}.tsv",
            ab=abus,
            min_ab=config["min_ab"],
        ),
    params:
        vocs=lambda wildcards, input: ",".join(config["vocs"]),
    threads: 12
    run:
        outdir = "/".join(str(output.preds[0]).split("/")[0:-2])
        shell("mkdir -p {outdir}")
        for ab in abus:
            shell(
                # "srun --ntasks=1 --cpus-per-task={threads} "
                "kallisto quant -t {threads} -b {config[bootstraps]} "
                "-i {input.idx} -o {outdir}/{wildcards.voc}_ab{ab} "
                "benchmarks/{wildcards.dataset}_contamination/wwsim_{wildcards.voc}_ab{ab}_1.fastq "
                "benchmarks/{wildcards.dataset}_contamination/wwsim_{wildcards.voc}_ab{ab}_2.fastq "
            )
            shell(
                # "srun --ntasks=1 --cpus-per-task=1 "
                "python pipeline/output_abundances.py "
                "-m {config[min_ab]} "
                "-o {outdir}/{wildcards.voc}_ab{ab}/predictions_m{config[min_ab]}.tsv "
                "--metadata {config[ref]}/metadata.tsv "
                "--voc {params.vocs} "
                "{outdir}/{wildcards.voc}_ab{ab}/abundance.tsv "
            )


rule create_figs_compare_error:
    input:
        expand(
            "benchmarks/{{dataset}}_{{format}}/out/{voc}_ab{ab}_er{er}/predictions_m{min_ab}.tsv",
            voc=pangolin,
            ab=abus,
            er=errors,
            min_ab=config["min_ab"],
        ),
    output:
        expand(
            "benchmarks/figs/{{dataset}}_{{format}}/{file}_{fonts}.{ext}",
            ext=config["plot_exts"],
            file=[
                "freq_error_plot",
                "freq_error_plot_logscale",
                "freq_scatter_loglog",
                "error_error_plot",
                "error_error_plot_logscale",
            ],
            fonts=config["plot_font_sizes"],
        ),
        snek="benchmarks/figs/{dataset}_{format}/snek",
        dir=directory("benchmarks/figs/{dataset}_{format}"),
    params:
        vocs=lambda wildcards, input: ",".join(config["vocs"]),
        exts=lambda wildcards, input: ",".join(config["plot_exts"]),
        json=lambda wildcards, input: json.dumps(config),
    run:
        for size in config["plot_font_sizes"]:
            shell(
                "python benchmarking/evaluate_error.py "
                "--voc {params.vocs} "
                "--plot_abundance_value {config[plot_abundance_value]} "
                "--plot_error_value {config[plot_error_value]} "
                "-o {output.dir} "
                "--output_format {params.exts} "
                "--font_size {size} "
                "--suffix _{size} "
                "benchmarks/{wildcards.dataset}_{wildcards.format}/out/*/predictions_m{config[min_ab]}.tsv "
                "&& echo {params.json} > {output.snek}"
                # "-m {config[min_ab]} "
            )


# rule create_figs_compare_designer:
#     input:
#         expand(
#             "benchmarks/des_{{ref}}_{{dataset}}_{{format}}/out/{voc}_er{er}/predictions_m{min_ab}.tsv",
#             voc=config["designer_vocs"],
#             ab=abus,
#             er=errors,
#             min_ab=config["min_ab"],
#         ),
#     output:
#         expand(
#             "benchmarks/figs/des_{{ref}}_{{dataset}}_{{format}}/{file}_{fonts}.{ext}",
#             ext=config["plot_exts"],
#             file=[
#                 "freq_error_plot",
#                 "freq_error_plot_logscale",
#                 "freq_scatter_loglog",
#                 "error_error_plot",
#                 "error_error_plot_logscale",
#             ],
#             fonts=config["plot_font_sizes"],
#         ),
#         snek="benchmarks/figs/des_{ref}_{dataset}_{format}/snek",
#         dir=directory("benchmarks/figs/des_{ref}_{dataset}_{format}"),
#     params:
#         vocs=lambda wildcards, input: ",".join(config["designer_vocs"]),
#         exts=lambda wildcards, input: ",".join(config["plot_exts"]),
#         json=lambda wildcards, input: json.dumps(config),
#     run:
#         for size in config["plot_font_sizes"]:
#             shell(
#                 "python benchmarking/evaluate_designer.py "
#                 "--voc {params.vocs} "
#                 "--abundance {config[designer_abu]} "
#                 "--plot_error_value {config[plot_error_value]} "
#                 "-o {output.dir} "
#                 "--output_format {params.exts} "
#                 "--font_size {size} "
#                 "--suffix _{size} "
#                 "benchmarks/des_{wildcards.ref}_{wildcards.dataset}_{wildcards.format}/out/*/predictions_m{config[min_ab]}.tsv "
#                 "&& echo {params.json} > {output.snek}"
#                 # "-m {config[min_ab]} "
#             )


rule create_figs_contamination:
    input:
        expand(
            "benchmarks/{{dataset}}_contamination/out/{voc}_ab{ab}/predictions_m{min_ab}.tsv",
            voc=pangolin,
            ab=abus,
            min_ab=config["min_ab"],
        ),
    output:
        expand(
            "benchmarks/figs/{{dataset}}_contamination/{file}_{fonts}.{ext}",
            ext=config["plot_exts"],
            file=[
                "freq_error_plot",
                "freq_error_plot_logscale",
                "freq_scatter_loglog",
            ],
            fonts=config["plot_font_sizes"],
        ),
        snek="benchmarks/figs/{dataset}_contamination/snek",
        dir=directory("benchmarks/figs/{dataset}_contamination"),
    params:
        vocs=lambda wildcards, input: ",".join(config["vocs"]),
        exts=lambda wildcards, input: ",".join(config["plot_exts"]),
        json=lambda wildcards, input: json.dumps(config),
    run:
        cir = ""
        if config["conts_in_ref"]:
            cir = "--conts_in_meta"
        for size in config["plot_font_sizes"]:
            shell(
                "python benchmarking/evaluate_contamination.py "
                "--voc {params.vocs} "
                "-o {output.dir} "
                "--output_format {params.exts} "
                "--font_size {size} "
                "--suffix _{size} "
                " {cir} "
                "--full_cont_count {config[conts_amount]} "
                "benchmarks/{wildcards.dataset}_contamination/out/*/predictions_m{config[min_ab]}.tsv "
                "&& echo {params.json} > {output.snek}"
            )
