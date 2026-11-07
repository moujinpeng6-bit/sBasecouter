nextflow.enable.dsl = 2

params.manifests = null
params.outdir = 'runs/recount'
params.cores = 8
params.memory_gb = 64
params.parallel = 1

def quoted(value) {
    "'" + value.toString().replace("'", "'\"'\"'") + "'"
}

process RECOUNT {
    tag { manifest.baseName }
    cpus (params.cores as int)
    memory "${params.memory_gb} GB"
    maxForks (params.parallel as int)
    errorStrategy 'ignore'
    // Raw assets are declared inside the manifest, outside the Nextflow cache key.
    cache false

    input:
    tuple path(manifest), val(destination)

    script:
    """
    python -m sbasecouter recount --manifest ${quoted(manifest)} --out ${quoted(destination)} \\
        --cores ${task.cpus} --memory-gb ${task.memory.toGiga()}
    """
}

workflow {
    if (!params.manifests) {
        error 'Provide --manifests pointing to a directory of sample JSON files'
    }
    def directory = file(params.manifests, checkIfExists: true)
    def output = file(params.outdir).toAbsolutePath()
    jobs = channel.fromPath("${directory}/*.json", checkIfExists: true)
        .map { manifest -> tuple(manifest, output.resolve(manifest.baseName).toString()) }
    RECOUNT(jobs)
}
