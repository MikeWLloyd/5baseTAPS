# JAX-GT/5baseTAPS: Contributing Guidelines

Hi there! Many thanks for taking an interest in improving JAX-GT/5baseTAPS.

We manage tasks and bug reports using [GitHub Issues](https://github.com/TheJacksonLaboratory/5baseTAPS/issues). Please use the pre-filled issue templates to save time.

## Contribution workflow

1. Check that there isn't already an issue about your idea in the [issue tracker](https://github.com/TheJacksonLaboratory/5baseTAPS/issues).
2. [Fork](https://help.github.com/en/github/getting-started-with-github/fork-a-repo) the [5baseTAPS repository](https://github.com/TheJacksonLaboratory/5baseTAPS) to your GitHub account.
3. Make the necessary changes within your forked repository following the [Pipeline conventions](#pipeline-contribution-conventions) below.
4. Submit a Pull Request against the `dev` branch and wait for review.

If you're not used to this workflow, see GitHub's [collaborating with pull requests](https://help.github.com/en/github/collaborating-with-issues-and-pull-requests) guide.

## Tests

Validate the pipeline config syntax locally before submitting:

```bash
nextflow config .
```

Full testing (stub or real run) requires HPC access with reference genome and real input data. GitHub Actions CI runs config validation automatically on each PR.

## Patch releases

If a bug is discovered after a release:

- On your fork, create a `patch` branch based on `master`.
- Fix the bug and bump the patch version (X.Y.Z+1).
- Open a PR from `patch` to `master`.

## Pipeline contribution conventions

### Adding a new step

1. Define the input channel from the expected previous process.
2. Write the process block.
3. Define the output channel if needed.
4. Add any new parameters to `nextflow.config` with a default value.
5. Add any new parameters to `nextflow_schema.json` with help text.
6. Add sanity checks for all relevant parameters.
7. Validate locally before submitting.
8. If applicable, add a test in `.github/workflows/ci.yml`.
9. Update `assets/multiqc_config.yml` for any new MultiQC modules.
10. Document new output files in `docs/output.md`.

### Default values

Parameters should be initialised in the `params` scope in `nextflow.config`.

### Resource requirements

Default process resources (CPUs / memory / time) are defined in `conf/base.config` using `withLabel:` selectors. The process resources are passed to tools dynamically via `${task.cpus}` and `${task.memory}` in the `script:` block.

### Naming schemes

- Initial process channels: `ch_output_from_<process>`
- Intermediate/terminal channels: `ch_<previousprocess>_for_<nextprocess>`
