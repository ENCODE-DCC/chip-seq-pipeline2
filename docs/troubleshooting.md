
## Trimmomatic errors

### `Error: Unable to detect quality encoding`

Take a look at your FASTQs first. If there are multiple technical replicates for one biological replicate, then pipeline merge all technical replicates' FASTQs first and then crop the merged one with Trimmomatic. If you mix up FASTQs with different base encoding then you will see this error.

Add `phred33` or `phred64` to your input JSON. It's `auto` by default.

```javascript
{
	"chip.trimmomatic_phred_score_format": "phred33"
}
```

## Conda environment

It takes long (>30 minutes) to resolve pipeline's Conda environment since the pipeline uses lots of dependencies including `bowtie2`, `samtools`, `bedtools` and math libraries like `numpy`.

If it takes too long (>hours) then try with a different method other than Conda. Install docker and try with  `--docker` for `caper run` or `caper submit`. You don't need to define docker/singularity images and their versions (image's tag) since they are already defined in pipeline' WDL and caper take it automatically.
