# CLI

PlotNado includes a template-driven CLI for users who prefer a file-based workflow.

The standard flow is:

1. `plotnado init` to infer a YAML template from track files
2. `plotnado validate` to catch missing files or bad group references
3. `plotnado plot` to render one or more regions from that template

## `plotnado init`

Generate a starting template from track files:

```bash
plotnado init sample1.bw sample2.bw peaks.narrowpeak --auto --output template.yaml
```

Useful options:

- `--output` / `-o`: where to write the YAML file
- `--genome` / `-g`: set a default genome such as `hg38` or `mm10`
- `--group-by`: group tracks by a predefined strategy or regex
- `--auto`: skip interactive prompts and use inferred defaults
- `--no-genes`: do not add a genes guide track by default

Examples:

```bash
plotnado init *.bw --auto
plotnado init sample1_H3K27ac.bw sample1_H3K4me3.bw sample2_H3K27ac.bw --group-by sample
plotnado init control_r1.bw control_r2.bw treat_r1.bw treat_r2.bw --group-by '([^_]+)_r[0-9]'
```

The generated template is plain YAML and intended to be edited. A typical file looks like:

```yaml
genome: hg38
guides:
  genes: true
tracks:
  - path: sample1.bw
    type: bigwig
    title: sample1
    group: sample
  - path: peaks.narrowpeak
    type: narrowpeak
    title: peaks
```

## `plotnado validate`

Validate a template before rendering:

```bash
plotnado validate template.yaml
```

Validation checks:

- template can be loaded as YAML
- local track files exist
- group references resolve against track `name` or `title`
- the template compiles into a render plan cleanly

## `plotnado plot`

Render a template for one or more regions:

```bash
plotnado plot template.yaml --region chr1:1,000,000-1,100,000 --output output.png
```

Useful options:

- `--region` / `-r`: genomic region or gene name; repeat for multiple outputs
- `--output` / `-o`: explicit output path for a single region
- `--format` / `-f`: output format such as `png`, `pdf`, `svg`, or `jpg`
- `--width` / `-w`: override template width in inches
- `--dpi`: output resolution

Examples:

```bash
plotnado plot template.yaml --region chr1:1M-2M
plotnado plot template.yaml --region GNAQ
plotnado plot template.yaml --region chr1:1M-2M --region chr2:5M-6M
plotnado plot template.yaml --region chr1:1,000,000-2,000,000 --output plot.pdf --dpi 300
```

Gene-name resolution requires `genome` to be defined in the template.

## Python bridge

The CLI template format is not separate from the Python API. The same template can be loaded directly in Python:

```python
from plotnado import GenomicFigure

fig = GenomicFigure.from_template("template.yaml")
fig.save("output.png", region="chr1:1,000,000-1,100,000")
```
