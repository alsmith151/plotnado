# PlotNado

<div class="pn-hero">
    <div class="pn-hero__copy">
        <p class="pn-eyebrow">Genome browser figures without the browser</p>
        <h1>Publication-ready genomic tracks in Python.</h1>
        <p class="pn-lead">PlotNado composes BigWig signals, genes, peaks, interactions, and overlays into clean genome-browser style figures with a high-level API and a file-driven CLI.</p>
        <div class="pn-actions">
            <a class="pn-button pn-button--primary" href="quickstart/">Start with Quick Start</a>
            <a class="pn-button pn-button--secondary" href="track_catalog/">Browse track catalog</a>
            <a class="pn-button pn-button--ghost" href="cli/">Explore the CLI</a>
        </div>
        <div class="pn-stats">
            <div class="pn-stat">
                <strong>Python API</strong>
                <span>Chain one track method at a time.</span>
            </div>
            <div class="pn-stat">
                <strong>CLI + YAML</strong>
                <span>Infer, edit, validate, and render templates.</span>
            </div>
            <div class="pn-stat">
                <strong>Track-first design</strong>
                <span>Signals, annotations, overlays, and matrices.</span>
            </div>
        </div>
    </div>
    <div class="pn-hero__visual">
        <img src="images/examples/basic_figure.png" alt="PlotNado example figure" />
        <p class="pn-caption">A composed PlotNado figure with signal, annotation, and structural tracks.</p>
    </div>
</div>

<div class="pn-section">
    <h2>Pick your entry point</h2>
    <p>Most readers need one of three routes: get a figure rendered quickly, understand which track type to use, or tune the figure until it is publication-ready.</p>
</div>

<div class="pn-card-grid">
    <div class="pn-card" markdown="1">
        <p class="pn-kicker">Start Here</p>
        <h3>Install and render your first figure</h3>
        <p>Use the Quick Start when you want a working figure with minimal setup and clear defaults.</p>
        <div class="pn-linklist">
            <a href="installation/">Installation</a>
            <a href="quickstart/">Quick Start</a>
            <a href="quickstart_tracks/">Track construction patterns</a>
        </div>
    </div>
    <div class="pn-card" markdown="1">
        <p class="pn-kicker">Build With Confidence</p>
        <h3>Choose the right track and input model</h3>
        <p>Learn which track to reach for, what each one expects, and where to find concrete examples.</p>
        <div class="pn-linklist">
            <a href="track_catalog/">Track Catalog</a>
            <a href="data_inputs/">Data Inputs</a>
            <a href="example_coverage/">Example Coverage</a>
        </div>
    </div>
    <div class="pn-card" markdown="1">
        <p class="pn-kicker">Tune The Result</p>
        <h3>Dial in styling, scaling, and layout</h3>
        <p>Use the aesthetics and recipes guides when the figure works but does not yet communicate clearly.</p>
        <div class="pn-linklist">
            <a href="aesthetics/">Aesthetics</a>
            <a href="recipes/">Recipes</a>
            <a href="troubleshooting/">Troubleshooting</a>
        </div>
    </div>
</div>

<div class="pn-section">
    <h2>Two workflows, one rendering model</h2>
    <p>PlotNado offers a fluent Python API for notebooks and pipelines, plus a CLI that turns genomic inputs into editable YAML templates.</p>
</div>

<div class="pn-card-grid">
    <div class="pn-card">
        <p class="pn-kicker">Python API</p>
        <h3>Compose tracks directly in code</h3>
        <p>Add tracks one method at a time and save the result in a single region-specific render step.</p>

```python
from plotnado import GenomicFigure

fig = GenomicFigure(theme="publication")
fig.scalebar()
fig.axis()
fig.genes("hg38")
fig.bigwig("signal.bw", title="ChIP signal", style="fill")
fig.save("output.png", region="chr1:1,010,000-1,080,000")
```
    </div>
    <div class="pn-card">
        <p class="pn-kicker">CLI + YAML</p>
        <h3>Infer a template, then edit and render it</h3>
        <p>This route is useful when the inputs already exist on disk and you want a reproducible, file-driven workflow.</p>

```bash
plotnado init *.bw peaks.narrowpeak --auto --output template.yaml
plotnado validate template.yaml
plotnado plot template.yaml --region chr1:1,000,000-1,100,000 --output out.png
```
    </div>
    <div class="pn-card">
        <p class="pn-kicker">Need Options Fast?</p>
        <h3>Inspect the runtime metadata</h3>
        <p>Every track exposes discoverable options, so you do not have to guess field names or dig through source first.</p>

```python
from plotnado import GenomicFigure

GenomicFigure.track_options("overlay")
GenomicFigure.track_options_markdown("bigwig")
```
    </div>
</div>

<div class="pn-section">
    <h2>Rendered examples</h2>
    <p>The docs include plotted outputs, not just code snippets, so you can see what each configuration actually produces.</p>
</div>

<div class="pn-showcase">
    <div class="pn-showcase__item">
        <img src="images/quickstart.png" alt="Quick start output" />
        <h3>First figure</h3>
        <p>Start from a minimal stack: scale bar, axis, genes, and one signal track.</p>
    </div>
    <div class="pn-showcase__item">
        <img src="images/aesthetics/overlay_autoscale.png" alt="Overlay autoscale output" />
        <h3>Overlay + autoscale</h3>
        <p>Overlay multiple signals in one panel while keeping y-scaling explicit and readable.</p>
    </div>
    <div class="pn-showcase__item">
        <img src="images/examples/recipe_theme_labels.png" alt="Theme and labels output" />
        <h3>Theme-driven polish</h3>
        <p>Refine color, labels, and spacing once the track composition is correct.</p>
    </div>
</div>

<div class="pn-section">
    <h2>Where to go next</h2>
    <p class="pn-inline-note">If the figure is blank or scaling looks wrong, go straight to <a href="troubleshooting/">Troubleshooting</a>. If you are deciding between track types, start with the <a href="track_catalog/">Track Catalog</a>. If you want exact parameters, jump to <a href="reference/">Reference</a> and <a href="api_reference/">API Reference</a>.</p>
</div>
