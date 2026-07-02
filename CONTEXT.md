# ESPI Analysis

Shared glossary for the ESPI single-cell analysis pipeline and paper-focused outputs.

## Language

**p27/CDKN1B**:
A cell-cycle regulator whose loss defines the p27CKO context for this analysis. Use `p27` in prose and `CDKN1B` when naming the gene.
_Avoid_: p27Kip1, Cdkn1b

**p27CKO**:
The p27 conditional knockout condition used as the control condition in this analysis.
_Avoid_: p27 knockout, p27 loss condition

**E-Stim**:
The electrical stimulation condition compared against p27CKO control samples.
_Avoid_: electrostimulation, stimulation

**BrdU**:
A proliferation readout added after E-Stim to mark cells that entered DNA synthesis during the labeling window.
_Avoid_: EdU, proliferation dye

**NPC**:
Neural progenitor cell identity used when interpreting marker evidence and focused proliferation hypotheses.
_Avoid_: neural stem cell, progenitor cell

**OPC**:
Oligodendrocyte precursor cell identity used when interpreting marker evidence.
_Avoid_: oligodendrocyte progenitor, pre-oligodendrocyte

**Interpreted identity**:
A biological cell identity assigned from marker evidence, such as NPC, neuron, oligodendrocyte, OPC, astrocyte, microglia, or endothelial cell.
_Avoid_: cell type label, cluster identity

**Normalization branch**:
A Seurat object produced from the same filtered cells and genes using one normalization method, either default log1p normalization or PFlog normalization. Each normalization branch owns its canonical `pca` reduction.
_Avoid_: normalization layer, normalization method

**Candidate clustering**:
A clustering result created to compare PC dimensions or Leiden resolution choices within a normalization branch.
_Avoid_: clustering trial, clustering run

**Chosen clustering**:
A candidate clustering selected for downstream interpretation within a normalization branch, including UMAP annotation, heatmaps, pseudobulk summaries, and focused tests.
_Avoid_: final clustering, selected clustering

**Pseudobulk sample**:
An aggregate expression profile whose replicate unit is a Mouse × Condition sample. Identity-stratified pseudobulk adds interpreted identity as a stratum, but the replicate unit remains Mouse × Condition.
_Avoid_: bulk sample, cell-level sample

**Focused test**:
A named condition-effect test chosen before running differential expression or enrichment. Focused tests are separate from exploratory clustering, marker discovery, and heatmap interpretation.
_Avoid_: targeted test, exploratory test
