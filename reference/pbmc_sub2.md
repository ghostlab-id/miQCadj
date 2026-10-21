# Example PBMC dataset for miQCadj

A subset of the `pbmcsca` dataset from
[SeuratData](https://github.com/satijalab/seurat-data) containing two
sample–technology combinations: pbmc1 generated using Smart-seq2 and
pbmc2 generated with CEL-Seq2.

## Format

A `Seurat` object containing 33,694 features and 526 cells. The RNA
assay contains raw counts, and the cell metadata include:

- `orig.ident`:

  Original sample identity.

- `nCount_RNA`:

  Total RNA counts per cell.

- `nFeature_RNA`:

  Number of detected RNA features per cell.

- `percent.mito`:

  Original proportion of counts assigned to mitochondrial genes.

- `percent.mt`:

  Recomputed proportion of counts assigned to mitochondrial genes.

- `Cluster`:

  Original cluster annotation.

- `CellType`:

  Original cell-type annotation.

- `sample`:

  Sample identifier.

- `tech`:

  Sequencing technology.

- `sample_tech`:

  Combined sample and technology identifier.

## Source

Derived from the `pbmcsca` dataset distributed through the
[SeuratData](https://github.com/satijalab/seurat-data) package.

## Details

The object contains raw RNA counts and cell-level metadata used to
demonstrate miQC-based quality-control modelling and filtering.
