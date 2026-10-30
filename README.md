# Bone marrow CODEX / H&E registration toolkit

All commands run on jhcpu1. The raw dataset is read-only; derived outputs go to the dataset's `registration/` directory unless an output path is supplied.

## Dataset audit

`python src/audit_dataset.py --root /jhcnas6/Spatial_Proteomics_Temp/Bone_Marrow_CODEX`

## Register one HE/SP pair

The generalized implementation takes explicit image paths. HE is the moving image; SP is fixed. The SP channel is used only for feature matching (default channel 0, usually DAPI for CODEX).

```bash
python src/register_pair.py \
  --he /path/to/sample_HE.qptiff \
  --sp /path/to/sample_SP.qptiff \
  --out-dir /path/to/output/sample \
  --sp-channel 0
```

Inputs: one H&E WSI and one spatial-proteomics WSI readable as pyramidal TIFF/QPTIFF. Outputs: registered tiled pyramidal BigTIFF in SP coordinates, JSON affine transform (HE pixels to SP pixels), and a thumbnail contact sheet/overlay for review. Existing registered output is protected; use `--overwrite` to replace it. The raw inputs are never changed. Requires OpenCV, tifffile, numpy, and the `vips` executable.

For this dataset, example paths are:

```bash
python /ssd2/chenyuxi/registration_sp/src/register_pair.py \
  --he /jhcnas6/Spatial_Proteomics_Temp/Bone_Marrow_CODEX/Normal_H\&E_Images/H10_HE.qptiff \
  --sp /jhcnas6/Spatial_Proteomics_Temp/Bone_Marrow_CODEX/Normal_CODEX_Images/H10_Scan1.qptiff \
  --out-dir /jhcnas6/Spatial_Proteomics_Temp/Bone_Marrow_CODEX/registration/manual/H10
```

## Quality control

`python src/tissue_qc.py --root DATASET --out registration/tissue_qc.csv` creates tissue coverage and blur/foreground metrics from pyramid thumbnails. Review thresholds before filtering; it never deletes or modifies source images.

`python src/nuclei_distance.py --help` compares centroids from supplied HE and DAPI label masks after applying a saved transform. Reports nearest-neighbor distances in pixels and microns when pixel sizes are provided. Segmentation masks must be in their respective image coordinate systems.

## Centroid distance QC without segmentation masks

When no SP segmentation mask is available, use `centroid_qc.py`. It reads the registered HE image and the original SP image, extracts DAPI nuclei from SP channel 0 by thresholding/connected components, extracts blue/purple nuclei from HE, and compares each HE centroid with its nearest DAPI centroid. Both images are analyzed at their coarsest pyramid level. This is a screening metric, not a replacement for manual overlay review.

```bash
ENV=/ssd2/chenyuxi/conda_envs/aslide_env/bin/python
ROOT=/jhcnas6/Spatial_Proteomics_Temp/Bone_Marrow_CODEX
$ENV src/centroid_qc.py \
  --he-registered $ROOT/registration/registered_h\&e/Normal/H10_HE_to_CODEX.ome.tif \
  --sp $ROOT/Normal_CODEX_Images/H10_Scan1.qptiff \
  --out-dir $ROOT/registration/centroid_qc/H10 \
  --sp-channel 0 \
  --min-area 8 --max-area 2000 --max-distance 25
```

Outputs are `centroid_pairs.csv` and `summary.json`. Tune `--min-area`, `--max-area`, and `--max-distance` after reviewing detection overlays.

## Visual review

`python src/view_registration.py --help` creates an interactive napari view for a WSI or a selected patch. Use `--sample H10 --level 3 --x 12000 --y 8000 --size 2048` for a patch. The registered H&E and CODEX DAPI are shown in fixed CODEX coordinates.

## Dependencies

Use the existing registration environment with OpenCV, tifffile, numpy, pandas, scikit-image, scipy, napari, and pyvips/libvips. The registration script requires the `vips` executable. See `docs/` for coordinate conventions and QC limitations.
