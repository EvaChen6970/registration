# Interpretation

The manifest is the source of truth for attempted pairs; a TIFF on disk is not proof of acceptance. Distinguish `registered`, `registered_fallback`, `failed_qc`, and `failed_visual_qc`. Review overlays and transform statistics together. SIFT inliers and thumbnail reprojection error are not calibrated physical errors.

Tissue coverage is a screening metric only. Thresholds vary by stain, scanner, and background. Keep the CSV and make filtering decisions explicitly; scripts never delete files.

Nucleus distances are nearest-neighbor diagnostics, not a registration score by themselves: segmentation quality, cell density, tissue composition, and cross-modal cell detection affect them. Use matched tissue masks and report distributions, not just the median. Transform direction is HE coordinates to CODEX coordinates.

## Nuclei instance overlap

Segment nuclei on the registered HE image, so the HE instance-label TIFF is already in the SP pixel grid. Compare it with the SP/CODEX instance-label TIFF using `nuclei_overlap.py`. The tool requires identical 2D shapes, treats label 0 as background, reports one-to-one object matches ranked by IoU, and writes per-pair CSV plus summary JSON. The IoU threshold is configurable; report sensitivity across thresholds (for example 0.1, 0.25, 0.5). This measures segmentation agreement as well as registration, so inspect overlays and segmentation quality too.
