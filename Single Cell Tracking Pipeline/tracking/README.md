# Single Cell Tracking

This repository provides a robust, ID-stable cell tracking algorithm designed for Cellpose-generated segmentation masks.  
The tracker is optimized for long-term single-cell imaging with shape jitter, merges/splits, moderate drift, and occasional detection loss.

The script supports:

- Stable “sticky” ID assignment  
- Hungarian matching with IoU-dominant cost  
- Constant-velocity prediction  
- Anti-takeover logic  
- Reviving recently lost tracks  
- Hard locking of prior owners  
- Optional overlay PNGs and MP4 videos  
- Optional per-frame TIFFs storing track ID labels  
- Full parameter control via command line arguments  

---

## 1. Folder Structure

Your data should look like:

```
FOV_XX/
├── Segmented/
│ ├── frame_000_cp_masks.tif
│ ├── frame_001_cp_masks.tif
│ └── ...
├── Donor/ (optional for FRET pipeline)
├── Acceptor/ (optional for FRET pipeline)
└── raw frames (frame_000.tif, frame_001.tif, ...)
```

This runs the tracker with safe default parameters.

## 2. Override parameters

Example: allow faster movement + higher IoU strictness:

```
python track_render_label_overlay.py \
    --mask-dir /path/to/FOV_40/Segmented \
    --search-range 90 \
    --beta 0.90
```
Example: enable track TIFF output:

```
python track_render_label_overlay.py \
    --mask-dir /path/to/FOV_12/Segmented \
    --write-track-tiffs
```

Show all available options:

```
python track_render_label_overlay.py --help
```
## 3. Outputs

The script generates:
```
tracks_<name>.csv
tracks_overlay_slow_<name>.mp4
overlays_<name>/   (per-frame PNG overlays)
tracks_tiff_<name>/ (optional TIFFs with track IDs)
```
The CSV includes:

- `frame`
- `track_id`
- `cx`, `cy` (centroid)
- `area`
- `label` (Cellpose label)

## 4. Parameters and Suggested Ranges (Compact)

Below is a concise summary of all tracking parameters with typical workable ranges. Adjust them depending on drift, segmentation noise, and cell behavior.

**Motion / Matching**

- `search-range: Max movement (px) per frame. Typical 20–100 (default 70)`

- `alpha: Distance weight. Typical 0.1–1.0`

- `beta: IoU weight (higher = more stable IDs). Typical 0.5–1.0`

- `gamma: Area-change penalty. Typical 0.01–0.2`

- `max-cost: Max allowed matching cost. Typical 1.0–3.0`

**Memory / Revive**

- `memory: Frames a track can vanish but stay active. Typical 3–15`

- `revive-gap: How far back ended tracks can be revived. Typical 20–100`

**Morphology**

- `dilate-radius: Dilation radius for IoU tolerance. Typical 5–15 px`

- `min-area: Minimum object size. Typical 30–200 px`

**Anti-Takeover & Revive Strictness**

- `revive-iou-min: Minimum IoU to revive/reuse ID. Typical 0.05–0.3`

- `no-takeover-recent: Protect last owner N frames. Typical 3–15`

- `mutual-iou-tol: Tolerance for mutual best-match. Typical 1e-6–1e-3`

- `no-enforce-mutual-iou: Disable mutual IoU rule (not recommended)`

**Hard-Lock (pre-assignment)**

- `lock-iou: IoU threshold for hard lock. Typical 0.02–0.1`

- `lock-dist: Max predicted distance for lock. Typical 30–100 px`

**Rendering**

- `fps: MP4 playback speed. Typical 1–10 FPS`

- `repeat: Repeat each frame (slower playback). Typical 1–5`

- `write-track-tiffs: Save per-frame track-ID TIFFs`

## 5. Recommended Presets

**Slow-moving cells:**
- `search-range 50`
- `beta 0.85`
- `lock-iou 0.08`
- `lock-dist 40`

**Fast-moving or noisy segmentation:**
- `search-range 100`
- `beta 0.65`
- `revive-gap 80`
- `no-takeover-recent 4`

**High-precision tracking (strict):**
- `beta 1.0`
- `revive-iou-min 0.15`
- `lock-iou 0.07`
- `max-cost 1.2`

## 6. Citation

If you use this code in a publication, please cite the repository or description in your methods section.
