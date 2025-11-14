# Single Cell FRET Extraction

This repository provides a robust, track-aware single-cell FRET extraction pipeline designed for Cellpose-segmented timelapse data.  
The extractor uses stable track IDs, per-frame donor/acceptor images, Sauvola gating (optional), ring-based background subtraction, and bleed-through correction to produce clean, time-aligned single-cell FRET traces.

The script supports:

- Ring-based background subtraction  
- Optional Sauvola foreground gating  
- Gaussian smoothing  
- Bleed-through correction  
- Track-aware FRET computation  
- Per-track Excel sheets  
- Per-track PNG plots  
- Full parameter control via command-line arguments  

---

## 1. Folder Structure

Your data should look like:
```
FOV_XX/
├── Segmented/
│ ├── frame_000_cp_masks.tif
│ ├── frame_001_cp_masks.tif
│ └── ...
├── Donor/
│ ├── frame_000.tif
│ └── ...
├── Acceptor/
│ ├── frame_000.tif
│ └── ...
└── tracks_v3.7_sticky_new.csv
```
