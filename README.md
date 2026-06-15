# RUSLE Soil Erosion Model for Rwanda (Google Earth Engine)

A ready-to-run **Revised Universal Soil Loss Equation (RUSLE)** implementation in Google Earth Engine that estimates annual soil loss for any district in Rwanda — or any FAO GAUL level-2 district worldwide — at **10 m resolution**.

> **Change one line — the district name — and run.** No other edits required.

`A = R × K × LS × C × P`  →  soil loss in tonnes per hectare per year (t ha⁻¹ yr⁻¹)

---
Run it live at "https://code.earthengine.google.com/218ca7590cfd028af8e75a9d3e378e31"
## Why this version is scientifically sound

Most quick RUSLE scripts badly **over-estimate** erosion on steep terrain — and Rwanda is "the land of a thousand hills." The usual culprit is the old Wischmeier & Smith slope-steepness factor, which was never validated above ~50% slope and explodes on hillsides. This implementation fixes that and uses defensible, peer-reviewed methods for every factor:

| Factor | Method | Source |
|--------|--------|--------|
| **R** – rainfall erosivity | Linear, `R = 38.5 + 0.35 × P` | Roose-type, calibrated for East Africa |
| **K** – soil erodibility | EPIC equation from sand/silt/clay/organic carbon | Williams (1995) |
| **LS** – topography | **Sine-based S factor** + RUSLE β slope-length | **McCool et al. (1987)** |
| **C** – cover management | `C = (1 − NDVI) / 2`, full-year composite | Durigon et al. (2014) |
| **P** – support practice | Land cover + slope-dependent terracing factor | ESA WorldCover 10 m |

The **LS fix is the key one**: McCool's formula uses the *sine* of the slope angle (not the tangent) specifically to stop steep-slope over-estimation. On a 30° hillside this is roughly one-third of the old polynomial value, which is the difference between forest reading "moderate" instead of an impossible "severe."

---

## Data sources (all free, in the GEE catalog)

- **Rainfall:** CHIRPS pentad precipitation
- **Terrain:** SRTM 30 m DEM
- **Vegetation:** Sentinel-2 Surface Reflectance (Harmonized)
- **Land cover:** ESA WorldCover v200 (10 m)
- **Soil:** OpenLandMap sand / clay / organic-carbon fractions
- **Boundaries:** FAO GAUL simplified, level 2

---

## How to use

1. Open the [Google Earth Engine Code Editor](https://code.earthengine.google.com/).
2. Paste the contents of [`RUSLE_Rwanda_SoilErosion.js`](RUSLE_Rwanda_SoilErosion.js).
3. Edit the parameter block at the top:

   ```javascript
   var COUNTRY  = 'Rwanda';
   var DISTRICT = 'Muhanga';   // <<< change to any district
   var YEAR     = 2024;
   ```

4. Press **Run**. Read the `>>> Soil loss percentiles` print in the Console to sanity-check before trusting the map.
5. Open the **Tasks** tab to export the soil-loss and severity-class GeoTIFFs to Google Drive.

A healthy result for a Rwandan district has a **median in the single digits to mid-teens** and varies smoothly across all six classes. If forest is still red, check the LS and C diagnostics printed in the console.

---

## Severity classes

The first break (10 t ha⁻¹ yr⁻¹) is the conventional **soil-loss tolerance** limit — above it, erosion is considered unsustainable. The high end is split so "severe" does not become a catch-all.

| Class | Range (t ha⁻¹ yr⁻¹) |
|-------|---------------------|
| Very slight | < 5 |
| Tolerable | 5 – 10 |
| Moderate | 10 – 20 |
| High | 20 – 40 |
| Severe | 40 – 80 |
| Extreme | > 80 |

---

## Example result — Muhanga District, 2024

| Class | Area (km²) | Share |
|-------|-----------|-------|
| Very slight (<5) | 154 | ~20% |
| Tolerable (5–10) | 151 | ~20% |
| Moderate (10–20) | 286 | ~37% |
| High (20–40) | 163 | ~21% |
| Severe (40–80) | 7 | ~1% |
| Extreme (>80) | 7 | ~1% |

About 40% of the district sits within tolerable limits, while the small **severe/extreme** pockets — steep cultivated hillsides — are the priority targets for terracing and conservation. These figures are consistent with careful published RUSLE assessments for Rwanda (e.g. Sebeya watershed, 14.7–34.7 t ha⁻¹ yr⁻¹).

---

## Screenshots

<img width="1920" height="991" alt="image" src="https://github.com/user-attachments/assets/7d4ccc07-1ca2-4ac7-a64f-0a3b358b27ec" />

<img width="1920" height="991" alt="image" src="https://github.com/user-attachments/assets/6c6cfa3a-354e-42fa-9e45-aec5bcd8d436" />
<img width="1920" height="991" alt="image" src="https://github.com/user-attachments/assets/4d839737-6521-4ec3-9db7-01bc65dd0585" />


| Soil loss class map | Factor layers | Class area chart |
|---------------------|---------------|------------------|
| `screenshots/soil_loss_class.png` | `screenshots/factors.png` | `screenshots/class_chart.png` |

---

## References

- Renard, K.G. et al. (1997). *Predicting Soil Erosion by Water: RUSLE.* USDA Agriculture Handbook 703.
- McCool, D.K. et al. (1987). Revised slope steepness factor for the USLE. *Transactions of the ASAE* 30(5): 1387–1396.
- Williams, J.R. (1995). The EPIC model. In *Computer Models of Watershed Hydrology.*
- Durigon, V.L. et al. (2014). NDVI-based C factor for tropical regions. *Int. J. Remote Sensing* 35(2): 441–453.
- Karamage, F. et al. (2016). Extent of Cropland and Related Soil Erosion Risk in Rwanda. *IJERPH* 13(5): 503.

---

## Licence

MIT — free to use, adapt and share. Attribution appreciated.

**Author:** Jean De Dieu Niyogisubizo & Furaha Albert Patrick · Rwanda
