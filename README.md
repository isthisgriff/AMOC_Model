# AMOC 3-Box Model — Interactive Simulation & Presentation Dashboard
 
A Python implementation of the **Rooth (1982) three-box model** of the Atlantic Meridional Overturning Circulation (AMOC), built as an interactive presentation tool for exploring how climate forcing factors — Greenland and Antarctic ice melt, atmospheric warming, and trade wind changes — drive the AMOC toward a potential tipping point and collapse.
 
---
 
## What Is the AMOC and Why Does It Matter?
 
The **Atlantic Meridional Overturning Circulation** is one of the most important heat-transport systems on Earth. Often called the "ocean conveyor belt," it moves warm surface water northward from the tropics into the North Atlantic, releasing heat into the atmosphere along the way. This heat transfer is responsible for keeping northwestern Europe significantly warmer than it would otherwise be at those latitudes — London, for example, sits at the same latitude as Calgary.
 
At the northern end of this system, the warm surface water cools, becomes denser, and sinks to the ocean floor, flowing southward as cold deep water before eventually upwelling again in the Southern Ocean. This continuous loop transports roughly **1.3 petawatts** of heat — comparable to the output of a million power stations.
 
**If the AMOC weakens or collapses:**
- Northern Europe experiences significantly colder winters
- Sea levels rise faster along the US East Coast as warm water no longer piles up offshore
- Tropical monsoon patterns shift, affecting rainfall across Africa, Asia, and South America
- Less oxygen reaches the deep ocean, threatening marine ecosystems
 
In recent decades, scientists have observed signs of AMOC weakening and have begun studying the conditions under which it could transition to a permanently weaker — or entirely collapsed — state.
 
---
 
## What This Model Does
 
This project implements the **Rooth (1982) three-box model**, one of the foundational conceptual tools for understanding AMOC dynamics. The model simplifies the Atlantic Ocean into three well-mixed boxes:
 
| Box | Region | Characteristics |
|-----|--------|-----------------|
| **North (N)** | High-latitude North Atlantic | Cold, salty, dense — the sinking region |
| **Tropics (T)** | Low-latitude Atlantic | Warm, salty — high evaporation, heat source |
| **South (S)** | High-latitude South Atlantic | Cold — upwelling region |
 
Each box is described by two state variables: **temperature** (T) and **salinity** (S). These determine each box's density using a linear equation of state:
 
```
ρ = ρ₀ × (1 - α×(T - T_ref) + β×(S - S_ref))
```
 
where α is the thermal expansion coefficient and β is the haline contraction coefficient.
 
The AMOC overturning flow rate `q` is then driven by the density difference between the North and South boxes:
 
```
q = k × (ρ_North - ρ_South)
```
 
When `q > 0`, the North is denser than the South — water sinks in the North and the AMOC runs normally. When `q` approaches 0, the circulation stalls. When `q < 0`, the circulation has reversed.
 
### Forcing Parameters
 
The simulation allows you to apply the following real-world climate forcings:
 
| Parameter | Physical Meaning | Effect on AMOC |
|-----------|-----------------|----------------|
| **Greenland meltwater flux** | Freshwater runoff from the Greenland ice sheet into the North Atlantic | Directly dilutes the North, reducing its density and sinking — fastest path to collapse |
| **Antarctic meltwater flux** | Freshwater runoff into the South Atlantic | Makes the South lighter, initially widening the density gap before indirect weakening takes over — slower acting than Greenland |
| **Atmospheric warming (ΔT)** | Global mean temperature rise above pre-industrial levels | Warms surface water, reducing thermal density contrast; acts alongside meltwater forcing |
| **Wind stress factor** | Strength of trade-wind-driven gyre circulation | Modulates the horizontal exchange of heat and salt between boxes independent of the overturning |
| **Forcing ramp duration** | How gradually the above forcings build up over time | Controls the rate of change — abrupt vs. gradual perturbation |
 
---
 
## Key Scientific Insights This Model Illustrates
 
### 1. Greenland Melt Is More Dangerous on Human Timescales
Greenland meltwater enters the North Atlantic directly at the sinking site. This immediately reduces the density of the water that drives the AMOC, causing rapid weakening. Antarctic meltwater enters the South Atlantic, which initially *increases* the density contrast before the signal slowly propagates around the loop — making it far slower to weaken the AMOC. This is one of Rooth's original key findings.
 
### 2. The AMOC Has a Tipping Point
The model demonstrates that the AMOC is not a system that weakens smoothly and proportionally. Once the density difference between North and South is eroded past a threshold, the circulation can cross into a collapsed state.
 
### 3. Hysteresis — The Circulation May Not Recover
Even if you reduce the forcing after the AMOC has weakened significantly, it does not simply recover. The system exhibits **hysteresis**: the forcing level required to collapse the AMOC is higher than the level required to keep it collapsed. This means that climate intervention after the fact may not restore the circulation to its original state.
 
### 4. Compound Forcing Accelerates Collapse
Atmospheric warming and freshwater forcing together push the system toward collapse faster than either does alone. This is because warming reduces the thermal component of the density difference at the same time that freshwater reduces the saline component.
 
---
 
## Installation & Usage
 
### Requirements
```
numpy
scipy
matplotlib
```
 
All three are standard scientific Python packages. Check if you already have them:
```bash
python -c "import numpy, scipy, matplotlib; print('All dependencies found')"
```
 
If any are missing:
```bash
pip install numpy scipy matplotlib
```
 
### Running the Dashboard
```bash
python amoc_presentation.py
```
 
### Running the Full Analysis (scenarios + hysteresis + dashboard)
```bash
python amoc_3box_model.py
```
> Note: the hysteresis computation runs ~60 sequential simulations and takes approximately 30 seconds.
 
---
 
## Model Limitations
 
This is a **conceptual model** designed for insight and visualization, not operational prediction. Key simplifications to be aware of:
 
- **No spatial structure** — each box is perfectly well-mixed; there are no horizontal or vertical gradients within a box
- **No wind-driven gyre dynamics** — the gyre is parameterized as a diffusive exchange term rather than resolved explicitly
- **No sea ice** — ice formation and melt affect both temperature and salinity in ways not captured here
- **No atmospheric feedbacks** — the atmosphere is represented only as a fixed restoring temperature; in reality, a weakening AMOC would change the atmospheric circulation
- **Linear equation of state** — real seawater density has nonlinear temperature and salinity dependence
- **Fixed box volumes** — the geometry of the ocean does not change with circulation state
 
Despite these simplifications, the Rooth 3-box model correctly captures the key qualitative behaviors of the AMOC: thermohaline-driven overturning, sensitivity to freshwater forcing, bistability, and hysteresis. It remains a standard pedagogical and conceptual tool in physical oceanography.
 
---
 
## References
 
- Rooth, C. (1982). Hydrology and ocean circulation. *Progress in Oceanography*, 11(2), 131–149.
- Rahmstorf, S. (1996). On the freshwater forcing and transport of the Atlantic thermohaline circulation. *Climate Dynamics*, 12(12), 799–811.
- Stommel, H. (1961). Thermohaline convection with two stable regimes of flow. *Tellus*, 13(2), 224–230.
- IPCC AR6 WGI (2021). Chapter 9: Ocean, Cryosphere and Sea Level Change. *FAQ 9.3: Will the Gulf Stream shut down?*
 
---
 
## License
 
This project was developed for academic and educational purposes.# AMOC_Model
