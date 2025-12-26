# RGB-UAV - Canopy Gaps Detection Pipeline 

[DOI]

**R pipeline for automated detection of canopy gaps using UAV DSM difference**  
*Simonetti et al. (2026)*


## 🎯 **What it does**

Detects forest canopy gaps through DSM change detection with these steps:
1. **DSM change** (imagedate_time2 - imagedate_time2) to identify height loss
2. **Focal median filter** (99×99 cells) to correct slow vertical warping
3. **Threshold** (< -5m height loss) to isolate gaps
4. **Geometric filtering** (area > 5m², area/perimeter ratio > 0.6)


## 🔗 **Related Work**

**Araujo, R. F. et al. (2021)**  
*Strong temporal variation in treefall and branchfall rates in a tropical forest is related to extreme rainfall: results from 5 years of monthly drone data for a 50 ha plot.*  
**Biogeosciences**, 18, 6517–6531.  
[https://doi.org/10.5194/bg-18-6517-2021](https://doi.org/10.5194/bg-18-6517-2021)

**Simonetti, A. et al. (2023)**  
*Canopy gaps and associated losses of biomass – combining UAV imagery and field data in a central Amazon forest.*  
**Biogeosciences**, 20, 3651–3666.  
[https://doi.org/10.5194/bg-20-3651-2023](https://doi.org/10.5194/bg-20-3651-2023)
