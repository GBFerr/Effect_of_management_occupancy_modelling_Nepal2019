# Effect_of_management_occupancy_modelling_Nepal2019
Code used in the paper 'Wildlife response to management regime and habitat loss in the Terai Arc Landscape of Nepal'; Ferreira et al 2022

Part of the Biome Health Project

PRE ANALYSIS
File: LandCoverAND... - extract % land cover from raster file to be used as covariate in the analysis; explore covariates - other variables were obtained from field forms and QGIS

File: building_spp_matrices... - generate site by date matrices of species occurrence. This is the detection/non-detection data used in the multi-species occupancy models


ANALYSIS:
File: MSOM_ManagementEffect... - implements the multi-species occupancy model to estimate the effect of management (Figs. 2-3)

File: MSOM_HabitatLoss... - implements the multi-species occupancy model to estimate the interaction between habitat loss caused by agriculture and management regime (Fig. 4)
