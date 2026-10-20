mkdir renew_20250301_{mouse,macaque,marmoset}
micromamba run -n R4.1 Rscript gen_mouse.renew.2025.r /cluster/home/majianyun/share/HD_data/Spatial_Hypothalamus/renewed_ST_data_250227/Mouse250313 renew_20250301_mouse /cluster/home/majianyun/share/HD_data/Spatial_Hypothalamus/renewed_ST_data_250227/SpeciesBrainRegion-XSJ_location_v5.mouse.csv
micromamba run -n R4.1 Rscript gen_macaque.renew.2025.r /cluster/home/majianyun/share/HD_data/Spatial_Hypothalamus/renewed_ST_data_250227/Macaque250327 renew_20250301_macaque /cluster/home/majianyun/share/HD_data/Spatial_Hypothalamus/renewed_ST_data_250227/SpeciesBrainRegion-XSJ_location_v5.macaque.csv
micromamba run -n R4.1 Rscript gen_marmoset.renew.2025.r /cluster/home/majianyun/share/HD_data/Spatial_Hypothalamus/renewed_ST_data_250227/Marmoset250220 renew_20250301_marmoset /cluster/home/majianyun/share/HD_data/Spatial_Hypothalamus/renewed_ST_data_250227/SpeciesBrainRegion-XSJ_location_v5.marmoset.csv

#### convert to SCT h5ad for transfer
micromamba run -n convert Rscript /cluster/home/majianyun/data/Hypothalamus/seu_to_h5ad.spatial_SCT.r /cluster/home/majianyun/share/HD_data/Spatial_Hypothalamus/mjy_script/renew_20250301_*/mseu.rds

