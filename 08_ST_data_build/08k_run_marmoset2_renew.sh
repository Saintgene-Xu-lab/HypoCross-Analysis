mkdir renew_20250301_marmoset2
micromamba run -n R4.1 Rscript gen_marmoset2.renew.2025.r /cluster/home/majianyun/share/HD_data/Spatial_Hypothalamus/renewed_ST_data_250227/marmoset2 renew_20250301_marmoset2 /cluster/home/majianyun/share/HD_data/Spatial_Hypothalamus/renewed_ST_data_250227/SpeciesBrainRegion-XSJ_location_v5.marmoset.csv

#### convert to SCT h5ad for transfer
micromamba run -n convert Rscript /cluster/home/majianyun/data/Hypothalamus/seu_to_h5ad.spatial_SCT.r /cluster/home/majianyun/share/HD_data/Spatial_Hypothalamus/mjy_script/renew_20250301_marmoset2/mseu.rds

