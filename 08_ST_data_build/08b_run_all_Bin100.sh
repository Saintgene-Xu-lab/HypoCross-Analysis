mkdir renew_Bin100_20250414_{mouse,macaque,marmoset}
micromamba run -n R4.1 Rscript /cluster/home/majianyun/share/HD_data/Spatial_Hypothalamus/mjy_script/gen_mouse_renew.bin50.r /cluster/home/majianyun/share/HD_data/Spatial_Hypothalamus/renewed_bin100_250227/mouse renew_Bin100_20250414_mouse /cluster/home/majianyun/share/HD_data/Spatial_Hypothalamus/renewed_ST_data_250227/SpeciesBrainRegion-XSJ_location_v5.mouse.csv
micromamba run -n R4.1 Rscript /cluster/home/majianyun/share/HD_data/Spatial_Hypothalamus/mjy_script/gen_macaque_renew.bin50.r /cluster/home/majianyun/share/HD_data/Spatial_Hypothalamus/renewed_bin100_250227/macaque renew_Bin100_20250414_macaque /cluster/home/majianyun/share/HD_data/Spatial_Hypothalamus/renewed_ST_data_250227/SpeciesBrainRegion-XSJ_location_v5.macaque.csv
micromamba run -n R4.1 Rscript gen_marmoset_renew.bin50.r /cluster/home/majianyun/share/HD_data/Spatial_Hypothalamus/renewed_bin100_250227/marmoset renew_Bin100_20250414_marmoset /cluster/home/majianyun/share/HD_data/Spatial_Hypothalamus/renewed_ST_data_250227/SpeciesBrainRegion-XSJ_location_v5.marmoset.csv

#### convert to SCT h5ad for transfer
# micromamba run -n convert Rscript /cluster/home/majianyun/data/Hypothalamus/seu_to_h5ad.spatial_SCT.r /cluster/home/majianyun/share/HD_data/Spatial_Hypothalamus/mjy_script/renew_Bin100_20250414_*/mseu.rds

