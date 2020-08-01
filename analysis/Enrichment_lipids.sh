## gwas_enrichment.ipynb is downloaded from https://github.com/cumc/bioworkflows/tree/master/fine-mapping/gwas_enrichment.ipynb

cd /home/gw/tmp/18-Jul-2020/
for dataset in hdl.zscore.gz ldl.zscore.gz tc.zscore.gz tg.zscore.gz; do

sos run /home/gw/GIT/bioworkflows/fine-mapping/gwas_enrichment.ipynb range2var_annotation --z-score $dataset --single-annot /home/gw/tmp/18-Jul-2020/lipids_bed.txt --annotation_dir /home/gw/Documents/annotations/bed/ --cwd /home/gw/tmp/18-Jul-2020/ -j 10
sos run /home/gw/GIT/bioworkflows/fine-mapping/gwas_enrichment.ipynb enrichment --z-score $dataset --single-annot lipids_bed.txt --annotation_dir /home/gw/Documents/annotations/bed/ --cwd /home/gw/tmp/18-Jul-2020/ -j 5

done

cd /home/gw/tmp/31-Jul-2020/
for dataset in original_z_mvp.ldl.zscore.gz original_z_mvp.tc.zscore.gz original_z_mvp.tg.zscore.gz original_z_mvp.hdl.zscore.gz; do

sos run /home/gw/GIT/bioworkflows/fine-mapping/gwas_enrichment.ipynb range2var_annotation --z-score $dataset --single-annot /home/gw/tmp/18-Jul-2020/lipids_bed.txt --annotation_dir /home/gw/Documents/annotations/bed/ --cwd /home/gw/tmp/31-Jul-2020/ -j 10
sos run /home/gw/GIT/bioworkflows/fine-mapping/gwas_enrichment.ipynb enrichment --z-score $dataset --single-annot /home/gw/tmp/18-Jul-2020/lipids_bed.txt --annotation_dir /home/gw/Documents/annotations/bed/ --cwd /home/gw/tmp/31-Jul-2020/ -j 5

done