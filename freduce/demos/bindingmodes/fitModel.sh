# Step 1a: train ntgaStcan type symmetric PSAM on original PBM binding data
java -Xmx1500M FeatureReduce bindingMode.init -c 0 3 2 -l zhu2009.9bp -unifiedSeeds ATGACTCAT -correspond "1 2 3 4 ~5 ~4 ~3 ~2 ~1" -i gcn4.dream.v11.clean.txt -savePsam xml results/zhu2009.9bp.gcn4.v11.psam.9nt.xml -displayMotifs No

# Step 1b: train ntgaCGtcan type symmetric PSAM on original PBM binding data
java -Xmx1500M FeatureReduce bindingMode.init -c 0 3 2 -l zhu2009.10bp -unifiedSeeds ATGACGTCAT -correspond "1 2 3 4 5 ~5 ~4 ~3 ~2 ~1" -i gcn4.dream.v11.clean.txt -savePsam xml results/zhu2009.10bp.gcn4.v11.psam.10nt.xml -displayMotifs No

# Step 2a: retrain first PSAM on residuals of second PSAM
java -Xmx1500M FeatureReduce bindingMode.init -c 0 3 2 -l zhu2009.9bp.resid -unifiedSeeds ATGACTCAT -correspond "1 2 3 4 ~5 ~4 ~3 ~2 ~1" -i gcn4.dream.v11.clean.txt -residuals xml results/zhu2009.10bp.gcn4.v11.psam.10nt.xml -savePsam xml results/zhu2009.9bp.gcn4.final.v11.psam.9nt.xml -displayMotifs No

# Step 2a: retrain second PSAM on residuals of first PSAM
java -Xmx1500M FeatureReduce bindingMode.init -c 0 3 2 -l zhu2009.10bp.resid -unifiedSeeds ATGACGTCAT -correspond "1 2 3 4 5 ~5 ~4 ~3 ~2 ~1" -i gcn4.dream.v11.clean.txt -residuals xml results/zhu2009.9bp.gcn4.v11.psam.9nt.xml  -savePsam xml results/zhu2009.10bp.gcn4.final.v11.psam.10nt.xml -displayMotifs No

# Step 3: fit model in which the PSAMs compete
java -Xmx1500M FeatureReduce bindingModesCompare.init -c 0 3 2 -l zhu2009.compare -compareMotifs results/zhu2009.9bp.gcn4.final.v11.psam.9nt.xml results/zhu2009.10bp.gcn4.final.v11.psam.10nt.xml  -comparePosWeights results/zhu2009.9bp.gcn4.psam.9nt.positionalBias.table results/zhu2009.10bp.gcn4.psam.10nt.positionalBias.table -i gcn4.dream.v11.clean.txt -displayMotifs No
