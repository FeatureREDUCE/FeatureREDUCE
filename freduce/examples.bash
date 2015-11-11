#!/bin/bash

### Make sure that Xvfb is running
### or start it here
# bash -v start.xvfb.bash

### Set the CLASSPATH and DISPLAY
# source .envrc

##################################################
# Training the models
##################################################

### Train all the models (TF_1 thru TF_66) and also train the All-Kmers model
# java -Xmx5000M FeatureReduce dream.init -kmer -l dream -i DREAM5_PBM_Data_Needed_For_Predictions.txt 2>&1 | tee dream.training.all.out


### Train the model for TF_1 including the All-Kmers model
java -Xmx5000M FeatureReduce dream.init -ids TF_1 -kmer -l dream -i DREAM5_PBM_Data_Needed_For_Predictions.txt 2>&1 | tee dream.training.TF_1.out


##################################################
# Predicting Probe intensities with the models
##################################################

### Predict all the probe intensities (TF_1 thru TF_66) and use the All-Kmers Models
# java -Xmx5000M FeatureReduce dream.init -kmer -l dream -a DREAM5_PBM_TeamName_Predictions.txt -o DREAM5_PBM_FeatureREDUCE_Predictions.all.txt  2>&1 | tee dream.predictions.all.out

### Predict the probe intensities for TF_1 and use the All-Kmers Model
java -Xmx5000M FeatureReduce dream.init -ids TF_1 -kmer -l dream -a DREAM5_PBM_TeamName_Predictions.txt -o DREAM5_PBM_FeatureREDUCE_Predictions.TF_1.txt  2>&1 | tee dream.predictions.TF_1.out

### Predict the probe intensities for TF_1 and DON'T use the All-Kmers Model
# java -Xmx5000M FeatureReduce dream.init -ids TF_1 -l dream -a DREAM5_PBM_TeamName_Predictions.txt -o DREAM5_PBM_FeatureREDUCE_Predictions.TF_1.txt  2>&1 | tee dream.predictions.TF_1.out

###################################################################################
# Training the models and Predicting Probe intensities with the models together
###################################################################################

### Train all the models (TF_1 thru TF_66) and also train the All-Kmers model, then predict all the probe intensities (TF_1 thru TF_66) and use the All-Kmers Models in the predictions
# java -Xmx5000M FeatureReduce dream.init -kmer -l dream -i DREAM5_PBM_Data_Needed_For_Predictions.txt -a DREAM5_PBM_TeamName_Predictions.txt -o DREAM5_PBM_FeatureREDUCE_Predictions.all.txt 2>&1 | tee dream.training.predicting.all.out


### Train the model for TF_1 including the All-Kmers model, then predict the probe intensities for TF_1 and use the All-Kmers Model in the predictions
# java -Xmx5000M FeatureReduce dream.init -ids TF_1 -kmer -l dream -i DREAM5_PBM_Data_Needed_For_Predictions.txt -a DREAM5_PBM_TeamName_Predictions.txt -o DREAM5_PBM_FeatureREDUCE_Predictions.TF_1.txt 2>&1 | tee dream.training.predicting.TF_1.out


##################################################
# Predicting (non-PBM) Sequence Affinities
##################################################

### Predict affinities for non-PBM DNA for all the IDs
# java -Xmx5000M FeatureReduce dream.init -noPosBias -l dream -a DREAM5_PBM_TeamName_Predictions.txt -o DREAM5_PBM_FeatureREDUCE_Predictions.all.txt  2>&1 | tee dream.predictions.all.out

### Predict affinities for non-PBM DNA for just TF_1
# java -Xmx5000M FeatureReduce dream.init -ids TF_1 -noPosBias -l dream -a DREAM5_PBM_TeamName_Predictions.txt -o DREAM5_PBM_FeatureREDUCE_Predictions.TF_1.txt  2>&1 | tee dream.predictions.TF_1.out

###################################################################################
# Getting the corresponding K-mer Affinities from an FSAM
###################################################################################

### Get the affinity sphere for all 10-mers with a relative affinity >= 0.01 of the highest affinity 10-mer.
# java -Xmx1000M FeatureReduce -fsam results.dream/dream.TF_1.fsam.10nt.ser -affinitySphere 0 0.01 TF_1.fsam.10nt.affinitySphere.01.table


###################################################################################
# Loading an FSAM model for viewing
###################################################################################

### Just load an FSAM model for viewing
# java FeatureReduce -fsam results.dream/dream.TF_1.fsam.10nt.ser

### Just dump the FSAM parameters without showing the Logo
# java FeatureReduce -fsam results.dream/dream.TF_1.fsam.10nt.ser -displayMotifs No

###################################################################################
# Loading a PSAM model for viewing
###################################################################################

### Just load an PSAM model for viewing the logo
# java FeatureReduce -psam ser results.dream/dream.TF_1.psam.10nt.ser
# java FeatureReduce -psam xml results.dream/dream.TF_1.psam.10nt.xml
