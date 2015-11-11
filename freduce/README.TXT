************************************************************************
                          FeatureREDUCE 1.0
                             Todd Riley
                           Bussemaker Lab
                         Columbia University
************************************************************************


Welcome to FeatureREDUCE! This tool builds protein-DNA affinity models
from occupancy assays like PBM, SELEX-seq, and MITOMI. For more info
about the rich functionality of FeatureREDUCE (and a pleasant read)
please refer to the User Manual.

*******************************************************
1.0 Installation
*******************************************************

FeatureREDUCE is a java application that interfaces to R using the
rJava/JRI libraries.  Some regressions are performed in R and others in
java.  FeatureREDUCE also makes use of the BioJava, GNU Trove,
jFreeChart, Apache Commons, and iText classes.

The Java dependencies are:
	JVM >= 1.5
	rJava/JRI 0.6-3
	GNU Trove 3.0.2
	Apache Commons
	BioJava 1.5
	jFreeChart 1.0.13
	iText 5.1.0

NOTE: ALL THE JAVA LIBRARIES ARE INCLUDED IN THIS INSTALL EXCEPT FOR THE JVM!
      THEY ARE LOCATED IN THE freduce/bin DIRECTORY

The R dependencies are:
	R  2.10.1
	Mass package
	Matrix package
	nnls package
	rJava package 0.6-3

**************************
1.1 - Installing rJava/JRI
**************************

Installing rJava/JRI can be a pain. In the R.packages directory are already
precompiled, linux 64-bit rJava/JRI binaries that you can unzip into your R library
directory, as well as the rJava/JRI source package. Read the README.TXT in the
R.packages directory for further info in how to install and test the rJava/JRI
package.

*****************************
1.2 - Setting the Environment
*****************************

Execute "source .envrc" in a bash shell to add environment variables and augment
the CLASSPATH with the necessary jar files for FeatureREDUCE. Note that you will
probably need to change the R environment variables to match your particular
installation:

Script .envrc:

# make sure that the necessary R environment variables are set
# CHANGE THE R ENVIRONMENT VARIABLES TO MATCH YOUR R INSTALL!!
R_HOME=/usr/lib/R
export R_HOME

RJAVA_HOME=/usr/local/lib/R
export RJAVA_HOME

R_SHARE_DIR=/usr/share/R/share
export R_SHARE_DIR

LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$RJAVA_HOME/site-library/rJava/jri
export LD_LIBRARY_PATH

# augment the CLASSPATH with all the FeatureREDUCE dependencies
CLASSPATH=$CLASSPATH:$HOME/freduce/bin:$HOME/freduce/bin/biojava.jar
:$HOME/freduce/bin/Biojava.suppl.v1.01.jar:$HOME/freduce/bin/bytecode.jar
:$HOME/freduce/bin/commons-cli.jar:$HOME/freduce/bin/commons-collections-2.1.jar
:$HOME/freduce/bin/commons-pool-1.1.jar:$HOME/freduce/bin/trove-3.0.2.jar
:$HOME/freduce/bin/jfreechart-1.0.13.jar:$HOME/freduce/bin/jcommon-1.0.16.jar
:$RJAVA_HOME/site-library/rJava/jri/JRI.jar:$HOME/freduce/bin/commons-lang-2.4.jar
:$HOME/freduce/bin/itextpdf-5.1.0.jar
export CLASSPATH

# Set the display if you're using Xvfb
DISPLAY=:1.0
export DISPLAY



*******************************************************
2.0 Logos and Graphics
*******************************************************

FeatureREDUCE generates sequence logos, FSAM logos, and plots of the positional
bias profiles using the Biojava and jFreeChart libraries.  The iText
libraries are used to save PDF versions of the figures.  Both the PNG
and PDF versions can be found in the [OUTPUT] -> Directory specified in the
INIT_FILE.  Use the optional parameter "-displayMotifs No" to turn off the
generation, display, and saving of these logos.

To generate and save the logos without displaying them, use Xvfb (X virtual
frame buffer). This tool can configure a virtual monitor that doesn't actually
display anything, but can still be used as an invisible sandbox to draw graphics.

Script start.xvfb.bash:
# run Xvfb to listen on ":1.0"
# to write graphics to this Xvfb instance set environment "DISPLAY=:1.0"
sudo Xvfb :1 -screen 0 1280x1024x24 &

*******************************************************
2.0 Execution
*******************************************************

start.xvfb.bash - a script that will start Xvfb (X virtual frame buffer) for
a server environment without a designated high resolution monitor.

.envrc - an environment script that sets needed R and FeatureREDUCE environment
variables. Execute "source .envrc" in a bash shell to set the environment.

examples.bash - a bash script that sets up an execution environment and contains
many commented examples for how to run FeatureREDUCE on the sequences from the
DREAM5 PBM Competition.

NOTE:  BEFORE YOU CAN RUN THE EXAMPLES IN THE EXAMPLES.BASH SCRIPT YOU MUST COPY
THE FOLLOWING DREAM5 COMPETITION FILES INTO THE freduce DIRECTORY (the same
directory where this README.TXT file is located):
    DREAM5_PBM_Data_Needed_For_Predictions.txt
    DREAM5_PBM_TeamName_Predictions.txt

It is strongly recommended to have a different label and output directory for
each probe intensities dataset (eg. with different microarray
preprocessing). For example, for the set of probe intensities that have been spatially
detrended use "-l dream.spatialDetrend" on the command line, and
"Directory = dream.spatialDetrend" in the init file. This will make it much
easier to keep track of your different models.

