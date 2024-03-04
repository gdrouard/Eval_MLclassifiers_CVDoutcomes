# Repository of: "Evaluation of machine learning strategies for prediction of cardiovascular risk factors from multi-omic data"

**Abbreviations:**
ML: machine learning. SSAE: semi-supervised autoencoder. USAE: Unsupervised autoencoder. X.train: training set. Y.train: considered outcome in the training set. X.test: test sample. Y.test: observed outcome in the test sample.

Any question can be addressed to: *gabin.drouard@helsinki.fi* 

**Script description:**

*Compare USAE vs. SSAE 1.R:* First script used to evaluate in which settings SSAE outperformed USAE. The plots were used in the supplementary material.

*Compare USAE vs. SSAE 2.R:* Second script to evaluate how often SSAE outperformed USAE.

*Data-adapted VCDN.R:* A script that uses single-omic predictions to generate multi-omic predictions. The script is adapted to the datasets used and an introduction to VCDN can be found at doi: 10.1038/s41467-021-23774-w.

*Dimension reduction using SSAE and USAE and variable importance investigation.R*: Script showing how autoencoders were constructed and used to reduce the dimensionality of transcriptomic and metabolomic data.

*Plot Figure 4.R*: Script showing how Figure 4, comparing F1 scores of ML classifiers, was generated.

*Training Classifiers.R*: A generic script showing how ML classifiers were applied to the data to predict 1-sd classes.

*Transfer learning pre-trained SSAE - transcriptomics.R*: Script showing how pre-trained SSAE were externalised to another cohort to predict hypertensive participants.
