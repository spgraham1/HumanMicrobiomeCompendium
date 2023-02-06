# HumanMicrobiomeCompendium
Code used in analysis of the Human Microbiome Compendium

The "Distance to next sample" is a metric used to understand how similar any given sample is to the rest of the compendium as a whole.
Since we're interested in regional variation in the microbiome, we're using this metric to understand how well any given region is 
represented in teh compendium as a whole. We randomly select n microbiome samples, and one additional microbiome sample from the region
of interest. We then calculate the distance between this new microbiome sample and all other microbiomes in the original n samples. 
We would expect regions that are better represented in the compendium to have a lower "Distance to next sample".
