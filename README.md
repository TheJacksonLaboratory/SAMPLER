# SAMPLER: Unsupervised representations for rapid analysis of whole slide tissue images
Here we provide the code associated with the SAMPLER manuscript"SAMPLER: Unsupervised representations for rapid analysis of whole slide tissue images" [1]. We provide exmaples how to use it and share SAMPLER represeantions of frozen and FFPE breast cancer (BRCA), non-small cell lung cancer (NSCLC), and renal cell carcinoma (RCC) whole slide images (WSIs) of TCGA, as well as BRCA and NSCLC WSIs of CPTAC.

![SAMPLER overview](https://github.com/TheJacksonLaboratory/SAMPLER/blob/main/mainfig1.png)

# Downloading SAMPLER representations
Use the following links to download SAMPLER representations of:

-FFPE TCGA-BRCA WSIs using an inceptionV3 backbone pretrained on imagenet: https://zenodo.org/records/10359808

-Frozen TCGA-BRCA WSIs using an inceptionV3 backbone pretrained on imagenet: https://zenodo.org/records/10360086

-FFPE and Frozen TCGA-LUAD WSIs using an inceptionV3 backbone pretrained on imagenet: https://zenodo.org/records/10360102

-FFPE and Frozen TCGA-LUSC WSIs using an inceptionV3 backbone pretrained on imagenet: https://zenodo.org/records/10360117

-FFPE and Frozen TCGA-KIRC WSIs using an inceptionV3 backbone pretrained on imagenet: https://zenodo.org/records/10360119

-FFPE and Frozen TCGA-KIRP WSIs using an inceptionV3 backbone pretrained on imagenet: https://zenodo.org/records/10360125

-FFPE and Frozen TCGA-KICH WSIs using an inceptionV3 backbone pretrained on imagenet: https://zenodo.org/records/10360131

-FFPE TCGA-lung (LUAD and LUSC) WSIs using the single-scale tile-level features provided in [2]: https://zenodo.org/records/10360135

-CPTAC-BRCA WSIs using an inceptionV3 backbone pretrained on imagenet: https://zenodo.org/records/10360143

-CPTAC-LUAD WSIs using an inceptionV3 backbone pretrained on imagenet: https://zenodo.org/records/10360245

-CPTAC-LSCC WSIs using an inceptionV3 backbone pretrained on imagenet: https://zenodo.org/records/10360256



# Using SAMPLER
The easiest way to use SAMPLER is to copy the "SAMPLER.py" file into the main folder of the python code and load it with "import". SAMPLER only uses basic functions in os, numpy, scipy, and pickle packages, and is not sensitive to the exact version of the isntalled package. Alternatively, clone this repository and run setup.py in the SAMPLER folder. Most of SAMPLER's helper functions are for creating attention maps. If only SAMPLER representations are needed, they can be created using a single line of code:

wsi_sampler=np.ravel(np.percentile(X,pers,axis=0)).

X is the WSI feature matrix where rows are tiles and columns are features. pers are the percentiles that comprise SAMPLER represenations. We used ten deciles, i.e., pers=[5,15,25,...,95], but higher/lower resolutions can be considered. If you have the SAMPLER package installed (or copied) SAMPLER.SAMPLERrep(X) computes the SAMPLER represenation.


We have provided a Jupyter notebook example for computing SAMPLER represeantions from tile features, building SAMPLER-based predictive models, and generating SAMPLER attention maps at https://zenodo.org/records/10346576.

# References

[1] https://www.biorxiv.org/content/10.1101/2023.08.01.551468v1

[2] https://github.com/hustvl/MMIL-Transformer

