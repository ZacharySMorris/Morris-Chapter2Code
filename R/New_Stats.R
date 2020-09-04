###New script for angular comparison stats###


##Perform Major Axis calculation and resampling for each subset
Ecomorph_Ontogeny.MajorAxis <- resample.major.axis(Crocodylia_Ecomorphs, method = "bootstrap", iter = 9999)
Non_Gavialid_Ontogeny.MajorAxis <- resample.major.axis(Crocodylia_Ecomorphs_Non_Extreme, method = "bootstrap", iter = 9999)

Clades1.MajorAxis <- resample.major.axis(Combined_Morphospace_Clades1, method = "bootstrap", iter = 9999)
Clades1_merged.MajorAxis <- resample.major.axis(Combined_Morphospace_Clade1_merged, method = "bootstrap", iter = 9999)

Clades2.MajorAxis <- resample.major.axis(Combined_Morphospace_Clades2, method = "bootstrap", iter = 9999)
Clades3.MajorAxis <- resample.major.axis(Combined_Morphospace_Clades3, method = "bootstrap", iter = 9999)
##

##Comparisons across groups##
Ontogenies_merged_v_Clades1_merged.MajorAxis <- posthoc_MASA_comp(Non_Gavialid_Ontogeny.MajorAxis,Clades1_merged.MajorAxis)
Ontogenies_merged_v_Clades2.MajorAxis <- posthoc_MASA_comp(Non_Gavialid_Ontogeny.MajorAxis,Clades2.MajorAxis)
Ontogenies_merged_v_Clades3.MajorAxis <- posthoc_MASA_comp(Non_Gavialid_Ontogeny.MajorAxis,Clades3.MajorAxis)

Ontogenies_v_Clades1_merged.MajorAxis <- posthoc_MASA_comp(Ecomorph_Ontogeny.MajorAxis,Clades1_merged.MajorAxis)
Ontogenies_v_Clades2.MajorAxis <- posthoc_MASA_comp(Ecomorph_Ontogeny.MajorAxis,Clades2.MajorAxis)
Ontogenies_v_Clades3.MajorAxis <- posthoc_MASA_comp(Ecomorph_Ontogeny.MajorAxis,Clades3.MajorAxis)
##

