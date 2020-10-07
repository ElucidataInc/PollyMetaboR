### Introduction

The PollyMetaboR can be used to process metabolomics data generated using untargeted approach. This package works on the peak detailed format of El-MAVEN output.  

### Steps to use this package


```R
library(PollyMetaboR)
```

#### Load Data


```R
data(demo_peak_detailed_elmaven, package = "PollyMetaboR")
data(KEGG_mzMass, package = "PollyMetaboR")
```

#### Make XCMS Object


```R
xcms_obj <- create_xcms_object_from_elmaven(maven_data)
```

#### Perform Annotation using CAMERA


```R
camera_output <- perform_annotation_by_camera(xcms_obj, polarity = "positive", ppm = 10, mzabs = 0.1)
```


#### Restructure CAMERA output


```R
restructure_camera <- restructure_camera_annotations(camera_output, polarity = "positive")
```

#### Get feature representative for each feature group


```R
representative_df <- get_feature_group_representative(restructure_camera$combined, polarity = "positive")

```


#### Perform Metabolite Identification


```R
identified_df <- perform_metabolite_identification(mz_data = representative_df, comp_data = KEGG_mzMass, mz_colname = 'basemass', mz_tolerence_unit= 'ppm', mz_tolerence = 20, numcores = 4)

```

#### Merge identified data with restrutured camera output


```R
metab_ident <- merge_identified_with_restructured_camera(identified_df, restructure_camera$combined)

```

#### Summary of Annotation


```R
p <- plot_hist_elements_frequency(restructure_camera$combined$pcgroup, frequency_type = "by_occurrence", 
                                  plot_title = "Number of features vs counts of pcgroup",
                                  xaxis_title = "Number of features",
                                  yaxis_title = "Counts of pcgroup")

```

#### Make group summary format from identification step output format


```R
metscape <- make_group_summary_from_metab_ident_format(identified_df)
```