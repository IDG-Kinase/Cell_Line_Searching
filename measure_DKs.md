Dark Kinase Expression in ARCHS4 Data
================

Cell Line Counts
================

``` r
cell_line = c("AsPC1","BxPC3","Capan1","Capan2","CFPAC-1","HPAF-II","Hs766T",
              "Miapaca","Panc-1","SW1990","A427","A549","H1299","H1395",
              "H1703","H1703","H2126","H2170","H2228","H23","H358","H441",
              "H520","H522","H727","SK-MES","SW900")

cell_line_data = as.data.frame(
  meta_data[c("Sample_characteristics_ch1",
              "Sample_description",
              "Sample_source_name_ch1",
              "Sample_title",
              "Sample_geo_accession")])

hit_sets = list()
hit_counts = list()

for (this_line in cell_line) {
  title_hits = cell_line_data %>%
    filter(grepl(this_line,Sample_title,ignore.case = T) |
           grepl(this_line,Sample_characteristics_ch1,ignore.case = T) |
           grepl(this_line,Sample_description,ignore.case = T) |
           grepl(this_line,Sample_source_name_ch1,ignore.case = T))
  
  control_hits = title_hits %>%
    filter(grepl('control|ctrl|DMSO',Sample_title,ignore.case = T) |
           grepl('control|ctrl|DMSO',Sample_characteristics_ch1,ignore.case = T) |
           grepl('control|ctrl|DMSO',Sample_description,ignore.case = T) |
           grepl('control|ctrl|DMSO',Sample_source_name_ch1,ignore.case = T))
  
  hit_sets$this_line = c(hit_sets$this_line,this_line)
  hit_sets$title = c(hit_sets$title,title_hits)
  hit_sets$control = c(hit_sets$control,control_hits)
  
  hit_counts$cell_line = c(hit_counts$cell_line, this_line)
  hit_counts$line_counts = c(hit_counts$line_counts, dim(title_hits)[1])
  hit_counts$control_counts = c(hit_counts$control_counts, dim(control_hits)[1])
}
hit_counts = as.data.frame(hit_counts)

hit_counts = hit_counts %>% select(cell_line,line_counts,control_counts) %>% arrange(desc(control_counts))
```

Cell Line Count Table
=====================

``` r
knitr::kable(hit_counts)
```

| cell\_line |  line\_counts|  control\_counts|
|:-----------|-------------:|----------------:|
| A549       |           514|              140|
| H358       |           135|               22|
| Panc-1     |           114|               20|
| H23        |            27|               15|
| H1299      |            75|               12|
| Miapaca    |            57|                9|
| BxPC3      |            18|                3|
| A427       |            33|                3|
| AsPC1      |            10|                2|
| H441       |            14|                2|
| Capan1     |             4|                1|
| SW900      |             2|                1|
| Capan2     |             2|                0|
| CFPAC-1    |            10|                0|
| HPAF-II    |             4|                0|
| Hs766T     |             0|                0|
| SW1990     |            10|                0|
| H1395      |             1|                0|
| H1703      |             0|                0|
| H1703      |             0|                0|
| H2126      |             0|                0|
| H2170      |             0|                0|
| H2228      |             3|                0|
| H520       |             1|                0|
| H522       |             1|                0|
| H727       |             0|                0|
| SK-MES     |             0|                0|
