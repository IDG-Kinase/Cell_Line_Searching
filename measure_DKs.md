Dark Kinase Expression in ARCHS4 Data
================

Cell Line Counts
================

Step one is to try to get a handle on how many experiments have data concerning each of our cell lines of interest, the following code attempt to extract two sets of data from the ARCHS4 data:

-   The experiments that appear to contain data from the cell line
-   Of those experiments, which appear to be some sort of "control" experiment

This is just a quick pass through the data using exact text matching for the cell line names, if these string happen to show up in an identification string that experiment will be marked. The "control" experiment extraction is a filtered sub-set for each cell line where the string control, ctrl or DMSO (all case-insensitive) appears.

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

control_hit_sets = list()
line_hit_sets = list()
hit_counts = list()

for (this_line in cell_line) {
  line_hits = cell_line_data %>%
    filter(grepl(this_line,Sample_title,ignore.case = T) |
           grepl(this_line,Sample_characteristics_ch1,ignore.case = T) |
           grepl(this_line,Sample_description,ignore.case = T) |
           grepl(this_line,Sample_source_name_ch1,ignore.case = T))
  
  control_hits = line_hits %>%
    filter(grepl('control|ctrl|DMSO',Sample_title,ignore.case = T) |
           grepl('control|ctrl|DMSO',Sample_characteristics_ch1,ignore.case = T) |
           grepl('control|ctrl|DMSO',Sample_description,ignore.case = T) |
           grepl('control|ctrl|DMSO',Sample_source_name_ch1,ignore.case = T))
  
  if (dim(line_hits)[1] > 0) {
    line_hit_sets[[this_line]] = line_hits  
  }
  
  if (dim(control_hits)[1] > 0) {
    control_hit_sets[[this_line]] = control_hits  
    
  }
  
  hit_counts$cell_line = c(hit_counts$cell_line, this_line)
  hit_counts$line_counts = c(hit_counts$line_counts, dim(line_hits)[1])
  hit_counts$control_counts = c(hit_counts$control_counts, dim(control_hits)[1])
}
hit_counts = as.data.frame(hit_counts)

hit_counts = hit_counts %>% 
  select(cell_line,line_counts,control_counts) %>% 
  arrange(desc(control_counts))
```

Cell Line Count Table
---------------------

For reference, there are 1035 hits to the list of cell lines, with 230 matching either 'control', 'ctrl' or 'DMSO'.

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

Kinase Expression
=================

So now we have a set of experiments and we want to look at the kinase expression for each cell line. To simplify the data I'm simply going to take the average expression reported for each cell line, ignoring any variation from experiment set to experiment set.

"Control" Experiment Data Sets
------------------------------

Let's start with the experiments marked as controls, we have 12 cell lines have at least one "control" experiment.

``` r
kinase_control_data = list()
for (this_line in names(control_hit_sets)) {
  this_line_expression = expression[,which(meta_data$Sample_geo_accession %in% control_hit_sets[[this_line]]$Sample_geo_accession)]
  
  kinase_control_data$name = c(kinase_control_data$name,rep(this_line,length(gene_names)))
  kinase_control_data$gene_names = c(kinase_control_data$gene_names,gene_names)
  kinase_control_data$hgnc_id = c(kinase_control_data$hgnc_id,hgnc_id)
  
  if (class(this_line_expression) == 'matrix') {
    kinase_control_data$expression = c(kinase_control_data$expression,rowMeans(this_line_expression))
    kinase_control_data$experiment_count = c(kinase_control_data$experiment_count,
                                        rep(dim(this_line_expression)[2],length(gene_names)))
  } else {
    kinase_control_data$expression = c(kinase_control_data$expression,this_line_expression)
    kinase_control_data$experiment_count = c(kinase_control_data$experiment_count,
                                        rep(1,length(gene_names)))

  }
}

kinase_control_data = as.data.frame(kinase_control_data)

kinase_control_data = left_join(kinase_control_data, all_kinases)
```

    ## Joining, by = "hgnc_id"

    ## Warning: Column `hgnc_id` joining factor and character vector, coercing
    ## into character vector

``` r
summary_kinase_expression = kinase_control_data %>% 
  #remove all the genes with fewer than 200 pseudo-counts
  filter(expression > 200) %>%
  group_by(name) %>% 
  summarise(light_kinase_count = sum(class == "Light"),
            dark_kinase_count = sum(class == "Dark"),
            experiment_count = mean(experiment_count)) %>%
  arrange(desc(dark_kinase_count))
```

Control Expression Table
------------------------

| name    |  light\_kinase\_count|  dark\_kinase\_count|  experiment\_count|
|:--------|---------------------:|--------------------:|------------------:|
| H358    |                   273|                  106|                 22|
| H441    |                   269|                  105|                  2|
| A549    |                   265|                  102|                140|
| H1299   |                   253|                   99|                 12|
| Panc-1  |                   235|                   92|                 20|
| Capan1  |                   246|                   91|                  1|
| SW900   |                   250|                   91|                  1|
| Miapaca |                   220|                   89|                  9|
| AsPC1   |                   246|                   87|                  2|
| BxPC3   |                   224|                   83|                  3|
| H23     |                   163|                   60|                 15|
| A427    |                   126|                   40|                  3|

All Hits Experiment Data Sets
-----------------------------

Let's also go ahead and look at all the cell line hits, without regard to passing the "control" filters (¯\_(ツ)\_/¯, yolo). This one should be taken with a grain of salt as who knows how these cells were treated (if you want to thumb through the experiments though, let me know).

For reference, 20 cell lines have at

Let's start with the experiments marked as controls, we have cell lines have at least one experiment.

``` r
kinase_all_data = list()
for (this_line in names(line_hit_sets)) {
  this_line_expression = expression[,which(meta_data$Sample_geo_accession %in% line_hit_sets[[this_line]]$Sample_geo_accession)]
  
  kinase_all_data$name = c(kinase_all_data$name,rep(this_line,length(gene_names)))
  kinase_all_data$gene_names = c(kinase_all_data$gene_names,gene_names)
  kinase_all_data$hgnc_id = c(kinase_all_data$hgnc_id,hgnc_id)
  
  if (class(this_line_expression) == 'matrix') {
    kinase_all_data$expression = c(kinase_all_data$expression,rowMeans(this_line_expression))
    kinase_all_data$experiment_count = c(kinase_all_data$experiment_count,
                                        rep(dim(this_line_expression)[2],length(gene_names)))
  } else {
    kinase_all_data$expression = c(kinase_all_data$expression,this_line_expression)
    kinase_all_data$experiment_count = c(kinase_all_data$experiment_count,
                                        rep(1,length(gene_names)))

  }
}

kinase_all_data = as.data.frame(kinase_all_data)

kinase_all_data = left_join(kinase_all_data, all_kinases)
```

    ## Joining, by = "hgnc_id"

    ## Warning: Column `hgnc_id` joining factor and character vector, coercing
    ## into character vector

``` r
summary_kinase_expression_all = kinase_all_data %>% 
  #remove all the genes with fewer than 200 pseudo-counts
  filter(expression > 200) %>%
  group_by(name) %>% 
  summarise(light_kinase_count = sum(class == "Light"),
            dark_kinase_count = sum(class == "Dark"),
            experiment_count = mean(experiment_count)) %>%
  arrange(desc(dark_kinase_count))
```

Control Expression Table
------------------------

| name    |  light\_kinase\_count|  dark\_kinase\_count|  experiment\_count|
|:--------|---------------------:|--------------------:|------------------:|
| H1395   |                   304|                  123|                  1|
| H2228   |                   281|                  114|                  3|
| Capan1  |                   271|                  112|                  4|
| Capan2  |                   274|                  112|                  2|
| BxPC3   |                   278|                  110|                 18|
| H23     |                   282|                  109|                 27|
| SW1990  |                   267|                  104|                 10|
| A549    |                   263|                  102|                514|
| CFPAC-1 |                   267|                  100|                 10|
| H441    |                   254|                   97|                 14|
| HPAF-II |                   270|                   96|                  4|
| SW900   |                   249|                   96|                  2|
| Panc-1  |                   240|                   93|                114|
| AsPC1   |                   257|                   92|                 10|
| Miapaca |                   238|                   92|                 57|
| H522    |                   227|                   91|                  1|
| H1299   |                   244|                   90|                 75|
| H358    |                   245|                   88|                135|
| A427    |                   132|                   47|                 33|
| H520    |                     1|                    0|                  1|

``` r
proc.time() - processing_start_time
```

    ##    user  system elapsed 
    ##  61.096   4.834  66.255
