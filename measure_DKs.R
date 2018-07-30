library(rhdf5)
library(tidyverse)


cell_line = c("AsPC1","BxPC3","Capan1","Capan2","CFPAC-1","HPAF-II","Hs766T",
              "Miapaca","Panc-1","SW1990","A427","A549","H1299","H1395",
              "H1703","NCI-H1703","H2126","H2170","H2228","H23","H358","H441",
              "H520","H522","H727","SK-MES","SW900")

# meta_data = h5read('human_matrix.h5','/meta')

cell_line_data = as.data.frame(
  meta_data[c("Sample_characteristics_ch1",
              "Sample_description",
              "Sample_source_name_ch1",
              "Sample_title")])

hit_sets = list()

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
  
  hit_sets$title = c(hit_sets$title,title_hits)
  hit_sets$control = c(hit_sets$title,control_hits)
  
  if (dim(control_hits)[1] == 0) {
    print(this_line)
    print(dim(title_hits))
    print(dim(control_hits))
    print("")
  }
  
}