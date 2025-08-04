setwd()

library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(caret)

####################################
### DIPHENHYDRAMINE STUDY ###
####################################

# Read data
annotations_DPH <- read.delim("merged_results_with_gnps_diphenhydramine.tsv")
metadata_DPH <- read.csv("diphenhydramine_metadata.csv", sep = ",")
metadata_druglib <- read.csv("GNPS_Drug_Library_Metadata_Drugs.csv", sep = ",")
metadata_analog_suspect <- read.csv("GNPS_Drug_Library_Metadata_Drug_Analogs_Updated.csv", sep = ",")

annotations_DPH <- annotations_DPH %>%
  dplyr::select(SpectrumID, SpectrumFile, Compound_Name, X.Scan., SharedPeaks, MQScore)

# Exclude annotations with shared peaks <=4 and cosine <0.9
annotations_DPH_filtered <- annotations_DPH %>% dplyr::filter(!(SharedPeaks <= 4 & MQScore < 0.9))

annotations_DPH_filtered$SpectrumFile <- gsub(".mzML", "", annotations_DPH_filtered$SpectrumFile)
annotations_DPH_filtered$SpectrumFile <- sub(".*/", "", annotations_DPH_filtered$SpectrumFile)

metadata_DPH$filename <- gsub(".mzXML", "", metadata_DPH$filename)

metadata_druglib <- metadata_druglib %>% dplyr::select(-smiles)
metadata_druglib$name_connected_compound <- NA
colnames(metadata_druglib)
metadata_analog_suspect <- metadata_analog_suspect %>% 
  dplyr::select(-delta_mass, -analog_mgf_scan, -analog_usi, -parent_drug_libid)
metadata_analog_suspect$name_parent_compound <- metadata_analog_suspect$name_connected_compound
colnames(metadata_analog_suspect)[colnames(metadata_analog_suspect) == "name_analog"] <- "name_compound"
colnames(metadata_analog_suspect)
colnames(metadata_analog_suspect)[colnames(metadata_analog_suspect) == "analog_libid"] <- "gnps_libid"
all_metadata <- rbind(metadata_druglib, metadata_analog_suspect)

# Combine annotations from GNPS with metadata for the GNPS drug library
DHP <- annotations_DPH_filtered %>%
  dplyr::filter(SpectrumID %in% all_metadata$gnps_libid) %>%
  left_join(all_metadata, by = c("SpectrumID" = "gnps_libid")) %>% 
  dplyr::filter(chemical_source %in% c("Medical","Drug metabolite","Drug_analog", "Drug_suspect"))

unique(metadata_DPH$Subject)

DPH_data <- metadata_DPH %>%
  dplyr::left_join(DHP, by = c("filename" = "SpectrumFile")) %>% 
  distinct(filename, name_parent_compound, .keep_all = TRUE) %>% 
  dplyr::filter(!(Subject %in% c("not applicable", "not collected")))

DPH_data_2 <- DPH_data %>%
  mutate(plasma = ifelse(str_detect(UBERONBodyPartName, "blood plasma"), 1, 0)) %>%
  mutate(skin = ifelse(str_detect(UBERONBodyPartName, "skin"), 1, 0))

DPH_data_2$timepoint_min <- as.numeric(DPH_data_2$timepoint_min)

only_DPH <- DPH_data_2 %>% 
  dplyr::filter(str_detect(name_parent_compound, "diphenhydramine"))

timepoints <- c(0, 30, 60, 90, 120, 240, 360, 480, 600, 720, 1440)

plasma_DPH <- only_DPH %>% 
  dplyr::filter(UBERONBodyPartName == "blood plasma")

collapsed_plasma_DPH <- plasma_DPH %>%
  dplyr::distinct(Subject, timepoint_min, name_parent_compound, .keep_all = TRUE)

skin_DPH <- only_DPH %>% 
  dplyr::filter(str_detect(UBERONBodyPartName, "skin"))

collapsed_skin_DPH <- skin_DPH %>%
  dplyr::distinct(Subject, timepoint_min, name_parent_compound, .keep_all = TRUE)

data_plasma_skin_DPH <- rbind(collapsed_plasma_DPH, collapsed_skin_DPH)

# Create a dataframe of all unique subjects and all timepoints
complete_timepoints <- DPH_data_2 %>%
  dplyr::distinct(Subject) %>%  
  tidyr::crossing(timepoint_min = timepoints)  

# Left join the complete_timepoints back to the original dataset
# This will fill missing timepoints with NA for the other columns
DPH_data_3 <- complete_timepoints %>%
  left_join(data_plasma_skin_DPH, by = c("Subject", "timepoint_min"))

DPH_data_3$timepoint_min <- as.numeric(DPH_data_3$timepoint_min)

DPH_data_4 <- DPH_data_3 %>%
  dplyr::mutate(plasma = replace_na(plasma, 0)) %>% 
  dplyr::mutate(skin = replace_na(skin, 0))

DPH_data_4$plasma_expected <- 1
DPH_data_4$skin_expected <- 1

detection_summary_DPH <- DPH_data_4 %>%
  group_by(timepoint_min) %>%
  summarize(
    total_samples = 10,
    plasma_detections = sum(plasma == plasma_expected),
    skin_detections = sum(skin == skin_expected),
    plasma_detection_percentage = (plasma_detections / total_samples) * 100,
    skin_detection_percentage = (skin_detections / total_samples) * 100,
  )

detection_summary_DPH <- detection_summary_DPH %>%
  mutate(timepoint_hours = timepoint_min / 60)

# Reshape the data to long format for plotting
detection_long_DPH <- detection_summary_DPH %>%
  pivot_longer(cols = starts_with("plasma_detection_percentage"):starts_with("skin_detection_percentage"),
               names_to = "Biofluid", 
               values_to = "detection_percentage") %>%
  mutate(Biofluid = case_when(
    Biofluid == "plasma_detection_percentage" ~ "Plasma",
    Biofluid == "skin_detection_percentage" ~ "Skin"
  ))

# Create the line plot
lineplot_DPH <- ggplot(detection_long_DPH, aes(x = timepoint_hours, y = detection_percentage, color = Biofluid)) +
  geom_line(size = 1.2) +     
  geom_point(size = 3) +        
  ylim(0, 100) +               
  scale_x_continuous(limits = c(0, 24), breaks = seq(0, 24, by = 3)) + 
  scale_color_manual(values = c("Plasma" = "#C7395F",   
                                "Skin" = "#EDCBD2")) +  
  labs(x = "Timepoint (hrs)",
       y = "Detection (%)",
       color = "Biofluid") +
  theme_classic() +
  theme(axis.ticks.length=unit(.25, "cm"),
        axis.line.x.bottom=element_line(size=1.2), axis.line.y.left = element_line(size=1.2),
        axis.ticks.x.bottom = element_line(size=1.2), axis.ticks.y.left = element_line(size=1.2),
        axis.text.y.left = element_text(size = 20), axis.text.x.bottom = element_text(size = 20),
        axis.title = element_text(size = 20))
print(lineplot_DPH)
#ggsave("Lineplot_validation_diphenhydramine_with_analogs.svg", plot = lineplot_DPH, width = 5.3, height = 4)

####################################
### COOPERSTOWN STUDY ###
####################################

# Read data
annotations_cooper_blood <- read.delim("merged_results_with_gnps_cooperstown_blood.tsv")
annotations_cooper_other <- read.delim("merged_results_with_gnps_cooperstown_other_biofluids.tsv")
metadata_cooper <- read.csv("cooperstown_metadata.csv", sep = ",")

annotations_cooper_blood <- annotations_cooper_blood %>%
  dplyr::select(SpectrumID, SpectrumFile, Compound_Name, X.Scan., SharedPeaks, MQScore)

annotations_cooper_other <- annotations_cooper_other %>%
  dplyr::select(SpectrumID, SpectrumFile, Compound_Name, X.Scan., SharedPeaks, MQScore)

# Exclude SharedPeaks <=4 and cosine <0.9
annotations_cooper_blood_filtered <- annotations_cooper_blood %>% dplyr::filter(!(SharedPeaks <= 4 & MQScore < 0.9))
annotations_cooper_blood_filtered$SpectrumFile <- gsub(".mzML", "", annotations_cooper_blood_filtered$SpectrumFile)
annotations_cooper_blood_filtered$SpectrumFile <- sub(".*/", "", annotations_cooper_blood_filtered$SpectrumFile)
annotations_cooper_other_filtered <- annotations_cooper_other %>% dplyr::filter(!(SharedPeaks <= 4 & MQScore < 0.9))
annotations_cooper_other_filtered$SpectrumFile <- gsub(".mzML", "", annotations_cooper_other_filtered$SpectrumFile)
annotations_cooper_other_filtered$SpectrumFile <- sub(".*/", "", annotations_cooper_other_filtered$SpectrumFile)

cooper_other_1 <- annotations_cooper_other_filtered %>%
  dplyr::filter(SpectrumID %in% all_metadata$gnps_libid) %>%
  left_join(all_metadata, by = c("SpectrumID" = "gnps_libid"))

metadata_cooper$filename <- gsub(".mzXML", "", metadata_cooper$filename)

cooper_other_2 <- metadata_cooper %>%
  dplyr::left_join(cooper_other_1, by = c("filename" = "SpectrumFile")) %>% 
  distinct(filename, name_parent_compound, .keep_all = TRUE) %>% 
  dplyr::filter(!Sample_Type == "not applicable") %>% 
  dplyr::filter(!Sample_Type == "blood_plasma")

cooper_blood_1 <- annotations_cooper_blood_filtered %>%
  dplyr::filter(SpectrumID %in% all_metadata$gnps_libid) %>%
  left_join(all_metadata, by = c("SpectrumID" = "gnps_libid")) 

cooper_blood_2 <- metadata_cooper %>%
  dplyr::right_join(cooper_blood_1, by = c("filename" = "SpectrumFile")) %>% 
  distinct(filename, name_parent_compound, .keep_all = TRUE) %>% 
  dplyr::filter(!Sample_Type == "not applicable")

cooper_data <- rbind(cooper_blood_2, cooper_other_2)

cooper_data$name_parent_compound <- ifelse(cooper_data$name_parent_compound == "omeprazole sulfide 5-carboxylic acid", "omeprazole", cooper_data$name_parent_compound)
cooper_data$name_parent_compound <- ifelse(cooper_data$name_parent_compound == "hydroxyomeprazole", "omeprazole", cooper_data$name_parent_compound)
cooper_data$name_parent_compound <- ifelse(cooper_data$name_parent_compound == "carboxyomeprazole", "omeprazole", cooper_data$name_parent_compound)
cooper_data$name_parent_compound <- ifelse(cooper_data$name_parent_compound == "esomeprazole", "omeprazole", cooper_data$name_parent_compound)
cooper_data$name_parent_compound <- ifelse(cooper_data$name_parent_compound == "omeprazole sulfide 5-carboxylic acid", "omeprazole", cooper_data$name_parent_compound)
cooper_data$name_parent_compound <- ifelse(cooper_data$name_parent_compound == "esomeprazole|omeprazole", "omeprazole", cooper_data$name_parent_compound)

#############
# PLASMA
#############
cooper_plasma <- cooper_data %>%
  dplyr::filter(Sample_Type == "blood_plasma") %>% 
  dplyr::filter(Study_Day == "1") 

unique_subjects <- cooper_plasma %>%
  distinct(Subject)

cooper_plasma_selected <- cooper_plasma %>%
  mutate(omeprazole = ifelse(str_detect(name_parent_compound, "omeprazole"), 1, 0)) %>%
  mutate(midazolam = ifelse(str_detect(name_parent_compound, "midazolam"), 1, 0)) %>%
  mutate(caffeine = ifelse(str_detect(name_parent_compound, "caffeine"), 1, 0)) 

cooper_plasma_selected <- cooper_plasma_selected %>% 
  dplyr::filter(str_detect(name_parent_compound, "omeprazole|midazolam|caffeine"))

timepoints <- c(0, 5, 30, 60, 120, 240, 300, 360, 480)

omeprazole <- cooper_plasma_selected %>% 
  filter(name_parent_compound == "omeprazole")

collapsed_omeprazole <- omeprazole %>%
  distinct(Subject, Time_Point_Mins, name_parent_compound, .keep_all = TRUE)

midazolam <- cooper_plasma_selected %>% 
  filter(name_parent_compound == "midazolam")

collapsed_midazolam <- midazolam %>%
  distinct(Subject, Time_Point_Mins, name_parent_compound, .keep_all = TRUE)

caffeine <- cooper_plasma_selected %>% 
  filter(name_parent_compound == "caffeine")

collapsed_caffeine <- caffeine %>%
  distinct(Subject, Time_Point_Mins, name_parent_compound, .keep_all = TRUE)

cooper_plasma_1 <- rbind(collapsed_omeprazole, collapsed_caffeine)
cooper_plasma_2 <- rbind(cooper_plasma_1, collapsed_midazolam)

complete_timepoints <- cooper_plasma %>%
  distinct(Subject) %>%  
  tidyr::crossing(Time_Point_Mins = timepoints)  

cooper_plasma_2$Time_Point_Mins <- as.numeric(cooper_plasma_2$Time_Point_Mins)

# Left join the complete_timepoints back to the original dataset
# This will fill missing timepoints with NA for the other columns
cooper_plasma_3 <- complete_timepoints %>%
  left_join(cooper_plasma_2, by = c("Subject", "Time_Point_Mins"))

cooper_plasma_3 <- cooper_plasma_3 %>% 
  dplyr::select(Subject, Time_Point_Mins, Study_Day, name_parent_compound, omeprazole, midazolam, caffeine) %>%
  dplyr::mutate(omeprazole = replace_na(omeprazole, 0)) %>% 
  dplyr::mutate(midazolam = replace_na(midazolam, 0)) %>%
  dplyr::mutate(caffeine = replace_na(caffeine, 0))

cooper_plasma_3$omeprazole_expected <- 1
cooper_plasma_3$midazolam_expected <- 1
cooper_plasma_3$caffeine_expected <- 1

detection_summary_cooper_plasma <- cooper_plasma_3 %>%
  group_by(Time_Point_Mins) %>%
  summarize(
    total_samples = 13,
    omeprazole_detections = sum(omeprazole == omeprazole_expected),
    midazolam_detections = sum(midazolam == midazolam_expected),
    caffeine_detections = sum(caffeine == caffeine_expected),
    omeprazole_detection_percentage = (omeprazole_detections / total_samples) * 100,
    midazolam_detection_percentage = (midazolam_detections / total_samples) * 100,
    caffeine_detection_percentage = (caffeine_detections / total_samples) * 100
  )

detection_summary_cooper_plasma <- detection_summary_cooper_plasma %>%
  mutate(timepoint_hours = Time_Point_Mins / 60)

# Reshape the data to long format for plotting
detection_long_cooper_plasma <- detection_summary_cooper_plasma %>%
  pivot_longer(cols = starts_with("omeprazole_detection_percentage"):starts_with("caffeine_detection_percentage"),
               names_to = "Drug", 
               values_to = "detection_percentage") %>%
  mutate(Drug = case_when(
    Drug == "omeprazole_detection_percentage" ~ "Omeprazole",
    Drug == "midazolam_detection_percentage" ~ "Midazolam",
    Drug == "caffeine_detection_percentage" ~ "Caffeine"
  ))

# Create the line plot for all drugs
lineplot_cooper_plasma <- ggplot(detection_long_cooper_plasma, aes(x = timepoint_hours, y = detection_percentage, color = Drug)) +
  geom_line(size = 1.2) +       
  geom_point(size = 3) +        
  ylim(0, 100) +                
  scale_x_continuous(limits = c(0, 8), breaks = seq(0, 8, by = 2)) + 
  scale_color_manual(values = c("Omeprazole" = "#3b5ba5",   
                                "Midazolam" = "#e87a5d",    
                                "Caffeine" = "#f3b941")) + 
  labs(x = "Timepoint (hrs)",
       y = "Detection (%)",
       color = "Compound name") +
  theme_classic() +
  theme(axis.ticks.length=unit(.25, "cm"),
        axis.line.x.bottom=element_line(size=1.2), axis.line.y.left = element_line(size=1.2),
        axis.ticks.x.bottom = element_line(size=1.2), axis.ticks.y.left = element_line(size=1.2),
        axis.text.y.left = element_text(size = 20), axis.text.x.bottom = element_text(size = 20),
        axis.title = element_text(size = 20))
print(lineplot_cooper_plasma)
#ggsave("Lineplot_validation_cooper_plasma_with_analog.svg", plot = lineplot_cooper_plasma, width = 5.75, height = 4)

#############
# FECES
#############
cooper_feces <- cooper_data %>%
  dplyr::filter(Sample_Type == "fecal")

unique_subjects <- cooper_feces %>%
  distinct(Subject)

cooper_feces_selected <- cooper_feces %>%
  mutate(omeprazole = ifelse(str_detect(name_parent_compound, "omeprazole"), 1, 0)) %>%
  mutate(midazolam = ifelse(str_detect(name_parent_compound, "midazolam"), 1, 0)) %>%
  mutate(caffeine = ifelse(str_detect(name_parent_compound, "caffeine"), 1, 0)) 

cooper_feces_selected <- cooper_feces_selected %>% 
  dplyr::filter(str_detect(name_parent_compound, "omeprazole|midazolam|caffeine"))

timepoints <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9)

caffeine <- cooper_feces_selected %>% 
  filter(name_parent_compound == "caffeine")

midazolam <- cooper_feces_selected %>% 
  filter(name_parent_compound == "midazolam")

omeprazole <- cooper_feces_selected %>% 
  filter(name_parent_compound == "omeprazole")

collapsed_omeprazole <- omeprazole %>%
  distinct(Subject, Study_Day, name_parent_compound, .keep_all = TRUE)

cooper_feces_1 <- rbind(collapsed_omeprazole, midazolam)

complete_timepoints <- cooper_feces %>%
  distinct(Subject) %>%  
  tidyr::crossing(Study_Day = timepoints )  

cooper_feces_1$Study_Day <- as.numeric(cooper_feces_1$Study_Day)

# Left join the complete_timepoints back to the original dataset
# This will fill missing timepoints with NA for the other columns
cooper_feces_2 <- complete_timepoints %>%
  left_join(cooper_feces_1, by = c("Subject", "Study_Day"))

cooper_feces_3 <- cooper_feces_2 %>% 
  dplyr::select(Subject, Study_Day, name_parent_compound, omeprazole, midazolam) %>%
  dplyr::mutate(omeprazole = replace_na(omeprazole, 0)) %>% 
  dplyr::mutate(midazolam = replace_na(midazolam, 0)) 

cooper_feces_3$omeprazole_expected <- 1
cooper_feces_3$midazolam_expected <- 1

detection_summary_cooper_feces <- cooper_feces_3 %>%
  group_by(Study_Day) %>%
  summarize(
    total_samples = 14,
    omeprazole_detections = sum(omeprazole == omeprazole_expected),
    midazolam_detections = sum(midazolam == midazolam_expected),
    omeprazole_detection_percentage = (omeprazole_detections / total_samples) * 100,
    midazolam_detection_percentage = (midazolam_detections / total_samples) * 100
  )

# Reshape the data to long format for plotting
detection_long_cooper_feces <- detection_summary_cooper_feces %>%
  pivot_longer(cols = starts_with("omeprazole_detection_percentage"):starts_with("midazolam_detection_percentage"),
               names_to = "Drug", 
               values_to = "detection_percentage") %>%
  mutate(Drug = case_when(
    Drug == "omeprazole_detection_percentage" ~ "Omeprazole",
    Drug == "midazolam_detection_percentage" ~ "Midazolam"
  ))

# Create the line plot for all drugs
lineplot_cooper_feces <- ggplot(detection_long_cooper_feces, aes(x = Study_Day, y = detection_percentage, color = Drug)) +
  geom_line(size = 1.2) +       
  geom_point(size = 3) +        
  ylim(0, 100) +                
  scale_x_continuous(limits = c(0, 9), breaks = seq(0, 9, by = 2)) + 
  scale_color_manual(values = c("Omeprazole" = "#3b5ba5",   
                                "Midazolam" = "#e87a5d")) + 
  labs(x = "Time (Day)",
       y = "Detection (%)",
       color = "Compound name") +
  theme_classic() +
  theme(axis.ticks.length=unit(.25, "cm"),
        axis.line.x.bottom=element_line(size=1.2), axis.line.y.left = element_line(size=1.2),
        axis.ticks.x.bottom = element_line(size=1.2), axis.ticks.y.left = element_line(size=1.2),
        axis.text.y.left = element_text(size = 20), axis.text.x.bottom = element_text(size = 20),
        axis.title = element_text(size = 20))
print(lineplot_cooper_feces)
#ggsave("Lineplot_validation_cooper_feces_with_analog.svg", plot = lineplot_cooper_feces, width = 5.75, height = 4)

