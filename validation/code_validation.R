setwd("add directory with data")

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
metadata_druglib <- read.csv("20240904_druglib_metadata_FINAL.csv", sep = ",")
metadata_analog_suspect <- read.csv("20240904_drug_analog_suspect_metadata_FINAL.csv", sep = ",")

annotations_DPH <- annotations_DPH %>%
  dplyr::select(SpectrumID, SpectrumFile, Compound_Name, X.Scan., SharedPeaks, MQScore)

# Exclude annotations with shared peaks <=4 and cosine <0.9
annotations_DPH_filtered <- annotations_DPH %>% dplyr::filter(!(SharedPeaks <= 4 & MQScore < 0.9))

annotations_DPH_filtered$SpectrumFile <- gsub(".mzML", "", annotations_DPH_filtered$SpectrumFile)
annotations_DPH_filtered$SpectrumFile <- sub(".*/", "", annotations_DPH_filtered$SpectrumFile)

metadata_DPH$filename <- gsub(".mzXML", "", metadata_DPH$filename)

# Combine metadata for GNPS drug library with analog and suspect library
metadata_druglib <- metadata_druglib %>% dplyr::select(-smiles)
metadata_druglib$name_connected_compound <- NA
metadata_analog_suspect$name_parent_compound <- metadata_analog_suspect$name_connected_compound
metadata_analog_suspect <- metadata_analog_suspect %>% dplyr::select(-gnps_libid, -analog_mgf_scan)
colnames(metadata_analog_suspect)[colnames(metadata_analog_suspect) == "analog_libid"] <- "gnps_libid"
colnames(metadata_analog_suspect)[colnames(metadata_analog_suspect) == "name_analog"] <- "name_compound"
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
  scale_color_manual(values = c("Plasma" = "#f9c58d",   
                                "Skin" = "#f492f0")) +  
  labs(x = "Timepoint (hrs)",
       y = "Detection (%)",
       color = "Biofluid") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10))
print(lineplot_DPH)
ggsave("Lineplot_validation_diphenhydramine.pdf", plot = lineplot_DPH, width = 6, height = 3.5, dpi = 900)
getwd()

# Create confusion matrix for plasma
DPH_data_plasma <- DPH_data_4 %>% 
  dplyr::filter(UBERONBodyPartName == "blood plasma") %>% 
  dplyr::select(Subject, timepoint_min, name_parent_compound,
                plasma, plasma_expected)

DPH_data_plasma_all <- complete_timepoints %>%
  left_join(DPH_data_plasma, by = c("Subject", "timepoint_min"))  %>%
  dplyr::mutate(plasma_expected = ifelse(timepoint_min == 0, 0, 1)) %>%
  dplyr::mutate(plasma = ifelse(is.na(plasma), 0, plasma))

DPH_predose_plasma <- DPH_data_plasma_all %>% 
  dplyr::filter(timepoint_min == 0)
table(DPH_predose_plasma$plasma)

DPH_postdose_plasma <- DPH_data_plasma_all %>%
  dplyr::filter(timepoint_min >= 30 & timepoint_min <= 1440) %>%
  group_by(Subject) %>%
  summarise(plasma = max(plasma)) %>%
  ungroup()
table(DPH_postdose_plasma$plasma)

conf_matrix_DPH_plasma <- matrix(c(10, 0, 0, 10), nrow = 2, byrow = TRUE)
colnames(conf_matrix_DPH_plasma) <- c("YES", "NO")
rownames(conf_matrix_DPH_plasma) <- c("YES", "NO")
conf_matrix_DPH_plasma_df <- as.data.frame(as.table(conf_matrix_DPH_plasma))

CM_DPH_plasma <- ggplot(conf_matrix_DPH_plasma_df, aes(Var2, Var1, fill = Freq)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "#f9c58d", high = "#F59C3D") +  
  geom_text(aes(label = Freq), color = "black", size = 14) +  
  theme_minimal() +
  labs(x = "Administered", y = "Detected", fill = "Count") +
  theme(axis.title = element_text(size = 16),        
        axis.text = element_text(size = 14),         
        legend.title = element_text(size = 14),      
        legend.text = element_text(size = 12)) 
ggsave("CM_validation_diphenhydramine_plasma.pdf", plot = CM_DPH_plasma, width = 6, height = 4.5, dpi = 900)
getwd()

# Generate confusion matrix for diphenhydramine detection in skin
DPH_data_skin <- DPH_data_4 %>% 
  dplyr::filter(str_detect(UBERONBodyPartName, "skin")) %>% 
  dplyr::select(Subject, timepoint_min, name_parent_compound,
                skin, skin_expected)

DPH_data_skin_all <- complete_timepoints %>%
  left_join(DPH_data_skin, by = c("Subject", "timepoint_min"))  %>%
  dplyr::mutate(skin_expected = ifelse(timepoint_min == 0, 0, 1)) %>%
  dplyr::mutate(skin = ifelse(is.na(skin), 0, skin))

DPH_predose_skin <- DPH_data_skin_all %>% 
  dplyr::filter(timepoint_min == 0)
table(DPH_predose_skin$skin)

DPH_postdose_skin <- DPH_data_skin_all %>%
  dplyr::filter(timepoint_min >= 30 & timepoint_min <= 1440) %>%
  group_by(Subject) %>%
  summarise(skin = max(skin)) %>%
  ungroup()
table(DPH_postdose_skin$skin)

conf_matrix_DPH_skin <- matrix(c(10, 0, 0, 10), nrow = 2, byrow = TRUE)
colnames(conf_matrix_DPH_skin) <- c("YES", "NO")
rownames(conf_matrix_DPH_skin) <- c("YES", "NO")
conf_matrix_DPH_skin_df <- as.data.frame(as.table(conf_matrix_DPH_skin))

CM_DPH_skin <- ggplot(conf_matrix_DPH_skin_df, aes(Var2, Var1, fill = Freq)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "#f9c58d", high = "#F59C3D") +  
  geom_text(aes(label = Freq), color = "black", size = 14) +  
  theme_minimal() +
  labs(x = "Administered", y = "Detected", fill = "Count") +
  theme(axis.title = element_text(size = 16),        
        axis.text = element_text(size = 14),         
        legend.title = element_text(size = 14),      
        legend.text = element_text(size = 12)) 
ggsave("CM_validation_diphenhydramine_skin.pdf", plot = CM_DPH_skin, width = 6, height = 4.5, dpi = 900)
getwd()

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
  scale_color_manual(values = c("Omeprazole" = "#f9c58d",   
                                "Midazolam" = "#f492f0",    
                                "Caffeine" = "#7B9E87")) + 
  labs(x = "Timepoint (hrs)",
       y = "Detection (%)",
       color = "Compound name") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10))
ggsave("Lineplot_validation_cooper_plasma.pdf", plot = lineplot_cooper_plasma, width = 6, height = 3.5, dpi = 900)
getwd()

cooper_plasma_4 <- cooper_plasma_3 %>%
  dplyr::mutate(omeprazole_expected = ifelse(Time_Point_Mins == 0, 0, 1)) %>%
  dplyr::mutate(caffeine_expected = ifelse(Time_Point_Mins == 0, 0, 1)) %>%
  dplyr::mutate(midazolam_expected = ifelse(Time_Point_Mins == 0, 0, 1)) 

# Generate confusion matrix for caffeine in plasma
caffeine_plasma_predose <- cooper_plasma_4 %>%
  dplyr::filter(Time_Point_Mins == 0)
table(caffeine_plasma_predose$caffeine)

caffeine_plasma_postdose <- cooper_plasma_4 %>%
  dplyr::filter(Time_Point_Mins >= 5 & Time_Point_Mins <= 480) %>%
  group_by(Subject) %>%
  summarise(caffeine = max(caffeine)) %>%
  ungroup()
table(caffeine_plasma_postdose$caffeine)

conf_matrix_caffeine_plasma <- matrix(c(13, 1, 0, 12), nrow = 2, byrow = TRUE)
colnames(conf_matrix_caffeine_plasma) <- c("YES", "NO")
rownames(conf_matrix_caffeine_plasma) <- c("YES", "NO")
conf_matrix_caffeine_plasma_df <- as.data.frame(as.table(conf_matrix_caffeine_plasma))

CM_caffeine_plasma <- ggplot(conf_matrix_caffeine_plasma_df, aes(Var2, Var1, fill = Freq)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "#f9c58d", high = "#F59C3D") +  
  geom_text(aes(label = Freq), color = "black", size = 14) +  
  theme_minimal() +
  labs(x = "Expected", y = "Detected", fill = "Count") +
  theme(axis.title = element_text(size = 16),        
        axis.text = element_text(size = 14),         
        legend.title = element_text(size = 14),      
        legend.text = element_text(size = 12)) 
ggsave("CM_validation_cooperstown_caffeine_plasma.pdf", plot = CM_caffeine_plasma, width = 6, height = 4.5, dpi = 900)
getwd()

# Generate confusion matrix for midazolam in plasma
midazolam_plasma_predose <- cooper_plasma_4 %>%
  dplyr::filter(Time_Point_Mins == 0)
table(midazolam_plasma_predose$midazolam)

midazolam_plasma_postdose <- cooper_plasma_4 %>%
  dplyr::filter(Time_Point_Mins >= 5 & Time_Point_Mins <= 480) %>%
  group_by(Subject) %>%
  summarise(midazolam = max(midazolam)) %>%
  ungroup()
table(midazolam_plasma_postdose$midazolam)

conf_matrix_midazolam_plasma <- matrix(c(9, 0, 4, 13), nrow = 2, byrow = TRUE)
colnames(conf_matrix_midazolam_plasma) <- c("YES", "NO")
rownames(conf_matrix_midazolam_plasma) <- c("YES", "NO")
conf_matrix_midazolam_plasma_df <- as.data.frame(as.table(conf_matrix_midazolam_plasma))

CM_midazolam_plasma <- ggplot(conf_matrix_midazolam_plasma_df, aes(Var2, Var1, fill = Freq)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "#f9c58d", high = "#F59C3D") +  
  geom_text(aes(label = Freq), color = "black", size = 14) +  
  theme_minimal() +
  labs(x = "Expected", y = "Detected", fill = "Count") +
  theme(axis.title = element_text(size = 16),        
        axis.text = element_text(size = 14),         
        legend.title = element_text(size = 14),      
        legend.text = element_text(size = 12)) 
ggsave("CM_validation_cooperstown_midazolam_plasma.pdf", plot = CM_midazolam_plasma, width = 6, height = 4.5, dpi = 900)
getwd()

# Generate confusion matrix for omeprazole in plasma
omeprazole_plasma_predose <- cooper_plasma_4 %>%
  dplyr::filter(Time_Point_Mins == 0)
table(omeprazole_plasma_predose$omeprazole)

omeprazole_plasma_postdose <- cooper_plasma_4 %>%
  dplyr::filter(Time_Point_Mins >= 5 & Time_Point_Mins <= 480) %>%
  group_by(Subject) %>%
  summarise(omeprazole = max(omeprazole)) %>%
  ungroup()
table(omeprazole_plasma_postdose$omeprazole)

conf_matrix_omeprazole_plasma <- matrix(c(13, 0, 0, 13), nrow = 2, byrow = TRUE)
colnames(conf_matrix_omeprazole_plasma) <- c("YES", "NO")
rownames(conf_matrix_omeprazole_plasma) <- c("YES", "NO")
conf_matrix_omeprazole_plasma_df <- as.data.frame(as.table(conf_matrix_omeprazole_plasma))

CM_omeprazole_plasma <- ggplot(conf_matrix_omeprazole_plasma_df, aes(Var2, Var1, fill = Freq)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "#f9c58d", high = "#F59C3D") + 
  geom_text(aes(label = Freq), color = "black", size = 14) +  
  theme_minimal() +
  labs(x = "Expected", y = "Detected", fill = "Count") +
  theme(axis.title = element_text(size = 16),        
        axis.text = element_text(size = 14),        
        legend.title = element_text(size = 14),      
        legend.text = element_text(size = 12)) 
ggsave("CM_validation_cooperstown_omeprazole_plasma.pdf", plot = CM_omeprazole_plasma, width = 6, height = 4.5, dpi = 900)
getwd()

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
  mutate(cefprozil = ifelse(str_detect(name_parent_compound, "cefprozil"), 1, 0)) %>%
  mutate(caffeine = ifelse(str_detect(name_parent_compound, "caffeine"), 1, 0)) 

cooper_feces_selected <- cooper_feces_selected %>% 
  dplyr::filter(str_detect(name_parent_compound, "omeprazole|midazolam|caffeine|cefprozil"))

timepoints <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9)

caffeine <- cooper_feces_selected %>% 
  filter(name_parent_compound == "caffeine")

cefprozil <- cooper_feces_selected %>% 
  filter(name_parent_compound == "cefprozil")

collapsed_cefprozil<- cefprozil %>%
  distinct(Subject, Study_Day, name_parent_compound, .keep_all = TRUE)

midazolam <- cooper_feces_selected %>% 
  filter(name_parent_compound == "midazolam")

omeprazole <- cooper_feces_selected %>% 
  filter(name_parent_compound == "omeprazole")

collapsed_omeprazole <- omeprazole %>%
  distinct(Subject, Study_Day, name_parent_compound, .keep_all = TRUE)

cooper_feces_1 <- rbind(collapsed_cefprozil, midazolam)
cooper_feces_2 <- rbind(cooper_feces_1, collapsed_omeprazole)

complete_timepoints <- cooper_feces %>%
  distinct(Subject) %>%  
  tidyr::crossing(Study_Day = timepoints )  

cooper_feces_2$Study_Day <- as.numeric(cooper_feces_2$Study_Day)

# Left join the complete_timepoints back to the original dataset
# This will fill missing timepoints with NA for the other columns
cooper_feces_3 <- complete_timepoints %>%
  left_join(cooper_feces_2, by = c("Subject", "Study_Day"))

cooper_feces_3 <- cooper_feces_3 %>% 
  dplyr::select(Subject, Study_Day, name_parent_compound, omeprazole, midazolam, cefprozil) %>%
  dplyr::mutate(omeprazole = replace_na(omeprazole, 0)) %>% 
  dplyr::mutate(midazolam = replace_na(midazolam, 0)) %>%
  dplyr::mutate(cefprozil = replace_na(cefprozil, 0))

cooper_feces_3$omeprazole_expected <- 1
cooper_feces_3$midazolam_expected <- 1
cooper_feces_3$cefprozil_expected <- 1

detection_summary_cooper_feces <- cooper_feces_3 %>%
  group_by(Study_Day) %>%
  summarize(
    total_samples = 14,
    omeprazole_detections = sum(omeprazole == omeprazole_expected),
    midazolam_detections = sum(midazolam == midazolam_expected),
    cefprozil_detections = sum(cefprozil == cefprozil_expected),
    omeprazole_detection_percentage = (omeprazole_detections / total_samples) * 100,
    midazolam_detection_percentage = (midazolam_detections / total_samples) * 100,
    cefprozil_detection_percentage = (cefprozil_detections / total_samples) * 100
  )

# Reshape the data to long format for plotting
detection_long_cooper_feces <- detection_summary_cooper_feces %>%
  pivot_longer(cols = starts_with("omeprazole_detection_percentage"):starts_with("cefprozil_detection_percentage"),
               names_to = "Drug", 
               values_to = "detection_percentage") %>%
  mutate(Drug = case_when(
    Drug == "omeprazole_detection_percentage" ~ "Omeprazole",
    Drug == "midazolam_detection_percentage" ~ "Midazolam",
    Drug == "cefprozil_detection_percentage" ~ "Cefprozil"
  ))

# Create the line plot for all drugs
lineplot_cooper_feces <- ggplot(detection_long_cooper_feces, aes(x = Study_Day, y = detection_percentage, color = Drug)) +
  geom_line(size = 1.2) +       
  geom_point(size = 3) +        
  ylim(0, 100) +                
  scale_x_continuous(limits = c(0, 9), breaks = seq(0, 9, by = 2)) + 
  scale_color_manual(values = c("Omeprazole" = "#f9c58d",   
                                "Midazolam" = "#f492f0",    
                                "Cefprozil" = "#7B9E87")) + 
  labs(x = "Time (Day)",
       y = "Detection (%)",
       color = "Compound name") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10))
ggsave("Lineplot_validation_cooper_feces.pdf", plot = lineplot_cooper_feces, width = 6, height = 3.5, dpi = 900)
getwd()

cooper_feces_4 <- cooper_feces_3 %>%
  dplyr::mutate(cefprozil_expected = ifelse(Study_Day == 0, 0, 1)) %>%
  dplyr::mutate(midazolam_expected = ifelse(Study_Day == 0, 0, 1)) %>%
  dplyr::mutate(omeprazole_expected = ifelse(Study_Day == 0, 0, 1)) %>%
  dplyr::mutate(cefprozil_expected = ifelse(Study_Day == 1, 0, cefprozil_expected)) 

# Generate confusion matrix for cefprozil in feces
cefprozil_feces_predose <- cooper_feces_4 %>%
  dplyr::filter(Study_Day %in% c(0, 1))
table(cefprozil_feces_predose$cefprozil)

cefprozil_feces_postdose <- cooper_feces_4 %>%
  dplyr::filter(Study_Day >= 2 & Study_Day <= 9) %>%
  group_by(Subject) %>%
  summarise(cefprozil = max(cefprozil)) %>%
  ungroup()
table(cefprozil_feces_postdose$cefprozil)

conf_matrix_cefprozil_feces <- matrix(c(10, 1, 4, 13), nrow = 2, byrow = TRUE)
colnames(conf_matrix_cefprozil_feces) <- c("YES", "NO")
rownames(conf_matrix_cefprozil_feces) <- c("YES", "NO")
conf_matrix_cefprozil_feces_df <- as.data.frame(as.table(conf_matrix_cefprozil_feces))

CM_cefprozil_feces <- ggplot(conf_matrix_cefprozil_feces_df, aes(Var2, Var1, fill = Freq)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "#f9c58d", high = "#F59C3D") +  
  geom_text(aes(label = Freq), color = "black", size = 14) +  
  theme_minimal() +
  labs(x = "Expected", y = "Detected", fill = "Count") +
  theme(axis.title = element_text(size = 16),        
        axis.text = element_text(size = 14),         
        legend.title = element_text(size = 14),      
        legend.text = element_text(size = 12)) 
ggsave("CM_validation_cooperstown_cefprozil_feces.pdf", plot = CM_cefprozil_feces, width = 6, height = 4.5, dpi = 900)
getwd()

# Generate confusion matrix for midazolam in feces
midazolam_feces_predose <- cooper_feces_4 %>%
  dplyr::filter(Study_Day %in% c(0))
table(midazolam_feces_predose$midazolam)

midazolam_feces_postdose <- cooper_feces_4 %>%
  dplyr::filter(Study_Day >= 1 & Study_Day <= 9) %>%
  group_by(Subject) %>%
  summarise(midazolam = max(midazolam)) %>%
  ungroup()
table(midazolam_feces_postdose$midazolam)

conf_matrix_midazolam_feces <- matrix(c(1, 0, 13, 14), nrow = 2, byrow = TRUE)
colnames(conf_matrix_midazolam_feces) <- c("YES", "NO")
rownames(conf_matrix_midazolam_feces) <- c("YES", "NO")
conf_matrix_midazolam_feces_df <- as.data.frame(as.table(conf_matrix_midazolam_feces))

CM_midazolam_feces <- ggplot(conf_matrix_midazolam_feces_df, aes(Var2, Var1, fill = Freq)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "#f9c58d", high = "#F59C3D") +  
  geom_text(aes(label = Freq), color = "black", size = 14) +  
  theme_minimal() +
  labs(x = "Expected", y = "Detected", fill = "Count") +
  theme(axis.title = element_text(size = 16),        
        axis.text = element_text(size = 14),         
        legend.title = element_text(size = 14),     
        legend.text = element_text(size = 12)) 
ggsave("CM_validation_cooperstown_midazolam_feces.pdf", plot = CM_midazolam_feces, width = 6, height = 4.5, dpi = 900)
getwd()

# Generate confusion matrix for omeprazole in feces
omeprazole_feces_predose <- cooper_feces_4 %>%
  dplyr::filter(Study_Day %in% c(0))
table(omeprazole_feces_predose$omeprazole)

omeprazole_feces_postdose <- cooper_feces_4 %>%
  dplyr::filter(Study_Day >= 1 & Study_Day <= 9) %>%
  group_by(Subject) %>%
  summarise(omeprazole = max(omeprazole)) %>%
  ungroup()
table(omeprazole_feces_postdose$omeprazole)

conf_matrix_omeprazole_feces <- matrix(c(4, 0, 10, 14), nrow = 2, byrow = TRUE)
colnames(conf_matrix_omeprazole_feces) <- c("YES", "NO")
rownames(conf_matrix_omeprazole_feces) <- c("YES", "NO")
conf_matrix_omeprazole_feces_df <- as.data.frame(as.table(conf_matrix_omeprazole_feces))

CM_omeprazole_feces <- ggplot(conf_matrix_omeprazole_feces_df, aes(Var2, Var1, fill = Freq)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "#f9c58d", high = "#F59C3D") +  
  geom_text(aes(label = Freq), color = "black", size = 14) +  
  theme_minimal() +
  labs(x = "Expected", y = "Detected", fill = "Count") +
  theme(axis.title = element_text(size = 16),        
        axis.text = element_text(size = 14),         
        legend.title = element_text(size = 14),      
        legend.text = element_text(size = 12)) 
ggsave("CM_validation_cooperstown_omeprazole_feces.pdf", plot = CM_omeprazole_feces, width = 6, height = 4.5, dpi = 900)
getwd()
