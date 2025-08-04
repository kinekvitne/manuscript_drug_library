setwd()

library(dplyr)
library(readxl)
library(tidyverse)
library(stringr)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(rworldmap)
library(sf)

# Read data
feature_table <- read.csv("20240507_AGP_NoGapFill_BlkRm_LODcheck.csv", sep = ",") 
annotations <- read.delim("merged_results_with_gnps.tsv", header = TRUE, stringsAsFactors = FALSE) 
ReDu_AGP <- read.csv("ReDu_metadata_AGP.csv") #MSV000080673 
data_10317_2136 <- read.csv("qiita_metadata_IDs_10317_2136.csv") 

feature_table_t <- feature_table %>% column_to_rownames("X") %>% t() %>% as.data.frame() %>% rownames_to_column("SampleID")
feature_table_t_pivot <- feature_table_t %>% pivot_longer(!SampleID, names_to = "Drug", values_to = "PeakArea") %>% dplyr::filter(PeakArea =="TRUE")

# Exclude annotations with shared peaks <=4 and cosine <0.9
annotations <- annotations %>% dplyr::filter(!(SharedPeaks <= 4 & MQScore < 0.9))

# Manual inspection (shared peaks <= 4 and cosine >= 0.9)
# Features to be excluded
annotations_curated <- annotations %>%
  dplyr::filter(!(X.Scan. %in% c("3799")))

AGP_annotations <- annotations_curated %>%
  dplyr::select(SpectrumID, Compound_Name, X.Scan.)
AGP_annotations$X.Scan. <- as.character(AGP_annotations$X.Scan.)

AGP_info <- feature_table_t_pivot %>% left_join(AGP_annotations, by = c("Drug" = "X.Scan."))
AGP_info$SampleID <- gsub(".mzML.Peak.area", "", AGP_info$SampleID)
AGP_info$SampleID <- substring(AGP_info$SampleID, 2)

ReDu_AGP$filename <- str_extract(ReDu_AGP$filename, "[^/]*$")
ReDu_AGP$filename <- gsub(".mzML", "", ReDu_AGP$filename)
ReDu_AGP$filename <- gsub(".mzXML", "", ReDu_AGP$filename)

ReDu_AGP_annotations <- AGP_info %>% left_join(ReDu_AGP, by = c("SampleID" = "filename"))

AGP_combined <- ReDu_AGP_annotations %>% 
  dplyr::mutate(anonymized_name = gsub("_.*", "", SampleID)) %>% 
  left_join(data_10317_2136)

# Look at samples collected from the ICU in AGP
AGP_2136 <- AGP_combined %>%
  dplyr::filter(Cohort == "ICU")
unique_IDs_2136 <- AGP_2136 %>% # n = 93
  dplyr::distinct(host_subject_id, .keep_all = TRUE)

# Exclude ICU samples
AGP_10317 <- AGP_combined %>%
  dplyr::filter(Cohort != "ICU") %>%
  filter(!is.na(ATTRIBUTE_DatasetAccession))

# Look if any IDs have multiple samples
ids_with_multiple_timestamps <- AGP_10317 %>%
  group_by(host_subject_id) %>%
  filter(n_distinct(collection_timestamp) > 1) %>%
  distinct(host_subject_id) # 47 individuals had multiple samples

timestamp_counts <- AGP_10317 %>%
  group_by(host_subject_id) %>%
  summarize(n_timestamps = n_distinct(collection_timestamp)) %>%
  filter(n_timestamps > 1)

# Keep only the first sample for individuals with multiple samples
filtered_data <- AGP_10317 %>%
  group_by(host_subject_id) %>%
  slice_min(order_by = collection_timestamp, with_ties = FALSE) %>%
  dplyr::select(host_subject_id, collection_timestamp)

merged_data <- AGP_10317 %>%
  inner_join(filtered_data, by = c("host_subject_id", "collection_timestamp"))

# Exclude countries with few samples
country_counts <- merged_data %>%
  group_by(Country) %>%
  summarize(unique_host_subjects = n_distinct(host_subject_id)) %>%
  arrange(desc(unique_host_subjects))

selected_regions <- merged_data %>%
  mutate(Region = case_when(
    Country == "United States of America" ~ "United States",
    Country %in% c("United Kingdom", "France", "Germany", "Norway", "Netherlands", 
                   "Ireland", "Italy", "Sweden", "Latvia", "Spain", "Austria", "Czech Republic",
                   "Finland", "Switzerland", "Slovakia") ~ "Europe",
    Country == "Australia" ~ "Australia",
    TRUE ~ NA_character_  
  )) %>%
  dplyr::filter(Region %in% c("Australia", "Europe", "United States"))

unique_IDs_AGP_10317 <- selected_regions %>% # n = 1,993
  dplyr::distinct(host_subject_id, .keep_all = TRUE)

# Save the final AGP data that will be used for analysis to a csv file to avoid rerunning the entire code
#write.csv(selected_regions, "AGP_10317_final.csv")

# Read final data AGP
AGP_10317 <- read.csv("AGP_10317_final.csv", sep = ",") %>%
  dplyr::select(-1)

# Look at demographic data such as age and sex
AGP_age_sex <- AGP_10317 %>% # Have info on age and sex for 1845 individuals
  dplyr::distinct(host_subject_id, .keep_all = TRUE) %>%
  dplyr::mutate(AgeInYears = round(as.numeric(AgeInYears))) %>%
  filter(!is.na(AgeInYears),
         !is.na(BiologicalSex),
         !AgeInYears %in% c("not collected", "not applicable", "not specified", ""),
         !BiologicalSex %in% c("not collected", "not applicable", "not specified", ""))

class_counts <- table(AGP_age_sex$BiologicalSex)
print(class_counts) 
mean(AGP_age_sex$AgeInYears)
sd(AGP_age_sex$AgeInYears) 
min(AGP_age_sex$AgeInYears) 
max(AGP_age_sex$AgeInYears) 

# Read data GNPS drug library
metadata_druglib <- read.csv("GNPS_Drug_Library_Metadata_Drugs.csv", sep = ",")
metadata_druglib <- metadata_druglib %>%
  dplyr::select(-smiles)
metadata_druglib$name_connected_compound <- NA
colnames(metadata_druglib)

# Read data GNPS analog/suspect library
metadata_analog_suspect <- read.csv("GNPS_Drug_Library_Metadata_Drug_Analogs_Updated.csv", sep = ",")
metadata_analog_suspect <- metadata_analog_suspect %>% 
  dplyr::select(-delta_mass, -analog_mgf_scan, -analog_usi, -parent_drug_libid)
metadata_analog_suspect$name_parent_compound <- metadata_analog_suspect$name_connected_compound
colnames(metadata_analog_suspect)[colnames(metadata_analog_suspect) == "name_analog"] <- "name_compound"
colnames(metadata_analog_suspect)
colnames(metadata_analog_suspect)[colnames(metadata_analog_suspect) == "analog_libid"] <- "gnps_libid"

# Combine metadata for the GNPS drug library and analog/suspect library
all_metadata <- rbind(metadata_druglib, metadata_analog_suspect)

# Combine data with the GNPS drug library metadata
# Remove duplicates
AGP_druglib <- AGP_10317 %>%
  dplyr::filter(SpectrumID %in% all_metadata$gnps_libid) %>%
  left_join(all_metadata, by = c("SpectrumID" = "gnps_libid")) %>% 
  distinct(SampleID, name_compound, .keep_all = TRUE) %>% 
  dplyr::filter(chemical_source %in% c("Medical","Drug metabolite","Drug_analog"))

# Generate a frequency table to summarize detected drugs
freq_table <- table(AGP_druglib$name_parent_compound) 
freq_df <- as.data.frame(freq_table) %>% arrange(desc(Freq))

# Standardize specific compound names for consistency before downstream analysis
AGP_druglib$name_parent_compound <- ifelse(AGP_druglib$name_parent_compound == "n-acetyl mesalazine", "mesalazine", AGP_druglib$name_parent_compound)
AGP_druglib$name_parent_compound <- ifelse(AGP_druglib$name_parent_compound == "omeprazole sulfide 5-carboxylic acid", "omeprazole", AGP_druglib$name_parent_compound)
AGP_druglib$name_parent_compound <- ifelse(AGP_druglib$name_parent_compound == "atenolol|metoprolol", "metoprolol", AGP_druglib$name_parent_compound)
AGP_druglib$name_parent_compound <- ifelse(AGP_druglib$name_parent_compound == "clarithromycin|erythromycin", "erythromycin", AGP_druglib$name_parent_compound)
AGP_druglib$name_parent_compound <- ifelse(AGP_druglib$name_parent_compound == "quetiapine sulfoxide ", "quetiapine", AGP_druglib$name_parent_compound)
AGP_druglib$name_parent_compound <- ifelse(AGP_druglib$name_parent_compound == "clindamycin n-oxide", "clindamycin", AGP_druglib$name_parent_compound)
AGP_druglib$name_parent_compound <- ifelse(AGP_druglib$name_parent_compound == "clindamycin sulfoxide", "clindamycin", AGP_druglib$name_parent_compound)
AGP_druglib$name_parent_compound <- ifelse(AGP_druglib$name_parent_compound == "albendazole ", "albendazole", AGP_druglib$name_parent_compound)
AGP_druglib$name_parent_compound <- ifelse(AGP_druglib$name_parent_compound == "losartan|losartan", "losartan", AGP_druglib$name_parent_compound)
AGP_druglib$name_parent_compound <- ifelse(AGP_druglib$name_parent_compound == "losartan ", "losartan", AGP_druglib$name_parent_compound)


# Remove duplicates again since for some individuals we detect both the parent compound and drug metabolites
AGP_druglib <- AGP_druglib %>%
  dplyr::filter(chemical_source %in% c("Medical","Drug metabolite","Drug_analog")) %>% 
  distinct(SampleID, name_parent_compound, .keep_all = TRUE)

# Save the final AGP data combined with the GNPS drug library metadata to a csv file to avoid rerunning the entire code
#write.csv(AGP_druglib, "AGP_10317_druglib_annotations.csv")

AGP_druglib <- read.csv("AGP_10317_druglib_annotations.csv", sep = ",") %>%
  dplyr::select(-1) 

# Count number of unique parent compounds detected
count_detected_drugs <- AGP_druglib %>% 
  dplyr::distinct(name_parent_compound, .keep_all = TRUE) # n = 89

# Count number of unique IDs where drug was detected
unique_IDs_drug <- AGP_druglib %>%
  dplyr::distinct(host_subject_id, .keep_all = TRUE)

###### FIGURE 1g - TOP DETECTED PHARMACOLOGIC CLASSES ######
# Remove duplicates in case some individuals use more than one drug within the same pharmacologic class
topdrugs <- AGP_druglib %>% 
  dplyr::distinct(host_subject_id, pharmacologic_class, .keep_all = TRUE) 
topdrugs <- topdrugs %>% dplyr::filter(pharmacologic_class != "no_match")

# Count occurrences of each drug
drug_counts <- topdrugs %>% 
  group_by(pharmacologic_class) %>% 
  summarise(drug_count = n())

# Calculate percentage
total_individuals <- n_distinct(AGP_10317$host_subject_id)
drug_percentages <- drug_counts %>%
  mutate(drug_percentage = (drug_count / total_individuals) * 100) %>%
  dplyr::arrange(desc(drug_percentage)) %>%
  head(22) %>% 
  dplyr::filter(!pharmacologic_class %in% c("antirheumatic|antimalarial", "antileishmanial drug"))


# Figure 1g
top20_pharmacologic_class <- ggplot(drug_percentages, aes(x = reorder(pharmacologic_class, drug_percentage), y = drug_percentage)) +
  geom_bar(stat = "identity", fill = "#990066", width = 0.7) +  
  geom_text(aes(label = paste0(round(drug_percentage, 1), "%")), vjust = 0.4, hjust = -0.1, size = 3) +  
  labs(title = "Top 20 Most Detected FDA Drug Classes - Feces",
       x = "",
       y = "Percentage (%)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_flip()
top20_pharmacologic_class
#ggsave("AGP_top20_pharmacologic_class.svg", plot = top20_pharmacologic_class, width = 6, height = 4, dpi = 900)
#getwd()

###### FIGURE 1h - DETECTED DRUGS ACROSS AGE AND SEX ######
# Create age bins to categorize individuals into age groups. Exclude individuals 90-100 yrs
bins <- seq(0, 100, by = 10)
AGP_AgeGroup <- AGP_age_sex %>%
  mutate(AgeGroup = cut(AgeInYears, breaks = bins, right = FALSE)) %>%
  dplyr::filter(!AgeGroup %in% c("[90,100)"))

# Calculate total sample count per age group and sex
total_sample_count <- AGP_AgeGroup %>%
  group_by(AgeGroup, BiologicalSex) %>%
  summarise(TotalSampleCount = n(), .groups = 'drop') 

# Filter to keep only individuals with available age and sex data
AGP_final_drug_data <- AGP_druglib %>%
  dplyr::mutate(AgeInYears = round(as.numeric(AgeInYears))) %>%
  filter(!is.na(AgeInYears),
         !is.na(BiologicalSex),
         !AgeInYears %in% c("not collected", "not applicable", "not specified", ""),
         !BiologicalSex %in% c("not collected", "not applicable", "not specified", ""))

# Define age bins to group into age groups
bins <- seq(0, 100, by = 10)
AGP_final_drug_data <- AGP_final_drug_data %>%
  mutate(AgeGroup = cut(AgeInYears, breaks = bins, right = FALSE)) %>%
  dplyr::filter(!AgeGroup %in% c("[90,100)"))

# CARDIOVASCULAR DRUGS
# Filter for cardiovascular drugs and remove duplicates
CV_HT_drugs <- AGP_final_drug_data %>%
  dplyr::filter(str_detect(therapeutic_area, "cardiolog")) %>%
  dplyr::filter(!str_detect(therapeutic_indication, "erectile dysfunction")) %>%
  dplyr::filter(!str_detect(name_compound, "phenylephr")) %>%
  distinct(host_subject_id, .keep_all = TRUE) 

cardiology_counts <- CV_HT_drugs %>%
  group_by(AgeGroup, BiologicalSex) %>%
  summarise(DrugCount = n(), .groups = 'drop')

# Merge with total sample counts and calculate normalized values
normalized_data <- cardiology_counts %>%
  left_join(total_sample_count, by = c("AgeGroup", "BiologicalSex")) %>%
  mutate(normalized = DrugCount / TotalSampleCount)

# Create a data frame with all combinations of age groups and biological sexes
age_groups <- c("[0,10)", "[10,20)", "[20,30)", "[30,40)", 
                "[40,50)", "[50,60)", "[60,70)", "[70,80)", "[80,90)")
biological_sexes <- c("female", "male")
all_combinations <- expand.grid(AgeGroup = age_groups, BiologicalSex = biological_sexes)

# Combine with original data
merged_data <- merge(all_combinations, normalized_data, by = c("AgeGroup", "BiologicalSex"), all = TRUE) %>%
  mutate(normalized = ifelse(is.na(normalized), 0, normalized))

colors <- c("female" = "#C9D3EC", "male" = "#FBB465")

cvd_normalized <- ggplot(merged_data, aes(x = AgeGroup, y = normalized, fill = BiologicalSex)) +
  geom_col(position = "dodge", alpha=0.5) +
  geom_smooth(aes(group = BiologicalSex, color = BiologicalSex), method = "loess", span = 0.5, se = F, linetype = "solid", size = 2) +
  labs(title = "Cardiovascular drugs",
       x = "Age (yrs)",
       y = "Normalized Count", 
       fill = "") +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),  
    axis.line = element_line(color = "black"),  
    axis.ticks = element_line(color = "black"),  
    axis.ticks.length = unit(0.2, "cm"),  
    axis.text = element_text(color = "black"),  
    axis.title = element_text(color = "black"),  
    plot.title = element_text(hjust = 0.5),  
    axis.text.x = element_text(angle = 45, hjust = 1))
cvd_normalized
#ggsave("AGP_cardiovascular_drugs.svg", plot = cvd_normalized, width = 6, height = 4, dpi = 900)
#getwd()

# DRUGS FOR ERECTILE DYSFUNCTION
# Filter for drugs used for erectile dysfunction and remove duplicates
male_drugs <- AGP_final_drug_data %>% dplyr::filter(str_detect(therapeutic_indication, "erectile")) %>%
  distinct(host_subject_id, .keep_all = TRUE)

male_counts <- male_drugs %>%
  group_by(AgeGroup, BiologicalSex) %>%
  summarise(DrugCount = n(), .groups = 'drop') 

normalized_data_2 <- male_counts %>%
  left_join(total_sample_count, by = c("AgeGroup", "BiologicalSex")) %>%
  mutate(normalized = DrugCount / TotalSampleCount)

# Combine with original data
merged_data_2 <- merge(all_combinations, normalized_data_2, by = c("AgeGroup", "BiologicalSex"), all = TRUE) %>%
  mutate(normalized = ifelse(is.na(normalized), 0, normalized))

male_normalized <- ggplot(merged_data_2, aes(x = AgeGroup, y = normalized, fill = BiologicalSex)) +
  geom_col(position = "dodge", alpha=0.5) +
  geom_smooth(aes(group = BiologicalSex, color = BiologicalSex), method = "loess", span = 0.5, se = F, linetype = "solid", size = 2, , span=0.5) +
  labs(title = "Drugs for erectile dysfunction",
       x = "Age (yrs)",
       y = "Normalized Count", 
       fill = "") +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),  
    axis.line = element_line(color = "black"),  
    axis.ticks = element_line(color = "black"),  
    axis.ticks.length = unit(0.2, "cm"),  
    axis.text = element_text(color = "black"),  
    axis.title = element_text(color = "black"),  
    plot.title = element_text(hjust = 0.5),  
    axis.text.x = element_text(angle = 45, hjust = 1))
male_normalized
#ggsave("AGP_erectile_dysfunction.svg", plot = male_normalized, width = 6, height = 4, dpi = 900)
#getwd()

# ANTIDEPRESSANTS
depression_drugs <- AGP_final_drug_data %>% dplyr::filter(str_detect(therapeutic_indication, "depress")) %>%
  distinct(host_subject_id, .keep_all = TRUE)

depression_counts <- depression_drugs %>%
  group_by(AgeGroup, BiologicalSex) %>%
  summarise(DrugCount = n(), .groups = 'drop') 

normalized_data_3 <- depression_counts %>%
  left_join(total_sample_count, by = c("AgeGroup", "BiologicalSex")) %>%
  mutate(normalized = DrugCount / TotalSampleCount)

merged_data_3 <- merge(all_combinations, normalized_data_3, by = c("AgeGroup", "BiologicalSex"), all = TRUE) %>%
  mutate(normalized = ifelse(is.na(normalized), 0, normalized))

depression_normalized <- ggplot(merged_data_3, aes(x = AgeGroup, y = normalized, fill = BiologicalSex)) +
  geom_col(position = "dodge", alpha=0.5) +
  geom_smooth(aes(group = BiologicalSex, color = BiologicalSex), method = "loess", span = 0.5, se = F, linetype = "solid", size = 2) +
  labs(title = "Antidepressants",
       x = "Age (yrs)",
       y = "Normalized Count", 
       fill = "") +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),  
    axis.line = element_line(color = "black"),  
    axis.ticks = element_line(color = "black"),  
    axis.ticks.length = unit(0.2, "cm"),  
    axis.text = element_text(color = "black"),  
    axis.title = element_text(color = "black"),  
    plot.title = element_text(hjust = 0.5),  
    axis.text.x = element_text(angle = 45, hjust = 1))
depression_normalized
#ggsave("AGP_antidepressants.svg", plot = depression_normalized, width = 6, height = 4, dpi = 900)
#getwd()

# PAINKILLERS
pain_drugs <- AGP_final_drug_data %>% dplyr::filter(str_detect(mechanism_of_action, "cyclooxygenase")) %>%
  dplyr::filter(!name_parent_compound == "mesalazine") %>% 
  distinct(host_subject_id, .keep_all = TRUE)

pain_counts <- pain_drugs %>%
  group_by(AgeGroup, BiologicalSex) %>%
  summarise(DrugCount = n(), .groups = 'drop') 

normalized_data_4 <- pain_counts %>%
  left_join(total_sample_count, by = c("AgeGroup", "BiologicalSex")) %>%
  mutate(normalized = DrugCount / TotalSampleCount)

merged_data_4 <- merge(all_combinations, normalized_data_4, by = c("AgeGroup", "BiologicalSex"), all = TRUE) %>%
  mutate(normalized = ifelse(is.na(normalized), 0, normalized))

pain_normalized <- ggplot(merged_data_4, aes(x = AgeGroup, y = normalized, fill = BiologicalSex)) +
  geom_col(position = "dodge", alpha=0.5) +
  geom_smooth(aes(group = BiologicalSex, color = BiologicalSex), method = "loess", span = 0.5, se = F, linetype = "solid", size = 2) +
  labs(title = "Analgesics",
       x = "Age (yrs)",
       y = "Normalized Count", 
       fill = "") +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),  
    axis.line = element_line(color = "black"),  
    axis.ticks = element_line(color = "black"),  
    axis.ticks.length = unit(0.2, "cm"),  
    axis.text = element_text(color = "black"),  
    axis.title = element_text(color = "black"),  
    plot.title = element_text(hjust = 0.5), 
    axis.text.x = element_text(angle = 45, hjust = 1))
pain_normalized
#ggsave("AGP_painkillers.svg", plot = pain_normalized, width = 6, height = 4, dpi = 900)
#getwd()

# Perform Chi-square test to investigate differences in painkiller detection by sex
AGP_sex <- AGP_10317 %>%
  dplyr::distinct(host_subject_id, .keep_all = TRUE) %>%
  dplyr::filter(!is.na(BiologicalSex),
                !BiologicalSex %in% c("not collected", "not applicable", "not specified", ""))

total_counts2 <- AGP_sex %>%
  group_by(BiologicalSex) %>%
  summarise(Count = n(), .groups = 'drop') 

pain_counts2 <- pain_drugs %>%
  group_by(BiologicalSex) %>%
  summarise(Count = n(), .groups = 'drop') 

# Create the contingency table and perform Chi-square test
contingency_table <- matrix(c(64, 977, 22, 895), nrow = 2, byrow = TRUE)
rownames(contingency_table) <- c("Female", "Male")
colnames(contingency_table) <- c("Used Painkillers", "Did Not Use Painkillers")
chisq.test(contingency_table)

# ALLERGY
allergy_drugs <- AGP_final_drug_data %>% dplyr::filter(str_detect(therapeutic_indication, "allerg")) %>%
  distinct(host_subject_id, .keep_all = TRUE)

allergy_counts <- allergy_drugs %>%
  group_by(AgeGroup, BiologicalSex) %>%
  summarise(DrugCount = n(), .groups = 'drop') 

normalized_data_5 <- allergy_counts %>%
  left_join(total_sample_count, by = c("AgeGroup", "BiologicalSex")) %>%
  mutate(normalized = DrugCount / TotalSampleCount)

merged_data_5 <- merge(all_combinations, normalized_data_5, by = c("AgeGroup", "BiologicalSex"), all = TRUE) %>%
  mutate(normalized = ifelse(is.na(normalized), 0, normalized))

allergy_normalized <- ggplot(merged_data_5, aes(x = AgeGroup, y = normalized, fill = BiologicalSex)) +
  geom_col(position = "dodge", alpha=0.5) +
  geom_smooth(aes(group = BiologicalSex, color = BiologicalSex), method = "loess", span = 0.5, se = F, linetype = "solid", size = 2) +
  labs(title = "Antihistamines",
       x = "Age (yrs)",
       y = "Normalized Count", 
       fill = "") +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),  
    axis.line = element_line(color = "black"),  
    axis.ticks = element_line(color = "black"), 
    axis.ticks.length = unit(0.2, "cm"),  
    axis.text = element_text(color = "black"),  
    axis.title = element_text(color = "black"),  
    plot.title = element_text(hjust = 0.5),  
    axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(limits = c(0.00, 0.14), breaks = c(0, 0.03, 0.06, 0.09, 0.12))
allergy_normalized
#ggsave("AGP_antihistamines.svg", plot = allergy_normalized, width = 6, height = 4, dpi = 900)
#getwd()

# ANTIBIOTICS
# Filter to keep only antibiotics
abx_drugs <- AGP_final_drug_data %>% dplyr::filter(str_detect(therapeutic_area, "infectious")) %>%
  filter(!(name_parent_compound %in% c("acyclovir", "amphotericin b", "elvitegravir","fluconazole", "hydroxychloroquine",
                                       "miltefosine","raltegravir",  
                                       "nystatin", "praziquantel", "pyrantel", "rimantadine", "quinine"))) %>%
  distinct(host_subject_id, .keep_all = TRUE)

abx_counts <- abx_drugs %>%
  group_by(AgeGroup, BiologicalSex) %>%
  summarise(DrugCount = n(), .groups = 'drop') 

normalized_data_6 <- abx_counts %>%
  left_join(total_sample_count, by = c("AgeGroup", "BiologicalSex")) %>%
  mutate(normalized = DrugCount / TotalSampleCount)

merged_data_6 <- merge(all_combinations, normalized_data_6, by = c("AgeGroup", "BiologicalSex"), all = TRUE) %>%
  mutate(normalized = ifelse(is.na(normalized), 0, normalized))

abx_normalized <- ggplot(merged_data_6, aes(x = AgeGroup, y = normalized, fill = BiologicalSex)) +
  geom_col(position = "dodge", alpha=0.5) +
  geom_smooth(aes(group = BiologicalSex, color = BiologicalSex), method = "loess", span = 0.5, se = F, linetype = "solid", size = 2) +
  labs(title = "Antibiotics",
       x = "Age (yrs)",
       y = "Normalized Count", 
       fill = "") +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),  
    axis.line = element_line(color = "black"),  
    axis.ticks = element_line(color = "black"),  
    axis.ticks.length = unit(0.2, "cm"),  
    axis.text = element_text(color = "black"),  
    axis.title = element_text(color = "black"),  
    plot.title = element_text(hjust = 0.5),  
    axis.text.x = element_text(angle = 45, hjust = 1))
abx_normalized
#ggsave("AGP_antibiotics.svg", plot = abx_normalized, width = 6, height = 4, dpi = 900)
#getwd()

###### SI FIGURE S4b - WORLDMAP ######
# Filter to keep only individuals with available geographic info
AGP_10317_with_region <- AGP_10317 %>%
  dplyr::filter(!(LatitudeandLongitude %in% c("not specified", "not specified|not specified", "0"))) %>%
  dplyr::select(SampleID, Drug, Country, LatitudeandLongitude, host_subject_id, Region)

country_counts <- AGP_10317 %>%
  group_by(Country) %>%
  summarize(unique_host_subjects = n_distinct(host_subject_id)) %>%
  arrange(desc(unique_host_subjects))

unique_IDs_AGP <- AGP_10317_with_region %>%
  dplyr::distinct(host_subject_id, .keep_all = TRUE) # n = 1903

# Look at data where drug was detected and geographic data is available
AGP_druglib_with_region <- AGP_druglib %>%
  dplyr::filter(!(LatitudeandLongitude %in% c("not specified", "not specified|not specified", "0"))) %>%
  dplyr::select(SampleID, name_parent_compound, Country, LatitudeandLongitude, host_subject_id, Region)

# Count of drugs per individual
AGP_druglib_with_region <- AGP_druglib_with_region %>%
  group_by(host_subject_id) %>%
  dplyr::mutate(count = n_distinct(name_parent_compound)) %>%
  ungroup() %>%
  dplyr::mutate(count = ifelse(count > 4, 4, count)) 

AGP_IDs_AGP_drug <- AGP_druglib_with_region %>%
  dplyr::distinct(host_subject_id, .keep_all = TRUE)  %>%
  dplyr::select(!name_parent_compound)

common_ids <- intersect(AGP_druglib_with_region$host_subject_id, AGP_10317_with_region$host_subject_id)

# Filter out rows from AGP_10317_with_region where the host_subject_id is present in both datasets
df1_filtered <- AGP_10317_with_region %>%
  dplyr::filter(!host_subject_id %in% common_ids)

df1_filtered$count <- 0
unique_df1_filtered <- df1_filtered %>%
  dplyr::distinct(host_subject_id, .keep_all = TRUE) %>%
  dplyr::select(!Drug)

unique_country <- rbind(unique_df1_filtered, AGP_IDs_AGP_drug)

# Define color palette and shapes for drug count
colors <- brewer.pal(5, "Paired")
print(colors)
colors2 <- c(colors[5], colors[1], colors[3], colors[2], colors[4])
print(colors)
shapes <- c(2, 17, 15, 7, 16) 

# Separate Latitude and Longitude
unique_country2 <- unique_country %>%
  separate(LatitudeandLongitude, into = c("lat", "long"), sep = "[|;,]", convert = TRUE)

world_map <- map_data("world")

# United States
worldmap_AGP_US <- ggplot() +
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group), fill = "#F0F0F0", color = "#F0F0F0") +
  geom_point(data = unique_country2, aes(x = long, y = lat, shape = factor(count), color = factor(count)), size = 2.5, alpha = 0.6) +
  scale_color_manual(values = colors2) +
  scale_shape_manual(values = shapes) +
  theme_minimal() +
  theme(
    axis.line = element_line(color = "black"),  
    axis.ticks = element_line(color = "black"), 
    axis.ticks.length = unit(0.2, "cm"),  
    axis.text = element_text(color = "black"), 
    axis.title = element_text(color = "black"),  
    plot.title = element_text(hjust = 0.5),  
    legend.position = "bottom") +
  labs(title = "Drug Count by Geographic Location", x = "Longitude", y = "Latitude", color = "Drug Count", shape = "Drug Count") +
  coord_sf(xlim = c(-125, -67), ylim = c(23, 49))
print(worldmap_AGP_US)
#ggsave("AGP_worldmap_US.svg", plot = worldmap_AGP_US, width = 6, height = 4, dpi = 900)

# United Kingdom
worldmap_AGP_UK <- ggplot() +
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group), fill = "#F0F0F0", color = "#F0F0F0") +
  geom_point(data = unique_country2, aes(x = long, y = lat, shape = factor(count), color = factor(count)), size = 2.5, alpha = 0.6) +
  scale_color_manual(values = colors2) +
  scale_shape_manual(values = shapes) +
  theme_minimal() +
  theme(
    axis.line = element_line(color = "black"),  
    axis.ticks = element_line(color = "black"),  
    axis.ticks.length = unit(0.2, "cm"),  
    axis.text = element_text(color = "black"),  
    axis.title = element_text(color = "black"),  
    plot.title = element_text(hjust = 0.5),  
    legend.position = "bottom") +
  labs(title = "Drug Count by Geographic Location", x = "Longitude", y = "Latitude", color = "Drug Count", shape = "Drug Count") +
  coord_sf(xlim = c(-7, 2), ylim = c(49, 60))
print(worldmap_AGP_UK)
#ggsave("AGP_worldmap_Europe.svg", plot = worldmap_AGP_UK, width = 6, height = 4, dpi = 900)

# Australia
worldmap_AGP_australia <- ggplot() +
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group), fill = "#F0F0F0", color = "#F0F0F0") +
  geom_point(data = unique_country2, aes(x = long, y = lat, shape = factor(count), color = factor(count)), size = 2.5, alpha = 0.6) +
  scale_color_manual(values = colors2) +
  scale_shape_manual(values = shapes) +
  theme_minimal() +
  theme(
    axis.line = element_line(color = "black"),  
    axis.ticks = element_line(color = "black"),  
    axis.ticks.length = unit(0.2, "cm"),  
    axis.text = element_text(color = "black"),  
    axis.title = element_text(color = "black"),  
    plot.title = element_text(hjust = 0.5),  
    legend.position = "bottom") +
  labs(title = "Drug Count by Geographic Location", x = "Longitude", y = "Latitude", color = "Drug Count", shape = "Drug Count") +
  coord_sf(xlim = c(112, 155), ylim = c(-39, -12))
print(worldmap_AGP_australia)
#ggsave("AGP_worldmap_Australia.svg", plot = worldmap_AGP_australia, width = 5, height = 3.2, dpi = 900)

# Calculate the percentage of individuals with each drug count for the world map 
drug_count_freq <- unique_country %>%
  group_by(count) %>%
  summarise(freq = n())

total_individuals_2 <- n_distinct(unique_country$host_subject_id)
drug_count_fre_2 <- drug_count_freq %>%
  mutate(percentage = (freq / total_individuals_2) * 100)

# Count by individual region
drug_count_freq_3 <- unique_country %>%
  group_by(count, Region) %>%
  summarise(freq = n())

data <- data.frame(
  Region = c("Australia", "Europe", "United States", "Australia", "Europe", "United States", 
             "Australia", "Europe", "United States", 
             "Australia", "Europe", "United States", "Australia", "Europe", "United States"),
  count = c(0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4),
  freq = c(92, 400, 771, 18, 94, 277, 9, 28, 109, 1, 7, 53, 5, 4, 35)
)

# Create a frequency matrix
freq_matrix <- xtabs(freq ~ Region + count, data = data)

# Perform chi-square test
chisq.test(freq_matrix)

# Calculate the percentage of individuals with each drug count for the US
US <- unique_country %>%
  dplyr::filter((Region %in% c("United States")))

drug_count_freq_US <- US %>%
  group_by(count) %>%
  summarise(freq = n())

drug_count_freq_US <- drug_count_freq_US %>%
  mutate(percentage = (freq / 1245) * 100)

# Calculate the percentage of individuals with each drug count for Europe
Europe <- unique_country %>%
  dplyr::filter((Region %in% c("Europe")))

drug_count_freq_Europe <- Europe %>%
  group_by(count) %>%
  summarise(freq = n())

drug_count_freq_Europe <- drug_count_freq_Europe %>%
  mutate(percentage = (freq / 533) * 100)

# Calculate the percentage of individuals with each drug count for Australia
Australia <- unique_country %>%
  dplyr::filter((Region %in% c("Australia")))

drug_count_freq_Australia <- Australia %>%
  group_by(count) %>%
  summarise(freq = n())

drug_count_freq_Australia <- drug_count_freq_Australia %>%
  mutate(percentage = (freq / 125) * 100)

#
UK <- unique_country %>%
  dplyr::filter((Country %in% c("United Kingdom"))) %>%
  dplyr::distinct(host_subject_id, .keep_all = TRUE)

