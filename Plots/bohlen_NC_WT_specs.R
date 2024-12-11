sample_names = c( #TOCHANGE
  "SRR10416856",  # 80S WT
  "SRR11553468"  # 40S WT
)

elements = c("cds_end", "cds_start", "stst", "u9", "u9_end")

sample80 <- read.csv("sample80_WT.csv")
sample40 <- read.csv("sample40_WT.csv")

.sample_80S = c("SRR10416856")
.sample_40S = c("SRR11553468")

count_fraction = 0.5

read_data_count <- function(ribo_type, element) {
  samples = switch (ribo_type,
                    "80S" = .sample_80S,
                    "40S" = .sample_40S
  )
  
  sample1 = read.table(paste0("transcript-count/", samples[1], "_read_to_", tolower(element), ".tsv"), header = T) %>%
    mutate(across(3:4, ~ replace_na(., 0)))
  
  sample = sample1
  
  return(sample)
}

read_data_transcript <- function(ribo_type, element) {
  samples = switch (ribo_type,
                    "80S" = .sample_80S,
                    "40S" = .sample_40S)
  
  sample1 = read.table(paste0("transcript-data/", samples[1], "_read_to_", tolower(element), ".tsv"), header = T) %>%
    mutate(across(3:4, ~ replace_na(., 0))) 
  
  sample = sample1 %>%
    group_by(distance2element, length) %>%
    summarise_all(sum) %>%
    ungroup()
  
  return(sample)
}

read_data_3 <- function(ribo_type, element) {
  samples = switch (ribo_type,
                    "80S" = .sample_80S,
                    "40S" = .sample_40S
  )
  
  sample1 = read.table(paste0("transcript-3prime/", samples[1], "_read_to_", tolower(element), ".tsv"), header = T) %>%
    mutate(across(3:4, ~ replace_na(., 0))) 
  
  sample = sample1 %>%
    group_by(distance2element, length) %>%
    summarise_all(sum) %>%
    ungroup()
  
  return(sample)
}


