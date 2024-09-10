# processing data from the ferrozine assay
# kfp, 2024-03-01

library(tidyverse)
library(googlesheets4)

ferrozine_weights = read_sheet("https://docs.google.com/spreadsheets/d/1bpRgvOgWwtvSYZS7CHiNt-2Ck3G5AVg4tGnq_zWsJ90/edit#gid=873884422", sheet = "iron-weights")
import_iron = function(FILEPATH = "1-data/iron-ferrozine"){
  
  # import map
  ferrozine_map = read_sheet("https://docs.google.com/spreadsheets/d/1bpRgvOgWwtvSYZS7CHiNt-2Ck3G5AVg4tGnq_zWsJ90/edit#gid=873884422", sheet = "iron-map", col_types = "c") %>% janitor::clean_names()
  
  # import data files (plate reader)
  filePaths_ferrozine <- list.files(path = FILEPATH, pattern = ".xlsx", full.names = TRUE, recursive = TRUE)
  ferrozine_data <- do.call(bind_rows, lapply(filePaths_ferrozine, function(path) {
    df <- readxl::read_excel(path, skip = 25) %>% mutate_all(as.character) %>% janitor::clean_names()
    df = df %>% mutate(source = basename(path))
    df}))
  
  list(ferrozine_map = ferrozine_map,
       ferrozine_data = ferrozine_data)
  
}
process_iron = function(ferrozine_map, ferrozine_data){
  
  # clean the map
  map_processed = 
    ferrozine_map %>% 
    mutate(date = ymd(date)) %>% 
    fill(date, tray, analysis) %>% 
    pivot_longer(-c(date, tray, analysis, dilution, letter, notes), names_to = "number", values_to = "sample_label") %>% 
    filter(!is.na(sample_label)) %>% 
    mutate(number = parse_number(number)) %>% 
    mutate_at(vars(c(dilution, number)), as.numeric) %>% 
    mutate(well_position = paste0(letter, number),
           dilution = if_else(is.na(dilution), 1, dilution)) %>% 
    arrange(date, tray, number, letter) %>% 
    mutate(sample_type = case_when(grepl("mM", sample_label) ~ "standard",
                                   TRUE ~ "sample")) %>% 
    dplyr::select(date, tray, analysis, dilution, well_position, sample_label, sample_type)
  
  # clean the data
  data_processed = 
    ferrozine_data %>% 
    mutate_all(na_if,"") %>% 
    fill(x1) %>% 
    filter(x14 == "562") %>% 
    rename(x = x1,
           x1 = x1_2) %>% 
    pivot_longer(-c(source, x), values_to = "absorbance_562") %>% 
    mutate(name = str_remove(name, "x"),
           well_position = paste0(x, name),
           date = str_extract(source, "[0-9]{4}-[0-9]{2}-[0-9]{2}"),
           date = ymd(date),
           tray = str_extract(source, "plate[1-9][a-z]?"),
           tray = parse_number(tray, "plate"),
           tray = as.character(tray),
           absorbance_562 = as.numeric(absorbance_562)) %>% 
    dplyr::select(date, tray, well_position, absorbance_562) %>% 
    right_join(map_processed, by = c("date", "tray", "well_position")) %>% 
    filter(!grepl("skip", sample_label))
  
  calibrate_ferrozine_data = function(data_processed){
    # now do the calibrations
    # standards are in mM
    # molecular formula for FAS = (NH₄)₂Fe(SO₄)₂·6H₂O
    # so 1 M FAS = 1M Fe
    # 1 mM FAS = 1 * 55.85 mg Fe in 1 L solution = 55.85 mg Fe in 1 L solution
    # therefore 1 mM = 55.85 mg/L or 55.85 ppm
    
    standards = 
      data_processed %>% 
      filter(sample_type == "standard") %>% 
      mutate(standard_mM = parse_number(sample_label),
             standard_type = str_extract(sample_label, "FAS|FeCl3"),
             standard_ppm =  case_when(standard_type == "FAS" ~ standard_mM * 55.85)) %>% 
      #dplyr::select(date, tray, absorbance_562, standard_ppm) %>% 
      mutate(standard_ppm = as.numeric(standard_ppm))
    
    # reduction efficiency
    reduction_efficiency = 
      standards %>% filter(standard_mM == 2) %>% 
      group_by(date, standard_type) %>% 
      dplyr::summarize(abs = mean(absorbance_562)) %>% 
      pivot_wider(names_from = "standard_type", values_from = "abs") %>% 
      mutate(red_eff = 100 * FeCl3/FAS) %>% 
      dplyr::select(-FAS, -FeCl3)
    # 99.6 % efficiency

    gg_calibration = 
      standards %>% 
      ggplot(aes(x = standard_ppm, y = absorbance_562, color = as.character(tray)))+
      geom_point()+
      geom_smooth(method = "lm", se = F)+
      facet_wrap(~date + tray)
    
    calibration_coef = 
      standards %>% 
      drop_na() %>% 
      dplyr::group_by(date) %>% 
      dplyr::summarize(slope = lm(absorbance_562 ~ standard_ppm)$coefficients["standard_ppm"], 
                       intercept = lm(absorbance_562 ~ standard_ppm)$coefficients["(Intercept)"])
    
    # y = mx + c
    # abs = m*ppm + c
    # ppm = abs-c/m
    
    data_calibrated = 
      data_processed %>% 
      left_join(calibration_coef, by = c("date")) %>% 
      mutate(ppm_calculated = ((absorbance_562 - intercept) / slope))
    
    list(calibration_coef = calibration_coef,
         data_calibrated = data_calibrated,
         gg_calibration = gg_calibration,
         reduction_efficiency = reduction_efficiency)
  }
  
  calibration_curves = calibrate_ferrozine_data(data_processed)$gg_calibration
  reduction = calibrate_ferrozine_data(data_processed)$reduction_efficiency
  
  samples = 
    calibrate_ferrozine_data(data_processed)$data_calibrated %>% 
    filter(sample_type == "sample") %>% 
    mutate(ppm_calculated = ppm_calculated * dilution) %>% 
    dplyr::select(date, sample_label, analysis, ppm_calculated) %>% 
    filter(!is.na(ppm_calculated)) %>% 
    group_by(date, sample_label) %>% 
    dplyr::summarise(Fe_total = mean(ppm_calculated)) %>% 
    left_join(reduction) %>% 
    mutate(Fe_total = Fe_total * 100/red_eff,  # correct for 99.6% reduction efficiency of ascorbic acid
           Fe_total = round(Fe_total, 2)) %>% 
    rename(FeTotal_ppm = Fe_total) %>% 
    ungroup() %>% 
    dplyr::select(-date, -red_eff)

  # convert to ug/g
  samples_ugg = 
    samples %>% 
    left_join(ferrozine_weights %>% dplyr::select(sample_label, wt_g, HCl_mL)) %>% 
    mutate(Fe_ugg = FeTotal_ppm * HCl_mL/wt_g,
           Fe_ugg = round(Fe_ugg, 2))  %>% 
    dplyr::select(sample_label, Fe_ugg)
}

#
# hog island ----
process_iron_hog = function(ferrozine_map, ferrozine_data){
  
  # clean the map
  map_processed = 
    ferrozine_map %>% 
    mutate(date = ymd(date)) %>% 
    fill(date, tray, analysis) %>% 
    pivot_longer(-c(date, tray, analysis, dilution, letter, notes), names_to = "number", values_to = "sample_label") %>% 
    filter(!is.na(sample_label)) %>% 
    mutate(number = parse_number(number)) %>% 
    mutate_at(vars(c(dilution, number)), as.numeric) %>% 
    mutate(well_position = paste0(letter, number),
           dilution = if_else(is.na(dilution), 1, dilution)) %>% 
    arrange(date, tray, number, letter) %>% 
    mutate(sample_type = case_when(grepl("mM", sample_label) ~ "standard",
                                   TRUE ~ "sample")) %>% 
    dplyr::select(date, tray, analysis, dilution, well_position, sample_label, sample_type)
  
  # clean the data
  data_processed = 
    ferrozine_data %>% 
    mutate_all(na_if,"") %>% 
    fill(x1) %>% 
    filter(x14 == "562") %>% 
    rename(x = x1,
           x1 = x1_2) %>% 
    pivot_longer(-c(source, x), values_to = "absorbance_562") %>% 
    mutate(name = str_remove(name, "x"),
           well_position = paste0(x, name),
           date = str_extract(source, "[0-9]{4}-[0-9]{2}-[0-9]{2}"),
           date = ymd(date),
           tray = str_extract(source, "plate[1-9][a-z]?"),
           tray = parse_number(tray, "plate"),
           tray = as.character(tray),
           absorbance_562 = as.numeric(absorbance_562)) %>% 
    dplyr::select(date, tray, well_position, absorbance_562) %>% 
    right_join(map_processed, by = c("date", "tray", "well_position")) %>% 
    filter(!grepl("skip", sample_label))
  
  calibrate_ferrozine_data = function(data_processed){
    # now do the calibrations
    # standards are in mM
    # molecular formula for FAS = (NH₄)₂Fe(SO₄)₂·6H₂O
    # so 1 M FAS = 1M Fe
    # 1 mM FAS = 1 * 55.85 mg Fe in 1 L solution = 55.85 mg Fe in 1 L solution
    # therefore 1 mM = 55.85 mg/L or 55.85 ppm
    
    standards = 
      data_processed %>% 
      filter(sample_type == "standard") %>% 
      mutate(standard_mM = parse_number(sample_label),
             standard_type = str_extract(sample_label, "FAS|FeCl3"),
             standard_ppm =  case_when(standard_type == "FAS" ~ standard_mM * 55.85)) %>% 
      #dplyr::select(date, tray, absorbance_562, standard_ppm) %>% 
      mutate(standard_ppm = as.numeric(standard_ppm))
    
    # reduction efficiency
    reduction_efficiency = 
      standards %>% filter(standard_mM == 2) %>% 
      group_by(date, standard_type) %>% 
      dplyr::summarize(abs = mean(absorbance_562)) %>% 
      pivot_wider(names_from = "standard_type", values_from = "abs") %>% 
      mutate(red_eff = 100 * FeCl3/FAS) %>% 
      dplyr::select(-FAS, -FeCl3)
    # 99.6 % efficiency
    
    gg_calibration = 
      standards %>% 
      ggplot(aes(x = standard_ppm, y = absorbance_562, color = as.character(tray)))+
      geom_point()+
      geom_smooth(method = "lm", se = F)+
      facet_wrap(~date + tray)
    
    calibration_coef = 
      standards %>% 
      drop_na() %>% 
      dplyr::group_by(date) %>% 
      dplyr::summarize(slope = lm(absorbance_562 ~ standard_ppm)$coefficients["standard_ppm"], 
                       intercept = lm(absorbance_562 ~ standard_ppm)$coefficients["(Intercept)"])
    
    # y = mx + c
    # abs = m*ppm + c
    # ppm = abs-c/m
    
    data_calibrated = 
      data_processed %>% 
      left_join(calibration_coef, by = c("date")) %>% 
      mutate(ppm_calculated = ((absorbance_562 - intercept) / slope))
    
    list(calibration_coef = calibration_coef,
         data_calibrated = data_calibrated,
         gg_calibration = gg_calibration,
         reduction_efficiency = reduction_efficiency)
  }
  
  calibration_curves = calibrate_ferrozine_data(data_processed)$gg_calibration
  reduction = calibrate_ferrozine_data(data_processed)$reduction_efficiency
  
  samples = 
    calibrate_ferrozine_data(data_processed)$data_calibrated %>% 
    filter(sample_type == "sample") %>% 
    mutate(ppm_calculated = ppm_calculated * dilution) %>% 
    dplyr::select(date, sample_label, analysis, ppm_calculated) %>% 
    filter(!is.na(ppm_calculated)) %>% 
    group_by(date, sample_label) %>% 
    dplyr::summarise(Fe_total = mean(ppm_calculated)) %>% 
    left_join(reduction) %>% 
    mutate(Fe_total = Fe_total * 100/red_eff,  # correct for 99.6% reduction efficiency of ascorbic acid
           Fe_total = round(Fe_total, 2)) %>% 
    rename(FeTotal_ppm = Fe_total) %>% 
    ungroup() %>% 
    dplyr::select(-date, -red_eff)
  
  # convert to ug/g
  samples_ugg = 
    samples %>% 
    left_join(ferrozine_weights %>% dplyr::select(sample_label, wt_g, HCl_mL)) %>% 
    mutate(Fe_ugg = FeTotal_ppm * HCl_mL/wt_g,
           Fe_ugg = round(Fe_ugg, 2))  %>% 
    dplyr::select(sample_label, Fe_ugg)
}

data_hog = import_iron(FILEPATH = "1-data/iron-ferrozine/2024-02-29_hog_island")$ferrozine_data
map_hog = import_iron(FILEPATH = "1-data/iron-ferrozine/2024-02-29_hog_island")$ferrozine_map

processed_hog = process_iron_hog(ferrozine_map = map_hog, ferrozine_data = data_hog)
processed_hog %>% write.csv("1-data/processed/hog_iron.csv", row.names = F, na = "")

# upland
process_iron_upland = function(ferrozine_map, ferrozine_data){
  
  # clean the map
  map_processed = 
    ferrozine_map %>% 
    mutate(date = ymd(date)) %>% 
    fill(date, tray, analysis) %>% 
    pivot_longer(-c(date, tray, analysis, dilution, letter, notes), names_to = "number", values_to = "sample_label") %>% 
    filter(!is.na(sample_label)) %>% 
    mutate(number = parse_number(number)) %>% 
    mutate_at(vars(c(dilution, number)), as.numeric) %>% 
    mutate(well_position = paste0(letter, number),
           dilution = if_else(is.na(dilution), 1, dilution)) %>% 
    arrange(date, tray, number, letter) %>% 
    mutate(sample_type = case_when(grepl("mM", sample_label) ~ "standard",
                                   TRUE ~ "sample")) %>% 
    dplyr::select(date, tray, analysis, dilution, well_position, sample_label, sample_type)
  
  # clean the data
  data_processed = 
    ferrozine_data %>% 
    mutate_all(na_if,"") %>% 
    fill(x1) %>% 
    filter(x14 == "562") %>% 
    rename(x = x1,
           x1 = x1_2) %>% 
    pivot_longer(-c(source, x), values_to = "absorbance_562") %>% 
    mutate(name = str_remove(name, "x"),
           well_position = paste0(x, name),
           date = str_extract(source, "[0-9]{4}-[0-9]{2}-[0-9]{2}"),
           date = ymd(date),
           tray = str_extract(source, "plate[1-9][a-z]?"),
           tray = parse_number(tray, "plate"),
           tray = as.character(tray),
           absorbance_562 = as.numeric(absorbance_562)) %>% 
    dplyr::select(date, tray, well_position, absorbance_562) %>% 
    right_join(map_processed, by = c("date", "tray", "well_position")) %>% 
    filter(!grepl("skip", sample_label))
  
  calibrate_ferrozine_data = function(data_processed){
    # now do the calibrations
    # standards are in mM
    # molecular formula for FAS = (NH₄)₂Fe(SO₄)₂·6H₂O
    # so 1 M FAS = 1M Fe
    # 1 mM FAS = 1 * 55.85 mg Fe in 1 L solution = 55.85 mg Fe in 1 L solution
    # therefore 1 mM = 55.85 mg/L or 55.85 ppm
    
    standards = 
      data_processed %>% 
      filter(sample_type == "standard") %>% 
      mutate(standard_mM = parse_number(sample_label),
             standard_type = str_extract(sample_label, "FAS|FeCl3"),
             standard_ppm =  case_when(standard_type == "FAS" ~ standard_mM * 55.85)) %>% 
      #dplyr::select(date, tray, absorbance_562, standard_ppm) %>% 
      mutate(standard_ppm = as.numeric(standard_ppm))
    
    # reduction efficiency
    reduction_efficiency = 
      standards %>% filter(standard_mM == 1) %>% 
      group_by(date, standard_type) %>% 
      dplyr::summarize(abs = mean(absorbance_562)) %>% 
      pivot_wider(names_from = "standard_type", values_from = "abs") %>% 
      mutate(red_eff = 100 * FeCl3/FAS) %>% 
      dplyr::select(-FAS, -FeCl3)
    # 99.6 % efficiency
    
    gg_calibration = 
      standards %>% 
      ggplot(aes(x = standard_ppm, y = absorbance_562, color = as.character(tray)))+
      geom_point()+
      geom_smooth(method = "lm", se = F)+
      facet_wrap(~date + tray)
    
    calibration_coef = 
      standards %>% 
      drop_na() %>% 
      dplyr::group_by(date) %>% 
      dplyr::summarize(slope = lm(absorbance_562 ~ standard_ppm)$coefficients["standard_ppm"], 
                       intercept = lm(absorbance_562 ~ standard_ppm)$coefficients["(Intercept)"])
    
    # y = mx + c
    # abs = m*ppm + c
    # ppm = abs-c/m
    
    data_calibrated = 
      data_processed %>% 
      left_join(calibration_coef, by = c("date")) %>% 
      mutate(ppm_calculated = ((absorbance_562 - intercept) / slope))
    
    list(calibration_coef = calibration_coef,
         data_calibrated = data_calibrated,
         gg_calibration = gg_calibration,
         reduction_efficiency = reduction_efficiency)
  }
  
  calibration_curves = calibrate_ferrozine_data(data_processed)$gg_calibration
  reduction = calibrate_ferrozine_data(data_processed)$reduction_efficiency
  
  samples = 
    calibrate_ferrozine_data(data_processed)$data_calibrated %>% 
    filter(sample_type == "sample") %>% 
    mutate(ppm_calculated = ppm_calculated * dilution) %>% 
    dplyr::select(date, sample_label, analysis, ppm_calculated) %>% 
    filter(!is.na(ppm_calculated)) %>% 
    group_by(date, sample_label) %>% 
    dplyr::summarise(Fe_total = mean(ppm_calculated)) %>% 
    left_join(reduction) %>% 
    mutate(Fe_total = Fe_total * 100/red_eff,  # correct for 99.6% reduction efficiency of ascorbic acid
           Fe_total = round(Fe_total, 2)) %>% 
    rename(FeTotal_ppm = Fe_total) %>% 
    ungroup() %>% 
    dplyr::select(-date, -red_eff)
  
  # convert to ug/g
  samples_ugg = 
    samples %>% 
    left_join(ferrozine_weights %>% dplyr::select(sample_label, wt_g, HCl_mL)) %>% 
    mutate(Fe_ugg = FeTotal_ppm * HCl_mL/wt_g,
           Fe_ugg = round(Fe_ugg, 2))  %>% 
    dplyr::select(sample_label, Fe_ugg)
}

data_upland = import_iron(FILEPATH = "1-data/iron-ferrozine/2024-07-22_upland")$ferrozine_data
map_upland = import_iron(FILEPATH = "1-data/iron-ferrozine/2024-07-22_upland")$ferrozine_map

processed_upland = process_iron_upland(ferrozine_map = map_upland, ferrozine_data = data_upland) %>% filter(!is.na(Fe_ugg))
processed_upland %>% write.csv("1-data/processed/upland_iron.csv", row.names = F, na = "")
