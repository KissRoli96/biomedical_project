# Biostatisztika final project - imputacios feladat
#
# Ez a script az oktato altal kert imputacios outputokhoz keszult:
# 1) Kolmogorov-Smirnov (KS) osszehasonlito grafikon
# 2) Mean Absolute Error (MAE) osszehasonlito grafikon
# 3) kivalasztott imputacios algoritmus
# 4) egy valtozo eredeti es imputalt eloszlasanak histogramja

# Ha valamelyik csomag hianyzik, RStudio Console-ban futtasd:
# install.packages(c("dplyr", "ggplot2", "naniar", "missCompare"))

library(dplyr)
library(ggplot2)
library(naniar)
library(missCompare)


# 1. Adat beolvasasa ---------------------------------------------------------

adult_data <- read.csv("adult22csv/adult22.csv", stringsAsFactors = FALSE)

phq_items <- paste0("PHQ8", 1:8, "_A")
gad_items <- paste0("GAD7", 1:7, "_A")

project_data <- adult_data[, c(
  phq_items,
  gad_items,
  "PHQCAT_A",
  "GADCAT_A",
  "AGEP_A",
  "SEX_A",
  "PHSTAT_A",
  "SMKCIGST_A",
  "DRKSTAT_A",
  "PA18_02R_A",
  "BMICAT_A",
  "EDUCP_A"
)]


# 2. Adattisztitas es specialis hianyzo kodok NA-ra alakitasa ---------------

recode_item <- function(x) {
  ifelse(x %in% 1:4, x - 1, NA)
}

analysis_data <- project_data

# PHQ-8 es GAD-7 itemek: 1-4 valaszbol 0-3 pont lesz, 7/8/9 -> NA
analysis_data[phq_items] <- lapply(analysis_data[phq_items], recode_item)
analysis_data[gad_items] <- lapply(analysis_data[gad_items], recode_item)

# Kategorias es kovarians valtozok: ervenyes kodok megtartasa, tobbi NA
analysis_data <- analysis_data %>%
  mutate(
    phq_category = ifelse(PHQCAT_A %in% 1:4, PHQCAT_A, NA),
    gad_category = ifelse(GADCAT_A %in% 1:4, GADCAT_A, NA),
    age = ifelse(AGEP_A >= 18 & AGEP_A <= 85, AGEP_A, NA),
    sex = ifelse(SEX_A %in% 1:2, SEX_A, NA),
    health_status = ifelse(PHSTAT_A %in% 1:5, PHSTAT_A, NA),
    smoking = ifelse(SMKCIGST_A %in% 1:4, SMKCIGST_A, NA),
    alcohol = ifelse(DRKSTAT_A %in% 1:8, DRKSTAT_A, NA),
    physical_activity = ifelse(PA18_02R_A %in% 1:3, PA18_02R_A, NA),
    bmi_cat = ifelse(BMICAT_A %in% 1:4, BMICAT_A, NA),
    education = ifelse(EDUCP_A %in% 1:10, EDUCP_A, NA)
  ) %>%
  select(
    all_of(phq_items),
    all_of(gad_items),
    phq_category,
    gad_category,
    age,
    sex,
    health_status,
    smoking,
    alcohol,
    physical_activity,
    bmi_cat,
    education
  )


# 3. Hianyzo adatok felt terkepezese ----------------------------------------

missing_table <- data.frame(
  variable = names(analysis_data),
  missing_n = colSums(is.na(analysis_data)),
  missing_percent = round(colMeans(is.na(analysis_data)) * 100, 2)
)

print(missing_table)

# Egyszeru hianyzo adat abra
gg_miss_var(analysis_data)

# Reszletesebb vizualizacio
vis_miss(analysis_data)


# 4. missCompare elokeszites ------------------------------------------------

# A missCompare::clean() eltavolithat nagyon sok hianyt tartalmazo valtozokat
# vagy sorokat. A kuszobertekeket a kurzusanyaghoz hasonloan allitjuk be.
cleaned <- missCompare::clean(
  analysis_data,
  var_removal_threshold = 0.5,
  ind_removal_threshold = 0.8
)

# MCAR teszt: azt vizsgalja, hogy a hianyzas teljesen veletlenszeru-e.
mcar_test(cleaned)

metadata <- missCompare::get_data(
  cleaned,
  matrixplot_sort = TRUE,
  plot_transform = TRUE
)

metadata$Complete_cases
metadata$NA_Correlation_plot
metadata$Fraction_missingnes
metadata$Fraction_missingness_per_variable
metadata$NA_per_variable
metadata$min_PDM_thresholds
metadata$Vars_above_half


# 5. Imputacios algoritmusok osszehasonlitasa -------------------------------

# Ez a lepes adja a KS es MAE osszehasonlito eredmenyeket/grafikonokat.
# Minel kisebb a KS es MAE, annal jobb az imputacios algoritmus.

simulated <- missCompare::simulate(
  rownum = metadata$Rows,
  colnum = metadata$Columns,
  cormat = metadata$Corr_matrix,
  meanval = 0,
  sdval = 5
)

imputation_comparison <- missCompare::impute_simulated(
  rownum = metadata$Rows,
  colnum = metadata$Columns,
  cormat = metadata$Corr_matrix,
  MD_pattern = metadata$MD_Pattern,
  NA_fraction = metadata$Fraction_missingness,
  min_PDM = 5,
  n.iter = 5,
  assumed_pattern = NA
)

# Az imputation_comparison objektum tartalmazza az algoritmusok
# osszehasonlitasahoz szukseges eredmenyeket. RStudio-ban kattints ra az
# Environment panelen, es nezd meg, milyen nevu elemekben vannak a KS es MAE
# grafikonok/táblák.
imputation_comparison


# 6. Kivalasztott algoritmus -------------------------------------------------

# A kurzusanyagban a missForest algoritmus szerepel peldakent.
# Ha a KS es MAE grafikonok alapjan mas algoritmus jobb, a sel_method erteket
# ahhoz kell igazitani.
#
# A missCompare dokumentacio szerint a method kodokat itt tudod megnezni:
# ?missCompare::impute_data
#
# Kiindulo valasztas:
# sel_method = c(14)  # missForest

selected_algorithm <- "missForest"
selected_method_code <- c(14)

imputed_missForest <- missCompare::impute_data(
  cleaned,
  scale = FALSE,
  n.iter = 10,
  sel_method = selected_method_code
)

# A kurzusanyag alapjan a missForest imputalt adat altalaban itt van:
imputed_data <- imputed_missForest$missForest_imputation[[1]]


# 7. Poszt-imputacios diagnosztika es histogram -----------------------------

diag <- missCompare::post_imp_diag(
  cleaned,
  imputed_data,
  scale = FALSE,
  n.boot = 5
)

# Pelda: egy PHQ item eredeti es imputalt eloszlasanak histogramja.
# Ha mas valtozot szeretnel bemutatni, csereld ki a valtozo nevet.
diag$Histograms$PHQ83_A

# Alternativ, sajat ggplot histogram egy valtozora:
histogram_variable <- "PHQ83_A"

hist_data <- data.frame(
  value = c(cleaned[[histogram_variable]], imputed_data[[histogram_variable]]),
  source = c(
    rep("Eredeti adat", nrow(cleaned)),
    rep("Imputalt adat", nrow(imputed_data))
  )
)

ggplot(hist_data, aes(x = value, fill = source)) +
  geom_histogram(position = "identity", alpha = 0.5, bins = 10) +
  labs(
    title = paste("Eredeti es imputalt eloszlas:", histogram_variable),
    x = histogram_variable,
    y = "Gyakorisag",
    fill = "Adat"
  ) +
  theme_minimal()


# 8. Rovid eredmenyszoveg sablon --------------------------------------------

cat("\nKivalasztott imputacios algoritmus:", selected_algorithm, "\n")
cat("A vegleges valasztast a KS es MAE osszehasonlito grafikonok alapjan kell indokolni.\n")
cat("A kisebb KS es MAE ertek jobb imputacios teljesitmenyt jelez.\n")

