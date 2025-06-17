# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# PROJECT: Process Persnet
# PURPOSE: Organize and analyze egonetwork data from the persnet instrument
# DIR:     ""
# INPUTS:  Fake data ("fake_data.csv")
# OUTPUTS: A .csv file of network measures, a grid of sociograms, and individual
# sociograms. Built on the legacy code of persnet analyses.
# Requirements: tidyverse (version 2.0.0), tidygraph (1.3.1), stats (4.1.1) and
# purrr (1.0.2); igraph (2.0.3); patchwork (1.3.0), ggraph, patchwork, cowplot,
# rstudioapi
# AUTHORS: Zachary Wehrwein, Liam McCafferty, Amar Dhand
# CREATED: 03/26/2025
# LATEST:  06/08/2025
# PSERIES: NA
# NSERIES: NA
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# This file will not run if the following occurs:
# 1. igraph, tidyverse, tidygraph go through a major update
# 2. the way in which persnet data is loaded changes (e.g. character vs numeric)
# 3. different variable names, but instructions to adjust are included below
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

################################# Code Setup ##################################

#Attaching packages
library(tidyverse)
library(igraph)
library(cowplot)
library(tidygraph)
library(purrr)
library(stats)
library(ggraph)
library(patchwork)
library(rstudioapi)

### Locate data
data_file <- "L:/Research Project Current/Social Connectedness/data/persnet/soccon_persnet_demographic_2025-06-10_fakenames.csv"
# data_file <- "fake_data.csv"
dict_file <- "L:/Research Project Current/Social Connectedness/data/persnet/soccon_datadict.csv"

out_dir <- "L:/Research Project Current/Social Connectedness/data/persnet/"
out_suffix <- "fakenames"

#Clearing Working Directory
rm(list = ls())

#Automatically set working directory, requires rstudioapi package. If not working set WD manually
#to do so manually go to Session -> Set Working Directory -> To Source File location
#copy paste result from console over following commented line of code:
#Example:
# setwd("~/Desktop/network")
setwd(dirname(getActiveDocumentContext()$path))

# Get funcs
source("funcs.R")


######################### ANALYSIS ##############################
#Importing data
df_input <- read.csv(
    data_file,
    stringsAsFactors = FALSE,
    colClasses = c(zip = "character")
) #this reads in zip retains 0s

df_dict <- read.csv(dict_file)


### record_id
record_id = df_input$record_id

### age
age = df_input$age

### sex
sex_map <- get_codebook_mapping(df_dict, "sex")
sex <- factor(df_input$sex, levels = sex_map$levels, labels = sex_map$labels)

### race
#purrr applies the same function to a list, in this case, each row of df_input
#the blank space between .x, , is selecting all columns of df_input
#while drop= FALSE forces R to return not a vector but a dataframe object
race_labels <- get_codebook_mapping(df_dict, "race")$labels
ego_race <- lapply(1:4, function(n) {
    {
        purrr::map_chr(
            1:nrow(df_input),
            ~ extract_attribute(
                df_input[.x, , drop = FALSE],
                n,
                race_labels,
                "race"
            )
        )
    } %>%
        factor(race_labels)
})

### education
# Sometimes this is "edu", sometimes "education". Since used twice, set as val here
edu_name <- "education"
edu_map <- get_codebook_mapping(df_dict, edu_name)
df_input$edu <- factor(
    df_input[, edu_name],
    levels = edu_map$levels,
    labels = edu_map$labels
)
education <- df_input$edu

### zip
zip = as.character(df_input$zip)

### smoking
smoke_map <- get_codebook_mapping(df_dict, "smoke")

df_input$smoke <- factor(
    df_input$smoke,
    levels = smoke_map$levels,
    labels = smoke_map$labels
)
smoke = df_input$smoke

### Network stats
network_size = calc_total_alters_df(df_input)
density = calc_egoless_density_df(df_input)

gra_list <- organize_list_tidygraphs(df_input)
constraint = sapply(gra_list, calc_node_constraint)
effsize = sapply(gra_list, calc_node_ens)
mean_degree = sapply(gra_list, calc_egoless_mean_degree)
max_degree = sapply(gra_list, calc_egoless_max_degree)

### Alters data
# All relationship mapping will be the same, so just take the first one
age_sd = purrr::map_dbl(
    1:nrow(df_input),
    ~ sd_age_alters_row(df_input[.x, , drop = FALSE])
)

gender_map <- codebook_to_dict(get_codebook_mapping(df_dict, "name1sex"))
iqv_sex <- purrr::map_dbl(
    1:nrow(df_input),
    ~ calc_attribute_iqv(df_input[.x, , drop = FALSE], "sex", gender_map)
)

race_map <- codebook_to_dict(get_codebook_mapping(df_dict, "name1race"))
iqv_race <- purrr::map_dbl(
    1:nrow(df_input),
    ~ calc_attribute_iqv(df_input[.x, , drop = FALSE], "race", race_map)
)

educ_map <- list(
    only_high_school = c(1, 2),
    some_college = c(3, 4),
    college_grad = c(5, 6),
    dont_know = c(99)
)
iqv_educ <- purrr::map_dbl(
    1:nrow(df_input),
    ~ calc_attribute_iqv(df_input[.x, , drop = FALSE], "educ", educ_map)
)

relat_map <- codebook_to_dict(get_codebook_mapping(df_dict, "name1relat"))
kin_prop <- purrr::map_dbl(
    1:nrow(df_input),
    ~ calc_prop_alters_multians(
        df_input[.x, , drop = FALSE],
        c("Spouse", "Family"),
        relat_map,
        "relat"
    )
)

speak_map <- codebook_to_dict(get_codebook_mapping(df_dict, "name1speak"))
weak_freq_prop <- purrr::map_dbl(
    1:nrow(df_input),
    ~ calc_prop_alters_singleans(
        df_input[.x, , drop = FALSE],
        c("Monthly", "Less often"),
        speak_map,
        "speak",
    )
)

len_map <- codebook_to_dict(get_codebook_mapping(df_dict, "name1length"))
weak_dur_prop <- purrr::map_dbl(
    1:nrow(df_input),
    ~ calc_prop_alters_singleans(
        df_input[.x, , drop = FALSE],
        c("Less than three", "Three to six"),
        len_map,
        "length"
    )
)

dist_map <- codebook_to_dict(get_codebook_mapping(df_dict, "name1dist"))
far_dist_prop <- purrr::map_dbl(
    1:nrow(df_input),
    ~ calc_prop_alters_singleans(
        df_input[.x, , drop = FALSE],
        c("16-50 miles", "50+ miles"),
        len_map,
        "dist"
    )
)

als_map <- codebook_to_dict(get_codebook_mapping(df_dict, "name1als"))
met_through_als_prop <- purrr::map_dbl(
    1:nrow(df_input),
    ~ calc_prop_alters_singleans(
        df_input[.x, , drop = FALSE],
        c("Yes"),
        als_map,
        "als"
    )
)


# heavy_drinkers_prop <- purrr::map_dbl(
#     1:nrow(df_input),
#     ~ prop_heavy_drinkers_row(df_input[.x, , drop = FALSE])
# )

# smoking_prop <- purrr::map_dbl(
#     1:nrow(df_input),
#     ~ calc_prop_alters_smoke(df_input[.x, , drop = FALSE])
# )

# no_exercise_prop <- purrr::map_dbl(
#     1:nrow(df_input),
#     ~ calc_prop_alters_exercise(df_input[.x, , drop = FALSE])
# )

# bad_diet_prop <- purrr::map_dbl(
#     1:nrow(df_input),
#     ~ calc_prop_alters_diet(df_input[.x, , drop = FALSE])
# )

# health_prob_prop <- purrr::map_dbl(
#     1:nrow(df_input),
#     ~ alter_health_problems_row(df_input[.x, , drop = FALSE])
# )

blau_gender <- purrr::map_dbl(
    1:nrow(df_input),
    ~ calc_blau_alter_heterophily(
        df_input[.x, , drop = FALSE],
        "sex",
        gender_map,
        FALSE
    )
)

blau_educ <- purrr::map_dbl(
    1:nrow(df_input),
    ~ calc_blau_alter_heterophily(
        df_input[.x, , drop = FALSE],
        attribute = "educ",
        educ_map,
        FALSE
    )
)

blau_distance <- purrr::map_dbl(
    1:nrow(df_input),
    ~ calc_blau_alter_heterophily(
        df_input[.x, , drop = FALSE],
        attribute = "dist",
        dist_map,
        FALSE
    )
)

blau_length <- purrr::map_dbl(
    1:nrow(df_input),
    ~ calc_blau_alter_heterophily(
        df_input[.x, , drop = FALSE],
        attribute = "length",
        len_map,
        FALSE
    )
)

blau_speak <- purrr::map_dbl(
    1:nrow(df_input),
    ~ calc_blau_alter_heterophily(
        df_input[.x, , drop = FALSE],
        attribute = "speak",
        speak_map,
        FALSE
    )
)


# If you have taken away variables or changed variables, you need to comment out
# variables not included in your dataset in the code section below.
df_clean = tibble::tibble(
    record_id = record_id,
    age = age,
    sex = sex,
    race1 = race1,
    race2 = race2,
    education = education,
    zip = zip,
    # employment = employment,
    # occupation = occupation,
    # income = income,
    # married = married,
    # live_alone = live_alone,
    # household_number = household_number,
    # ego_alcohol = alcohol,
    # ego_smoke = smoke,
    # ego_exercise = exercise,
    # ego_healty_diet = diet,
    # health_problems1 = health_prob1,
    # health_problems2 = health_prob2,
    # health_problems3 = health_prob3,
    # health_problems4 = health_prob4,
    network_size = network_size,
    density = density,
    constraint = constraint,
    effsize = effsize,
    max_degree = max_degree,
    mean_degree = mean_degree,
    kin_prop = kin_prop,
    age_sd = age_sd,
    IQV_sex = iqv_sex,
    # IQV_race = iqv_race,
    IQV_educ = iqv_educ,
    weak_freq_prop = weak_freq_prop,
    weak_dur_prop = weak_dur_prop,
    far_dist_prop = far_dist_prop,
    heavy_drinkers_prop = heavy_drinkers_prop,
    smoking_prop = smoking_prop,
    # no_exercise_prop = no_exercise_prop,
    bad_diet_prop = bad_diet_prop,
    # health_prob_prop = health_prob_prop,
    blau_gender = blau_gender,
    blau_educ = blau_educ,
    blau_distance = blau_distance,
    blau_length = blau_length,
    blau_speak = blau_speak
)

df_clean$zip <- as.character(df_clean$zip) #double check zip as string

######################### FIGURES ##############################
df_relat_names <- names_to_relat(df_input, relat_map)
list_tidygras <- organize_list_tidygraphs(df_relat_names)
list_network_plots_labels <- mapply(
    plot_single_network_node_labels,
    tidygra = list_tidygras,
    ego_name = record_id,
    SIMPLIFY = FALSE
)


#This section creates a montage of network graphs in a grid orientation from top
#  left to bottom right. These graphs are purposefully stripped of alter names
#  and simplified for a better visualization at smaller sizes. Output is a PDF by default.

#Default sizing/scaling for this section is set for the fake datasets of n = 6. You may need to
#  make adjustments for larger/smaller datasets and differently sized outputs. This section
#  only needs data import and functions in the "Establishing Network Graph Functions" section to work.

#To customize the montage easily, modify these three variables to create a desired output.
#  Note that at smaller output sizes and/or larger datasets you may need to reduce the
#  graph_scale, otherwise the edges/nodes will be too large for their respective graphs.


#Width in inches of the output PDF
output_width <- 7.5
#Height in inches of the output PDF
output_height <- 5
#Change the size of the edges/nodes for graphs within montage
graph_scale <- 4
## for ~100 networks, we recommend 8.5, 11, and 2

#Creating list of graphs
list_tidygras <- organize_list_tidygraphs(df_input)

#Gives how many graphs the montage will contain
graph_count <- length(list_tidygras)

#Area of the output PDF
workspace_area <- output_width * output_height
#Estimated area of each graph
est_graph_size <- sqrt(workspace_area / graph_count)
#Estimated number of columns to most evenly fill grid
est_column_count <- ceiling(output_width / est_graph_size)
#Changing scale of graphs by how large graphs will be
est_graph_scale <- graph_scale / est_graph_size

#Constructing graphs, note that this is where our edge/node size scaling actually occurs
#  and this section determines the ratio of edge size to node size.
list_network_plots <- lapply(
    list_tidygras,
    plot_single_network,
    edge_size = 0.25 * est_graph_scale,
    node_size = 1 * est_graph_scale
)
#Constructing grid of plots, note this is where we set how many columns the grid will contain,
#  therefore setting the size and aspect ratio of each plot.
wrap_list_network_plots <- wrap_plots(
    plotlist = list_network_plots,
    ncol = est_column_count
)


### Writing files
write.csv(
    df_clean,
    paste(
        out_dir,
        "Clean_Data",
        ifelse(out_suffix != "", paste("_", out_suffix, sep = ""), ""),
        ".csv",
        sep = ""
    ),
    row.names = FALSE
)

#Creating output PDF. Note that if you'd like to change the output to something other than a PDF,
#  change the extention within the output file name and uncomment the last two lines of code.
ggsave(
    paste(
        out_dir,
        "Social_Network_Grid",
        ifelse(out_suffix != "", paste("_", out_suffix, sep = ""), ""),
        ".pdf",
        sep = ""
    ),
    wrap_list_network_plots,
    width = output_width,
    height = output_height,
    #This code is only needed in case the user wants to change the output to a raster output (not PDF)
    # units = "in",
    # dpi = 300
)

pdf(
    paste(
        out_dir,
        "Single_Networks",
        ifelse(out_suffix != "", paste("_", out_suffix, sep = ""), ""),
        ".pdf",
        sep = ""
    ),
    width = 10,
    height = 10
)
# each print() automatically starts a new page when mfrow = c(1,1)
par(mfrow = c(1, 1))

for (plt in list_network_plots_labels) {
    print(plt)
}

dev.off()

### Employment
#### TODO: Update this so that it works for employment___X
# df_input$employment <- factor(df_input$employment,
#                               levels = c(1, 2, 3, 4, 5, 6, 7, 0),
#                               labels = c("Employed for wages", "Self-employed",
#                                          "Out of work and looking for work",
#                                          "Student", "Retired",
#                                          "Unable to work", "Prefer not to answer",
#                                          "Out of work but not currently looking for work"))
# employment = df_input$employment
#

### Occupation
# df_input$occupation <- factor(df_input$occupation,
#                               levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "77"),
#                               labels = c("Executive, manager",
#                                          "Sales or clerical worker",
#                                          "Mechanic, electrician, skilled worker",
#                                          "Machine operator, inspector, bus/cab driver",
#                                          "Service worker",
#                                          "Professional", "Business owner",
#                                          "Laborer, unskilled worker", "Farming",
#                                          "Military", "Other"))
# occupation = df_input$occupation

### Income
# df_input$income <- factor(df_input$income,
#                           levels = c("1", "2", "3", "4", "5"),
#                           labels = c("less than $5,000", "$5,000 to $49,000",
#                                      "$50,000 to $169,000", "$170,000 to $490,000",
#                                      "more than $500,000"))
# income = df_input$income

### Married
# df_input$married <- factor(df_input$married,
#                            levels = c(0, 1),
#                            labels = c("Not married", "Married"))
# married = df_input$married

### Live Alone
# df_input$live_alone <- factor(df_input$live_alone, levels = c(0, 1),
#                               labels = c("No", "Yes"))
# live_alone = df_input$live_alone

### Household Number
# household_number = df_input$household_number
#

### Alcohol
# df_input$alcohol <- factor(df_input$alcohol,
#                            levels = c(0, 1, 9),
#                            labels = c("No", "Yes", "I do not drink heavily"))
# alcohol = df_input$alcohol

### Exercise
# df_input$exercise <- factor(df_input$exercise,
#                             levels = c(0, 1),
#                             labels = c("No", "Yes"))
# exercise = df_input$exercise

### Diet
# df_input$diet <- factor(df_input$diet,
#                         levels = c(0, 1),
#                         labels = c("No", "Yes"))
# diet = df_input$diet

### Health Problems
# health_prob1 = purrr::map_chr(
#     1:nrow(df_input),
#     ~ extract_health_problems(df_input[.x, , drop = FALSE], health_numeric = 1)
# ) %>%
#     factor(c(
#         "General",
#         "Pain",
#         "Cognitive_MentalHealth",
#         "Cardiac",
#         "NoProblems"
#     ))
# health_prob1

# health_prob2 = purrr::map_chr(
#     1:nrow(df_input),
#     ~ extract_health_problems(df_input[.x, , drop = FALSE], health_numeric = 2)
# ) %>%
#     factor(c(
#         "General",
#         "Pain",
#         "Cognitive_MentalHealth",
#         "Cardiac",
#         "NoProblems"
#     ))
# health_prob2

# health_prob3 = purrr::map_chr(
#     1:nrow(df_input),
#     ~ extract_health_problems(df_input[.x, , drop = FALSE], health_numeric = 3)
# ) %>%
#     factor(c(
#         "General",
#         "Pain",
#         "Cognitive_MentalHealth",
#         "Cardiac",
#         "NoProblems"
#     ))
# health_prob3

# health_prob4 = purrr::map_chr(
#     1:nrow(df_input),
#     ~ extract_health_problems(df_input[.x, , drop = FALSE], health_numeric = 4)
# ) %>%
#     factor(c(
#         "General",
#         "Pain",
#         "Cognitive_MentalHealth",
#         "Cardiac",
#         "NoProblems"
#     ))
# health_prob4

############################# Trouble-shooting ################################

# - If you have difficulties, run line-by-line to find the error

# - Look at your raw data and clean or remove rows with missing data. This code
# does not clean your data. It's processes data assuming the survey was filled
# out correctly.

# - If still not working, then email adhand@bwh.harvard.edu.
