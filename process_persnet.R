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

#Clearing Working Directory
rm(list = ls())

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

#Automatically set working directory, requires rstudioapi package. If not working set WD manually
#to do so manually go to Session -> Set Working Directory -> To Source File location
#copy paste result from console over following commented line of code:
#Example:
# setwd("~/Desktop/network")
setwd(dirname(getActiveDocumentContext()$path))

#Importing data
df_input <- read.csv("fake_data.csv",
                     stringsAsFactors = FALSE,
                     colClasses = c(zip = "character")) #this reads in zip retains 0s

################### Establishing Network Graph Functions ######################

remove_ego_from_igraph <- function(tg_graph) {
  # # # # # # # #
  # Function: Removes the ego node from a tidygraph object, returning  
  #           a graph without the ego.
  # Inputs:  
  #   tg_graph = A tidygraph object representing a personal network
  # Outputs:  
  #   A tidygraph object without the ego node
  # # # # # # # #
  
  tg_graph_egoless <- tg_graph %>%
    tidygraph::activate(nodes) %>%
    dplyr::filter(name != "ego")
  
  return(tg_graph_egoless)
}

organize_row_to_tidygraph <- function(df_row_input) {
  # # # # # # # #
  # Function: Converts a single row of a personal network data frame into  
  #           a tidygraph object, representing network ties
  # Inputs: df_row_input = A single row of a personal network data frame
  # Outputs: A tidygraph object representing the network structure
  # # # # # # # #
  
  # Extract personal network (psn) tie values from the input row
  psn_edges_variables <- dplyr::select(df_row_input, tie1:a_tie105)
  psn_edges_values <- as.integer(psn_edges_variables)  # convert to integer
  
  # Initialize an adjacency matrix for the ego network
  psn_mat <- matrix(NA, 16, 16)  # create a 16x16 matrix (max 15 alters + ego)
  diag(psn_mat) <- 0  # Set diagonal to 0 (no self-loops)
  
  # Populate the lower triangular part of the adjacency matrix with tie values
  psn_mat[lower.tri(psn_mat)] <- psn_edges_values
  
  # Make the matrix symmetric by copying the lower triangle to the upper triangle
  psn_mat <- t(psn_mat)
  psn_mat[lower.tri(psn_mat)] <- psn_edges_values
  
  # Extract column names for alters from the input row
  psn_name_cols <- dplyr::select(df_row_input, c(
    name1, name2, name3, name4, name5,
    name6, name7, name8, name9, name10,
    name11, name12, name13, name14, name15
  ))
  
  # Determine which alters to keep based on a "keep" column (binary indicator)
  psn_keep_name_cols <- dplyr::select(df_row_input, name_1:name_15)
  psn_names_to_keep <- c("ego", ifelse(psn_keep_name_cols == 1, psn_name_cols, NA))
  
  # Flatten the names to a vector, ignoring names marked as NA
  psn_names_for_mat <- as.vector(unlist(psn_names_to_keep, use.names = FALSE))
  colnames(psn_mat) <- rownames(psn_mat) <- psn_names_for_mat  # Set matrix row/column names
  
  # Remove rows and columns corresponding to NAs in the matrix
  psn_mat <- psn_mat[!is.na(rownames(psn_mat)), !is.na(colnames(psn_mat))]
  
  # Identify edges (non-zero entries in the matrix)
  psn_edges <- which(psn_mat != 0, arr.ind = TRUE)
  
  # Check if there are edges (non-isolate)
  if (!is.null(psn_edges) && is.matrix(psn_edges) && nrow(psn_edges) > 0) {
    psn_el <- data.frame(
      from = rownames(psn_mat)[psn_edges[, 1]],
      to = colnames(psn_mat)[psn_edges[, 2]],
      weight = psn_mat[psn_edges]
    )
    
    # Ensure the edge list is undirected (from < to condition)
    psn_el <- psn_el[psn_el$from < psn_el$to, ]
    
    # convert to tidygraph object
    tgra <- tidygraph::as_tbl_graph(psn_el, directed = FALSE) %>%
      dplyr::mutate(record_id = df_row_input$record_id,  # Add record ID
                    node_id = dplyr::row_number())  # Create a unique node ID
    
    return(tgra)
  } else {
    # if no edges exist, create an isolate 
    tgra <- tidygraph::tbl_graph(
      nodes = tibble::tibble(name = c("ego"))) %>%
      tidygraph::activate(nodes) %>%
      dplyr::mutate(record_id = df_row_input$record_id,
                    node_id = dplyr::row_number())
    
    return(tgra)
  }
}

organize_list_tidygraphs <- function(persnet_df) {
  # # # # # # # #
  # Function: Converts a personal network data frame into a list of tidygraph 
  #           objects, with each row representing an individual network.
  # Inputs: persnet_df = A personal network data frame
  # Outputs: A list of tidygraph objects
  # # # # # # # #
  
  # Split the data frame into a list of individual rows
  df_as_list <- persnet_df %>%
    dplyr::mutate(index = 1:dplyr::n()) %>%  
    dplyr::group_split(index)  # Split by the index
  
  # Convert each row into a tidygraph object using row_to_tidygraph()
  tidygraph_list <- lapply(df_as_list, organize_row_to_tidygraph)
  
  return(tidygraph_list)
}

######################### Ego Based Calculations ##############################

############################# Record ID #######################################

record_id = df_input$record_id
record_id

################################ Age ##########################################

age = df_input$age
age

################################ Sex ##########################################

sex = factor(df_input$sex, levels = c(0, 1, 2),
             labels = c("Female", "Male", "other"))
sex

################################ Race #########################################

extract_race_identity <- function(persnet_row, race_numeric) {
  # race_numeric is which of 1st, 2nd, 3rd, or 4th race 
  # race designation labels (with capitalization and white spaces)
  race_labels <- c("Black", "White", "American Indian", "Asian",
                   "Hawaiian and Pacific Islander", "Other", NA)
  #these specific race columns can be edited depending on column names
  race_columns <- c("race___1", "race___2", "race___3",
                    "race___4", "race___5", "race___77", "race___88")
  # identify which race columns are marked as 1
  selected_races <- race_labels[which(persnet_row[race_columns] == 1)]
  #select which item in list of race columns == 1
  return(selected_races[race_numeric])
}

#purrr applies the same function to a list, in this case, each row of df_input
#the blank space between .x, , is selecting all columns of df_input
#while drop= FALSE forces R to return not a vector but a dataframe object
race1 = purrr::map_chr(1:nrow(df_input),
                       ~ extract_race_identity(df_input[.x, , drop = FALSE], race_numeric = 1)) %>%
  factor(c("Black", "White", "American Indian", "Asian",
           "Hawaiian and Pacific Islander", "Other", NA))
race1

race2 = purrr::map_chr(1:nrow(df_input),
                       ~ extract_race_identity(df_input[.x, , drop = FALSE], race_numeric = 2)) %>%
  factor(c("Black", "White", "American Indian", "Asian",
           "Hawaiian and Pacific Islander", "Other", NA))
race2

race3 = purrr::map_chr(1:nrow(df_input),
                       ~ extract_race_identity(df_input[.x, , drop = FALSE], race_numeric = 3)) %>%
  factor(c("Black", "White", "American Indian", "Asian",
           "Hawaiian and Pacific Islander", "Other", NA))
race3

race4 = purrr::map_chr(1:nrow(df_input),
                       ~ extract_race_identity(df_input[.x, , drop = FALSE], race_numeric = 4)) %>%
  factor(c("Black", "White", "American Indian", "Asian",
           "Hawaiian and Pacific Islander", "Other", NA))
race4

################################## Education ##################################

df_input$edu <- factor(df_input$edu, 
                       levels = c(1, 2, 3, 4, 5, 6, 88),
                       labels = c("Some High School", "High School Grad",
                                  "Some College", "Associate Degree", 
                                  "Bachelor's Degree", "Graduate Degree", 
                                  "No Answer"))
education <- df_input$edu
education

################################## Zip Code ###################################

zip = as.character(df_input$zip)
zip

################################# Employment ##################################

df_input$employment <- factor(df_input$employment, 
                              levels = c(1, 2, 3, 4, 5, 6, 7, 0),
                              labels = c("Employed for wages", "Self-employed",
                                         "Out of work and looking for work", 
                                         "Student", "Retired", 
                                         "Unable to work", "Prefer not to answer",
                                         "Out of work but not currently looking for work"))
employment = df_input$employment
employment

############################## Occupation #####################################

df_input$occupation <- factor(df_input$occupation, 
                              levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "77"),
                              labels = c("Executive, manager", 
                                         "Sales or clerical worker", 
                                         "Mechanic, electrician, skilled worker", 
                                         "Machine operator, inspector, bus/cab driver", 
                                         "Service worker", 
                                         "Professional", "Business owner", 
                                         "Laborer, unskilled worker", "Farming", 
                                         "Military", "Other"))
occupation = df_input$occupation
occupation

############################## Income #########################################

df_input$income <- factor(df_input$income, 
                          levels = c("1", "2", "3", "4", "5"),
                          labels = c("less than $5,000", "$5,000 to $49,000", 
                                     "$50,000 to $169,000", "$170,000 to $490,000", 
                                     "more than $500,000"))
income = df_input$income
income

########################### Marital Status ####################################

df_input$married <- factor(df_input$married, 
                           levels = c(0, 1),
                           labels = c("Not married", "Married"))
married = df_input$married
married

############################ Live Alone #######################################

df_input$live_alone <- factor(df_input$live_alone, levels = c(0, 1), 
                              labels = c("No", "Yes"))
live_alone = df_input$live_alone
live_alone

###################### Number of Household Members ############################

household_number = df_input$household_number
household_number

################################ Alcohol ######################################

df_input$alcohol <- factor(df_input$alcohol, 
                           levels = c(0, 1, 9),
                           labels = c("No", "Yes", "I do not drink heavily"))
alcohol = df_input$alcohol
alcohol

################################ Smoking ######################################

df_input$smoke <- factor(df_input$smoke, 
                         levels = c(0, 1, 9),
                         labels = c("No", "Yes", "I do not smoke"))
smoke = df_input$smoke
smoke

############################### Exercise #####################################

df_input$exercise <- factor(df_input$exercise, 
                            levels = c(0, 1),
                            labels = c("No", "Yes"))
exercise = df_input$exercise
exercise

############################### Diet ##########################################

df_input$diet <- factor(df_input$diet, 
                        levels = c(0, 1),
                        labels = c("No", "Yes"))
diet = df_input$diet
diet

############################## Health Problems ################################

extract_health_problems <- function(persnet_row, health_numeric) {
  # # # # # # # #
  # Function: Extracts a specific health problem label from a single row of a personal
  #           network dataframe. The function checks a set of predefined health problem
  #           columns and returns the health label corresponding to the nth (health_numeric)
  #           problem that is flagged (i.e., set to 1).
  #
  # Inputs:
  #   persnet_row  - A single row of a personal network dataframe. This row should include
  #                  a set of columns (e.g., "health___1", "health___2", etc.) representing
  #                  different health problems, where a value of 1 indicates that the health
  #                  problem is present.
  #   health_numeric - A numeric value (1, 2, 3, …) indicating which flagged health problem
  #                  to extract. For example, 1 returns the first flagged health problem, 2
  #                  returns the second, and so on.
  #
  # Process:
  #   1. A vector 'health_labels' is defined containing the labels for each health problem.
  #   2. A corresponding vector 'health_columns' is set with the names of the dataframe columns.
  #   3. The function checks the specified health columns in the given row. For any column with
  #      a value of 1, it selects the corresponding label from 'health_labels'.
  #   4. It then checks if the number of flagged health problems is at least as large as 'health_numeric'.
  #
  # Outputs:
  #   - If there are enough flagged health problems, the function returns the label for the
  #     health problem at the specified position (e.g., the 1st, 2nd, or 3rd problem).
  #   - If there are fewer flagged health problems than requested, the function returns NA.
  # # # # # # # #
  health_labels <- c("General","Pain","Cognitive_MentalHealth",
                     "Cardiac","NoProblems")
  #tweak name of health columns here if need be
  health_columns <- c("health___1", "health___2", "health___3",
                      "health___4", "health___0")
  selected_health <- health_labels[which(persnet_row[health_columns] == 1)]
  # Return the requested health problem number (1st, 2nd, 3rd, 4th) or NA
  if (length(selected_health) >= health_numeric) {
    return(selected_health[health_numeric])
  } else {
    return(NA)
  }
}
health_prob1 = purrr::map_chr(1:nrow(df_input),
                              ~ extract_health_problems(df_input[.x, , drop = FALSE], health_numeric = 1)) %>%
  factor(c("General","Pain","Cognitive_MentalHealth", "Cardiac","NoProblems"))
health_prob1

health_prob2 = purrr::map_chr(1:nrow(df_input),
                              ~ extract_health_problems(df_input[.x, , drop = FALSE], health_numeric = 2)) %>%
  factor(c("General","Pain","Cognitive_MentalHealth", "Cardiac","NoProblems"))
health_prob2

health_prob3 = purrr::map_chr(1:nrow(df_input),
                              ~ extract_health_problems(df_input[.x, , drop = FALSE], health_numeric = 3)) %>%
  factor(c("General","Pain","Cognitive_MentalHealth", "Cardiac","NoProblems"))
health_prob3

health_prob4 = purrr::map_chr(1:nrow(df_input),
                              ~ extract_health_problems(df_input[.x, , drop = FALSE], health_numeric = 4)) %>%
  factor(c("General","Pain","Cognitive_MentalHealth", "Cardiac","NoProblems"))
health_prob4

######################## Network Structure Measures ###########################

###################### Total Network Size (Unique Alters) #####################

calc_total_alters_row <- function(persnet_row) {
  # # # # # # # #
  # Function: Computes the total number of unique alters in a given  
  #           personal network row, combining named alters and additional  
  #           names listed in "more_names" columns.
  # Inputs:  
  #   persnet_row = A single row of a personal network data frame
  # Outputs:  
  #   Total count of unique alters in the network
  # # # # # # # #
  
  # check if inout is valid
  if (is.null(persnet_row) || !"data.frame" %in% class(persnet_row)) {
    stop("The input must be a single row of a personal network dataframe.")
  }
  
  # Convert the row into a tidygraph object
  tc_tidygraph <- organize_row_to_tidygraph(persnet_row)
  
  # Extract unique alter names from the tidygraph, excluding the ego node
  unique_names_tidygra <- tc_tidygraph %>%
    tidygraph::activate(nodes) %>%       
    dplyr::as_tibble() %>%           
    dplyr::filter(name != "ego") %>% 
    dplyr::pull(name) %>%           
    unique()
  
  # Extract and clean additional alter names from "more_names" columns
  all_more_names <- paste(persnet_row$more_names_1, 
                          persnet_row$more_names_2, 
                          persnet_row$more_names_3, 
                          sep = ", ") %>%
    stringr::str_split(",\\s*") %>%  # split by commas
    unlist() %>%            # conver to vector
    .[. != ""] %>%          # semove empty strings     
    unique() 
  
  # Return the count of unique alters from both sources
  return(length(union(unique_names_tidygra, all_more_names)))
}

calc_total_alters_df <- function(persnet_df) {
  #applies total alters row to each dataframe
  tryCatch({
    #apply to each row the function of calculate total alters 
    vector_count_total_alters <- apply(persnet_df, 1, function(row) 
      calc_total_alters_row(persnet_df %>% dplyr::filter(record_id == row["record_id"])))
    return(vector_count_total_alters)  
  }, error = function(e) { 
    # If an error occurs, stop with a message indicating that the input is not a valid persnet dataframe
    warning(" error: there was an error calculating number of unqiue alters, returning NA")  
    return(NA)
  })
}

network_size = calc_total_alters_df(df_input)
network_size

########################### Network Density ###################################

calc_egoless_density <- function(tg_graph) {
  # # # # # # # #
  # Function: Computes the density of a personal network graph after  
  #           removing the ego node.
  # Inputs:  
  #   tg_graph = A tidygraph object representing a personal network
  # Outputs:  
  #   Density of the network without the ego node (rounded to 2 decimal places)
  # # # # # # # #
  
  return(round(igraph::edge_density(remove_ego_from_igraph(tg_graph)), 2))
}

calc_egoless_density_df <- function(persnet_df) {
  tryCatch({
    # Organize a list of tidygraphs from the persnet dataframe
    gra_list <- organize_list_tidygraphs(persnet_df)
    # Calculate the density of each network after removing the ego node
    vector_egoless_density <- sapply(gra_list, calc_egoless_density)
    return(vector_egoless_density)  
  }, error = function(e) { 
    # If an error occurs, stop with a message indicating that the input is not a valid persnet dataframe
    stop("Error: not a valid persnet dataframe")  
  })
}

density = calc_egoless_density_df(df_input)
density

######################### Network Constraint ##################################

calc_node_constraint <- function(tidygra, node_index = NULL) {
  # # # # # # # #
  # Function: Computes Burt's Constraint measure for a specified node  
  #           in a tidygraph object, quantifying structural holes.
  # Inputs:  
  #   tidygra    = A tidygraph object representing a personal network
  #   node_index = (Optional) The node for which constraint is calculated 
  # (default = "ego")
  # Outputs:  
  #   Burt's Constraint value for the specified node
  # # # # # # # #
  
  # test if the input is a valid tidygraph object
  if (is.null(tidygra) || !"tbl_graph" %in% class(tidygra)) {
    stop("The input graph must be a tidygraph graph object.")
  }
  
  # check if the graph is an isolate (no edges), return NA otherwise
  if (dim(igraph::as_adjacency_matrix(tidygra))[1] < 2) return(NA)
  
  # check if edge attribute weight is missing
  edge_attrs <- tidygra %>% tidygraph::activate(edges) %>% tibble::as_tibble()
  if (!"weight" %in% colnames(edge_attrs)) {
    warning("Network has no weights; assigning all tie weights as 1.")
    tidygra <- tidygra %>%
      tidygraph::activate(edges) %>%
      dplyr::mutate(weight = 1)
  }
  
  # Convert to adjacency matrix (undirected) with weight
  adj_matrix <- igraph::as_adjacency_matrix(tidygra, attr = "weight", sparse = FALSE)
  node_names <- igraph::V(tidygra)$name
  
  # identify ego node index
  if (is.null(node_index)) {
    if (!"ego" %in% node_names) {
      stop("No node labeled 'ego' exists in the graph.")
    }
    ego_index <- which(node_names == "ego")
  } else {
    if (!node_index %in% node_names) {
      stop("The specified node does not exist.")
    }
    ego_index <- which(node_names == node_index)
  }
  
  num_nodes <- nrow(adj_matrix)
  
  # total tie strength for each node (undirected)
  sum_tie_strengths <- rowSums(adj_matrix)
  
  # compute proportional tie strengths for undirected networks
  matrix_proportions <- matrix(0, nrow = num_nodes, ncol = num_nodes)
  for (i in 1:num_nodes) {
    if (sum_tie_strengths[i] == 0) next
    for (j in setdiff(1:num_nodes, i)) {
      matrix_proportions[i, j] <- adj_matrix[i, j] / sum_tie_strengths[i]
    }
  }
  
  # Calculate Burt's constraint for ego node
  constraint <- 0
  for (j in setdiff(1:num_nodes, ego_index)) {
    redundancy <- sum(matrix_proportions[ego_index, -c(ego_index, j)] * matrix_proportions[-c(ego_index, j), j])
    constraint <- constraint + (matrix_proportions[ego_index, j] + redundancy)^2
  }
  
  return(round(constraint, 2))
}

gra_list <- organize_list_tidygraphs(df_input)
constraint = sapply(gra_list, calc_node_constraint)
constraint

######################### Effective Network Size ##############################

calc_node_ens <- function(tidygra, node_index = NULL) {
  # # # # # # # #
  # Function: Computes the Effective Network Size (ENS) for a specified node 
  #           in a tidygraph object, considering tie strengths.
  # Inputs:  
  #   tidygra    = A tidygraph object representing a personal network
  #   node_index = (Optional) The node for which ENS is calculated (default = "ego")
  # Outputs:  
  #   Effective Network Size (ENS) value for the specified node
  # # # # # # # #
  
  # Test if the input is valid
  if (is.null(tidygra) || !"tbl_graph" %in% class(tidygra)) {
    stop("The input graph must be a tidygraph graph object.")
  }
  
  # Check if the network is an isolate and return NA if so
  if (dim(igraph::as_adjacency_matrix(tidygra))[1] < 2) return(NA)
  
  # Check if the edge attribute 'weight' is missing
  edge_attrs <- tidygra %>% tidygraph::activate(edges) %>% tibble::as_tibble()
  if (!"weight" %in% colnames(edge_attrs)) {
    warning("Warning: Network has no weights. Effective Network Size typically includes weight. Calculating with all tie weights set to 1.")
    tidygra <- tidygra %>%
      tidygraph::activate(edges) %>%
      dplyr::mutate(weight = 1)  # Add default weight of 1
  }
  
  # convert to adjacency matrix
  adj_matrix <- igraph::as_adjacency_matrix(tidygra, attr = "weight", sparse = FALSE)
  node_names <- igraph::V(tidygra)$name  # Extract node names from graph
  
  # determine the node index for calculation
  if (is.null(node_index)) {
    if (!"ego" %in% node_names) {
      stop("No node labeled 'ego' exists in the graph.")
    }
    ego_index <- which(node_names == "ego")
  } else {
    # validate that the specified node exists
    if (!node_index %in% node_names) {
      stop("The specified node does not exist.")
    }
    ego_index <- which(node_names == node_index)
  }
  
  # number of nodes in the network
  num_nodes <- dim(adj_matrix)[1]
  
  # compute total tie strength for ego
  total_tie_strength <- sum(adj_matrix[ego_index, ])  
  
  # initialize effective size result
  effective_size_result <- 0  
  
  # loop through nodes excluding the specified ego node
  for (current_node in setdiff(1:num_nodes, ego_index)) {  
    max_tie_strength <- max(adj_matrix[current_node, ], na.rm = TRUE)  # Max tie strength for current node
    if (max_tie_strength == 0) next  # Skip if no ties exist
    
    redundancy_sum <- 0  # Initialize redundancy sum for current node
    for (other_neighbor in setdiff(1:num_nodes, c(ego_index, current_node))) {
      # proportion of ego's total ties to other_neighbor
      ego_proportion <- adj_matrix[ego_index, other_neighbor] / total_tie_strength
      
      # proportion of current_node's strongest tie to other_neighbor
      node_redundancy <- adj_matrix[current_node, other_neighbor] / max_tie_strength
      
      redundancy_sum <- redundancy_sum + (ego_proportion * node_redundancy)
    }
    
    # update effective size result
    effective_size_result <- effective_size_result + (1 - redundancy_sum)
  }
  return(round(effective_size_result,2)) 
}
gra_list <- organize_list_tidygraphs(df_input)
effsize = sapply(gra_list, calc_node_ens)
effsize

############################ Mean Degree ######################################

calc_egoless_mean_degree <- function(tg_graph) {
  # # # # # # # #
  # Function: Computes the mean degree of nodes in a personal  
  #           network after removing the ego node.
  # Inputs:  
  #   tg_graph = A tidygraph object representing a personal network
  # Outputs:  
  #   Mean degree of nodes in the network (excluding ego)
  # # # # # # # #
  
  tryCatch({
    # Remove ego
    egoless_graph <- remove_ego_from_igraph(tg_graph)
    
    # Check if the graph is empty (no nodes)
    if (igraph::vcount(egoless_graph) == 0) {
      return(0)
    }
    
    # Else calculate mean degree
    return(round(mean(igraph::degree(egoless_graph)), 2))
  }, error = function(e) {
    # If an error occurs, return NA
    return(NA)
  })
}

gra_list <- organize_list_tidygraphs(df_input)
mean_degree = sapply(gra_list, calc_egoless_mean_degree)
mean_degree

########################## Max Degree #########################################

calc_egoless_max_degree <- function(tg_graph) {
  # # # # # # # #
  # Function: Computes the maximum degree of nodes in a personal  
  #           network after removing the ego node.
  # Inputs:  
  #   tg_graph = A tidygraph object representing a personal network
  # Outputs:  
  #   Maximum degree of nodes in the network (excluding ego)
  # # # # # # # #
  
  tryCatch({
    # Remove ego
    
    egoless_graph <- remove_ego_from_igraph(tg_graph)
    
    # Check if the graph is empty, if so return 0
    if (igraph::vcount(egoless_graph) == 0) {
      return(0)
    }
    
    # Else calculate max degree
    return(max(igraph::degree(egoless_graph)))
  }, error = function(e) {
    # If an error occurs, return NA
    return(NA)
  })
}
gra_list <- organize_list_tidygraphs(df_input)
max_degree = sapply(gra_list, calc_egoless_max_degree)
max_degree

####################### Network Composition Measures ##########################

########################### Relationship Type #################################

calc_prop_alters_relationship <- function(persnet_row, relationship_type) {
  # # # # # # # #
  # Function: Computes the proportion of alters with a specified relationship  
  #           type in a personal network row.
  # Inputs:  
  #   persnet_row       = A single row of a personal network data frame
  #   relationship_type = One of 'spouse', 'family', 'friend', 'advice',  
  #                      'coworker', or 'other'
  # Outputs:  
  #   Proportion of alters with the specified relationship type
  # # # # # # # #
  
  # Relationship types mapped to numeric codes
  relationship_map <- c(
    spouse = 1,
    family = 2,
    friend = 3,
    advice = 4,
    coworker = 5,
    other = 77
  )
  
  # Validate the relationship_type input
  if (!(relationship_type %in% names(relationship_map))) {
    stop("Error: Choose one of 'spouse', 'family', 'friend', 'advice', 'coworker', or 'other'.")
  }
  
  # Convert row to a tidygraph object and check if the network is an isolate
  tg_graph <- organize_row_to_tidygraph(persnet_row)
  if (igraph::vcount(tg_graph) == 1) {
    return(NA)
  }
  
  # Identify relationship-related columns
  relationship_columns <- grep("^name\\d+relat___\\d+$", colnames(persnet_row), value = TRUE)
  
  # Handle missing relationship columns
  if (length(relationship_columns) == 0) {
    return("Error: Relationship columns not present in this version of persnet.")
  }
  
  # Filter for the relevant relationship type
  relevant_columns <- grep(paste0("___", relationship_map[[relationship_type]]), relationship_columns, value = TRUE)
  if (length(relevant_columns) == 0) {
    return(0)  # If no relevant columns exist, proportion is 0
  }
  
  # Select the relevant relationship type columns
  relationship_selected <- persnet_row %>% dplyr::select(dplyr::all_of(relevant_columns))
  
  # Calculate and return the proportion of alters with the selected relationship type
  return(
    length(which(relationship_selected == 1)) /
      sum(persnet_row %>% dplyr::select(tie1:tie15) != 0, na.rm = TRUE)
  )
}

prop_kin_persnet_row <- function(persnet_row) {
  prop_kin <- calc_prop_alters_relationship(persnet_row, "spouse") +
    calc_prop_alters_relationship(persnet_row, "family")
  return(round(prop_kin,2))
}

#this function applies the same function to each row
kin_prop = purrr::map_dbl(1:nrow(df_input), ~ prop_kin_persnet_row(df_input[.x, , 
                                                                            drop = FALSE]))
kin_prop

########################### Age Standard Deviation ############################

sd_age_alters_row <- function(persnet_row) {
  # Compute the standard deviation of the ages of the alters in the egonetwork
  return(
    round(stats::sd(persnet_row %>% dplyr::select(name1age:name15age) %>% unlist(), 
                    na.rm = TRUE), 2)
  )
}
age_sd = purrr::map_dbl(1:nrow(df_input), 
                        ~ sd_age_alters_row(df_input[.x, , drop = FALSE]))
age_sd

####################### Diversity Calculation Functions #######################

calc_blau_alter_heterophily <- function(persnet_row, attribute = NULL) {
  # # # # # # # #
  # Function: Computes the Blau heterophily index for a specified alter attribute.  
  #           The index measures diversity within a personal network.
  # Inputs:  
  #   persnet_row = A single row of a personal network data frame
  #   attribute   = One of 'gender', 'race', 'educ', 'support' (types of support),  
  #                 'distance' (how far alters live from ego), 'length'  
  #                 (how long alters have known ego), or 'speak'  
  #                 (how often alters speak to ego)
  # Outputs:  
  #   Blau heterophily index for the specified attribute
  # # # # # # # #
  
  valid_attributes <- c(
    "gender", "educ", 
    "distance", "length", "speak",
    "support", "relationships_family", "race"
  ) #
  
  if (is.null(attribute) || is.na(attribute) || !(attribute %in% valid_attributes)) {
    warning("Error: Choose one of 'gender', 'educ', 
    'support' (types of support), 'distance' (how far alters live from ego),
    'length' (how long alters have known ego), or 'speak' (how often alters speak to ego).")
    return(NA)
  } #'relationships_family' (family vs non-family),
  
  # Calculate and return Blau heterophily index for gender; this is meaningless?
  if (attribute == 'gender') {
    gender_cols <- persnet_row %>% dplyr::select(name1sex:name15sex)
    prop_men <- sum(gender_cols, na.rm = TRUE) / sum(persnet_row %>% dplyr::select(tie1:tie15) != 0, na.rm = TRUE)
    prop_women <- length(which(gender_cols == 0)) / sum(persnet_row %>% dplyr::select(tie1:tie15) != 0, na.rm = TRUE)
    blau_index_gender <- 1 - (prop_men^2 + prop_women^2)
    return(round(blau_index_gender,2))                               
  }
  
  # Calculate and return Blau heterophily index for race; commented out as by default
  ## alters can be of multiple races, and Blau hterophily is not defined for overlapping categories
  if (attribute == 'race') {
    prop_black <- calc_prop_alters_race_identity(persnet_row,"black") 
    prop_white <- calc_prop_alters_race_identity(persnet_row,"white")  
    prop_native_american <- calc_prop_alters_race_identity(persnet_row,"native_american")  
    prop_asian <- calc_prop_alters_race_identity(persnet_row,"asian")
    prop_pacific_islanders <- calc_prop_alters_race_identity(persnet_row,"pacific_islander")  
    prop_unknown <- calc_prop_alters_race_identity(persnet_row,"unknown") 
    blau_index_race <- 1 - sum(prop_black^2,
                               prop_white^2,
                               prop_native_american^2,
                               prop_asian^2,
                               prop_pacific_islanders^2,
                               prop_unknown^2
    )
    return(round(blau_index_race,2))                               
  }
  
  # Calculate and return Blau heterophily index for education
  if (attribute == 'educ') {
    prop_only_highschool <- calc_prop_alters_education_level(persnet_row,'only_high_school')  
    prop_some_college <- calc_prop_alters_education_level(persnet_row,'some_college')  
    prop_college_grad <- calc_prop_alters_education_level(persnet_row,'college_grad')  
    prop_dont_know_edu <- calc_prop_alters_education_level(persnet_row,'dont_know')
    blau_index_educ <- 1 - sum(prop_only_highschool^2, prop_some_college^2, prop_college_grad^2, prop_dont_know_edu^2)
    return(round(blau_index_educ,2))                               
  }
  
  # Blau heterophily assumes mutually exclusive categories,
  # but alters can be in multiple relationships and offer multiple types of support.
  
  # Calculate and return Blau heterophily index for distance
  if (attribute == 'distance') {
    blau_index_distance <- 1 - sum(
      calc_prop_alters_distance_away(persnet_row,"same_house")^2,
      calc_prop_alters_distance_away(persnet_row,"1_5miles")^2,
      calc_prop_alters_distance_away(persnet_row,"6_15miles")^2,
      calc_prop_alters_distance_away(persnet_row,"16_50miles")^2,
      calc_prop_alters_distance_away(persnet_row,"more50miles")^2
    )
    return(round(blau_index_distance,2))
  }
  
  # Calculate and return Blau heterophily index for frequency of speaking
  if (attribute == 'speak') {
    blau_index_speak <- 1 - sum(
      calc_prop_alters_freq_speak(persnet_row,'daily')^2,
      calc_prop_alters_freq_speak(persnet_row,'weekly')^2,
      calc_prop_alters_freq_speak(persnet_row,'monthly')^2,
      calc_prop_alters_freq_speak(persnet_row,'more_monthly')^2,
      calc_prop_alters_freq_speak(persnet_row,'dont_know')^2
    )
    return(round(blau_index_speak,2))
  }
  
  # Calculate and return Blau heterophily index for length of knowing ego
  if (attribute == 'length') {
    blau_index_length <- 1 - sum(
      calc_prop_alters_known_length(persnet_row,"less_than_3years")^2,
      calc_prop_alters_known_length(persnet_row,"three_to_6years")^2,
      calc_prop_alters_known_length(persnet_row,"more_than_6years")^2,
      calc_prop_alters_known_length(persnet_row,"unknown")^2
    )
    return(round(blau_index_length,2))
  }
}

calc_attribute_iqv <- function(persnet_row, attribute = NULL) {
  # # # # # # # #
  # Function: Computes the Index of Qualitative Variation (IQV) for a specified
  #           alter attribute. The IQV measures diversity and variation within
  #           a personal network.
  # Inputs:
  #   persnet_row = A single row of a personal network data frame
  #   attribute   = One of 'gender', 'educ', 'relationships' (family vs non-family),
  #                 or 'support' (types of support)
  # Outputs:
  #   IQV value for the specified attribute
  # # # # # # # #
  
  # define valid attributes for IQV calculation
  valid_attributes <- c(
    "race",
    "gender",
    "educ",
    "relationships",
    "support"
  )
  
  # validate the attribute input
  if (is.null(attribute) || is.na(attribute) || !(attribute %in% valid_attributes)) {
    warning("Error: Choose one of gender' or 'educ'")
    return(NA)
  }
  
  # calculate IQV using Blau heterophily index and corresponding normalization factor
  if (attribute == 'race') {
    return(
      round(calc_blau_alter_heterophily(persnet_row, 'race') / (1 - 1 / 5),2)
    )
  }
  if (attribute == 'gender') {
    return(
      round(calc_blau_alter_heterophily(persnet_row, 'gender') / (1 - 1 / 2),2)
    )
  }
  if (attribute == 'educ') {
    return(
      round(calc_blau_alter_heterophily(persnet_row, 'educ') / (1 - 1 / 4),2)
    )
  }
  if (attribute == 'relationships') {
    return(
      round(calc_blau_alter_heterophily(persnet_row, 'relationships') / (1 - 1 / 6),2)
    )
  }
  if (attribute == 'support') {
    return(
      round(calc_blau_alter_heterophily(persnet_row, 'support') / (1 - 1 / 5),2)
    )
  }
}

############################# Diversity of Sex ################################

#this command applies the function "calc_attribute_iqv" to each row of 
##df_input (the .x is the index), drop=FALSE forces R to return a dataframe
#rather than vector
# IQV is essentially a normalized version of Blau's heterogenity index, and hence 
#the IQV function relies on the blau alter heterophily
# Both of these can be tweaked for specific variable names
iqv_sex <- purrr::map_dbl(1:nrow(df_input), ~ calc_attribute_iqv(df_input[.x, , drop = FALSE],
                                                                 attribute = "gender"))
iqv_sex

############################## Race Diversity #################################

calc_prop_alters_race_identity <- function(persnet_row, race_category) {
  # # # # # # # #
  # Function: Computes the proportion of alters who belong to a specified  
  #           racial category.
  # Inputs:  
  #   persnet_row  = A single row of a personal network data frame
  #   race_category = One of 'black', 'white', 'native_american',  
  #                   'asian', 'pacific_islander', or 'unknown'
  # Outputs:  
  #   Proportion of alters within the specified racial category
  # # # # # # # #
  
  # Map race categories to numeric codes
  race_map <- c(
    black = 1,
    white = 2,
    native_american = 3,
    asian = 4,
    pacific_islander = 5,
    unknown = 77
  )
  
  # Validate the race_category input
  if (!(race_category %in% names(race_map))) {
    stop("Error: Choose one of 'black', 'white', 'native_american', 'asian', 
         'pacific_islander', or 'unknown'.")
  }
  
  # Identify race-related columns
  race_cols_string <- grep("^name\\d+race$", colnames(persnet_row), value = TRUE)
  
  # Handle missing columns
  if (length(race_cols_string) == 0) {
    stop("Error: Cannot find columns related to race.")
  }
  
  # Select race columns
  race_cols <- persnet_row %>% dplyr::select(name1race:name15race)
  
  # Calculate and return the proportion of alters within the specified race category
  return(
    length(which(race_cols == race_map[[race_category]])) /
      sum(persnet_row %>% dplyr::select(tie1:tie15) != 0, na.rm = TRUE)
  )
}
iqv_race <- purrr::map_dbl(1:nrow(df_input), ~ calc_attribute_iqv(df_input[.x, , drop = FALSE], 
                                                                  attribute = "race"))
iqv_race

######################### Diversity of Education ##############################

calc_prop_alters_education_level <- function(persnet_row, education_level) {
  # # # # # # # #
  # Function: Computes the proportion of alters who fall into a specified  
  #           education level category.
  # Inputs:  
  #   persnet_row     = A single row of a personal network data frame
  #   education_level = One of 'only_high_school', 'some_college',  
  #                     'college_grad', or 'dont_know'
  # Outputs:  
  #   Proportion of alters within the specified education level
  # # # # # # # #
  
  # Map education levels to numeric codes
  education_map <- list(
    only_high_school = c(1, 2),
    some_college = c(3, 4),
    college_grad = c(5, 6),
    dont_know = c(99)
  )
  
  # Validate the education_level input
  if (!(education_level %in% names(education_map))) {
    stop("Error: Choose one of 'only_high_school', 'some_college', 'college_grad', or 'dont_know'.")
  }
  
  # Identify education-related columns
  edu_cols_string <- grep("^name\\d+educ$", colnames(persnet_row), value = TRUE)
  
  # Handle missing columns
  if (length(edu_cols_string) == 0) {
    stop("Error: Cannot find columns related to education.")
  }
  
  # Select education columns
  edu_cols <- persnet_row %>% dplyr::select(name1educ:name15educ)
  
  # Calculate and return the proportion of alters within the specified education level
  return(
    length(which(edu_cols %in% education_map[[education_level]])) /
      sum(persnet_row %>% dplyr::select(tie1:tie15) != 0, na.rm = TRUE)
  )
} 
iqv_educ <- purrr::map_dbl(1:nrow(df_input), ~ calc_attribute_iqv(df_input[.x, , drop = FALSE], 
                                                                  attribute = "educ"))
iqv_educ

####################### Frequency of Interaction ##############################

calc_prop_alters_freq_speak <- function(persnet_row, frequency) {
  # # # # # # # #
  # Function: Computes the proportion of alters who speak to the ego  
  #           at a specified frequency.
  # Inputs:  
  #   persnet_row = A single row of a personal network data frame
  #   frequency   = One of 'daily', 'weekly', 'monthly', 'more_monthly', or 'dont_know'
  # Outputs:  
  #   Proportion of alters who speak to the ego at the specified frequency
  # # # # # # # #
  
  # Map frequency labels to numeric codes
  frequency_map <- c(
    daily = 1,
    weekly = 2,
    monthly = 3,
    more_monthly = 4,
    dont_know = 99
  )
  
  # Validate the frequency input
  if (!(frequency %in% names(frequency_map))) {
    stop("Error: Choose one of 'daily', 'weekly', 'monthly', 'more_monthly', or 'dont_know'.")
  }
  
  # Identify speaking frequency-related columns
  speak_cols_string <- grep("^name\\d+speak$", colnames(persnet_row), value = TRUE)
  
  # Handle missing columns
  if (length(speak_cols_string) == 0) {
    stop("Error: Cannot find columns related to how frequently alters speak to the patient.")
  }
  
  # Select speaking frequency columns
  speak_cols <- persnet_row %>% dplyr::select(name1speak:name15speak)
  
  # Calculate and return the proportion of alters with the specified frequency
  return(
    length(which(speak_cols == frequency_map[[frequency]])) /
      sum(persnet_row %>% dplyr::select(tie1:tie15) != 0, na.rm = TRUE)
  )
}
weak_freq_prop_row <- function(persnet_row) {
  # Proportion of alters who contact ego monthly or less frequently
  return(
    round(
      calc_prop_alters_freq_speak(persnet_row, "monthly") +
        calc_prop_alters_freq_speak(persnet_row, "more_monthly"),
      2
    )
  )
}
weak_freq_prop <- purrr::map_dbl(1:nrow(df_input), ~ weak_freq_prop_row(df_input[.x, , drop = FALSE]))
weak_freq_prop

########################## Duration of Relationship ###########################

calc_prop_alters_known_length <- function(persnet_row, time_category) {
  # # # # # # # #
  # Function: Computes the proportion of alters who have known the ego  
  #           for a specified length of time.
  # Inputs:  
  #   persnet_row   = A single row of a personal network data frame
  #   time_category = One of 'less_than_3years', 'three_to_6years',  
  #                   'more_than_6years', or 'unknown'
  # Outputs:  
  #   Proportion of alters within the specified time category
  # # # # # # # #
  
  # Map time categories to numeric codes
  time_map <- c(
    less_than_3years = 1,
    three_to_6years = 2,
    more_than_6years = 3,
    unknown = 99
  )
  
  # Validate the time_category input
  if (!(time_category %in% names(time_map))) {
    stop("Error: Choose one of 'less_than_3years', 'three_to_6years', 'more_than_6years', or 'unknown'.")
  }
  
  # Identify length-of-relationship-related columns
  length_cols_string <- grep("^name\\d+length$", colnames(persnet_row), value = TRUE)
  
  # Handle missing columns
  if (length(length_cols_string) == 0) {
    stop("Error: Cannot find columns related to how long alters have known the ego.")
  }
  
  # Select length-of-relationship columns
  length_cols <- persnet_row %>% dplyr::select(name1length:name15length)
  
  # Calculate and return the proportion of alters within the specified time category
  return(
    length(which(length_cols == time_map[[time_category]])) /
      sum(persnet_row %>% dplyr::select(tie1:tie15) != 0, na.rm = TRUE)
  )
}
weak_dur_prop_row <- function(persnet_row) {
  # Proportion of alters who have known the ego for less than 6 years
  return(
    round(
      calc_prop_alters_known_length(persnet_row, "less_than_3years") +
        calc_prop_alters_known_length(persnet_row, "three_to_6years"),
      2
    )
  )
}
weak_dur_prop <- purrr::map_dbl(1:nrow(df_input), ~ weak_dur_prop_row(df_input[.x, , drop = FALSE]))
weak_dur_prop

############################ Distance from Alter ##############################

calc_prop_alters_distance_away <- function(persnet_row, distance_category) {
  # # # # # # # #
  # Function: Computes the proportion of alters who live within a specified  
  #           distance category from the ego.
  # Inputs:  
  #   persnet_row       = A single row of a personal network data frame
  #   distance_category = One of 'same_house', '1_5miles', '6_15miles',  
  #                       '16_50miles', or 'more50miles'
  # Outputs:  
  #   Proportion of alters within the specified distance category
  # # # # # # # #
  
  # Map distance categories to numeric codes
  distance_map <- c(
    same_house = 1,
    "1_5miles" = 2,
    "6_15miles" = 3,
    "16_50miles" = 4,
    more50miles = 5
  )
  
  # Validate the distance_category input
  if (!(distance_category %in% names(distance_map))) {
    stop("Error: Choose one of 'same_house', '1_5miles', '6_15miles', '16_50miles', or 'more50miles'.")
  }
  
  # Identify distance-related columns
  dist_cols_string <- grep("^name\\d+dist$", colnames(persnet_row), value = TRUE)
  
  # Handle missing columns
  if (length(dist_cols_string) == 0) {
    stop("Error: Cannot find columns related to how far away alters live.")
  }
  
  # Select distance columns
  dist_cols <- persnet_row %>% dplyr::select(name1dist:name15dist)
  
  # Calculate and return the proportion of alters within the specified distance
  return(
    length(which(dist_cols == distance_map[[distance_category]])) /
      sum(persnet_row %>% dplyr::select(tie1:tie15) != 0, na.rm = TRUE)
  )
}
far_dist_prop_row <- function(persnet_row) {
  # Proportion of alters who live > 15 miles away
  return(
    round(
      calc_prop_alters_distance_away(persnet_row, "16_50miles") +
        calc_prop_alters_distance_away(persnet_row, "more50miles"),
      2
    )
  )
}
far_dist_prop <- purrr::map_dbl(1:nrow(df_input), ~ far_dist_prop_row(df_input[.x, , drop = FALSE]))
far_dist_prop

############################## Alters' Alcohol #################################

prop_heavy_drinkers_row <- function(persnet_row) {
  # # # # # # # #
  # Function: Computes the proportion of alters who are 'heavy drinkers', defined
  # as those that alters say have or have not count down on drinking in the past
  # six months
  # Inputs:  
  #   persnet_row       = A single row of a personal network data frame
  # Outputs:  
  #   Proportion of alters ego states have or have not cut down on drinking. those 
  # who do not drink heavily or do not drink are considered 'not heavy drinkers'
  # # # # # # # #
  alcohol_cols_string <- grep("^name\\d+alcohol$", colnames(persnet_row), value = TRUE)
  if (length(alcohol_cols_string) == 0) {
    warning("Error: cannot find variables related to whether alters drink alcohol regularly. Returning NA")
    return(NA)
  } else {
    alcohol_cols <- persnet_row %>% dplyr::select(name1alcohol:name15alcohol)
    return(
      round(  
        (length(which(alcohol_cols == 0)) + length(which(alcohol_cols == 1))) /
          sum(persnet_row %>% dplyr::select(tie1:tie15) != 0, na.rm = TRUE),
        2)
    )
  }
}
heavy_drinkers_prop <- purrr::map_dbl(1:nrow(df_input), ~ prop_heavy_drinkers_row(df_input[.x, , drop = FALSE]))
heavy_drinkers_prop

############################## Alters' Smoking #################################

calc_prop_alters_smoke <- function(persnet_row) {
  # # # # # # # #
  # Function: Computes the proportion of alters who do smoke  
  # Inputs:  
  #   persnet_row = A single row of a personal network data frame
  # Outputs:  
  #   Proportion of alters for whom alter answered 1 or 0 to smoking questions
  # # # # # # # #
  
  # Identify smoking-related columns
  smoke_cols_string <- grep("^name\\d+smoke$", colnames(persnet_row), value = TRUE)
  
  # Test if any of smoking columns are missing
  if (length(smoke_cols_string) == 0) {
    warning("Error: Cannot find variables related to whether alters smoke.")
    return(NA)
  }
  
  # Select smoking-related columns
  smoke_cols <- persnet_row %>% dplyr::select(name1smoke:name15smoke)
  
  # Calculate and return the proportion of alters who do smoke
  prop_smoke_value <- sum(length(which(smoke_cols == 0)),length(which(smoke_cols == 1))) /
    sum(persnet_row %>% dplyr::select(tie1:tie15) != 0, na.rm = TRUE)
  return(
    round(prop_smoke_value,2)
  )
}
smoking_prop <- purrr::map_dbl(1:nrow(df_input), ~ calc_prop_alters_smoke(df_input[.x, , drop = FALSE]))
smoking_prop

############################ Alters' Exercise ##################################

calc_prop_alters_exercise <- function(persnet_row) {
  # # # # # # # #
  # Function: Computes the proportion of alters who do not exercise regularly  
  # Inputs:  
  #   persnet_row = A single row of a personal network data frame
  # Outputs:  
  #   Proportion of alters who do not exercise regularly
  # # # # # # # #
  
  # Identify exercise-related columns
  exer_cols_string <- grep("^name\\d+exer$", colnames(persnet_row), value = TRUE)
  
  # Handle missing columns
  if (length(exer_cols_string) == 0) {
    stop("Error: Cannot find variables related to whether alters exercised regularly.")
  }
  
  # Select exercise-related columns
  exer_cols <- persnet_row %>% dplyr::select(name1exer:name15exer)
  
  # Calculate the proportion of alters who do not exercise regularly
  return(
    round((length(which(exer_cols == 0)) /
             sum(persnet_row %>% dplyr::select(tie1:tie15) != 0, na.rm = TRUE)),2)
  )
}
no_exercise_prop <- purrr::map_dbl(1:nrow(df_input), ~ calc_prop_alters_exercise(df_input[.x, , drop = FALSE]))
no_exercise_prop

############################### Alters' Diet ###################################

calc_prop_alters_diet <- function(persnet_row) {
  # # # # # #
  # Function: Computes the proportion of alters who have a bad diet  
  # Inputs:  
  #   persnet_row = A single row of a personal network data frame
  # Outputs:  
  #   Proportion of alters ego knew ``did not eat a healthy diet regularly over the past 3 months'
  # # # # # #
  
  # Identify diet-related columns
  diet_cols_string <- grep("^name\\d+diet$", colnames(persnet_row), value = TRUE)
  
  # Handle missing columns
  if (length(diet_cols_string) == 0) {
    warning("Error: Cannot find variables related to whether alters had a good diet.")
    return(NA)
  }
  
  # Select diet-related columns
  diet_cols <- persnet_row %>% dplyr::select(name1diet:name15diet)
  
  # Calculate and return the proportion of alters with a good diet
  diet_prop_value <- length(which(diet_cols == 0)) /
    sum(persnet_row %>% dplyr::select(tie1:tie15) != 0, na.rm = TRUE)
  
  return(
    round(diet_prop_value,2)
  )
}
bad_diet_prop <- purrr::map_dbl(1:nrow(df_input), ~ calc_prop_alters_diet(df_input[.x, , drop = FALSE]))
bad_diet_prop

######################### Alters' Health Problems ##############################

alter_health_problems_row <- function(persnet_row) {
  alter_health_cols_string <- grep("^name\\d+health___\\d+$", colnames(persnet_row), value = TRUE)
  if (length(alter_health_cols_string) == 0) {
    stop("Error: cannot find variables related to alter health problems.")
  } else {
    alter_health_problem_cols <- persnet_row %>% 
      dplyr::select(name1health___1:name15health___99) %>%
      dplyr::select(!dplyr::contains("___99")) %>% 
      dplyr::select(!dplyr::contains("___0"))
    
    return(
      round(
        sum(alter_health_problem_cols) /
          sum(persnet_row %>% dplyr::select(tie1:tie15) != 0, na.rm = TRUE),
        2
      )
    )
  }
}
health_prob_prop <- purrr::map_dbl(1:nrow(df_input), ~ alter_health_problems_row(df_input[.x, , drop = FALSE]))
health_prob_prop

######## Composition Diversity (gender, edu, dist, length, frequency) #########

blau_gender <- purrr::map_dbl(1:nrow(df_input), ~ calc_blau_alter_heterophily(df_input[.x, , drop = FALSE], attribute = "gender"))
blau_educ <- purrr::map_dbl(1:nrow(df_input), ~ calc_blau_alter_heterophily(df_input[.x, , drop = FALSE], attribute = "educ"))
blau_distance <- purrr::map_dbl(1:nrow(df_input), ~ calc_blau_alter_heterophily(df_input[.x, , drop = FALSE], attribute = "distance"))
blau_length <- purrr::map_dbl(1:nrow(df_input), ~ calc_blau_alter_heterophily(df_input[.x, , drop = FALSE], attribute = "length"))
blau_speak <- purrr::map_dbl(1:nrow(df_input), ~ calc_blau_alter_heterophily(df_input[.x, , drop = FALSE], attribute = "speak"))

blau_gender
blau_educ
blau_distance
blau_length
blau_speak

##################### Constructing/Exporting Data Frame #######################

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
  employment = employment,
  occupation = occupation,
  income = income,
  married = married,
  live_alone = live_alone,
  household_number = household_number,
  ego_alcohol = alcohol,
  ego_smoke = smoke,
  ego_exercise = exercise,
  ego_healty_diet = diet,
  health_problems1 = health_prob1,
  health_problems2 = health_prob2,
  health_problems3 = health_prob3,
  health_problems4 = health_prob4,
  network_size = network_size,
  density = density,
  constraint = constraint,
  effsize = effsize,
  max_degree = max_degree,
  mean_degree = mean_degree,
  kin_prop = kin_prop,
  age_sd = age_sd,
  IQV_sex = iqv_sex,
  IQV_race = iqv_race,
  IQV_educ = iqv_educ,
  weak_freq_prop = weak_freq_prop,
  weak_dur_prop	= weak_dur_prop,
  far_dist_prop	= far_dist_prop,
  heavy_drinkers_prop = heavy_drinkers_prop,
  smoking_prop = smoking_prop,
  no_exercise_prop = no_exercise_prop,
  bad_diet_prop	= bad_diet_prop,
  health_prob_prop = health_prob_prop,
  blau_gender	= blau_gender,
  blau_educ	= blau_educ,
  blau_distance	= blau_distance,
  blau_length	= blau_length,
  blau_speak = blau_speak
)

df_clean$zip <- as.character(df_clean$zip) #double check zip as string

write.csv(df_clean, "Clean_Data.csv", row.names = FALSE)

###################### Single Network Visualizations ###########################

#Note: if outputting individual network graphs through pdf() function begins giving
#  errors of invalid fonts, this may be caused by ggraph's fonts not being used.
#  You may need to change font family under theme_graph() to a font listed at this link
#  https://stat.ethz.ch/R-manual/R-devel/library/grDevices/html/postscriptFonts.html

plot_single_network_node_labels <- function(tidygra) {
  # # # # #
  # Function: Plots a personal network graph using ggraph, distinguishing 
  #           between strong and weak ties and highlighting the ego node.
  # Inputs: tidygra = A tidygraph object representing a personal network
  # Outputs: A ggplot object visualizing the network structure
  # # # # #
  
  # Test if valid tidygra input
  if (is.null(tidygra) || !"tbl_graph" %in% class(tidygra)) {
    warning("Warning: Graph is missing or invalid. Skipping.")
    return(NA)  # Skip invalid graphs
  }
  
  # Test whether the network is an isolate (no weight edge attribute)
  edge_attributes <- tidygra %>% tidygraph::activate(edges) %>% tibble::as_tibble()
  if (!"weight" %in% colnames(edge_attributes)) {
    # Plot isolate (ego only)
    tg_plot <- ggraph(tidygra, layout = "fr") +
      geom_node_point(size = 4, color = 'black', show.legend = FALSE) +
      theme_graph(base_family = "Helvetica-Narrow")
    
    return(tg_plot)
  } else {
    
    # Check if ego is present in the graph, and if so, focus layout on ego
    node_names <- unique(tidygra %N>% dplyr::pull(name))
    
    # Transform tie strength into strong/weak strings and create an alter dummy variable
    tidygra <- tidygra %>%
      tidygraph::activate(edges) %>%
      dplyr::mutate(strength_of_tie = ifelse(weight == 1, "weak", "strong")) %>%
      tidygraph::activate(nodes) %>%
      dplyr::mutate(alter_dummy = ifelse(name != 'ego', 1, 0)) %>%
      dplyr::mutate(node_fill = ifelse(name == "ego", "black", "white")) %>%
      dplyr::mutate(node_text_color = ifelse(name == "ego", "white", "black"))
    
    if ("ego" %in% node_names) {
      focus_index <- which(node_names == "ego")
      record_id <- tidygra %N>% dplyr::pull(record_id) %>% unique() %>% .[1]
      
      # Plot network with ego as focal point
      tg_plot <- ggraph(tidygra, layout = "focus", focus = focus_index) +
        geom_edge_link(aes(color = strength_of_tie, linetype = strength_of_tie),
                       edge_width = 0.75, show.legend = FALSE) +
        scale_edge_linetype_manual(values = c("weak" = "dashed", "strong" = "solid")) +
        scale_edge_color_manual(values = c("weak" = "#5B8FA8FF", "strong" = "#800000FF")) +
        geom_node_point(aes(color = factor(alter_dummy)), size = 4, show.legend = FALSE) +
        scale_colour_manual(values = c('black', 'grey66')) +
        geom_node_label(aes(label = name, fill = node_fill, color = node_text_color),
                        size = 4, label.padding = unit(0.25, "lines"),
                        label.size = 0.5, show.legend = FALSE) +
        scale_color_identity() +
        scale_fill_identity() +
        ggtitle(paste("Record ID:", record_id)) +
        theme_graph(base_family = "Helvetica-Narrow")
      
      return(tg_plot)
    } else {
      
      record_id <- tidygra %N>% dplyr::pull(record_id) %>% unique() %>% .[1]
      
      # Plot network without ego as focal point
      tg_plot <- ggraph(tidygra, layout = "fr") +
        geom_edge_link(aes(color = strength_of_tie, linetype = strength_of_tie),
                       edge_width = 0.75, show.legend = FALSE) +
        scale_edge_linetype_manual(values = c("weak" = "dashed", "strong" = "solid")) +
        scale_edge_color_manual(values = c("weak" = "#5B8FA8FF", "strong" = "#800000FF")) +
        geom_node_point(size = 4, color = 'grey66', show.legend = FALSE) +
        geom_node_label(aes(label = name, fill = node_fill, color = node_text_color),
                        size = 4, label.padding = unit(0.25, "lines"),
                        label.size = 0.5, show.legend = FALSE) +
        scale_color_identity() +
        scale_fill_identity() +
        ggtitle(paste("Record ID:", record_id)) +
        theme_graph(base_family = "Helvetica-Narrow")
      
      return(tg_plot)
    }
  }
}

list_tidygras <- organize_list_tidygraphs(df_input)
list_network_plots_labels <- lapply(list_tidygras, plot_single_network_node_labels)

pdf("Single_Networks.pdf", width = 10, height = 10)
# each print() automatically starts a new page when mfrow = c(1,1)
par(mfrow = c(1,1))

for (plt in list_network_plots_labels) {
  print(plt)
}

dev.off()

###################### Network Montage Visualization ##########################

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

plot_single_network <- function(tidygra, edge_size = 0.5, node_size = 2) {
  # # # # # # # #
  # Function: Plots a personal network graph using ggraph, distinguishing 
  #           between strong and weak ties and highlighting the ego node.
  # Inputs: tidygra = A tidygraph object representing a personal network
  # Outputs: A ggplot object visualizing the network structure
  # # # # # # # #
  
  # Test if valid tidygra input
  if (is.null(tidygra) || !"tbl_graph" %in% class(tidygra)) {
    warning("Warning: Graph is missing or invalid. Skipping.")
    return(NA)  # Skip invalid graphs
  }
  
  # Test whether the network is an isolate (no weight edge attribute)
  edge_attributes <- tidygra %>% tidygraph::activate(edges) %>% tibble::as_tibble()
  if (!"weight" %in% colnames(edge_attributes)) {
    # Plot isolate (ego only)
    tg_plot <- ggraph(tidygra, layout = "fr") +
      geom_node_point(size = node_size, color = 'black', show.legend = FALSE) +
      theme_graph(plot_margin = unit(c(0, 0, 0, 0), "mm"))
    
    return(tg_plot)
  } else {
    
    # Check if ego is present in the graph, and if so, focus layout on ego
    node_names <- unique(tidygra %N>% dplyr::pull(name))
    
    # Transform tie strength into strong/weak strings and create an alter dummy variable
    tidygra <- tidygra %>%
      tidygraph::activate(edges) %>%
      dplyr::mutate(strength_of_tie = ifelse(weight == 1, "weak", "strong")) %>%
      tidygraph::activate(nodes) %>%
      dplyr::mutate(alter_dummy = ifelse(name != 'ego', 1, 0))
    
    if ("ego" %in% node_names) {
      focus_index <- which(node_names == "ego")
      
      # Plot network with ego as focal point
      tg_plot <- ggraph(tidygra, layout = "focus", focus = focus_index) +
        geom_edge_link(aes(color = strength_of_tie, linetype = strength_of_tie),
                       edge_width = edge_size, show.legend = FALSE) +
        scale_edge_linetype_manual(values = c("weak" = "solid", "strong" = "solid")) +
        scale_edge_colour_manual(values = c("weak" = "#5B8FA8FF", "strong" = "#800000FF")) +
        geom_node_point(aes(color = factor(alter_dummy)), size = node_size, show.legend = FALSE) +
        scale_colour_manual(values = c('black', 'grey66')) +
        theme_graph(plot_margin = unit(c(0, 0, 0, 0), "mm"))
      
      return(tg_plot)
    } else {
      
      # Plot network without ego as focal point
      tg_plot <- ggraph(tidygra, layout = "fr") +
        geom_edge_link(aes(color = strength_of_tie, linetype = strength_of_tie),
                       edge_width = edge_size, show.legend = FALSE) +
        scale_edge_linetype_manual(values = c("weak" = "solid", "strong" = "solid")) +
        scale_edge_colour_manual(values = c("weak" = "#5B8FA8FF", "strong" = "#800000FF")) +
        geom_node_point(size = node_size, color = 'grey66', show.legend = FALSE) +
        theme_graph(plot_margin = unit(c(0, 0, 0, 0), "mm"))
      
      return(tg_plot)
    }
  }
}

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
list_network_plots <- lapply(list_tidygras, plot_single_network,
                             edge_size = 0.25 * est_graph_scale,
                             node_size = 1 * est_graph_scale)
#Constructing grid of plots, note this is where we set how many columns the grid will contain,
#  therefore setting the size and aspect ratio of each plot.
wrap_list_network_plots <- wrap_plots(plotlist = list_network_plots, 
                                      ncol = est_column_count)

#Creating output PDF. Note that if you'd like to change the output to something other than a PDF,
#  change the extention within the output file name and uncomment the last two lines of code.
ggsave("Social_Network_Grid.pdf", wrap_list_network_plots, 
       width = output_width,
       height = output_height,
       #This code is only needed in case the user wants to change the output to a raster output (not PDF)
       # units = "in",
       # dpi = 300
       )

############################# Trouble-shooting ################################

# - If you have difficulties, run line-by-line to find the error

# - Look at your raw data and clean or remove rows with missing data. This code
# does not clean your data. It's processes data assuming the survey was filled
# out correctly.

# - If still not working, then email adhand@bwh.harvard.edu.
