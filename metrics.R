##### Functions related to statistics calculations

###################### Utility Functions ######################

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
    psn_edges_values <- as.integer(psn_edges_variables) # convert to integer

    # Initialize an adjacency matrix for the ego network
    psn_mat <- matrix(NA, 16, 16) # create a 16x16 matrix (max 15 alters + ego)
    diag(psn_mat) <- 0 # Set diagonal to 0 (no self-loops)

    # Populate the lower triangular part of the adjacency matrix with tie values
    psn_mat[lower.tri(psn_mat)] <- psn_edges_values

    # Make the matrix symmetric by copying the lower triangle to the upper triangle
    psn_mat <- t(psn_mat)
    psn_mat[lower.tri(psn_mat)] <- psn_edges_values

    # Extract column names for alters from the input row
    psn_name_cols <- dplyr::select(
        df_row_input,
        c(sprintf("name%s", 1:15))
    )

    # Determine which alters to keep based on a "keep" column (binary indicator)
    psn_keep_name_cols <- dplyr::select(df_row_input, name_1:name_15)
    psn_names_to_keep <- c(
        "ego",
        ifelse(psn_keep_name_cols == 1, psn_name_cols, NA)
    )

    # Flatten the names to a vector, ignoring names marked as NA
    psn_names_for_mat <- as.vector(unlist(psn_names_to_keep, use.names = FALSE))
    colnames(psn_mat) <- rownames(psn_mat) <- psn_names_for_mat # Set matrix row/column names

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
            dplyr::mutate(
                record_id = df_row_input$record_id, # Add record ID
                node_id = dplyr::row_number()
            ) # Create a unique node ID

        return(tgra)
    } else {
        # if no edges exist, create an isolate
        tgra <- tidygraph::tbl_graph(
            nodes = tibble::tibble(name = c("ego"))
        ) %>%
            tidygraph::activate(nodes) %>%
            dplyr::mutate(
                record_id = df_row_input$record_id,
                node_id = dplyr::row_number()
            )

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
        dplyr::group_split(index) # Split by the index

    # Convert each row into a tidygraph object using row_to_tidygraph()
    tidygraph_list <- lapply(df_as_list, organize_row_to_tidygraph)

    return(tidygraph_list)
}


n_alters <- function(persnet_row) {
    # # # # # # # #
    # Function: Counts the number of alters in a network
    # Inputs: persnet_row = one row of a persnet dataframe
    # Outputs: `double` the nubmer of alters in the persnet row's network
    # # # # # # # #

    return(sum(
        persnet_row %>% dplyr::select(tie1:tie15) != 0,
        na.rm = TRUE
    ))
}

extract_ego_multi_attributes <- function(
    persnet_row,
    attribute,
    mapping
) {
    # # # # # # #
    # Function: Extract the labels for all values of the attribute in a checkbox/multi answer style field (e.g. race___1, race___2, etc.)
    # Inputs:
    #   persnet_row = A single row of a personal network data frame
    #   attribute = Character name of attribute in question (e.g. "race")
    #   mapping = `named double` with labels as names and integers as values
    # Outputs:
    #   `Character` containing the labels from mapping for each value in the attribute columns
    # # # # # # #

    # Get attribute column names
    cols <- grep(
        sprintf("^%s_*\\d+$", attribute),
        names(persnet_row),
        value = TRUE
    )

    # Get NUM in attribute___NUM for all columns where the value == 1, get label from mapping
    return(
        mapping[
            sub(".*_", "", cols[which(persnet_row[cols] == 1)]) %>%
                as.numeric()
        ] %>%
            names()
    )
}

calc_total_alters_row <- function(persnet_row) {
    ########## NOTE: This has more_names column too, which doesn't align with the graphed network

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
    all_more_names <- paste(
        persnet_row$more_names_1,
        persnet_row$more_names_2,
        persnet_row$more_names_3,
        sep = ", "
    ) %>%
        stringr::str_split(",\\s*") %>% # split by commas
        unlist() %>% # conver to vector
        .[. != ""] %>% # semove empty strings
        unique()

    # Return the count of unique alters from both sources
    return(length(union(unique_names_tidygra, all_more_names)))
}

calc_total_alters_df <- function(persnet_df) {
    # applies total alters row to each dataframe
    tryCatch(
        {
            # apply to each row the function of calculate total alters
            vector_count_total_alters <- apply(persnet_df, 1, function(row) {
                calc_total_alters_row(
                    persnet_df %>% dplyr::filter(record_id == row["record_id"])
                )
            })
            return(vector_count_total_alters)
        },
        error = function(e) {
            # If an error occurs, stop with a message indicating that the input is not a valid persnet dataframe
            warning(
                " error: there was an error calculating number of unqiue alters, returning NA"
            )
            return(NA)
        }
    )
}

############################ Metric Functions ###################################

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
    tryCatch(
        {
            # Organize a list of tidygraphs from the persnet dataframe
            gra_list <- organize_list_tidygraphs(persnet_df)
            # Calculate the density of each network after removing the ego node
            vector_egoless_density <- sapply(gra_list, calc_egoless_density)
            return(vector_egoless_density)
        },
        error = function(e) {
            # If an error occurs, stop with a message indicating that the input is not a valid persnet dataframe
            stop("Error: not a valid persnet dataframe")
        }
    )
}

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
    if (dim(igraph::as_adjacency_matrix(tidygra))[1] < 2) {
        return(NA)
    }

    # check if edge attribute weight is missing
    edge_attrs <- tidygra %>% tidygraph::activate(edges) %>% tibble::as_tibble()
    if (!"weight" %in% colnames(edge_attrs)) {
        warning("Network has no weights; assigning all tie weights as 1.")
        tidygra <- tidygra %>%
            tidygraph::activate(edges) %>%
            dplyr::mutate(weight = 1)
    }

    # Convert to adjacency matrix (undirected) with weight
    adj_matrix <- igraph::as_adjacency_matrix(
        tidygra,
        attr = "weight",
        sparse = FALSE
    )
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
        if (sum_tie_strengths[i] == 0) {
            next
        }
        for (j in setdiff(1:num_nodes, i)) {
            matrix_proportions[i, j] <- adj_matrix[i, j] / sum_tie_strengths[i]
        }
    }

    # Calculate Burt's constraint for ego node
    constraint <- 0
    for (j in setdiff(1:num_nodes, ego_index)) {
        redundancy <- sum(
            matrix_proportions[ego_index, -c(ego_index, j)] *
                matrix_proportions[-c(ego_index, j), j]
        )
        constraint <- constraint +
            (matrix_proportions[ego_index, j] + redundancy)^2
    }

    return(round(constraint, 2))
}

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
    if (dim(igraph::as_adjacency_matrix(tidygra))[1] < 2) {
        return(NA)
    }

    # Check if the edge attribute 'weight' is missing
    edge_attrs <- tidygra %>% tidygraph::activate(edges) %>% tibble::as_tibble()
    if (!"weight" %in% colnames(edge_attrs)) {
        warning(
            "Warning: Network has no weights. Effective Network Size typically includes weight. Calculating with all tie weights set to 1."
        )
        tidygra <- tidygra %>%
            tidygraph::activate(edges) %>%
            dplyr::mutate(weight = 1) # Add default weight of 1
    }

    # convert to adjacency matrix
    adj_matrix <- igraph::as_adjacency_matrix(
        tidygra,
        attr = "weight",
        sparse = FALSE
    )
    node_names <- igraph::V(tidygra)$name # Extract node names from graph

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
        max_tie_strength <- max(adj_matrix[current_node, ], na.rm = TRUE) # Max tie strength for current node
        if (max_tie_strength == 0) {
            next
        } # Skip if no ties exist

        redundancy_sum <- 0 # Initialize redundancy sum for current node
        for (other_neighbor in setdiff(
            1:num_nodes,
            c(ego_index, current_node)
        )) {
            # proportion of ego's total ties to other_neighbor
            ego_proportion <- adj_matrix[ego_index, other_neighbor] /
                total_tie_strength

            # proportion of current_node's strongest tie to other_neighbor
            node_redundancy <- adj_matrix[current_node, other_neighbor] /
                max_tie_strength

            redundancy_sum <- redundancy_sum +
                (ego_proportion * node_redundancy)
        }

        # update effective size result
        effective_size_result <- effective_size_result + (1 - redundancy_sum)
    }
    return(round(effective_size_result, 2))
}

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

    tryCatch(
        {
            # Remove ego
            egoless_graph <- remove_ego_from_igraph(tg_graph)

            # Check if the graph is empty (no nodes)
            if (igraph::vcount(egoless_graph) == 0) {
                return(0)
            }

            # Else calculate mean degree
            return(round(mean(igraph::degree(egoless_graph)), 2))
        },
        error = function(e) {
            # If an error occurs, return NA
            return(NA)
        }
    )
}

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

    tryCatch(
        {
            # Remove ego

            egoless_graph <- remove_ego_from_igraph(tg_graph)

            # Check if the graph is empty, if so return 0
            if (igraph::vcount(egoless_graph) == 0) {
                return(0)
            }

            # Else calculate max degree
            return(max(igraph::degree(egoless_graph)))
        },
        error = function(e) {
            # If an error occurs, return NA
            return(NA)
        }
    )
}

####################### Network Composition Measures ##########################

########################### Age Standard Deviation ############################

sd_age_alters_row <- function(persnet_row) {
    # Compute the standard deviation of the ages of the alters in the egonetwork
    return(
        round(
            stats::sd(
                persnet_row %>% dplyr::select(name1age:name15age) %>% unlist(),
                na.rm = TRUE
            ),
            2
        )
    )
}

####################### Diversity Calculation Functions #######################

calc_blau_alter_heterophily <- function(
    persnet_row,
    attribute,
    mapping
) {
    # # # # # # # #
    # Function: Computes the Blau heterophily index for a specified alter attribute.
    #           The index measures diversity within a personal network.
    # Inputs:
    #   persnet_row = A single row of a personal network data frame
    #   attribute   = `Character` attribute on which to calculate metric
    #   mapping     = `named double` corresponding labels (names) and values (numeric) for the attribute
    # Outputs:
    #   Blau heterophily index for the specified attribute
    # # # # # # # #

    # Blau heterophily assumes mutually exclusive categories, so only use 'singleans'

    # Apply calc_prop_alters_singleans to all names in mapping, returning a list, then get the proportion
    return(
        lapply(names(mapping), function(x) {
            calc_prop_alters_singleans(
                persnet_row,
                x,
                mapping,
                attribute
            )^2
        }) %>%
            unlist() %>%
            {
                function(x) 1 - sum(x)
            }() %>%
            round(digits = 2)
    )
}


calc_attribute_iqv <- function(
    persnet_row,
    attribute,
    mapping
) {
    # # # # # # # #
    # Function: Computes the Index of Qualitative Variation (IQV) for a specified
    #           alter attribute. The IQV measures diversity and variation within
    #           a personal network.
    # Inputs:
    #   persnet_row = A single row of a personal network data frame
    #   attribute   = `Character` attribute on which to calculate metric
    #   mapping     = `named double` corresponding labels (names) and values (numeric) for the attribute
    # Outputs:
    #   IQV value for the specified attribute
    # # # # # # # #

    # calculate IQV using Blau heterophily index and corresponding normalization factor
    return(round(
        calc_blau_alter_heterophily(
            persnet_row,
            attribute,
            mapping
        ) /
            1 -
            (1 / length(mapping[[1]])),
        digits = 2
    ))
}

calc_prop_alters_multians <- function(
    persnet_row,
    categories,
    mapping,
    attribute
) {
    # # # # # # # #
    # Function: Computes the proportion of alters that have a given category or categories
    #   for attributes that may have multiple answers (such as "relationship").
    #   Does not double count in the case where one alter is assigned multiple values in `categories`.
    #
    # Inputs:
    #   persnet_row = A single row of a personal network data frame
    #   categories  = `Character` containing one or multiple categories from which to calcluate proportion (e.g., c("Friend", "Family") OR "Friend")
    #   mapping     = `named double` corresponding labels (names) and values (numeric) for the attribute
    #   attribute   = `Character` attribute in which categories can be found (e.g. "relat")
    # Outputs:
    #   `double` Proportion of alters that have the category/one of the categories in `categories`.
    #   Returns NA if cannot find columns related to `attribute`
    # # # # # # # #

    # Validate the category input
    if (!all((categories %in% names(mapping)))) {
        stop(sprintf(
            "Error: categories `%s` is not in the provided mapping. Must be in: %s",
            paste(categories, collapse = "`, `"),
            paste(names(mapping), collapse = ", ")
        ))
    }

    # Convert row to a tidygraph object and check if the network is an isolate
    if (igraph::vcount(organize_row_to_tidygraph(persnet_row)) == 1) {
        return(NA)
    }

    # Find all columns pertaining to everything in categories
    relevant_cols <- purrr::map(categories, function(x) {
        persnet_row %>%
            dplyr::select(dplyr::matches(sprintf(
                "^name%s%s_*%s$",
                1:15,
                attribute,
                mapping[[x]]
            )))
    }) %>%
        unlist()

    # How many unique names == 1 (avoids double counting categories)
    if (is.null(relevant_cols)) {
        warning(sprintf(
            "Cannot find columns related to `%s`. Returning NA...",
            attribute
        ))
        return(NA)
    } else {
        return(
            length(
                gsub("_(.*)", "", names(which(relevant_cols == 1))) %>% unique()
            ) /
                n_alters(persnet_row)
        )
    }
}

calc_prop_alters_singleans <- function(
    persnet_row,
    categories,
    mapping,
    attribute
) {
    # # # # # # # #
    # Function: Computes the proportion of alters that have a given category or categories
    #   for attributes that may have only have one answer (such as "speak" and "length").
    #
    # Inputs:
    #   persnet_row = A single row of a personal network data frame
    #   categories  = `Character` containing one or multiple categories from which to calcluate proportion (e.g., c("Friend", "Family") OR "Friend")
    #   mapping     = `named double` corresponding labels (names) and values (numeric) for the attribute
    #   attribute   = `Character` attribute in which categories can be found (e.g. "relat")
    # Outputs:
    #   `double` Proportion of alters that have the category/one of the categories in `categories`.
    #   Returns NA if cannot find columns related to `attribute`
    # # # # # # # #

    # Validate the category input
    if (!all((categories %in% names(mapping)))) {
        stop(sprintf(
            "Error: categories `%s` is not in the provided mapping. Must be in: %s",
            paste(categories, collapse = "`, `"),
            paste(names(mapping), collapse = ", ")
        ))
    }

    # Convert row to a tidygraph object and check if the network is an isolate
    if (igraph::vcount(organize_row_to_tidygraph(persnet_row)) == 1) {
        return(NA)
    }

    # Select attribute columns of the proper form
    selected_cols <- persnet_row %>%
        dplyr::select(matches(c(sprintf("^name%s%s$", 1:15, attribute))))

    # Calculate and return the proportion of alters within the specified category
    if (dim(selected_cols)[2] == 0) {
        warning(sprintf(
            "Cannot find columns related to `%s`. Returning NA...",
            attribute
        ))
        return(NA)
    } else {
        return(
            length(which(selected_cols %in% mapping[categories])) /
                n_alters(persnet_row)
        )
    }
}
