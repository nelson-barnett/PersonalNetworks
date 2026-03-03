###### Functions related to data visualization

###################### Single Network Visualizations ###########################

#Note: if outputting individual network graphs through pdf() function begins giving
#  errors of invalid fonts, this may be caused by ggraph's fonts not being used.
#  You may need to change font family under theme_graph() to a font listed at this link
#  https://stat.ethz.ch/R-manual/R-devel/library/grDevices/html/postscriptFonts.html

source("metrics.R")

plot_single_network_node_labels <- function(
    tidygra,
    ego_name = NULL,
    fig_title = NULL,
    node_size = 6,
    friend_fill = "#89b53c",
    family_fill = "#007080",
    friend_txt = "white",
    family_txt = "white",
    weak_color = "#236782",
    strong_color = "#c26c21",
    repel = FALSE
) {
    # # # # #
    # Function: Plots a personal network graph using ggraph, distinguishing
    #           between strong and weak ties and highlighting the ego node.
    # Inputs:
    #   tidygra      = A tidygraph object representing a personal network
    #   ego_name     = Character to use instead of "ego" in the graph. Default: NULL (does not replace "ego")
    #   fig_title    = Character to use as the title of the graph. Default: NULL (title is "Record ID: record_id", where "record_id" is extracted from from tidygra)
    #   friend_fill  = Character color with which to fill nodes that contain the text "friend". Default: "#89b53c"
    #   family_fill  = Character color with which to fill nodes that contain the text "family". Default: "#007080"
    #   friend_txt   = Character color with which to write text in nodes that contain the text "friend". Default: "white"
    #   family_txt   = Character color with which to write text in nodes that contain the text "family". Default: "white"
    #   weak_color   = Character color of line indicating weak connections between nodes. Default: "#236782"
    #   strong_color = Character color of line indicating strong connections between nodes. Default: "#c26c21"
    # Outputs: A ggplot object visualizing the network structure
    # # # # #

    # Test if valid tidygra input
    if (is.null(tidygra) || !"tbl_graph" %in% class(tidygra)) {
        warning("Warning: Graph is missing or invalid. Skipping.")
        return(NA) # Skip invalid graphs
    }

    # Test whether the network is an isolate (no weight edge attribute)
    edge_attributes <- tidygra %>%
        tidygraph::activate(edges) %>%
        tibble::as_tibble()
    if (!"weight" %in% colnames(edge_attributes)) {
        # Plot isolate (ego only)
        tg_plot <- ggraph(tidygra, layout = "fr") +
            geom_node_point(size = 4, color = 'black', show.legend = FALSE) +
            theme_graph(base_family = "Helvetica-Narrow")

        return(tg_plot)
    }

    # Check if ego is present in the graph, and if so, focus layout on ego
    node_names <- unique(tidygra %N>% dplyr::pull(name))

    # Transform tie strength into strong/weak strings and create an alter dummy variable
    tidygra <- tidygra %>%
        tidygraph::activate(edges) %>%
        dplyr::mutate(
            strength_of_tie = ifelse(weight == 1, "weak", "strong")
        ) %>%
        tidygraph::activate(nodes) %>%
        dplyr::mutate(
            alter_dummy = ifelse(name != "ego", 1, 0),
            node_fill = dplyr::case_when(
                name == "ego" ~ "black",
                grepl("^Friend", name) ~ friend_fill,
                grepl("^Family", name) ~ family_fill,
                .default = "white"
            ),
            node_text_color = dplyr::case_when(
                name == "ego" ~ "white",
                grepl("^Friend", name) ~ friend_txt,
                grepl("^Family", name) ~ family_txt,
                .default = "black"
            )
        )

    record_id <- tidygra %N>% dplyr::pull(record_id) %>% unique() %>% .[1]
    ttl <- ifelse(is.null(fig_title), paste("Record ID:", record_id), fig_title)

    if ("ego" %in% node_names) {
        focus_index <- which(node_names == "ego")

        # Change "ego" name if a new one was passed
        if (!is.null(ego_name)) {
            tidygra <- tidygra %N>%
                dplyr::mutate(
                    name = ifelse(row_number() == focus_index, ego_name, name)
                )
        }

        # Plot network with ego as focal point
        tg_plot <- ggraph(tidygra, layout = "focus", focus = focus_index) +
            geom_edge_link(
                aes(color = strength_of_tie, linetype = strength_of_tie),
                edge_width = 0.75,
                show.legend = FALSE
            ) +
            scale_edge_linetype_manual(
                values = c("weak" = "dashed", "strong" = "solid")
            ) +
            scale_edge_color_manual(
                values = c("weak" = weak_color, "strong" = strong_color)
            ) +
            geom_node_point(
                aes(color = factor(alter_dummy)),
                size = node_size,
                show.legend = FALSE,
            ) +
            # scale_colour_manual(values = c('black', 'grey66')) +
            geom_node_label(
                aes(label = name, fill = node_fill, color = node_text_color),
                size = node_size,
                label.padding = unit(0.25, "lines"),
                label.size = 0.5,
                show.legend = FALSE,
                repel = repel
            ) +
            scale_color_identity() +
            scale_fill_identity() +
            ggtitle(ttl) +
            theme_graph(base_family = "Helvetica-Narrow") +
            scale_x_continuous(expand = expansion(.15)) +
            scale_y_continuous(expand = expansion(.25))
    } else {
        # Plot network without ego as focal point
        tg_plot <- ggraph(tidygra, layout = "fr") +
            geom_edge_link(
                aes(color = strength_of_tie, linetype = strength_of_tie),
                edge_width = 0.75,
                show.legend = FALSE
            ) +
            scale_edge_linetype_manual(
                values = c("weak" = "dashed", "strong" = "solid")
            ) +
            scale_edge_color_manual(
                values = c("weak" = weak_color, "strong" = strong_color)
            ) +
            geom_node_point(
                size = node_size,
                color = 'grey66',
                show.legend = FALSE
            ) +
            geom_node_label(
                aes(label = name, fill = node_fill, color = node_text_color),
                size = node_size,
                label.padding = unit(0.25, "lines"),
                label.size = 0.5,
                show.legend = FALSE
            ) +
            scale_color_identity() +
            scale_fill_identity() +
            ggtitle(ttl) +
            theme_graph(base_family = "Helvetica-Narrow")
    }
    return(tg_plot)
}

###################### Network Montage Visualization ##########################

plot_single_network <- function(
    tidygra,
    edge_size = 0.5,
    node_size = 2,
    friend_fill = "#89b53c",
    family_fill = "#007080",
    weak_color = "#236782",
    strong_color = "#c26c21"
) {
    # # # # # # # #
    # Function: Plots a personal network graph using ggraph, distinguishing
    #           between strong and weak ties and highlighting the ego node.
    # Inputs:
    #   tidygra      = A tidygraph object representing a personal network
    #   edge_size    = Numeric width of the edges. Default: 0.5
    #   node_size    = Numeric size of the nodes. Default: 2
    #   weak_color   = Character color of line indicating weak connections between nodes. Default: "#236782"
    #   strong_color = Character color of line indicating strong connections between nodes. Default: "#c26c21"
    # Outputs: A ggplot object visualizing the network structure
    # # # # # # # #

    # Test if valid tidygra input
    if (is.null(tidygra) || !"tbl_graph" %in% class(tidygra)) {
        warning("Warning: Graph is missing or invalid. Skipping.")
        return(NA) # Skip invalid graphs
    }

    # Test whether the network is an isolate (no weight edge attribute)
    edge_attributes <- tidygra %>%
        tidygraph::activate(edges) %>%
        tibble::as_tibble()
    if (!"weight" %in% colnames(edge_attributes)) {
        # Plot isolate (ego only)
        tg_plot <- ggraph(tidygra, layout = "fr") +
            geom_node_point(
                size = node_size,
                color = 'black',
                show.legend = FALSE
            ) +
            theme_graph(plot_margin = unit(c(0, 0, 0, 0), "mm"))
    } else {
        # Check if ego is present in the graph, and if so, focus layout on ego
        node_names <- unique(tidygra %N>% dplyr::pull(name))

        # Transform tie strength into strong/weak strings and create an alter dummy variable
        tidygra <- tidygra %>%
            tidygraph::activate(edges) %>%
            dplyr::mutate(
                strength_of_tie = ifelse(weight == 1, "weak", "strong")
            ) %>%
            tidygraph::activate(nodes) %>%
            dplyr::mutate(
                alter_dummy = ifelse(name != 'ego', 1, 0),
                node_fill = dplyr::case_when(
                    name == "ego" ~ "black",
                    grepl("^Friend", name) ~ friend_fill,
                    grepl("^Family", name) ~ family_fill,
                    .default = "grey"
                ),
            )

        if ("ego" %in% node_names) {
            focus_index <- which(node_names == "ego")

            # Plot network with ego as focal point
            tg_plot <- ggraph(tidygra, layout = "focus", focus = focus_index) +
                geom_edge_link(
                    aes(color = strength_of_tie, linetype = strength_of_tie),
                    edge_width = edge_size,
                    show.legend = FALSE
                ) +
                scale_edge_linetype_manual(
                    values = c("weak" = "dashed", "strong" = "solid")
                ) +
                scale_edge_color_manual(
                    values = c("weak" = weak_color, "strong" = strong_color)
                ) +
                geom_node_point(
                    aes(color = node_fill),
                    size = node_size,
                    show.legend = FALSE
                ) +
                scale_colour_identity() +
                theme_graph(plot_margin = unit(c(0.01, 0.01, 0.01, 0.01), "mm"))
        } else {
            # Plot network without ego as focal point
            tg_plot <- ggraph(tidygra, layout = "fr") +
                geom_edge_link(
                    aes(color = strength_of_tie, linetype = strength_of_tie),
                    edge_width = edge_size,
                    show.legend = FALSE
                ) +
                scale_edge_linetype_manual(
                    values = c("weak" = "solid", "strong" = "solid")
                ) +
                scale_edge_colour_manual(
                    values = c("weak" = weak_color, "strong" = strong_color)
                ) +
                geom_node_point(
                    size = node_size,
                    color = 'grey66',
                    show.legend = FALSE
                ) +
                theme_graph(plot_margin = unit(c(0, 0, 0, 0), "mm"))
        }
    }
    return(tg_plot)
}

plot_prop_singleans <- function(
    df,
    mapping,
    attribute,
    attribute_only_cols,
    key_name,
    palette = "Dark2",
    plot_type = "grouped_bar"
) {
    n_alters <- n_alters_with_data(df)
    p_names <- sprintf("Participant %s -- n = %s", 1:nrow(df), n_alters)

    # Get proportions for every mapping value
    props <- purrr::map(names(mapping), \(category) {
        prop_alters_singleans(
            df,
            category,
            mapping,
            attribute,
            attribute_only_cols
        )
    }) %>%
        setNames(names(mapping)) %>%
        dplyr::bind_cols() %>%
        t() %>%
        as.data.frame() %>%
        setNames(p_names) %>%
        tibble::rownames_to_column(key_name)

    # Make plot
    plt <- switch(
        plot_type,
        "pie" = lapply(
            p_names,
            \(column) {
                ggplot(
                    props,
                    aes(
                        x = "",
                        y = !!ensym(column),
                        fill = !!ensym(key_name)
                    )
                ) +
                    geom_bar(stat = "identity", width = 1, color = "white") +
                    coord_polar("y", start = 0) +
                    scale_fill_brewer(palette = palette) +
                    theme_void() +
                    ggtitle(column)
            }
        ) %>%
            gridExtra::grid.arrange(
                grobs = .,
                nrow = ceiling(length(p_names) / 2)
            ),
        "grouped_bar" = props %>%
            pivot_longer(
                cols = -any_of(key_name),
                names_to = "Participant",
                values_to = "Proportion"
            ) %>%
            ggplot(
                aes(fill = Participant, y = Proportion, x = !!ensym(key_name))
            ) +
            geom_bar(position = "dodge", stat = "identity") +
            scale_fill_brewer(palette = palette),
        stop(sprintf("Unknown plot_type argument: %s", plot_type))
    )

    return(plt)
}

plot_prop_multians_piecharts <- function(df, mapping, attribute, key_name) {
    n_alters <- n_alters_with_data(df)
    p_names <- sprintf("Participant %s -- n = %s", 1:nrow(df), n_alters)

    # For each alter
    out_list <- vector(mode = "list", length = 15)
    for (alter_num in 1:15) {
        # Get all of the answers for this attribute
        m <- df %>%
            select(matches(sprintf(
                "^name%s%s_+\\d+$",
                alter_num,
                attribute
            ))) %>%
            dplyr::rename_with(
                ~ names(mapping[mapping == as.integer(gsub("(.*_+)", "", .x))])
            )
        # Get proportion and adjust table
        m <- (m / rowSums(m, na.rm = TRUE)) %>%
            t() %>%
            as.data.frame() %>%
            setNames(p_names) %>%
            tibble::rownames_to_column(key_name)

        # Make piechart
        out_list[[alter_num]] <- lapply(
            p_names,
            \(column) {
                ggplot(
                    m,
                    aes(
                        x = "",
                        y = !!ensym(column),
                        fill = !!ensym(key_name)
                    )
                ) +
                    geom_bar(
                        stat = "identity",
                        width = 1,
                        color = "white"
                    ) +
                    coord_polar("y", start = 0) +
                    scale_fill_brewer(palette = palette) +
                    theme_void() +
                    ggtitle(column)
            }
        )
    }
    return(out_list)

    # gridExtra::grid.arrange(grobs = rev(out_list[[1]]), ncol = 1)
}


plot_prop_multians_scatter <- function(df, mapping, attribute, y_name) {
    plots <- df %>%
        select(matches(sprintf(
            "^name\\d+%s_+\\d+$",
            attribute
        ))) %>%
        tibble::rownames_to_column(var = "Participant") %>%
        pivot_longer(
            -Participant,
            names_to = c("Alter", y_name),
            values_to = "TF",
            names_pattern = "name(\\d+).*_+(\\d+)",
            names_transform = list(y_name = as.integer)
        ) %>%
        mutate(
            Alter = paste("Alter", Alter),
            !!ensym(y_name) := dplyr::recode(
                !!ensym(y_name),
                !!!setNames(names(mapping), mapping)
            )
        ) %>%
        filter(TF == 1) %>%
        group_by(Participant) %>%
        group_map(
            ~ ggplot(., aes(x = Alter, y = !!ensym(y_name))) +
                geom_point(size = 5, color = "Blue") +
                theme_bw() +
                ggtitle(paste("Participant", .y[[1]])) +
                theme(
                    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
                )
        ) %>%
        gridExtra::grid.arrange(grobs = ., nrow = ceiling(nrow(df) / 2))

    return(plots)
}
