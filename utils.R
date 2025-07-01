##### Helpful utility functions

get_codebook_mapping <- function(df, field_name) {
    # # # # # # # #
    # Function: Extracts labels and values from redcap codebook
    # Inputs:
    #   df = Dataframe codebook (typically read in via read.csv())
    #   field_name = `Character` the value in the column "Variable / Field Name" from which to extract keys and values
    # Outputs: `named double` with character labels as names and integers as values
    # # # # # # # #

    this_row <- dplyr::filter(df, Variable...Field.Name == field_name)

    # Handle yes/no questions. REDCap standard is no=0, yes=1
    if (purrr::is_empty(this_row$Field.Type)) {
        stop(sprintf("Cannot find field_name `%s` in df", field_name))
    } else if (this_row$Field.Type == "yesno") {
        vals <- c(0, 1)
        keys <- c("no", "yes")
    } else {
        # REDCap standard is to separate options with | and keys/values with ,
        items_list <- this_row %>%
            dplyr::pull(Choices..Calculations..OR.Slider.Labels) %>%
            strsplit(, split = "|", fixed = TRUE) %>%
            unlist() %>%
            strsplit(split = ",", fixed = TRUE) %>%
            lapply(trimws)

        # Separate keys and values
        vals <- purrr::map(items_list, 1) %>% unlist() %>% as.numeric()
        keys <- purrr::map(items_list, 2) %>% unlist()
    }
    return(setNames(vals, keys))
}

names_to_relat <- function(df, mapping) {
    # # # # # # # #
    # Function: Returns a new dataframe with the names replaced by unique relationship labels
    # Inputs: 
    #   df = persnet `dataframe`
    #   mapping = `named double` with relationship labels as names and integers as values (typically from get_codebook_mapping)
    # Outputs: 
    #   `dataframe` the original dataframe with alter names replaced by their relationship label. 
    #   If there are multiple relationships for the same alter, the first one is used.
    # # # # # # # #
   
    F <- function(persnet_row) {
        # Go through each "nameNUM" column
        for (i in 1:15) {
            # Find corresponding nameNUMrelat___NUM columns
            relat_cols <- dplyr::select(
                persnet_row,
                dplyr::matches(sprintf("^name%srelat_*\\d+$", i))
            )
            if (!any(relat_cols)) {
                next
            }
            # Take the first one that == 1
            label <- names(mapping)[
                mapping ==
                    sub(".*_", "", names(relat_cols[which(relat_cols == 1)][1]))
            ]

            # Check if label exists
            col_with_label <- grep(sprintf("^%s$", label), persnet_row)
            # If label exists, add a 1 suffix to it
            if (length(col_with_label) > 0L) {
                persnet_row[[col_with_label]] <- sprintf("%s_1", label)
            }

            # If there are already numbered labels, make one with a unique suffix
            if (length(grep(sprintf("^%s_", label), persnet_row)) > 0L) {
                # Loop until new_name has a unique suffix
                suffix <- 2
                new_name <- sprintf("%s_%s", label, suffix)
                while (
                    new_name %in%
                        dplyr::select(persnet_row, c(sprintf("name%s", 1:15)))
                ) {
                    suffix <- suffix + 1
                    new_name <- sprintf("%s_%s", label, suffix)
                }
            } else {
                # This is first time seeing this label, so make that the name
                new_name <- label
            }

            # Replace the value of this name
            persnet_row[[sprintf("name%s", i)]] <- new_name
        }
        return(persnet_row)
    }

    # Apply the function to each row of df and reconstruct the table
    return(df %>% rowwise() %>% summarize(F(across()))) #### TODO: Update F(across()) (depreciated now)
}
