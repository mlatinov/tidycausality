
explore_causal <- function(data, outcome = NULL, treatment = NULL, covariate = NULL, type){

  # Flags
  has_outcome   <- !is.null(outcome)
  has_treatment <- !is.null(treatment)
  has_covariate <- !is.null(covariate)

  # Case assignment
  case <- case_when(
    has_covariate  &  !has_treatment   &  !has_outcome ~ "covariate_only",
    has_covariate  &   has_treatment   &  !has_outcome ~ "covariates_treatment",
    has_covariate  &  !has_treatment   &   has_outcome ~ "covariate_outcome",
    has_treatment  &  !has_covariate   &   has_outcome ~ "treatment_outcome",
    has_covariate  &   has_treatment   &   has_outcome ~ "full_causal"
  )

  # Define a global theme and color scheme
  causal_theme <- function() {
    theme_minimal(base_size = 14, base_family = "Helvetica") +
      theme(
        plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
        plot.subtitle = element_text(size = 12, hjust = 0.5, color = "grey40"),
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 11),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "grey85", linewidth = 0.3),
        strip.text = element_text(face = "bold", size = 13, color = "grey25")
      )
  }
  #### Distribution plot Function ####
  distribution_plot <- function(){

    if (case == "covariate_only") {

      # Density plot
      if (is.numeric(data[[covariate]])) {
        plot <-
          ggplot(data, aes(x = .data[[covariate]])) +
          geom_density(fill = "steelblue") +
          labs(
            title = paste0("Distribution of ", covariate ," "),
            y = "Density",
            x = covariate
            )+
          causal_theme()+
          theme(legend.position = "none")
        # Barplot
      } else if (is.factor(data[[covariate]])) {
        plot <-
          ggplot(data, aes(x = fct_infreq(.data[[covariate]]),fill = .data[[covariate]])) +
          geom_bar()+
          labs(
            title = paste0("Bar plot of " , covariate ," "),
            x = covariate,
            y = "Count"
          )+
          scale_fill_viridis_d(option = "A")+
          causal_theme()
      } else {
        stop("Covariate must be numeric or factor")
      }

    } else if (case == "covariates_treatment") {

      # Boxplot
      if (is.numeric(data[[covariate]])) {
        plot <- ggplot(data, aes(x = .data[[covariate]], fill = .data[[treatment]])) +
          geom_boxplot()+
          labs(
            title = paste0("Boxplot " ,covariate , " by ", treatment ," " )
          )+
          coord_flip()+
          scale_fill_viridis_d(option = "A")+
          causal_theme()

        # Stacked Bar plot
      } else if (is.factor(data[[covariate]])) {
        plot <-
          ggplot(data, aes(x = .data[[covariate]], fill = .data[[treatment]])) +
          geom_bar(position = "fill") +
          labs(
            title = paste0("Stacked Bar Plot " ,covariate, "+" , treatment ," "),
            x = covariate,
            y = "Percent",
            fill = treatment
          )+
          scale_y_continuous(labels = scales::percent_format())+
          scale_fill_viridis_d(option = "A")+
          causal_theme()
      }
      # Density plot
    } else if (case == "covariate_outcome") {
      plot <-
        ggplot(data, aes(x = .data[[covariate]], fill = .data[[outcome]])) +
        geom_density()+
        labs(
          title = paste0("Density plot of " ,covariate, " by " , outcome )
        )+
        scale_fill_viridis_d(option = "A")+
        causal_theme()
      # Boxplot
    } else if (case == "full_causal") {
      plot <-
        ggplot(data, aes(x = .data[[covariate]], fill = .data[[outcome]])) +
        geom_boxplot() +
        facet_wrap(vars(.data[[treatment]]))+
        coord_flip()+
        scale_fill_viridis_d(option = "A")+
        labs(
          title = paste0("Boxplot " , covariate , " by " , outcome)
        )+
        causal_theme()

    } else {
      stop("Unsupported case")
    }

    return(plot)
  }

  #### Relationship Function ####
  relationship <- function(){

    if (case == "covariates_treatment") {

      # Scatterplot covariate vs treatment
      if (is.numeric(data[[covariate]])) {
        plot <-
          ggplot(data, aes(x = .data[[covariate]],y = .data[[treatment]],color = .data[[covariate]]))+
          geom_jitter()+
          scale_color_viridis_c(option = "A")+
          labs(
            title = paste0("Jitter plot " , covariate , " by " , treatment)
          )+
          causal_theme()

        # Mosaic
      }else if (is.factor(data[[covariate]])) {
        x_var <- covariate
        fill_var <- treatment

        plot <-
        ggplot(data) +
          geom_mosaic(aes(x = product(!!sym(x_var)), fill = !!sym(fill_var))) +
          theme_mosaic() +
          scale_fill_viridis_d(option = "A")+
          labs(
            title = paste0("Mosaic Plot " ,covariate , " by " , treatment)
          )+
          causal_theme()
      }

    }else if (case == "covariate_outcome"){

      # Scatterplot covariate vs outcome with smooth line
      if (is.numeric(data[[covariate]])) {
        plot <-
          ggplot(data,aes(x = .data[[covariate]],y = .data[[outcome]]))+
          geom_point()+
          scale_color_viridis_c(option = "A")+
          labs(
            title = paste0("Scatter plot " , covariate , " by " , outcome)
          )+
          causal_theme()
        # Boxplot of outcome by covariate levels
      }else if (is.factor(data[[covariate]])) {
        plot <-
          ggplot(data,aes(x = .data[[outcome]],fill = .data[[covariate]]))+
          geom_boxplot()+
          coord_flip()+
          scale_fill_viridis_d(option = "A")+
          labs(
            title = paste0("Boxplot " , covariate , " by " , outcome)
          )+
          causal_theme()

      }else if (case == "full_causal") {

        # Covariate vs outcome scatterplot, colored by treatment
        if (is.numeric(data[[covariate]] & is.numeric(data[[outcome]]))) {
          plot <-
            ggplot(data,aes(x = .data[[covariate]],y = .data[[outcome]],color = .data[[treatment]]))+
            geom_point()+
            scale_color_viridis_c(option = "A")+
            labs(
              title = paste0("Scatter plot " , covariate , " by " , outcome)
            )+
            causal_theme()
        # Boxplot of outcome by covariate, faceted/colored by treatment
        }else if (is.factor(data[[covariate]] & is.numeric(data[[outcome]]))) {
          plot <-
            ggplot(data,aes(x = .data[[outcome]], fill = .data[[covariate]]))+
            geom_boxplot()+
            facet_wrap(vars(.data[[treatment]]))+
            scale_fill_viridis_d(option = "A")+
            labs(
              title = paste0("Boxplot " , outcome , " by " , covariate)
            )+
            causal_theme()

          # Boxplot covariate by outcome, faceted by treatment
        }else if (is.numeric(data[[covariate]] & is.factor(data[[outcome]]))) {
          plot <-
            ggplot(data,aes(x = .data[[covariate]], fill = .data[[outcome]]))+
            geom_boxplot()+
            facet_wrap(vars(.data[[treatment]]))+
            scale_fill_viridis_d(option = "A")+
            labs(
              title = paste0("Boxplot " , outcome , " by " , covariate)
            )+
            causal_theme()

          # Stacked bar plot outcome vs covariate, faceted by treatment
        }else if (is.factor(.data[[covariate]] & is.factor(outcome))) {
          plot <-
            ggplot(data, aes(x = .data[[outcome]], fill = .data[[covariate]])) +
            geom_bar(position = "fill") +
            facet_wrap(vars(.data[[treatment]]))+
            labs(
              title = paste0("Stacked Bar Plot " ,outcome, "+" , covariate ," "),
              x = covariate,
              y = "Percent",
              fill = treatment
            )+
            scale_y_continuous(labels = scales::percent_format())+
            scale_fill_viridis_d(option = "A")+
            causal_theme()
        }
      }
    }

    # Return Plot
    return(plot)
  }

  #### DAG Function ####
  dag <- function() {

    # Implement Hill Climimg Automated DAG
    bn <- hc(data, score = "bic-cg")

    # Convert bnlearn DAG to igraph
    edges <- as.data.frame(arcs(bn))
    g <- graph_from_data_frame(edges, directed = TRUE)

    # Annotate edge types
    edges$edge_type <- "Other"
    edges$edge_type[edges$from == treatment & edges$to == outcome] <- "Treatment → Outcome"
    edges$edge_type[edges$to == outcome & edges$from != treatment] <- "Covariate → Outcome"
    edges$edge_type[edges$to == treatment & edges$from != outcome] <- "Covariate → Treatment"

    # Add edge_type as edge attribute to igraph
    E(g)$edge_type <- edges$edge_type

    # Annotate node types
    node_types <- rep("Covariate", vcount(g))
    names(node_types) <- V(g)$name
    node_types[treatment] <- "Treatment"
    node_types[outcome] <- "Outcome"

    # Detect potential confounders: nodes affecting both treatment and outcome
    confounders <- V(g)$name[
      V(g)$name %in% edges$from[edges$to == treatment] &
        V(g)$name %in% edges$from[edges$to == outcome]
    ]

    # Update node type for confounders
    node_types[confounders] <- "Confounder"
    V(g)$type <- node_types

    # Colors for edges and nodes
    edge_colors <- c(
      "Treatment → Outcome" = "red",
      "Covariate → Outcome" = "orange",
      "Covariate → Treatment" = "purple",
      "Other" = "gray50"
    )

    node_colors <- c(
      "Treatment" = "red",
      "Outcome" = "green",
      "Covariate" = "steelblue",
      "Confounder" = "orange"
    )

    # Plot DAG
    plot <-
      ggraph(g, layout = "tree") +
      geom_edge_link(aes(color = edge_type),
                     arrow = arrow(length = unit(4, 'mm')),
                     end_cap = circle(3, 'mm'),
                     width = 0.8) +
      geom_node_point(aes(color = type), size = 8, alpha = 0.9) +
      geom_node_text(aes(label = name), vjust = -1.2, size = 4, fontface = "bold") +
      scale_color_manual(values = node_colors) +
      scale_edge_color_manual(values = edge_colors) +
      ggtitle("Causal DAG with Potential Confounders") +
      theme_void() +
      theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))

    # Return
    return(plot)
  }
  # Call function depending on type
  plot <- switch(type,
                 "distribution" = distribution_plot(),
                 "relationship" = relationship(),
                 "dag"          = dag(),
                 stop("Invalid type")
  )

  return(plot)
}










# Plot Checked

#### Distributions ####

# 1 Bar plot Covariate Factor
explore_causal(data = data,covariate = "browser" ,type = "distribution")
# 2  Density plot Covariate
explore_causal(data = data,covariate = "platform_os" ,type = "distribution")
# 3  Stacted Bar plot Covariate + treatment factor
explore_causal(data = data,covariate = "browser",treatment = "treatment",type = "distribution")
# 4 Boxplot Covariate + treatment int
explore_causal(data = data,covariate = "hour",treatment = "treatment",type = "distribution")
# 5 Covariate + treatment + outcome
explore_causal(data = data,covariate = "hour",treatment = "treatment",outcome = "outcome",type = "distribution")

#### Relationship ####

# 1 Covarite + Treatment numeric
explore_causal(data = data,covariate = "hour",treatment = "treatment",type = "relationship")
# 2 Covariate + Treatment Factor
explore_causal(data = data,covariate = "browser",treatment = "treatment",type = "relationship")

#### DAG ####
explore_causal(data = data,covariate = "browser",outcome = "outcome",treatment = "treatment",type = "dag")



