
#' Explore Potential Causal Relationships in a Dataset
#'
#' This function provides visual exploration tools to investigate causal relationships between covariates, treatment, and outcome. It can produce:
#' - Distribution plots of a single covariate
#' - Relationship plots between covariates and treatment or outcome
#' - Automated DAGs to detect potential confounders using Hill-Climbing structure learning
#'
#' @param data A data.frame containing the variables of interest.
#' @param outcome Optional. Name of the outcome variable.
#' @param treatment Optional. Name of the treatment variable.
#' @param covariate Optional. Name of the covariate variable to explore.
#' @param type A character string specifying the type of exploration:
#' - "distribution": distribution or bar plot of a covariate
#' - "relationship": relationship between covariate and treatment/outcome
#' - "dag": generates a causal DAG and detects potential confounders
#'
#' @return A ggplot2 object (for "distribution" and "relationship") or a ggraph DAG plot (for "dag").
#'
#' @details
#' - Distribution plots:
#' - Numeric covariates: density plots
#' - Factor covariates: bar plots
#' - Relationship plots:
#' - Covariate vs treatment or outcome, with boxplots, scatterplots, or mosaic plots depending on variable types
#' - DAG plotting:
#' - Uses bnlearn::hc() (Hill-Climbing) to learn the DAG structure from the data
#' - Converts DAG to igraph and plots using ggraph
#' - Annotates nodes as Treatment, Outcome, Covariate, or Confounder (nodes affecting both treatment and outcome)
#' - Highlights edges with colors according to type (Treatment → Outcome, Covariate → Outcome, Covariate → Treatment)
#'
#' @section Notes:
#' - The function can automatically identify potential confounders when type = "dag"; these can be directly fed into adjust_confounders() as the confounders argument.
#' - Supports numeric and factor variables for covariates, treatment, and outcome.
#' - Requires packages: ggplot2, ggraph, igraph, bnlearn, forcats, viridis, and ggmosaic.
#'
#' @examples
#' \dontrun{
#' # Distribution of a numeric covariate
#' explore_causal(data = df, covariate = "age", type = "distribution")
#'
#' # Relationship between covariate and treatment
#' explore_causal(data = df, covariate = "income", treatment = "treat", type = "relationship")
#'
#' # DAG with automatic detection of potential confounders
#' dag_plot <- explore_causal(data = df, outcome = "Y", treatment = "treat", covariate = "age", type = "dag")
#' print(dag_plot)
#'
#' # Use detected confounders for adjustment
#' confounders <- explore_causal(data = df, outcome = "Y", treatment = "treat", type = "dag")$confounders
#' adjusted_df <- adjust_confounders(data = df, confounders = confounders, outcome = "Y", treatment = "treat", method = "regression")
#' }
#'
#' @seealso \code{\link{adjust_confounders}}, \code{\link{ggplot2}}, \code{\link{bnlearn}}, \code{\link{ggraph}}
#'
#' @export
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

    # Make sure that the package bnlearn is installed
    if (!requireNamespace("bnlearn", quietly = TRUE)) {
      stop("The 'bnlearn' package is required for method = 'entropy'. Please install it with install.packages('bnlearn').")
    }
    # Make sure that the package igraph is installed
    if (!requireNamespace("igraph", quietly = TRUE)) {
      stop("The 'igraph' package is required for method = 'entropy'. Please install it with install.packages('igraph').")
    }
    # Make sure that the package ggraph is installed
    if (!requireNamespace("ggraph", quietly = TRUE)) {
      stop("The 'ggraph' package is required for method = 'entropy'. Please install it with install.packages('ggraph').")
    }

    # Learn DAG using hybrid algorithm (MMHC)
    bn <- bnlearn::mmhc(data)

    # Bootstrap for edge stability
    boot <- bnlearn::boot.strength(data = data,R = 100,algorithm = "mmhc")
    avg_dag <- bnlearn::averaged.network(boot, threshold = 0.5)

    #  Detect confounders
    detect_confounders <- function(data, treatment, outcome, p_threshold = 0.05) {

      confounders <- c()
      covariates <- setdiff(names(data), c(treatment, outcome))

      # For every covariate run a stat test and extract the p values
      for (cov in covariates) {
        x <- data[[cov]]
        tr <- data[[treatment]]
        y <- data[[outcome]]

        # Association with treatment
        assoc_tr <- NA
        if (is.numeric(x) & is.numeric(tr)) assoc_tr <- cor.test(x, tr)$p.value
        else if (is.numeric(x) & is.factor(tr)) assoc_tr <- summary(aov(x ~ tr))[[1]][["Pr(>F)"]][1]
        else if (is.factor(x) & is.numeric(tr)) assoc_tr <- summary(aov(tr ~ x))[[1]][["Pr(>F)"]][1]
        else if (is.factor(x) & is.factor(tr)) assoc_tr <- chisq.test(table(x, tr))$p.value

        # Association with outcome
        assoc_y <- NA
        if (is.numeric(x) & is.numeric(y)) assoc_y <- cor.test(x, y)$p.value
        else if (is.numeric(x) & is.factor(y)) assoc_y <- summary(aov(x ~ y))[[1]][["Pr(>F)"]][1]
        else if (is.factor(x) & is.numeric(y)) assoc_y <- summary(aov(y ~ x))[[1]][["Pr(>F)"]][1]
        else if (is.factor(x) & is.factor(y)) assoc_y <- chisq.test(table(x, y))$p.value

        if (!is.na(assoc_tr) & !is.na(assoc_y)) {
          if (assoc_tr < p_threshold & assoc_y < p_threshold) confounders <- c(confounders, cov)
        }
      }

      return(confounders)

    }

    confounders_assoc <- detect_confounders(data, treatment, outcome, p_threshold = 0.05)

    # DAG-based confounders (parents of both treatment and outcome)
    edges <- as.data.frame(bnlearn::arcs(avg_dag))
    confounders_dag <- intersect(edges$from[edges$to == treatment],
                                 edges$from[edges$to == outcome])

    # Robust confounders: union of DAG + association
    confounders <- unique(c(confounders_dag, confounders_assoc))

    # Step 4: Build igraph for plotting
    all_nodes <- colnames(data)
    g <- igraph::graph_from_data_frame(edges, directed = TRUE, vertices = all_nodes)

    node_types <- rep("Covariate", length(all_nodes))
    names(node_types) <- all_nodes
    node_types[treatment] <- "Treatment"
    node_types[outcome] <- "Outcome"
    node_types[confounders] <- "Confounder"
    igraph::V(g)$type <- node_types

    edge_colors <- c("Treatment → Outcome" = "red",
                     "Covariate → Outcome" = "orange",
                     "Covariate → Treatment" = "purple",
                     "Other" = "gray50")

    node_colors <- c("Treatment" = "red",
                     "Outcome" = "green",
                     "Covariate" = "steelblue",
                     "Confounder" = "orange")

    plot <- ggraph::ggraph(g, layout = "tree") +
      ggraph::geom_edge_link(aes(color = edge_colors),
                             arrow = arrow(length = unit(4, 'mm')),
                             end_cap = ggraph::circle(3, 'mm')) +
      ggraph::geom_node_point(aes(color = type), size = 8) +
      ggraph::geom_node_text(aes(label = name), vjust = -1.2) +
      scale_color_manual(values = node_colors) +
      theme_void() +
      ggtitle("Causal DAG with Potential Confounders")

    return(list(plot = plot, confounders = confounders, bn = avg_dag))
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




