
explore_causal <- function(data,outcome = NULL,treatment = NULL,covariate = NULL,type){

  # Extract the outcome treatment and covariate
  outcome_col <- data[[outcome]]
  treatment_col <- data[[treatment]]
  covariate_col <- data[[covariate]]

  # Check which one is null and break it into cases
  has_outcome   <- !is.null(outcome)
  has_treatment <- !is.null(treatment)
  has_covariate <- !is.null(covariate)

  case <- case_when(
    has_covariate  &  !has_treatment   &  !has_outcome ~ "covariate_only",
    has_covariate  &   has_treatment   &  !has_outcome ~ "covariates_treatment",
    has_covariate  &  !has_treatment   &   has_outcome ~ "covariate_outcome"
    has_treatment  &  !has_covariate   &   has_outcome ~ "treatment_outcome",
    has_covariate  &   has_treatment   &   has_outcome ~ "full_causal",
  )

  # Function for ploting Bar plots
  distribution_plot <- function(data,outcome = outcome,treatment = treatment,covariate = covariate,case = case){

    # If only covariate specified
    if (case == "covariate_only") {

      # if the type is numeric
      if (is.numeric(covariate)) {

        # Density plot
        density <-
          ggplot(data = data,aes(x = covariate))+
          geom_density()

        # Boxplot
        boxplot <-
          ggplot(data = data,aes(x = covariate))+
          geom_boxplot()

        # Combine the two plots
        plot <- density + boxplot

        # If the type is factor
      }else if(is.factor(covariate)){

        # Barplot
        bar_plot <-
          ggplot(data = data,aes(x = fct_infreq(covariate)))+
          geom_bar()

        # Return barplot
        plot <- bar_plot

      }

      # If only covariate and treatment specified
    } else if (case == "covariates_treatment") {

      # if the covariate is numeric
      if (is.numeric(covariate)) {

        # Boxplot of covariate by treatment
        plot <-
          ggplot(data = data,aes(x = covariate))+
          geom_boxplot()+
          facet_wrap(~treatment)

        # If the covariate is factor
      }else if (is.factor(covarite)) {

        # Stacked bar covariate by treatment
        plot <-
          ggplot(data = data,aes(x = covariate, fill = treatment))+
          geom_bar(position = "fill", alpha = 0.8) +
          scale_y_continuous(labels = scales::percent_format()) +
          labs(
            x = "Covariate",
            y = "Proportion",
            title = "Proportion of Categorical Outcome by Treatment"
          ) +
          scale_fill_manual(values = c("skyblue", "salmon")) +
          theme_minimal()

        # If covariate is not num or factor get error
      }else{
        stop("Covariate must be numeric or factor")
      }

      # If only covariate and outcome specified
    }else if (case == "covariate_outcome") {




    }

  }

    # If the outcome is numeric plot density plot

    p <-
      ggplot(data = data,aes(x = outcome,fill = treatment))+
      geom_density()



    # If the outcome is factor and covariate is factor plot Proportion bar plot

    p <-
      ggplot(data = data,aes(x = covariate,fill = outcome))+
      geom_bar(position = "fill",alpha = 0.8)+
      scale_y_continuous(labels = scales::percent_format())+
      labs(
        x = "Covatiate",
        y = "Proportion",
        title = "Proportion of Categorical Outcome by Covariate"
      ) +
      scale_fill_manual(values = c("skyblue", "salmon")) +
      theme_minimal()








}

explore_causal(data = data,outcome = "outcome",treatment = "treatment",type = "distribution")


data %>%



ggplot(data = data,aes(x = treatment, fill = outcome))+
  geom_bar(position = "fill", alpha = 0.8) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(
    x = "Treatment",
    y = "Proportion",
    title = "Proportion of Categorical Outcome by Treatment"
  ) +
  scale_fill_manual(values = c("skyblue", "salmon")) +
  theme_minimal()


