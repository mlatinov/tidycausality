
# Function to clean the as
clean_ad_smart <- function(data){
  # Clean the data
  data <- data %>%
    mutate(
      # Drop rare levels from the device
      device_make = as.factor(
        case_when(
          str_detect(device_make,pattern = "Samsung") ~ "Samsung",
          str_detect(device_make,pattern = "XiaoMi") ~ "XiaoMi",
          str_detect(device_make,pattern = "Nokia") ~ "Nokia",
          str_detect(device_make,pattern = "Pixel ") ~ "Pixel ",
          str_detect(device_make,pattern = "OnePlus") ~ "OnePlus",
          str_detect(device_make,pattern = "iPhone") ~   "iPhone",
          .default = "other"
        )
      ),
      # Drop rare levels from the browser
      browser = as.factor(
        case_when(
          str_detect(browser,pattern = "Chrome") ~ "Chrome",
          str_detect(browser,pattern = "Facebook") ~ "Facebook",
          str_detect(browser,pattern = "Safari") ~ "Safari",
          str_detect(browser,pattern = "Samsung Internet") ~ " Samsung_Internet",
          .default = "other"
        )
      ),
      # Create one column for outcome Engaged vs Not Engaged
      outcome = as.factor(if_else(condition = yes == 0 & no == 0,true = "0",false = "1")),

      # Extract the day of the week from the date and encode it as factor
      day = as.factor(wday(ymd(date),label = TRUE)),

      # Change the experiment to factor
      experiment = as.factor(experiment)

    ) %>%
    # Rename
    rename(
      treatment = experiment,
      device = device_make,
    ) %>%
    # Remove the redundant columns
    select(-date,-yes,-no)
  # Return
  return(data)
}
