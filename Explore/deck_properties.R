# Assuming your data is in a dataframe called 'formatted_data'

# Create a results dataframe to store the stats for each deck
deck_stats <- data.frame(
  deck = 1:4,
  mean_outcome = numeric(4),
  prop_negative = numeric(4),
  mean_positive = numeric(4),
  mean_loss = numeric(4)
)

# Loop through each deck
for (d in 1:4) {
  # Subset data for this deck where choice is 1 (selected)
  deck_data <- subset(formatted_data, shown == d & choice == 1)
  
  # Mean outcome (only when choice is 1)
  deck_stats$mean_outcome[d] <- mean(deck_data$outcome)
  
  # Proportion of negative outcomes
  # Negative outcomes are when outcome < 0
  deck_stats$prop_negative[d] <- mean(deck_data$outcome < 0)
  
  # Mean positive outcomes (outcome > 0)
  positive_data <- subset(deck_data, outcome > 0)
  deck_stats$mean_positive[d] <- ifelse(nrow(positive_data) > 0, 
                                        median(positive_data$outcome), 
                                        NA)
  
  # Mean loss outcomes (outcome < 0)
  loss_data <- subset(deck_data, outcome < 0)
  deck_stats$mean_loss[d] <- ifelse(nrow(loss_data) > 0, 
                                    median(loss_data$outcome), 
                                    NA)
}

# Print the results
print(deck_stats)
