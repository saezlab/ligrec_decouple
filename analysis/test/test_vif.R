bp <- read.table("~/Downloads/bp_example.txt", header = TRUE, row.names = 1)
bp

fit <- lm(BP ~ ., data = bp)
summary(fit)

colinear <- car::vif(fit) %>%
    broom::tidy() %>%
    filter(x > 5) %>%
    pluck("names")

# if any has vif > 5 -> iteratively remove those

refit <- lm(BP ~ ., data = bp %>% select(-colinear))
summary(refit)
car::vif(refit)
