load("selfquiz_answers.rda")

# this loads some data and four answers (ans1-ans4)

ls()

# fill in the '...' below with R code using the data in 'df'

head(df)

# 1. what's the mean score of the rows where type == A

ans1 == ...

# 2. what are the IDs of the top 20 rows by score (higher to lower)

all(ans2 == ...)

# 3. There are some IDs in a variable 'my.ids'. What are the scores for these?

head(my.ids)

all(ans3 == ...)

# 4. Make a new score column equal to the original scores, but multiplied by -1 
# for the rows where type == A

df$new.score <- ...
all(ans4 == df$new.score)

