combos <- expand.grid(1:2, 3:5)
test <- data.frame(combos[1, ], combos[2, ],
    edges = c(0.4, 0.4, 0, 0, 0, 0)
)
test2 <- test[which(test$edges == 1), ]
cliques <- max_cliques(graph_from_data_frame(test2, directed = FALSE))
non_overlapping_cliques(cliques)
# for (clique in cliques) {
#     print(clique)
#     for (marker in clique) {
#         print(marker)
#     }
# }

set.seed(1)
test <- matrix(runif(1000000, 0, 1), 1000, 1000,
  dimnames = list(1000:2000, 1000:2000))
row(test)
row(test[1:10, 1:10])
?do.call
as.vector(matrix(rep(1:10, 10), 10, 10))
as.vector(matrix(rep(1:10, 10), 10, 10, byrow = TRUE))
?matrix
colnames(test)
# test <- ifelse(test > 0.2 & test < 0.5, 1, 0)
set.seed(1)
pos <- sort(runif(1000, 0, 1000))
test2 <- matrix(rep(0, 4), 2, 2, dimnames = list(1:2, 1:2))
set.seed(1)
test3 <- matrix(runif(25, 0, 1), 5, 5, dimnames = list(2:6, 2:6))
test3 <- ifelse(test3 > 0.2 & test3 < 0.5, 1, 0)
test4 <- matrix(as.vector(test3), 5, 5)
test3[1:3, 1:3]
n <- 20
max_window <- 10
ld_floor <- 0.2
ld_ceiling <- 0.5
adj_mat_range <- matrix(rep(0, length(test)), nrow = nrow(test), 
  ncol = ncol(test))
for (i in 1:nrow(test)) {
  j <- i
  while (
    j < ncol(test) && length(i:j) <= n &&
    pos[j] - pos[i] < max_window) {
      j <- j + 1
  }
  if (i < ncol(test)) {
    i_row <- vector()
    for (k in (i + 1):j) {
      dist <- pos[k] - pos[i]
      ld <- test[i, k]
      ld_scaled <- ld_floor + ((ld_ceiling - ld_floor) *
        (1 - (dist / max_window) ^ 2))
      i_row <- c(i_row, ifelse(ld >= ld_floor && ld <= ld_ceiling, 1, 0))
    }
    adj_mat_range[i, (i + 1):j] <- i_row
  }
}
test[1:10, 1:10]
adj_mat_range[1:10, 1:10]


row(test)
col(test)

length((graph_from_adjacency_matrix(test) %>% max_cliques())[[1]])
cliques_frame <- tibble()
nrow(cliques_frame)
ld_floor <- 0.2
ld_ceiling <- 0.5
ld_diff <- ld_ceiling - ld_floor
max_window <- 10
adjacency_mat <- function (x, i, j) {
  distance <- abs(pos[j] - pos[i])
  if (distance < max_window) {
    ld_scaled <- ld_floor + (ld_diff * (1 - distance / max_window))
    ifelse(x > ld_scaled && x < ld_ceiling, 1, 0)
  } else {
    ifelse(x > ld_floor && x < ld_ceiling, 1, 0)
  }
}
test <- matrix(
  mapply(adjacency_mat, test, row(test), col(test)), nrow = nrow(test))
test

mat <- matrix(c(1, 2, 3, 4), 2, 2)
matrix(mapply(function(x, i, j) x + i + j, mat, row(mat), col(mat)),
  nrow = nrow(mat))

frame_test <- data.frame(x = runif(10, 1, 10), y = runif(10, 15, 30))

frame_test[order()]
order(frame_test)

frame_test[, order(names(frame_test))]

frame_test %>% arrange(names(frame_test))

what <- function () {
  k <- 3
  tryt <- function (n) {
    k <<- k + n
  }
  tryt(1)
  print(k)
}
what()

sub_test <- test[1:10, 1:10]
sub_test[1, 1]
ld_floor <- 0.2
ld_ceiling <- 0.5
pos <- sort(runif(10, 0, 10))
adj_vec_maker <- function (ld, row, col) {
  dist <- abs(pos[row] - pos[col])
  ld_scaled <- ld_ceiling - ((ld_ceiling - ld_floor) * ((1 / dist) ^ 2))
  ifelse(ld >= ld_scaled && ld <= ld_ceiling, 1, 0)
}
adj_vec <- mapply(adj_vec_maker, sub_test, row(sub_test), col(sub_test))
# test this
adj_mat <- matrix(adj_vec, 10, 10)
library(tidyverse)
library(igraph)
cliques <- graph_from_adjacency_matrix(adj_mat, mode = "undirected") %>%
    largest_cliques()
unlist(cliques[[1]])
length(cliques[[1]])

test4 <- tibble(c(1, 2), c(3, 4))
colnames(test4) <- c("x", "y")
test4[1, ]$x