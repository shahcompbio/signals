x <- cbind(c(0,0, 0.1,0.05, 0.1,0.05, 0,0, 0,0 , 0,0),
      c(0,0, 0.05,0.1, 0.05,0.1, 0,0, 0,0 , 0,0),
      c(0.1,0.05, 0,0, 0.1,0.05, 0,0, 0,0, 0,0),
      c(0.05,0.1, 0,0, 0.05,0.1, 0,0, 0,0, 0,0),
      c(0.1,0.05, 0.1,0.05, 0,0, 0.01,0, 0,0, 0,0),
      c(0.05,0.1, 0.05,0.1, 0,0, 0,0.01, 0,0, 0,0),
      c(0,0, 0,0, 0.01,0, 0,0, 0.1,0.05, 0.1,0.05),
      c(0,0, 0,0, 0,0.01, 0,0, 0.05,0.1, 0.05,0.1),
      c(0,0, 0,0, 0,0, 0.1,0.05, 0,0, 0.1,0.05),
      c(0,0, 0,0, 0,0, 0.05,0.1, 0,0, 0.05,0.1),
      c(0,0, 0,0, 0,0, 0.1,0.05, 0.1,0.05, 0,0),
      c(0,0, 0,0, 0,0, 0.05,0.1, 0.05,0.1, 0,0))

L = diag(colSums(x)) - x
spec <- eigen(L)
fiedler <- spec$vectors[,rank(spec$values)[2]]
fiedler < 0

A <- c(2/3, 1/3, 1/3, 2/3)
B <- c(1/3, 2/3, 2/3, 1/3)
v[fiedler > 0]

create_adj_matrix <- function(A, B){
  v <- c(rbind(A, B))
  m <- matrix(0, nrow = length(v), ncol = length(v))
  
  for (i in 1:length(v)){
    for (j in 1:length(v)){
         m[i,j] <- 1 - abs(v[i] - v[j])
    }
  }
  
  chks <- split(1:length(v), ceiling(seq_along(1:length(v)) / 2))
  for (i in 1:length(chks)){
    m[chks[[i]][1], chks[[i]][1]] <- 0
    m[chks[[i]][2], chks[[i]][2]] <- 0
    m[chks[[i]][1], chks[[i]][2]] <- 0
    m[chks[[i]][2], chks[[i]][1]] <- 0
  }
  
  return(m)
}

get_fiedler_vec <- function(x){
  L = diag(colSums(x)) - x
  spec <- eigen(L)
  fiedler <- spec$vectors[,rank(spec$values)[2]]
  sign <- fiedler > 0
  
  chks2 <- split(sign, ceiling(seq_along(1:length(sign)) / 2))
  myswitch <- unlist(lapply(chks2, which)) == 1
  
  as.vector(myswitch)
}

chks <- split(1:length(v), ceiling(seq_along(1:length(v)) / 2))
for (i in 1:(length(chks) - 1)){
  m[chks[[i]][1], chks[[i + 1]][1]] <- 1 - abs(v[chks[[i]][1]] - v[chks[[i + 1]][1]])
  m[chks[[i]][2], chks[[i+1]][2]] <- 1 - abs(v[chks[[i]][1]] - v[chks[[i + 1]][2]])
  m[chks[[i]][1], chks[[i + 1]][2]] <- 1 - abs(v[chks[[i]][1]] - v[chks[[i + 1]][2]])
  m[chks[[i]][2], chks[[i+1]][1]] <- 1 - abs(v[chks[[i]][2]] - v[chks[[i + 1]][2]])
}

copy_upper <- function(m) {
  m[lower.tri(m)] <- t(m)[lower.tri(m)]
  m
}

m <- copy_upper(m)
random.stuff <- matrix(runif(prod(dim(m)), min = -0.001, max = 0.001), nrow = nrow(m))
x <- abs(m + random.stuff)
x <- copy_upper(x)
x <- m

