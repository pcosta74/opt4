#
# I/O Constants
#
.__FILE__ <- (function() {
  attr(body(sys.function()), "srcfile")
})()$filename

PWD   <- dirname(.__FILE__)
FILE1 <-'cuttingstock.csv'
# FILE2 <-'mem1.csv'
# FILE3 <-'mem2.csv'
# FILE4 <-'neighbours.csv'
# FILE5 <-'output.csv'

#
# Problem constants
#

CH.Sbest <- c(1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,1,1,1)
Attr.S0  <- c(0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,1,1,1,1)
Skew.L   <- c(1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,1,0,1,1,0,0,0)
Skew.M   <- c(0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0)
Skew.R   <- c(0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,1,1,1,1)

AR   <- 1500
DATA <- transform(read.csv2(file.path(PWD,FILE1)),
                  vi=as.numeric(as.character(vi)))

NCOL <- length(Attr.S0)
ECOL <- 2


#
# Component Swap neighbourhood
#
component.swap <- function(Scurr) {
  
  bits.in   <- which(Scurr == 1)
  bits.out  <- which(Scurr == 0)
  
  Acurr <- Scurr %*% DATA$ai
  Vcurr <- Scurr %*% DATA$vi
  Nghd.Scurr <- matrix(c(Scurr, Acurr, Vcurr), ncol=NCOL+ECOL)
  
  for(i in bits.in) {
    for(o in bits.out) {
      Sngh <- Scurr
      Sngh[i] <- 0
      Sngh[o] <- 1
      Angh <- Sngh %*% DATA$ai
      Vngh <- Sngh %*% DATA$vi
      if(Angh <= AR) { 
        Nghd.Scurr <- rbind(Nghd.Scurr, c(Sngh,Angh,Vngh))
      }
    }
  }
  
  df.nghd <- as.data.frame(unique(Nghd.Scurr))
  colnames(df.nghd) <- c(1:length(Scurr),'A','OFV')
  
  return(df.nghd)
}

#
# Bit-flip neighbourhood
#
bit.flip <- function(Scurr, strategy=c('add','sub')) {
  stopifnot(strategy %in% c('add','sub'))
  
  bits.in  <- which(Scurr == 1)
  bits.out <- which(Scurr == 0)
  
  if(strategy == 'add') {
    # if strategy=='add' then bit.set = bits.out ...
    bit.set <- bits.out
    S <- .bf.move.sub
  } else {
    # ... else bit.set = bits.in
    bit.set <- bits.in
    S <- .bf.move.add
  }

  Acurr <- Scurr %*% DATA$ai
  Vcurr <- Scurr %*% DATA$vi
  
  Nghd.Scurr <- matrix(c(Scurr, Acurr, Vcurr), ncol=NCOL+ECOL)
  for(b in bit.set) {
    Sngh    <- Scurr    # create neighbour
    Sngh[b] <- !Sngh[b] # bit-flip
    Nghd.Scurr <- rbind(Nghd.Scurr, S(Sngh, c(b))) # execute strategy
  }

  df.nghd <- as.data.frame(unique(Nghd.Scurr))
  colnames(df.nghd) <- c(1:length(Scurr),'A','OFV')
  
  return(df.nghd)
}

# Subtract strategy
.bf.move.sub <- function(Scurr, except) {
  Nghd.Scurr <- matrix(integer(0), ncol=NCOL+ECOL)
  
  bits.in <- which(Scurr == 1)
  for(b in setdiff(bits.in,except)) {
    Sngh    <- Scurr    # create neighbour
    Sngh[b] <- 0        # bit-flip
    Angh    <- Sngh %*% DATA$ai
    Vngh    <- Sngh %*% DATA$vi
    
    if(Angh <= AR) { 
      # obeys restriction, add to neighbourhood
      Nghd.Scurr <- rbind(Nghd.Scurr, c(Sngh, Angh, Vngh))
    } else {
      # keep on removing
      Nghd.Scurr <- rbind(Nghd.Scurr, .bf.move.sub(Sngh, except))
    }
  }
  
  return(Nghd.Scurr)
}

# Add strategy
.bf.move.add <- function(Scurr, except) {
  Nghd.Scurr <- matrix(integer(0), ncol=NCOL+ECOL)
  
  is.bound <- FALSE
  bits.out <- which(Scurr == 0)
  for(b in setdiff(bits.out,except)) {
    Sngh    <- Scurr    # create neighbour
    Sngh[b] <- 1        # bit-flip
    Angh    <- Sngh  %*% DATA$ai
    Acurr   <- Scurr %*% DATA$ai
    Vcurr   <- Scurr %*% DATA$vi
    
    if(Angh <= AR) {
      # keep on ading
      Nghd.Scurr <- rbind(Nghd.Scurr, .bf.move.add(Sngh, except))
    } else if(!is.bound) {
      # obeys restriction, add to neighbourhood
      Nghd.Scurr <- rbind(Nghd.Scurr, c(Scurr, Acurr, Vcurr))
      is.bound <- TRUE
    }
  }
  
  return(Nghd.Scurr)
}

#
# Find best neighbour
#
best.neighbour <- function(Nghd, T, Scurr, Vbest, `%op%`) {
  if(length(T)>0) {
    L <- list()
    cols <- seq(NCOL)
    for(i in 1:nrow(Nghd)) {
      l <- unique(which(Nghd[i,cols] %op% Scurr))
      if(USE.TYPE) { l <- DATA$ti[l] }
      L <- append(L, list(l))        # flipped bits in the neighbourhood 
    }
    
    R <- c()
    for(i in seq_along(L)) {
      if(any(L[[i]] %in% T)) {
        R <- c(R,i)                  # tabu list exclusion
      }
    }

    AC <- R[Nghd$OFV[R] > Vbest]     # Overriden by Asp. Criterion (AC)
    R  <- R[!(R %in% AC)]
    Nghd$OFV[R] <- -Inf
  }
  which.max(Nghd$OFV[-1]) + 1        # first neighbout == Scurr
}

#
# main
#
ITER     <- 200
USE.TYPE <- FALSE

Sbest <- c()
Abest <- 0
Vbest <- 0

T       <- list()
Scurr   <- Skew.R
Acurr   <- as.integer(Scurr %*% DATA$ai)
Vcurr   <- as.integer(Scurr %*% DATA$vi)

iter <- ITER
mem1 <- matrix(NA,nrow=iter,ncol=NCOL,dimnames=list(NULL,1:NCOL))
mem2 <- matrix(0, nrow=NCOL,ncol=NCOL,dimnames=list(NULL,1:NCOL))

while(iter > 0) {
  iter  <- iter-1

  #Nghd <- component.swap(Scurr)
  #`%op%` <- `!=` 

  #Nghd <- bit.flip(Scurr,strategy='sub')
  #`%op%` <- `>`
  
  Nghd <- bit.flip(Scurr,strategy='add')
  `%op%` <- `<` 
  #write.csv2(Nghd, file.path(PWD,FILE4))

  n <- best.neighbour(Nghd, T, Scurr, Vbest, `%op%`=`%op%`)
  Sngh <- as.numeric(Nghd[n, 1:25])
  Vngh <- as.numeric(Nghd$OFV[n])
  Angh <- as.numeric(Nghd$A[n])
  
  n <- which(Scurr %op% Sngh)
  if(USE.TYPE) { n <- DATA$ti[n] }
  attr <- list(n)

  if(Vngh > Vbest) {
    Sbest <- Sngh
    Abest <- Angh
    Vbest <- Vngh
  }

  T <- append(T,attr)
  if(length(T)==8) { T <- T[2:8]} # 7-item long tabu list

  add <- which(Scurr < Sngh)
  sub <- which(Scurr > Sngh)
  mem2[add, sub] <- mem2[add, sub] + 1
  # write.csv2(mem2, file.path(PWD,FILE3))
  
  Scurr <- Sngh
  Acurr <- Angh
  Vcurr <- Vngh
  
  mem1[ITER - iter,] <- Scurr
  # write.csv2(mem1, file.path(PWD,FILE2))
}

