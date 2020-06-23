## ----setup, include=FALSE--------------------------------------------------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, fig.align='center')


## ---- eval=F, warning=F, message=F-----------------------------------------------------------------------------------------------------------------------------------
## d0 <- read_csv('../input/hashcode-photo-slideshow/d_pet_pictures.txt')
## kable(str(d0))
## kable(colnames(d0))


## ---- warning=F, message=F-------------------------------------------------------------------------------------------------------------------------------------------
# install packages if necessary
if(!require(corrr))       install.packages("corrr")
if(!require(usedist))     install.packages("usedist")
# if(!require(Biomanager))  install.packages("BiocManager")
if(!require(ComplexHeatmap)) install("ComplexHeatmap")
if(!require(DataCombine)) install.packages("DataCombine")
if(!require(coop))       install.packages("coop")
if(!require(knitr))      install.packages("knittr")
if(!require(kableExtra)) install.packages("kableExtra")

# import library
library(tidyverse)
library(corrr)
library(stringr)
library(lubridate)
library(dplyr)
library(DataCombine)
library(factoextra)
library(dendextend)
library(cluster)
library(pheatmap)
library(usedist)
library(ComplexHeatmap)
library(DataCombine)
library(data.table)
library(coop)
library(knitr)
library(kableExtra)


## ---- warning=F, message=F-------------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------
# Clear memory / specify input file / record time of start
rm(list=ls())
in_file <- '../input/hashcode-photo-slideshow/d_pet_pictures.txt'
sys1 <- Sys.time()


## ---- warning=F, message=F-------------------------------------------------------------------------------------------------------------------------------------------
n <- 128   # preliminary execution
# n <- 1024
# n <- 2048  # preliminary
# n <- 8192
# n <- 90000 # full version
k1 <- 32  # searching margin for vertial pairing
k2 <- 128 # searching margin for whole sequence

## ---- warning=F, message=F-------------------------------------------------------------------------------------------------------------------------------------------
d1 <- read_csv(in_file, skip=1, col_names=F)
d1 <- head(d1, n) 
d1 <- d1 %>% mutate(H = str_extract(X1, '^[H|V]'), 
                    X2=str_extract(X1, '^[H|V]\\s(\\d*)\\s')) %>% 
  mutate(X2=str_replace(X2, '^[H|V] ', '')) %>% 
  mutate(X2=str_replace(X2, '\\s', ''))  %>% 
  mutate(X2=as.integer(X2))  %>% 
  mutate(X1=str_extract(X1,'\\st.*$')) %>%
  mutate(X1=str_replace(X1,'^\\s', ''))  %>% 
  mutate(X0=1:n) %>%
  subset(select=c(4,2,3,1))
knitr::kable(head(d1, 3))


## ---- warning=F, message=F-------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------
# house keeping parameters
#--------------------------------------------
maxs <- max(d1$X2) + 3L
tags <- str_split(d1$X1, '\\s')
tags <- unique(unlist(tags))
tags <- tags[order(tags)]
ntags <- length(tags)
# write_csv(d1, 'b1x.txt')


## ---- warning=F, message=F-------------------------------------------------------------------------------------------------------------------------------------------
#===========================================
# create a big table  => skip usually
#--------------------------------------------
d2 <- sapply(1:(nrow(d1)), function(i){
  y <- str_split(d1$X1[i], '\\s', simplify = T)
  s <- rep(0L, length(tags))
  for(j in 1:(d1$X2[i])){
    yx <- paste('^', y[j], '$', sep='')
    s <- s + ifelse(str_detect(tags, yx), 1L, 0L)
  }
  return(as.integer(s))
})
d2 <- transpose(data.frame(d2))
colnames(d2) <- tags
d2 <- as.tibble(d2)


## ---- warning=F, message=F-------------------------------------------------------------------------------------------------------------------------------------------
# write_csv(d2, 'b2.txt', append = F)
#===========================================
# d2 <- read_csv('b2.txt', col_types=cols(.default = col_integer()))
d2 <- head(d2, n)
knitr::kable(head(d2, 3))


## ---- warning=F, message=F-------------------------------------------------------------------------------------------------------------------------------------------
d3 <- bind_cols(X0=d1$X0, X0V=0, H=d1$H, X2=d1$X2, d2) 
knitr::kable(head(d3,3))


## ---- warning=F, message=F-------------------------------------------------------------------------------------------------------------------------------------------
#============================================
# scores : acutal and max / 
#          individual and cumulative
#============================================
s_sc <- function(dx, total=T){
#############################################
  s1 <- 0
  s2 <- numeric()
  dd <- dx[,tags]

  nr1 <- nrow(dx) -1
  sc2 <-0
  sm2 <-0
  for (i in 1:nr1){

    x <- unlist(dd[i, ])
    y <- unlist(dd[i+1, ])

    j1 <- sum((x - y) == 1)
    j2 <- sum(x * y)
    j3 <- sum((y - x) == 0)

    sc1 <-  min(c(j1,j2,j3))
    sm1 <-  floor(sum(x)/2)
    sc2 <-  sc2 + sc1
    sm2 <-  sm2 + sm1

    s1  <- s1 +sc1
    s2   <- rbind(s2, c(i, sc1, sm1, sc2, sm2))
  }
 ifelse(total, return(s1), return(s2))
}
#------------------------------------------


## ---- warning=F, messag=F--------------------------------------------------------------------------------------------------------------------------------------------
#============================================
# distance of two slides
# x: y: one line of photos
s_di <- function(x,y){
#============================================
  j1 <- sum((x - y) > 0)
  j2 <- sum(x * y)
  j3 <- sum((y - x) > 0)
  #  ntags/(min(c(j1,j2,j3)) +1)
  return(min(c(j1,j2,j3)))
}
#------------------------------------------


## ---- warning=F, message=F-------------------------------------------------------------------------------------------------------------------------------------------
#===========================================
# Reorder/cluster pictures with teh stacks of
# same number of tags
v_cor <- function(d5a){
#------------------------------------------
# order by correlation

  uq <- unique(d5a$X2)
  
  d5o <- d5a[F,]
  for (u in uq){
    # progress report # surpressed for the report
    print('')
    print(Sys.time())
    print(paste('u:', u, sep='')  )
    flush.console()
    d5i <- d5a %>% filter(X2==u)
    if (nrow(d5i)==1) {
      d5o <- bind_rows(d5o, d5i)
      next
    }
    print('transpose...') # progress report # surpresed in report
    flush.console()
    c1 <- transpose(d5i[,tags])
    colnames(c1) <- d5i$X0

    print('correlate...')
    flush.console()
    c2 <- tcrossprod(data.matrix(c1)) %>% as_cordf() %>% rearrange()

    print('mutate...')
    flush.console()
    c3 <- c2 %>% mutate(XR = rownames(c2)) %>%
      rename(X0=rowname) %>% select(X0, XR) %>%
      mutate(X0=as.integer(as.numeric(X0)), XR=as.integer(as.numeric(XR))) %>%
      arrange(X0)

    c4 <- left_join(d5i, c3, by='X0') %>% arrange(XR) %>% select(-XR)
    d5o <- bind_rows(d5o, c4)
  }
  return(d5o)
}


## ---- warning=F, message=F-------------------------------------------------------------------------------------------------------------------------------------------
#===========================================
# Pairng vertical pictures. 
# In return data.fram, vertical pictures are
# counted as one horizontal picture
v_pa <- function(d3){
#------------------------------------------
# First we sort out the vertical photos only
# The horizontal photos d5H is added back to the 
# stack at the end of the function
  
  d5H <- d3 %>% filter(H=='H')
  d5V <- d3 %>% filter(H=='V')

  # Number of vertical photos
  n_d5 <- nrow(d5V)

  # All vertical photos are split into two stacks, 
  # one (na) with more tags (ntags=19...10), 
  # and the other(nb) with less tags (ntags=10...3)
  
  na <- floor(n_d5 / 2L)
  nb <- na

  # Sort the cards in the descending order of number of tags (X2).
  # Correlate cards inside the stack
  # Second ordering key (X0) is necessary so that d5a and d5b
  # does not have an overlap
  
  d5a <- d5V %>% arrange(desc(X2),desc(X0)) %>% head(na)
  d5a <- v_cor(d5a)

  # Sort the cards in the ascending order of number of tags (X2).
  # Correlate cards inside the stack

  d5b <- d5V %>% arrange(X2, X0) %>% head(nb)       
  d5b <- map_df(v_cor(map_df(d5b, rev)), rev)

  # av is the first stack. Only tag part is extracted

  av <- d5a[,tags] 

  print('')
  print('initial preparation finished.')
  flush.console()

  #-------------------------------------
  # Loop for searchign pairs
  # b is a buffer to store the pair
  # i : index for av
  # j : index for bv
  
  d5bi <- d5b
  b <- numeric()
  av0s <- 0

  for (i in 1:nrow(d5a)) {

    # Second stack
    bv <- d5bi[,tags]

    # Staring minimu. We will try to make x1_min = 0
    
    x1_min <- maxs 
    j_min <- 1

    # Look at only first k1 cards.
    # In case that the remaining cards are less than k1
    # stop loop there. 
    
    nbv <- min( c(k1, nrow(d5bi)  ))
    
    # initial values    
    j  <- 1
    j1 <- j
    
    # unlist makes the calculation between av and bv significantly faster
    aa <- unlist(av[i,])
    bb <- unlist(bv[j,])

    # While looking for a match in bv, 
    # if aa changes from even to odd, or vice versa
    # report it
    
    if (av0s != sum(aa)){

        while ((sum(aa) + sum(bb)) %%2 !=0 ){

        j <- j + 1
        bb <- unlist(bv[j,])

        if(j >= nrow(bv) ){
          j <- nrow(bv)
          break
        }
      }
      j1 <- j
    }
    nbv <- ifelse( (j1 + nbv) <= nrow(bv), nbv, nrow(bv) - j1)

    #------------------------------
    # Start looking for good bb
    while(j <= nbv) {
      
      # Candidate of match.
      bb <- unlist(bv[j,])
      
      # If the nubmer of unique tags in aa and bb 
      # does not make up an even number, 
      # just increment j (got to next candidate)
      # wihtout further processing. 
      # Incement nbv as well not to hit the limit. 
      
      while (sum(bb) %% 2 != (av0s %% 2) & j < (j1 + nbv)){
        j <- j + 1
        bb <- unlist(bv[j,])
        nbv <- ifelse((nbv +1 + j) <= nrow(bv), nbv+1, nbv)
        
      }
      # Now the retultant pair should have even number of 
      # unique tags
      
      bb <- unlist(bv[j,])

      # Is there any common tags? 
      x1 <- sum(aa * bb)

      if (x1 < x1_min){
        # If we find better match, replace the current choice

        x1_min <- x1
        j_min <- j
      }
      if (x1==0){
        # If we hit the best score, stop looking for match 
        
        x1_min <- x1
        j_min <- j

        break
      }
      j <- j + 1

    } # j loop (bv) finished

    bb_min <- unlist(bv[j_min,])

    # Store the pair in the buffer
    b <- append(b, c(d5a$X0[i], d5bi$X0[j_min], x1_min, sum(aa)+ sum(bb_min)))

    # Show progress report # but surpressed for report
        print(paste('i:', i, '/', na, ' paired:', d5a$X0[i], '/', d5bi$X0[j_min], ' | ',
                    'j1:', j1, ' nbv:', nbv,' ',format(Sys.time(), "%H:%M:%S"),' | ', 
                    sum(aa), '/', sum(bb_min), ' x1_min: ', x1_min, sep=''))
      flush.console()

    # prepare for the next round
    av0s <- sum(aa)
    
    # Remove the mahing pair from the stack
    d5bi <- d5bi[-j_min,]

  }
  # Shape the buffer b so that one can 
  # easily see the pair
  b <- matrix(b, na, 4, byrow=T)

  #-------------------------------------------
  # Now we know ID pairs which vertical photo
  #  makes a good pair with other vertical photo. 
  # Here actual reordering of pictures.

  for(i in 1:nrow(b)){

    # index in d5a and d5b
    # that mathces the IDs in b
    # b: photo index 
    # iax: index in d5a and d5b
    
    ia <- which(d5a$X0==b[i,1]) 
    ib <- which(d5b$X0==b[i,2])
    
    ax <- unlist(d5a[ia,tags])
    bx <- unlist(d5b[ib,tags])

    cx <- ax + bx
    cx[which(cx >= 1)] <- 1

    # Store peir information in a new column X0V
    d5a[ia,]$X0V <- b[i,2]

    # browser()        
    # Update the number of tags (X2)
    # as the sum of tags in aa and bb
  
    if(i %% 100 == 0){
    print(paste(i, 'cx transpose...,', sep=' ')) # progress report # surpresed in report
    flush.console()
    }
    d5a[ia,tags] <- t(cx)
    d5a[ia,]$X2  <- sum(cx)

  }
  
  # Attach the marged data.set d5a
  # back to the list of the horizontal photos
  # that we separted at the beginning 
  # of this function
  
  d3x <- bind_rows(d5H, d5a)
  return(d3x)
}


## ---- warning=F, message=F-------------------------------------------------------------------------------------------------------------------------------------------
# Pair
d3x <- v_pa(d3)

# Store result in a text file
# so that we can restart teh program from here
# write_csv(d3x, 'b3x.txt')


## ---- warning=F, message=F-------------------------------------------------------------------------------------------------------------------------------------------
#===========================================
# Greedy search
#-------------------------------------------
# How long did it take until here?
sys2 <- Sys.time()
knitr::kable(sys2-sys1)

#-------------------------------------------
# In case we resume from the point 
# that the vertical pairs are created 
# d3x <- read_csv('b3x.txt',col_types=cols(.default = col_integer(), H = col_character()))

#-------------------------------------------
# Shuffle horizontal photos and vertical pairs

set.seed(1978, sample.kind = 'Rounding')
i_rows <- sample(nrow(d3x))
d3x <- d3x[i_rows,]

#-------------------------------------------
# Sort pictures in the descending order of 
# number of tags
d3x <- d3x %>% arrange(desc(X2))

#-------------------------------------------
# Cluster the photos so that 
# we can find a mathching counterpart 
# with in short margin of search.

d3x <- v_cor(d3x)

#-------------------------------------------
# Create placeholders

d3xi <- d3x
d4x <- d3x[F,]

#-------------------------------------------
# This is the starting picture
# We will start the picture 
# that come up the top of the stack
# when we sort the cards according to 
# the number of tags

ax  <- unlist(d3xi[1,tags])

# Remove it from the stack now
d3xi <- d3xi[-1,]

# Store it in submission buffer
d4x[1,] <- d3xi[1,]

#-------------------------------------------
# Initialize cummulative score
s1 <- 0
#-------------------------------------------
# Initialize cummulative maximum (optimal) score
s2 <- 0

#-------------------------------------------
# Here starts the search
# ii is the index of the second cards of 
# transitions

k3 <- k2

for(ii in 2:nrow(d3xi)){

  # As it becomes progressively harder
  # when the number of tags becomes less, 
  # we will extend the margin of search 
  # for pictures that have less tags. 
  
  if(sum(ax) <= 4 ){k3 <- k2 *4}
  if(sum(ax) <= 10){k3 <- k2* 2}

  # The search stops either when we hit
  # we the limit k3, or the stadk of 
  # the card (= pictures) exhausted. 
  
  nr <- min(c(k3, nrow(d3xi)))  # look first 1/8
  
  # Initial distance befor the search
  # We will make a_dist as large as possible
  
  a_dist <-0

  # ix is the index of the second card 
  
  for(ix in 1:nr){
    a_distx <- s_di(ax, unlist(d3xi[ix,tags]))

    # If we find a better match, replace it. 
    if(a_distx > a_dist){
      a_dist <- a_distx
      ix0 <- ix
    }
    
    # If we hit the best possible score
    # stop the search
    if(a_distx == floor(sum(ax)/2)) {
      a_dist <- a_distx
      ix0 <- ix
      break
    }
  }


  # Update the scores
  # s1: actutal score
  # s2: theoretical maximum 
  
  s1 <- s1 + a_dist
  s2 <- s2 + floor(sum(ax)/2)

  # The best match is stored in the buffer d4x  
  d4x <- bind_rows(d4x, d3xi[ix0,] )

  # Set the second card of the transition found this round 
  # as the first card of the next round 

  ax <- unlist(d3xi[ix0,tags])
  
  # Progress report at every 100 trials. 
  # surpressed for report
   if(ii %% 100 == 0){
     print(paste('ii:', ii, '/', nrow(d3xi), ' score:', a_dist, '/', 
                  d3xi[ix0,]$X2, ' | ',  s1, '/', s2, ' || ',      
                format(Sys.time(), "%H:%M:%S"), sep=''))
    flush.console()
   }
  
  # Remove the card already taken 
  # from the stack
  
  d3xi <- d3xi[-ix0,]

  # If there is no card left in the stack,
  # get out of the loop
  
  if (nrow(d3xi) <= 1 | any(is.na(d3xi)) ) {
    print('break...')
    flush.console()
    break
  }
}
# Store the result here, so that we can resume the program 
# from here

# write_csv(d4x, 'b4x.txt')

#--------------------------
sys3 <- Sys.time()
knitr::kable(sys3 - sys2)
#=============================================================
# Prepare submission file 
c4 <- character()
c4 <- toString(nrow(d4x))
for (i in 1:nrow(d4x)){
  c <- ifelse(d4x[i,]$X0V==0, as.character(d4x[i,]$X0-1), 
              paste(as.character(d4x[i,]$X0-1), as.character(d4x[i,]$X0V-1), sep=' '))
  c4 <- append(c4, c)
}
knitr::kable(head(c4,10))
# write_csv(data.frame(c4), 'c4.txt', col_names=F)
#=============================================================
# END
##############################################################


##---- warning=F, message=F-------------------------------------------------------------------------------------------------------------------------------------------
#=============================================================
# The total score
knitr::kable(s_sc(d4x, total=T),  col.names='Final score')



## ---- warning=F, message=F-------------------------------------------------------------------------------------------------------------------------------------------
#=============================================================
# How the score increases
#-------------------------------------------------------------
c5 <-  s_sc(d4x, total=F)
colnames(c5) <- c('i', 'score', 'max_score', 'cum', 'max_cum')
rownames(c5) <- 1:nrow(c5)
c5 <- data.frame(c5)

c5 %>% ggplot() + 
       geom_line(aes(x=i, y=max_cum, color='coral')) + 
       geom_line(aes(x=i, y=cum, color='darkcyan')) + 
       xlab('Number of pictures') + ylab('Scores')+   
       scale_color_discrete(name = "", labels = c("Maximum", "Actual"))


## ---- warning=F, message=F-------------------------------------------------------------------------------------------------------------------------------------------
#=============================================================
# Take a look at the record where the critical loss happen
#-------------------------------------------------------------

c5 %>% ggplot(aes(x=i, y=(max_score - score))) + 
       geom_line(color='coral') + 
       xlab('Number of pictures') + ylab('Loss')

c5 %>% ggplot(aes(x=i, y=(max_cum-cum)/max_cum)) + 
       geom_line(color='coral') + 
       xlab('Number of pictures') + ylab('Fractional loss')  

