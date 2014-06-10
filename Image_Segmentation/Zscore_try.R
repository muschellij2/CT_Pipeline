rm(list=ls())
library(matrixStats)
library(microbenchmark)
source('~/Dropbox/CTR/DHanley/CT_Registration/programs/Zscore.R')

dim = c(512, 512, 35)
img = array(rnorm(prod(dim), mean=4, sd=4), 
            dim=dim)

truth2 = img
for (i in 1:dim(img)[2]) {
  truth2[,i,] = (truth2[,i,]- mean(truth2[,i,]))/sd(truth2[,i,])
}

truth1 = img
for (i in 1:dim(img)[1]) {
  truth1[i,,] = (truth1[i,,]- mean(truth1[i,,]))/sd(truth1[i,,])
}

truth3 = img
for (i in 1:dim(img)[3]) {
  truth3[,,i] = (truth3[,,i]- mean(truth3[,,i]))/sd(truth3[,,i])
}


microbenchmark({try3 = zscore3(img, margin=3)}, unit="ms", times=5)
microbenchmark({try3 = zscore(img, margin=3)}, unit="ms", times=5)

system.time({
  try3 = zscore3(img, margin=3)
  print(all.equal(try3, truth3))
})

system.time({
  try2 = zscore3(img, margin=2)
  print(all.equal(try2, truth2))
})

system.time({
  try1 = zscore3(img, margin=1)
  print(all.equal(try1, truth1))
})
system.time({print(all.equal(zscore(img, margin=3), truth3))})
system.time({print(all.equal(zscore(img, margin=2), truth2))})
system.time({print(all.equal(zscore(img, margin=1), truth1))})

all.equal(try3, truth3)
try2 = zscore3(img, margin=2)
all.equal(try2, truth2)

try1 = zscore3(img, margin=1)
all.equal(try1, truth1)


mat = array(1:(3*5*2), dim= c(3, 5, 2))
margin = 1
if (margin == 3){
  perm = 1:3
  revperm = 1:3
}
if (margin == 2){
  perm = c(1, 3, 2)
  revperm = c(1, 3, 2)
}  
if (margin == 1){
  perm = c(2, 3, 1)
  revperm = c(3, 1, 2)
}
img = aperm(mat, perm)

vec = matrix(img, ncol=dim(mat)[margin])
m = colMeans(vec)
s = colSds(vec)

vecc = (t(vec) - m)/s
vecc = t(vecc)
imgc = array(vecc, 
             dim = dim(img))
imgc = aperm(imgc, revperm)
