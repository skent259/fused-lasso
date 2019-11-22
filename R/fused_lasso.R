library(flsa)
library(png)

## Read in images
lena_noi <- readPNG("./BM3D_images/Lena512_noi_s40.png")
lena_est <- readPNG("./BM3D_images/Lena512_est_s40.png")
lena_true <- readPNG("./BM3D_images/Lena512.png")

house <- readPNG("./BM3D_images/house.png")
house_est <- readPNG("./BM3D_images/house_est_s40.png")

## Set SD
sd <- 40 ## Based on [0,255] grey-scale

sed.seed(8)
house_noi <- house + (sd/255)*matrix(rnorm(prod(dim(house))), nrow = dim(house)[1])
house_noi_2 <- house + noise
writePNG(house_noi, "./BM3d_images/house_noi_s40.png")

## Fused Lasso
lambda2= seq(0.1,0.1, length.out = 10)

system.time(
  fl.house.1 <- flsa(y = house_noi, lambda2 = 0.1)
) ## 10.42 s

system.time(
  fl.house.2 <- flsa(y = house_noi, lambda2 = 0.2)
)

writePNG(fl.house.1[1,,], "./BM3d_images/house_est_s40.png")


psnr <- function(y_true, y_est) {
  # 10*log10(1/mean((y(:)-y_est(:)).^2))
  10*log(base = 10, x = ( 1/mean((y_true-y_est)^2) ) )
}

psnr(house, fl.house.1[1,,] )
psnr(house, house_est)

psnr(lena_true, lena_est)
## Testing
# y2Dim = matrix(rnorm(100), ncol=10)
# #res <- flsa(y, lambda2=lambda2)
# res2Dim <- flsa(y2Dim, lambda2=lambda2)
# 
# 
# resSolObj2Dim <- flsa(y2Dim, lambda2=NULL)
# res2Dim2 <- flsaGetSolution(resSolObj2Dim, lambda2=lambda2)


