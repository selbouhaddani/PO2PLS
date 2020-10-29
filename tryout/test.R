library(OmicsPLS)

N=20;p=10;q=20
r=2
rx=1
ry=1

x <- matrix(rnorm(p*N), N,p)
y <- matrix(rnorm(q*N), N,q)
z <- rnorm(N)

params_true <- generate_params(x, y, z, r, rx, ry, alpha = 0.01, type='o2m')
params_true$a <- 0*t(3:2)
params_true$b <- 0*t(2:1)
params_true$sig2G <- 0*0.001
dat <- generate_data(N, params_true)

dim(dat)

X <- dat[,1:p]
Y <- dat[,(p+1):(p+q)]
Z <- dat[,p+q+1]

rm(x,y,z)

fit <- Su_PO2PLS(X, Y, Z, r, rx, ry, steps = 1000, init_param = params_true)
#fit$debug %>% sapply(function(e) e$a[1]) %>% c(params_true$a[1],.) %>% plot
#abline(h=params_true$a[1])

fit_o2 <- PO2PLS::PO2PLS(X, Y, r, rx, ry, steps = 1000, init_param = params_true)


E_step(X, Y, Z, params_true) %>% str
PO2PLS::E_step(X, Y, params_true)[-(1:2)] %>% str

# Shh different in E steps
