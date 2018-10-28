library(ggpubr)
require(graphics)

theta <- seq(0, pi, 0.01)
eps <- seq(0, 1, 0.01)
pa <- c("black", colorRampPalette(c("white", "orange"), space="rgb")(101))

dgdc <- function(The, Eps) (3*sin(The)^2 + 2*Eps)/(4*(1 + Eps))

domains <- expand.grid(eps, theta)
names(domains) <- c("eps", "theta")
images <- dgdc(domains['theta'], domains['eps'])

df <- data.frame(domains, images)
names(df) <- c("eps", "theta", "val")

df <- within(df, {
  co <- cos(theta)
})

p <- ggline(df, x = "co", y = "val", color = "eps", group = "eps", plot_type = "l", size=0.8)
p <- p + scale_color_gradientn(colors = pa)
p