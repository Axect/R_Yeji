# Title     : TODO
# Objective : TODO
# Created by: kavis
# Created on: 18. 10. 24

library(ggplot2)
library(latex2exp)

y <- function(m) (3*sin(m['th'])^2 + 2*m['eps'])/(4*(1 + m['eps']))
y2 <- function(x) (3*sin(x)^2)/(4)

theta <- seq(0, pi, 0.01);
eps <- seq(0, 1, 0.01);

domains <- expand.grid(eps, theta)
names(domains) <- c("eps", "th")
images  <- y(domains)

df <- data.frame(domains, images)
dg <- data.frame(cos(theta), y2(theta))
names(df) <- c("eps", "th", "val")
names(dg) <- c("th", "val")
df <- within(df, {
  co <- cos(th)
})


#=====================================#
#Plotting
#=====================================#
png("Plot2.png", units="in", width=8, height=6,res=600)
p <- ggplot() + geom_line(data=df, aes(x = co, y = val, color = eps, group = eps))
p <- p + geom_line(data=dg, aes(x=th, y=val, color="black"));
p <- p + scale_colour_gradient(name="$\\epsilon$", low="white", high="orange")
p <- p + xlab(TeX('$\\cos(\\theta)$')) + ylab(TeX('$\\frac{1}{\\Gamma}$ $\\frac{d\\Gamma}{d\\cos(\\theta)}$'))
p <- p + theme(plot.title = element_text(hjust = 0.5))
p <- p + theme_bw()
p
dev.off()
