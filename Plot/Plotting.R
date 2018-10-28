library(ggplot2)
library(latex2exp)
library(tikzDevice)

The <- seq(0, pi, 0.1)
Eps <- seq(0, 1, 0.1)

y <- function(The, Eps) (3*sin(The)^2 + 2*Eps)/(4*(1 + Eps))

#col <- colorRampPalette(c("orange", "white"))

cc <- scales::seq_gradient_pal("orange", "white", "Lab")(seq(0,1,length.out=10))

domains <- expand.grid(Eps, The) # expand.grid(fixed value for 1 sweep, sweeped value) 
names(domains) <- c("Eps", "The")
images <- y(domains[,2],domains[,1])
  
df <- data.frame(domains, images, y(domains[,2], domains[1,1]))
names(df) <- c("Eps", "The", "Ima", "One")
df <- within(df, {co <- cos(The)})


tikz(file = "plot_test.tex", width = 6, height = 4)
#=====================================#
#Plotting
#=====================================#
#png("Plot2.png", units="in", width=8, height=6,res=600)
p <- ggplot(data=df, aes(x = co))
p <- p + geom_line(aes(y = Ima, color = Eps, group = Eps)) + scale_colour_continuous(name="$\\epsilon$", low="white", high="orange")# guide_colorbar(barheight = unit(4,"cm"))) 
p <- p + geom_line(aes(y = One), color = "black", size =1.2)
#p <- p + guides(guide_legend(keyheight = unit(10,"cm")))
#p <- p + theme(legend.key.height(3,"cm"))
p <- p + xlab('$\\cos(\\theta)$') + ylab('$\\frac{1}{\\Gamma}$ $\\frac{d\\Gamma}{d\\cos(\\theta)}$') 
p <- p + theme(plot.title = element_text(hjust = 0.5))
p <- p + theme_bw()
p
dev.off()
