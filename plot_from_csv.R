library(ggplot2)
library(grid)
library(plyr)
library(readr)

vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)

plot_csv <- function(csv, plot_name) {
  
csv$mut_mod <- factor(csv$mut_mod, levels=c("1","0.8", "0.6", "0.4"))
meds <- ddply(csv, .(mut_mod), summarize, med = median(accuracy))

plot_med <- ggplot(csv, aes(x=mut_mod, y=accuracy)) + 
  coord_cartesian(ylim = c(0, 0.5)) +
  geom_violin(trim=TRUE) +
  geom_boxplot(width=0.1, fill="white")+
  labs(title=plot_name, x="Clone limit (%)", y = "Error") +
  theme(plot.title = element_text(hjust = 1))

return(plot_med)
}

png(filename="plot_experiment.png")

perfect <- read_csv("res_exp_perfect.csv")
persistent_kfull <- read_csv("res_exp_persistent_full.csv")
dollo <- read_csv("res_exp_dollo.csv")
caminsokal <- read_csv("res_exp_caminsokal.csv")


plot_perfect <- plot_csv(perfect, "A. Perfect Phylogeny")
plot_pers_full <- plot_csv(persistent_kfull, "B. Persistent Phylogeny")
plot_dollo <- plot_csv(dollo, "C. Dollo Phylogeny")
plot_cs <- plot_csv(caminsokal, "D. Camin-Sokal Phylogeny")

grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 2)))
print(plot_perfect, vp = vplayout(1, 1))
print(plot_pers_full, vp = vplayout(1, 2))
print(plot_dollo, vp = vplayout(2, 1))
print(plot_cs, vp = vplayout(2, 2))

dev.off()
