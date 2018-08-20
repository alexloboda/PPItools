library(ggplot2)
library(ggthemes)
library(RColorBrewer)

exps <- 2:8

data <- data.frame()
for (exp in exps) {
  curr <- scan(paste0("res", exp))
  data <- rbind(data, data.frame(exp = exp, value = curr))
}

baseline <- data.frame(exp = exps, val = sapply(exps, function(x) 10 / x))

g <- ggplot(data) +
  geom_boxplot(mapping = aes(exp, value, group = exp), fill = exps, alpha = 0.2, width = 0.2) +
  scale_fill_brewer(palette = "BuPu") +
  theme_bw() +
  xlab("") + ylab("") +
  theme(axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16)) +
  geom_point(data = baseline, shape = 4, size = 2, stroke = 2, mapping = aes(exp, val), color = "darkred")
g
