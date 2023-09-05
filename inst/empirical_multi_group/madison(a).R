

int = c(-2.44, -0.38, -3.17, -1.86, -1.32, -0.91, -4.62, -0.99, -1.03, -2.08, -1.34, -2.12,
  -0.77, -3.38, -1.17, -0.80, -4.33,-0.82, -3.87, -1.97, -2.38)

me = c(2.67, 2.18, 2.49, 2.18, 3.88, 3.14, 4.93, 1.46, 2.50, 4.59, 
  4.99, 3.89, 2.93, 4.20, 2.23, 2.04, 3.11, 4.58, 4.42, 2.25, 2.22)

se_int = c(.17, .10, .25, .13, .14, .09, .51, .08, .11, .28, .20, .20, 
           .12, .27, .12, .11, .38, .22, 1.43, .18, .16)
se_me = c(.18, .15, .26, .15, .24, .24, .52, .12, .15, .27, .29, .20, .18, .29, .19, .18, .39, .42, 1.43, .18, .20)

plot(int, ylim = c(-5,5), xlab = "Item Number", ylab = "Estimated Value")
points(me, pch = 20)


plotdf = data.frame(int,me,i=1:21)
library(ggplot2)
ggplot(plotdf)+
  geom_point(aes(x=i,y=int,col='Intercepts'), shape = 1, size = 2)+
  geom_point(aes(x=i,y=me,col='Main Effects'), shape = 15, size = 2)+
  geom_errorbar(aes(x=i,ymin=int-2*se_int, ymax=int+2*se_int,col='Intercepts')) +
  geom_errorbar(aes(x=i,ymin=me-2*se_me, ymax=me+2*se_me,col='Main Effects')) +
  coord_cartesian(ylim = c(-6.5, 7)) +
  theme(legend.position = "none") +
  xlab("Item Numbers") +
  ylab(expression(Beta ~ "Estimated Value"))

  