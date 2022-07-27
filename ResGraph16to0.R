library(ggplot2)
library(patchwork)

loog500=readRDS("individuals_16to0_m45_g50_LoogPro_res.RData")

dim(loog500)
loog500[,,500]
yearBP=loog500[,2,1]-2000
sfa=loog500[,5,1]
loog=array(data=NA,dim=c(121,4,1))
loog=as.data.frame(loog)
colnames(loog)=c("yearBP","sfa","ec","size")
loog$yearBP=yearBP
loog$sfa=sfa
loog$size=loog500[,4,1]

for(i in 1:121) {
  loog[i,3]=sum(loog500[i,14,seq(1,500,1)]>=loog500[i,7,1])/500
}
range(loog$ec)




theme_set(theme_bw())

p=ggplot(data=loog,aes(x=yearBP,y=sfa))
p=p+geom_line(colour="black",)
#p=p+geom_line(aes(y=V10),colour="darkred",linetype="dashed")
#p=p+geom_vline(xintercept=3500,linetype="dotted",colour="darkred")
p=p+labs(y=expression(paste("Scaling Factor Angle ",alpha)),x=element_blank(),
         title="Mobility in European ancient DNA from 14000 to 2000 years BP")
p=p+theme(axis.title.x=element_text(margin=margin(t=10),size=15),
          axis.title.y=element_text(margin=margin(r=10),size=20),
          plot.title=element_text(margin=margin(r=10,0,10,0),size=34,hjust=0.5),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          plot.margin=margin(t=10,r=10,b=0,l=10),
          axis.text.y=element_text(size=15),
          panel.grid.minor=element_blank())
p=p+scale_x_reverse(breaks=seq(2000,14000,1000),limits=c(14000,2000),expand=expansion(mult=c(.025,.025)))
p=p+scale_y_continuous(breaks=seq(0,60,10),limits=c(0,60))


p1=ggplot(data=loog,aes(x=yearBP,y=ec))
p1=p1+geom_line(colour="black",)
p1=p1+geom_line(aes(y=0.05),colour="darkred",linetype="solid")
p1=p1+labs(x=element_blank(),y="P-value")
p1=p1+theme(axis.title.x=element_text(margin=margin(t=10),size=15),
            axis.title.y=element_text(margin=margin(r=10),size=20),
            plot.title=element_text(margin=margin(r=10,0,0,0),size=15,hjust=0.5),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            plot.margin=margin(t=0,r=10,b=0,l=10),
            axis.text.y=element_text(size=15),
            panel.grid.minor=element_blank())
p1=p1+scale_x_reverse(breaks=seq(2000,14000,1000),limits=c(14000,2000),expand=expansion(mult=c(.025,.025)))
p1=p1+scale_y_continuous(breaks=seq(0,0.6,0.1),limits=c(0,0.6))


p2=ggplot(data=loog,aes(x=yearBP,y=size))
p2=p2+geom_line(colour="black",)
#p2=p2+geom_line(aes(y=0.05),colour="darkred",linetype="dotted")
p2=p2+labs(x="Years BP",y="Sample size")
p2=p2+theme(axis.title.x=element_text(margin=margin(t=10),size=25),
            axis.title.y=element_text(margin=margin(r=10),size=20),
            plot.title=element_text(margin=margin(r=10,0,0,0),size=20,hjust=0.5),
            plot.margin=margin(t=0,r=10,b=10,l=10),
            axis.text.x=element_text(size=15),
            axis.text.y=element_text(size=15),
            panel.grid.minor=element_blank())
p2=p2+scale_x_reverse(breaks=seq(2000,14000,1000),limits=c(14000,2000),expand=expansion(mult=c(.025,.025)))
p2=p2+scale_y_continuous(breaks=seq(0,750,150),limits=c(0,750))


p3=p/p1/p2
p3=p3+plot_layout(heights=c(3,1,1))
print(p3)


### X vs Auto###

auto=readRDS("ind_m45_g50_chr1-22_res_1000.RData")
xchr=readRDS("ind_m45_g50_XnoY_res_1000.RData")

auto[,2,1]
xchr[,2,1]



yearBP=auto[seq(1,41,1),2,1]-1500
sfa=auto[seq(1,41,1),5,1]
autores=array(data=NA,dim=c(41,4,1))
autores=as.data.frame(autores)
colnames(autores)=c("yearBP","sfa","ec","size")
autores$yearBP=yearBP
autores$sfa=sfa
autores$size=auto[seq(1,41,1),4,1]

for(i in 1:41) {
  autores[i,3]=sum(auto[1,14,seq(1,1000,1)]>=auto[i,7,1])/1000
}
head(autores)
autores$ec=1-autores$ec
range(autores$ec)

auto[,7,1]


p=ggplot(data=autores,aes(x=yearBP,y=sfa))
p=p+geom_line(colour="black",)
#p=p+geom_line(aes(y=V10),colour="darkred",linetype="dashed")
#p=p+geom_vline(xintercept=3500,linetype="dotted",colour="darkred")
p=p+labs(y=expression(paste("Scaling Factor Angle ",alpha)),x=element_blank(),
         title="Mobility of autosomes in European ancient DNA")
p=p+theme(axis.title.x=element_text(margin=margin(t=10),size=15),
          axis.title.y=element_text(margin=margin(r=10),size=20),
          plot.title=element_text(margin=margin(r=10,0,10,0),size=35,hjust=0.5),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          plot.margin=margin(t=10,r=10,b=0,l=10),
          axis.text.y=element_text(size=15),
          panel.grid.minor=element_blank())
p=p+scale_x_reverse(breaks=seq(3500,7500,500),limits=c(7500,3500),expand=expansion(mult=c(.025,.025)))
p=p+scale_y_continuous(breaks=seq(0,60,10),limits=c(0,60))


p1=ggplot(data=autores,aes(x=yearBP,y=ec))
p1=p1+geom_line(colour="black",)
p1=p1+geom_line(aes(y=0.05),colour="darkred",linetype="solid")
p1=p1+labs(x=element_blank(),y="P-value")
p1=p1+theme(axis.title.x=element_text(margin=margin(t=10),size=15),
            axis.title.y=element_text(margin=margin(r=10),size=20),
            plot.title=element_text(margin=margin(r=10,0,0,0),size=15,hjust=0.5),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            plot.margin=margin(t=0,r=10,b=0,l=10),
            axis.text.y=element_text(size=15),
            panel.grid.minor=element_blank())
p1=p1+scale_x_reverse(breaks=seq(3500,7500,500),limits=c(7500,3500),expand=expansion(mult=c(.025,.025)))
p1=p1+scale_y_continuous(breaks=seq(0,0.6,0.1),limits=c(0,0.6))


p2=ggplot(data=autores,aes(x=yearBP,y=size))
p2=p2+geom_line(colour="black",)
#p2=p2+geom_line(aes(y=0.05),colour="darkred",linetype="dotted")
p2=p2+labs(x="Years BP",y="Sample size")
p2=p2+theme(axis.title.x=element_text(margin=margin(t=10),size=25),
            axis.title.y=element_text(margin=margin(r=10),size=20),
            plot.title=element_text(margin=margin(r=10,0,0,0),size=20,hjust=0.5),
            plot.margin=margin(t=0,r=10,b=10,l=10),
            axis.text.x=element_text(size=15),
            axis.text.y=element_text(size=15),
            panel.grid.minor=element_blank())
p2=p2+scale_x_reverse(breaks=seq(3500,7500,500),limits=c(7500,3500),expand=expansion(mult=c(.025,.025)))
p2=p2+scale_y_continuous(breaks=seq(0,750,150),limits=c(0,750))


p3=p/p1/p2
p3=p3+plot_layout(heights=c(3,1,1))
print(p3)



yearBP=xchr[seq(1,41,1),2,1]-1500
sfa=xchr[seq(1,41,1),5,1]
xchrres=array(data=NA,dim=c(41,4,1))
xchrres=as.data.frame(xchrres)
colnames(xchrres)=c("yearBP","sfa","ec","size")
xchrres$yearBP=yearBP
xchrres$sfa=sfa
xchrres$size=xchr[seq(1,41,1),4,1]

for(i in 1:41) {
  xchrres[i,3]=sum(xchr[1,14,seq(1,1000,1)]>=xchr[i,7,1])/1000
}
head(xchrres)

range(xchrres$ec)




p=ggplot(data=xchrres,aes(x=yearBP,y=sfa))
p=p+geom_line(colour="black",)
#p=p+geom_line(aes(y=V10),colour="darkred",linetype="dashed")
#p=p+geom_vline(xintercept=3500,linetype="dotted",colour="darkred")
p=p+labs(y=expression(paste("Scaling Factor Angle ",alpha)),x=element_blank(),
         title="Mobility of the X chromosome in European ancient DNA")
p=p+theme(axis.title.x=element_text(margin=margin(t=10),size=15),
          axis.title.y=element_text(margin=margin(r=10),size=20),
          plot.title=element_text(margin=margin(r=10,0,10,0),size=35,hjust=0.5),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          plot.margin=margin(t=10,r=10,b=0,l=10),
          axis.text.y=element_text(size=15),
          panel.grid.minor=element_blank())
p=p+scale_x_reverse(breaks=seq(3500,7500,500),limits=c(7500,3500),expand=expansion(mult=c(.025,.025)))
p=p+scale_y_continuous(breaks=seq(0,60,10),limits=c(0,60))


p1=ggplot(data=xchrres,aes(x=yearBP,y=ec))
p1=p1+geom_line(colour="black",)
p1=p1+geom_line(aes(y=0.05),colour="darkred",linetype="solid")
p1=p1+labs(x=element_blank(),y="P-value")
p1=p1+theme(axis.title.x=element_text(margin=margin(t=10),size=15),
            axis.title.y=element_text(margin=margin(r=10),size=20),
            plot.title=element_text(margin=margin(r=10,0,0,0),size=15,hjust=0.5),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            plot.margin=margin(t=0,r=10,b=0,l=10),
            axis.text.y=element_text(size=15),
            panel.grid.minor=element_blank())
p1=p1+scale_x_reverse(breaks=seq(3500,7500,500),limits=c(7500,3500),expand=expansion(mult=c(.025,.025)))
p1=p1+scale_y_continuous(breaks=seq(0,0.6,0.1),limits=c(0,0.6))


p2=ggplot(data=xchrres,aes(x=yearBP,y=size))
p2=p2+geom_line(colour="black",)
#p2=p2+geom_line(aes(y=0.05),colour="darkred",linetype="dotted")
p2=p2+labs(x="Years BP",y="Sample size")
p2=p2+theme(axis.title.x=element_text(margin=margin(t=10),size=25),
            axis.title.y=element_text(margin=margin(r=10),size=20),
            plot.title=element_text(margin=margin(r=10,0,0,0),size=20,hjust=0.5),
            plot.margin=margin(t=0,r=10,b=10,l=10),
            axis.text.x=element_text(size=15),
            axis.text.y=element_text(size=15),
            panel.grid.minor=element_blank())
p2=p2+scale_x_reverse(breaks=seq(3500,7500,500),limits=c(7500,3500),expand=expansion(mult=c(.025,.025)))
p2=p2+scale_y_continuous(breaks=seq(0,750,150),limits=c(0,750))


p3=p/p1/p2
p3=p3+plot_layout(heights=c(3,1,1))
print(p3)

