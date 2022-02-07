# read 
n=lst_2022_01_13$n
V=lst_2022_01_13$V
f=lst_2022_01_13$f
WEuclid=lst_2022_01_13$WEuclid
Wgraph=lst_2022_01_13$Wgraph
Whst=lst_2022_01_13$Whst
W=lst_2022_01_13$W

# 
library(tidyverse)
library(devtools)
source_url("https://raw.githubusercontent.com/miruetoto/yechan/master/_notebooks/heavysnow.R")

# 유클리드와 그래프도메인의 정보를 플랏.
df<-tibble(V)
df$f<-as.vector(f) 
df$divlink<-0
df$conlink<-0
df$link<-0
for(i in 1:23){
  df$divlink[i]<-sum(W[i,])
  df$conlink[i]<-sum(W[,i])
  df$link[i]<-(sum(W[i,])+sum(W[,i]))/2
}
df$type<-rep('solo',23)
df[V %in% c("The Avengers" ,"Avengers: Age of Ultron","Captain America: Civil War","Avengers: Infinity War","Avengers: Endgame"),'type']='multi'
df

# 
p1<-ggplot(df, aes(x=conlink, y=divlink,z=f))+
  geom_point(aes(x=conlink, y=divlink,size=f**2,shape=type))+
  scale_size(range=c(2,20))+
  scale_alpha_manual(values=c(1,0.6))+
  scale_shape_manual(values=c(9,1))+
  theme_classic()+
  guides(size="none")+
  geom_text_repel(aes(x=conlink,y=divlink,label=V),col="gray40")+
  #geom_density2d(aes(x=conlink,y=divlink,fill=tatf,color=f)),+
  ylab("Out-degree")+xlab("In-degree")+
  ggtitle("")+theme(plot.title=element_text(face="bold.italic",size=rel(2)))+
  theme(axis.title.x=element_text(face=4,size=rel(1)))+theme(axis.title.y=element_text(face=4,size=rel(1)))
#theme(legend.position = "none")
p1
#ggsave(plot=p1,"./2022-01-13-p1.pdf",width=10,height=7)

# 
g1<-gfft(f,WEuclid)
g2<-gfft(f,Wgraph)
g3<-gfft(f,Whst)


#
library(gridExtra)
e1<-eigenplot(g1)+theme(plot.title=element_text(face="bold.italic",size=rel(1)))
e2<-eigenplot(g2)+theme(plot.title=element_text(face="bold.italic",size=rel(1)))
e3<-eigenplot(g3)+theme(plot.title=element_text(face="bold.italic",size=rel(1)))
s1<-specplot(g1)
s2<-specplot(g2)
s3<-specplot(g3)
p2<-grid.arrange(e1,e2,e3,s1,s2,s3,nrow=2)
#ggsave(plot=p2,"./2022-01-13-p2.pdf",width=5,height=3)

#
d1<-decompose(f,WEuclid,V=V) # 0, 35000, 60000, 80000
d2<-decompose(f,Wgraph,V=V) # 0, 35000, 60000, 80000
d3<-decompose(f,Whst,V=V) # 0, 35000, 60000, 80000

d1$case<-"Euclid"
d2$case<-"Graph"
d3$case<-"HST"
df2_ <-rbind(d1,d2,d3)
df2_ %>% group_by(case,eigenvectorindex) #%>% mutate(textsize= 10*(abs(fhat)>70))
df2=merge(df, df2_ %>% group_by(case,eigenvectorindex) %>% mutate(textsize= 10*(abs(fhat)>70))) %>% as_tibble
complst_ = c() 
for(i in 1:23) complst_[i] = sum(filter(df2,eigenvectorindex==i,case=='HST')$fhat**2)/sum(filter(df2,eigenvectorindex==i,case=='HST')$f**2)
print(which(complst_ > 0.01))
p3<-ggplot(data=filter(df2,eigenvectorindex %in% which(complst_ > 0.01), case=='HST'),aes(x=V,y=fhat))+
  geom_col(aes(fill=fhat>0),width=0.7)+
  #geom_text_repel(aes(x=Vindex,y=fhat,label=V,size=textsize),col=1,fontface=4,alpha=0.8,segment.size=0.2,segment.color="gray60",min.segment.length=5,hjust=0.1)+
  scale_radius(range = c(0,1.8))+
  guides(size='none')+
  facet_wrap("eigenvectorindex",ncol=1,scales='free')+
  geom_hline(aes(yintercept=0),col="gray60",lty=2)+
  xlab("")+ylab("")+guides(fill='none')+
  theme_classic()+
  #theme(strip.text.x = element_text(angle=45,size = 10, color = "black", face = "bold.italic"))+
  #theme(strip.text.y = element_text(size = 10, color = "black", face = "bold.italic"))+
  theme(axis.text.x=element_text(angle=35, hjust=1, vjust=1,size=rel(0.8)))+
  theme(plot.title=element_text(face="bold.italic",size=rel(1.5)))+
  scale_x_discrete(limits=V)

df2_ = df2 %>% mutate(fhatsquare= fhat**2) 
df3_ = df2_ %>% group_by(eigenvectorindex,case) %>% summarise(fhatsquaresum= sum(fhatsquare))
df3_ = merge(df2_,df3_,key=eigenvectorindex) %>% as_tibble
df3 = df3_ %>% mutate(fhatprop= fhatsquare/fhatsquaresum)
df3


## 1. define friendship 
library(igraph)
frnd_ship<-friendship(W)

## 2. define relations
relations<-expand.grid(from=V, to=V)
relations<-cbind(relations,frnd_ship)

## 3. make gdf and weight
gdf<-graph_from_data_frame(relations[frnd_ship>0,],directed=TRUE,vertices=V)
wght<-frnd_ship[frnd_ship>0]
library(ggraph)
p5<-ggraph(gdf)+theme_void()+
  geom_edge_arc(curvature=0.1,aes(alpha=frnd_ship),edge_width=0.2,arrow=arrow(length=unit(2,'mm'),type='closed',angle=10),linewidth=20)+
  geom_node_point(size=3+f/100,alpha=0.7,color=abs(f)>500)+
  #geom_node_text(aes(label = V,size=f), color="black", repel=T,fontface=4)+
  guides(size=FALSE) + guides(edge_alpha=FALSE) + 
  ggtitle("")+ theme(plot.title=element_text(hjust=0.09))+theme(plot.title=element_text(face="bold.italic",size=rel(1.5)))
p5
