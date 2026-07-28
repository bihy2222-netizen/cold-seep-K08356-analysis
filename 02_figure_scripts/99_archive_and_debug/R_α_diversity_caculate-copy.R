#设置文件路径
setwd("/Users/catherine/Downloads/metaphlan 作图/class 级别作图-柱状图/") 
rm(list = ls())

#没有的话就安装这些安装包
install.packages("reshape2")
install.packages("dplyr")
install.packages("stringr")
install.packages("vegan")
install.packages("ggplot2")
install.packages("tidyr")
install.packages("crayon")


#导入需要的包
library(reshape2)
library(dplyr)
library(stringr)
library(vegan)
library(tidyr)
library(crayon)
library(ggplot2)


#读取物种数据和分组数据
class <- read.delim('class34.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
group <- read.delim("class——group.txt",header = T)

###计算多样性指数
#Shannon指数
Shannon <- diversity(order, index = "shannon", MARGIN = 2, base = exp(1))   #MARGIN：1用列数据算，2用行数据算
#simpson指数
Simpson <- diversity(order, index = "simpson", MARGIN = 2, base = exp(1))
#物种丰富度
Richness <- specnumber(order, MARGIN = 2) 

#将以上三个指数先统计成表
index <- as.data.frame(cbind(Shannon, Simpson, Richness))

#计算obs，chao，ace指数
tgenus<-ceiling(as.data.frame(t(order)))  #转置物种数据之后取整数（向上取整）
obs_chao_ace <- t(estimateR(tgenus))  #estimateR获取obs，chao，ace指数

#将obs，chao，ace指数与前面指数计算结果合并
index$Chao <- obs_chao_ace[,2]
index$Ace <- obs_chao_ace[,4]
index$obs <- obs_chao_ace[,1]

#计算Pielou指数
index$Pielou <- Shannon / log(Richness)  #这里使用自然对数ln，部分文献会使用log10和log2

#加一列sample以便后续匹配
index$sample<- c(rownames(index))

#合并分组信息与多样性指数
data <- merge(index,group,by = 'sample')

#导出alpha多样性指数
write.csv(data,"group1_alpha_data_order数据.csv")

###以Shannon指数为例画箱线图（即上文图一），需要其他指数替换Shannon即可
###多余颜色代码,'#0571b0','#0D99D4','#f4a582','#89CBB8','#FBAEA6','#FFCBAC','#81CCB6','#A3C39D','#F86A44','#99D1E6','#228EC7'
Shannon <- ggplot(data,aes(x=group,y=Shannon,color=group))+  
  stat_boxplot(geom = "errorbar", linewidth=2)+  #添加误差线
  geom_boxplot(linewidth=2)+  #箱线图
  scale_color_manual(values =c('#ca0020','#A39AF3','#F3C051','#92c5de','#0571b0','#0D99D4'))+ 
  geom_jitter(width = 0.1,alpha =0.5,size=2)+   #添加抖动点
  labs(title = "Shannon diversity index")+
  theme_bw()+
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),  
        plot.title = element_text(hjust = 0.5,size=30,face = "bold"),
        axis.title.x =element_blank(),
        axis.title.y = element_text(size=25,colour ="black",face = "bold"),
        axis.text = element_text(size=25,colour ="black",face = "bold"),
        legend.position = "none")

Shannon

#导出
ggsave("group2-Shannon-order.png",plot = Shannon,device = "png",height =9,width =16)

###所有指数画分面图（上文图二）
data_long <- tidyr::gather(data, variable, value, -sample, -group) #先转换数据框为长格式

#画图
total<-ggplot(data_long, aes(x = group, y = value,color=group)) +
  stat_boxplot(geom = "errorbar", linewidth=2)+
  geom_boxplot(linewidth=2)+
  scale_color_manual(values =c('#ca0020','#A39AF3','#F3C051','#92c5de','#0571b0','#0D99D4'))+
  geom_jitter(width = 0.1,alpha =0.5,size=2)+   
  facet_grid(variable~., scales = "free_y") + #按variable纵向分面
  labs(title = "alpha diversity", x = "", y = "")+
  theme_bw()+
  theme(panel.background = element_blank(),
        panel.grid = element_blank(), 
        strip.text.y = element_text(size = 15, colour = "black",face = "bold"), #分面标题文字设置
        plot.title = element_text(hjust = 0.5,size=30,face = "bold"),
        axis.text.x = element_text(size=15,colour ="black",face = "bold"),
        axis.text.y = element_text(size=10,colour ="black"),
        legend.position = "none")

total

#导出
ggsave("group4-total_alpha_class.png",plot = total,device = "png",height =16,width =9)



