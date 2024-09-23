# Analysis of bacterial genomes evolved in the absence of phage
# Reena Debray
# Feb 1, 2022
# Modified: Evelia Coss (03/Sep/2024)

#TODO: Este parrafo no se entiende, no encuentro la función fit_growthmodel() en otras partes del codigo. 

## This function uses the fit_growthmodel() function of the growthrates package to calculate growth rates from microplate data
## It takes data in the following form: each row is the measurement of a single well at a single time point, with additional columns "Treatment" and "OD".
## A logistic curve is fit to each well separately using the OD over time
## I used the following parameters: lower=c(y0=0.000001,mumax=0,K=0),upper=c(y0=0.05,mumax=5,K=1.5)),p=c(y0=0.01,mumax=0.2,K=0.1)

#TODO: Creo que faltaría la descripción del script.

#TODO: Una sugerencia sería crear un proyecto de R para que se ajusten automaticamente el working directory
#TODO: También, en el espiritu de la reproducibilidad, hay un paquete que se llama renv que guarda un ambiente con los paquetes que se utilizan dentro del proyecto (https://rstudio.github.io/renv/articles/renv.html). Segun entiendo empiezas con un ambiente vacío y jala los paquetes que ya tienes instalados al ambiente pero guarda todas las versiones (incluso me parece que la de R también).

# To enhance reproducibility we set up an environment for the following project.
# --- set environment ---
library(renv)

# This command initializes the R environment it should only be used once. 
# renv::init()

# To ensure you have the proper versions and packages needed for this project run the following command.
renv::restore()


# --- packages ----
library(growthrates)
library(lmerTest)
library(readxl)
library(tidyverse)
library(reshape2)


#---- GC_long_to_growthrates() -------
#' Calculate Growth Rates from Microplate Data
#'
#' This function calculates growth rates from microplate data using a logistic growth model.
#' It fits a logistic curve to each well's optical density (OD) over time and extracts the growth rate.
#'
#' @param GC_long A data frame containing microplate data. It must have columns `Well`, `Time`, `OD`, and `Treatment`.
#' @param lower A numeric vector specifying the lower bounds for the parameters of the logistic growth model.
#' @param upper A numeric vector specifying the upper bounds for the parameters of the logistic growth model.
#' @param p A numeric vector specifying the initial parameter values for the logistic growth model.
#'
#' @return A data frame with columns:
#' \item{well}{The identifier of the well.}
#' \item{treatment}{The treatment associated with the well.}
#' \item{r}{The growth rate obtained from the fitted logistic model.}
#'
#' @import growthrates
#' @import lmerTest
#'
#' @export
#'
#' @examples
#' # Example data
#' GC_long <- data.frame(
#'   Well = rep(c("A1", "A2"), each = 10),
#'   Time = rep(1:10, 2),
#'   OD = runif(20),
#'   Treatment = rep(c("Control", "Treatment"), each = 10)
#' )
#'
#' # Parameters for the logistic model
#' lower <- c(y0 = 0.000001, mumax = 0, K = 0)
#' upper <- c(y0 = 0.05, mumax = 5, K = 1.5)
#' p <- c(y0 = 0.01, mumax = 0.2, K = 0.1)
#'
#' # Calculate growth rates
#' growth_rates <- GC_long_to_growthrates(GC_long, lower, upper, p)
#' print(growth_rates)

GC_long_to_growthrates<-function(GC_long,lower,upper,p){
  ###initialize data frame
  growthrates<-data.frame(matrix(nrow=0,ncol=3))
  #populate with model fit
  for (well in unique(GC_long$Well)){
    treatment<-GC_long[GC_long$Well==well,"Treatment"][1]
    r<-as.numeric(coef(fit_growthmodel(FUN = grow_logistic, GC_long[GC_long$Well==well,"Time"], GC_long[GC_long$Well==well,"OD"],p=p,lower=lower,upper=upper))[2])
    growthrates<-rbind(growthrates,c(well,treatment,r))
  }
  ###return output
  colnames(growthrates)=c("well","treatment","r")
  growthrates$r<-as.numeric(growthrates$r)
  return(growthrates)
}

## ----- Figure 1: Resistance is costly -------

### Read in "Costs_of_Res.xlsx" as costs_of_res
costs_of_res <- read_excel("Costs_of_Res.xlsx")


### Calculate fitted values controlling for phage resistance and plate layout
### Express fitted values as a percentage of wild-type fitness
costs_of_res$fitted<-fitted.values(lm(r~(Population=="ancDC3000")+Column+Plate,costs_of_res))
costs_of_res$fitted_percWT<-costs_of_res$fitted/mean(unlist(costs_of_res[costs_of_res$Population=="ancDC3000","fitted"]))*100

#TODO: ¿Esta bien que algunos genes tengan valores faltantes?
table(costs_of_res$Gene)


### Reorder the isolates by resistance gene and growth rate
costs_of_res<-costs_of_res[order(costs_of_res$Gene,-costs_of_res$r),]

costs_of_res[costs_of_res$Population=="MS15","Gene"]="Z" # place the isolate with no detected genetic differences on the right-hand side of the graph
costs_of_res[costs_of_res$Population=="MS15","Annotation"]="Z" # place the isolate with no detected genetic differences on the right-hand side of the graph

costs_of_res$order<-seq(1,nrow(costs_of_res))

### Resistant bacteria grow more slowly than their sensitive ancestor
costs_of_res_R<-costs_of_res[costs_of_res$Population!="ancDC3000",]

costs_agg<-aggregate(costs_of_res_R$fitted_percWT,by=list(costs_of_res_R$Population,costs_of_res_R$Gene),FUN=mean)

colnames(costs_agg)<-c("Population","Gene","fitted_percWT")
t.test(x = costs_agg$fitted_percWT,mu = 100,alternative = "less")

### Variation in growth rates is not explained by resistance gene (exclude population with no detected mutations)
anova(lm(fitted_percWT~Gene,costs_agg[costs_agg$Population!="MS15",]))

### Code for Figure 1B

#TODO: Los titulos no estan bien puestos. Les faltaba color y quitar el último.

library(ggh4x) # para ponerle color a los titulos
library(RColorBrewer)

# Asignar colores a los titulos de los facets
strip <- strip_themed(text_x = elem_list_text(color = brewer.pal(5, "Dark2")))


#TODO: faltan los titulos y los ticks del eje x en un angulo de 90.

# Guardar las etiquetas que van en los títulos de los facets
strip_label <- distinct(costs_of_res[,c("Gene","Annotation")])[-c(1:2),] %>% mutate(label = paste(Gene,Annotation, sep = "\n"))


costs_of_res[costs_of_res$Population!="ancDC3000",] %>% 
  mutate(Population = ifelse(Population == "MS15","no detected\nmutation",Population),
         Gene = ifelse(Gene == "Z","",Gene),
         Annotation = ifelse(Annotation == "Z","",Annotation),
         label = paste(Gene,Annotation, sep = "\n")) %>% 
  mutate(label = factor(label, levels = c(strip_label$label,"\n"))) %>% 
ggplot()+
  stat_summary(aes(x = reorder(Population,-fitted_percWT),
                   y = fitted_percWT,
                   group=Population,
                   fill=Gene),
               geom="pointrange",
               shape=21,
               size=0.8)+
  geom_hline(yintercept=100,linetype="dashed")+
  scale_fill_brewer(palette="Dark2")+
  guides(fill=F)+
  ylab("Fitness (% of wild-type)")+
  facet_wrap2(~label,
             scales="free_x",
             #space="free_x", 
             strip.position = "bottom",
             nrow = 1,
             strip = strip
             )+
  theme_classic(base_size=18)+
  theme(axis.title.x=element_blank(),
                axis.text.x=element_text(angle = 90),
                #axis.ticks.x=element_blank(),
        strip.background = element_blank(),
        strip.placement = "outside",
        strip.text = element_text(size = 12))


## ------- Figure 2A: Growth over time ---------

### Read in "Fitness_over_Time.xlsx" as fitness_over_time
fitness_over_time <- read_excel("Fitness_over_Time.xlsx")

### Calculate fitted values controlling for phage resistance, passage, and plate layout
fitness_over_time$fitted<-fitted.values(lm(r~Passage*Type+Column,fitness_over_time))

### Express fitted values as a percentage of wild-type fitness

for (passage in seq(0,12,2)){
  anc_means<-mean(unlist(fitness_over_time[fitness_over_time$Passage==passage & fitness_over_time$Type=="ANC","fitted"]))
  fitness_over_time[fitness_over_time$Passage==passage & fitness_over_time$Type=="PR","fitted_percWT"]<-100*fitness_over_time[fitness_over_time$Passage==passage & fitness_over_time$Type=="PR","fitted"]/anc_means
}

### Code for Figure 2A

df<-data.frame(
  aggregate(unlist(fitness_over_time[fitness_over_time$Type=="PR","fitted_percWT"]),
            by=list(unlist(fitness_over_time[fitness_over_time$Type=="PR","Passage"])),FUN=mean),
  aggregate(unlist(fitness_over_time[fitness_over_time$Type=="PR","fitted_percWT"]),
            by=list(unlist(fitness_over_time[fitness_over_time$Type=="PR","Passage"])),FUN=sd)[,2])

colnames(df)=c("passage","mean","SD")

#TODO: Faltaban las lineas de atrás

ggplot()+
  geom_point(data=df,aes(x=passage*3,y=mean),size=3)+
  geom_errorbar(data=df,aes(x=passage*3,ymin=mean-SD,ymax=mean+SD),width=0,linewidth=1)+
  geom_line(data=df,aes(x=passage*3,y=mean))+
  geom_line(data = fitness_over_time, aes(x = Passage*3, y = fitted_percWT, group = Population), linewidth = 0.2, alpha = 0.4)+
  theme_classic(base_size=16)+
  geom_hline(yintercept=100,linetype="dashed")+
  scale_x_continuous(breaks=seq(0,36,6))+
  ylab("Fitness (% of wild-type)")+
  xlab("Day of experimental evolution")+
  coord_cartesian(ylim=c(60,115))


## ------ Figure 2B: Resistance over time -------

#TODO: Read in "Resistance_over_Time.xlsx" as res_over_time
res_over_time <- read_excel("Resistance_over_Time.xlsx")
res_long<-melt(res_over_time,id.vars = "isolate",measure.vars = c("Prop_S_P0","Prop_S_P4","Prop_S_P8","Prop_S_P12"))
res_long$variable<-as.character(res_long$variable)
res_long$passage<-as.numeric(substr(res_long$variable,9,nchar(res_long$variable)))

### Increase in fitness over time; no main or interaction effect of populations that reverted to sensitivity
fitness_over_time[fitness_over_time$Population%in%c("SNK7","QAC5","VCM4"),"Outcome"]<-"S"
fitness_over_time[!fitness_over_time$Population%in%c("SNK7","QAC5","VCM4"),"Outcome"]<-"R"
anova(lmer(fitted_percWT~Passage*Outcome+(1|Population),fitness_over_time))

### No effect of costs of resistance on reversion to sensitivity
costs_agg[costs_agg$Population%in%c("SNK7","QAC5","VCM4"),"Outcome"]<-"S"
costs_agg[!costs_agg$Population%in%c("SNK7","QAC5","VCM4"),"Outcome"]<-"R"
t.test(fitted_percWT~Outcome,costs_agg)


### ----- Code for Figure 2B -----
# It could have been an additional transformation of res_over_time or a similar dataset.
LOR_res_long <- melt(res_over_time, id.vars = "isolate", measure.vars = c("Prop_S_P0","Prop_S_P4","Prop_S_P8","Prop_S_P12"))
LOR_res_long$variable <- as.character(LOR_res_long$variable)
LOR_res_long$passage <- as.numeric(substr(LOR_res_long$variable,9,nchar(LOR_res_long$variable)))
LOR_res_long <- LOR_res_long[!LOR_res_long$isolate %in% c("FMS6","SNK6"),]  # Insulation filtering

### Excludes populations with colonies that appeared sensitive on plates but were not affected by phage in growth curves
ggplot(LOR_res_long[!LOR_res_long$isolate%in%c("FMS6","SNK6"),],
       aes(x = passage*3,
           y = value*100,
           group=isolate,
           color=(isolate%in%c("SNK7","QAC5","VCM4"))))+
  geom_point(size=4.5,shape=18)+geom_line(size=0.8)+
  theme_classic()+
  guides(color=F)+
  scale_color_brewer(palette="Dark2")+
  scale_x_continuous(breaks=c(0,12,24,36))+
  theme_classic(base_size=16)+ylab("Proportion of phage sensitivity in population (%)")+
  xlab("Day of experimental evolution")+
  theme(axis.title.y=element_text(size=14))

## ----- Figure 3: Replay experiment -------

### Read in "Replay.xlsx" as replay
replay <- read_excel("Replay.xlsx", sheet = "Data")
replay$Prop_S_P12 <- as.numeric(as.character(replay$Prop_S_P12)) # Convert as numeric
### Populations derived from a founder that re-evolved sensitivity were more likely to re-evolve sensitivity as well
t.test(Prop_S_P12~Founder_outcome,data=replay)

### Code for Figure 3C
agg<-data.frame(aggregate(100*replay$Prop_S_P12,
                          by=list(replay$Founder_outcome),
                          FUN=function(x){mean(x,na.rm=T)}),
                aggregate(100*replay$Prop_S_P12,
                          by=list(replay$Founder_outcome),
                          FUN=function(x){sd(x,na.rm=T)/sqrt(length(x))}))

colnames(agg)=c("Founder_outcome","mean","Founder_2","SE")

ggplot()+
  geom_errorbar(data=agg,
                aes(Founder_outcome,ymax=mean+SE,ymin=mean-SE),
                alpha=0.4,
                width=0.6,
                size=0.8)+
  stat_summary(data=agg,
               aes(Founder_outcome,mean,fill=Founder_outcome),
               geom="bar",
               fun=mean,
               width=0.8,
               color="black")+
  geom_jitter(data=replay,
              aes(Founder_outcome,Prop_S_P12*100),
              size=3,
              alpha=0.4,
              width=0.05,
              height=1)+
  theme_classic(base_size=16)+
  scale_fill_brewer(palette="Dark2")+
  scale_x_discrete(labels = c("R" = "Remained\nphage-resistant","S"="Re-evolved\nphage sensitivity"))+
  guides(fill=F)+
  xlab("Outcome of founder population")+
  ylab("Proportion of phage sensitivity in replay populations (%)")


### supplement

#diminishing returns
cast<-dcast(fitness_over_time[fitness_over_time$Sample!="ANC1_P0",],Population~Passage,FUN=mean,value.var = "fitted_percWT")
colnames(cast)=c("Population","P0","P2","P4","P6","P8","P10","P12")
anova(lm((P12-P0)~P0,cast))
ggplot(cast,aes(P0,P12-P0))+stat_smooth(method="lm",color="gray30")+geom_jitter(width=0.8,height=0.8,size=2,alpha=0.8)+theme_classic(base_size=16)+xlab("Initial fitness\n(% wild-type, day 0)")+ylab("Change in fitness\n(% wild-type, day 0 - day 36)")
