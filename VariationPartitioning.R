# Variation partitioning


fallVar<-cca(field.com[,-c(7,13,22)]~
               `Average of STDM (g.m-2)`+
               Benthic_CHLA_MG.M2+Discharge.cms+Dvar+
               HUC12num, data=field.env)
fallVarR<-RsquareAdj(fallVar)$r.squared

#Field Data
showvarparts(3, bg=2:6)
FVar<-varpart(field.com[,-c(7,13,22)], 
              ~`Average of STDM (g.m-2)`, 
              ~Benthic_CHLA_MG.M2+Discharge.cms+Dvar,
              ~HUC12num,
              data=field.env, chisquare = T, 
              permutations=1000)
FVar
#fmusVar = a+f+d+g or a+al+ae+ael
fmusVar<-cca(field.com[,-c(7,13,22)]~
               `Average of STDM (g.m-2)`, data=field.env)
fmusVarR<-RsquareAdj(fmusVar)$r.squared
#fenvVar = b+d+g+e or b+ae+el+ael
fenvVar<-cca(field.com[,-c(7,13,22)]~
               Benthic_CHLA_MG.M2+Discharge.cms+Dvar, data=field.env)
fenvVarR<-RsquareAdj(fenvVar)$r.squared
#flocVar = c+f+g+e or l+al+el+ael
flocVar<-cca(field.com[,-c(7,13,22)]~
               HUC12num, data=field.env)
flocVarR<-RsquareAdj(flocVar)$r.squared

#fonlymus = a
fonlymus<-cca(field.com[,-c(7,13,22)]~
               `Average of STDM (g.m-2)`+
               Condition(Benthic_CHLA_MG.M2+Discharge.cms+Dvar+
               HUC12num), data=field.env)
RsquareAdj(fonlymus)
#fonlyenv= b or e
fonlyenv<-cca(field.com[,-c(7,13,22)]~
               Condition(`Average of STDM (g.m-2)`+HUC12num)+
               Benthic_CHLA_MG.M2+Discharge.cms+Dvar, data=field.env)
RsquareAdj(fonlyenv)
#fonlyloc = c or l
fonlyloc<-cca(field.com[,-c(7,13,22)]~
                Condition(`Average of STDM (g.m-2)`+
                Benthic_CHLA_MG.M2+Discharge.cms+Dvar)+
                HUC12num, data=field.env)
RsquareAdj(fonlyloc)

#fae.l = a+b+d or a+e+ae
fae.l<-cca(field.com[,-c(7,13,22)]~
             `Average of STDM (g.m-2)`+
             Benthic_CHLA_MG.M2+Discharge.cms+Dvar+
             Condition(HUC12num), data=field.env)
#fal.e = a+c+f or a+l+al
fal.e<-cca(field.com[,-c(7,13,22)]~
             `Average of STDM (g.m-2)`+ HUC12num +
             Condition(Benthic_CHLA_MG.M2+Discharge.cms+Dvar),
           data=field.env)
#fel.a = c+e+b or e+l+el
fel.a<-cca(field.com[,-c(7,13,22)]~
             Condition(`Average of STDM (g.m-2)`)+
             Benthic_CHLA_MG.M2+Discharge.cms+Dvar+
             HUC12num, data=field.env)
a<-RsquareAdj(fonlymus)$r.squared #a
b<-RsquareAdj(fonlyenv)$r.squared #b
c<-RsquareAdj(fonlyloc)$r.squared #c
d<-RsquareAdj(fae.l)$r.squared - RsquareAdj(fonlymus)$r.squared - 
  RsquareAdj(fonlyenv)$r.squared #d 
f<-RsquareAdj(fal.e)$r.squared - RsquareAdj(fonlymus)$r.squared - 
  RsquareAdj(fonlyloc)$r.squared #f 
e<-RsquareAdj(fel.a)$r.squared - RsquareAdj(fonlyenv)$r.squared - 
  RsquareAdj(fonlyloc)$r.squared #e 
g<-RsquareAdj(fallVar)$r.squared-a-b-c-d-e-f
FVar$part$indfract
data.frame(FVar$part$indfract[3],man=c(a,b,c,d,e,f,g,NA))
sum(a,b,c,d,e,f,g)
RsquareAdj(fallVar)


library(venneuler)
fven <- venneuler(c(Mussels=a*100, Environment=b*100,
                    Location=c*100,
                    "Mussels&Environment"=d*100,
                    "Environment&Location"=e*100,
                    "Mussels&Location"=f*100,
                    "Mussels&Location&Environment"=g))
plot(fven)
library(eulerr)
fvarr2<-c(c(Mussels=round(a*100,1), 
        Environment=round(b*100,1),
        Location=round(c*100,1),
        "Mussels&Environment"=round(d*100,1),
        "Environment&Location"=round(e*100,1),
        "Mussels&Location"=round(f*100,1),
        "Mussels&Location&Environment"=0))
plot(venn(fvarr2), quantities=T)

# Enclosures

EVar
eb<-EVar$part$fract$R.squared[1]+EVar$part$fract$R.squared[2]-
  EVar$part$fract$R.squared[3]
ea<-EVar$part$fract$R.squared[1]-eb
ec<-EVar$part$fract$R.squared[2]-eb
evarr2<-c(Mussels=round(ea*100,1), 
            Environment=round(ec*100,1),
            "Mussels&Environment"=round(eb*100,1))
plot(venn(evarr2), quantities=T)


EVar;plot(EVar)
even <- venneuler(c(Mussels=0.9, Environment=0.2, 
                    "Mussels&Environment"=0.3))
plot(even)

SVar

sb<-SVar$part$fract$R.squared[1]+SVar$part$fract$R.squared[2]-
  SVar$part$fract$R.squared[3]
sa<-SVar$part$fract$R.squared[1]-sb
sc<-SVar$part$fract$R.squared[2]-sb
sum(SVar$part$fract$R.squared[1:3]) #~0.1447
svarr2<-c(Mussels=round(sa*100,1), 
          Environment=round(sc*100,1),
          "Mussels&Environment"=round(sb*100,1))
plot(venn(svarr2), quantities=T)
sven <- venneuler(c(Mussels=3.1, Chlorophylla=1.4, 
                    "Mussels&Chlorophylla"=0.8))
plot(sven)


CCAres<-data.frame(scale=rep(c("Field","Enc","Shell"),each=3),
                   Type=c("Location","Explained","Unexplained"),
                   Value=c(6.276, 35.050, 71.226, 
                           NA, 2.512, 97.488,
                           NA, 5.151, 94.849))

ggplot(CCAres[CCAres$scale=="Shell",],
       aes(x = "", y = Value, fill = Type,label=Type)) +
  geom_bar(stat = "identity", color = "white") +
  geom_text_repel()+
  coord_polar("y", start = 0)+theme_void()+
  scale_fill_grey()+
  theme(legend.position = "none")



