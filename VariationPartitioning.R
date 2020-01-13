# Variation partitioning

#Field Data
showvarparts(3, bg=2:6)
FVar<-varpart(field.com[,-c(7,13,22)], 
              ~`Average of STDM (g.m-2)`, 
              ~Benthic_CHLA_MG.M2+Discharge.cms+Dvar,
              ~HUC12num,
              data=field.env, chisquare = T, 
              permutations=1000)
FVar
fmusVar<-cca(field.com[,-c(7,13,22)]~
               `Average of STDM (g.m-2)`, data=field.env)
fmusVarR<-RsquareAdj(fmusVar)$r.squared
fenvVar<-cca(field.com[,-c(7,13,22)]~
               Benthic_CHLA_MG.M2+Discharge.cms+Dvar, data=field.env)
fenvVarR<-RsquareAdj(fenvVar)$r.squared
flocVar<-cca(field.com[,-c(7,13,22)]~
               HUC12num, data=field.env)
flocVarR<-RsquareAdj(flocVar)$r.squared
fallVar<-cca(field.com[,-c(7,13,22)]~
               `Average of STDM (g.m-2)`+
               Benthic_CHLA_MG.M2+Discharge.cms+Dvar+
               HUC12num, data=field.env)
fallVarR<-RsquareAdj(fallVar)$r.squared

fallVar<-cca(field.com[,-c(7,13,22)]~
               `Average of STDM (g.m-2)`+
               Benthic_CHLA_MG.M2+Discharge.cms+Dvar+
               HUC12num, data=field.env)
fonlymus<-cca(field.com[,-c(7,13,22)]~
               `Average of STDM (g.m-2)`+
               Condition(Benthic_CHLA_MG.M2+Discharge.cms+Dvar+
               HUC12num), data=field.env)
RsquareAdj(fonlymus)
fonlyenv<-cca(field.com[,-c(7,13,22)]~
               Condition(`Average of STDM (g.m-2)`+HUC12num)+
               Benthic_CHLA_MG.M2+Discharge.cms+Dvar, data=field.env)
RsquareAdj(fonlyenv)
fonlyloc<-cca(field.com[,-c(7,13,22)]~
                Condition(`Average of STDM (g.m-2)`+
                Benthic_CHLA_MG.M2+Discharge.cms+Dvar)+
                HUC12num, data=field.env)
RsquareAdj(fonlyloc)
fallVar<-cca(field.com[,-c(7,13,22)]~
               `Average of STDM (g.m-2)`+
               Benthic_CHLA_MG.M2+Discharge.cms+Dvar+
               HUC12num, data=field.env)
