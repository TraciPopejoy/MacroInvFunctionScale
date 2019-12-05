# 3_statistical_analyses
library(MuMIn);library(quantreg);library(ggplot2)

# read in binned size spectra data
binned <- read.csv("binned_size_spectra.csv",
                   stringsAsFactors = FALSE)
# make a data frame with the appropriate info: Mussel biomass, year, season, huc12
SizeRegData<-binned %>% left_join(fieldSiteKey, by="Sample") %>% 
  left_join(FieldSpData@data) %>% 
  dplyr::select(-SamplingSeason,-Quadrat.n, -Reach,-`USGS gage Identifier`,
                -SiteID,-HUC12,-HUC8)

# global models ####
# full quadratic
global.quadratic <- lm(log_count_corrected~ 
                         log_mids_center *Mussel.bm.g +
                         I(log_mids_center^2) + 
                         I(log_mids_center^2):Mussel.bm.g,
                       data=SizeRegData,
                       na.action = "na.fail")
# full linear model
global.linear <- lm(log_count_corrected ~ 
                      log_mids_center*Mussel.bm.g,
                    data=SizeRegData,
                    na.action = "na.fail")
# compare quadratic and linear models
AIC(global.quadratic, global.linear)
# quadratic much better than linear
# quad AIC = 4265.491
# linear AIC = 5017.710
# move forward with quad model

# systematically test simplified quadratic models using MuMIn::dredge 
dredge.models <- dredge(global.quadratic,
                        beta = "none",
                        extra = "R^2")

# table 2 for MS #### 
# table with results for all simplified models
sink("table_2.txt")
dredge.models
sink()

# pick best model based on AICc using MuMIn::get.models
top.models <- get.models(dredge.models,
                         subset = delta < 2)
# single top model selected

# save top.model for making figure in script 4
saveRDS(top.models[[1]],
        "top_quadratic_model.RDS")


# Mrange ####
# Mrange ####
# calculate range of M values
mrange <- SizeRegData %>% 
  group_by(site, Mussel.bm.g) %>%
  summarise(mrange = max(log_mids_center) - min(log_mids_center))
# save Mrange as csv
write.csv(mrange, "results/mrange_data.CSV",
          row.names = FALSE)

# mrange ~ gradient linear model
m.mod <- lm(mrange ~ pca1, data = mrange)
summary(m.mod)

# Quantile regression
M.quant <- rq(log_mids~pca1, data = binned, tau = c(0.05, 0.95))
summary(M.quant)


