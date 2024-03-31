library(surveyPrev)
library(rdhs)
library(tidyverse)
library(sf)

country <- 'Malawi'
survey_year <- 2015
indicator <- 'wasting'
#wasting is rare
# stunting is not

load(file = paste0('/Users/alanamcgovern/Desktop/Research/UN_Estimates/UN-Subnational-Estimates/Info/',paste0(country, "_general_info.Rdata"))) # load the country info
countryId <- dhs_countries()[dhs_countries()$ISO3_CountryCode==toupper(gadm.abbrev),]
potential_surveys <- dhs_datasets(countryIds = countryId$DHS_CountryCode, surveyYearStart = 2010, surveyType = 'DHS') %>% 
  dplyr::filter((FileType == 'Births Recode' & FileFormat=='Stata dataset (.dta)'))
survey_years <- as.numeric(potential_surveys$SurveyYear)

poly.path <- paste0("/Users/alanamcgovern/Desktop/Research/UN_Estimates/UN-Subnational-Estimates/Data/shapeFiles_gadm/gadm41_",gadm.abbrev,"_shp")
# other cases
if(!dir.exists(poly.path)){
  poly.path <- paste0("/Users/alanamcgovern/Desktop/Research/UN_Estimates/UN-Subnational-Estimates/Data/shapeFiles_alt/",country)
}

FulldhsData <- surveyPrev::getDHSdata(country = country, indicator = indicator, year = survey_year)
dhsdata <- getDHSindicator(FulldhsData, indicator = indicator)
geo <- getDHSgeo(country = country, year = survey_year)

MWIAdm1 <- read_sf(dsn = file.path(paste0(poly.path,'/',poly.layer.adm1,'.shp')))
MWIAdm1 <- as_Spatial(MWIAdm1)
MWIAdm1$NAME_1 <- MWIAdm1$DHSREGEN

MWIAdm2 <- read_sf(dsn = file.path(paste0(poly.path,'/',poly.layer.adm2,'.shp')))
MWIAdm2 <- as_Spatial(MWIAdm2)
MWIAdm2$NAME_2 <- MWIAdm2$NAME_1

admin1.data <- adminInfo(geo=MWIAdm1,admin=1)
admin1.data$admin.info$admin1 <- 1:nrow(admin1.data$admin.info)
admin2.data <- adminInfo(geo=MWIAdm2,admin=2)
admin2.data$admin.info$admin2 <- 1:nrow(admin2.data$admin.info)
colnames(admin2.data$admin.mat) <- rownames(admin2.data$admin.mat) <- admin2.data$admin.info$DistrictName

cluster.info <- clusterInfo(geo = geo, poly.adm1 = MWIAdm1, poly.adm2 =MWIAdm2)

#direct estimate
dir1 <- directEST(data = dhsdata, cluster.info = cluster.info, admin = 1)
dir1

#negbinomial estimate




