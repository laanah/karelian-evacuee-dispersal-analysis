library(dplyr)
library(gamlss)
library(MASS)
library(Rcpp)
library(lme4)
library(MuMIn)
library(lsmeans)

# Data preparation
######################################################################

## Filter data
######################################################################

data <- readRDS("~/person_data.rds")
data <- data %>% filter(birthregion == "karelia" & primaryperson == 1)
data <- data[ which(data$birthyear < 1925), ]


## Generate initial location in Karelia before evacuation
######################################################################

living <- readRDS("~/living_record.rds")
drive <- read.csv("~/location_sheet.csv")

newdata <- left_join(drive, living, by = c("placename" = "name"))
newdata$fromkarelia <- with(newdata, ifelse(movedOut < 41 & movedOut > 38 & region.x == "karelia", 1, 0))

fromK <- newdata[ which(newdata$fromkarelia == 1), ]
fromK <- fromK %>% dplyr::select("placename", "Population", "Latitude", "Longitude", "personId", "movedOut")
fromK <- fromK[ order(fromK$movedOut), ]
fromK <- fromK[ !duplicated (fromK$personId), ]


## Rename columns
######################################################################

data <- left_join(data, fromK, by = c("id" = "personId"))
names(data)[ names(data)== "Population" ] <- "frompop"
names(data)[ names(data)== "Latitude" ] <- "fromlat"
names(data)[ names(data)== "Longitude" ] <- "fromlon"
names(data)[ names(data)== "movedOut" ] <- "fromyear"
names(data)[ names(data)== "placeId" ] <- "fromplace"


## Log populations and scale them to floats between 0 and 1
######################################################################

data$Birth_pop_log <- log(data$birthpopulation)
data$from_pop_log <- log(data$frompop)
data$FDF_pop_log <- log(data$fdf_population)

data$FDF_pop_log <- data$FDF_pop_log - min(data$FDF_pop_log, na.rm=TRUE)
data$FDF_pop_log <- data$FDF_pop_log / max(data$FDF_pop_log, na.rm=TRUE)

data$fdf_latitude1 <- data$fdf_latitude - min(data$fdf_latitude, na.rm=TRUE)
data$fdf_latitude1 <- data$fdf_latitude / max(data$fdf_latitude, na.rm=TRUE)

data$fdf_longitude1 <- data$fdf_longitude - min(data$fdf_longitude, na.rm=TRUE)
data$fdf_longitude1 <- data$fdf_longitude / max(data$fdf_longitude, na.rm=TRUE)

data$Birth_pop_log <- data$Birth_pop_log - min(data$Birth_pop_log, na.rm=TRUE)
data$Birth_pop_log <- data$Birth_pop_log / max(data$Birth_pop_log, na.rm=TRUE)

# Prepare evacuee data for models and maps

all_evacuees <- as.data.frame(table(data$fdf_name))
names(all_evacuees)[ names(all_evacuees) == "Var1" ] <- "fdf_name"
names(all_evacuees)[ names(all_evacuees) == "Freq" ] <- "freqAll"


# Models
######################################################################

## FDF Population model
######################################################################

first_dest_finland_model <- data %>% dplyr::select("sex", "age_1970",
                                                   "fdf_latitude1", "fdf_longitude1",
                                                   "FDF_pop_log", "returnedkarelia")
first_dest_fin_complete <- first_dest_finland_model[ complete.cases(first_dest_finland_model), ]

FDFpopmodel <- glm(
  returnedkarelia ~ sex + age_1970 + fdf_latitude1 + fdf_longitude1 + FDF_pop_log,
  data=first_dest_finland_model, family = binomial(link = "logit")
)

summary(FDFpopmodel)


## Birth model
######################################################################

birth_pop_model <- data %>% dplyr::select("sex", "age_1970", "returnedkarelia",
                                          "Birth_pop_log")
birth_pop_complete <- birth_pop_model[complete.cases(birth_pop_model),]

birthmodel <- glm(
  returnedkarelia ~ sex + age_1970 + Birth_pop_log,
  data = birth_pop_model, family = binomial(link = "logit")
)

summary(birthmodel)


## FDF Population model with bedlan data
######################################################################
bedlan <- read.csv2("~/bedlan.csv")
bedlan$SwePROSENT <- with(bedlan, RUOTSINK / YHTEENSA * 100)
dataB <- left_join(data, bedlan, by = c("fdf_name" = "NIMI"))

dataS <- left_join(data, all_evacuees, by = c("fdf_name" = "fdf_name"))
dataS$SCprec <- with(dataS, freqAll /fdf_population  * 100)
dataS$freqAll_log <- log(dataS$freqAll)

dataS <- dataS %>% dplyr::select("id","SCprec")
dataB <- left_join(dataB, dataS, by = c("id" = "id"))

dataB$MEANTEMP1 <- dataB$MEANTEMP - min(dataB$MEANTEMP, na.rm=TRUE)
dataB$MEANTEMP1 <- dataB$MEANTEMP / max(dataB$MEANTEMP, na.rm=TRUE)

dataB$RAINYDAYS1 <- dataB$RAINYDAYS - min(dataB$RAINYDAYS, na.rm=TRUE)
dataB$RAINYDAYS1 <- dataB$RAINYDAYS / max(dataB$RAINYDAYS, na.rm=TRUE)

dataB$RAINFALL1 <- dataB$RAINFALL - min(dataB$RAINFALL, na.rm=TRUE)
dataB$RAINFALL1 <- dataB$RAINFALL / max(dataB$RAINFALL, na.rm=TRUE)

dataB$SNOWDEPTH1 <- dataB$SNOWDEPTH - min(dataB$SNOWDEPTH, na.rm=TRUE)
dataB$SNOWDEPTH1 <- dataB$SNOWDEPTH / max(dataB$SNOWDEPTH, na.rm=TRUE)

dataB$SNOWYDAYS1 <- dataB$SNOWYDAYS - min(dataB$SNOWYDAYS, na.rm=TRUE)
dataB$SNOWYDAYS1 <- dataB$SNOWYDAYS / max(dataB$SNOWYDAYS, na.rm=TRUE)

dataB$SwePROSENT1 <- dataB$SwePROSENT - min(dataB$SwePROSENT, na.rm=TRUE)
dataB$SwePROSENT1 <- dataB$SwePROSENT / max(dataB$SwePROSENT, na.rm=TRUE)

dataB$SCprec1 <- dataB$SCprec - min(dataB$SCprec, na.rm=TRUE)
dataB$SCprec1 <- dataB$SCprec / max(dataB$SCprec, na.rm=TRUE)

dataBcomp <- dataB[ which(dataB$MEANTEMP1 < Inf), ]
mean(dataBcomp$MEANTEMP1)

dataMT<- dataB %>% dplyr::select(sex,age_1970,fdf_latitude1,fdf_longitude1,
                                 FDF_pop_log,returnedkarelia,MEANTEMP1,RAINFALL1,RAINYDAYS1,
                                 SNOWYDAYS1,SNOWDEPTH1,SwePROSENT1,SCprec1)
dataMT <- dataMT[complete.cases(dataMT),]
# make a correlation matrix out of dataframe m
cor_mat <- cor(dataMT)

# correlation threshold
cor_threshold <- 0.5

# set correlations to 0 if absolute value bigger than threshold
cor_mat <- ifelse(abs(cor_mat) > cor_threshold, 0 , 1)

# generate the lower-triangular matrix with correlations
cor_mat[upper.tri(cor_mat, diag = TRUE)] <- NA

options(na.action = "na.fail")

modeltest <- glm(
  returnedkarelia ~ sex + age_1970 + fdf_latitude1 + fdf_longitude1 + FDF_pop_log + MEANTEMP1 +
    RAINYDAYS1 + SNOWDEPTH1 + SwePROSENT1 + SCprec1,
  data=dataMT, family = binomial(link = "logit")
)

summary(modeltest)


#Get top model
modelset<-dredge(modeltest, rank = AICc, trace=FALSE,subset = cor_mat )

#modelset
summary(modelset)


write.csv(modelset, file = "~/modelset.csv")

################## AVERAGE MODELS WITHIN 2 AICC POINTS ###############
avgmodel<-model.avg(modelset, subset = delta < 2 )

topmodel<-get.models(modelset, subset = 1) [[1]]
summary(topmodel, type = "response")


# Maps
######################################################################

## Color coded return to Karelia from Finland map
######################################################################

returnees_only <- data[ which(data$returnedkarelia==1), ]

returnees_with_fdf <- as.data.frame(table(returnees_only$fdf_name))
names(returnees_with_fdf)[ names(returnees_with_fdf) == "Var1" ] <- "fdf_name"
names(returnees_with_fdf)[ names(returnees_with_fdf) == "Freq" ] <- "freqReturn"

all_evacuees <- left_join(all_evacuees, returnees_with_fdf, by = c("fdf_name" = "fdf_name"))
all_evacuees$percentage <- with(all_evacuees, freqReturn / freqAll * 100)

drive_map1 <- left_join(drive, all_evacuees, by = c ("placename"="fdf_name"))
drive_map1 <- drive_map1[ which(drive_map1$freqAll > 0), ]
drive_map1 <- drive_map1[ !duplicated(drive_map1$Latitude), ]

drive_map1$popR <- cut(drive_map1$Population, breaks = c(-Inf,5000,10000,15000,20000,Inf),
                       labels = c("2", "4", "6", "8", "10"))

map50 <- drive_map1[ which (drive_map1$freqAll > 49), ]
write.csv(map50, "~/map50.csv")


## Birthplace return rate map
######################################################################

evacuee_birthplaces <- as.data.frame(table(data$birthplace))
names(evacuee_birthplaces)[ names(evacuee_birthplaces) == "Var1" ] <- "birthplace"
names(evacuee_birthplaces)[ names(evacuee_birthplaces) == "Freq" ] <- "freqAll"

returnee_rdk_names <- as.data.frame(table(data$rdk_name))
names(returnee_rdk_names)[names(returnee_rdk_names)== "Var1"] <- "rdk_name"
names(returnee_rdk_names)[names(returnee_rdk_names)== "Freq"] <- "freqReturn"

evacuee_birthplaces <- left_join(evacuee_birthplaces, returnee_rdk_names, by = c("birthplace"="rdk_name"))
evacuee_birthplaces$percentage <- with(evacuee_birthplaces, freqReturn/freqAll * 100)

kareliaplaces <- left_join(drive, evacuee_birthplaces, by = c("placename"="birthplace"))

drive_map2 <- kareliaplaces[ which(kareliaplaces$freqAll > 0), ]
drive_map2 <- drive_map2[ !duplicated(drive_map2$Latitude), ]

mapBA100 <- drive_map2[ which(drive_map2$freqAll > 99), ]
write.csv(mapBA100, "~/mapBA100.csv")
