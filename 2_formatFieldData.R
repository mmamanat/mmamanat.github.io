
## Load required packages
require(rgdal)
require(raster)
require(ggplot2)
require(reshape2)

## Set main working directory 
main.dir <- "C:/Users/mmama/Downloads/GIS714/"
setwd(main.dir)

## Read in lakes shapefile 
lakes <- readOGR(paste0(main.dir,"data/lakes"),"updatedValidLakes")
## Subset to Lake Jesup
jesup.shp <- subset(lakes, GNIS_NAM_1 == "Jesup, Lake")



## Read in field data 
## Data from: http://seminole.wateratlas.usf.edu/lake/waterquality.asp?wbodyid=7589&wbodyatlas=lake
wq.files <- lapply(list.files(paste0(main.dir,"data/water_quality"), "*.csv$", full.names = T), function(x){read.csv(x, header = T)})
wq.xy <- lapply(wq.files, function(x){x <- x[,c("Actual_Longitude","Actual_Latitude")]; return(x)})
wq.spdf <- lapply(1:length(wq.files), function(x){SpatialPointsDataFrame(coords = wq.xy[[x]], data = wq.files[[x]], proj4string =  CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))})
wq.data <- lapply(wq.spdf, function(x){spTransform(x, projection(jesup.shp))})
## Subset to data from 2018
wq.date <- lapply(wq.data, function(x){x$SampleDate <- as.Date(x$SampleDate, "%m/%d/%Y"); return(x)})
wq.subset.date <- lapply(wq.date, function(x){subset(x, substr(SampleDate,1,4) == "2018")})
## Subset to just include the columns of interest (location, date, value)
wq.subset.cols <- lapply(wq.subset.date, function(x){x[c("Actual_Latitude","Actual_Longitude","SampleDate","Parameter","Result_Value")]})
## Aggregate to average result value for each location and for each date 
wq.agg <- lapply(wq.subset.cols, function(x){aggregate(Result_Value ~ Actual_Latitude + Actual_Longitude + SampleDate, x, FUN = max)})
wq.agg.param <- lapply(1:length(wq.agg), function(x){wq.agg[[x]]$Parameter <- unique(wq.subset.cols[[x]]$Parameter); return(wq.agg[[x]])})
## Create a column for the month of the year for each date 
wq.monthly <- lapply(wq.agg.param, function(x){x$MONTH <- strftime(x$SampleDate, format = "%m"); x$SampleDate <- NULL; return(x)})
## Combine into a single dataframe
all.wq <- do.call(rbind, wq.monthly)
## Check for normality 
wq.normal <- lapply(as.list(as.character(unique(all.wq$Parameter))), function(x){shapiro.test(subset(all.wq, Parameter == x)$Result_Value)$p.value})
## Transform using the box-cox transformation 
#all.wq$Result_Value <- log(all.wq$Result_Value)
## Make a new dataframe for each unique month 
all.wq.split <- split(all.wq, f = all.wq$MONTH)
all.wq.split <- lapply(all.wq.split, function(x){x$MONTH <- NULL; return(x)})
## Reshape so each parameter is a column header 
all.wq.dcast <- lapply(all.wq.split, function(x){dcast(x, Actual_Latitude + Actual_Longitude ~ Parameter, value.var = "Result_Value", mean)})
## Remove any rows with a latitude greater than 28.76 (not actually in Jesup)
all.wq.rm <- all.wq.dcast



## Lakewatch data 
lakewatch <- lapply(list.files(paste0(main.dir,"data/lakewatch"), "*.csv$", full.names = T), function(x){read.csv(x, header = T)})
## Keep columns of interest 
lakewatch.subset <- lapply(lakewatch, function(x){x[c("Month","Latitude","Longitude","TP..µg.L.","TN..µg.L.","SECCHI..ft.")]})
lakewatch.all <- do.call(rbind, lakewatch.subset)
## Make into long format
lakewatch.melt <- melt(lakewatch.all, id.vars = c("Month","Latitude","Longitude"))
## Log transform data 
#lakewatch.melt$value <- log(lakewatch.melt$value)
## Make a unique dataframe for each month
lakewatch.split <- split(lakewatch.melt, f = lakewatch.melt$Month)
## Reshape
lakewatch.dcast <- lapply(lakewatch.split, function(x){dcast(x, Latitude + Longitude ~ variable, value.var = "value", mean)})



## SJRWMD data 
sjrwmd <- lapply(list.files(paste0(main.dir,"data/sjrwmd"), "*.csv$", full.names = T), function(x){read.csv(x, header = T)})
## Subset to year of interest and columns of interest 
sjrwmd.date <- lapply(sjrwmd, function(x){x$Date <- as.Date(x$Date, "%m/%d/%Y"); return(x)})
sjrwmd.subset <- lapply(sjrwmd.date, function(x){subset(x, substr(as.Date(x$Date, "%m/%d/%Y"),1,4) == "2018")})
sjrwmd.subset <- lapply(sjrwmd.subset, function(x){x[c("Date","Latitude","Longitude","Parameter","Measured_Value")]})
sjrwmd.subset <- lapply(sjrwmd.subset, function(x){x[which(x$Parameter %in% c("DO","Secchi","Turbidity")),]})
## Convert secchi from meters to feet 
sjrwmd.subset <- lapply(sjrwmd.subset, function(x){x$Measured_Value[which(x$Parameter == "Secchi")] <- x$Measured_Value[which(x$Parameter == "Secchi")] / 0.3048; return(x)})
## Make a column for month 
sjrwmd.monthly <- lapply(sjrwmd.subset, function(x){x$MONTH <- strftime(as.character(x$Date), format = "%m"); x$Date <- NULL; return(x)})
## Combine into single dataframe
sjrwmd.all <- do.call(rbind, sjrwmd.monthly)
## Log transform data 
#sjrwmd.all$Measured_Value <- log(sjrwmd.all$Measured_Value)
## Make a unique dataframe for each month
sjrwmd.split <- split(sjrwmd.all, f = sjrwmd.all$MONTH)
## Reshape
sjrwmd.dcast <- lapply(sjrwmd.split, function(x){dcast(x, Latitude + Longitude ~ Parameter, value.var = "Measured_Value", mean)})



## For each list of dataframes, give them the same column names 
## sjrwmd.dcast
sjrwmd.final <- lapply(sjrwmd.dcast, function(x){names(x) <- c("LAT","LON","DO","SECCHI","TURBIDITY"); return(x)})
sjrwmd.final <- lapply(sjrwmd.final, function(x){x$TN <- NA; x$TP <- NA; return(x)})
sjrwmd.final <- lapply(sjrwmd.final, function(x){x[c("LAT","LON","DO","TN","TP","SECCHI","TURBIDITY")]})
names(sjrwmd.final) <- c(1:12)
## lakewatch.dcast
lakewatch.final <- lapply(lakewatch.dcast, function(x){names(x) <- c("LAT","LON","TP","TN","SECCHI"); return(x)})
lakewatch.final <- lapply(lakewatch.final, function(x){x$DO <- NA; x$TURBIDITY <- NA; return(x)})
lakewatch.final <- lapply(lakewatch.final, function(x){x[c("LAT","LON","DO","TN","TP","SECCHI","TURBIDITY")]})
names(lakewatch.final) <- c(2:11)
## all.wq.rm
wq.final <- lapply(all.wq.rm, function(x){x$Chla_ugl <- NULL; return(x)})
wq.final <- lapply(wq.final, function(x){names(x) <- c("LAT","LON","DO","TN","TP","SECCHI","TURBIDITY"); return(x)})
wq.final <- lapply(wq.final, function(x){x[c("LAT","LON","DO","TN","TP","SECCHI","TURBIDITY")]})
names(wq.final) <- c(1:12)

## Rbind data from each source into data frame for each month 
output.data <- list()
output.data[[1]] <- rbind(sjrwmd.final[[1]], wq.final[[1]])
output.data[[2]] <- rbind(sjrwmd.final[[2]], wq.final[[2]], lakewatch.final[[1]])
output.data[[3]] <- rbind(sjrwmd.final[[3]], wq.final[[3]], lakewatch.final[[2]])
output.data[[4]] <- rbind(sjrwmd.final[[4]], wq.final[[4]], lakewatch.final[[3]])
output.data[[5]] <- rbind(sjrwmd.final[[5]], wq.final[[5]], lakewatch.final[[4]])
output.data[[6]] <- rbind(sjrwmd.final[[6]], wq.final[[6]], lakewatch.final[[5]])
output.data[[7]] <- rbind(sjrwmd.final[[7]], wq.final[[7]], lakewatch.final[[6]])
output.data[[8]] <- rbind(sjrwmd.final[[8]], wq.final[[8]], lakewatch.final[[7]])
output.data[[9]] <- rbind(sjrwmd.final[[9]], wq.final[[9]], lakewatch.final[[8]])
output.data[[10]] <- rbind(sjrwmd.final[[10]], wq.final[[10]], lakewatch.final[[9]])
output.data[[11]] <- rbind(sjrwmd.final[[11]], wq.final[[11]], lakewatch.final[[10]])
output.data[[12]] <- rbind(sjrwmd.final[[12]], wq.final[[12]])

## Combine to a single value per date 
output.zeros <- lapply(output.data, function(x){x[which(is.na(x),arr.ind = T)] <- 0; return(x)})
output.zeros <- lapply(output.zeros, function(x){round(x, 4)})
output.agg <- lapply(output.zeros, function(x){aggregate(. ~ LAT + LON, x, FUN = max)})
output.agg <- lapply(output.agg, function(x){x[which(x == 0,arr.ind = T)] <- NA; return(x)})


months <- c("JAN","FEB","MAR","APR","MAY","JUN","JUL","AUG","SEP","OCT","NOV","DEC")


## To SPDF
require(raster)
output.spdf <- lapply(output.agg, function(x){SpatialPointsDataFrame(x[,c("LON","LAT")], x, proj4string = CRS(projection(lakes)))})
save.spdf <- lapply(1:length(output.spdf), function(x){writeOGR(output.spdf[[x]], "spdfs", paste0(months[x],"_shp"), driver = "ESRI Shapefile")})


## Save results 

csv.data <- lapply(1:length(output.agg), function(x){write.csv(output.agg[[x]], paste0(main.dir, "csvs/", months[x],"_data.csv"), row.names = F)})



  