## Homogenize each daily series with respect to the previous one in the given order

library(dataresqc)
library(lubridate)

inpath <- "../Zurich/SEF/5_monthly"
dailypath <- "../Zurich/SEF/4_daily"
refpath <- "../Bern/SEF/refs"
outpath <- "../Zurich/SEF/6_hom"

target_series <- c("Zurich_Observatory", "GCOS_Winterthur", "Bruegger", "Kuesnacht",
                   "NGZ_1842", "NGZ_1836", "Hofmeister", "SNG_1831", "Horner_1823",
                   "Feer", "Escher_181606", "Escher_181601", "Horner_1812",
                   "Hirzel_1795", "Muralt_1770", "Hirzel_1761", "Muralt_1760", 
                   "Meyer", "Ott") # NB: the first series will not be homogenized

ref_series <- c("Delemont", "Solothurn", "Weissenstein", "StGallen",
                "Aarau_Zschogge", "Uetliberg", "Nufenen", "Neuchatel",
                "Annone", "Schaffhausen", "Wolleb", "Ryhiner",
                "Herisau_Nef", "Herisau_Merz", "Luzern_Ineichen", "Huber_1789",
                "Merian_1828", "Merian_1833", "Merian_1837",
                "Aarau_Bronner", "Zug", "Gotthard", "Waldenburg",
                "Chur_Buol", "Chur_Herold", "Winterthur_Unknown",
                "Tegerfelden", "Gais", "Ellikon", "Mulhouse", "Socin",
                "Bern_Observatory", "Reinhard", "Burgdorf_1851", "Koch", 
                "Wolf", "Benoit", "Trechsel_1848", "Trechsel_1826", 
                "Studer_1803", "Studer_1801", "Sutz", "Vevey",
                "Studer_1797", "Bueren", "Lombach_1785", "Lombach_1777", 
                "Studer_1779", "Tavel", target_series)

nmin <- 60 # minimum overlap in months above which reference series are not used

breakpoints <- read.csv("../Metadata/Zurich_breakpoints_ta.csv", sep=";", stringsAsFactors=FALSE)
breakpoints_ref <- read.csv("../Metadata/Refs_breakpoints_ta.csv", sep=";", stringsAsFactors=FALSE)

target_files <- list.files(inpath, pattern="_ta_", full.names=TRUE)
target_files_daily <- list.files(dailypath, pattern="_ta_", full.names=TRUE)
ref_files <- list.files(refpath, pattern="_ta_", full.names=TRUE)


## Remove all files in output folder
file.remove(list.files(outpath,full.names=TRUE))


## Read first target series (assumed homogeneous)
x <- list() # list of homogeneous segments
x[[1]] <- read_sef(target_files[grep(target_series[1],target_files)])
x[[1]] <- x[[1]][which(!is.na(x[[1]]$Value)), ]
x[[1]]$dates <- as.Date(paste(x[[1]]$Year,x[[1]]$Month,1,sep="-"), format="%Y-%m-%d")
x[[1]]$station <- target_series[1]
x[[1]]$start <- x[[1]]$dates[1]
x[[1]]$end <- rev(x[[1]]$dates)[1]
x[[1]]$Value_hom <- x[[1]]$Value
file.copy(target_files_daily[grep(target_series[1],target_files_daily)], outpath)


## Read in series to be homogenized (fill in list x)
j <- 1
for (i in 2:length(target_series)) {
  xtmp <- read_sef(target_files[grep(target_series[i],target_files)])
  xtmp <- xtmp[which(!is.na(xtmp$Value)), ]
  xtmp$dates <- as.Date(paste(xtmp$Year,xtmp$Month,1,sep="-"), format="%Y-%m-%d")
  xtmp$station <- target_series[i]
  dailyfile <- target_files_daily[grep(target_series[i],target_files_daily)]
  if (length(dailyfile) == 0 | target_series[i]=="Tavel") {
    xtmp_daily <- xtmp
  } else {
    xtmp_daily <- read_sef(dailyfile)
    xtmp_daily$dates <- as.Date(paste(xtmp_daily$Year,xtmp_daily$Month,xtmp_daily$Day,sep="-"), 
                                format="%Y-%m-%d")
  }
  
  ## Read breakpoints
  id <- strsplit(target_series[i], "_")[[1]][1]
  if (id == "GCOS") id <- strsplit(target_series[i], "_")[[1]][2]
  bps <- breakpoints$Date[breakpoints$Station==id]
  if (length(bps) > 0) {
    bps <- as.Date(bps, format="%d.%m.%Y")
    bps <- c(min(xtmp_daily$dates), bps, max(xtmp_daily$dates))
    bps <- bps[bps>=bps[1] & bps<=bps[length(bps)] & !duplicated(bps)]
  } else {
    bps <- c(min(xtmp_daily$dates), max(xtmp_daily$dates))
  }
  if (length(dailyfile) == 0 | target_series[i]=="Tavel") {
    bps[length(bps)] <- bps[length(bps)] + days_in_month(bps[length(bps)]) - 1
  }
  
  ## Loop on breakpoints
  if (length(bps) > 2) {
    for (i_bp in length(bps):2) {
      j <- j + 1
      if (day(bps[i_bp]) >= 28) { # use breakpoint month only if breakpoint is at the end
        x[[j]] <- xtmp[xtmp$dates<=bps[i_bp]&xtmp$dates>=bps[i_bp-1], ]
      } else {
        x[[j]] <- xtmp[xtmp$dates<=(bps[i_bp]-28)&xtmp$dates>=bps[i_bp-1], ]
      }
      x[[j]]$start <- bps[i_bp-1]
      x[[j]]$end <- bps[i_bp]
    }
  } else {
    j <- j + 1
    x[[j]] <- xtmp
    x[[j]]$start <- bps[1]
    x[[j]]$end <- bps[2]
  }
  if (target_series[i] == "Escher_181601") ie <- j
}


## Read in reference series (fill in list y)
y <- list()
j <- 0
for (i in 1:length(ref_series)) {
  ytmp <- read_sef(ref_files[grep(ref_series[i],ref_files)])
  ytmp <- ytmp[which(!is.na(ytmp$Value)), ]
  ytmp$dates <- as.Date(paste(ytmp$Year,ytmp$Month,1,sep="-"), format="%Y-%m-%d")
  ytmp$station <- ref_series[i]
  
  ## Read breakpoints
  id <- strsplit(ref_series[i], "_")[[1]]
  if (length(id) > 1) {
    if (id[1] == "GCOS") {
      id <- id[2]
    } else if (substr(id[2],1,2) %in% c("17","18")) {
      id <- id[1]
    } else {
      id <- ref_series[i]
    }
  } else {
    id <- ref_series[i]
  }
  bps <- breakpoints_ref$Date[breakpoints_ref$Station==id]
  if (length(bps) > 0) {
    bps <- as.Date(bps, format="%d.%m.%Y")
    bps <- c(min(ytmp$dates), bps, max(ytmp$dates))
    bps <- bps[bps>=bps[1] & bps<=bps[length(bps)]]
  } else {
    bps <- c(min(ytmp$dates), max(ytmp$dates))
  }
  bps[length(bps)] <- bps[length(bps)] + days_in_month(bps[length(bps)]) - 1
  
  ## Loop on breakpoints
  if (length(bps) > 2) {
    for (i_bp in length(bps):2) {
      j <- j + 1
      if (day(bps[i_bp]) >= 28) { # use breakpoint month only if breakpoint is at the end
        y[[j]] <- ytmp[ytmp$dates<=bps[i_bp]&ytmp$dates>=bps[i_bp-1], ]
      } else {
        y[[j]] <- ytmp[ytmp$dates<=(bps[i_bp]-28)&ytmp$dates>=bps[i_bp-1], ]
      }
      if (nrow(y[[j]]) == 0) {
        y[[j]] <- NULL
        j <- j - 1
      }
    }
  } else {
    j <- j + 1
    y[[j]] <- ytmp
  }
}


# Homogenization
last_station <- "0"
months <- 1:12
days <- 1:365
trigon <- cbind(cos(2*pi*(months/12)), sin(2*pi*(months/12)), 
                cos(4*pi*(months/12)), sin(4*pi*(months/12)))
trigon_daily <- cbind(cos(2*pi*(days/365)), sin(2*pi*(days/365)), 
                      cos(4*pi*(days/365)), sin(4*pi*(days/365)))
for (i in 2:length(x)) {
  
  if (x[[i]]$station[1]!=last_station) { # write homogenized daily file and read in the next to homogenize
    if (!last_station %in% c("0","Escher_181601")) {
      filename <- target_files_daily[grep(last_station,target_files_daily)]
      filename <- rev(strsplit(filename,"/")[[1]])[1]
      out <- x_daily[,c("Year","Month","Day","Hour","Minute","Value")]
      out$Value <- round(out$Value, 1)
      x_daily$Meta <- sub("^[|]", "", x_daily$Meta)
      write_sef(out, outpath, meta["var"], meta["id"], meta["name"], meta["lat"], 
                meta["lon"], meta["alt"], meta["source"], meta["link"], meta["units"],
                meta["stat"], meta["meta"], x_daily$Meta, period="day",
                outfile=filename, keep_na=TRUE)
    }
    last_station <- x[[i]]$station[1]
    if (last_station != "Escher_181601") {
      x_daily <- read_sef(target_files_daily[grep(last_station,target_files_daily)], all=TRUE)
      meta <- read_meta(target_files_daily[grep(last_station,target_files_daily)])
      x_daily <- x_daily[!grepl("qc=",x_daily$Meta), ]
      x_daily$dates <- as.Date(paste(x_daily$Year,x_daily$Month,x_daily$Day,sep="-"), 
                               format="%Y-%m-%d")
      x_daily$julian <- as.integer(format(x_daily$dates,"%j"))
    }
  }
  
  message("WORKING ON ", toupper(x[[i]]$station[1]), ": ", 
          min(x[[i]]$Year), "-", max(x[[i]]$Year))
  x[[i]]$Value_hom <- NA
  
  ## Merge all segments homogenized so far into one
  xhom <- Reduce(rbind, x[1:(i-1)])
  xhom <- xhom[!duplicated(xhom$dates), ]
  xhom <- xhom[order(xhom$dates), ]
  n <- sum(x[[i]]$dates %in% xhom$dates)
  adjustments <- rep(NA, 12)
  
  if (n >= 12) { # homogenized data as reference
    message("Using homogenized data as reference (", n, " months of overlap)")
    x1 <- xhom[xhom$dates%in%x[[i]]$dates, ] # homogeneous series
    x2 <- x[[i]][x[[i]]$dates%in%xhom$dates, ] # series to homogenize
    for (im in 1:12) {
      adjustments[im] <- mean(x1$Value_hom[x1$Month==im]-x2$Value[x2$Month==im], na.rm=TRUE)
    }
  }
  
  if (n < nmin) { # reference stations as reference
    n1 <- n2 <- c()
    for (j in 1:length(y)) {
      if (y[[j]]$station[1] != x[[i]]$station[1]) {
        n1 <- append(n1, sum(xhom$dates%in%y[[j]]$dates & xhom$dates>max(x[[i]]$dates)))
        n2 <- append(n2, sum(x[[i]]$dates%in%y[[j]]$dates))
        if (n1[j]>=12 & n2[j]>=12) {
          message("Using ", y[[j]]$station[1], " as reference (", 
                  n1[j],"+",n2[j], " months of overlap)")
        }
      } else {
        n1 <- append(n1, 0)
        n2 <- append(n2, 0)
      }
    }
    k <- which(n1>=12 & n2>=12)
    if (length(k) == 0 | (x[[i]]$station[1]=="Muralt"&x[[i]]$Year[1]==1770)) { # Muralt has data only for Jan-Apr
      k <- which(n1>0 & n2>0)
      ref_adjustments <- array(dim=c(length(k),12))
      for (ik in 1:length(k)) {
        xref1 <- y[[k[ik]]][y[[k[ik]]]$dates%in%xhom$dates & y[[k[ik]]]$dates>max(x[[i]]$dates), ]
        xref2 <- y[[k[ik]]][y[[k[ik]]]$dates%in%x[[i]]$dates, ]
        x1 <- xhom[xhom$dates%in%xref1$dates, ] # series to be homogenized against
        x2 <- x[[i]][x[[i]]$dates%in%xref2$dates, ] # series to homogenize
        message("Using ", xref1$station[1], " as reference (", 
                n1[k[ik]],"+",n2[k[ik]], " months of overlap)")
        for (im in months) {
          ref_adjustments[ik,im] <- mean(x1$Value_hom[x1$Month==im]-xref1$Value[xref1$Month==im],na.rm=TRUE) -
            mean(x2$Value[x2$Month==im]-xref2$Value[xref2$Month==im],na.rm=TRUE)
        }
        ref_adjustments[ik,] <- rep(mean(ref_adjustments[ik,],na.rm=TRUE), 12)
      }
      adjustments <- rbind(adjustments, ref_adjustments)
      warning(paste("Constant correction applied to", x[[i]]$station[1]))
    } else {
      ref_adjustments <- array(dim=c(length(k),12))
      for (ik in 1:length(k)) {
        xref1 <- y[[k[ik]]][y[[k[ik]]]$dates%in%xhom$dates & y[[k[ik]]]$dates>max(x[[i]]$dates), ]
        xref2 <- y[[k[ik]]][y[[k[ik]]]$dates%in%x[[i]]$dates, ]
        x1 <- xhom[xhom$dates%in%xref1$dates, ] # series to be homogenized against
        x2 <- x[[i]][x[[i]]$dates%in%xref2$dates, ] # series to homogenize
        for (im in months) {
          ref_adjustments[ik,im] <- mean(x1$Value_hom[x1$Month==im]-xref1$Value[xref1$Month==im],na.rm=TRUE) -
            mean(x2$Value[x2$Month==im]-xref2$Value[xref2$Month==im],na.rm=TRUE)
        }
      }
      adjustments <- rbind(adjustments, ref_adjustments)
    }
  }
  
  ## Take median and smooth adjustments
  print(adjustments)
  if (!is.null(dim(adjustments))) {
    adjustments <- apply(adjustments, 2, median, na.rm=TRUE)
  }
  params <- lm(adjustments ~ trigon)$coefficients
  adjustments <- params[1] + params[2]*trigon[,1] + params[3]*trigon[,2] + 
    params[4]*trigon[,3] + params[5]*trigon[,4]
  print(adjustments)
  adjustments_daily <- params[1] + params[2]*trigon_daily[,1] + params[3]*trigon_daily[,2] + 
    params[4]*trigon_daily[,3] + params[5]*trigon_daily[,4]
  for (im in months) {
    x[[i]]$Value_hom[x[[i]]$Month==im] <- round(x[[i]]$Value[x[[i]]$Month==im] + adjustments[im], 1)
  }
  
  ## Apply daily corrections
  adjustments_daily[366] <- adjustments_daily[365]
  for (j in 1:366) {
    k <- which(x_daily$julian==j & x_daily$date>=x[[i]]$start[1] & x_daily$date<=x[[i]]$end[1])
    x_daily$Value[k] <- x_daily$Value[k] + adjustments_daily[j]
    x_daily$Meta[k] <- paste0(x_daily$Meta[k], "|correction=", round(adjustments_daily[j],2))
  }
  
}

## Write last daily file
filename <- target_files_daily[grep(last_station,target_files_daily)]
filename <- rev(strsplit(filename,"/")[[1]])[1]
out <- x_daily[,c("Year","Month","Day","Hour","Minute","Value")]
out$Value <- round(out$Value, 1)
x_daily$Meta <- sub("^[|]", "", x_daily$Meta)
write_sef(out, outpath, meta["var"], meta["id"], meta["name"], meta["lat"], 
          meta["lon"], meta["alt"], meta["source"], meta["link"], meta["units"],
          meta["stat"], meta["meta"], x_daily$Meta, period="day",
          outfile=filename, keep_na=TRUE)


## Write Escher monthly file
out <- x[[ie]]
out$Hour <- out$Minute <- NA
out <- out[,c("Year","Month","Day","Hour","Minute","Value_hom")]
meta <- read_meta(target_files[grep(x[[ie]]$station[1],target_files)])
write_sef(out, outpath, meta["var"], meta["id"], meta["name"], meta["lat"], 
          meta["lon"], meta["alt"], meta["source"], meta["link"], meta["units"],
          meta["stat"], meta["meta"], period="month", note="monthly", keep_na=TRUE)
