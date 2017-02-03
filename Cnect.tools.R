##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title Function to create the zones of release ans sink for Ichthyop
##' @param lon : matrix witht the values of longitude
##' @param lat : Matrix with the values of latitude
##' @param rel.pol number of the polygon of realease. this correspond to the row of the matris lat and lon
##' @param rec.pol number of the polygon of recruitment. this correspond to the row of the matris lat and lon. The default value is the same than rel.pol
##' @param rel.depth Depth of recruitment
##' @param rec.depth Depth of larval release. if null use the default values (0,50) meters
##' @param Atlantis If the information came from the atlantis,  is necesary to considere that the polygons start at 0
##' @param pol.names Vector with the names of the polygons
##' @param name.file Name of the output file
##' @param z.names Name of the zones. useful for ichthyops
##' @param rec_z.names Name of the zones of release
##' @return the file the zones of ichthyops
##' @author Demiurgo
zones.xml <- function(lon, lat, rel.pol, rec.pol = NULL, rel.depth = NULL, rec.depth = NULL,
                      Atlantis = TRUE, pol.names, name.file = NULL){
    if(is.null(rec.pol))   rec.pol    <- rel.pol
    if(is.null(rel.depth)) rel.depth  <- c(0, 50) # default values
    if(is.null(rec.depth)) rec.depth  <- c(0, 1200)
    if(is.null(name.file)) name.file  <- 'NO_NAME'
    if(isTRUE(Atlantis)){
        rel.pol <- rel.pol + 1
        rec.pol <- rec.pol + 1
    }
    z.names     <- pol.names[rel.pol]
    rec_z.names <- pol.names[rec.pol]


    ## creating the release areas
    sink(paste(name.file, ".xml", sep = ''))
    cat('<?xml version="1.0" encoding="UTF-8"?>\n')
    cat('<zones>\n')
    for( i in 1 : length(rel.pol)){
        cat('<zone>\n')
        cat('<key>', paste('rel', z.names[i],sep='_'),'</key>\n')
        cat('<enabled>true</enabled>\n')
        cat('<type>release</type>\n')
        cat('<polygon>\n')
        for(j in 1 : length(na.omit(lat[rel.pol[i], ]))){
            cat('<point>\n')
            cat('<index>', j - 1, '</index>\n', sep = '')
            cat('<lon>', lon[rel.pol[i], j], '</lon>\n', sep = '')
            cat('<lat>', lat[rel.pol[i], j], '</lat>\n', sep = '')
            cat('</point>\n')
        }
        cat('</polygon>\n')
        cat('<bathy_mask>\n')
        cat('<enabled>true</enabled>\n')
        cat('<line_inshore>0.0</line_inshore>\n')
        cat('<line_offshore>12000.0</line_offshore>\n')
        cat('</bathy_mask>\n')
        cat('<thickness>\n')
        cat('<enabled>true</enabled>\n')
        cat('<upper_depth>', rel.depth[1], '</upper_depth>\n', sep = '')
        cat('<lower_depth>', rel.depth[2], '</lower_depth>\n', sep = '')
        cat('</thickness>\n')
        cat('<color>[r=255,g=255,b=255]</color>\n')
        cat('</zone>\n')
    }
    ## creating the sin areas
    for( i in 1 : length(rec.pol)){
        cat('<zone>\n')
        cat('<key>', paste('rec', rec_z.names[i], sep = '_'),'</key>\n')
        cat('<enabled>true</enabled>\n')
        cat('<type>recruitment</type>\n')
        cat('<polygon>\n')
        for(j in 1 : length(na.omit(lat[rec.pol[i], ]))){
            cat('<point>\n')
            cat('<index>', j - 1, '</index>\n', sep = '')
            cat('<lon>', lon[rec.pol[i], j], '</lon>\n', sep = '')
            cat('<lat>', lat[rec.pol[i], j], '</lat>\n', sep = '')
            cat('</point>\n')
        }
        cat('</polygon>\n')
        cat('<bathy_mask>\n')
        cat('<enabled>true</enabled>\n')
        cat('<line_inshore>0.0</line_inshore>\n')
        cat('<line_offshore>12000.0</line_offshore>\n')
        cat('</bathy_mask>\n')
        cat('<thickness>\n')
        cat('<enabled>true</enabled>\n')
        cat('<upper_depth>', rec.depth[1], '</upper_depth>\n', sep = '')
        cat('<lower_depth>', rec.depth[2], '</lower_depth>\n', sep = '')
        cat('</thickness>\n')
        cat('<color>[r=0,g=255,b=153]</color>\n')
        cat('</zone>\n')
    }
    cat('</zones>')

    sink()
}

##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title Creation of the connectivity matrix from ichthyos output
##' @param dirpath Path to the directory with the netcdf files
##' @param prop Default TRUE,  calculate the value
##' @param real
##' @return A array with the structure A*B*C were A is the source, B is the sink and C in the nuber of years
##' @author Demiurgo
connect.matrix<- function(dirpath, prop = TRUE, real = TRUE){
    library('RNetCDF')
    library(XML)
    ## The directory that contains the series of netcdf input files
    filenames   <- list.files(path = dirpath, full.names = TRUE)
    nc          <- open.nc(filenames) ## open the file
    filezone    <- att.get.nc(nc,'NC_GLOBAL', 'release.zone.zone_file')
    simu        <- length(strsplit(att.get.nc(nc,'NC_GLOBAL', 'release.schedule.events'), '\\\"year')[[1]]) - 1
    t.zones     <- att.get.nc(nc,'NC_GLOBAL', 'nb_zones')
    time.end    <- dim.inq.nc(nc, 'time')$length
    rec.zones   <- dim.inq.nc(nc, 'recruitment_zone')$length
    rel.zones   <- t.zones - rec.zones
    n.part      <- dim.inq.nc(nc, 'drifter')$length
    prt.per.sim <- n.part / simu
    ## Particule relaese zone
    pr.zone     <- rep(sort(rep(seq( 1 : rel.zones), (prt.per.sim / rel.zones))), simu)
    ## matching the number of release particles!!
    pr.zone       <- c(pr.zone, rep(NA, n.part - length(pr.zone)))
    connect.mat   <- array(NA, c(rel.zones, rec.zones, simu))
    step.r        <- seq(from = 0, to = n.part, by = prt.per.sim)

    ## Main loop to bould the connectivity matrix by year (or simulation)
    for( z in 1 : rec.zones){
        recruited <- var.get.nc(nc, 'recruited_zone', c(z, 1, time.end), c(1, NA, 1))
        by.zone   <- recruited * pr.zone
        rec.ini   <- rep(0, rel.zones)
        i=1
        for(i in 1 : simu){
            ##recruited <- var.get.nc(nc, 'recruited_zone', c(z, 1, time.rec[i]), c(1, NA, 1))
            by.zone.h <- hist(by.zone[(step.r[i] + 1) : step.r[i + 1]], breaks = seq( 0, rel.zones + 1) - 0.5,
                              plot=FALSE)$counts[2 :(rel.zones + 1)]
            ## it's necesary to remove the particules for the previouys simulation
            connect.mat[, z, i] <- by.zone.h

        }
    }
    ### Calculatin the proportions
    if(isTRUE(prop)){
        if(isTRUE(real)){
            ## Asumming that you have looses out of the system
            for(year in 1 : simu){
                connect.mat[, , year] <- connect.mat[, , year] / prt.per.sim
            }
        } else {
            ## The connectivity matrix sum to 1
            connect.mat[, , year] <-  t(apply(connect.mat[, , year], 1, function(x) x / sum(x)))
        }
    }
    ## Zone's names
    zones   <- paste('zone', c(0 : (t.zones -1)),sep='')
    zones.n <- sapply(zones,function(x)  att.get.nc(nc, x,'long_name'))
    rownames(connect.mat) <- zones.n[1 : rel.zones]
    colnames(connect.mat) <- zones.n[(rel.zones+1) : length(zones.n)]
    return(connect.mat)
}

##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title Creation of the netcdf based on the connectivity matrix thought time
##' @param connect.l List of matrix witht the connectivity
##' @param names Names of the fucntional groups
##' @return A netcdf file ready to use for atlantis
##' @author Demiurgo
create.netcdf <- function(connect.l, names){
    ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
    ## ~            Conectivity recruitment         ~ ##
    ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##

    ## Libraries
    library(ncdf4)
    ## Read file
    ## -  -  -  -  -  -  -  -
    ## Sure that we are dealing with a list
    if(!is.list(connect.l)) connect.l <- list(connect.l)
    ## Create connectivity array
    time     <- 1 : dim(connect.l[[1]])[3]  ## all the matrix and output should have the same time duration
    n.pol    <- 0 : 50
    n.poly   <- c('Poly_0', 'Poly_1', 'Poly_2', 'Poly_3', 'Poly_4', 'Poly_5', 'Poly_6', 'Poly_7', 'Poly_8', 'Poly_9', 'Poly_10',
                  'Poly_11', 'Poly_12', 'Poly_13', 'Poly_14', 'Poly_15', 'Poly_16', 'Poly_17', 'Poly_18', 'Poly_19', 'Poly_20',
                  'Poly_21', 'AS_S_out', 'AS_N_out', 'JF5_N_out', 'RCSC_N_out', 'RCSC_S_out', 'JF5_S_out', 'JF4', 'JF3', 'JF2',
                  'JF1_E', 'JF1_M', 'JF1_NO', 'JF1_SO', 'RCSC_NE', 'RCSC_SE', 'RCSC_NO', 'RCSC_SO', 'JF6_S', 'JF6_N', 'AS_NE',
                  'AS_SE', 'AS_NO', 'AS_SO', 'RCSC', 'JF5', 'JF6', 'AS', 'BO2', 'BO1')

    t     <- ncdim_def(name     = 't',
                       units    = 'year since 1990-01-01 00:00:00 +10',
                       longname = 'Time',
                       unlim    = TRUE,
                       vals     = as.double(time))
    d      <- ncdim_def(name     = 'd',
                        longname = 'Donor',
                        vals     = n.pol, unit = 'prop')

    r      <- ncdim_def(name     = 'r',
                        longname = 'Recip',
                        vals     = n.pol, unit = 'prop')

    matrices <- list(array(0, c(length(n.pol), length(n.pol), length(time))))
    matrices <- rep(matrices, length(connect.l))
    for( fg in 1 : length(connect.l)){
        assign(paste('con.mat.', names[fg], sep = ''),
               ncvar_def(name     = paste(names[fg],'_Connnectivity', sep = ''),
                         longname = paste('Connectivity of ', names[fg], sep = ''),
                         units    = 'proportion',
                         dim      = list(r, d, t),
                         missval  = -9999,
                         prec     = 'double'))
        loc.col  <- match(gsub('rec_', '', colnames(connect.l[[fg]])), n.poly)
        loc.row  <- match(gsub('rel_', '', rownames(connect.l[[fg]])), n.poly)
        ## put the values in the new matrix
        if(dim(connect.l[[fg]])[3] !=  length(time)){
            ## this arregement is to avoid the problem of two spawning perior per year
            cat('\n WARNINGS -  - ', names[fg], ',  has different number of simulations.  check that\n')
            arre <- seq(from = 1, to = dim(connect.l[[fg]])[3], by = 2)
            matrices[[fg]][loc.col, loc.row, time] <- connect.l[[fg]][, , arre]
        } else {
            matrices[[fg]][loc.col, loc.row, time] <- connect.l[[fg]]
        }
    }
    ## list witht the varialbes
    vars   <- mget(ls(pattern = 'con.mat.*'))
    o.file <- 'FG_connect.nc'
    nc     <- nc_create(o.file, vars)
    ## putting the variables
    for(fg in 1 : length(connect.l)){
        ncvar_put(nc, vars[[fg]], matrices[[fg]])
    }


    nc_close(nc)
}
