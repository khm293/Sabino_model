#script to make slopes and dem in to a netCDF file

#load packages
library(abind)
library(ncdf4)
library(lattice)
library(fields)

setwd("/Users/katie/Dropbox/ParFlow Domain/slopes/")
fundir=("/Users/katie/Dropbox/ParFlow Domain/slopes/functions/")
outdir=("/Users/katie/Dropbox/netCDF stuff/netCDF_in_R/")

dem=matrix(scan("dem.txt"), ncol=246, byrow=T)
slopex=matrix(scan("slopex.txt"), ncol=246, byrow=T)
slopey=matrix(scan("slopey.txt"), ncol=246, byrow= T)

ny=nrow(dem)
nx=ncol(dem)
demT=t(dem[ny:1,])
image.plot(demT)

# write to netCDF file
#create path and outfile name
outpath <- outdir
outname <- paste("elevslopes",".nc", sep="")
outfname <- paste(outpath, outname, "", sep="")

x=seq(from=0, to=245, by=1)
y=seq(from=0,to=177, by=1)


# create and write the netCDF file -- ncdf4 version
# define dimensions
xdim <- ncdim_def(name="lon", units='', longname = '', vals= x) 
ydim <- ncdim_def(name="lat", units='', longname='', vals= y) 


# define variables
fillvalue <- 1e32
dlname <- "elevation"
elev_def <- ncvar_def(name="Elev",units="m",dim=list(xdim,ydim),missval=fillvalue,longname=dlname, prec="double")
dlname <- "slopes in x-direction"
x_def <- ncvar_def(name="slopex",units="m/m",dim=list(xdim,ydim),missval=fillvalue,longname=dlname, prec="double")
dlname <- "slopes in y-direction"
y_def <- ncvar_def(name="slopey",units="m/m",dim=list(xdim,ydim),missval=fillvalue,longname=dlname, prec="double")

vars=list(elev_def, x_def, y_def)

# create netCDF file and put arrays
cdf=nc_create(outfname, vars, force_v4=TRUE, verbose=TRUE)

# put variables
ncvar_put(cdf,elev_def,demT)
ncvar_put(cdf,x_def,slopex)
ncvar_put(cdf,y_def,slopey)


# Get a summary of the created file:
cdf
nc_close(cdf)
