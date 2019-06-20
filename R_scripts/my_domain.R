# Workflow Example 1: Rectangular domain no river mask
# This example walks through the simplest case where
# you have a rectangular domain and no river network specified a-priori
# in this case the only input required is a DEM

setwd("/Users/katie/Dropbox/ParFlow Domain/slopes/")
fundir=("/Users/katie/Dropbox/ParFlow Domain/slopes/")
outdir=("/Users/katie/Dropbox/ParFlow Domain/slopes/")

##########################################
#Source Libraries and functions
rm(list=ls())
library('fields')
source("./functions/D4_Traverse.R")
source("./functions/Init_Queue.R")
source("./functions/Stream_Traverse.R")
source("./functions/Find_Orphan.R")
source("./functions/drainage_area.R")
source("./functions/Slope_Calc_Upwind.R")
source("./functions/Get_Border.R")
source("./functions/Define_Subbasins.R")
source("./functions/Write_Raster.R")

##########################################
#Settings
#Settings for the slope processing to change

#DEM processing
ep=0.01 #The epsilon value applied to flat cells

#Slope scaling
maxslope=0.5	#maximum slope (slopes with absolute value greater than this will be set to minslope), set to -1 if you don't want to set a threshold
minslope=1e-5	#minimum slope (slopes with absolute value less than this will be set to minslope), set to -1 if you don't want to set a threshold
scale=-1 #The maximum ratio of secondary to primary flow directions (set to -1 if you don't want the secondary slopes to be scaled, set to 0 if you want only primary flow directios)

#River and subbasin size for slope calculations
sub_th=100 #area threshold (number of grid cells) to use for subbasin delineation
riv_th=50 #optional additional area threshold (number of grid cells) to use for the river mask for slope processing. See notes below if you want to change this to be different from the subbasin threshold
riv_method=2 #method for processing river cellls (0=treat river cells the same as the rest of the domain, 1=set secondary slopes along the river to zero, 2=apply subbasin average slopes to river cells, 3=apply subbasin average river slopes of the river cells)
mrg_th=10	#Threshold number of grid cells for merging small subbasins in the subbasin analysis

#Grid dimensions for slopes
dx=90 #grid cell size for slope calculations
dy=90 #grid cell size for slope calcualtions

#Runname - for output writing
runname='tucson'

##########################################
# Read in Inputs
# The DEM should be formated as a matrix with the same
# dimensions as the domain
dem=matrix(scan("dem.txt"), ncol=246, byrow=T)
ny=nrow(dem)
nx=ncol(dem)
demT=t(dem[ny:1,]) #transforming the dem so it is indexed as [x,y] for functions

#check that you aren't upside down and backwards somehow...
#if you've formatted your input correctly this should look like your DEM
# without any additional transforming on the matrix
image.plot(demT)
tiny_dem=matrix(NA, ncol=10, nrow=7)
tiny_dem=demT[56:66,156:163]
image.plot(tiny_dem)

write.table(demT, 'demT.txt', row.names = FALSE, col.names = FALSE)

##########################################
# Process the DEM
#1. initialize the queue with all the rectangular border cells
#2. Process the DEM so that all cells drain to the boundaries

#1. Initialize queue
init=InitQueue(demT) #using the rectangular boundary
#2. Process DEM
travHS=D4TraverseB(demT, init$queue, init$marked, basins=init$basins, epsilon=ep)

#Look at the outputs
ls(travHS) #to see a list of everything that comes out of the processing
image.plot(travHS$dem) # the processed DEM
image(travHS$marked) # Mask of cells the processing algorithm covered (should be everything for this approach)
image(travHS$step) # The step at which each cell was processed in the algorithm
image(travHS$basins) # The resulting drainage basins


##########################################
# Calculate the slopes
# Note this step also fixes the directions of the borders because
# directions are not provided when the queue is initialized

### Option 1: just calcualte the slopes for the entire domain with no distinction between river and hillslope cells
#In this example secondary slope scaling is turned on and the secondary
#Slopes in the secondary direction are set to a maximum of 0.1*primary flow direction
#To calculate only slopes in the primary flow direction set the secondaryTH to 0
#Additionally primary slopes are limited by min slope and max slope thresholds
#slopesUW=SlopeCalcUP(dem=travHS$dem, direction=travHS$direction, dx=dx, dy=dy,  secondaryTH=scale, maxslope=maxslope, minslope=minslope)


### Option 2: If you would like to handle river cells differently from the rest of the domain

#do a preliminary slope calc just to get the flow directions on the boundary fixed
slopesUW=SlopeCalcUP(dem=travHS$dem, direction=travHS$direction, dx=dx, dy=dy,  secondaryTH=scale, maxslope=maxslope, minslope=minslope)

# Calculate the drainage area
area=drainageArea(slopesUW$direction, printflag=F)

# Define subbasins for calcualting river reach slopes
# the riv_th here is the drainage area threshold for splitting the river network branches
# when you do this you can still end up with subbasins with drainage areas less than the riv_th
# when multiple branches come together in a small area.
# To fix this you can set a merge threshold (merge_th) so that subbains with areas < merge_th autmoatically get merged with their downstream neighbor
subbasin=CalcSubbasins(slopesUW$direction, area, riv_th=sub_th, merge_th=mrg_th)
#plot the resulting subbasins and rivers
temp=subbasin$RiverMask
temp[temp==0]=NA
maskcol=colorRampPalette(c('black', 'black'))
#maskcol=colorRampPalette(c('white', 'white'))
image.plot(subbasin$subbasins)
image.plot((temp*2), add=T, col=maskcol(2), legend=F)

#Calculate the slopes
# The "river_method' flag here determines how the river cells will be handeled (e.g. using subbasin averages along reaches). Refer to the top of this script or the function for details.
slopesUW=SlopeCalcUP(dem=travHS$dem, direction=travHS$direction, dx=dx, dy=dy, secondaryTH=scale, maxslope=maxslope, minslope=minslope, river_method=riv_method, rivermask=subbasin$RiverMask, subbasin=subbasin$subbasins)

### Option 2b: Alternate more advanced approach: Define a river mask separate from the subbasin river mask and use this for the slope calculations. If you do this the average slopes will still be calculated
#by subbasin using the sub_th, but you can apply those average sloeps to more river cells by setting a lower threshold here. This is the 'riv_th' set at the top
#if you set riv_th=sub_th at the top this will have the same effect as just running the slope calc with the subbasin$RiverMask
rivers=area
rivers[area<riv_th]=0
rivers[area>=riv_th]=1

#plot the subbasins with the new river mask to check that the threshold is good
temp=rivers
temp[temp==0]=NA
maskcol=colorRampPalette(c('black', 'black'))
maskcol=colorRampPalette(c('white', 'white'))
image.plot(subbasin$subbasins)
image.plot((temp*2), add=T, col=maskcol(2), legend.only=FALSE)

slopesUW=SlopeCalcUP(dem=travHS$dem, direction=travHS$direction, dx=dx, dy=dy, secondaryTH=scale, maxslope=maxslope, minslope=minslope, river_method=riv_method, rivermask=rivers, subbasin=subbasin$subbasins)

#Look at the slopes and directions
image(slopesUW$slopex)
image(slopesUW$slopey)
image.plot(slopesUW$direction)

# slopesUW$slopex[156,58]=-1e-4
# slopesUW$slopex[157,58]=-1e-4
# slopesUW$slopex[158,58]=-1e-4

write.table(slopesUW$slopex, 'slopex.txt', row.names = FALSE, col.names = FALSE)
write.table(slopesUW$slopey, 'slopey.txt', row.names = FALSE, col.names = FALSE)

##########################################
# Calculate the drainage area - if you went with option 1 for slopes and you didn't do this already
area=drainageArea(slopesUW$direction, printflag=T) #rectangular boundary
image.plot(area)

##########################################
#Write the slopes out in PF format
slopeUWx=slopeUWy=rep(0, nx*ny)
jj=1
for(j in 1:ny){
	for(i in 1:nx){
		slopeUWx[jj]=slopesUW$slopex[i,j]
		slopeUWy[jj]=slopesUW$slopey[i,j]
		jj=jj+1
	}
}

fout=paste(runname, ".slopex.sa", sep="")
write.table( t(c(nx,ny,1)), fout, append=F, row.names=F, col.names=F)
write.table(slopeUWx, fout, append=T, row.names=F, col.names=F)
fout=paste(runname,".slopey.sa", sep="")
write.table( t(c(nx,ny,1)), fout, append=F, row.names=F, col.names=F)
write.table(slopeUWy, fout, append=T, row.names=F, col.names=F)

##########################################
#Example writing out other variables as matrices
#write.table( t(slopesUW$direction[,ny:1]) ,paste(runname, ".direction.out.txt", sep=""), row.names=F, col.names=F)
#write.table( t(travHS$dem[,ny:1]) ,paste(runname, ".dem.out.txt", sep=""), row.names=F, col.names=F)
#write.table( t(area[,ny:1]) , paste(runname, ".area.out.txt", sep=""), row.names=F, col.names=F)
#write.table( t(subbasin$subbasins[,ny:1]) , paste(runname, ".subbasins.out.txt", sep=""), row.names=F, col.names=F)
#write.table( t(subbasin$segments[,ny:1]) , paste(runname, ".subbasin_streams.out.txt", sep=""), row.names=F, col.names=F)
write.table( t(rivers[,ny:1]) , paste(runname,"rivermask.out.", runname,".txt", sep=""), row.names=F, col.names=F)


#Example writing out a variable as a raster
#write.raster( t(subbasin$segments[,ny:1]) , paste(runname, ".subbasin_streams.out.asc", sep=""), xllcorner=0.0, yllcorner=0.0, dx=dx, naval=-999)

## write out the subbasin summary information
#write.table(subbasin$summary, paste(runname, ".Subbasin_Summary.txt", sep=""), row.names=F)

##########################################
#Plotting
outdir=("/Users/katie/Desktop/ParFlow Domain/slopes/")
rivermask=rivers #This should be the river mask you input to the slope processing. could also be 'rivers' if you calcualted your own river mask for that

#Find all the headwater cells to use as starting points
#calculate the number of river cells draining to any cell
d4=c(1,2,3,4)
down=up=left=right=matrix(0, nrow=nx, ncol=ny) 
down[which(slopesUW$direction==d4[1])]=1
left[which(slopesUW$direction==d4[2])]=1
up[which(slopesUW$direction==d4[3])]=1
right[which(slopesUW$direction==d4[4])]=1
draincount=matrix(0, nrow=nx, ncol=ny)
draincount[,1:(ny-1)]=draincount[,1:(ny-1)]+down[,2:ny]*rivermask[,2:ny]
draincount[,2:ny]=draincount[,2:ny]+up[,1:(ny-1)]*rivermask[,1:(ny-1)]
draincount[1:(nx-1),]=draincount[1:(nx-1),]+left[2:nx,]*rivermask[2:nx,]
draincount[2:nx, ]=draincount[2:nx,]+right[1:(nx-1),]*rivermask[1:(nx-1),]

#Identify all the headwater cells
headwater=matrix(0, nrow=nx, ncol=ny)
headwater[which(draincount==0 & rivermask==1)]=1
#image.plot(headwater)
headlist=which(headwater==1, arr.ind=T)
colnames(headlist)=c("x", "y")
#headlist
nhead=nrow(headlist)

#Pick a point and walk downriver to the outlet
for(h in 1:nhead){
  #h=40 #picking one of the headwater cells
  x=headlist[h,1]
  y=headlist[h,2]
  
  #Plot to show you where you picked
  image.plot(rivermask) 
  points(x/nx, y/ny, pch="*", col="white")
  
  #Walk downriver
  path_mat=matrix(0, nrow=nx, ncol=ny)
  pathlist=NULL
  #pathlist=c(1, x, y,area[x,y], slopesUW$direction[x,y], slopesUW$slopex[x,y], slopesUW$slopey[x,y],subbasin$segments[x,y], subbasin$subbasins[x,y])
  
  #walk down the flow path
  step=1
  while(x<=nx & x>0 & y<=ny & y>0){
    path_mat[x,y]=1
    if(step==1){
      pathlist=c(1, x, y,area[x,y], slopesUW$direction[x,y], slopesUW$slopex[x,y], slopesUW$slopey[x,y],subbasin$segments[x,y], subbasin$subbasins[x,y])
    }else{
      pathlist=rbind(pathlist, c(step, x, y,area[x,y], slopesUW$direction[x,y], slopesUW$slopex[x,y], slopesUW$slopey[x,y], subbasin$segments[x,y], subbasin$subbasins[x,y]))
    }
    dir=slopesUW$direction[x,y]
    if(dir==1){y=y-1}
    if(dir==2){x=x-1}
    if(dir==3){y=y+1}
    if(dir==4){x=x+1}
    step=step+1
    #print(step)
  }
  #colnames(pathlist)=c("Step", "X", "Y", "Area", "Direction", "Slopex", "Slopey", "Subbbasin_Segment", "Subbasin_Subbasin")
  #image.plot(path_mat)
  
  #plot(pathlist[,1], pathlist[,4], type="l", xlab="step", ylab="area")
  
  #Plot
  nstep=length(pathlist)/9  ##CHANGE IF YOU ADD VARIBLES TO THE PATHLIST!!
  #only plot if thre is more than one step and you are starting in a numbered subbasin
  if(nstep>1 & subbasin$subbasins[headlist[h,1],headlist[h,2]]>0){
    slopetemp=rep(0, nstep)
    for(i in 1:nstep){slopetemp[i]=max(abs(pathlist[i,6]), abs(pathlist[i,7]))}
    
    fout=paste(outdir, runname, "Stream_Trace", h, ".pdf", sep="")
    pdf(fout, width=14, height=6)
    par(mfrow=c(1,3))
    plot(pathlist[,1], slopetemp, type="l", xlab="step", ylab="slope")
    plot(pathlist[,1], pathlist[,8], type="l", xlab="step", ylab="Subbasin Number")
    
    temp=rivermask  
    temp[temp==0]=NA
    temp=temp+path_mat
    maskcol=colorRampPalette(c('black', 'white'))
    #maskcol=colorRampPalette(c('white', 'white'))
    image.plot(subbasin$subbasins)
    image.plot((temp*2), add=T, col=maskcol(2), legend=F)
    dev.off()
  }
}

# #Testing find and replace
# nxt=nyt=10000
# test=matrix(floor(10*abs(rnorm(nx*ny, 0,1))), nrow=nx)
# search=16
# replace=34
# ilist=which(test==search)
# search[ilist]=replace


# #Identify nodes where area jumpes more than 2
# par(mfrow=c(2,3))
# xrange=60:75
# yrange=60:75
# xpoint=7
# ypoint=9
# segtemp=subbasin$segments[xrange,yrange]
# segtemp[which(segtemp==0)]=NA
# image.plot(rivermask[xrange,yrange])
# points(xpoint/15,ypoint/15, pch='*', col="white")
# image.plot(subbasin$subbasins[xrange,yrange],zlim=c(69,75))
# points(xpoint/15,ypoint/15, pch='*', col="white")
# image.plot(segtemp,zlim=c(69,75))
# points(xpoint/15,ypoint/15, pch='*', col="white")
# image.plot(slopesUW$direction[xrange,yrange])
# points(xpoint/15,ypoint/15, pch='*', col="white")
# sxtemp=slopesUW$slopex[xrange,yrange]
# sxtemp[which(sxtemp==0)]=NA
# image.plot(sxtemp, zlim=c(-0.3, 0.3))
# points(xpoint/15,ypoint/15, pch='*', col="white")
# sytemp=slopesUW$slopey[xrange,yrange]
# sytemp[which(sytemp==0)]=NA
# image.plot(sytemp, zlim=c(-0.3, 0.3))
# points(xpoint/15,ypoint/15, pch='*', col="black")

### making Sabino mask

#mask=subbasin$subbasins
#mask[mask==(1)]=0
#mask[mask==(20)]=1
#mask[mask==(19)]=1
#mask[mask==(16)]=1
#mask[mask==(15)]=1
#mask[mask==(14)]=1
#mask[mask==(13)]=1
#mask[mask==(12)]=1
#mask[mask==(11)]=1
#mask[mask==(10)]=1
#mask[mask==(7)]=1
#mask[mask==(6)]=1
#mask[mask==(21)]=0
#mask[mask==(18)]=0
#mask[mask==(17)]=0
#mask[mask==(9)]=0
#mask[mask==(8)]=0
#mask[mask==(5)]=0
#mask[mask==(4)]=0
#mask[mask==(3)]=0
#mask[mask==(2)]=0
#image.plot(mask)

#######plotting subbasin averaged slopes, actual river elevations, and the stream trace#######
outdir=("/Users/katie/Dropbox/ParFlow Domain/slopes/")
rivermask=rivers #This should be the river mask you input to the slope processing. could also be 'rivers' if you calcualted your own river mask for that

#Find all the headwater cells to use as starting points
#calculate the number of river cells draining to any cell
d4=c(1,2,3,4)
down=up=left=right=matrix(0, nrow=nx, ncol=ny) 
down[which(slopesUW$direction==d4[1])]=1
left[which(slopesUW$direction==d4[2])]=1
up[which(slopesUW$direction==d4[3])]=1
right[which(slopesUW$direction==d4[4])]=1
draincount=matrix(0, nrow=nx, ncol=ny)
draincount[,1:(ny-1)]=draincount[,1:(ny-1)]+down[,2:ny]*rivermask[,2:ny]
draincount[,2:ny]=draincount[,2:ny]+up[,1:(ny-1)]*rivermask[,1:(ny-1)]
draincount[1:(nx-1),]=draincount[1:(nx-1),]+left[2:nx,]*rivermask[2:nx,]
draincount[2:nx, ]=draincount[2:nx,]+right[1:(nx-1),]*rivermask[1:(nx-1),]

#Identify all the headwater cells
headwater=matrix(0, nrow=nx, ncol=ny)
headwater[which(draincount==0 & rivermask==1)]=1
#image.plot(headwater)
headlist=which(headwater==1, arr.ind=T)
colnames(headlist)=c("x", "y")
#headlist
nhead=nrow(headlist)

#Pick a point and walk downriver to the outlet
for(h in 1:nhead){
  #h=40 #picking one of the headwater cells
  x=headlist[h,1]
  y=headlist[h,2]
  
  #Plot to show you where you picked
  image.plot(rivermask) 
  points(x/nx, y/ny, pch="*", col="white")
  
  #Walk downriver
  path_mat=matrix(0, nrow=nx, ncol=ny)
  pathlist=NULL
  #pathlist=c(1, x, y,area[x,y], slopesUW$direction[x,y], slopesUW$slopex[x,y], slopesUW$slopey[x,y],subbasin$segments[x,y], subbasin$subbasins[x,y])
  
  #walk down the flow path
  step=1
  while(x<=nx & x>0 & y<=ny & y>0){
    path_mat[x,y]=1
    if(step==1){
      pathlist=c(1, x, y, slopesUW$direction[x,y], slopesUW$slopex[x,y], slopesUW$slopey[x,y],subbasin$segments[x,y], subbasin$subbasins[x,y], demT[x,y])
    }else{
      pathlist=rbind(pathlist, c(step, x, y, slopesUW$direction[x,y], slopesUW$slopex[x,y], slopesUW$slopey[x,y], subbasin$segments[x,y], subbasin$subbasins[x,y],demT[x,y]))
    }
    dir=slopesUW$direction[x,y]
    if(dir==1){y=y-1}
    if(dir==2){x=x-1}
    if(dir==3){y=y+1}
    if(dir==4){x=x+1}
    step=step+1
    #print(step)
  }
  #colnames(pathlist)=c("Step", "X", "Y", "Area", "Direction", "Slopex", "Slopey", "Subbbasin_Segment", "Subbasin_Subbasin")
  #image.plot(path_mat)
  
  #plot(pathlist[,1], pathlist[,4], type="l", xlab="step", ylab="area")
  
  #Plot
  nstep=length(pathlist)/9  ##CHANGE IF YOU ADD VARIBLES TO THE PATHLIST!!
  #only plot if thre is more than one step and you are starting in a numbered subbasin
  if(nstep>1 & subbasin$subbasins[headlist[h,1],headlist[h,2]]>0){
    slopetemp=rep(0, nstep)
    for(i in 1:nstep){slopetemp[i]=max(abs(pathlist[i,5]), abs(pathlist[i,6]))}
    
    fout=paste(outdir, runname, "Stream_Trace", h, ".pdf", sep="")
    pdf(fout, width=14, height=6)
    par(mfrow=c(1,3))
    plot(pathlist[,1], slopetemp, type="l", xlab="step", ylab="slope")
    plot(pathlist[,1], pathlist[,9], type="l", xlab="step", ylab="elevation")
    
    temp=rivermask  
    temp[temp==0]=NA
    temp=temp+path_mat
    maskcol=colorRampPalette(c('black', 'white'))
    #maskcol=colorRampPalette(c('white', 'white'))
    image.plot(subbasin$subbasins)
    image.plot((temp*2), add=T, col=maskcol(2), legend=F)
    dev.off()
  }
}


