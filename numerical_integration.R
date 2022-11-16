###Numerical integration###

########################## Functions ##############################################


grid <- function(r, center){
  #A function taking number of grid points, r, returning the points on the real line
  #that are spaced according to the rule 2/3 within three standard deviations from
  #the mean and the rest spaced logarithmically in the tails.
  
  c(center -3 -4*log(r/(1:(r-1))),   #i=1,...,r-1
    center -3 +3*(0:(4*r))/(2*r),    #i=r,...,5r
    center +3 +4*log(r/((r-1):1)))}  #i=5r+1,...,6r-1



grid_bounded <- function(points, bounds){
  #This function takes a set of grid points found using the function "grid" and
  #removes all points outside the bounds. The bounds are added into the vector to
  #become the first and last points.
  
  x <- max(points[1],bounds[1]) #first element is either first point or lower bound
  x <- c(x,points[which(points>x & points<bounds[2])])  #add in all points between bounds
  if(bounds[2]<points[length(points)])
    {x <- c(x, bounds[2])} #last element is either last point or upper bound
  return(x)
}



grid_z <- function(points){
  #This function calculates the z values - taking z1,z3,...,zm as the points in the 
  #bounded grid and the even z values, z2,z4,...,zm-1 as the midpoints of these.
  
  m <- 2*length(points)-1 #length of the vector to be created
  z1 <- rep(0, m)
  z1[seq(1,m,by=2)] <- points #z1,z3,...,zm
  z1[seq(2,m-1,by=2)] <- points[-length(points)] + diff(points)/2 #z2,z4,...,zm-1
  return(z1)}



weights <- function(points){
  #This function applies weights to the vector found using the grid functions. This
  #is done using Simposon's Rule.
  m <- length(points) #length of the vector
  w1 <- rep(0,m) #empty vector of weights
  if(length(points)==3){ #account for very small number of points where below does not work.
    h<-diff(points)[1]
    w1 <- c(h/6,4*h/6,h/6)}
  else{
    w1[1] <- (points[3]-points[1])/6 #w1
    w1[seq(3,m-2,by=2)] <- (points[seq(5,m,by=2)]-points[seq(1,m-4,by=2)])/6 #w3,w5,...,w(m-2)
    w1[seq(2,m-1,by=2)] <- 4*(points[seq(3,m,by=2)]-points[seq(1,m-2,by=2)])/6 #w2,w4,...,w(m-1)
    w1[m] <- (points[m]-points[m-2])/6 #wm
  }
  return(w1)}



