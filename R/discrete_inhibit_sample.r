##' @title Spatially discrete sampling
##' @description Draw a spatially discrete sample from a specified set of spatial locations within a polygonal sampling region according to an \bold{"inhibitory plus close pairs"} specification.
##' @param obj a \code{sf} or \code{sp} object where each line corresponds to a spatial location containing values of two-dimensional coordinates and, optionally, the values of one or more associated values, typically an outcome of interest and any associated covariates.
##' @param size a non-negative integer giving the total number of locations to be sampled.
##' @param delta minimum permissible distance between any two locations in preliminary sample. This can be allowed to vary with number of \code{'close pairs'} if a \bold{simple inhibitory} design is compared to one of the \bold{inhibitory plus close pairs} design.
##' @param delta.fix 'logical' specifies whether \code{'delta'} is fixed or allowed to vary with number of close pairs \eqn{k}. Default is \code{delta.fix = FALSE}.
##' @param k number of close-pair locations in the sample. Must be an integer between 0 and \code{size}/2.
##' @param cp.criterion criterion for choosing close pairs \eqn{k}. The \code{"cp.zeta"} criterion chooses locations not included in the initial sample, from the uniform distribution of a disk with radius \code{'zeta'} (NB: \code{zeta} argument must be provided for this criterion). The \code{"cp.neighb"} criterion chooses nearest neighbours amongst locations not included in the initial sample (\code{'zeta'} becomes trivial for \code{'cp.neighb'} criterion).
##' @param zeta maximum permissible distance (radius of a disk with center \eqn{x^{*}_{j}, j = 1, \ldots, k}) within which a close-pair point is placed. See \bold{Details}.
##' @param ntries number of rejected proposals after which the algorithm terminates.
##' @param poly 'optional', a \code{sf} or \code{sp} polygon object in which the design sits. The default is the bounding box of points given by \code{obj}.
##' @param plotit 'logical' specifying if graphical output is required. Default is \code{plotit = TRUE}.
##'
##' @details To draw a sample of size \eqn{n} from a population of spatial locations \eqn{X_{i} : i  = 1,\ldots,N}, with the property that the distance between any two sampled locations is at least \eqn{\delta}, the function implements the following algorithm.
##' \itemize{
##' \item{Step 1.} Draw an initial sample of size \eqn{n}  completely at random and call this \eqn{x_{i}  : i  = 1,\dots, n}.
##' \item{Step 2.} Set \eqn{i  = 1}.
##' \item{Step 3.} Calculate the smallest distance, \eqn{d_{\min}}, from \eqn{x_{i}}  to all other \eqn{x_{j}}  in the initial sample.
##' \item{Step 4.} If \eqn{d_{\min} \ge \delta}, increase \eqn{i}  by 1 and return to step 2 if \eqn{i \le n}, otherwise stop.
##' \item{Step 5.} If \eqn{d_{\min} < \delta}, draw an integer \eqn{j}  at random from \eqn{1,  2,\ldots,N}, set \eqn{x_{i}  = X_{j}}  and return to step 3.}
##'
##' Samples generated in this way exhibit  more regular spatial arrangements than would random samples of the same size. The degree of regularity achievable will be influenced by the spatial arrangement of the population \eqn{X_{i}  : i  = 1,\ldots,N}, the specified value of \eqn{\delta}  and the sample size \eqn{n}. For any given population, if \eqn{n}  and/or \eqn{\delta} is too large, a sample of the required size with the distance between any two sampled locations at least \eqn{\delta} will not be achievable; the algorithm will then find \eqn{n_{s} < n} points that can be placed for the given parameters.
##'
##' \bold{Sampling close pairs of points.}
##'
##' For some purposes, typically when using the same sample for parameter estimation and spatial prediction, it is desirable that a spatial sampling scheme include pairs of closely spaced points \eqn{x}. The function offers two ways of specifying close pairs, either as the closest available unsampled point to an existing sampled point \code{(cp.critetrion = cp.neighb)}, or as a random choice from amongst all available unsampled points within distance \eqn{zeta} of an existing sampled point \code{(cp.criterion = cp.zeta)}.
##' The algorithm proceeds as follows.
##'
##' Let \eqn{k} be the required number of close pairs.
##' \itemize{
##' \item{Step 1.} Construct a simple inhibitory design \bold{SI}\eqn{(n - k, \delta)}.
##' \item{Step 2.} Sample \eqn{k} from \eqn{x_{1}, \ldots, x_{n - k}} without replacement and call this set \eqn{x_{j} : j = 1, \ldots, k}.
##' \item{Step 3.} For each \eqn{x_{j}: j = 1, \ldots, k}, select a close pair \eqn{x_{n-k+j}} according to the specified criterion.}
##'
##' \bold{Note:} Depending on the spatial configuration of potential sampling locations and, when the selection criterion \code{cp.criterion = cp.zeta}, the specified value of \eqn{zeta}, it is possible that one or more of the selected points  \eqn{x_{j}} in Step 2 will not have an eligible ``close pair''. In this case, the algorithm will try  find an alternative \eqn{x_{j}} and report a warning if it fails to do so.
##'
##' @return a list with the following four components:
##' @return \code{unique.locs:} the number of unique sampled locations.
##' @return \code{delta:} the value of \eqn{\delta} after taking into account the number of close pairs \eqn{k}. If \code{delta.fix = TRUE}, this will be \eqn{\delta} input by the user.
##' @return \eqn{k:} the number of close pairs included in the sample (for \bold{inhibitory plus close pairs} design).
##' @return \code{sample.locs:} a \code{sf} or \code{sp} object containing the final sampled locations and any associated values.
##'
##' @note If \code{'delta'} is set to 0, a completely random sample is generated. In this case, \emph{'close pairs'} are not permitted and \code{'zeta'} becomes trivial.
##'
##' @references Chipeta  M G, Terlouw D J, Phiri K S and Diggle P J. (2016). Inhibitory geostatistical designs for spatial prediction taking account of uncertain covariance structure, \emph{Enviromentrics}, pp. 1-11.
##' @references Diggle P J. (2014). \emph{Statistical Analysis of Spatial and Spatio-Temporal Point Patterns.} 3rd ed., Boca Raton: CRC Press
##' @references Diggle P J and Lophaven S. (2006). Bayesian geostatistical design, \emph{Scandinavian Journal of Statistics} \bold{33}(1) pp. 53 - 64.
##'
##' @author Michael G. Chipeta \email{mchipeta@@mlw.mw}
##' @author Peter J. Diggle \email{p.diggle@@lancaster.ac.uk}
##'
##'
##' @examples
##' library("sf")
##' set.seed(1234)
##' x <- 0.015+0.03*(1:33)
##' xall <- rep(x,33)
##' yall <- c(t(matrix(xall,33,33)))
##' xy <- cbind(xall,yall)+matrix(-0.0075+0.015*runif(33*33*2),33*33,2)
##'
##'
##' # Convert to SF object
##' xy <- xy %>%
##'   as.data.frame %>%
##'   sf::st_as_sf(coords = c(1,2))
##'
##'
##' # Plot the points
##' plot(st_geometry(xy),pch=19,cex=0.25,xlab="longitude",ylab="latitude",
##'      cex.lab=1,cex.axis=1,cex.main=1, axes = TRUE)
##'
##'
##' # Generate spatially random sample
##' set.seed(15892)
##' xy.sample1 <- xy[sample(1:dim(xy)[1],50,replace=FALSE),]
##' plot(xy.sample1, pch = 19, col = 'black', add = TRUE)
##'
##'
# Generate spatially simple inhibitory sample
##' set.seed(15892)
##' xy.sample2 <- discrete.inhibit.sample(obj=xy,size = 100,
##'                                       delta = 0.08,plotit = TRUE)
##' plot(st_geometry(xy),pch=19, cex = 0.25, col="black", add = TRUE)
##'
##'
##' # Generate spatially inhibitory sample
##' # with close pairs (cp.zeta criterion):
##' set.seed(15892)
##' xy.sample3 <- discrete.inhibit.sample(obj=xy, size = 100,delta = 0.065,
##'                                      k = 25,cp.criterion = "cp.zeta",
##'                                      zeta = 0.025, plotit = TRUE)
##' plot(st_geometry(xy),pch=19, cex = 0.25, col="black", add = TRUE)
##'
##'
##' # Generate spatially inhibitory sample
##' # with close pairs (cp.neighb criterion):
##' set.seed(15892)
##' xy.sample4 <- discrete.inhibit.sample(obj=xy,size = 100,
##'                                       delta = 0.065, k = 25,cp.criterion = "cp.neighb",
##'                                       plotit = TRUE)
##' plot(st_geometry(xy),pch=19, cex = 0.25, col="black", add = TRUE)
##'
##'
##' # Generate spatially inhibitory sample
##' # with close pairs (cp.zeta criterion):
##' set.seed(15892)
##' xy.sample5 <- discrete.inhibit.sample(obj=xy,size = 100,
##'                                       delta = 0.065, cp.criterion = "cp.zeta",
##'                                       zeta = 0.025, delta.fix = TRUE,
##'                                       k = 25, plotit = TRUE)
##' plot(st_geometry(xy),pch=19, cex = 0.25, col="black", add = TRUE)
##'
##'
##' # Generate simple inhibitory sample from a regular grid
##' library("PrevMap")
##' data("sim.data")
##' set.seed(15892)
##' xy.sample6 <- discrete.inhibit.sample(obj = sim.data,
##'                                       size = 50, delta = 0.08,plotit = TRUE)
##' plot(st_geometry(sim.data),pch=19,col="black", cex = 0.25, add = TRUE)
##'
##'
##' # Generate inhibitory plus close pairs sample from a regular grid
##' set.seed(15892)
##' xy.sample7 <- discrete.inhibit.sample(obj = sim.data,
##'                                       cp.criterion = "cp.neighb", size = 50,
##'                                       delta = 0.1, k = 5, plotit =TRUE)
##' plot(st_geometry(sim.data),pch=19,col="black", cex = 0.25, add = TRUE)
##' @import sf
##' @import sp
##' @export

discrete.inhibit.sample  <- function(obj, size, delta, delta.fix = FALSE, k = 0, cp.criterion = NULL,
                                      zeta, ntries = 10000, poly = NULL, plotit = TRUE)
{

  obj.origin <- obj
  if(!inherits(obj, 'SpatialPointsDataFrame')){
    if(!inherits(obj, 'SpatialPoints')){
      if(!inherits(obj,"sf") & !inherits(obj, "data.frame")){
        stop("\n 'obj' must be of class 'sp' or 'sf'")
      }
    }
  }
  if(inherits(obj, 'Spatial')){
    obj <- sf::st_as_sf(obj)
  }
  if (any(!is.numeric(sf::st_coordinates(obj))))
    stop("\n non-numerical values in the coordinates")
  if(any(is.na(sf::st_geometry(obj)))){
    warning("\n NA's not allowed in 'obj' coordinates")
    obj <- obj[complete.cases(obj), , drop = FALSE]
    warning("\n eliminating rows with NA's")
  }
  if(!is.null(poly)){
    poly.origin <- poly
    if(!inherits(poly, 'SpatialPolygonsDataFrame'))
    if(!inherits(poly, 'SpatialPolygons'))
    if(!inherits(poly, 'Polygons'))
    if(!inherits(poly, 'Polygon'))
    if(!inherits(poly, 'sfc_POLYGON'))
    if(!inherits(poly, 'sfc'))
      if(!inherits(poly, 'sf'))
        stop("\n 'poly' must be of class 'sp' or 'sf'")
  }
  if (inherits(poly, 'Spatial')){
    poly <- st_as_sf(poly)
  } else{
    poly <- poly
    if (!identical(st_crs(obj), st_crs(poly)))
      stop("\n 'obj' and 'poly' are not in the same coordinate system")
  }
  if(length(size) > 0){
    if(!is.numeric(size) | size <= 0)
      stop("\n 'size' must be a positive integer")
    else
      orig.size <- size
  }
  if(length(k) > 0){
    if(k > 0){
      if(!is.numeric(k) | k < 0)
        stop("\n 'k' must be a positive integer >= 0")
      if (k > size/2)
        stop("\n 'k' must be between 0 and size/2")
      if(is.null(cp.criterion))
        stop("\n Close pairs selection criterion 'cp.criterion' must be provided")
      if (cp.criterion != "cp.zeta" & cp.criterion != "cp.neighb")
        stop("\n 'cp.criterion' must be either 'cp.neighb' or 'cp.zeta'")
    }
  }
  if(length(delta) > 0){
    if(!is.numeric(delta) | delta < 0)
      stop("\n 'delta' must be a positive integer >= 0")
  }
  if(delta == 0){
    if(k>0){
      stop("\n close pairs not allowed for completely random sample")
    } else {
      res1 <- as.data.frame(unique(st_coordinates(obj)))
      N   <- dim(res1)[1]
      index  <- 1:N
      index.sample  <- sample(index, size, replace = FALSE)
      xy.sample  <- res1[index.sample,]; dim(xy.sample) #to remove this
    }
  } else {
    delta.orig <- delta
    if (delta.fix == TRUE){
      delta = delta
    } else {
      delta <- delta * sqrt(size/(size - k));delta
    }
    dsq  <- delta*delta; dsq
    if (is.null(poly)){
      poly.shape <- sf::st_convex_hull(st_union(obj))
    } else {
      poly.shape <- poly
    }
    if(!is.infinite(size) && (size * pi * dsq/4 > as.numeric(sf::st_area(poly.shape))))
      stop("\n Polygon is too small to fit ", size, " points, with 'k' = ", k, " close pairs,",
           " at minimum separation ", round(delta, digits = 4))

    # "XnotinF" Function
    xnotiny  <- function(a1,a2)
    {
      a1.vec <- apply(a1, 1, paste, collapse = "")
      a2.vec <- apply(a2, 1, paste, collapse = "")
      a1.without.a2.rows <- as.data.frame(a1[!a1.vec %in% a2.vec,])
      return(a1.without.a2.rows)
    }

    res1 <- as.data.frame(unique(sf::st_coordinates(obj)))
    N   <- dim(res1)[1]
    index  <- 1:N
    index.sample  <- sample(index, 1, replace = FALSE)
    xy.sample  <- res1[index.sample,]
    for (i in 2:size){
      dmin  <- 0
      iter <- 1
      while (dmin < dsq){
        take <- sample(index, 1)
        iter <- iter+1
        dvec<-(res1[take,1]-xy.sample[,1])^2+(res1[take,2]-xy.sample[,2])^2
        dmin<-min(dvec)
        if(iter == ntries)
          break
      }
      xy.sample[i,]  <- res1[take,]
      num <- dim(xy.sample)[1]
      if(iter == ntries && dim(xy.sample)[1] < size){
        warning("\n For the given 'delta' and 'size', only ", num,
                " inhibitory sample locations placed out of ", size,
                ". Consider revising 'delta' and/or 'size'")
        break
      }
    }
  }

  if (k>0) {
    k.origin <- k
    size <- dim(unique(xy.sample))[1]
    reduction <- ((orig.size - size)/orig.size)
    if (k > size/2){
      k <- floor(k*(1-reduction))
      warning("\n For the given parameters, only ", k,
              " close pairs could be placed out of ", k.origin)
    }
    take<-matrix(sample(1:size,2*k,replace=FALSE),k,2)
    xy.sample <- unique(xy.sample)
    if(cp.criterion == "cp.neighb"){
      for (j in 1:k) {
        take1<-take[j,1]; take2<-take[j,2]
        xy1<-as.numeric(c(xy.sample[take1,]))
        dvec<-(res1[,1]-xy1[1])^2+(res1[,2]-xy1[2])^2
        neighbour<-order(dvec)[2]
        xy.sample[take2,]<-res1[neighbour,]
      }
    }
    if(cp.criterion == "cp.zeta"){
      if(!is.numeric(zeta) | zeta < 0)
        stop("\n 'zeta' must be between > 0 and 'delta'/2")
      if(zeta < delta.orig*0.005)
        stop("\n 'zeta' too small.")
      if(zeta > delta.orig/2)
        stop("\n 'zeta' must be between > 0 and 'delta'/2")
      for (j in 1:k){
        take1<-take[j,1]; take2<-take[j,2]
        xy1<-as.numeric(c(xy.sample[take1,]))
        dvec<-(res1[,1]-xy1[1])^2+(res1[,2]-xy1[2])^2
        z.vec <- which(dvec > 0 & dvec <= zeta*0.25)
        z.vec.pts <- (1:dim(res1)[1])[z.vec]
        avail.locs <- xnotiny(res1[z.vec,], xy.sample)
        if (nrow(avail.locs) > 0) {
            rep.loc <- sample(1:dim(avail.locs)[1],1,replace = F)
            xy.sample[take2,]<-avail.locs[rep.loc,]
          } else {
              warning("\n One or more locations do not have
                      eligible 'close pairs'")
              break
          }
        }
    }
  }

  xy.sample <- sf::st_as_sf(xy.sample, coords = c("X", "Y"))
  st_crs(xy.sample) <- st_crs(obj)
  xy.sample <- obj[xy.sample, ]

  if(plotit==TRUE){
    par(oma=c(5, 5, 5, 5.5), mar=c(5.5, 5.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)
    plot(st_geometry(xy.sample),pch=19,col=1,axes = TRUE,
         xlab="longitude",ylab="lattitude", font.main = 3,
         cex.main = 1.2, col.main = "blue",
         main = paste("Inhibitory sampling design,", k,
                      "close pairs", sep = " "),
         xlim = c(range(st_coordinates(poly.shape)[,1])),
         ylim = c(range(st_coordinates(poly.shape)[,2])))
    if (!is.null(poly)){
      plot(st_geometry(poly.shape), add= TRUE)
    }
  }

  res <- list()
  res$unique.locs <- num
  res$delta = delta
  res$k <- k
  res$sample.locs = xy.sample
  st_crs(res$sample.locs) <- st_crs(obj)
  if (class(xy.sample)[1] != class(obj.origin)[1]){
    res$sample.locs <- sf::as_Spatial(res$sample.locs, "Spatial")
  }

  return(res)
}



