##' @title Spatially continuous sampling
##' @description Draws a spatially continous sample of locations within a polygonal sampling region according to an \bold{"inhibitory plus close pairs"} specification.
##' @param poly a \code{sf} or \code{sp} polygon in which to generate the design.
##' @param size a non-negative integer giving the total number of locations to be sampled.
##' @param delta minimum permissible distance between any two locations in preliminary sample. This can be allowed to vary with the number of \code{'close pairs'} if a \bold{simple inhibitory} design is compared to one of the \bold{inhibitory plus close pairs} design.
##' @param delta.fix 'logical' specifies whether \code{delta} is fixed or allowed to vary with number of close pairs \eqn{k}. Default is \code{delta.fix = FALSE}.
##' @param k number of locations in preliminary sample to be replaced by near neighbours of other preliminary sample locations to form \code{close pairs} (integer between 0 and \code{size/2}). A \bold{simple inhibitory} deisgn is generated when \eqn{k = 0}.
##' @param rho maximum distance between the two locations in a \code{'close-pair'}.
##' @param ntries number of rejected proposals after which the algorithm will terminate.
##' @param plotit 'logical' specifying if graphical output is required. Default is \code{plotit = TRUE}.
##'
##' @details  To draw a simple inhibitory (\bold{SI}) sample of size \code{n}  from a spatially continuous region \eqn{A}, with the property that the distance between any two sampled locations is at least \code{delta}, the following algorithm is used.
##' \itemize{
##' \item{Step 1.} Set \eqn{i  = 1} and generate a point \eqn{x_{1}}  uniformly distributed on \eqn{{\cal D}}.
##' \item{Step 2.} Generate a point \eqn{x}  uniformly distributed on \eqn{{\cal D}} and calculate the minimum, \eqn{d_{\min}}, of the distances from \eqn{x_{i}} to all \eqn{x_{j}: j \leq i }.
##' \item{Step 3.} If \eqn{d_{\min} \ge \delta}, increase \eqn{i}  by 1, set \eqn{x_{i} = x} and return to step 2 if \eqn{i \le n}, otherwise stop;
##' \item{Step 4.} If \eqn{d_{\min} < \delta}, return to step 2 without increasing \eqn{i}.
##' }
##'
##' \bold{Sampling close pairs of points.}
##'
##' For some purposes, it is desirable that a spatial sampling scheme include pairs of closely spaced points, resulting in an inhibitory plus close pairs (\bold{ICP}) design. In this case, the above algorithm requires the following additional steps to be taken.
##' Let \eqn{k}  be the required number of close pairs. Choose a value \code{rho}  such that a close pair  of points will be a pair of points separated by a distance of at most \code{rho}.
##' \itemize{
##' \item{Step 5.} Set \eqn{j  = 1} and draw a random sample of size 2 from integers \eqn{1, 2, \ldots, n}, say \eqn{(i_1, i_2)};
##' \item{Step 6.} Replace \eqn{x_{i_{1}}} by \eqn{x_{i_{2}} + u} , where \eqn{u}  is uniformly distributed on the disc with centre \eqn{x_{i_{2}}} and radius \code{rho}, increase \eqn{i} by 1 and return to step 5 if \eqn{i \le k}, otherwise stop.
##' }
##'
##' When comparing a \bold{SI} design to one of the \bold{ICP} designs, the inhibitory components should have the same degree of spatial regularity.
##' This requires \eqn{\delta} to become a function of \eqn{k} namely \deqn{\delta_{k} = \delta_{0}\sqrt{n/(n - k)}} with \eqn{\delta_{0}} held fixed.
##'
##' @return a list with the following four components:
##' @return \code{size:} the total number of sampled locations.
##' @return \code{delta:} the value of \eqn{\delta} after taking into account the number of close pairs \eqn{k}. If \code{delta.fix = TRUE}, this will be \eqn{\delta} input by the user.
##' @return \eqn{k:} the number of close pairs included in the sample (for \bold{inhibitory plus close pairs} design).
##' @return \code{sample.locs:} a \code{sf} or \code{sp} object containing coordinates of dimension \code{n} by 2 containing the sampled locations.
##'
##' @note If \code{'delta'} is set to 0, a completely random sample is generated. In this case, \code{'close pairs'} are not permitted and \code{rho} is irrelevant.
##'
##' @seealso \code{\link[geosample:random.sample]{random.sample}} and \code{\link[geosample:discrete.inhibit.sample]{discrete.inhibit.sample}}
##'
##' @references Chipeta  M G, Terlouw D J, Phiri K S and Diggle P J. (2016b). Inhibitory geostatistical designs for spatial prediction taking account of uncertain covariance structure, \emph{Enviromentrics}, pp. 1-11.
##'
##' @author Michael G. Chipeta \email{mchipeta@@mlw.mw}
##' @author Peter J. Diggle \email{p.diggle@@lancaster.ac.uk}
##'
##' @examples
##' library("geoR")
##' library("sf")
##' data("parana")
##' poly <- parana$borders
##' poly <- matrix(c(poly[,1],poly[,2]),dim(poly)[1],2,byrow=FALSE)
##' #convert matrix to polygon
##' poly <- st_sf(st_sfc(st_polygon(list(as.matrix(poly)))))
##' #poly <- as(poly, "Spatial")
##' poly
##'
##' # Generate spatially regular sample
##' set.seed(5871121)
##' xy.sample1 <- contin.inhibit.sample(poly=poly,size = 100, delta = 30, plotit = TRUE)
##'
##'
##' # Generate spatially regular sample with 10 close pairs
##' set.seed(5871122)
##' xy.sample2 <- contin.inhibit.sample(poly,size = 100, delta = 30,
##'                                     k = 5, rho = 15, plotit = TRUE)
##'
##' # Generate spatially regular sample with 10 close pairs
##' set.seed(5871123)
##' xy.sample3 <- contin.inhibit.sample(poly,size = 100, delta = 30, delta.fix = TRUE,
##'                                     k = 10, rho = 15, plotit = TRUE)
##'
##' @import sp
##' @import sf
##' @importFrom splancs csr
##' @export

contin.inhibit.sample<-function(poly,size,delta, delta.fix = FALSE,
                                k=0,rho=NULL, ntries = 10000, plotit = TRUE) {

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
    poly <- sf::st_as_sf(poly)
  } else {
    poly <- poly
  }
  if (length(size) > 0){
    if (!is.numeric(size) | size <= 0)
      stop("\n 'size' must be a positive integer")
    else
      orig.size <- size
  }
  if(length(delta) > 0){
    if(!is.numeric(delta) | delta < 0)
      stop("\n 'delta' must be a positive integer >= 0")
    if(delta == 0 && k > 0)
        stop("\n Close pairs not allowed for completely
             random sample (i.e. when 'delta' = 0)")
    if(delta == 0 && k == 0)
      rho = NULL
  }
  if(length(k) > 0){
    if(!is.numeric(k) | k < 0)
      stop("\n 'k' must be a positive integer >= 0")
    if (k > size/2)
      stop("\n 'k' must be between 0 and 'size'/2")
    if(k>0 && is.null(rho)){
        stop("\n 'rho' must be provided if 'k' > 0")
    }
    if(k>0 && rho <= 0){
      stop("\n 'rho' must be positive,
           between > 0 and 'delta'/2")
    }
  }
  if(length(rho) > 0){
    if(!is.numeric(rho) | rho < 0)
      stop("\n 'rho' must be positive")
    if(rho > delta/2)
      stop("\n 'rho' must be between > 0
           and 'delta'/2")
  }

  st.poly <- sf::st_coordinates(poly)[,c(1:2)]
  xy.sample <- matrix(csr(st.poly,1),1,2)
  if(delta == 0){
    for (i in 2:size) {
        xy.try <- c(csr(st.poly,1))
        xy.sample <- rbind(xy.sample, xy.try)
      }
  } else {
      if (delta.fix == TRUE){
        delta = delta
      } else {
        delta <- delta * sqrt(size/(size - k))
      }
    dsq  <- delta*delta
    if (!is.infinite(size) && (size * pi * dsq/4 > st_area(poly)))
      stop("\n Polygon is too small to fit ", size,
           " points, with 'k' = ", k, " close pairs,",
           " at minimum separation ", round(delta, digits = 4))
    while (dim(xy.sample)[1] < size) {
      dmin<-0
      iter <- 1
      while (dmin<dsq) {
        xy.try<-c(csr(st.poly,1))
        dmin<-min((xy.sample[,1]-xy.try[1])^2+(xy.sample[,2]-xy.try[2])^2)
        iter <- iter + 1
        if(iter == ntries)
          break
      }
      xy.sample<-rbind(xy.sample,xy.try)
      if(iter == ntries && dim(xy.sample)[1] < size){
        warning("\n For the given 'delta' and 'size', only ", dim(xy.sample)[1],
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
    for (j in 1:k) {
      take1<-take[j,1]; take2<-take[j,2]
      xy1<-c(xy.sample[take1,])
      angle<-2*pi*runif(1)
      radius<-rho*sqrt(runif(1))
      xy.sample[take2,]<-xy1+radius*c(cos(angle),sin(angle))
    }
  }
  xy.sample <- xy.sample %>%
    as.data.frame %>%
    sf::st_as_sf(coords = c(1,2))
  sample.locs <- sf::st_as_sf(xy.sample)
  if(class(sample.locs)[1] != class(poly.origin)[1]){
    sample.locs <- sf::as_Spatial(xy.sample, "Spatial")
  }

  if(plotit==TRUE){
    par(oma=c(5, 5, 5, 5.5), mar=c(5.5, 5.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)
    plot(st_geometry(xy.sample),pch=19,col=1,axes = TRUE,
         xlab="longitude",ylab="lattitude", font.main = 3,
         cex.main = 1.2, col.main = "blue",
         main = paste("Continuous sampling design,", k,
                      "close pairs", sep = " "),
         xlim = c(range(st.poly[,1])),
         ylim = c(range(st.poly[,2])))
    plot(st_geometry(poly), add= TRUE)
  }
  res <- list()
  res$size <- dim(unique(xy.sample))[1]
  res$delta = delta
  res$k <- k
  res$sample.locs = sample.locs

  return(res)
}
