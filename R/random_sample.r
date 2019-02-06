##' @title Spatially random sample
##' @description This function draws a spatially random sample from a discrete set of units located over some defined geographical region or generate completely spatially random points within a polygon.
##' @param obj a \code{sf} or \code{sp} object (with \eqn{N \geq \code{size}}) where each line corresponds to one spatial location. It should contain values of 2D coordinates, data and, optionally, covariate(s) value(s) at the locations. This argument must be provided when sampling from a \code{"discrete"} set of points, see \code{'type'} below for details.
##' @param type random sampling, a choice of either \code{"discrete"}, from a set of \eqn{N} potential sampling points or \code{"continuum"} from independent, compeletely random points.
##' @param size a non-negative integer giving the total number of locations to be sampled.
##' @param poly 'optional' a \code{sf} or \code{sp} polygon in which to generate the design. The default is the bounding box of points given by \code{obj}. When sampling from a \code{"continuum"}, the argument \code{'poly'} must be provided.
##' @param plotit 'logical' specifying if graphical output is required. Default is \code{plotit = TRUE}.
##'
##' @return a \code{sf} or \code{sp} object of dimension \eqn{n} by \code{p=dim(obj)[2]} containing the final sampled locations and any associated values, if sampling from a \code{"discrete"} set of points. A matrix of \eqn{n} by 2 containing sampled locations, if sampling from a \code{"continuum"}.
##'
##' @examples
##' # 1. Sampling from a discrete set of points.
##' library("dplyr")
##' x <- 0.015+0.03*(1:33)
##' xall <- rep(x,33)
##' yall <- c(t(matrix(xall,33,33)))
##' xy <- cbind(xall,yall)+matrix(-0.0075+0.015*runif(33*33*2),33*33,2)
##' colnames(xy) <- c('X','Y')
##'
##' # Convert to SF
##' xy <- xy %>%
##'   as.data.frame %>%
##'   sf::st_as_sf(coords = c(1,2))
##' xy <- sf::st_as_sf(xy, coords = c('X', 'Y'))
##'
##'
##' # Sampling from a discrete set.
##' set.seed(15892)
##' xy.sample <- random.sample(obj = xy, size = 100, type = "discrete", plotit = TRUE)
##'
##'
##' # Sampling from a continuum.
##' library("geoR")
##' data("parana")
##' poly <- parana$borders
##' poly <- matrix(c(poly[,1],poly[,2]),dim(poly)[1],2,byrow=FALSE)
##' # Convert matrix to polygon
##' poly <- st_sf(st_sfc(st_polygon(list(as.matrix(poly)))))
##'
##' set.seed(15892)
##' xy.sample <- random.sample(poly = poly,size = 100, type = "continuum", plotit = TRUE)
##'
##'
##'
##'
##' @author Michael G. Chipeta \email{mchipeta@@mlw.mw}
##' @author Peter J. Diggle \email{p.diggle@@lancaster.ac.uk}
##'
##' @references Rowlingson, B. and Diggle, P. 1993 Splancs: spatial point pattern analysis code in S-Plus. Computers and Geosciences, 19, 627-655
##'
##' @import sp
##' @import sf
##' @importFrom splancs csr
##' @export

random.sample <- function(obj = NULL, poly = NULL, type, size, plotit = TRUE)

{
  if (is.null(type)){
    stop("\n 'type' must be provided")
  }
  if (type != "discrete" & type != "continuum")
    stop("'type' must be either 'discrete' or 'continuum'")

  if (type == "discrete"){
    obj.origin <- obj
  if (is.null(obj))
    stop("\n'obj' must be provided")
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
      stop("\n non-numerical values in 'obj' coordinates")
    if(any(is.na(sf::st_coordinates(obj)))){
      warning("\n NA's not allowed in 'obj' coordinates")
      obj <- obj[complete.cases(st_coordinates(obj)), , drop = FALSE]
      warning("\n eliminating rows with NA's")
    }
    if (is.null(poly)){
      poly <- sf::st_convex_hull(sf::st_union(obj))
    }
    if (length(size) > 0){
      if (!is.numeric(size) | size <= 0)
        stop("\n 'size' must be a positive integer")
    }
    if (size >= dim(obj)[1])
      stop("\n 'size' must be less than the total
           number of locations to sample from")
    if(size == 1){
      xy.sample <- obj[sample(1:dim(obj)[1], size, replace = FALSE), ]
      xy.sample <- xy.sample[which(!duplicated(xy.sample$geometry)), ]
    } else {
      ctr <- 1
      while(ctr < size){
        xy.sample <- obj[sample(1:dim(obj)[1], size, replace = FALSE), ]
        xy.sample <- xy.sample[which(!duplicated(xy.sample$geometry)), ]
        ctr <- dim(xy.sample)[1]
      }
    }
    res <- xy.sample
    if(class(xy.sample)[1] != class(obj.origin)[1]){
      res <- sf::as_Spatial(xy.sample, "Spatial")
    }
  }


  if (type == "continuum") {
    if (is.null(poly)){
      stop("\n Provide polygon in which to generate sample points")
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
    if(inherits(poly, 'Spatial')){
      plot.poly <- sf::st_as_sf(poly)
    } else {
      plot.poly <- poly
    }

    st.poly <- sf::st_coordinates(plot.poly)[,c(1:2)]
    xy.sample <- matrix(csr(st.poly,1),1,2)
    for (i in 2:size) {
      xy.try <- c(csr(st.poly,1))
      xy.sample <- rbind(xy.sample, xy.try)
    }
    xy.sample <- xy.sample %>%
      as.data.frame %>%
      sf::st_as_sf(coords = c(1,2))
    res <- xy.sample <- sf::st_as_sf(xy.sample)
    if(class(xy.sample)[1] != class(poly.origin)[1]){
      res <- sf::as_Spatial(xy.sample, "Spatial")
    }
  }

  if(plotit==TRUE){
    par(oma=c(5, 5, 5, 5.5), mar=c(5.5, 5.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)
    if (type == "discrete"){
      plot(st_geometry(xy.sample), pch = 19, col = 1, axes = TRUE,
           xlab = "longitude", ylab = "lattitude", font.main = 3,
           cex.main = 1.2, col.main = "blue",
           main = paste("Random sampling design,", size, "points", sep = " "))
      if (class(obj.origin)[1] == "sf"){
        plot(st_geometry(obj.origin),pch=19, cex = 0.25, col="black", add = TRUE)
      }else{
        plot(obj.origin,pch=19, cex = 0.25, col="black", add = TRUE)
      }
    } else{
      plot(st_geometry(xy.sample),pch=19,col=1,axes = TRUE,
           xlab="longitude",ylab="lattitude", font.main = 3, cex.main = 1.2, col.main = "blue",
           main = paste("Random sampling design,", size, "points", sep = " "),
           xlim = c(range(st.poly[,1])),
           ylim = c(range(st.poly[,2])))
      plot(st_geometry(plot.poly), add= TRUE)
    }
  }

  return(res)
}




