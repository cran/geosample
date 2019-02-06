## ----options, echo=FALSE----------------------------------------------------
options(prompt="R> ", continue="+  ", width=78, useFancyQuotes=FALSE, digits=4, message=FALSE)

## ----simpleinhibitory, include=TRUE, eval=TRUE, echo=TRUE, fig.height=9, fig.width=9, fig.cap="\\label{fig:simpleinhibitory} Simple inhibitory (discrete) design with $\\delta$ = 0.04 and $n_0$ = 30.", message = FALSE, fig.asp=1----
library("geosample")
library("viridisLite")
data(sim.data)
head(sim.data, n = 6L, addrownums = TRUE)
set.seed(123)
my.sample <- discrete.inhibit.sample(obj = sim.data, size = 30,
                                     delta = 0.04, plotit = TRUE)

## ----simestim, include=TRUE, eval=TRUE, echo=TRUE, fig.keep='none', message=FALSE----
library("PrevMap")
knots <- as.matrix(expand.grid(seq(-0.2, 1.2, length = 15),
                               seq(-0.2, 1.2, length = 15)))
mcmc.ctr <- control.mcmc.MCML(n.sim = 5500, burnin=500, thin = 5)
dat <- my.sample[[4]]
par0 <- c(0.001, 1, 0.4)
model.fit <-
  binomial.logistic.MCML(y ~ 1, units.m = ~units.m, data = dat, par0 = par0,
                         coords=~st_coordinates(dat), fixed.rel.nugget = 0,
                         start.cov.pars = par0[3], control.mcmc = mcmc.ctr,
                         low.rank = TRUE, knots = knots, kappa = 1.5,
                         method = "BFGS", messages = FALSE,
                         plot.correlogram = FALSE)

summary(model.fit, log.cov.pars = FALSE)

## ----simpred, include=TRUE, eval=TRUE, echo=TRUE----------------------------
model.pred <-
  spatial.pred.binomial.MCML(object = model.fit, type = "joint",
                             control.mcmc = mcmc.ctr,thresholds = 0.45,
                             grid.pred = st_coordinates(sim.data),
                             scale.predictions = "prevalence",
                             scale.thresholds ="prevalence",
                             standard.errors = TRUE, messages = FALSE,
                             plot.correlogram = FALSE)

## ----predvis, include = TRUE, eval=TRUE, echo=TRUE, fig.height=7, fig.width=7, fig.cap = "\\label{fig:predvis} Spatial prediction visualisation. Spatial predictions on the LHS and exceedance probabilities $P(x; 0.45)$ = $P$($prev$ > 0.45 at location $x$) on the RHS.", fig.asp=1----
par(mfrow = c(1,2))
plot(model.pred, type = "prevalence", col = viridis(256, direction = -1),
     summary = "predictions", zlim = c(0, 1))
contour(model.pred, type="prevalence", summary="predictions", zlim = c(0, 1),
        levels = seq(0.1,0.9, 0.1), add = TRUE)
plot(model.pred,summary="exceedance.prob",zlim=c(0,1),
     col = viridis(256, direction = -1))
contour(model.pred, summary = "exceedance.prob",zlim = c(0, 1),
        levels = seq(0.1,0.3, 0.1), add = TRUE)
par(mfrow = c(1,1))

## ----adaptivesample, include=TRUE, eval=TRUE, echo=TRUE, fig.height=9, fig.width=9, fig.cap = "\\label{fig:adaptivesample} Adaptive sampling design with $\\delta$ = 0.1 and $b = 10$, Dark blue dots ($n_0$ = 30) are the initial sampling locations. Red dots ($n_a$ = 10) are adaptive sampling locations added after analysing data from the initial design."----
obj.1 <- as.data.frame(cbind(model.pred$grid,
         c(model.pred$prevalence$standard.errors)^2))
colnames(obj.1) <- c("coord1", "coord2", "pred.var")
obj.1 <- sf::st_as_sf(obj.1, coords = c('coord1', 'coord2'))
adapt.sample.pv <-
  adaptive.sample(obj1 = obj.1, obj2 = dat,
                  pred.var.col = 1,criterion = "predvar",
                  delta = 0.1, batch.size = 10, poly = NULL,
                  plotit = TRUE)


## ----avail, include=TRUE, eval=TRUE, echo=TRUE, fig.height=9, fig.width=9,fig.cap="\\label{fig:majlocs} All potential household sampling locations in Majete.", fig.asp=1----
data("border")
data("majete")
plot(st_geometry(majete), pch = 19, cex = 0.5,
     xlim=range(st_coordinates(border)[,1]),
     ylim=range(st_coordinates(border)[,2]),
     axes = TRUE, xlab="longitude", ylab="latitude")
plot(border, lwd = 2, add= TRUE)

## ----applic1, include=TRUE, eval=TRUE, echo=TRUE, fig.height=10, fig.width=10, fig.cap="\\label{fig:majeteSI} Simple inhibitory (discrete) design with $\\delta$ = 400 meters and $n_0$ = 60 households (black dots) in Majete.", fig.asp=1----
set.seed(1234)
init.sample <-
  discrete.inhibit.sample(obj = majete, size = 60, delta = 0.4,
                          k = 0,  delta.fix = FALSE,
                          poly = border, plotit = TRUE)

## ----applicestim1, include=TRUE, eval=TRUE, echo=TRUE, fig.keep='none'------
mrdt <- init.sample[[4]]
glmfit <- glm(rdt~1, data = mrdt,
              family = binomial(link = logit))
ID.coords <- create.ID.coords(data=as.data.frame(mrdt),
                              coords=~st_coordinates(mrdt))
mrdt$units.m <- rep(1,nrow(mrdt))

par0 <- c(coef(glmfit),cov.pars=c(0.93171,3.9549))
control.mcmc <- control.mcmc.MCML(n.sim = 5500,burnin=500,thin=5)
model.fit <-
  binomial.logistic.MCML(rdt~1, units.m=~units.m,par0=par0,
                         coords=~st_coordinates(mrdt),data=mrdt,
                         ID.coords = ID.coords, kappa=0.5,
                         control.mcmc=control.mcmc, method="BFGS",
                         fixed.rel.nugget = 0, start.cov.pars=c(par0[3]),
                         messages = FALSE, plot.correlogram = FALSE)

summary(model.fit, log.cov.pars = FALSE)

## ----applicadapt1, include=TRUE, eval=TRUE, echo=TRUE, fig.height=10, fig.width=10, fig.cap = "\\label{fig:applicadapt1} Adaptive sampling design with $\\delta$ = 150 meters and $b = 40$, Blue dots ($n_0$ = 60) are the initial sampling households. Red dots ($n_a$ = 40) are adaptive samples added after analysing data from the initial design.", fig.asp=1----
avail.locs <- majete[!(majete$geometry) %in% (mrdt$geometry),]
model.pred <-
  spatial.pred.binomial.MCML(model.fit,
                             grid.pred=unique(st_coordinates(avail.locs)),
                             control.mcmc=control.mcmc, type = "marginal",
                             scale.predictions = "prevalence",
                             standard.errors = TRUE, thresholds = 0.15,
                             scale.thresholds = "prevalence",
                             messages = FALSE, plot.correlogram = FALSE)

pred.vars <- as.data.frame(cbind(model.pred$grid,
                                c(model.pred$prevalence$standard.errors)^2))

colnames(pred.vars)<- c("coord1", "coord2", "pred.var")
pred.vars <- sf::st_as_sf(pred.vars, coords = c('coord1', 'coord2'))
st_crs(pred.vars) <- st_crs(mrdt)
adapt.sample.pv <-
   adaptive.sample(obj1 = pred.vars, obj2 = mrdt,
                  pred.var.col = 1, criterion = "predvar",
                  delta = 0.15, batch.size = 40,
                  poly = border, plotit = TRUE)

## ----applicestim2, include=TRUE, eval=TRUE, echo=TRUE, fig.keep='none'------
mrdt <- majete[(majete$geometry) %in%
                 (adapt.sample.pv$sample.locs$curr.sample$geometry),]

ID.coords <- create.ID.coords(data=as.data.frame(mrdt),
                              coords=~st_coordinates(mrdt))
mrdt$units.m <- rep(1,nrow(mrdt))

par0 <- c(coef(model.fit))
model.fit <-
  binomial.logistic.MCML(rdt~1, units.m=~units.m,par0=par0,
                         coords=~st_coordinates(mrdt),data=mrdt,
                         ID.coords = ID.coords,
                         control.mcmc=control.mcmc, kappa=0.5,
                         fixed.rel.nugget = 0,
                         start.cov.pars=c(par0[3]),
                         method="BFGS", messages = FALSE,
                         plot.correlogram = FALSE)

summary(model.fit, log.cov.pars = FALSE)

## ----applicpred2, include=TRUE, eval=TRUE, echo=TRUE, fig.height=9, fig.width=9, fig.cap = "\\label{fig:applicpred2} (a) Malaria prevalence in Majete. (b) Exceedance probabilities $P(x; 0.15)$ for the predictions. $P(x; 0.15)$ = $P$($prev > 0.15$ at location $x$). (c) Standard errors of predictions.", fig.asp=1----
library(splancs)
pred.poly <- as_Spatial(border)@polygons[[1]]@Polygons[[1]]@coords
grid.pred <- gridpts(pred.poly, xs=0.05, ys=0.05)

model.pred <-
  spatial.pred.binomial.MCML(model.fit, grid.pred=grid.pred,
                             control.mcmc=control.mcmc,
                             type = "marginal",
                             scale.predictions = "prevalence",
                             standard.errors = TRUE, thresholds = 0.15,
                             scale.thresholds = "prevalence",
                             messages = FALSE, plot.correlogram = FALSE)
##1. Prevalence predictions
prevpred <-
  rasterFromXYZ(cbind(model.pred$grid[,1],
                      model.pred$grid[,2],
                      model.pred$prevalence$predictions))
prevpred <- raster::disaggregate(prevpred, fact = 10,
                                 method = "bilinear")

##2. Std error
stderror <-
  rasterFromXYZ(cbind(model.pred$grid[,1],
                      model.pred$grid[,2],
                      model.pred$prevalence$standard.errors))
stderror <- raster::disaggregate(stderror, fact = 10,
                                 method = "bilinear")

##3. Exceedance probablities
exceed <-
  rasterFromXYZ(cbind(model.pred$grid[,1],
                      model.pred$grid[,2],
                      model.pred$exceedance.prob))
exceed <- raster::disaggregate(exceed, fact = 10,
                                 method = "bilinear")
par(mfrow = c(2,2))
plot(prevpred, main = "(a)", col = viridis(256, direction = -1))
plot(exceed, main="(b)", zlim = c(0,1), col = viridis(256, direction = -1))
plot(stderror, main = "(c)", col = viridis(256, direction = -1))
par(mfrow = c(1,1))

