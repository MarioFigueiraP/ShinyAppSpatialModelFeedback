ss_userAction.Log <- periscope:::fw_get_user_log()

# Define server logic required to draw a histogram
shinyServer(function(input, output, session) {
    observeEvent(input$parametersMapSim|input$parametersSampleSim, {
        updateBox("resultbox",
                  action = "update",
                  options = list(width=ifelse(input$parametersMapSim|input$parametersSampleSim,9,12)))
    })
    
    observeEvent(input$LoadingPanel, {
        updateBox("TablePanelRaster",
                  action = "update",
                  options = list(width=ifelse(input$LoadingPanel,8,12)))
    })
    
    observeEvent(input$FitRandomPanel, {
        updateBox("plotstablesFitRandom",
                  action = "update",
                  options = list(width=ifelse(input$FitRandomPanel,9,12)))
    })
    
    # Server_Auxiliar_Functions_Section ----
    
    RW1 <- function(N, x0, mu=0, stdev){
        z<-cumsum(rnorm(n=N, mean=0, sd=stdev))
        t<-1:N
        x<-x0+t*mu+z
        return(x)
    }

    RW2 <- function(N, x0, mu=0, stdev){
        x <- vector(mode="numeric", length=N)
        x[1] <- x0+rnorm(1,mean=mu, sd=stdev)
        x[2] <- 2*x[1]-x0+rnorm(1,mean=mu, sd=stdev)
        for(i in 3:N){
            x[i] <- 2*x[i-1]+rnorm(1,mean=mu,sd=stdev)-x[i-2]
        }
        return(x)
    }
    
    lin.interpolate <- function(x,xseq,yseq){
        a <- (xseq-x)**2
        seq1 <- which.min(a); seq2 <- which.min(a[-seq1])
        seq2 <- which.min(a[-seq1])+abs(seq1-seq2-1)
        y <- (x-xseq[seq2])/(xseq[seq2]-xseq[seq1])*(yseq[seq2]-yseq[seq1])+yseq[seq2]
        return(y)
    }
    
    ProjectionFunctionRegularGrid <- function(lim1=c(0,1), lim2=c(0,1), loc, z, proj.grid){
        lattice <- fm_lattice_2d(seq(lim1[1],lim1[2],length.out=2),seq(lim2[1],lim2[2],length.out=2))
        mesh <- inla.mesh.create(loc=loc, boundary=lattice$segm, refine=list(max.edge=0.08))
        df <- data.frame(x=loc[,1], y=loc[,2], z=z)
        indx <- prodlim::row.match(as.data.frame(cbind(df$x,df$y)), mesh$loc[,1:2]) # library(prodlim)
        # df[order(indx),1:2]==mesh$loc
        df <- df[order(indx),]
        A <- inla.spde.make.A(mesh=mesh,loc=proj.grid)
        zproj <- A%*%df$z
        return(zproj)
    }
    
    InterpolateIrrGrid <- function(z, loc, gridInter){
        mesh <- inla.mesh.create(loc=loc)
        match.ind <- prodlim::row.match(as.data.frame(loc), as.data.frame(mesh$loc[,1:2]))
        z.ext <- matrix(data=NA, nrow=nrow(mesh$loc))
        z.ext[match.ind] <- z
        A.inter <- inla.spde.make.A(mesh=mesh,loc=as.matrix(gridInter))
        z.inter <- as.vector(A.inter%*%z.ext)
        InterpolateGgplot <- ggplot(data=data.frame(x=gridInter[,1], y=gridInter[,2],z=z.inter)) +
            geom_tile(aes(x=x,y=y,fill=z)) + scale_fill_viridis_c(option = "turbo")
        result <- list(Plot=InterpolateGgplot, DataInter=data.frame(x=gridInter[,1],y=gridInter[,2],z=z.inter))
        return(result)
    }
    
    RefineInterMesh <- function(z, loc, dimmap=200){
        mesh <- inla.mesh.create(loc=loc)
        match.ind <- prodlim::row.match(as.data.frame(loc), as.data.frame(mesh$loc[,1:2]))
        z.ext <- matrix(data=NA, nrow=nrow(mesh$loc))
        z.ext[match.ind] <- z
        grid.pred <- as.matrix(expand.grid(x=seq(min(mesh$loc[,1]),max(mesh$loc[,1]),length.out=dimmap), 
                                           y=seq(min(mesh$loc[,2]),max(mesh$loc[,2]),length.out=dimmap)))
        poly.mesh <- Polygon(unique(mesh$loc[mesh$segm$bnd$idx,1:2]))
        polys.mesh <- Polygons(list(poly.mesh),"poly.mesh")
        sppoly.mesh <- SpatialPolygons(list(polys.mesh), proj4string=CRS("+proj=longlat"))
        indx <- !is.na(over(SpatialPoints(grid.pred, sppoly.mesh@proj4string), sppoly.mesh))
        gridInter <- as.matrix(grid.pred)[indx,]
        A.inter <- inla.spde.make.A(mesh=mesh,loc=gridInter)
        z.inter <- as.vector(A.inter%*%z.ext)
        InterpolateGgplot <- ggplot(data=data.frame(x=gridInter[,1], y=gridInter[,2], z=z.inter)) +
            geom_tile(aes(x=x,y=y,fill=z)) + scale_fill_viridis_c(option = "turbo")
        result <- list(Plot=InterpolateGgplot, DataInter=data.frame(x=gridInter[,1],y=gridInter[,2], z=z.inter))
        return(result)
    }
    
    logdunif <- function(x,lim){
      z <-  -log(2) + log(x)
      if(x>lim[1]&x<lim[2]){
        log(1/diff(lim)) + z
      } else{-5000}
    }
    
    mesh.dual <- function(mesh){
        if (mesh$manifold=='R2'){
            ce <- t(sapply(1:nrow(mesh$graph$tv), function(i)
                colMeans(mesh$loc[mesh$graph$tv[i, ], 1:2])))
            library(parallel)
            pls <- mclapply(1:mesh$n, function(i) {
                p <- unique(Reduce('rbind', lapply(1:3, function(k) {
                    j <- which(mesh$graph$tv[,k]==i)
                    if (length(j)>0)
                        return(rbind(ce[j, , drop=FALSE],
                                     cbind(mesh$loc[mesh$graph$tv[j, k], 1] +
                                               mesh$loc[mesh$graph$tv[j, c(2:3,1)[k]], 1],
                                           mesh$loc[mesh$graph$tv[j, k], 2] +
                                               mesh$loc[mesh$graph$tv[j, c(2:3,1)[k]], 2])/2))
                    else return(ce[j, , drop=FALSE])
                })))
                j1 <- which(mesh$segm$bnd$idx[,1]==i)
                j2 <- which(mesh$segm$bnd$idx[,2]==i)
                if ((length(j1)>0) | (length(j2)>0)) {
                    p <- unique(rbind(mesh$loc[i, 1:2], p,
                                      mesh$loc[mesh$segm$bnd$idx[j1, 1], 1:2]/2 +
                                          mesh$loc[mesh$segm$bnd$idx[j1, 2], 1:2]/2,
                                      mesh$loc[mesh$segm$bnd$idx[j2, 1], 1:2]/2 +
                                          mesh$loc[mesh$segm$bnd$idx[j2, 2], 1:2]/2))
                    yy <- p[,2]-mean(p[,2])/2-mesh$loc[i, 2]/2
                    xx <- p[,1]-mean(p[,1])/2-mesh$loc[i, 1]/2
                }
                else {
                    yy <- p[,2]-mesh$loc[i, 2]
                    xx <- p[,1]-mesh$loc[i, 1]
                }
                Polygon(p[order(atan2(yy,xx)), ])
            })
            return(SpatialPolygons(lapply(1:mesh$n, function(i)
                Polygons(list(pls[[i]]), i))))
        }
        else stop("It only works for R2!")
    }
    
    # Server_Presentation_Section ----
    
    src = "https://raw.githubusercontent.com/MarioFigueiraP/ShinyAppSpatialModelFeedback/main/www/Simulation_Response_Sampling.png"
    output$Presentation1 <- renderText({c('<img src="',src,'">')})
    
    # Server_Simulation_Section ----
    
    meshSim <- eventReactive(input$makeSim,{
        xlattice0 <- seq(input$limlattice[1], input$limlattice[2], length.out=input$lengthlattice) 
        ylattice0 <- seq(input$limlattice[1], input$limlattice[2], length.out=input$lengthlattice)
        lattice0 <-  fm_lattice_2d(xlattice0, ylattice0)
        mesh <-  inla.mesh.create(lattice=lattice0, refine=list(max.edge=0.08*abs(diff(input$limlattice))))
        result <- list(mesh=mesh, lattice0=lattice0)
        return(result)
    })

    dataggplotMeshSim <- function() {
        data.frame(Longitude=meshSim()$mesh$loc[,1],Latitude=meshSim()$mesh$loc[,2])
    }

    ggplotMeshSim <- function(){
        ggplot() + gg(meshSim()$mesh)+theme_bw() + xlab("Latitude") + ylab("Longitude") +
            ggtitle("Mesh over the study region") +
            theme(plot.title=element_text(color = "black", size = 16,
                                          face = "bold", hjust = 0.5))
    }

    downloadFile(
        id = "ggplotMeshSim",
        logger = ss_userAction.Log,
        filenameroot = "ggplotMeshSim",
        aspectratio  = 1,
        downloadfxns = list(png  = ggplotMeshSim,
                            txt  = dataggplotMeshSim,
                            csv  = dataggplotMeshSim)
    )

    downloadablePlot(
        id = "ggplotMeshSim",
        logger = ss_userAction.Log,
        filenameroot = "ggplotMeshSim",
        aspectratio  = 1,
        downloadfxns = list(png  = ggplotMeshSim,
                            txt  = dataggplotMeshSim,
                            csv  = dataggplotMeshSim),
        visibleplot  = ggplotMeshSim)

    SimMap <- eventReactive(input$makeSim, {
        if(!is.na(input$seedGlobal)){set.seed(input$seedGlobal, kind="Mersenne-Twister", normal.kind="Inversion")}
        lattice0 <- meshSim()$lattice0
        mesh <- meshSim()$mesh
        lattice <- list(loc=mesh$loc, x=lattice0$x, y=lattice0$y)
        sigma0 <- input$sigma.matern; range0 <- input$range.matern
        kappa0 <- sqrt(8)/range0; tau0 <- 1/(sqrt(4*pi)*kappa0*sigma0)
        spde <- inla.spde2.matern(mesh, B.tau=cbind(log(tau0),-1,+1), B.kappa=cbind(log(kappa0),0,-1),
                                  theta.prior.mean=c(0,0), theta.prior.prec=c(0.1,1))
        Qu <-  inla.spde.precision(spde, theta=c(0, 0))

        if(is.na(input$seedSP)){u <- inla.qsample(n=1, Q=Qu)} else{
            u <- inla.qsample(n=1, Q=Qu, seed = input$seedSP)}
        u <-  as.vector(u[,1])

        #Now we must build the bathymetry and its effect on the predictor
        bat.fun <- function(x,y,s){
            z <- eval(parse(text=s))
        }
        bathymetry <- bat.fun(x=lattice$loc[,1], y=lattice$loc[,2], s=input$bathymetry.formula)
        beta <- as.numeric(unlist(strsplit(input$beta, ",")))

        if(input$typeBathyeff=="lin"){
            bat.eff <- beta[2]*bathymetry
        } else if(input$typeBathyeff=="rw1"){
            min.eff0 <- as.numeric(unlist(strsplit(input$minmax.effect.rw1, split=",")))[1]
            max.eff0 <- as.numeric(unlist(strsplit(input$minmax.effect.rw1, split=",")))[2]
            min.eff <- min.eff0-1; max.eff <- max.eff0+1
            tprewhile <- Sys.time();k <- 1
            while(min.eff0>min.eff|max.eff0<max.eff){ #be aweare about the relation between precision and number of knots 
                bat.eff.seq <- RW1(N=input$nknots.rw1, x0=input$init.rw1, mu=0, stdev=1/(input$prec.rw1)**0.5)
                min.eff <- min(bat.eff.seq);max.eff <- max(bat.eff.seq)
                tinwhile <- Sys.time()
                if (as.numeric(difftime(tinwhile,tprewhile, units="secs"))>20*k){
                    k <- k+1
                    warningMessage <- paste0("Simulate the random walk (rw1) is taking a while (more than ", 
                                   round(as.numeric(difftime(tinwhile,tprewhile, units="secs")), digits=0),
                           " secs). Maybe would be better restart the process with a wider interval between 
                           extreme values or with a higher precision.")
                    showNotification(ui=warningMessage, duration=10, closeButton=TRUE, type="warning")
                }
            }
            bat.seq <- seq(min(bathymetry), max(bathymetry), length.out=input$nknots.rw1)
            bat.eff <- sapply(X=bathymetry, FUN=lin.interpolate, xseq=bat.seq, yseq=bat.eff.seq)
        } else if(input$typeBathyeff=="rw2"){
            min.eff0 <- as.numeric(unlist(strsplit(input$minmax.effect.rw2, split=",")))[1]
            max.eff0 <- as.numeric(unlist(strsplit(input$minmax.effect.rw2, split=",")))[2]
            min.eff <- min.eff0-1; max.eff <- max.eff0+1
            tprewhile <- Sys.time();k <- 1
            while(min.eff0>min.eff|max.eff0<max.eff){ #be aweare about the relation between precision and number of knots 
                bat.eff.seq <- RW2(N=input$nknots.rw2, x0=input$init.rw2, mu=0, stdev=1/(input$prec.rw2)**0.5)
                min.eff <- min(bat.eff.seq);max.eff <- max(bat.eff.seq)
                tinwhile <- Sys.time()
                if (as.numeric(difftime(tinwhile,tprewhile, units="secs"))>20*k){
                    k <- k+1
                    warningMessage <- paste0("Simulate the random walk (rw1) is taking a while (more than ", 
                                             round(as.numeric(difftime(tinwhile,tprewhile, units="secs")), digits=0),
                                             " secs). Maybe would be better restart the process with a wider interval between 
                           extreme values or with a higher precision.")
                    showNotification(ui=warningMessage, duration=10, closeButton=TRUE, type="warning")
                }
            }
            bat.seq <- seq(min(bathymetry), max(bathymetry), length.out=input$nknots.rw2)
            bat.eff <- sapply(X=bathymetry, FUN=lin.interpolate, xseq=bat.seq, yseq=bat.eff.seq)
        } else if(input$typeBathyeff=="customfunction"){
            bat.effect.function <- function(b,s){
                bat.eff <- eval(parse(text=s))
            }
            bat.eff <- bat.effect.function(b=bathymetry, s=input$formula.eff.bathymetry)
        }

        pred <- beta[1] + bat.eff + u

        if(input$datadistributionSim=="gaussian"){
            sigma.epsilon=sqrt(input$var)
            abun <- rnorm(length(pred), (pred), sigma.epsilon)
        } else if(input$datadistributionSim=="gamma"){
            phi.gamma=input$var
            abun <- rgamma(length(pred), exp(pred)**2/phi.gamma, exp(pred)/phi.gamma)
        } else if(input$datadistributionSim=="bernoulli"){
            pres_abs <- rbinom(n=length(pred),size=1, prob=exp(pred)/(1+exp(pred))) 
        } else {print("Distribution must be \"gaussian\" or \"gamma\".")}
        
        grid <- as.matrix(expand.grid(seq(input$limlattice[1], input$limlattice[2],
                                          length.out=input$dimmap),
                                      seq(input$limlattice[1], input$limlattice[2],
                                          length.out=input$dimmap)))
        A <- inla.spde.make.A(mesh=mesh, loc=grid)

        spatial.effect <- as.vector(A %*% u)
        abundance.field <- as.vector(A %*% abun)
        bathymetry.field <- as.vector(A %*% bathymetry)
        bathymetry.effect.field <- as.vector(A %*% bat.eff)

        BathymetryEffect <- list(bathymetry=bathymetry, bathymetry.effect=bat.eff)
        
        DataSim <- list(gridx = grid[,1], gridy = grid[,2], abundance.field = abundance.field,
                        bathymetry.field = bathymetry.field, bathymetry.effect.field = bathymetry.effect.field,
                        spatial.effect = spatial.effect)
        result <- list(DataSim=list(DataSim=DataSim), BathymetryEffect=list(BathymetryEffect=BathymetryEffect))
        return(result)
    })
    
    ## Server_Simulation_Results_Subsection ====
    
    dataggplotSpSim <- function(){
        data.frame(Latitude=SimMap()$DataSim$DataSim$gridx, Longitude=SimMap()$DataSim$DataSim$gridy,
                   Spatial.effect=SimMap()$DataSim$DataSim$spatial.effect)
    }

    ggplotSpSim <- function(){
        DFM <- data.frame(Latitude=SimMap()$DataSim$DataSim$gridx, Longitude=SimMap()$DataSim$DataSim$gridy,
                          Spatial.effect=SimMap()$DataSim$DataSim$spatial.effect)
        ggplot <- ggplot(DFM) + geom_tile(aes(Latitude, Longitude, fill = Spatial.effect)) + 
            scale_fill_viridis_c(option = "turbo") + theme_bw() +
            xlab("Latitude") + ylab("Longitude") + ggtitle("Spatial effect") + labs(fill="Values") +
            theme(plot.title=element_text(color = "black", size = 16, face = "bold", hjust = 0.5))
        return(ggplot)
    }

    downloadFile(
        id = "ggplotSpSim",
        logger = ss_userAction.Log,
        filenameroot = "ggplotSpSim",
        aspectratio  = 1,
        downloadfxns = list(png  = ggplotSpSim,
                            txt  = dataggplotSpSim,
                            csv  = dataggplotSpSim)
    )

    downloadablePlot(
        id = "ggplotSpSim",
        logger = ss_userAction.Log,
        filenameroot = "ggplotSpSim",
        aspectratio  = 1,
        downloadfxns = list(png  = ggplotSpSim,
                            txt  = dataggplotSpSim,
                            csv  = dataggplotSpSim),
        visibleplot  = ggplotSpSim)

    # Data and bathymetric chart
    
    dataggplotBatSim <- function(){
        data.frame(Latitude=SimMap()$DataSim$DataSim$gridx, Longitude=SimMap()$DataSim$DataSim$gridy,
                   Bathymetry.field=SimMap()$DataSim$DataSim$bathymetry.field)
    }
    
    ggplotBatSim <- function(){
        DFM <- data.frame(Latitude=SimMap()$DataSim$DataSim$gridx, Longitude=SimMap()$DataSim$DataSim$gridy,
                          Bathymetry.field=SimMap()$DataSim$DataSim$bathymetry.field)
        ggplot <- ggplot(DFM) + geom_tile(aes(Latitude, Longitude, fill = Bathymetry.field)) + 
            scale_fill_viridis_c(option = "turbo") + theme_bw() + labs(fill="Values") +
            xlab("Latitude") + ylab("Longitude") + ggtitle("Covariate map") +
            theme(plot.title=element_text(color = "black", size = 16, face = "bold", hjust = 0.5))
        return(ggplot)
    }
    
    downloadFile(
        id = "ggplotBatChartSim",
        logger = ss_userAction.Log,
        filenameroot = "ggplotBatChartSim",
        aspectratio  = 1,
        downloadfxns = list(png  = ggplotBatSim,
                            txt  = dataggplotBatSim,
                            csv  = dataggplotBatSim)
    )
    
    downloadablePlot(
        id = "ggplotBatChartSim",
        logger = ss_userAction.Log,
        filenameroot = "ggplotBatChartSim",
        aspectratio  = 1,
        downloadfxns = list(png  = ggplotBatSim,
                            txt  = dataggplotBatSim,
                            csv  = dataggplotBatSim),
        visibleplot  = ggplotBatSim)
    
    # Data and relation between the covariate and its effect
    
    dataggplotBateffSim <- function(){
        data.frame(Bathymetry=SimMap()$BathymetryEffect$BathymetryEffect$bathymetry, 
                   Bathimetric.effect=SimMap()$BathymetryEffect$BathymetryEffect$bathymetry.effect)
    }
    
    ggplotBateffSim <- function(){
        DFM <- data.frame(Bathymetry=SimMap()$BathymetryEffect$BathymetryEffect$bathymetry, 
                          Bathimetric.effect=SimMap()$BathymetryEffect$BathymetryEffect$bathymetry.effect)
        ggplot <- ggplot(DFM, aes(x=Bathymetry, y=Bathimetric.effect)) + geom_point() + theme_bw() +
            xlab("Bathymetry") + ylab("Covariate effect") + ggtitle("Covariate effect") + labs(fill="Values") +
            theme(plot.title=element_text(color = "black", size = 16, face = "bold", hjust = 0.5))
        return(ggplot)
    }
    
    downloadFile(
        id = "ggplotBatEffSim",
        logger = ss_userAction.Log,
        filenameroot = "ggplotBatEfftSim",
        aspectratio  = 1,
        downloadfxns = list(png  = ggplotBateffSim,
                            txt  = dataggplotBateffSim,
                            csv  = dataggplotBateffSim)
    )
    
    downloadablePlot(
        id = "ggplotBatEffSim",
        logger = ss_userAction.Log,
        filenameroot = "ggplotBatEffSim",
        aspectratio  = 1,
        downloadfxns = list(png  = ggplotBateffSim,
                            txt  = dataggplotBateffSim,
                            csv  = dataggplotBateffSim),
        visibleplot  = ggplotBateffSim)
    
    # Covariate effect chart
    
    dataggplotBatEffChartSim <- function(){
        data.frame(Latitude=SimMap()$DataSim$DataSim$gridx, Longitude=SimMap()$DataSim$DataSim$gridy,
                   Bathymetry.field=SimMap()$DataSim$DataSim$bathymetry.effect.field)
    }
    
    ggplotBatEffChartSim <- function(){
        DFM <- data.frame(Latitude=SimMap()$DataSim$DataSim$gridx, Longitude=SimMap()$DataSim$DataSim$gridy,
                          Bathymetry.effect.field=SimMap()$DataSim$DataSim$bathymetry.effect.field)
        ggplot <- ggplot(DFM) + geom_tile(aes(Latitude, Longitude, fill = Bathymetry.effect.field)) + 
            scale_fill_viridis_c(option = "turbo") + theme_bw() +
            xlab("Latitude") + ylab("Longitude") + ggtitle("Covariate effect map") + labs(fill="Values") +
            theme(plot.title=element_text(color = "black", size = 16, face = "bold", hjust = 0.5))
        return(ggplot)
    }
    
    downloadFile(
        id = "ggplotBatEffChartSim",
        logger = ss_userAction.Log,
        filenameroot = "ggplotBatEffChartSim",
        aspectratio  = 1,
        downloadfxns = list(png  = ggplotBatEffChartSim,
                            txt  = dataggplotBatEffChartSim,
                            csv  = dataggplotBatEffChartSim)
    )
    
    downloadablePlot(
        id = "ggplotBatEffChartSim",
        logger = ss_userAction.Log,
        filenameroot = "ggplotBatEffChartSim",
        aspectratio  = 1,
        downloadfxns = list(png  = ggplotBatEffChartSim,
                            txt  = dataggplotBatEffChartSim,
                            csv  = dataggplotBatEffChartSim),
        visibleplot  = ggplotBatEffChartSim)
    
    # Variable response map
    
    dataggplotAbundanceSim <- function(){
        data.frame(Latitude=SimMap()$DataSim$DataSim$gridx, Longitude=SimMap()$DataSim$DataSim$gridy,
                   Abundance=SimMap()$DataSim$DataSim$abundance.field)
    }
    
    ggplotAbundanceSim <- function(){
        DFM <- data.frame(Latitude=SimMap()$DataSim$DataSim$gridx, Longitude=SimMap()$DataSim$DataSim$gridy,
                          Abundance.field=SimMap()$DataSim$DataSim$abundance.field)
        ggplot <- ggplot(DFM) + geom_tile(aes(Latitude, Longitude, fill = Abundance.field)) + 
            scale_fill_viridis_c(option = "turbo") + theme_bw() +
            xlab("Latitude") + ylab("Longitude") + ggtitle("Response variable map") + labs(fill="Values") +
            theme(plot.title=element_text(color = "black", size = 16, face = "bold", hjust = 0.5))
        return(ggplot)
    }
    
    downloadFile(
        id = "ggplotAbundanceSim",
        logger = ss_userAction.Log,
        filenameroot = "ggplotAbundanceSim",
        aspectratio  = 1,
        downloadfxns = list(png  = ggplotAbundanceSim,
                            txt  = dataggplotAbundanceSim,
                            csv  = dataggplotAbundanceSim)
    )
    
    downloadablePlot(
        id = "ggplotAbundanceSim",
        logger = ss_userAction.Log,
        filenameroot = "ggplotAbundanceSim",
        aspectratio  = 1,
        downloadfxns = list(png  = ggplotAbundanceSim,
                            txt  = dataggplotAbundanceSim,
                            csv  = dataggplotAbundanceSim),
        visibleplot  = ggplotAbundanceSim)
    
    ### Simulation of the samples
    
    Ind.sampling <- eventReactive(input$makeSample, {
        if(!is.na(input$seedSampleR)){set.seed(input$seedSampleR, kind="Mersenne-Twister", normal.kind="Inversion")}
        indx <- sample(1:length(SimMap()$DataSim$DataSim$gridx), size=input$niid.samples)
        RandomSample <- list(Latitude=SimMap()$DataSim$DataSim$gridx[indx], 
                                   Longitude=SimMap()$DataSim$DataSim$gridy[indx], 
                                   Abundance=SimMap()$DataSim$DataSim$abundance.field[indx],
                                   Bathymetry=SimMap()$DataSim$DataSim$bathymetry.field[indx])
        return(RandomSample)
    })
    
    dataggplotRandomSampleSim <- function(){
        data.frame(Latitude=Ind.sampling()$Latitude, Longitude=Ind.sampling()$Longitude,
                   Bathymetry=Ind.sampling()$Bathymetry, Abundance=Ind.sampling()$Abundance)
    }
    
    ggplotRandomSampleSim <- function(){
        DFM <- data.frame(Latitude=SimMap()$DataSim$DataSim$gridx, Longitude=SimMap()$DataSim$DataSim$gridy,
                          Abundance.field=SimMap()$DataSim$DataSim$abundance.field)
        DFS <- data.frame(Latitude=Ind.sampling()$Latitude, Longitude=Ind.sampling()$Longitude)
        ggplot <- ggplot(DFM) + geom_tile(aes(Latitude, Longitude, fill = Abundance.field)) + 
            scale_fill_viridis_c(option = "turbo") + theme_bw() +
            xlab("Latitude") + ylab("Longitude") + ggtitle("Independent and random sampling process") + 
            geom_point(data=DFS, aes(x=Latitude,y=Longitude)) + labs(fill="Values") +
            theme(plot.title=element_text(color = "black", size = 16, face = "bold", hjust = 0.5))
        return(ggplot)
    }
    
    downloadFile(
        id = "ggplotRandomSampleSim",
        logger = ss_userAction.Log,
        filenameroot = "ggplotRandomSampleSim",
        aspectratio  = 1,
        downloadfxns = list(png  = ggplotRandomSampleSim,
                            txt  = dataggplotRandomSampleSim,
                            csv  = dataggplotRandomSampleSim)
    )
    
    downloadablePlot(
        id = "ggplotRandomSampleSim",
        logger = ss_userAction.Log,
        filenameroot = "ggplotRandomSampleSim",
        aspectratio  = 1,
        downloadfxns = list(png  = ggplotRandomSampleSim,
                            txt  = dataggplotRandomSampleSim,
                            csv  = dataggplotRandomSampleSim),
        visibleplot  = ggplotRandomSampleSim)
    
    
    Pref.sampling <- eventReactive(input$makeSample, {
        if(!is.na(input$seedSampleP)){set.seed(input$seedSampleP, kind="Mersenne-Twister", normal.kind="Inversion")}
        indx <- sample(1:length(SimMap()$DataSim$DataSim$gridx), size=input$nps.samples,
                       prob=exp(input$r.scale*SimMap()$DataSim$DataSim$abundance.field/max(SimMap()$DataSim$DataSim$abundance.field)))
        PreferentialSample <- list(Latitude=SimMap()$DataSim$DataSim$gridx[indx],
                             Longitude=SimMap()$DataSim$DataSim$gridy[indx],
                             Abundance=SimMap()$DataSim$DataSim$abundance.field[indx],
                             Bathymetry=SimMap()$DataSim$DataSim$bathymetry.field[indx])
        return(PreferentialSample)
    })
    
    dataggplotPrefSampleSim <- function(){
        data.frame(Latitude=Pref.sampling()$Latitude, Longitude=Pref.sampling()$Longitude,
                   Bathymetry=Pref.sampling()$Bathymetry, Abundance=Pref.sampling()$Abundance)
    }
    
    ggplotPrefSampleSim <- function(){
        DFM <- data.frame(Latitude=SimMap()$DataSim$DataSim$gridx, Longitude=SimMap()$DataSim$DataSim$gridy,
                          Abundance.field=SimMap()$DataSim$DataSim$abundance.field)
        DFS <- data.frame(Latitude=Pref.sampling()$Latitude, Longitude=Pref.sampling()$Longitude)
        ggplot <- ggplot(DFM) + geom_tile(aes(Latitude, Longitude, fill = Abundance.field)) + 
            scale_fill_viridis_c(option = "turbo") + theme_bw() +
            xlab("Latitude") + ylab("Longitude") + ggtitle("Preferential sampling process") + 
            geom_point(data=DFS, aes(x=Latitude,y=Longitude)) + labs(fill="Values") +
            theme(plot.title=element_text(color = "black", size = 16, face = "bold", hjust = 0.5))
        return(ggplot)
    }
    
    downloadFile(
        id = "ggplotPrefSampleSim",
        logger = ss_userAction.Log,
        filenameroot = "ggplotPrefSampleSim",
        aspectratio  = 1,
        downloadfxns = list(png  = ggplotPrefSampleSim,
                            txt  = dataggplotPrefSampleSim,
                            csv  = dataggplotPrefSampleSim)
    )
    
    downloadablePlot(
        id = "ggplotPrefSampleSim",
        logger = ss_userAction.Log,
        filenameroot = "ggplotPrefSampleSim",
        aspectratio  = 1,
        downloadfxns = list(png  = ggplotPrefSampleSim,
                            txt  = dataggplotPrefSampleSim,
                            csv  = dataggplotPrefSampleSim),
        visibleplot  = ggplotPrefSampleSim)
    
    Lgcp.sampling <- eventReactive(input$makeSample, {
      if(!is.na(input$seedSampleP)){set.seed(input$seedSampleP, kind="Mersenne-Twister", normal.kind="Inversion")}
      indx <- sample(1:length(SimMap()$DataSim$DataSim$gridx), size=input$nps.samples,
                     prob=exp(input$r.scale*SimMap()$DataSim$DataSim$abundance.field/max(SimMap()$DataSim$DataSim$abundance.field)))
      PreferentialSample <- list(Latitude=SimMap()$DataSim$DataSim$gridx[indx],
                                 Longitude=SimMap()$DataSim$DataSim$gridy[indx],
                                 Precense=rep(1, times=length(SimMap()$DataSim$DataSim$abundance.field[indx])),
                                 Bathymetry=SimMap()$DataSim$DataSim$bathymetry.field[indx])
      return(PreferentialSample)
    })
    
    dataggplotLgcpSampleSim <- function(){
      data.frame(Latitude=Lgcp.sampling()$Latitude, Longitude=Lgcp.sampling()$Longitude,
                 Bathymetry=Lgcp.sampling()$Bathymetry, Abundance=Lgcp.sampling()$Abundance)
    }

    Mixture.sampling <- eventReactive(input$makeSampleMixture,{
      if(!is.na(input$seedSample.mixture1)){set.seed(input$seedSample.mixture1, kind="Mersenne-Twister", normal.kind="Inversion")}
      indx <- sample(1:length(SimMap()$DataSim$DataSim$gridx), size=input$n.mixture1.samples,
                     prob=exp(input$r.scale.mixture1*SimMap()$DataSim$DataSim$abundance.field/max(SimMap()$DataSim$DataSim$abundance.field)))
      MixtureSample1 <- list(Latitude=SimMap()$DataSim$DataSim$gridx[indx],
                             Longitude=SimMap()$DataSim$DataSim$gridy[indx],
                             Bathymetry=SimMap()$DataSim$DataSim$bathymetry.field[indx],
                             Abundance=SimMap()$DataSim$DataSim$abundance.field[indx],
                             group=rep("Mix_1",length(indx)))
      if(!is.na(input$seedSample.mixture2)){set.seed(input$seedSample.mixture2, kind="Mersenne-Twister", normal.kind="Inversion")}
      indx2 <- sample(1:length(SimMap()$DataSim$DataSim$gridx), size=input$n.mixture2.samples,
                     prob=exp(input$r.scale.mixture2*SimMap()$DataSim$DataSim$abundance.field/max(SimMap()$DataSim$DataSim$abundance.field)))
      MixtureSample2 <- list(Latitude=SimMap()$DataSim$DataSim$gridx[indx2],
                             Longitude=SimMap()$DataSim$DataSim$gridy[indx2],
                             Bathymetry=SimMap()$DataSim$DataSim$bathymetry.field[indx2],
                             Abundance=SimMap()$DataSim$DataSim$abundance.field[indx2],
                             group=rep("Mix_2",length(indx2)))
      MixtureSample <- data.frame(Latitude=c(MixtureSample1$Latitude, MixtureSample2$Latitude), 
                                  Longitude=c(MixtureSample1$Longitude, MixtureSample2$Longitude),
                                  Abundance=c(MixtureSample1$Abundance, MixtureSample2$Abundance),
                                  Bathymetry=c(MixtureSample1$Bathymetry, MixtureSample2$Bathymetry),
                                  Mixture=c(MixtureSample1$group, MixtureSample2$group))
      return(list(MixtureSample=MixtureSample))
    })
    
    dataggplotMixtureSampleSim <- function(){
      Mixture.sampling()$MixtureSample
    }
    
    ggplotMixtureSampleSim <- function(){
      DFM <- data.frame(Latitude=SimMap()$DataSim$DataSim$gridx, Longitude=SimMap()$DataSim$DataSim$gridy,
                        Abundance.field=SimMap()$DataSim$DataSim$abundance.field)
      DFS <- data.frame(Latitude=c(Mixture.sampling()$MixtureSample$Latitude), 
                        Longitude=c(Mixture.sampling()$MixtureSample$Longitude),
                        Mixture=c(Mixture.sampling()$MixtureSample$Mixture))
      ggplot <- ggplot(DFM) + geom_tile(aes(Latitude, Longitude, fill = Abundance.field)) + 
        scale_fill_viridis_c(option = "turbo") + theme_bw() + coord_equal(ratio = 1) +
        xlab("Latitude") + ylab("Longitude") + ggtitle("Mixture sampling process") + 
        geom_point(data=DFS, aes(x=Latitude,y=Longitude, color=as.factor(Mixture))) + labs(fill="Values", color="Mixture") +
        theme(plot.title=element_text(color = "black", size = 16, face = "bold", hjust = 0.5))
      return(ggplot)
    }
    
    downloadFile(
      id = "ggplotMixtureSampleSim",
      logger = ss_userAction.Log,
      filenameroot = "ggplotMixtureSampleSim",
      aspectratio  = 1,
      downloadfxns = list(png  = ggplotRandomSampleSim,
                          txt  = dataggplotMixtureSampleSim,
                          csv  = dataggplotMixtureSampleSim)
    )
    
    downloadablePlot(
      id = "ggplotMixtureSampleSim",
      logger = ss_userAction.Log,
      filenameroot = "ggplotMixtureSampleSim",
      aspectratio  = 1,
      downloadfxns = list(png  = ggplotMixtureSampleSim,
                          txt  = dataggplotMixtureSampleSim,
                          csv  = dataggplotRandomSampleSim),
      visibleplot  = ggplotMixtureSampleSim)
    
    # Server_UploadingData_Section ----
    
    file.sample.read <- reactive({
        if(is.null(input$file.uploadData)) return(NULL)
        input$file.uploadData
    })
    
    output$fileUploadedSample <- reactive({
        return(!is.null(file.sample.read()))
    })
    outputOptions(output, 'fileUploadedSample', suspendWhenHidden=FALSE)
    
    datareadSample <- function(){
        file <- file.sample.read()
        ext <- tools::file_ext(file$datapath)
        req(file)
        validate(need(ext == c("csv","rds"), "Please upload a csv or rds file"))
        if(ext=="csv"){file.read <- read.csv(file$datapath)}
        else file.read <- readRDS(file$datapath)
        return(file.read)
        # return(as.data.frame(file.read))
    }
    
    quiltplotSampleReadRaster <- function(){
      DF <- as.data.frame(datareadSample())
      gl <- list()
      for(i in 3:ncol(DF)){
        # options(repr.plot.width = 3, repr.plot.height =3)
        if(is.numeric(DF[,i])){
          assign(paste0("g",i-2),
                 ggplot(DF) + geom_point(aes_string(x=names(DF)[1],y=names(DF)[2],colour=names(DF)[i])) +
                   scale_colour_viridis_c(option = "turbo") + labs(colour=paste(names(DF)[i],"\nValues")) +
                   theme_bw() + xlab(names(DF)[1]) + ylab(names(DF)[2]) +
                   ggtitle(names(DF)[i]) + theme(plot.title=element_text(color = "black", size = 16, face = "bold", hjust = 0.5))
          )
        } else{
          assign(paste0("g",i-2),
                 ggplot(DF) + geom_point(aes_string(x=names(DF)[1],y=names(DF)[2],colour=names(DF)[i])) +
                   theme_bw() + xlab(names(DF)[1]) + ylab(names(DF)[2]) +
                   ggtitle(names(DF)[i]) + theme(plot.title=element_text(color = "black", size = 16, face = "bold", hjust = 0.5))
          )
        }
        gl[[i-2]] <- eval(parse(text=paste0("g",i-2)))
      }
      
      if(floor((ncol(DF)-2)/2)==((ncol(DF)-2)/2)){
        LayoutMatrix <- matrix(c(1:(ncol(DF)-2)), byrow=TRUE, ncol=2)
      } else{
        LayoutMatrix <- matrix(c(1:(ncol(DF)-2), (ncol(DF)-2)), byrow=TRUE, ncol=2)
      }
      
      gt <- grid.arrange(grobs=gl, layout_matrix = LayoutMatrix)
      return(gt)
    }
    
    downloadableTable("table.read.sample",
                      logger=ss_userAction.Log,
                      filenameroot="table.read.sample",
                      downloaddatafxns=list(csv=datareadSample,
                                            tsv=datareadSample),
                      tabledata=datareadSample,
                      caption="Sample Data"
    )
    
    downloadFile(
        id = "quilt.plot.sample",
        logger = ss_userAction.Log,
        filenameroot = "quilt.plot.sample",
        aspectratio  = 1,
        downloadfxns = list(png = quiltplotSampleReadRaster)
    )

    downloadablePlot(
        id = "quilt.plot.sample",
        logger = ss_userAction.Log,
        filenameroot = "quilt.plot.sample",
        aspectratio  = 0.125,
        downloadfxns = list(png=quiltplotSampleReadRaster),
        visibleplot  = quiltplotSampleReadRaster
    )
    
    hegihtScaleSample <- function(){
      return(400*floor((ncol(as.data.frame(datareadSample()))-2)/2))
      }
    
    observe({
      output$main_plot2 <- renderPlot({plot(quiltplotSampleReadRaster())}, height=hegihtScaleSample())
    })
    
    #### Start Raster
    
    file.Raster.read <- reactive({
      if(is.null(input$file.uploadDataRaster)) return(NULL)
      return(input$file.uploadDataRaster)
    })
    
    output$fileUploadedRaster <- reactive({
      return(!is.null(file.Raster.read()))
    })
    outputOptions(output, 'fileUploadedRaster', suspendWhenHidden=FALSE)
    
    datareadRaster <- function(){
      if(is.null(file.Raster.read())){
        file <- file.Raster.read()
        ext <- tools::file_ext(file$datapath)
        req(file)
        validate(need(ext == c("csv","rds"), "Please upload a csv or rds file"))
        if(ext=="csv"){file.read <- read.csv(file$datapath)}
        else file.read <- readRDS(file$datapath)
        return(file.read)
      } else{return(NULL)}
    }
    
    quiltplotRasterReadRaster <- function(){
      DF <- as.data.frame(datareadRaster())
      
      if(!sum( colnames(DF) == make.unique(colnames(DF), sep="") ) == ncol(DF)){
        colnames(DF) <- make.unique(colnames(DF), sep="")
      }
      
      ind.x <- "x1"
      k <- 1
      pos.x <- c()
      while(ind.x %in% names(DF)){
        pos.x[k] <- which(ind.x==names(DF))
        k <- k+1
        ind.x <- paste("x", k, sep="")
      }
      
      pos.ext.x <- pos.x
      pos.ext.x[length(pos.x)+1] <- ncol(DF)+1
      
      gl <- list()
      for(j in pos.x){
        k <- which(j==pos.x)
        for(i in (j+2):(pos.ext.x[which(j==pos.x)+1]-1)){
          if(is.numeric(DF[,i])){
            assign(paste0("g",i-k*2),
                   ggplot(DF) + geom_point(aes_string(x=names(DF)[j],y=names(DF)[j+1],colour=names(DF)[i])) +
                     scale_colour_viridis_c(option = "turbo") + labs(colour=paste(names(DF)[i],"\nValues")) +
                     theme_bw() + xlab(names(DF)[j]) + ylab(names(DF)[j+1]) +
                     ggtitle(names(DF)[i]) + theme(plot.title=element_text(color = "black", size = 16, face = "bold", hjust = 0.5))
            )
          } else{
            assign(paste0("g",i-k*2),
                   ggplot(DF) + geom_point(aes_string(x=names(DF)[j],y=names(DF)[j+1],colour=names(DF)[i])) +
                     theme_bw() + xlab(names(DF)[j]) + ylab(names(DF)[j+1]) +
                     ggtitle(names(DF)[i]) + theme(plot.title=element_text(color = "black", size = 16, face = "bold", hjust = 0.5))
            )
          }
          gl[[i-k*2]] <- eval(parse(text=paste0("g",i-k*2)))
        }
      }
      
      if(floor((ncol(DF)-2*k)/2)==((ncol(DF)-2*k)/2)){
        LayoutMatrix <- matrix(c(1:(ncol(DF)-2*k)), byrow=TRUE, ncol=2)
      } else{
        LayoutMatrix <- matrix(c(1:(ncol(DF)-2*k), (ncol(DF)-2*k)), byrow=TRUE, ncol=2)
      }

      gt <- grid.arrange(grobs=gl, layout_matrix = LayoutMatrix)
      return(gt)
    }
    
    val <- reactiveVal(FALSE)

    observeEvent(input$EnhancementIrrGrid, {
        val(!val())
    })
    
    downloadableTable("table.read.raster",
                      logger=ss_userAction.Log,
                      filenameroot="table.read.raster",
                      downloaddatafxns=list(csv=datareadRaster,
                                            tsv=datareadRaster),
                      tabledata=datareadSample,
                      caption="Raster Data"
    )
    
    downloadFile(
      id = "quilt.plot.raster",
      logger = ss_userAction.Log,
      filenameroot = "quilt.plot.raster",
      aspectratio  = 1,
      downloadfxns = list(png  = quiltplotRasterReadRaster)
    )
    
    downloadablePlot(
      id = "quilt.plot.raster",
      logger = ss_userAction.Log,
      filenameroot = "quilt.plot.raster",
      aspectratio  = 0.125,
      downloadfxns = list(png=quiltplotRasterReadRaster),
      visibleplot  = quiltplotRasterReadRaster
    )
    
    hegihtScaleRaster <- function(){
      return(400*floor((ncol(as.data.frame(datareadRaster()))-2)/2))
    }
    
    observe({
      output$main_plot3 <- renderPlot({plot(quiltplotRasterReadRaster())}, height=hegihtScaleRaster())
    })
    
    # Server_IndependentModelling_Section ----
    
    IndCheckBoxNames <- function(){
      if(input$IndDataSimulatedLoaded=="load"){
        DF <- as.data.frame(datareadSample())
        if(input$SelectIndFamily=="binomial"){DFnames <- names(DF)[c(5:ncol(DF))]}
        else{DFnames <- names(DF)[c(4:ncol(DF))]}
      } else if(input$IndDataSimulatedLoaded=="sim"){
        DF <- Ind.sampling()
        DFnames <- names(DF)[c(4)]
      }
      return(DFnames)
    }
    
    observe({
      output$checkBoxIndDataFrame <- renderUI({
        tagList(
          checkboxGroupInput(inputId="UserComponentsInd",
                             label="User Defined Components",
                             choices=IndCheckBoxNames(),
                             selected=c()
          )
        )
      })
    })
    
    observe({
      if(input$IndDataSimulatedLoaded=="load"){
        output$SelectIndEffectCov <- renderUI({
          
          if(length(input$UserComponentsInd)>0){
            box(id="IndCovariateEffects", width=12, title="Covariate Effects",
                lapply(seq_along(input$UserComponentsInd), function(i){
                  list(
                    selectInput(inputId=paste0("IndEffectCov",i), label=tags$span(style="color: blue; font-weight: bold;", paste(input$UserComponentsInd[i],"Effect")) ,
                                choices = list("Linear (or reference level)" = "linear", "Random Walk 1" = "rw1", "Random Walk 2" = "rw2", "SPDE 1" = "spde1", "IID"="iid"),
                                selected = "linear"),
                    conditionalPanel(condition=paste0("input.IndEffectCov",i,"=='rw1'||","input.IndEffectCov",i,"=='rw2'||","input.IndEffectCov",i,"=='spde1'"),
                                     numericInput(inputId=paste0("IndEffectCovNodes",i),
                                                  label="Number of nodes",
                                                  value=10, min=1, step=1
                                     )
                    ),
                    radioGroupButtons(
                      inputId = paste0("IndEffectCustomPrior",i),
                      label = "Custom Prior",
                      choices = list("Auto" = "auto", "Custom" = "custom"),
                      status = "primary"
                    ),
                    conditionalPanel(condition=paste0("input.IndEffectCustomPrior",i,"=='custom'"),
                                     list(
                                       selectInput(inputId = paste0("IndEffectCovKindPrior",i),
                                                   label = "Prior distribution",
                                                   choices = list("Base" = "base", "PC prior" = "pc", "Uniform" = "unif", "Flat Uniform" = "flatunif")),
                                       conditionalPanel(condition=paste0("input.IndEffectCovKindPrior",i,"=='base'"),
                                                        textInput(inputId = paste0("IndEffectCovPriorBaseValues",i),
                                                                  label = "Base Prior Values",
                                                                  value = "0, 1e-3")),
                                       conditionalPanel(condition=paste0("input.IndEffectCovKindPrior",i,"=='pc'"),
                                                        textInput(inputId = paste0("IndEffectCovPriorPCValues",i),
                                                                  label = "PC-Prior Values",
                                                                  value = "0, 1e-3")),
                                       conditionalPanel(condition=paste0("input.IndEffectCovKindPrior",i,"=='unif'"),
                                                        textInput(inputId = paste0("IndEffectCovPriorUnif",i),
                                                                  label = "Uniform Lower and upper values",
                                                                  value = "0, 10",
                                                                  placeholder = "Lower and upper values: 'U(a,b)'"))
                                     )
                    ),
                    conditionalPanel(condition=paste0("input.IndEffectCov",i,"=='rw1'||",
                                                      "input.IndEffectCov",i,"=='rw2'||",
                                                      "input.IndEffectCov",i,"=='iid'"),
                                     radioGroupButtons(
                                       inputId = paste0("IndEffectCovConstr",i),
                                       label = "Sum to zero",
                                       choices = list("TRUE" = "TRUE", "FALSE" = "FALSE"),
                                       status = "primary"
                                     ))
                  )
                }))
          } else{}
        }) } else if(input$IndDataSimulatedLoaded=="sim"){
          output$SelectIndEffectCov <- renderUI({
            if(length(input$UserComponentsInd)>0){
              box(id="IndCovariateEffects", width=12, title="Covariate Effects",
                  lapply(seq_along(input$UserComponentsInd), function(i){
                    list(
                      selectInput(inputId=paste0("IndEffectCov",i), label=tags$span(style="color: blue; font-weight: bold;", paste(input$UserComponentsInd[i],"Effect")) ,
                                  choices = list("Linear (or reference level)" = "linear", "Random Walk 1" = "rw1", "Random Walk 2" = "rw2", "SPDE 1" = "spde1", "IID"="iid"),
                                  selected = "linear"),
                      conditionalPanel(condition=paste0("input.IndEffectCov",i,"=='rw1'||","input.IndEffectCov",i,"=='rw2'||","input.IndEffectCov",i,"=='spde1'"),
                                       numericInput(inputId=paste0("IndEffectCovNodes",i),
                                                    label="Number of nodes",
                                                    value=10, min=1, step=1
                                       )
                      ),
                      radioGroupButtons(
                        inputId = paste0("IndEffectCustomPrior",i),
                        label = "Custom Prior",
                        choices = list("Auto" = "auto", "Custom" = "custom"),
                        status = "primary"
                      ),
                      conditionalPanel(condition=paste0("input.IndEffectCustomPrior",i,"=='custom'"),
                                       list(
                                         selectInput(inputId = paste0("IndEffectCovKindPrior",i),
                                                     label = "Prior distribution",
                                                     choices = list("Base" = "base", "PC prior" = "pc", "Uniform" = "unif", "Flat Uniform" = "unifflat")),
                                         conditionalPanel(condition=paste0("input.IndEffectCovKindPrior",i,"=='base'"),
                                                          textInput(inputId = paste0("IndEffectCovPriorBaseValues",i),
                                                                    label = "Base Prior Values",
                                                                    value = "0, 1e-3")),
                                         conditionalPanel(condition=paste0("input.IndEffectCovKindPrior",i,"=='pc'"),
                                                          textInput(inputId = paste0("IndEffectCovPriorPCValues",i),
                                                                    label = "PC-Prior Values",
                                                                    value = "0, 1e-3")),
                                         conditionalPanel(condition=paste0("input.IndEffectCovKindPrior",i,"=='unif'"),
                                                          textInput(inputId = paste0("IndEffectCovPriorUnif",i),
                                                                    label = "Uniform Lower and upper values",
                                                                    value = "0, 10",
                                                                    placeholder = "Lower and upper values: 'U(a,b)'"))
                                       )
                      ),
                      conditionalPanel(condition=paste0("input.IndEffectCov",i,"=='rw1'||",
                                                        "input.IndEffectCov",i,"=='rw2'||",
                                                        "input.IndEffectCov",i,"=='iid'"),
                                       radioGroupButtons(
                                         inputId = paste0("IndEffectCovConstr",i),
                                         label = "Sum to zero",
                                         choices = list("TRUE" = "TRUE", "FALSE" = "FALSE"),
                                         status = "primary"
                                       ))
                    )
                  }))
            } else{}
          }) }
    })
    
    observe({
      if(input$IndDataSimulatedLoaded=="load"){
        DF <- as.data.frame(datareadSample())
        if(length(input$UserComponentsInd)>0){
          IndUserComponent <- input$UserComponentsInd
          DF2 <- select(DF, IndUserComponent[!as.vector(unlist(lapply(X=select(DF,IndUserComponent), FUN=is.numeric)))])
          output$SelectLoadIndEffectCovFactorPred <- renderUI({
            if(ncol(DF2)>0){
              box(id="IndPredFactorLevel", width=12, title="Prediction Factor Level",
                  lapply(seq_along(names(DF2)), function(i){
                    choices <- unique(DF2[[names(DF2)[i]]])
                    list(
                      radioGroupButtons(inputId = paste0("IndKindPredictionFactorLevel",i), label = tags$span(style="color: blue; font-weight: bold;", paste(names(DF2)[i], "(prediction protocol)")),
                                        choices = c("Reference level" = "reference", "Nearest level" = "nearest"), status = "success", justified = TRUE),
                      conditionalPanel(
                        condition=paste0("input.IndKindPredictionFactorLevel",i,"=='reference'"),
                        selectInput(inputId=paste0("IndEffectCovFactorPred",i), label=paste(names(DF2)[i],"Reference Factor"),
                                    choices = choices, selected = choices[1])
                      )
                    )
                  }))
            }
          })
        } else{
          output$SelectLoadIndEffectCovFactorPred <- renderUI({})
        }
      }
    })
    
    JointPriorPreviewInd <- eventReactive(input$Indpreviewpriordistributions,{
      d <- 2; nu <- 1
      rhoinputPC <- as.numeric(unlist(strsplit(input$Indrangepcpriorprev, ",")))
      sigmainputPC <- as.numeric(unlist(strsplit(input$Indsigmapcpriorprev, ",")))
      rhoPC <- seq(rhoinputPC[1], rhoinputPC[2], length.out=rhoinputPC[3])
      sigmaPC <- seq(sigmainputPC[1], sigmainputPC[2], length.out=sigmainputPC[3])
      rhosigmaPC <- expand.grid(rho=rhoPC, sigma=sigmaPC)
      rho0PC <- rhoinputPC[4]; alpha1 <- rhoinputPC[5]
      sigma0PC <- sigmainputPC[4]; alpha2 <- sigmainputPC[5]
      
      lambdarhoPC <- -log(alpha1)*rho0PC**(d/2) #-(rho0/sqrt(8*nu))**(d/2)*log(alpha1)
      lambdasigmaPC <- -log(alpha2)/sigma0PC #-(sqrt(8*nu)/rhosigma$rho)**(-nu)*sqrt(gamma(nu)/(gamma(nu+d/2)*(4*pi)**(d/2)))*log(alpha2)/sigma0
      pirhosigmaPCM <- d/2*lambdarhoPC*rhosigmaPC$rho**(-1-d/2)*exp(-lambdarhoPC*rhosigmaPC$rho**(-d/2))*lambdasigmaPC*exp(-lambdasigmaPC*rhosigmaPC$sigma)
      
      probIntPC <- diff(range(rhoPC))*diff(range(sigmaPC))/length(pirhosigmaPCM)*sum(pirhosigmaPCM)
      
      rhoinputBase <- as.numeric(unlist(strsplit(input$Indrangebasepriorprev, ",")))
      sigmainputBase <- as.numeric(unlist(strsplit(input$Indsigmabasepriorprev, ",")))
      rho <- seq(rhoinputBase[1], rhoinputBase[2], length.out=rhoinputBase[3])
      sigma <- seq(sigmainputBase[1], sigmainputBase[2], length.out=sigmainputBase[3])
      rhosigma <- expand.grid(rho=rho, sigma=sigma)
      rho0 <- rhoinputBase[4]; sigma0 <- sigmainputBase[4]
      meanrho <- rhoinputBase[5]; meansigma <- sigmainputBase[5]
      sdrho <- rhoinputBase[6]; sdsigma <- sigmainputBase[6]
      pirho <- (sqrt(2*pi*sdrho**2)*rhosigma$rho)**(-1)*exp(-(log(rhosigma$rho/rho0)**2-2*log(rhosigma$rho/rho0)*meanrho+meanrho**2)/(2*sdrho**2))
      pisigma <- (sqrt(2*pi*sdsigma**2)*rhosigma$sigma)**(-1)*exp(-(log(rhosigma$sigma/sigma0)**2-2*log(rhosigma$sigma/sigma0)*meansigma+meansigma**2)/(2*sdsigma**2))
      pirhosigmaM <- pirho*pisigma
      
      probInt <- sum(pirhosigmaM)*diff(range(rhosigma$rho))*diff(range(rhosigma$sigma))/length(pirhosigmaM)
      
      ListPreview <- list(rhosigmaPC=rhosigmaPC, pirhosigmaPCM=pirhosigmaPCM, probIntPC=probIntPC,
                          rhosigma=rhosigma, pirhosigmaM=pirhosigmaM, probInt=probInt)
      
      return(ListPreview)
    })
    
    IndPreviewJointPriorPlot <- function(){
      ggplotPCprior <- ggplot() +
        geom_tile(data=data.frame(rho=JointPriorPreviewInd()$rhosigmaPC$rho, sigma=JointPriorPreviewInd()$rhosigmaPC$sigma, pi=JointPriorPreviewInd()$pirhosigmaPCM),
                  mapping=aes(x=rho, y=sigma, fill=pi)) +
        labs(title=paste("Joint PC-Prior. Cumulative Prob.=", round(JointPriorPreviewInd()$probIntPC, digits=2)),
             x=expression(rho), y=expression(sigma)) +
        scale_fill_viridis_c(option="turbo") + theme_bw() + theme(plot.title=element_text(face='bold'))
      
      ggplotENpriorAnalytic <- ggplot() +
        geom_tile(data=data.frame(rho=JointPriorPreviewInd()$rhosigma$rho, sigma=JointPriorPreviewInd()$rhosigma$sigma, pi=JointPriorPreviewInd()$pirhosigmaM),
                  mapping=aes(x=rho, y=sigma, fill=pi)) +
        scale_fill_viridis_c(option="turbo") +
        labs(title=paste("Joint Base Prior. Cumulative Prob.=", round(JointPriorPreviewInd()$probInt, digits=2)),
             x=expression(rho), y=expression(sigma)) +
        theme_bw() + theme(plot.title=element_text(face='bold'))
      
      ggplotJointPriorPreview <- ggplotPCprior + ggplotENpriorAnalytic
      return(ggplotJointPriorPreview)
    }
    
    downloadFile(
      id = "IndPreviewJointPriorPlot",
      logger = ss_userAction.Log,
      filenameroot = "IndPreviewJointPriorPlot",
      aspectratio  = 1,
      downloadfxns = list(png  = IndPreviewJointPriorPlot)
    )
    
    downloadablePlot(id = "IndPreviewJointPriorPlot",
                     logger = ss_userAction.Log,
                     filenameroot = "IndPreviewJointPriorPlot",
                     aspectratio  = 1,
                     downloadfxns = list(png=IndPreviewJointPriorPlot),
                     visibleplot  = IndPreviewJointPriorPlot)
    
    ## Mesh construction for IndependentModel ====
    
    IndMeshBase <- reactive({
      if(input$IndDataSimulatedLoaded=="sim"){
        DFsample <- as.data.frame(Ind.sampling())
        convexhull <- chull(DFsample[,1:2])
        convexhull <- c(convexhull, convexhull[1])
        qloc <- quantile(as.vector(dist(DFsample[sample(1:nrow(DFsample),size=min(c(50,nrow(DFsample)))),1:2])),probs=c(0.03,0.3))
        mesh <- fm_mesh_2d_inla(loc.domain=DFsample[convexhull,1:2], cutoff = qloc[1]/2, offset=c(-0.1, -0.2),
                             max.edge=c(qloc[1], qloc[2]))
        sample <- DFsample
      } else if(input$IndDataSimulatedLoaded=="load"){
        DFsample <- datareadSample()
        if(input$IndRasterSPDE=="raster"){
          rasterSample <- datareadRaster()[sample(1:nrow(datareadRaster()), min(c(50,nrow(datareadRaster())))),1:2]
          qloc <- quantile(as.vector(dist(rasterSample)),probs=c(0.03,0.3))
          convexhull <- chull(datareadRaster()[,1:2])
          convexhull <- c(convexhull, convexhull[1])
          mesh <- fm_mesh_2d_inla(loc.domain=datareadRaster()[convexhull,1:2], cutoff = qloc[1]/2, offset=c(-0.1, -0.2),
                               max.edge=c(qloc[1], qloc[2]))
          sample <- rasterSample
        } else if(input$IndRasterSPDE=="solvecov"){
          convexhull <- chull(rbind(DFsample[,1:2],datareadRaster()[,1:2]))
          convexhull <- c(convexhull, convexhull[1])
          qloc <- quantile(as.vector(dist(DFsample[sample(1:nrow(DFsample),size=min(c(50,nrow(DFsample)))),1:2])),probs=c(0.03,0.3))
          mesh <- fm_mesh_2d_inla(loc.domain=rbind(DFsample[,1:2],datareadRaster()[,1:2])[convexhull,], cutoff = qloc[1]/2, offset=c(-0.1, -0.2),
                               max.edge=c(qloc[1], qloc[2]))
          sample <- DFsample
        }
      }
      result <- list(mesh=mesh, qloc=qloc, Sample=sample)
      return(result)
    })
    
    IndMesh <- eventReactive(input$buildIndMesh, {
      if(input$IndDataSimulatedLoaded=="sim"){
        DFsample <- as.data.frame(Ind.sampling())
      } else if(input$IndDataSimulatedLoaded=="load"){
        DFsample <- as.data.frame(datareadSample())
        DFraster <- try(datareadRaster(), silent=TRUE)
        if(class(DFraster)!="try-error"){
          DFsample <- rbind(DFsample[,1:2], DFraster[,1:2])
        }
      }
      
      interiorNonConvex <- function(x, condition, convex, resolution, file.read){
        if(condition=="interiorIndMeshnonconvex"){
          interior <- inla.nonconvex.hull(points=x, convex=convex, resolution=resolution)
        } else if(condition=="interiorIndMeshcustomboundary"){
          ext <- tools::file_ext(file.read$datapath)
          validate(need(ext == c("csv","rds"), "Please upload a csv or rds file"))
          if(ext=="csv"){coords <- read.csv(file.read$datapath)}
          else coords <- readRDS(file.read$datapath)
          innerBorderMesh <- SpatialPolygons(Srl=list(Polygons(srl=list(Polygon(coords=coords)), ID="interiorBoundary")))
          interior <- inla.sp2segment(sp=innerBorderMesh)
        } else{interior <- NULL}
        return(interior)
      }
      
      boundaryNonConvex <- function(x, condition, convex, resolution, file.read){
        if(condition=="IndMeshnonconvex"){
          boundary <- inla.nonconvex.hull(points=x, convex=convex, resolution=resolution)
        } else if(condition=="IndMeshcustomboundary"){
          ext <- tools::file_ext(file.read$datapath)
          validate(need(ext == c("csv","rds"), "Please upload a csv or rds file"))
          if(ext=="csv"){coords <- read.csv(file.read$datapath)}
          else coords <- readRDS(file.read$datapath)
          outerBorderMesh <- SpatialPolygons(Srl=list(Polygons(srl=list(Polygon(coords=coords)), ID="externalBoundary")))
          boundary <- inla.sp2segment(sp=outerBorderMesh)
        } else{boundary <- NULL}
        return(boundary)
      }
      
      if(input$selectionIndMesh=="qlocation"){
        if(input$IndRasterSPDE=="raster"){
          rasterSample <- datareadRaster()[sample(1:nrow(datareadRaster()), min(c(50,nrow(datareadRaster())))),1:2]
          qloc <- quantile(as.vector(dist(rasterSample)),probs=c(ifelse(input$IndCustomMesh,input$IndMeshQloc,0.03),0.3))
          convexhull <- chull(DFsample[,1:2])
          convexhull <- c(convexhull, convexhull[1])
          mesh <- fm_mesh_2d_inla(loc.domain=DFsample[convexhull,1:2], cutoff = qloc[1]/2, offset=c(-0.1, -0.2), max.edge=c(qloc[1], qloc[2]),
                               boundary=list(interiorNonConvex(x=as.matrix(rasterSample[,1:2]), condition=input$interiorIndMesh, convex=input$interiorcurvatureIndMesh, resolution=input$interiorresolutionIndMesh, file.read=input$interiorshapefileIndMesh),
                                             boundaryNonConvex(x=as.matrix(rasterSample[,1:2]), condition=input$boundaryIndMesh, convex=input$curvatureIndMesh, resolution=input$resolutionIndMesh, file.read=input$shapefileIndMesh)))
          sample <- rasterSample
        } else if(input$IndRasterSPDE=="solvecov"){
          qloc <- quantile(as.vector(dist(DFsample[sample(1:nrow(DFsample),size=min(c(50,nrow(DFsample)))),1:2])),probs=c(ifelse(input$IndCustomMesh,input$IndMeshQloc,0.03),0.3))
          convexhull <- chull(DFsample[,1:2])
          convexhull <- c(convexhull, convexhull[1])
          mesh <- fm_mesh_2d_inla(loc.domain=DFsample[convexhull,1:2], cutoff = qloc[1]/2, offset=c(-0.1, -0.2), max.edge=c(qloc[1], qloc[2]),
                               boundary=list(interiorNonConvex(x=as.matrix(DFsample[,1:2]), condition=input$interiorIndMesh, convex=input$interiorcurvatureIndMesh, resolution=input$interiorresolutionIndMesh, file.read=input$interiorshapefileIndMesh),
                                             boundaryNonConvex(x=as.matrix(DFsample[,1:2]), condition=input$boundaryIndMesh, convex=input$curvatureIndMesh, resolution=input$resolutionIndMesh, file.read=input$shapefileIndMesh)))
          sample <- DFsample
        }
      } else if(input$selectionIndMesh=="edgelength"){
        if(input$IndRasterSPDE=="raster"){
          rasterSample <- datareadRaster()
          qloc <- input$EdgeLengthIndMesh
          convexhull <- chull(DFsample[,1:2])
          convexhull <- c(convexhull, convexhull[1])
          mesh <- fm_mesh_2d_inla(loc.domain=DFsample[convexhull,1:2], max.edge=c(1,1.5)*input$EdgeLengthIndMesh, cutoff=input$EdgeLengthIndMesh/5, offset=c(-0.1, -0.2), max.edge=c(qloc[1], qloc[2]),
                               boundary=list(interiorNonConvex(x=as.matrix(rasterSample[,1:2]), condition=input$interiorIndMesh, convex=input$interiorcurvatureIndMesh, resolution=input$interiorresolutionIndMesh, file.read=input$interiorshapefileIndMesh),
                                             boundaryNonConvex(x=as.matrix(rasterSample[,1:2]), condition=input$boundaryIndMesh, convex=input$curvatureIndMesh, resolution=input$resolutionIndMesh, file.read=input$shapefileIndMesh)))
          sample <- rasterSample
        } else if(input$IndRasterSPDE=="solvecov"){
          qloc <- input$EdgeLengthIndMesh
          convexhull <- chull(DFsample[,1:2])
          convexhull <- c(convexhull, convexhull[1])
          mesh <- fm_mesh_2d_inla(loc.domain=DFsample[convexhull,1:2], max.edge=c(1,1.5)*input$EdgeLengthIndMesh,  cutoff=input$EdgeLengthIndMesh/5, offset=c(-0.1, -0.2),
                               boundary=list(interiorNonConvex(x=as.matrix(DFsample[,1:2]), condition=input$interiorIndMesh, convex=input$interiorcurvatureIndMesh, resolution=input$interiorresolutionIndMesh, file.read=input$interiorshapefileIndMesh),
                                             boundaryNonConvex(x=as.matrix(DFsample[,1:2]), condition=input$boundaryIndMesh, convex=input$curvatureIndMesh, resolution=input$resolutionIndMesh, file.read=input$shapefileIndMesh)))
          sample <- DFsample
        }
      }
      
      result <- list(mesh=mesh, qloc=qloc, Sample=sample)
      return(result)
    })
    
    ggplotIndMesh <- function(){
      if(input$buildIndMesh==0){
        ggplot()+ gg(IndMeshBase()$mesh)+ theme_bw() + xlab("Latitude") + ylab("Longitude") +
          ggtitle("Mesh over the study region") +
          geom_point(data=IndMeshBase()$Sample,
                     aes(x=IndMeshBase()$Sample[,1],y=IndMeshBase()$Sample[,2]), size=1) +
          theme(plot.title=element_text(color = "black", size = 16, face = "bold", hjust = 0.5))
      } else{
        ggplot()+ gg(IndMesh()$mesh)+ theme_bw() + xlab("Latitude") + ylab("Longitude") +
          ggtitle("Mesh over the study region") +
          geom_point(data=IndMesh()$Sample,
                     aes(x=IndMesh()$Sample[,1],y=IndMesh()$Sample[,2]), size=1) +
          theme(plot.title=element_text(color = "black", size = 16, face = "bold", hjust = 0.5))
      }
      
    }
    
    downloadFile(
      id = "ggplotIndMesh",
      logger = ss_userAction.Log,
      filenameroot = "ggplotIndMesh",
      aspectratio  = 1,
      downloadfxns = list(png  = ggplotIndMesh)
    )
    
    downloadablePlot(
      id = "ggplotIndMesh",
      logger = ss_userAction.Log,
      filenameroot = "ggplotIndMesh",
      aspectratio  = 1,
      downloadfxns = list(png  = ggplotIndMesh),
      visibleplot  = ggplotIndMesh)
    
    ## Independent modelling code ====
    
    IndModelFit <- eventReactive(input$fitInd, {
      showNotification(ui=paste("Fitting the data."), duration = NULL)
      t1 <- Sys.time()
      #taking the data from simulation or from the loading tab
      if(input$IndDataSimulatedLoaded=="sim"){
        DFsample <- as.data.frame(Ind.sampling())
      } else if(input$IndDataSimulatedLoaded=="load"){
        DFsample <- as.data.frame(datareadSample())
        if(input$IndRasterSPDE=="raster"){
          DFraster <- as.data.frame(datareadRaster())
        }
      }
      
      variablesChosenDefault <- input$DefaultComponentsInd
      variablesChosenUser <- input$UserComponentsInd
      variablesChosen <- c(variablesChosenDefault, variablesChosenUser)
      
      if(input$buildIndMesh==0){
        server_mesh <- IndMeshBase()
        mesh <- server_mesh$mesh
      } else{
        server_mesh <- IndMesh()
        mesh <- server_mesh$mesh
      }
      
      if(input$optionIndS=="auto"){
        if(input$KindPriorSpatialEffectInd=="PC.prior"){
          prior.range <- c(mean(c(diff(range(mesh$loc[mesh$segm$int$idx[,1],1])),diff(range(mesh$loc[mesh$segm$int$idx[,1],2]))))/5, 0.5)
          prior.sigma <- c(1,0.5)
          spde <- inla.spde2.pcmatern(mesh, prior.range = prior.range, prior.sigma = prior.sigma, alpha=2)
        } else{
          prior.range <- mean(c(diff(range(mesh$loc[mesh$segm$int$idx[,1],1])),diff(range(mesh$loc[mesh$segm$int$idx[,1],2]))))/5
          prior.sigma <- 1
          alpha <- 2; d <- 2
          nu <-  alpha - d/2
          kappa0 <-  log(8*nu)/2 -log(prior.range)
          tau0 <-  0.5*(lgamma(nu) - lgamma(nu + d/2) - d/2*log(4*pi)) - nu*kappa0 - log(prior.sigma)
          spde <-  inla.spde2.matern(mesh = mesh, B.tau = cbind(tau0, nu, -1), B.kappa = cbind(kappa0, -1, 0),
                                     theta.prior.mean = c(0.5,0.5), theta.prior.prec = c(1,1))
        }
      } else if(input$optionIndS=="custom"){
        if(input$KindPriorSpatialEffectInd=="PC.prior"){
          prior.range <- as.numeric(unlist(strsplit(input$IndPriorRangePC, split=",")))
          prior.sigma <- as.numeric(unlist(strsplit(input$IndPriorStdevPC, split=",")))
          spde <- inla.spde2.pcmatern(mesh, prior.range = prior.range, prior.sigma = prior.sigma, alpha=2)
        } else{
          prior.range <- as.numeric(unlist(strsplit(input$IndPriorRangeBase, split=",")))
          prior.sigma <- as.numeric(unlist(strsplit(input$IndPriorStdevBase, split=",")))
          alpha <- 2; d <- 2
          nu <-  alpha - d/2
          kappa0 <-  log(8*nu)/2 -log(prior.range[1])
          tau0 <-  0.5*(lgamma(nu) - lgamma(nu + d/2) - d/2*log(4*pi)) - nu*kappa0 - log(prior.sigma[1])
          spde <-  inla.spde2.matern(mesh = mesh, B.tau = cbind(tau0, nu, -1), B.kappa = cbind(kappa0, -1, 0),
                                     theta.prior.mean = c(prior.range[2],prior.sigma[2]), theta.prior.prec = c(prior.range[3],prior.sigma[3]))
        }
      }
      
      spde.index <- inla.spde.make.index(name="Spatial", n.spde = spde$n.spde)
      
      # LGCP mesh operations
      
      n <- nrow(DFsample)
      lmat <- inla.spde.make.A(mesh, as.matrix(DFsample[,1:2]))
      A.inf <- lmat
      
      ### Prediction of covariates ====
      
      List.covariates.inf <- list()
      List.covariates.pred <- list()
      
      
      if((input$IndDataSimulatedLoaded=="load"&input$IndRasterSPDE=="solvecov")|input$IndDataSimulatedLoaded=="sim"|(input$IndDataSimulatedLoaded=="load"&input$IndRasterSPDE=="raster"|input$IndRasterPred=="SPDEraster")){
        x.pred <- seq(range(mesh$loc[mesh$segm$int$idx[,1], 1])[1], range(mesh$loc[mesh$segm$int$idx[,1], 1])[2], length.out=input$IndSPDErasterInput)
        y.pred <- seq(range(mesh$loc[mesh$segm$int$idx[,1], 2])[1], range(mesh$loc[mesh$segm$int$idx[,1], 2])[2], length.out=input$IndSPDErasterInput)
        xy.pred <- expand.grid(x=x.pred,y=y.pred)
        xy.pred <- xy.pred[which(!is.na(over(SpatialPoints(coords=xy.pred),SpatialPolygons(Srl=list(Polygons(srl=list(Polygon(coords=mesh$loc[mesh$segm$int$idx[,1], 1:2])), ID="int")))))),1:2]
        A.geo.pred <- inla.spde.make.A(mesh=mesh, loc=as.matrix(xy.pred))
      } else{
        xy.pred <- DFraster
        A.geo.pred <- inla.spde.make.A(mesh=mesh, loc=as.matrix(xy.pred))
      }
      
      for(i in seq_along(variablesChosenUser)){
        if(!is.character(DFsample[,variablesChosenUser[i]])){
          prior.range.cov <- c(mean(c(diff(range(mesh$loc[mesh$segm$int$idx[,1],1])),diff(range(mesh$loc[mesh$segm$int$idx[,1],2]))))/5, 0.5)
          prior.sigma.cov <- c(1,0.5)
          spde.cov <- inla.spde2.pcmatern(mesh, prior.range = prior.range.cov, prior.sigma = prior.sigma.cov, alpha=2)
          spde.cov.index <- inla.spde.make.index(name="spatial.cov", n.spde = spde.cov$n.spde)
          formula.cov <- y ~ -1 + Intercept + f(spatial.cov, model=spde.cov)
          
          if((input$IndDataSimulatedLoaded=="load"&input$IndRasterSPDE=="solvecov")|input$IndDataSimulatedLoaded=="sim"){
            x.pred <- seq(range(mesh$loc[mesh$segm$int$idx[,1], 1])[1], range(mesh$loc[mesh$segm$int$idx[,1], 1])[2], length.out=input$IndSPDErasterInput)
            y.pred <- seq(range(mesh$loc[mesh$segm$int$idx[,1], 2])[1], range(mesh$loc[mesh$segm$int$idx[,1], 2])[2], length.out=input$IndSPDErasterInput)
            xy.pred <- expand.grid(x=x.pred,y=y.pred)
            xy.pred <- xy.pred[which(!is.na(over(SpatialPoints(coords=xy.pred),SpatialPolygons(Srl=list(Polygons(srl=list(Polygon(coords=mesh$loc[mesh$segm$int$idx[,1], 1:2])), ID="int")))))),1:2]
            A.geo.pred <- inla.spde.make.A(mesh=mesh, loc=as.matrix(xy.pred))
            Inf.stack.cov <- inla.stack(data=list(y=DFsample[,variablesChosenUser[i]]),
                                        A=list(lmat, 1),
                                        effects=list(list(spatial.cov=spde.cov.index$spatial.cov),
                                                     list(Intercept=rep(1,nrow(DFsample)))
                                        ),
                                        tag="Inference.cov")
            Pred.stack.cov <- inla.stack(data=list(y=rep(NA,nrow(xy.pred))),
                                         A=list(A.geo.pred, 1),
                                         effects=list(list(spatial.cov=spde.cov.index$spatial.cov),
                                                      list(Intercept=rep(1,nrow(xy.pred)))
                                         ),
                                         tag="Prediction.cov")
            Total.stack.cov <- inla.stack(Inf.stack.cov, Pred.stack.cov)
            mod.cov <- inla(formula=formula.cov, data=inla.stack.data(Total.stack.cov), family="gaussian",
                            control.predictor = list(compute=FALSE, A=inla.stack.A(Total.stack.cov)))
            indx.inf <- inla.stack.index(Total.stack.cov, tag="Inference.cov")$data
            indx.pred <- inla.stack.index(Total.stack.cov, tag="Prediction.cov")$data
            List.covariates.inf[[variablesChosenUser[i]]] <- as.numeric(!is.na(DFsample[,variablesChosenUser[i]]))*DFsample[,variablesChosenUser[i]] + as.numeric(is.na(DFsample[,variablesChosenUser[i]]))*mod.cov$summary.fitted.values[indx.inf,"mean"]
            List.covariates.pred[[variablesChosenUser[i]]] <- mod.cov$summary.fitted.values[indx.pred,"mean"]
          } else if(input$IndDataSimulatedLoaded=="load"&input$IndRasterSPDE=="raster"|input$IndRasterPred=="SPDEraster"){
            x.pred <- seq(range(mesh$loc[mesh$segm$int$idx[,1], 1])[1], range(mesh$loc[mesh$segm$int$idx[,1], 1])[2], length.out=input$IndSPDErasterInput)
            y.pred <- seq(range(mesh$loc[mesh$segm$int$idx[,1], 2])[1], range(mesh$loc[mesh$segm$int$idx[,1], 2])[2], length.out=input$IndSPDErasterInput)
            xy.pred <- expand.grid(x=x.pred,y=y.pred)
            xy.pred <- xy.pred[which(!is.na(over(SpatialPoints(coords=xy.pred),SpatialPolygons(Srl=list(Polygons(srl=list(Polygon(coords=mesh$loc[mesh$segm$int$idx[,1], 1:2])), ID="int")))))),1:2]
            A.geo.pred <- inla.spde.make.A(mesh=mesh, loc=as.matrix(xy.pred))
            if(variablesChosenUser[i] %in% colnames(DFraster)){
              Inf.stack.cov <- inla.stack(data=list(y=c(DFsample[,variablesChosenUser[i]], DFraster[!is.na(DFraster[,variablesChosenUser[i]]),variablesChosenUser[i]])),
                                          A=list(inla.spde.make.A(mesh=mesh, loc=as.matrix(rbind(DFsample[,1:2],DFraster[!is.na(DFraster[,variablesChosenUser[i]]),1:2]))), 1),
                                          effects=list(list(spatial.cov=spde.cov.index$spatial.cov),
                                                       list(Intercept=rep(1,nrow(DFsample)+nrow(DFraster)))
                                          ),
                                          tag="Inference.cov")
            } else{
              Inf.stack.cov <- inla.stack(data=list(y=DFsample[,variablesChosenUser[i]]),
                                          A=list(inla.spde.make.A(mesh=mesh, loc=as.matrix(DFsample[,1:2])), 1),
                                          effects=list(list(spatial.cov=spde.cov.index$spatial.cov),
                                                       list(Intercept=rep(1,nrow(DFsample)))
                                          ),
                                          tag="Inference.cov")
            }
            Pred.stack.cov <- inla.stack(data=list(y=rep(NA,nrow(xy.pred))),
                                         A=list(A.geo.pred, 1),
                                         effects=list(list(spatial.cov=spde.cov.index$spatial.cov),
                                                      list(Intercept=rep(1,nrow(xy.pred)))
                                         ),
                                         tag="Prediction.cov")
            Total.stack.cov <- inla.stack(Inf.stack.cov, Pred.mesh.stack.cov, Pred.stack.cov)
            mod.cov <- inla(formula=formula.cov, data=inla.stack.data(Total.stack.cov), family="gaussian",
                            control.predictor = list(compute=FALSE, A=inla.stack.A(Total.stack.cov)))
            indx.inf <- inla.stack.index(Total.stack.cov, tag="Inference.cov")$data
            indx.pred <- inla.stack.index(Total.stack.cov, tag="Prediction.cov")$data
            List.covariates.inf[[variablesChosenUser[i]]] <- as.numeric(!is.na(DFsample[,variablesChosenUser[i]]))*DFsample[,variablesChosenUser[i]] + (as.numeric(is.na(DFsample[,variablesChosenUser[i]]))*(mod.cov$summary.fitted.values[indx.inf,"mean"])[1:nrow(DFsample)])
            List.covariates.pred[[variablesChosenUser[i]]] <- mod.cov$summary.fitted.values[indx.pred,"mean"]
          } else if(input$IndDataSimulatedLoaded=="load"&input$IndRasterSPDE=="raster"|input$IndRasterPred=="rasterpred"){
            xy.pred <- DFraste[,1:2]
            A.geo.pred <- inla.spde.make.A(mesh=mesh, loc=as.matrix(xy.pred))
            
            Inf.stack.cov <- inla.stack(data=list(y=c(DFsample[,variablesChosenUser[i]], DFraster[,variablesChosenUser[i]])),
                                        A=list(inla.spde.make.A(mesh=mesh, loc=as.matrix(rbind(DFsample[,1:2],DFraster[,1:2]))), 1),
                                        effects=list(list(spatial.cov=spde.cov.index$spatial.cov),
                                                     list(Intercept=rep(1,nrow(DFsample)+nrow(DFraster)))
                                        ),
                                        tag="Inference.cov")
            Total.stack.cov <- inla.stack(Inf.stack.cov, Pred.mesh.stack.cov)
            mod.cov <- inla(formula=formula.cov, data=inla.stack.data(Total.stack.cov), family="gaussian",
                            control.predictor = list(compute=FALSE, A=inla.stack.A(Total.stack.cov)))

            List.covariates.inf[[variablesChosenUser[i]]] <- DFsample[,variablesChosenUser[i]]
            List.covariates.pred[[variablesChosenUser[i]]] <- DFraster[,variablesChosenUser[i]]
          }
        } else{
          if((input$IndDataSimulatedLoaded=="load"&input$IndRasterSPDE=="solvecov")|input$IndDataSimulatedLoaded=="sim"){
            cov.inf.indx <- as.vector(unlist(lapply(X=1:nrow(DFsample), FUN=function(Y){which.min(apply(X=(DFsample[!is.na(DFsample[,variablesChosenUser[i]]),1:2]-matrix(DFsample[Y,1:2],ncol=2))**2, MARGIN=1, FUN=sum))})))
            List.covariates.inf[[variablesChosenUser[i]]] <- DFsample[cov.inf.indx, variablesChosenUser[i]]

            IndFactorModelDF <- data.frame(y=1, Ind=c(List.covariates.inf[[variablesChosenUser[i]]]))
            colnames(IndFactorModelDF) <- c("y", variablesChosenUser[i])
            idx.factor <- which(variablesChosenUser[i]==names(DFsample[!as.vector(unlist(lapply(X=DFsample, FUN=is.numeric)))])[names(DFsample[!as.vector(unlist(lapply(X=DFsample, FUN=is.numeric)))])%in%c(variablesChosenUser)])
            
            if(eval(parse(text=paste0("input$IndKindPredictionFactorLevel",idx.factor)))=="nearest"){
              cov.pred.ind <- as.vector(unlist(lapply(X=1:nrow(xy.pred), FUN=function(Y){which.min(apply(X=(DFsample[!is.na(DFsample[,variablesChosenUser[i]]),1:2]-matrix(xy.pred[Y,1:2],ncol=2))**2, MARGIN=1, FUN=sum))})))
              List.covariates.pred[[variablesChosenUser[i]]] <- DFsample[cov.pred.ind, variablesChosenUser[i]]
            } else{
              List.covariates.pred[[variablesChosenUser[i]]] <- rep(NA, times=nrow(xy.pred))
            }
            
          } else {
            DFsampleraster <- rbind(DFsample[,c(1:2,which(variablesChosenUser[i]==colnames(DFsample)))], DFraster[,c(1:2, which(variablesChosenUser[i]==colnames(DFraster)))])
            cov.inf.indx <- as.vector(unlist(lapply(X=1:nrow(DFsample), FUN=function(Y){which.min(apply(X=(DFsampleraster[!is.na(DFsampleraster[,variablesChosenUser[i]]),1:2]-matrix(DFsample[Y,1:2],ncol=2))**2, MARGIN=1, FUN=sum))})))
            List.covariates.inf[[variablesChosenUser[i]]] <- DFsample[cov.inf.indx, variablesChosenUser[i]]
            
            IndFactorModelDF <- data.frame(y=1, Ind=c(List.covariates.inf[[variablesChosenUser[i]]]))
            colnames(IndFactorModelDF) <- c("y", variablesChosenUser[i])
            idx.factor <- which(variablesChosenUser[i]==names(DFsample[!as.vector(unlist(lapply(X=DFsample, FUN=is.numeric)))])[names(DFsample[!as.vector(unlist(lapply(X=DFsample, FUN=is.numeric)))])%in%c(variablesChosenUser)])
            
            if(eval(parse(text=paste0("input$IndKindPredictionFactorLevel",idx.factor)))=="nearest"){
              cov.pred.ind <- as.vector(unlist(lapply(X=1:nrow(xy.pred), FUN=function(Y){which.min(apply(X=(DFsampleraster[!is.na(DFsampleraster[,variablesChosenUser[i]]),1:2]-matrix(xy.pred[Y,1:2],ncol=2))**2, MARGIN=1, FUN=sum))})))
              List.covariates.pred[[variablesChosenUser[i]]] <- DFsample[cov.pred.ind, variablesChosenUser[i]]
            } else{
              List.covariates.pred[[variablesChosenUser[i]]] <- rep(NA, times=nrow(xy.pred))
            }
            
          }
        }
      }
      
      ### Building the main model and stacks structure ====
      
      Inf.geo.effects.list <- list(
        list(),
        list()
      )
      
      Pred.geo.effects.list <- list(
        list(),
        list()
      )
      
      formula_mod <- c("y ~ -1")
      
      A_Inf.spde1 <- list()
      A_Pred.spde1 <- list()
      
      for(i in seq_along(variablesChosen)){
        if(variablesChosen[i]=="Intercept"){
          formula_mod <- paste(formula_mod, "f(Intercept, model='linear')", sep=" + ")
          Inf.geo.effects.list[[2]][["Intercept"]] <- c(rep(1, times=n))
          Pred.geo.effects.list[[2]][["Intercept"]] <- rep(1, times=nrow(A.geo.pred))
          
        } else if(variablesChosen[i]=="Spatial Effect"){
          formula_mod <- paste(formula_mod, "f(Spatial, model=spde)", sep=" + ")
          Inf.geo.effects.list[[1]][["Spatial"]] <- spde.index$Spatial
          Pred.geo.effects.list[[1]][["Spatial"]] <- spde.index$Spatial
          
        } else{
          j <- which(variablesChosen[i]==variablesChosenUser)
          if(!is.character(DFsample[,variablesChosenUser[j]])){
            if(eval(parse(text=paste0("input$IndEffectCov",j)))=="rw1"|eval(parse(text=paste0("input$IndEffectCov",j)))=="rw2"){
              if(eval(parse(text=paste0("input$IndEffectCustomPrior",j)))=="custom"){
                if(eval(parse(text=paste0("input$IndEffectCovKindPrior",j)))=="pc"){
                  assign(paste0("pc.values", variablesChosenUser[j]), c( as.numeric(unlist(strsplit(eval(parse(text=paste0("input$IndEffectCovPriorPCValues",j))), ","))) ))
                  formula_mod <- paste(formula_mod, paste0("f(",variablesChosenUser[j],", model='",eval(parse(text=paste0("input$IndEffectCov",j))),"', hyper=list(prec=list(prior='pc.prec', param=", eval(parse(text=paste0("pc.values", variablesChosenUser[j]))), ")), constr=",eval(parse(text=paste0("input$IndEffectCovConstr",j))),", scale.model=TRUE)"), sep=" + ")
                } else if(eval(parse(text=paste0("input$IndEffectCovKindPrior",j)))=="unif"){
                  lim <- as.numeric(unlist(strsplit(eval(parse(text=paste0("input$IndEffectCovPriorUnif",j))), ",")))
                  sigma <- seq(0, lim[2]*3, length.out=1E5)
                  theta <- -2*log(sigma)
                  logdens <- sapply(X=sigma, FUN=logdunif, lim=lim)
                  unif.prior <- list(theta=list(prior=paste0("table: ", paste(c(theta, logdens), collapse=" "))))
                  assign(paste0("hyper", variablesChosenUser[j]), unif.prior )
                  formula_mod <- paste(formula_mod, paste0("f(",variablesChosenUser[j],", model='",eval(parse(text=paste0("input$IndEffectCov",j))),"', hyper=",paste0("hyper",variablesChosenUser[j]),  ", constr=",eval(parse(text=paste0("input$IndEffectCovConstr",j))),", scale.model=TRUE)"), sep=" + ")
                } else if(eval(parse(text=paste0("input$IndEffectCovKindPrior",j)))=="unifflat"){
                  assign(paste0("unifflat.prior",variablesChosenUser[j]), "expression:
                  log_dens = 0 - log(2) - theta/2
                  ")
                  formula_mod <- paste(formula_mod, paste0("f(",variablesChosenUser[j],", model='",eval(parse(text=paste0("input$IndEffectCov",j))),"', hyper=list(prec=list(prior=",paste0("unifflat.prior",variablesChosenUser[j]),")), constr=",eval(parse(text=paste0("input$IndEffectCovConstr",j))),", scale.model=TRUE)"), sep=" + ")
                } else{
                  assign(paste0("hyper",variablesChosenUser[j]), list(prec=list(prior="loggamma",param=c( as.numeric(unlist(strsplit(eval(parse(text=paste0("input$IndEffectCovPriorBaseValues",j))), ","))) ))) )
                  formula_mod <- paste(formula_mod, paste0("f(",variablesChosenUser[j],", model='",eval(parse(text=paste0("input$IndEffectCov",j))),"', hyper=",paste0("hyper",variablesChosenUser[j]),  ", constr=",eval(parse(text=paste0("input$IndEffectCovConstr",j))),", scale.model=TRUE)"), sep=" + ")
                }
              } else{
                formula_mod <- paste(formula_mod, paste0("f(",variablesChosenUser[j],", model='",eval(parse(text=paste0("input$IndEffectCov",j))),  "', constr=",eval(parse(text=paste0("input$IndEffectCovConstr",j))),", scale.model=TRUE)"), sep=" + ")
              }
              group.cov <- inla.group(x=c(List.covariates.inf[[variablesChosenUser[j]]], List.covariates.pred[[variablesChosenUser[j]]]), n=eval(parse(text=paste0("input$IndEffectCovNodes",j))), method="cut")
              
              Inf.geo.effects.list[[2]][[variablesChosenUser[j]]] <- group.cov[seq_len(n)]
              Pred.geo.effects.list[[2]][[variablesChosenUser[j]]] <- group.cov[-seq_len(n)]
              
            } else if(eval(parse(text=paste0("input$IndEffectCov",j)))=="spde1"){
              Tot_cov <- c(List.covariates.inf[[variablesChosenUser[j]]], List.covariates.pred[[variablesChosenUser[j]]])
              spde1_nodes <- seq(min(Tot_cov), max(Tot_cov), length.out=eval(parse(text=paste0("input$IndEffectCovNodes",j))))
              mesh1d <- fm_mesh_1d(loc=spde1_nodes)
              
              if(eval(parse(text=paste0("input$IndEffectCustomPrior",j)))=="custom"){
                if(eval(parse(text=paste0("input$IndEffectCovKindPrior",j)))=="pc"){
                  hyper <- as.numeric(unlist(strsplit(eval(parse(text=paste0("input$IndEffectCovPriorPCValues",i))), ",")))
                  assign(paste0("spde1d_",variablesChosenUser[j]) ,inla.spde2.pcmatern(mesh=mesh1d, prior.range=c(abs(diff(range(Tot_cov)))/5, 0.5), prior.sigma=c(1,0.5)))
                } else{
                  prior.range <- as.numeric(unlist(strsplit(eval(parse(text=paste0("input$IndEffectCovPriorBaseValues",i))), ",")))[1:3]
                  prior.sigma <- as.numeric(unlist(strsplit(eval(parse(text=paste0("input$IndEffectCovPriorBaseValues",i))), ",")))[4:6]
                  alpha <- 2; d <- 2
                  nu <-  alpha - d/2
                  kappa0 <-  log(8*nu)/2 -log(prior.range[1])
                  tau0 <-  0.5*(lgamma(nu) - lgamma(nu + d/2) - d/2*log(4*pi)) - nu*kappa0 - log(prior.sigma[1])
                  assign(paste0("spde1d_",variablesChosenUser[j]), inla.spde2.matern(mesh = mesh1d, B.tau = cbind(tau0, nu, -1), B.kappa = cbind(kappa0, -1, 0),
                                                                                     theta.prior.mean = c(prior.range[2],prior.sigma[2]), theta.prior.prec = c(prior.range[3],prior.sigma[3])))
                }
              } else {
                assign(paste0("spde1d_",variablesChosenUser[j]), inla.spde2.pcmatern(mesh=mesh1d, prior.range=c(abs(diff(range(Tot_cov)))/5, 0.5), prior.sigma=c(1,0.5)))
              }
              
              spde1d.index <- inla.spde.make.index(name=paste0(variablesChosenUser[j]), n.spde=eval(parse(text=paste0("spde1d_",variablesChosenUser[j])))$n.spde)
              formula_mod <- paste(formula_mod, paste0("f(",variablesChosenUser[j],", model=", paste0("spde1d_",variablesChosenUser[j]),  ")"), sep=" + ")

              Inf.geo.effects.list[[length(A_Inf.spde1)+3]] <- list()
              Pred.geo.effects.list[[length(A_Inf.spde1)+3]] <- list()
              Inf.geo.effects.list[[length(A_Inf.spde1)+3]][[variablesChosenUser[j]]] <- spde1d.index[[1]]
              Pred.geo.effects.list[[length(A_Inf.spde1)+3]][[variablesChosenUser[j]]] <- spde1d.index[[1]]
              
              A_Inf.spde1[[length(A_Inf.spde1)+1]] <- inla.spde.make.A(mesh=mesh1d, loc=Tot_cov[seq_len(n)])
              A_Pred.spde1[[length(A_Pred.spde1)+1]] <- inla.spde.make.A(mesh=mesh1d, loc=Tot_cov[-seq_len(n)])
              
            } else if(eval(parse(text=paste0("input$IndEffectCov",j)))=="linear"){
              if(eval(parse(text=paste0("input$IndEffectCustomPrior",j)))=="custom"){
                hyper <- as.numeric(unlist(strsplit(eval(parse(text=paste0("input$IndEffectCovPriorBaseValues",j))), ",")))
                formula_mod <- paste(formula_mod, paste0("f(",variablesChosenUser[j],", model='",eval(parse(text=paste0("input$IndEffectCov",j))),"', mean.linear=", hyper[1], ",prec.linear=", hyper[2],")"), sep=" + ")
              } else{
                formula_mod <- paste(formula_mod, paste0("f(",variablesChosenUser[j],", model='",eval(parse(text=paste0("input$IndEffectCov",j))),"')"), sep=" + ")
              }
              
              Inf.geo.effects.list[[2]][[variablesChosenUser[j]]] <- c(List.covariates.inf[[variablesChosenUser[j]]])
              Pred.geo.effects.list[[2]][[variablesChosenUser[j]]] <- c(List.covariates.pred[[variablesChosenUser[j]]])
            } else{
              showNotification(ui=paste("The effect of numerical covariates cannot possess an independent and identically distributed (iid) structure. If this is required, the variable values should be recorded as text, not as numerical input.."), duration = NULL)
            }
          } else{ # Section for factor variables
            if(eval(parse(text=paste0("input$IndEffectCov",j)))=="iid"){
              IndFactorModelDF <- data.frame(y=1, Ind=c(List.covariates.inf[[variablesChosenUser[j]]]))
              colnames(IndFactorModelDF) <- c("y", variablesChosenUser[j])
              idx.factor <- which(variablesChosenUser[j]==names(DFsample[!as.vector(unlist(lapply(X=DFsample, FUN=is.numeric)))])[names(DFsample[!as.vector(unlist(lapply(X=DFsample, FUN=is.numeric)))])%in%c(variablesChosenUser)])
              
              Inf.geo.effects.list[[2]][[variablesChosenUser[j]]] <- c(List.covariates.inf[[variablesChosenUser[j]]])
              
              Pred.geo.effects.list[[2]][[variablesChosenUser[j]]] <- if(eval(parse(text=paste0("input$IndKindPredictionFactorLevel",idx.factor)))=="reference"){
                rep( eval(parse(text=paste0("input$IndEffectCovFactorPred",idx.factor))), times=length(List.covariates.pred[[variablesChosenUser[j]]]))
              } else{List.covariates.pred[[variablesChosenUser[j]]]}
              if(eval(parse(text=paste0("input$IndEffectCustomPrior",j)))=="custom"){
                if(eval(parse(text=paste0("input$IndEffectCovKindPrior",j)))=="pc"){
                  hyper <- list(prec=list(prior="pc.prior",param=c( as.numeric(unlist(strsplit(eval(parse(text=paste0("input$IndEffectCovPriorPCValues",j))), ","))) )))
                  formula_mod <- paste(formula_mod, paste0("f(",variablesChosenUser[j],", model='",eval(parse(text=paste0("input$IndEffectCov",j))),"', hyper=",hyper, ", constr=",eval(parse(text=paste0("input$IndEffectCovConstr",j))),")"), sep=" + ")
                } else if(eval(parse(text=paste0("input$IndEffectCovKindPrior",j)))=="base"){
                  hyper <- list(prec=list(prior="loggamma",param=c( as.numeric(unlist(strsplit(eval(parse(text=paste0("input$IndEffectCovPriorBaseValues",j))), ","))) )))
                  formula_mod <- paste(formula_mod, paste0("f(",variablesChosenUser[j],", model='iid', hyper=",hyper,  ", constr=",eval(parse(text=paste0("input$IndEffectCovConstr",j))),")"), sep=" + ")
                } else if(eval(parse(text=paste0("input$IndEffectCovKindPrior",j)))=="unif"){
                  lim <- as.numeric(unlist(strsplit(eval(parse(text=paste0("input$IndEffectCovPriorUnif",j))), ",")))
                  sigma <- seq(lim[1]+1E-5, lim[2]*3, length.out=1E5)
                  theta <- -2*log(sigma)
                  logdens <- sapply(X=sigma, FUN=logdunif, lim=lim)
                  unif.prior <- list(theta=list(prior=paste0("table: ", paste(c(theta, logdens), collapse=" "))))
                  assign(paste0("hyper", variablesChosenUser[j]), unif.prior )
                  formula_mod <- paste(formula_mod, paste0("f(",variablesChosenUser[j],", model='iid', hyper=",paste0("hyper",variablesChosenUser[j]),  ", constr=",eval(parse(text=paste0("input$IndEffectCovConstr",j))),", scale.model=TRUE)"), sep=" + ")
                } else{ # flatunif
                  unifflat.prior= "expression:
                  log_dens = 0 - log(2) - theta/2
                  "
                  formula_mod <- paste(formula_mod, paste0("f(",variablesChosenUser[j],", model='iid', hyper=list(prec=list(prior=",unifflat.prior,")), constr=",eval(parse(text=paste0("input$IndEffectCovConstr",j))),")"), sep=" + ")
                }
              } else{
                formula_mod <- paste(formula_mod, paste0("f(", variablesChosenUser[j], ", model='iid'",  ", constr=",eval(parse(text=paste0("input$IndEffectCovConstr",j))),")"), sep=" + ")
              }
            } else{
              IndFactorModelDF <- data.frame(y=1, Ind=c(List.covariates.inf[[variablesChosenUser[j]]]))
              colnames(IndFactorModelDF) <- c("y", variablesChosenUser[j])
              FactorVariables <- data.frame( model.matrix(object=as.formula(paste0("y~-1+", variablesChosenUser[j])), IndFactorModelDF ) )
              idx.factor <- which(variablesChosenUser[j]==names(DFsample[!as.vector(unlist(lapply(X=DFsample, FUN=is.numeric)))])[names(DFsample[!as.vector(unlist(lapply(X=DFsample, FUN=is.numeric)))])%in%c(variablesChosenUser)])
              for(l in setdiff(colnames(FactorVariables),paste0(variablesChosenUser[j], eval(parse(text=paste0("input$IndEffectCovFactorPred",idx.factor))) ))){
                ll <- l
                Inf.geo.effects.list[[2]][[l]] <- FactorVariables[,l]
                Pred.geo.effects.list[[2]][[l]] <- rep(0, length(List.covariates.pred[[variablesChosenUser[j]]]))
                if(eval(parse(text=paste0("input$IndEffectCustomPrior",j)))=="custom"|eval(parse(text=paste0("input$IndEffectCov",j)))=="linear"){
                  hyper <- as.numeric(unlist(strsplit(eval(parse(text=paste0("input$IndEffectCovPriorBaseValues",j))), ",")))
                  formula_mod <- paste(formula_mod, paste0("f(", l, ", model='linear', mean.linear=", hyper[1], ",prec.linear=", hyper[2],")"), sep=" + ")
                } else{
                  formula_mod <- paste(formula_mod, paste0("f(", l, ", model='linear')"), sep=" + ")
                }
              }
            }
          }
        }
      }
      
      ### Stacks of the geostatistical and prediction layers ====
      
      ResponseVariable <- DFsample[,3]
      
      A_inf_tot <- c(A.inf,1)
      if(length(A_Inf.spde1)>0){
        for(i in seq_along(A_Inf.spde1)){
          A_inf_tot[[2+i]] <- A_Inf.spde1[[i]]
        }
      }
      
      A_pred_tot <- c(A.geo.pred,1)
      if(length(A_Inf.spde1)>0){
        for(i in seq_along(A_Inf.spde1)){
          A_pred_tot[[2+i]] <- A_Pred.spde1[[i]]
        }
      }
      
      Inf.geo.stack <- inla.stack(data=list(y=ResponseVariable),
                                   A=A_inf_tot,
                                   effects=Inf.geo.effects.list,
                                   tag="Inference_geo")
      
      Pred.geo.stack <- inla.stack(data=list(y=matrix(NA, nrow=nrow(A.geo.pred), ncol=1)),
                                    A=A_pred_tot,
                                    effects=Pred.geo.effects.list,
                                    tag="Prediction_geo")
      
      Total.stack <- inla.stack(Inf.geo.stack, Pred.geo.stack)
      status.stacks <- "ok"
      
      ### INLA model ====
      
      formula_inla <- as.formula(formula_mod)
      
      if(input$autocustomIndFamily=='custom'){
        if(input$IndFamilyPriorKind=="pc"){
          family.pcprec <- as.numeric(unlist(strsplit(input$IndFamilyHyper,",")))
          controlFamily <- list(list(hyper = list(prec = list(prior="pc.prec", param=family.pcprec))), list())
        } else if(input$IndFamilyPriorKind=="unif"){
          lim <- as.numeric(unlist(strsplit(input$IndFamilyHyper,",")))
          sigma <- seq(0, lim[2]*3, length.out=1E5)
          theta <- -2*log(sigma)
          logdens <- sapply(X=sigma, FUN=logdunif, lim=lim)
          family.unif.prior <- list(theta=list(prior=paste0("table: ", paste(c(theta, logdens), collapse=" "))))
          controlFamily <- list(list(hyper = list(theta = list(prior=family.unif.prior))), list())
        } else if(input$IndFamilyPriorKind=="unifflat"){
          unifflat.prior= "expression:
                  log_dens = 0 - log(2) - theta/2
                  "
          controlFamily <- list(hyper = list(prec = list(prior=unifflat.prior)))
        } else{
          controlFamily <- list(hyper = list(prec = list(prior="loggamma", param=as.numeric(unlist(strsplit(input$IndFamilyHyper,","))))))
        }
      } else{
        controlFamily <- list()
      }
      
      if(input$autocustomIndMode=='custom'){
        controlModeTheta <- list(theta=as.numeric(unlist(strsplit(input$IndModeHyper,","))), restart=TRUE)
      } else{controlModeTheta <- inla.set.control.mode.default()}
      
      if(input$INLAModeInd=="classic"){
        controlINLA <- list(strategy=input$strategyapproxINLAInd,
                            int.strategy=input$strategyintINLAInd)
      } else{
        controlINLA <- list()
      }
      
      Ind.model <- inla(formula=formula_inla, family = input$SelectIndFamily,
                          data = inla.stack.data(Total.stack),
                          control.inla = controlINLA,
                          control.predictor = list(A = inla.stack.A(Total.stack), compute = TRUE, link = 1),
                          control.family = controlFamily,
                          control.mode = controlModeTheta,
                          control.compute = list(cpo = TRUE, dic = TRUE, waic = TRUE, config = TRUE),
                          inla.mode=input$INLAModeInd,
                          verbose=FALSE)
      
      index.pred <- inla.stack.index(Total.stack, "Prediction_geo")$data
      DFpred <- data.frame(Latitude=xy.pred[,1], Longitude=xy.pred[,2])
      DFpred$Abundance.mean <- Ind.model$summary.fitted.values[index.pred, "mean"]
      DFpred$Abundance.median <- Ind.model$summary.fitted.values[index.pred, "0.5quant"]
      DFpred$Abundance.sd <- Ind.model$summary.fitted.values[index.pred, "sd"]
      
      colnames(DFpred)[1:2] <- colnames(DFsample)[1:2]
      
      DFpredictorMeanMedianStdev <- data.frame(Latitude=xy.pred[,1], Longitude=xy.pred[,2],
                                               Predictor.mean=Ind.model$summary.linear.predictor[index.pred, "mean"],
                                               Predictor.median=Ind.model$summary.linear.predictor[index.pred, "0.5quant"],
                                               Predictor.stdev=Ind.model$summary.linear.predictor[index.pred, "sd"])
      
      colnames(DFpredictorMeanMedianStdev)[1:2] <- colnames(DFsample)[1:2]
      
      gridSpatial <- expand.grid(x=seq(min(DFpred[,1]), max(DFpred[,1]),length.out=input$dimIndmap),
                                 y=seq(min(DFpred[,2]), max(DFpred[,2]), length.out=input$dimIndmap))
      gridSpatial <- gridSpatial[which(!is.na(over(SpatialPoints(coords=gridSpatial),SpatialPolygons(Srl=list(Polygons(srl=list(Polygon(coords=mesh$loc[mesh$segm$int$idx[,1], 1:2])), ID="int")))))),1:2]
      
      A.spatial <- inla.spde.make.A(mesh=mesh, loc=as.matrix(gridSpatial))
      
      DFspatialMeanMedianStdev <- data.frame(Latitude=as.vector(gridSpatial[,1]), Longitude=as.vector(gridSpatial[,2]),
                                              Spatial.mean=as.vector(A.spatial%*%Ind.model$summary.random$Spatial$mean),
                                              Spatial.median=as.vector(A.spatial%*%Ind.model$summary.random$Spatial$`0.5quant`),
                                              Spatial.stdev=as.vector(A.spatial%*%Ind.model$summary.random$Spatial$sd))
      
      colnames(DFspatialMeanMedianStdev)[1:2] <- colnames(DFsample)[1:2]
      
      result <- list(DFpredAbunMeanMedianStdev=list(DFpredAbunMeanMedianStdev=DFpred),
                     DFpredictorMeanMedianStdev=list(DFpredictorMeanMedianStdev=DFpredictorMeanMedianStdev),
                     DFspatialMeanMedianStdev=list(DFspatialMeanMedianStdev=DFspatialMeanMedianStdev),
                     IndModel=list(IndModel=Ind.model),
                     DFPostFixed=list(DFPostFixed=Ind.model$marginals.fixed),
                     DFPostHyperpar=list(DFPostHyperpar=Ind.model$marginals.hyperpar),
                     Summary.fixed=list(Summary.fixed=Ind.model$summary.fixed),
                     Summary.hyperpar=list(Summary.hyperpar=Ind.model$summary.hyperpar),
                     SummaryInternalHyper=list(SummaryInternalHyper=Ind.model$internal.summary.hyperpar),
                     SummaryCPO=list(SummaryCPO=na.omit(Ind.model$cpo$cpo)),
                     DICmodel=list(DICmodel=data.frame(DIC=Ind.model$dic$family.dic, row.names="Geostatistical")),
                     WAICmodel=list(WAICmodel=data.frame(WAIC=cbind(unlist(lapply(na.omit(unique(Ind.model$dic$family)), function(i){sum(Ind.model$waic$local.waic[which(Ind.model$dic$family==i)])}))), row.names="Geostatistical"))
                     )
      
      t2 <- Sys.time()
      difftime(t2,t1, units="secs")
      showNotification(ui=paste("The model has been fitted:", as.numeric(round(Ind.model$cpu.used[4])),
                                "(Geostatistical model) and", as.numeric(round(difftime(t2,t1, units="secs"))),
                                "(overall process) secs." ), duration = NULL)
      # showNotification(ui=paste("The model's DIC is", Inderential.model$dic$dic), duration = NULL)
      return(result)
    })
    
    
    dataggplotIndAbundanceMeanMedianStdevFit <- function(){
      DF <- IndModelFit()$DFpredAbunMeanMedianStdev$DFpredAbunMeanMedianStdev
      return(DF)
    }
    
    DFIndAbundance <- reactive({
      DF <- IndModelFit()$DFpredAbunMeanMedianStdev$DFpredAbunMeanMedianStdev
      condition <- !(length(unique(DF[,1]))*length(unique(DF[,2]))==nrow(DF))&val()
      if(condition){
        grid <- expand.grid(seq(min(DF[,1]),max(DF[,1]), length.out=200),seq(min(DF[,2]),max(DF[,2]), length.out=200))
        DFInter <- data.frame(Latitude=grid[,1],Longitude=grid[,2])
        DFInter$Abundance.mean <- InterpolateIrrGrid(z=DF$Abundance.mean,loc=DF[,1:2], gridInter=grid)$DataInter$z
        DFInter$Abundance.median <- InterpolateIrrGrid(z=DF$Abundance.median,loc=DF[,1:2], gridInter=grid)$DataInter$z
        DFInter$Abundance.sd <- InterpolateIrrGrid(z=DF$Abundance.sd,loc=DF[,1:2], gridInter=grid)$DataInter$z
        colnames(DFInter)[1:2] <- colnames(DF)[1:2]
        DF <- DFInter
      }
      
      return(DF)
    })
    
    ggplotIndAbundanceMeanMedianStdevFit <- function(){
      DF <- DFIndAbundance()
      g1 <- ggplot(DF) + geom_tile(aes(x=DF[,1], y=DF[,2], fill=Abundance.mean)) +
        scale_fill_viridis_c(option = "turbo") + theme_bw() + coord_fixed() + labs(fill="Values") +
        xlab(colnames(DF)[1]) + ylab(colnames(DF)[1]) + ggtitle("Response predicted (mean)") +
        theme(plot.title=element_text(color = "black", size = 16, face = "bold", hjust = 0.5))
      
      g2 <- ggplot(DF) + geom_tile(aes(x=DF[,1], y=DF[,2], fill=Abundance.median)) +
        scale_fill_viridis_c(option = "turbo") + theme_bw() + coord_fixed() + labs(fill="Values") +
        xlab(colnames(DF)[1]) + ylab(colnames(DF)[1]) + ggtitle("Response predicted (median)") +
        theme(plot.title=element_text(color = "black", size = 16, face = "bold", hjust = 0.5))
      
      g3 <- ggplot(DF) + geom_tile(aes(x=DF[,1], y=DF[,2], fill=Abundance.sd))+
        scale_fill_viridis_c(option = "turbo") + theme_bw() + coord_fixed() + labs(fill="Values") +
        xlab(colnames(DF)[1]) + ylab(colnames(DF)[1]) + ggtitle("Response predicted (stdev.)") +
        theme(plot.title=element_text(color = "black", size = 16, face = "bold", hjust = 0.5))
      
      gt <- g1+g2+g3
      
      return(gt)
    }
    
    downloadFile(
      id = "ggplotIndAbundanceMeanMedianStdevFit",
      logger = ss_userAction.Log,
      filenameroot = "ggplotIndAbundanceMeanMedianStdevFit",
      aspectratio  = 1,
      downloadfxns = list(png  = ggplotIndAbundanceMeanMedianStdevFit,
                          csv = dataggplotIndAbundanceMeanMedianStdevFit,
                          txt = dataggplotIndAbundanceMeanMedianStdevFit)
    )
    
    downloadablePlot(id = "ggplotIndAbundanceMeanMedianStdevFit",
                     logger = ss_userAction.Log,
                     filenameroot = "ggplotIndAbundanceMeanMedianStdevFit",
                     aspectratio  = 1,
                     downloadfxns = list(png = ggplotIndAbundanceMeanMedianStdevFit,
                                         csv = dataggplotIndAbundanceMeanMedianStdevFit,
                                         txt = dataggplotIndAbundanceMeanMedianStdevFit),
                     visibleplot  = ggplotIndAbundanceMeanMedianStdevFit)
    
    
    dataggplotIndSpatialMeanMedianStdev <- function(){
      DF <- IndModelFit()$DFspatialMeanMedianStdev$DFspatialMeanMedianStdev
      return(DF)
    }
    
    ggplotIndSpatialMeanMedianStdev <- function(){
      DF <- IndModelFit()$DFspatialMeanMedianStdev$DFspatialMeanMedianStdev
      g1 <- ggplot(DF) + geom_tile(aes(x=DF[,1], y=DF[,2], fill=Spatial.mean)) +
        scale_fill_viridis_c(option = "turbo") + theme_bw() + coord_fixed() + labs(fill="Values") +
        xlab(colnames(DF)[1]) + ylab(colnames(DF)[1]) + ggtitle("Posterior spatial (mean)") +
        theme(plot.title=element_text(color = "black", size = 16, face = "bold", hjust = 0.5))
      
      g2 <- ggplot(DF) + geom_tile(aes(x=DF[,1], y=DF[,2], fill=Spatial.median)) +
        scale_fill_viridis_c(option = "turbo") + theme_bw() + coord_fixed() + labs(fill="Values") +
        xlab(colnames(DF)[1]) + ylab(colnames(DF)[1]) + ggtitle("Posterior spatial (median)") +
        theme(plot.title=element_text(color = "black", size = 16, face = "bold", hjust = 0.5))
      
      g3 <- ggplot(DF) + geom_tile(aes(x=DF[,1], y=DF[,2], fill=Spatial.stdev))+
        scale_fill_viridis_c(option = "turbo") + theme_bw() + coord_fixed() + labs(fill="Values") +
        xlab(colnames(DF)[1]) + ylab(colnames(DF)[1]) + ggtitle("Posterior spatial (stdev.)") +
        theme(plot.title=element_text(color = "black", size = 16, face = "bold", hjust = 0.5))
      
      gt <- g1+g2+g3
      
      return(gt)
    }
    
    downloadFile(
      id = "ggplotIndSpatialMeanMedianStdev",
      logger = ss_userAction.Log,
      filenameroot = "ggplotIndSpatialMeanMedianStdev",
      aspectratio  = 1,
      downloadfxns = list(png  = ggplotIndSpatialMeanMedianStdev,
                          csv = dataggplotIndSpatialMeanMedianStdev,
                          txt = dataggplotIndSpatialMeanMedianStdev)
    )
    
    downloadablePlot(id = "ggplotIndSpatialMeanMedianStdev",
                     logger = ss_userAction.Log,
                     filenameroot = "ggplotIndSpatialMeanMedianStdev",
                     aspectratio  = 1,
                     downloadfxns = list(png=ggplotIndSpatialMeanMedianStdev,
                                         csv = dataggplotIndSpatialMeanMedianStdev,
                                         txt = dataggplotIndSpatialMeanMedianStdev),
                     visibleplot  = ggplotIndSpatialMeanMedianStdev)
    
    dataggplotIndPredictorMeanMedianStdevFit <- function(){
      DF <- IndModelFit()$DFpredictorMeanMedianStdev$DFpredictorMeanMedianStdev
      return(DF)
    }
    
    DFIndPredictor <- reactive({
      DF <- IndModelFit()$DFpredictorMeanMedianStdev$DFpredictorMeanMedianStdev
      condition <- !(length(unique(DF[,1]))*length(unique(DF[,2]))==nrow(DF))&val()
      if(condition){
        grid <- expand.grid(seq(min(DF[,1]),max(DF[,1]), length.out=200),seq(min(DF[,2]),max(DF[,2]), length.out=200))
        DFInter <- data.frame(Latitude=grid[,1],Longitude=grid[,2])
        DFInter$Predictor.mean <- InterpolateIrrGrid(z=DF$Predictor.mean,loc=DF[,1:2], gridInter=grid)$DataInter$z
        DFInter$Predictor.median <- InterpolateIrrGrid(z=DF$Predictor.median,loc=DF[,1:2], gridInter=grid)$DataInter$z
        DFInter$Predictor.stdev <- InterpolateIrrGrid(z=DF$Predictor.stdev,loc=DF[,1:2], gridInter=grid)$DataInter$z
        colnames(DFInter)[1:2] <- colnames(DF)[1:2]
        DF <- DFInter
      }
      
      return(DF)
    })
    
    ggplotIndPredictorMeanMedianStdevFit <- function(){
      DF <- DFIndPredictor()
      g1 <- ggplot(DF) + geom_tile(aes(x=DF[,1], y=DF[,2], fill=Predictor.mean)) +
        scale_fill_viridis_c(option = "turbo") + theme_bw() + coord_fixed() + labs(fill="Values") +
        xlab(colnames(DF)[1]) + ylab(colnames(DF)[1]) + ggtitle("Linear predictor (mean)") +
        theme(plot.title=element_text(color = "black", size = 16, face = "bold", hjust = 0.5))
      
      g2 <- ggplot(DF) + geom_tile(aes(x=DF[,1], y=DF[,2], fill=Predictor.median)) +
        scale_fill_viridis_c(option = "turbo") + theme_bw() + coord_fixed() + labs(fill="Values") +
        xlab(colnames(DF)[1]) + ylab(colnames(DF)[1]) + ggtitle("Linear predictor (median)") +
        theme(plot.title=element_text(color = "black", size = 16, face = "bold", hjust = 0.5))
      
      g3 <- ggplot(DF) + geom_tile(aes(x=DF[,1], y=DF[,2], fill=Predictor.stdev))+
        scale_fill_viridis_c(option = "turbo") + theme_bw() + coord_fixed() + labs(fill="Values") +
        xlab(colnames(DF)[1]) + ylab(colnames(DF)[1]) + ggtitle("Linear predictor (stdev.)") +
        theme(plot.title=element_text(color = "black", size = 16, face = "bold", hjust = 0.5))
      
      gt <- g1+g2+g3
      
      return(gt)
    }
    
    downloadFile(
      id = "ggplotIndPredictorMeanMedianStdevFit",
      logger = ss_userAction.Log,
      filenameroot = "ggplotIndPredictorMeanMedianStdevFit",
      aspectratio  = 1,
      downloadfxns = list(png  = ggplotIndPredictorMeanMedianStdevFit,
                          csv = dataggplotIndPredictorMeanMedianStdevFit,
                          txt = dataggplotIndPredictorMeanMedianStdevFit)
    )
    
    downloadablePlot(id = "ggplotIndPredictorMeanMedianStdevFit",
                     logger = ss_userAction.Log,
                     filenameroot = "ggplotIndPredictorMeanMedianStdevFit",
                     aspectratio  = 1,
                     downloadfxns = list(png = ggplotIndPredictorMeanMedianStdevFit,
                                         csv = dataggplotIndPredictorMeanMedianStdevFit,
                                         txt = dataggplotIndPredictorMeanMedianStdevFit),
                     visibleplot  = ggplotIndPredictorMeanMedianStdevFit)
    
    dataggplotIndFixParamFit <- function(){
      DF <- as.data.frame(IndModelFit()$DFPostFixed$DFPostFixed)
      return(DF)
    }
    
    ggplotIndFixParamFit <- function(){
      DF <- IndModelFit()$DFPostFixed$DFPostFixed
      gl <- c()
      for(i in 1:length(DF)){
        assign(paste0("g",i),
               ggplot(data=data.frame(x=DF[[i]][,1], y=DF[[i]][,2]), aes(x=x,y=y)) + geom_line() +
                 theme_bw() + xlab(names(DF)[i]) + ylab(HTML(paste("Density f(",names(DF)[i],")"))) +
                 ggtitle(names(DF)[i]) + theme(plot.title=element_text(color = "black", size = 16, face = "bold", hjust = 0.5))
        )
        gl[i] <- paste0("g",i)
      }
      gt <- eval(parse(text=paste(gl, collapse="+")))
      return(gt)
    }
    
    downloadFile(
      id = "ggplotIndFixParamFit",
      logger = ss_userAction.Log,
      filenameroot = "ggplotIndFixParamFit",
      aspectratio  = 1,
      downloadfxns = list(png  = ggplotIndFixParamFit,
                          csv = dataggplotIndFixParamFit,
                          txt = dataggplotIndFixParamFit)
    )
    
    downloadablePlot(id = "ggplotIndFixParamFit",
                     logger = ss_userAction.Log,
                     filenameroot = "ggplotIndFixParamFit",
                     aspectratio  = 1,
                     downloadfxns = list(png = ggplotIndFixParamFit,
                                         csv = dataggplotIndFixParamFit,
                                         txt = dataggplotIndFixParamFit),
                     visibleplot  = ggplotIndFixParamFit)
    
    
    dataggplotIndHyperParamFit <- function(){
      DF <- as.data.frame(IndModelFit()$DFPostHyperpar$DFPostHyperpar)
      return(DF)
    }
    
    ggplotIndHyperParamFit <- function(){
      DF <- IndModelFit()$DFPostHyperpar$DFPostHyperpar
      gl <- c()
      for(i in 1:length(DF)){
        nm <- strsplit(names(DF)[i], " ")[[1]]
        if(nm[1]=="Precision"){
          namesDFold <- names(DF)[i]
          names(DF)[i] <- paste("Stdev.", paste(nm[-1], collapse=" "), collapse=" ")
          assign(paste0("g",i), try(
            ggplot(data=data.frame(x=inla.tmarginal(function(x) 1/x**0.5, DF[[i]])[,1],
                                   y=inla.tmarginal(function(x) 1/x**0.5, DF[[i]])[,2]), aes(x=x,y=y)) + geom_line() +
              theme_bw()+ xlab(names(DF)[i]) +
              ylab(HTML(paste("Density f(",names(DF)[i],")"))) + ggtitle(names(DF)[i]) +
              theme(plot.title=element_text(color = "black", size = 16, face = "bold", hjust = 0.5)), silent=TRUE)
          )
          if(sum(class(eval(parse(text=paste0("g",i))))=="try-error")==1){
            assign(paste0("g",i),
                   ggplot(data=data.frame(x=inla.smarginal(marginal=DF[[i]])[[1]],
                                          y=inla.smarginal(marginal=DF[[i]])[[2]]), aes(x=x,y=y)) + geom_line() +
                     theme_bw()+ xlab(names(DF)[i]) +
                     ylab(HTML(paste("Density f(",namesDFold,")"))) + ggtitle(namesDFold) +
                     theme(plot.title=element_text(color = "black", size = 16, face = "bold", hjust = 0.5))
            )
          }
        } else{
          assign(paste0("g",i),
                 ggplot(data=data.frame(x=DF[[i]][,1], y=DF[[i]][,2]), aes(x=x,y=y)) + geom_line() +
                   theme_bw()+ xlab(names(DF)[i]) + xlim(quantile(DF[[i]][,1], probs = c(0.025,0.975))) +
                   ylab(HTML(paste("Density f(",names(DF)[i],")"))) + ggtitle(names(DF)[i]) +
                   theme(plot.title=element_text(color = "black", size = 16, face = "bold", hjust = 0.5))
          )
        }
        gl[i] <- paste0("g",i)
      }
      gt <- eval(parse(text=paste(gl, collapse="+")))
      return(gt)
    }
    
    downloadFile(
      id = "ggplotIndHyperParamFit",
      logger = ss_userAction.Log,
      filenameroot = "ggplotIndHyperParamFit",
      aspectratio  = 1,
      downloadfxns = list(png  = ggplotIndHyperParamFit,
                          csv = dataggplotIndHyperParamFit,
                          txt = dataggplotIndHyperParamFit)
    )
    
    downloadablePlot(id = "ggplotIndHyperParamFit",
                     logger = ss_userAction.Log,
                     filenameroot = "ggplotIndHyperParamFit",
                     aspectratio  = 1,
                     downloadfxns = list(png=ggplotIndHyperParamFit,
                                         csv = dataggplotIndHyperParamFit,
                                         txt = dataggplotIndHyperParamFit),
                     visibleplot  = ggplotIndHyperParamFit)
    
    
    tableIndModelFixedPar <- function(){
      DF <- IndModelFit()$Summary.fixed$Summary.fixed %>%
        mutate(across(where(is.numeric), round, digits = 2))
    }
    
    downloadableTable("tableIndModelFixedPar",
                      logger=ss_userAction.Log,
                      filenameroot="tableIndModelFixedPar",
                      downloaddatafxns=list(csv=tableIndModelFixedPar,
                                            tsv=tableIndModelFixedPar),
                      tabledata=tableIndModelFixedPar, rownames = TRUE,
                      caption="Summary fixed parameters")
    
    tableIndModelHyperPar <- function(){
      DF <- IndModelFit()$Summary.hyperpar$Summary.hyperpar %>%
        mutate(across(where(is.numeric), round, digits = 2))
    }
    
    downloadableTable("tableIndModelHyperPar",
                      logger=ss_userAction.Log,
                      filenameroot="tableIndModelHyperPar",
                      downloaddatafxns=list(csv=tableIndModelHyperPar,
                                            tsv=tableIndModelHyperPar),
                      tabledata=tableIndModelHyperPar, rownames = TRUE,
                      caption="Summary hyperparameters")
    
    tableIndModelInternalHyperPar <- function(){
      DF <- IndModelFit()$SummaryInternalHyper$SummaryInternalHyper %>%
        mutate(across(where(is.numeric), round, digits = 2))
    }
    
    downloadableTable("tableIndModelInternalHyperPar",
                      logger=ss_userAction.Log,
                      filenameroot="tableIndModelInternalHyperPar",
                      downloaddatafxns=list(csv=tableIndModelInternalHyperPar,
                                            tsv=tableIndModelInternalHyperPar),
                      tabledata=tableIndModelInternalHyperPar, rownames = TRUE,
                      caption="Summary internal hyperparameters")
    
    dataIndDICtable <- function(){
      DF <- IndModelFit()$DICmodel$DICmodel %>%
        mutate(across(where(is.numeric), round, digits = 2))
      return(DF)
    }
    
    downloadableTable("dataIndDICtable",
                      logger=ss_userAction.Log,
                      filenameroot="dataIndDICtable",
                      downloaddatafxns=list(csv=dataIndDICtable,
                                            tsv=dataIndDICtable),
                      tabledata=dataIndDICtable, rownames = TRUE,
                      caption="Model DIC")
    
    dataIndWAICtable <- function(){
      DF <- IndModelFit()$WAICmodel$WAICmodel %>%
        mutate(across(where(is.numeric), round, digits = 2))
      return(DF)
    }
    
    downloadableTable("dataIndWAICtable",
                      logger=ss_userAction.Log,
                      filenameroot="dataIndWAICtable",
                      downloaddatafxns=list(csv=dataIndWAICtable,
                                            tsv=dataIndWAICtable),
                      tabledata=dataIndWAICtable, rownames = TRUE,
                      caption="Model WAIC")
    
    dataIndCPOtable <- function(){
      CPO <- IndModelFit()$SummaryCPO$SummaryCPO
      DF <- data.frame(n=length(CPO), mean=mean(CPO), median=median(CPO), stdev=sd(CPO),
                       quantile2.5=quantile(CPO,probs=c(0.025,0.975))[1],
                       quantile97.5=quantile(CPO,probs=c(0.025,0.975))[2]) %>%
        mutate(across(where(is.numeric), round, digits = 2))
      return(DF)
    }
    
    downloadableTable("dataIndCPOtable",
                      logger=ss_userAction.Log,
                      filenameroot="dataIndCPOtable",
                      downloaddatafxns=list(csv=dataIndCPOtable,
                                            tsv=dataIndCPOtable),
                      tabledata=dataIndCPOtable, rownames = FALSE,
                      caption="Summary CPO")
    
    # Server_LGCPModelling_Section ----
    
    LgcpCheckBoxNames <- function(){
      if(input$LgcpDataSimulatedLoaded=="load"){
        DF <- as.data.frame(datareadSample())
        if(input$SelectLgcpFamily=="binomial"){DFnames <- names(DF)[c(5:ncol(DF))]}
        else{DFnames <- names(DF)[c(4:ncol(DF))]}
      } else if(input$LgcpDataSimulatedLoaded=="sim"){
        DF <- Lgcp.sampling()
        DFnames <- names(DF)[c(4)]
      }
      return(DFnames)
    }
    
    observe({
      output$checkBoxLgcpDataFrame <- renderUI({
        tagList(
          checkboxGroupInput(inputId="UserComponentsLgcp",
                             label="User Defined Components",
                             choices=LgcpCheckBoxNames(),
                             selected=c()
          )
        )
      })
    })
    
    observe({
      if(input$LgcpDataSimulatedLoaded=="sim"){DF <- Lgcp.sampling()}
      else{DF <- as.data.frame(datareadSample())}
      output$checkBoxLgcpSharing <- renderUI({
        box(id="LgcpCovariateSharingTerms", width=12, title="Sharing Components",
            checkboxGroupInput(inputId=paste0("UserComponentsLgcpSharing"),
                               label=paste("User Defined Sharing Components:"),
                               choices=c(input$DefaultComponentsLgcp,input$UserComponentsLgcp),
                               selected=c(input$DefaultComponentsLgcp,input$UserComponentsLgcp)
            )
        )
      })
    })
    
    observe({
      if(input$LgcpDataSimulatedLoaded=="load"){
        output$SelectLgcpEffectCov <- renderUI({
          
          if(length(input$UserComponentsLgcp)>0){
            box(id="LgcpCovariateEffects", width=12, title="Covariate Effects",
                lapply(seq_along(input$UserComponentsLgcp), function(i){
                  list(
                    selectInput(inputId=paste0("LgcpEffectCov",i), label=tags$span(style="color: blue; font-weight: bold;", paste(input$UserComponentsLgcp[i],"Effect")) ,
                                choices = list("Linear (or reference level)" = "linear", "Random Walk 1" = "rw1", "Random Walk 2" = "rw2", "SPDE 1" = "spde1", "IID"="iid"),
                                selected = "linear"),
                    conditionalPanel(condition=paste0("input.LgcpEffectCov",i,"=='rw1'||","input.LgcpEffectCov",i,"=='rw2'||","input.LgcpEffectCov",i,"=='spde1'"),
                                     numericInput(inputId=paste0("LgcpEffectCovNodes",i),
                                                  label="Number of nodes",
                                                  value=10, min=1, step=1
                                     )
                    ),
                    radioGroupButtons(
                      inputId = paste0("LgcpEffectCustomPrior",i),
                      label = "Custom Prior",
                      choices = list("Auto" = "auto", "Custom" = "custom"),
                      status = "primary"
                    ),
                    conditionalPanel(condition=paste0("input.LgcpEffectCustomPrior",i,"=='custom'"),
                                     list(
                                       selectInput(inputId = paste0("LgcpEffectCovKindPrior",i),
                                                   label = "Prior distribution",
                                                   choices = list("Base" = "base", "PC prior" = "pc", "Uniform" = "unif", "Flat Uniform" = "flatunif")),
                                       conditionalPanel(condition=paste0("input.LgcpEffectCovKindPrior",i,"=='base'"),
                                                        textInput(inputId = paste0("LgcpEffectCovPriorBaseValues",i),
                                                                  label = "Base Prior Values",
                                                                  value = "0, 1e-3")),
                                       conditionalPanel(condition=paste0("input.LgcpEffectCovKindPrior",i,"=='pc'"),
                                                        textInput(inputId = paste0("LgcpEffectCovPriorPCValues",i),
                                                                  label = "PC-Prior Values",
                                                                  value = "0, 1e-3")),
                                       conditionalPanel(condition=paste0("input.LgcpEffectCovKindPrior",i,"=='unif'"),
                                                        textInput(inputId = paste0("LgcpEffectCovPriorUnif",i),
                                                                  label = "Uniform Lower and upper values",
                                                                  value = "0, 10",
                                                                  placeholder = "Lower and upper values: 'U(a,b)'"))
                                     )
                    ),
                    conditionalPanel(condition=paste0("input.LgcpEffectCov",i,"=='rw1'||",
                                                      "input.LgcpEffectCov",i,"=='rw2'||",
                                                      "input.LgcpEffectCov",i,"=='iid'"),
                                     radioGroupButtons(
                                       inputId = paste0("LgcpEffectCovConstr",i),
                                       label = "Sum to zero",
                                       choices = list("TRUE" = "TRUE", "FALSE" = "FALSE"),
                                       status = "primary"
                                     ))
                  )
                }))
          } else{}
        }) } else if(input$LgcpDataSimulatedLoaded=="sim"){
          output$SelectLgcpEffectCov <- renderUI({
            if(length(input$UserComponentsLgcp)>0){
              box(id="LgcpCovariateEffects", width=12, title="Covariate Effects",
                  lapply(seq_along(input$UserComponentsLgcp), function(i){
                    list(
                      selectInput(inputId=paste0("LgcpEffectCov",i), label=tags$span(style="color: blue; font-weight: bold;", paste(input$UserComponentsLgcp[i],"Effect")) ,
                                  choices = list("Linear (or reference level)" = "linear", "Random Walk 1" = "rw1", "Random Walk 2" = "rw2", "SPDE 1" = "spde1", "IID"="iid"),
                                  selected = "linear"),
                      conditionalPanel(condition=paste0("input.LgcpEffectCov",i,"=='rw1'||","input.LgcpEffectCov",i,"=='rw2'||","input.LgcpEffectCov",i,"=='spde1'"),
                                       numericInput(inputId=paste0("LgcpEffectCovNodes",i),
                                                    label="Number of nodes",
                                                    value=10, min=1, step=1
                                       )
                      ),
                      radioGroupButtons(
                        inputId = paste0("LgcpEffectCustomPrior",i),
                        label = "Custom Prior",
                        choices = list("Auto" = "auto", "Custom" = "custom"),
                        status = "primary"
                      ),
                      conditionalPanel(condition=paste0("input.LgcpEffectCustomPrior",i,"=='custom'"),
                                       list(
                                         selectInput(inputId = paste0("LgcpEffectCovKindPrior",i),
                                                     label = "Prior distribution",
                                                     choices = list("Base" = "base", "PC prior" = "pc", "Uniform" = "unif", "Flat Uniform" = "unifflat")),
                                         conditionalPanel(condition=paste0("input.LgcpEffectCovKindPrior",i,"=='base'"),
                                                          textInput(inputId = paste0("LgcpEffectCovPriorBaseValues",i),
                                                                    label = "Base Prior Values",
                                                                    value = "0, 1e-3")),
                                         conditionalPanel(condition=paste0("input.LgcpEffectCovKindPrior",i,"=='pc'"),
                                                          textInput(inputId = paste0("LgcpEffectCovPriorPCValues",i),
                                                                    label = "PC-Prior Values",
                                                                    value = "0, 1e-3")),
                                         conditionalPanel(condition=paste0("input.LgcpEffectCovKindPrior",i,"=='unif'"),
                                                          textInput(inputId = paste0("LgcpEffectCovPriorUnif",i),
                                                                    label = "Uniform Lower and upper values",
                                                                    value = "0, 10",
                                                                    placeholder = "Lower and upper values: 'U(a,b)'"))
                                       )
                      ),
                      conditionalPanel(condition=paste0("input.LgcpEffectCov",i,"=='rw1'||",
                                                        "input.LgcpEffectCov",i,"=='rw2'||",
                                                        "input.LgcpEffectCov",i,"=='iid'"),
                                       radioGroupButtons(
                                         inputId = paste0("LgcpEffectCovConstr",i),
                                         label = "Sum to zero",
                                         choices = list("TRUE" = "TRUE", "FALSE" = "FALSE"),
                                         status = "primary"
                                       ))
                    )
                  }))
            } else{}
          }) }
    })
    
    observe({
      if(input$LgcpDataSimulatedLoaded=="load"){
        DF <- as.data.frame(datareadSample())
        if(length(input$UserComponentsLgcp)>0){
          LgcpUserComponent <- input$UserComponentsLgcp
          DF2 <- select(DF, LgcpUserComponent[!as.vector(unlist(lapply(X=select(DF,LgcpUserComponent), FUN=is.numeric)))])
          output$SelectLoadLgcpEffectCovFactorPred <- renderUI({
            if(ncol(DF2)>0){
              box(id="LgcpPredFactorLevel", width=12, title="Prediction Factor Level",
                  lapply(seq_along(names(DF2)), function(i){
                    choices <- unique(DF2[[names(DF2)[i]]])
                    list(
                      radioGroupButtons(inputId = paste0("LgcpKindPredictionFactorLevel",i), label = tags$span(style="color: blue; font-weight: bold;", paste(names(DF2)[i], "(prediction protocol)")),
                                        choices = c("Reference level" = "reference", "Nearest level" = "nearest"), status = "success", justified = TRUE),
                      conditionalPanel(
                        condition=paste0("input.LgcpKindPredictionFactorLevel",i,"=='reference'"),
                        selectInput(inputId=paste0("LgcpEffectCovFactorPred",i), label=paste(names(DF2)[i],"Reference Factor"),
                                    choices = choices, selected = choices[1])
                      )
                    )
                  }))
            }
          })
        } else{
          output$SelectLoadLgcpEffectCovFactorPred <- renderUI({})
        }
      }
    })
    
    JointPriorPreviewLgcp <- eventReactive(input$Lgcppreviewpriordistributions,{
      d <- 2; nu <- 1
      rhoinputPC <- as.numeric(unlist(strsplit(input$Lgcprangepcpriorprev, ",")))
      sigmainputPC <- as.numeric(unlist(strsplit(input$Lgcpsigmapcpriorprev, ",")))
      rhoPC <- seq(rhoinputPC[1], rhoinputPC[2], length.out=rhoinputPC[3])
      sigmaPC <- seq(sigmainputPC[1], sigmainputPC[2], length.out=sigmainputPC[3])
      rhosigmaPC <- expand.grid(rho=rhoPC, sigma=sigmaPC)
      rho0PC <- rhoinputPC[4]; alpha1 <- rhoinputPC[5]
      sigma0PC <- sigmainputPC[4]; alpha2 <- sigmainputPC[5]
      
      lambdarhoPC <- -log(alpha1)*rho0PC**(d/2) #-(rho0/sqrt(8*nu))**(d/2)*log(alpha1)
      lambdasigmaPC <- -log(alpha2)/sigma0PC #-(sqrt(8*nu)/rhosigma$rho)**(-nu)*sqrt(gamma(nu)/(gamma(nu+d/2)*(4*pi)**(d/2)))*log(alpha2)/sigma0
      pirhosigmaPCM <- d/2*lambdarhoPC*rhosigmaPC$rho**(-1-d/2)*exp(-lambdarhoPC*rhosigmaPC$rho**(-d/2))*lambdasigmaPC*exp(-lambdasigmaPC*rhosigmaPC$sigma)
      
      probIntPC <- diff(range(rhoPC))*diff(range(sigmaPC))/length(pirhosigmaPCM)*sum(pirhosigmaPCM)
      
      rhoinputBase <- as.numeric(unlist(strsplit(input$Lgcprangebasepriorprev, ",")))
      sigmainputBase <- as.numeric(unlist(strsplit(input$Lgcpsigmabasepriorprev, ",")))
      rho <- seq(rhoinputBase[1], rhoinputBase[2], length.out=rhoinputBase[3])
      sigma <- seq(sigmainputBase[1], sigmainputBase[2], length.out=sigmainputBase[3])
      rhosigma <- expand.grid(rho=rho, sigma=sigma)
      rho0 <- rhoinputBase[4]; sigma0 <- sigmainputBase[4]
      meanrho <- rhoinputBase[5]; meansigma <- sigmainputBase[5]
      sdrho <- rhoinputBase[6]; sdsigma <- sigmainputBase[6]
      pirho <- (sqrt(2*pi*sdrho**2)*rhosigma$rho)**(-1)*exp(-(log(rhosigma$rho/rho0)**2-2*log(rhosigma$rho/rho0)*meanrho+meanrho**2)/(2*sdrho**2))
      pisigma <- (sqrt(2*pi*sdsigma**2)*rhosigma$sigma)**(-1)*exp(-(log(rhosigma$sigma/sigma0)**2-2*log(rhosigma$sigma/sigma0)*meansigma+meansigma**2)/(2*sdsigma**2))
      pirhosigmaM <- pirho*pisigma
      
      probInt <- sum(pirhosigmaM)*diff(range(rhosigma$rho))*diff(range(rhosigma$sigma))/length(pirhosigmaM)
      
      ListPreview <- list(rhosigmaPC=rhosigmaPC, pirhosigmaPCM=pirhosigmaPCM, probIntPC=probIntPC,
                          rhosigma=rhosigma, pirhosigmaM=pirhosigmaM, probInt=probInt)
      
      return(ListPreview)
    })
    
    LgcpPreviewJointPriorPlot <- function(){
      ggplotPCprior <- ggplot() +
        geom_tile(data=data.frame(rho=JointPriorPreviewLgcp()$rhosigmaPC$rho, sigma=JointPriorPreviewLgcp()$rhosigmaPC$sigma, pi=JointPriorPreviewLgcp()$pirhosigmaPCM),
                  mapping=aes(x=rho, y=sigma, fill=pi)) +
        labs(title=paste("Joint PC-Prior. Cumulative Prob.=", round(JointPriorPreviewLgcp()$probIntPC, digits=2)),
             x=expression(rho), y=expression(sigma)) +
        scale_fill_viridis_c(option="turbo") + theme_bw() + theme(plot.title=element_text(face='bold'))
      
      ggplotENpriorAnalytic <- ggplot() +
        geom_tile(data=data.frame(rho=JointPriorPreviewLgcp()$rhosigma$rho, sigma=JointPriorPreviewLgcp()$rhosigma$sigma, pi=JointPriorPreviewLgcp()$pirhosigmaM),
                  mapping=aes(x=rho, y=sigma, fill=pi)) +
        scale_fill_viridis_c(option="turbo") +
        labs(title=paste("Joint Base Prior. Cumulative Prob.=", round(JointPriorPreviewLgcp()$probInt, digits=2)),
             x=expression(rho), y=expression(sigma)) +
        theme_bw() + theme(plot.title=element_text(face='bold'))
      
      ggplotJointPriorPreview <- ggplotPCprior + ggplotENpriorAnalytic
      return(ggplotJointPriorPreview)
    }
    
    downloadFile(
      id = "LgcpPreviewJointPriorPlot",
      logger = ss_userAction.Log,
      filenameroot = "LgcpPreviewJointPriorPlot",
      aspectratio  = 1,
      downloadfxns = list(png  = LgcpPreviewJointPriorPlot)
    )
    
    downloadablePlot(id = "LgcpPreviewJointPriorPlot",
                     logger = ss_userAction.Log,
                     filenameroot = "LgcpPreviewJointPriorPlot",
                     aspectratio  = 1,
                     downloadfxns = list(png=LgcpPreviewJointPriorPlot),
                     visibleplot  = LgcpPreviewJointPriorPlot)
    
    ## Mesh construction for LgcperentialModel ====
    
    LgcpMeshBase <- reactive({
      if(input$LgcpDataSimulatedLoaded=="sim"){
        DFsample <- as.data.frame(Lgcp.sampling())
        convexhull <- chull(DFsample[,1:2])
        convexhull <- c(convexhull, convexhull[1])
        qloc <- quantile(as.vector(dist(DFsample[sample(1:nrow(DFsample),size=min(c(50,nrow(DFsample)))),1:2])),probs=c(0.03,0.3))
        mesh <- fm_mesh_2d_inla(loc.domain=DFsample[convexhull,1:2], cutoff = qloc[1]/2, offset=c(-0.1, -0.2),
                             max.edge=c(qloc[1], qloc[2]))
        sample <- DFsample
      } else if(input$LgcpDataSimulatedLoaded=="load"){
        DFsample <- datareadSample()
        if(input$LgcpRasterSPDE=="raster"){
          rasterSample <- datareadRaster()[sample(1:nrow(datareadRaster()), min(c(50,nrow(datareadRaster())))),1:2]
          qloc <- quantile(as.vector(dist(rasterSample)),probs=c(0.03,0.3))
          convexhull <- chull(datareadRaster()[,1:2])
          convexhull <- c(convexhull, convexhull[1])
          mesh <- fm_mesh_2d_inla(loc.domain=datareadRaster()[convexhull,1:2], cutoff = qloc[1]/2, offset=c(-0.1, -0.2),
                               max.edge=c(qloc[1], qloc[2]))
          sample <- rasterSample
        } else if(input$LgcpRasterSPDE=="solvecov"){
          convexhull <- chull(rbind(DFsample[,1:2],datareadRaster()[,1:2]))
          convexhull <- c(convexhull, convexhull[1])
          qloc <- quantile(as.vector(dist(DFsample[sample(1:nrow(DFsample),size=min(c(50,nrow(DFsample)))),1:2])),probs=c(0.03,0.3))
          mesh <- fm_mesh_2d_inla(loc.domain=rbind(DFsample[,1:2],datareadRaster()[,1:2])[convexhull,], cutoff = qloc[1]/2, offset=c(-0.1, -0.2),
                               max.edge=c(qloc[1], qloc[2]))
          sample <- DFsample
        }
      }
      result <- list(mesh=mesh, qloc=qloc, Sample=sample)
      return(result)
    })
    
    LgcpMesh <- eventReactive(input$buildLgcpMesh, {
      if(input$LgcpDataSimulatedLoaded=="sim"){
        DFsample <- as.data.frame(Lgcp.sampling())
      } else if(input$LgcpDataSimulatedLoaded=="load"){
        DFsample <- as.data.frame(datareadSample())
        DFraster <- try(datareadRaster(), silent=TRUE)
        if(class(DFraster)!="try-error"){
          DFsample <- rbind(DFsample[,1:2], DFraster[,1:2])
        }
      }
      
      interiorNonConvex <- function(x, condition, convex, resolution, file.read){
        if(condition=="interiorLgcpMeshnonconvex"){
          interior <- inla.nonconvex.hull(points=x, convex=convex, resolution=resolution)
        } else if(condition=="interiorLgcpMeshcustomboundary"){
          ext <- tools::file_ext(file.read$datapath)
          validate(need(ext == c("csv","rds"), "Please upload a csv or rds file"))
          if(ext=="csv"){coords <- read.csv(file.read$datapath)}
          else coords <- readRDS(file.read$datapath)
          innerBorderMesh <- SpatialPolygons(Srl=list(Polygons(srl=list(Polygon(coords=coords)), ID="interiorBoundary")))
          interior <- inla.sp2segment(sp=innerBorderMesh)
        } else{interior <- NULL}
        return(interior)
      }
      
      boundaryNonConvex <- function(x, condition, convex, resolution, file.read){
        if(condition=="LgcpMeshnonconvex"){
          boundary <- inla.nonconvex.hull(points=x, convex=convex, resolution=resolution)
        } else if(condition=="LgcpMeshcustomboundary"){
          ext <- tools::file_ext(file.read$datapath)
          validate(need(ext == c("csv","rds"), "Please upload a csv or rds file"))
          if(ext=="csv"){coords <- read.csv(file.read$datapath)}
          else coords <- readRDS(file.read$datapath)
          outerBorderMesh <- SpatialPolygons(Srl=list(Polygons(srl=list(Polygon(coords=coords)), ID="externalBoundary")))
          boundary <- inla.sp2segment(sp=outerBorderMesh)
        } else{boundary <- NULL}
        return(boundary)
      }
      
      if(input$selectionLgcpMesh=="qlocation"){
        if(input$LgcpRasterSPDE=="raster"){
          rasterSample <- datareadRaster()[sample(1:nrow(datareadRaster()), min(c(50,nrow(datareadRaster())))),1:2]
          qloc <- quantile(as.vector(dist(rasterSample)),probs=c(ifelse(input$LgcpCustomMesh,input$LgcpMeshQloc,0.03),0.3))
          convexhull <- chull(DFsample[,1:2])
          convexhull <- c(convexhull, convexhull[1])
          mesh <- fm_mesh_2d_inla(loc.domain=DFsample[convexhull,1:2], cutoff = qloc[1]/2, offset=c(-0.1, -0.2), max.edge=c(qloc[1], qloc[2]),
                               boundary=list(interiorNonConvex(x=as.matrix(rasterSample[,1:2]), condition=input$interiorLgcpMesh, convex=input$interiorcurvatureLgcpMesh, resolution=input$interiorresolutionLgcpMesh, file.read=input$interiorshapefileLgcpMesh),
                                             boundaryNonConvex(x=as.matrix(rasterSample[,1:2]), condition=input$boundaryLgcpMesh, convex=input$curvatureLgcpMesh, resolution=input$resolutionLgcpMesh, file.read=input$shapefileLgcpMesh)))
          sample <- rasterSample
        } else if(input$LgcpRasterSPDE=="solvecov"){
          qloc <- quantile(as.vector(dist(DFsample[sample(1:nrow(DFsample),size=min(c(50,nrow(DFsample)))),1:2])),probs=c(ifelse(input$LgcpCustomMesh,input$LgcpMeshQloc,0.03),0.3))
          convexhull <- chull(DFsample[,1:2])
          convexhull <- c(convexhull, convexhull[1])
          mesh <- fm_mesh_2d_inla(loc.domain=DFsample[convexhull,1:2], cutoff = qloc[1]/2, offset=c(-0.1, -0.2), max.edge=c(qloc[1], qloc[2]),
                               boundary=list(interiorNonConvex(x=as.matrix(DFsample[,1:2]), condition=input$interiorLgcpMesh, convex=input$interiorcurvatureLgcpMesh, resolution=input$interiorresolutionLgcpMesh, file.read=input$interiorshapefileLgcpMesh),
                                             boundaryNonConvex(x=as.matrix(DFsample[,1:2]), condition=input$boundaryLgcpMesh, convex=input$curvatureLgcpMesh, resolution=input$resolutionLgcpMesh, file.read=input$shapefileLgcpMesh)))
          sample <- DFsample
        }
      } else if(input$selectionLgcpMesh=="edgelength"){
        if(input$LgcpRasterSPDE=="raster"){
          rasterSample <- datareadRaster()
          qloc <- input$EdgeLengthLgcpMesh
          convexhull <- chull(DFsample[,1:2])
          convexhull <- c(convexhull, convexhull[1])
          mesh <- fm_mesh_2d_inla(loc.domain=DFsample[convexhull,1:2], max.edge=c(1,1.5)*input$EdgeLengthLgcpMesh, cutoff=input$EdgeLengthLgcpMesh/5, offset=c(-0.1, -0.2), max.edge=c(qloc[1], qloc[2]),
                               boundary=list(interiorNonConvex(x=as.matrix(rasterSample[,1:2]), condition=input$interiorLgcpMesh, convex=input$interiorcurvatureLgcpMesh, resolution=input$interiorresolutionLgcpMesh, file.read=input$interiorshapefileLgcpMesh),
                                             boundaryNonConvex(x=as.matrix(rasterSample[,1:2]), condition=input$boundaryLgcpMesh, convex=input$curvatureLgcpMesh, resolution=input$resolutionLgcpMesh, file.read=input$shapefileLgcpMesh)))
          sample <- rasterSample
        } else if(input$LgcpRasterSPDE=="solvecov"){
          qloc <- input$EdgeLengthLgcpMesh
          convexhull <- chull(DFsample[,1:2])
          convexhull <- c(convexhull, convexhull[1])
          mesh <- fm_mesh_2d_inla(loc.domain=DFsample[convexhull,1:2], max.edge=c(1,1.5)*input$EdgeLengthLgcpMesh,  cutoff=input$EdgeLengthLgcpMesh/5, offset=c(-0.1, -0.2),
                               boundary=list(interiorNonConvex(x=as.matrix(DFsample[,1:2]), condition=input$interiorLgcpMesh, convex=input$interiorcurvatureLgcpMesh, resolution=input$interiorresolutionLgcpMesh, file.read=input$interiorshapefileLgcpMesh),
                                             boundaryNonConvex(x=as.matrix(DFsample[,1:2]), condition=input$boundaryLgcpMesh, convex=input$curvatureLgcpMesh, resolution=input$resolutionLgcpMesh, file.read=input$shapefileLgcpMesh)))
          sample <- DFsample
        }
      }
      
      result <- list(mesh=mesh, qloc=qloc, Sample=sample)
      return(result)
    })
    
    ggplotLgcpMesh <- function(){
      if(input$buildLgcpMesh==0){
        ggplot()+ gg(LgcpMeshBase()$mesh)+ theme_bw() + xlab("Latitude") + ylab("Longitude") +
          ggtitle("Mesh over the study region") +
          geom_point(data=LgcpMeshBase()$Sample,
                     aes(x=LgcpMeshBase()$Sample[,1],y=LgcpMeshBase()$Sample[,2]), size=1) +
          theme(plot.title=element_text(color = "black", size = 16, face = "bold", hjust = 0.5))
      } else{
        ggplot()+ gg(LgcpMesh()$mesh)+ theme_bw() + xlab("Latitude") + ylab("Longitude") +
          ggtitle("Mesh over the study region") +
          geom_point(data=LgcpMesh()$Sample,
                     aes(x=LgcpMesh()$Sample[,1],y=LgcpMesh()$Sample[,2]), size=1) +
          theme(plot.title=element_text(color = "black", size = 16, face = "bold", hjust = 0.5))
      }
      
    }
    
    downloadFile(
      id = "ggplotLgcpMesh",
      logger = ss_userAction.Log,
      filenameroot = "ggplotLgcpMesh",
      aspectratio  = 1,
      downloadfxns = list(png  = ggplotLgcpMesh)
    )
    
    downloadablePlot(
      id = "ggplotLgcpMesh",
      logger = ss_userAction.Log,
      filenameroot = "ggplotLgcpMesh",
      aspectratio  = 1,
      downloadfxns = list(png  = ggplotLgcpMesh),
      visibleplot  = ggplotLgcpMesh)
    
    ## Lgcp modelling code ====
    
    LgcpModelFit <- eventReactive(input$fitLgcp, {
      showNotification(ui=paste("Fitting the data."), duration = NULL)
      t1 <- Sys.time()
      #taking the data from simulation or from the loading tab
      if(input$LgcpDataSimulatedLoaded=="sim"){
        DFsample <- as.data.frame(Lgcp.sampling())
      } else if(input$LgcpDataSimulatedLoaded=="load"){
        DFsample <- as.data.frame(datareadSample())
        if(input$LgcpRasterSPDE=="raster"){
          DFraster <- as.data.frame(datareadRaster())
        }
      }
      
      variablesChosenDefault <- input$DefaultComponentsLgcp
      variablesChosenUser <- input$UserComponentsLgcp
      variablesChosen <- c(variablesChosenDefault, variablesChosenUser)
      
      if(input$buildLgcpMesh==0){
        server_mesh <- LgcpMeshBase()
        mesh <- server_mesh$mesh
      } else{
        server_mesh <- LgcpMesh()
        mesh <- server_mesh$mesh
      }
      
      if(input$optionLgcpS=="auto"){
        if(input$KindPriorSpatialEffectLgcp=="PC.prior"){
          prior.range <- c(mean(c(diff(range(mesh$loc[mesh$segm$int$idx[,1],1])),diff(range(mesh$loc[mesh$segm$int$idx[,1],2]))))/5, 0.5)
          prior.sigma <- c(1,0.5)
          spde <- inla.spde2.pcmatern(mesh, prior.range = prior.range, prior.sigma = prior.sigma, alpha=2)
        } else{
          prior.range <- mean(c(diff(range(mesh$loc[mesh$segm$int$idx[,1],1])),diff(range(mesh$loc[mesh$segm$int$idx[,1],2]))))/5
          prior.sigma <- 1
          alpha <- 2; d <- 2
          nu <-  alpha - d/2
          kappa0 <-  log(8*nu)/2 -log(prior.range)
          tau0 <-  0.5*(lgamma(nu) - lgamma(nu + d/2) - d/2*log(4*pi)) - nu*kappa0 - log(prior.sigma)
          spde <-  inla.spde2.matern(mesh = mesh, B.tau = cbind(tau0, nu, -1), B.kappa = cbind(kappa0, -1, 0),
                                     theta.prior.mean = c(0.5,0.5), theta.prior.prec = c(1,1))
        }
      } else if(input$optionLgcpS=="custom"){
        if(input$KindPriorSpatialEffectLgcp=="PC.prior"){
          prior.range <- as.numeric(unlist(strsplit(input$LgcpPriorRangePC, split=",")))
          prior.sigma <- as.numeric(unlist(strsplit(input$LgcpPriorStdevPC, split=",")))
          spde <- inla.spde2.pcmatern(mesh, prior.range = prior.range, prior.sigma = prior.sigma, alpha=2)
        } else{
          prior.range <- as.numeric(unlist(strsplit(input$LgcpPriorRangeBase, split=",")))
          prior.sigma <- as.numeric(unlist(strsplit(input$LgcpPriorStdevBase, split=",")))
          alpha <- 2; d <- 2
          nu <-  alpha - d/2
          kappa0 <-  log(8*nu)/2 -log(prior.range[1])
          tau0 <-  0.5*(lgamma(nu) - lgamma(nu + d/2) - d/2*log(4*pi)) - nu*kappa0 - log(prior.sigma[1])
          spde <-  inla.spde2.matern(mesh = mesh, B.tau = cbind(tau0, nu, -1), B.kappa = cbind(kappa0, -1, 0),
                                     theta.prior.mean = c(prior.range[2],prior.sigma[2]), theta.prior.prec = c(prior.range[3],prior.sigma[3]))
        }
      }
      
      spde.index <- inla.spde.make.index(name="Spatial", n.spde = spde$n.spde)
      
      # LGCP mesh operations
      ldomain <- unique(mesh$loc[mesh$segm$int$idx,1:2])
      dmesh <- mesh.dual(mesh = mesh)
      domain.polys <- Polygons(list(Polygon(ldomain)), '0')
      domainSP <- SpatialPolygons(list(domain.polys))
      w <- sapply(1:length(dmesh), function(i) {
        if (gIntersects(dmesh[i, ], domainSP))
          return(gArea(gIntersection(dmesh[i, ], domainSP)))
        else return(0)
      })
      
      n <- nrow(DFsample)
      nv <- mesh$n
      imat <- Diagonal(nv, rep(1, nv))
      lmat <- inla.spde.make.A(mesh, as.matrix(DFsample[,1:2]))
      A.inf <- rbind(imat, lmat)
      
      ### Prediction of covariates ====
      
      List.covariates.inf <- list()
      List.covariates.mesh <- list()
      List.covariates.pred <- list()
      
      
      if((input$LgcpDataSimulatedLoaded=="load"&input$LgcpRasterSPDE=="solvecov")|input$LgcpDataSimulatedLoaded=="sim"|(input$LgcpDataSimulatedLoaded=="load"&input$LgcpRasterSPDE=="raster"|input$LgcpRasterPred=="SPDEraster")){
        x.pred <- seq(range(mesh$loc[mesh$segm$int$idx[,1], 1])[1], range(mesh$loc[mesh$segm$int$idx[,1], 1])[2], length.out=input$LgcpSPDErasterInput)
        y.pred <- seq(range(mesh$loc[mesh$segm$int$idx[,1], 2])[1], range(mesh$loc[mesh$segm$int$idx[,1], 2])[2], length.out=input$LgcpSPDErasterInput)
        xy.pred <- expand.grid(x=x.pred,y=y.pred)
        xy.pred <- xy.pred[which(!is.na(over(SpatialPoints(coords=xy.pred),SpatialPolygons(Srl=list(Polygons(srl=list(Polygon(coords=mesh$loc[mesh$segm$int$idx[,1], 1:2])), ID="int")))))),1:2]
        A.geo.pred <- inla.spde.make.A(mesh=mesh, loc=as.matrix(xy.pred))
      } else{ #if(input$LgcpDataSimulatedLoaded=="load"&input$LgcpRasterSPDE=="raster"|input$LgcpRasterPred=="rasterpred"){
        xy.pred <- DFraster
        A.geo.pred <- inla.spde.make.A(mesh=mesh, loc=as.matrix(xy.pred))
      }
      
      for(i in seq_along(variablesChosenUser)){
        if(!is.character(DFsample[,variablesChosenUser[i]])){
          prior.range.cov <- c(mean(c(diff(range(mesh$loc[mesh$segm$int$idx[,1],1])),diff(range(mesh$loc[mesh$segm$int$idx[,1],2]))))/5, 0.5)
          prior.sigma.cov <- c(1,0.5)
          spde.cov <- inla.spde2.pcmatern(mesh, prior.range = prior.range.cov, prior.sigma = prior.sigma.cov, alpha=2)
          spde.cov.index <- inla.spde.make.index(name="spatial.cov", n.spde = spde.cov$n.spde)
          formula.cov <- y ~ -1 + Intercept + f(spatial.cov, model=spde.cov)
          
          if((input$LgcpDataSimulatedLoaded=="load"&input$LgcpRasterSPDE=="solvecov")|input$LgcpDataSimulatedLoaded=="sim"){
            x.pred <- seq(range(mesh$loc[mesh$segm$int$idx[,1], 1])[1], range(mesh$loc[mesh$segm$int$idx[,1], 1])[2], length.out=input$LgcpSPDErasterInput)
            y.pred <- seq(range(mesh$loc[mesh$segm$int$idx[,1], 2])[1], range(mesh$loc[mesh$segm$int$idx[,1], 2])[2], length.out=input$LgcpSPDErasterInput)
            xy.pred <- expand.grid(x=x.pred,y=y.pred)
            xy.pred <- xy.pred[which(!is.na(over(SpatialPoints(coords=xy.pred),SpatialPolygons(Srl=list(Polygons(srl=list(Polygon(coords=mesh$loc[mesh$segm$int$idx[,1], 1:2])), ID="int")))))),1:2]
            A.geo.pred <- inla.spde.make.A(mesh=mesh, loc=as.matrix(xy.pred))
            Inf.stack.cov <- inla.stack(data=list(y=DFsample[,variablesChosenUser[i]]),
                                        A=list(lmat, 1),
                                        effects=list(list(spatial.cov=spde.cov.index$spatial.cov),
                                                     list(Intercept=rep(1,nrow(DFsample)))
                                        ),
                                        tag="Inference.cov")
            Pred.mesh.stack.cov <- inla.stack(data=list(y=rep(NA, mesh$n)),
                                              A=list(diag(1,nrow=mesh$n), 1),
                                              effects=list(list(spatial.cov=spde.cov.index$spatial.cov),
                                                           list(Intercept=rep(1,mesh$n))
                                              ),
                                              tag="Mesh.cov")
            Pred.stack.cov <- inla.stack(data=list(y=rep(NA,nrow(xy.pred))),
                                         A=list(A.geo.pred, 1),
                                         effects=list(list(spatial.cov=spde.cov.index$spatial.cov),
                                                      list(Intercept=rep(1,nrow(xy.pred)))
                                         ),
                                         tag="Prediction.cov")
            Total.stack.cov <- inla.stack(Inf.stack.cov, Pred.mesh.stack.cov, Pred.stack.cov)
            mod.cov <- inla(formula=formula.cov, data=inla.stack.data(Total.stack.cov), family="gaussian",
                            control.predictor = list(compute=FALSE, A=inla.stack.A(Total.stack.cov)))
            indx.inf <- inla.stack.index(Total.stack.cov, tag="Inference.cov")$data
            indx.mesh <- inla.stack.index(Total.stack.cov, tag="Mesh.cov")$data
            indx.pred <- inla.stack.index(Total.stack.cov, tag="Prediction.cov")$data
            List.covariates.inf[[variablesChosenUser[i]]] <- as.numeric(!is.na(DFsample[,variablesChosenUser[i]]))*DFsample[,variablesChosenUser[i]] + as.numeric(is.na(DFsample[,variablesChosenUser[i]]))*mod.cov$summary.fitted.values[indx.inf,"mean"]
            List.covariates.mesh[[variablesChosenUser[i]]] <- mod.cov$summary.fitted.values[indx.mesh,"mean"]
            List.covariates.pred[[variablesChosenUser[i]]] <- mod.cov$summary.fitted.values[indx.pred,"mean"]
          } else if(input$LgcpDataSimulatedLoaded=="load"&input$LgcpRasterSPDE=="raster"|input$LgcpRasterPred=="SPDEraster"){
            x.pred <- seq(range(mesh$loc[mesh$segm$int$idx[,1], 1])[1], range(mesh$loc[mesh$segm$int$idx[,1], 1])[2], length.out=input$LgcpSPDErasterInput)
            y.pred <- seq(range(mesh$loc[mesh$segm$int$idx[,1], 2])[1], range(mesh$loc[mesh$segm$int$idx[,1], 2])[2], length.out=input$LgcpSPDErasterInput)
            xy.pred <- expand.grid(x=x.pred,y=y.pred)
            xy.pred <- xy.pred[which(!is.na(over(SpatialPoints(coords=xy.pred),SpatialPolygons(Srl=list(Polygons(srl=list(Polygon(coords=mesh$loc[mesh$segm$int$idx[,1], 1:2])), ID="int")))))),1:2]
            A.geo.pred <- inla.spde.make.A(mesh=mesh, loc=as.matrix(xy.pred))
            if(variablesChosenUser[i] %in% colnames(DFraster)){
              Inf.stack.cov <- inla.stack(data=list(y=c(DFsample[,variablesChosenUser[i]], DFraster[!is.na(DFraster[,variablesChosenUser[i]]),variablesChosenUser[i]])),
                                          A=list(inla.spde.make.A(mesh=mesh, loc=as.matrix(rbind(DFsample[,1:2],DFraster[!is.na(DFraster[,variablesChosenUser[i]]),1:2]))), 1),
                                          effects=list(list(spatial.cov=spde.cov.index$spatial.cov),
                                                       list(Intercept=rep(1,nrow(DFsample)+nrow(DFraster)))
                                          ),
                                          tag="Inference.cov")
            } else{
              Inf.stack.cov <- inla.stack(data=list(y=DFsample[,variablesChosenUser[i]]),
                                          A=list(inla.spde.make.A(mesh=mesh, loc=as.matrix(DFsample[,1:2])), 1),
                                          effects=list(list(spatial.cov=spde.cov.index$spatial.cov),
                                                       list(Intercept=rep(1,nrow(DFsample)))
                                          ),
                                          tag="Inference.cov")
            }
            Pred.mesh.stack.cov <- inla.stack(data=list(y=rep(NA, mesh$n)),
                                              A=list(diag(1,nrow=mesh$n), 1),
                                              effects=list(list(spatial.cov=spde.cov.index$spatial.cov),
                                                           list(Intercept=rep(1,mesh$n))
                                              ),
                                              tag="Mesh.cov")
            Pred.stack.cov <- inla.stack(data=list(y=rep(NA,nrow(xy.pred))),
                                         A=list(A.geo.pred, 1),
                                         effects=list(list(spatial.cov=spde.cov.index$spatial.cov),
                                                      list(Intercept=rep(1,nrow(xy.pred)))
                                         ),
                                         tag="Prediction.cov")
            Total.stack.cov <- inla.stack(Inf.stack.cov, Pred.mesh.stack.cov, Pred.stack.cov)
            mod.cov <- inla(formula=formula.cov, data=inla.stack.data(Total.stack.cov), family="gaussian",
                            control.predictor = list(compute=FALSE, A=inla.stack.A(Total.stack.cov)))
            indx.inf <- inla.stack.index(Total.stack.cov, tag="Inference.cov")$data
            indx.mesh <- inla.stack.index(Total.stack.cov, tag="Mesh.cov")$data
            indx.pred <- inla.stack.index(Total.stack.cov, tag="Prediction.cov")$data
            List.covariates.inf[[variablesChosenUser[i]]] <- as.numeric(!is.na(DFsample[,variablesChosenUser[i]]))*DFsample[,variablesChosenUser[i]] + (as.numeric(is.na(DFsample[,variablesChosenUser[i]]))*(mod.cov$summary.fitted.values[indx.inf,"mean"])[1:nrow(DFsample)])
            List.covariates.mesh[[variablesChosenUser[i]]] <- mod.cov$summary.fitted.values[indx.mesh,"mean"]
            List.covariates.pred[[variablesChosenUser[i]]] <- mod.cov$summary.fitted.values[indx.pred,"mean"]
          } else if(input$LgcpDataSimulatedLoaded=="load"&input$LgcpRasterSPDE=="raster"|input$LgcpRasterPred=="rasterpred"){
            xy.pred <- DFraste[,1:2]
            A.geo.pred <- inla.spde.make.A(mesh=mesh, loc=as.matrix(xy.pred))
            
            Inf.stack.cov <- inla.stack(data=list(y=c(DFsample[,variablesChosenUser[i]], DFraster[,variablesChosenUser[i]])),
                                        A=list(inla.spde.make.A(mesh=mesh, loc=as.matrix(rbind(DFsample[,1:2],DFraster[,1:2]))), 1),
                                        effects=list(list(spatial.cov=spde.cov.index$spatial.cov),
                                                     list(Intercept=rep(1,nrow(DFsample)+nrow(DFraster)))
                                        ),
                                        tag="Inference.cov")
            Pred.mesh.stack.cov <- inla.stack(data=list(y=rep(NA, mesh$n)),
                                              A=list(diag(1,nrow=mesh$n), 1),
                                              effects=list(list(spatial.cov=spde.cov.index$spatial.cov),
                                                           list(Intercept=rep(1,mesh$n))
                                              ),
                                              tag="Mesh.cov")
            Total.stack.cov <- inla.stack(Inf.stack.cov, Pred.mesh.stack.cov)
            mod.cov <- inla(formula=formula.cov, data=inla.stack.data(Total.stack.cov), family="gaussian",
                            control.predictor = list(compute=FALSE, A=inla.stack.A(Total.stack.cov)))
            indx.mesh <- inla.stack.index(Total.stack.cov, tag="Mesh.cov")$data
            List.covariates.mesh[[variablesChosenUser[i]]] <- mod.cov$summary.fitted.values[indx.mesh,"mean"]
            
            List.covariates.inf[[variablesChosenUser[i]]] <- DFsample[,variablesChosenUser[i]]
            List.covariates.pred[[variablesChosenUser[i]]] <- DFraster[,variablesChosenUser[i]]
          }
        } else{
          if((input$LgcpDataSimulatedLoaded=="load"&input$LgcpRasterSPDE=="solvecov")|input$LgcpDataSimulatedLoaded=="sim"){
            cov.inf.indx <- as.vector(unlist(lapply(X=1:nrow(DFsample), FUN=function(Y){which.min(apply(X=(DFsample[!is.na(DFsample[,variablesChosenUser[i]]),1:2]-matrix(DFsample[Y,1:2],ncol=2))**2, MARGIN=1, FUN=sum))})))
            List.covariates.inf[[variablesChosenUser[i]]] <- DFsample[cov.inf.indx, variablesChosenUser[i]]
            mesh.indx <- as.vector(unlist(lapply(X=1:nrow(mesh$loc), FUN=function(Y){which.min(apply(X=(DFsample[!is.na(DFsample[,variablesChosenUser[i]]),1:2]-matrix(mesh$loc[Y,1:2],ncol=2))**2, MARGIN=1, FUN=sum))})))
            List.covariates.mesh[[variablesChosenUser[i]]] <- DFsample[mesh.indx, variablesChosenUser[i]]
            
            LgcpFactorModelDF <- data.frame(y=1, Lgcp=c(List.covariates.mesh[[variablesChosenUser[i]]], List.covariates.inf[[variablesChosenUser[i]]]))
            colnames(LgcpFactorModelDF) <- c("y", variablesChosenUser[i])
            idx.factor <- which(variablesChosenUser[i]==names(DFsample[!as.vector(unlist(lapply(X=DFsample, FUN=is.numeric)))])[names(DFsample[!as.vector(unlist(lapply(X=DFsample, FUN=is.numeric)))])%in%c(variablesChosenUser)])
            
            if(eval(parse(text=paste0("input$LgcpKindPredictionFactorLevel",idx.factor)))=="nearest"){
              cov.pred.ind <- as.vector(unlist(lapply(X=1:nrow(xy.pred), FUN=function(Y){which.min(apply(X=(DFsample[!is.na(DFsample[,variablesChosenUser[i]]),1:2]-matrix(xy.pred[Y,1:2],ncol=2))**2, MARGIN=1, FUN=sum))})))
              List.covariates.pred[[variablesChosenUser[i]]] <- DFsample[cov.pred.ind, variablesChosenUser[i]]
            } else{
              List.covariates.pred[[variablesChosenUser[i]]] <- rep(NA, times=nrow(xy.pred))
            }
            
          } else {
            DFsampleraster <- rbind(DFsample[,c(1:2,which(variablesChosenUser[i]==colnames(DFsample)))], DFraster[,c(1:2, which(variablesChosenUser[i]==colnames(DFraster)))])
            cov.inf.indx <- as.vector(unlist(lapply(X=1:nrow(DFsample), FUN=function(Y){which.min(apply(X=(DFsampleraster[!is.na(DFsampleraster[,variablesChosenUser[i]]),1:2]-matrix(DFsample[Y,1:2],ncol=2))**2, MARGIN=1, FUN=sum))})))
            List.covariates.inf[[variablesChosenUser[i]]] <- DFsample[cov.inf.indx, variablesChosenUser[i]]
            mesh.indx <- as.vector(unlist(lapply(X=1:nrow(mesh$loc), FUN=function(Y){which.min(apply(X=(DFsampleraster[!is.na(DFsampleraster[,variablesChosenUser[i]]),1:2]-matrix(mesh$loc[Y,1:2],ncol=2))**2, MARGIN=1, FUN=sum))})))
            List.covariates.mesh[[variablesChosenUser[i]]] <- DFsample[mesh.indx, variablesChosenUser[i]]
            
            LgcpFactorModelDF <- data.frame(y=1, Lgcp=c(List.covariates.mesh[[variablesChosenUser[i]]], List.covariates.inf[[variablesChosenUser[i]]]))
            colnames(LgcpFactorModelDF) <- c("y", variablesChosenUser[i])
            idx.factor <- which(variablesChosenUser[i]==names(DFsample[!as.vector(unlist(lapply(X=DFsample, FUN=is.numeric)))])[names(DFsample[!as.vector(unlist(lapply(X=DFsample, FUN=is.numeric)))])%in%c(variablesChosenUser)])
            
            if(eval(parse(text=paste0("input$LgcpKindPredictionFactorLevel",idx.factor)))=="nearest"){
              cov.pred.ind <- as.vector(unlist(lapply(X=1:nrow(xy.pred), FUN=function(Y){which.min(apply(X=(DFsampleraster[!is.na(DFsampleraster[,variablesChosenUser[i]]),1:2]-matrix(xy.pred[Y,1:2],ncol=2))**2, MARGIN=1, FUN=sum))})))
              List.covariates.pred[[variablesChosenUser[i]]] <- DFsample[cov.pred.ind, variablesChosenUser[i]]
            } else{
              List.covariates.pred[[variablesChosenUser[i]]] <- rep(NA, times=nrow(xy.pred))
            }
            
          }
        }
      }
      
      ### Building the main model and stacks structure ====
      
      Inf.lgcp.effects.list <- list(
        list(),
        list()
      )
      
      Pred.lgcp.effects.list <- list(
        list(),
        list()
      )
      
      formula_mod <- c("y ~ -1")
      
      A_Inf.spde1 <- list()
      A_Pred.spde1 <- list()
      
      for(i in seq_along(variablesChosen)){
        if(variablesChosen[i]=="Intercept"){
          formula_mod <- paste(formula_mod, "f(Intercept, model='linear')", sep=" + ")
          Inf.lgcp.effects.list[[2]][["Intercept"]] <- c(rep(1, times=nv), rep(1, times=n))
          Pred.lgcp.effects.list[[2]][["Intercept"]] <- rep(1, times=nrow(A.geo.pred))
          
        } else if(variablesChosen[i]=="Spatial Effect"){
          formula_mod <- paste(formula_mod, "f(Spatial, model=spde)", sep=" + ")
          Inf.lgcp.effects.list[[1]][["Spatial"]] <- spde.index$Spatial
          Pred.lgcp.effects.list[[1]][["Spatial"]] <- spde.index$Spatial
        } else{
          j <- which(variablesChosen[i]==variablesChosenUser)
          if(!is.character(DFsample[,variablesChosenUser[j]])){
            if(eval(parse(text=paste0("input$LgcpEffectCov",j)))=="rw1"|eval(parse(text=paste0("input$LgcpEffectCov",j)))=="rw2"){
              if(eval(parse(text=paste0("input$LgcpEffectCustomPrior",j)))=="custom"){
                if(eval(parse(text=paste0("input$LgcpEffectCovKindPrior",j)))=="pc"){
                  assign(paste0("pc.values", variablesChosenUser[j]), c( as.numeric(unlist(strsplit(eval(parse(text=paste0("input$LgcpEffectCovPriorPCValues",j))), ","))) ))
                  formula_mod <- paste(formula_mod, paste0("f(",variablesChosenUser[j],", model='",eval(parse(text=paste0("input$LgcpEffectCov",j))),"', hyper=list(prec=list(prior='pc.prec', param=", eval(parse(text=paste0("pc.values", variablesChosenUser[j]))), ")), constr=",eval(parse(text=paste0("input$LgcpEffectCovConstr",j))),", scale.model=TRUE)"), sep=" + ")
                } else if(eval(parse(text=paste0("input$LgcpEffectCovKindPrior",j)))=="unif"){
                  lim <- as.numeric(unlist(strsplit(eval(parse(text=paste0("input$LgcpEffectCovPriorUnif",j))), ",")))
                  sigma <- seq(0, lim[2]*3, length.out=1E5)
                  theta <- -2*log(sigma)
                  logdens <- sapply(X=sigma, FUN=logdunif, lim=lim)
                  unif.prior <- list(theta=list(prior=paste0("table: ", paste(c(theta, logdens), collapse=" "))))
                  assign(paste0("hyper", variablesChosenUser[j]), unif.prior )
                  formula_mod <- paste(formula_mod, paste0("f(",variablesChosenUser[j],", model='",eval(parse(text=paste0("input$LgcpEffectCov",j))),"', hyper=",paste0("hyper",variablesChosenUser[j]),  ", constr=",eval(parse(text=paste0("input$LgcpEffectCovConstr",j))),", scale.model=TRUE)"), sep=" + ")
                } else if(eval(parse(text=paste0("input$LgcpEffectCovKindPrior",j)))=="unifflat"){
                  assign(paste0("unifflat.prior",variablesChosenUser[j]), "expression:
                  log_dens = 0 - log(2) - theta/2
                  ")
                  formula_mod <- paste(formula_mod, paste0("f(",variablesChosenUser[j],", model='",eval(parse(text=paste0("input$LgcpEffectCov",j))),"', hyper=list(prec=list(prior=",paste0("unifflat.prior",variablesChosenUser[j]),")), constr=",eval(parse(text=paste0("input$LgcpEffectCovConstr",j))),", scale.model=TRUE)"), sep=" + ")
                } else{
                  assign(paste0("hyper",variablesChosenUser[j]), list(prec=list(prior="loggamma",param=c( as.numeric(unlist(strsplit(eval(parse(text=paste0("input$LgcpEffectCovPriorBaseValues",j))), ","))) ))) )
                  formula_mod <- paste(formula_mod, paste0("f(",variablesChosenUser[j],", model='",eval(parse(text=paste0("input$LgcpEffectCov",j))),"', hyper=",paste0("hyper",variablesChosenUser[j]),  ", constr=",eval(parse(text=paste0("input$LgcpEffectCovConstr",j))),", scale.model=TRUE)"), sep=" + ")
                }
              } else{
                formula_mod <- paste(formula_mod, paste0("f(",variablesChosenUser[j],", model='",eval(parse(text=paste0("input$LgcpEffectCov",j))),  "', constr=",eval(parse(text=paste0("input$LgcpEffectCovConstr",j))),", scale.model=TRUE)"), sep=" + ")
              }
              group.cov <- inla.group(x=c(List.covariates.mesh[[variablesChosenUser[j]]], List.covariates.inf[[variablesChosenUser[j]]], List.covariates.pred[[variablesChosenUser[j]]]), n=eval(parse(text=paste0("input$LgcpEffectCovNodes",j))), method="cut")
              
              Inf.lgcp.effects.list[[2]][[variablesChosenUser[j]]] <- group.cov[seq_len(n + nv)]
              Pred.lgcp.effects.list[[2]][[variablesChosenUser[j]]] <- group.cov[-seq_len(n + nv)]
              
            } else if(eval(parse(text=paste0("input$LgcpEffectCov",j)))=="spde1"){
              Tot_cov <- c(List.covariates.mesh[[variablesChosenUser[j]]], List.covariates.inf[[variablesChosenUser[j]]], List.covariates.pred[[variablesChosenUser[j]]])
              spde1_nodes <- seq(min(Tot_cov), max(Tot_cov), length.out=eval(parse(text=paste0("input$LgcpEffectCovNodes",j))))
              mesh1d <- fm_mesh_1d(loc=spde1_nodes)
              
              if(eval(parse(text=paste0("input$LgcpEffectCustomPrior",j)))=="custom"){
                if(eval(parse(text=paste0("input$LgcpEffectCovKindPrior",j)))=="pc"){
                  hyper <- as.numeric(unlist(strsplit(eval(parse(text=paste0("input$LgcpEffectCovPriorPCValues",i))), ",")))
                  assign(paste0("spde1d_",variablesChosenUser[j]) ,inla.spde2.pcmatern(mesh=mesh1d, prior.range=c(abs(diff(range(Tot_cov)))/5, 0.5), prior.sigma=c(1,0.5)))
                } else{
                  prior.range <- as.numeric(unlist(strsplit(eval(parse(text=paste0("input$LgcpEffectCovPriorBaseValues",i))), ",")))[1:3]
                  prior.sigma <- as.numeric(unlist(strsplit(eval(parse(text=paste0("input$LgcpEffectCovPriorBaseValues",i))), ",")))[4:6]
                  alpha <- 2; d <- 2
                  nu <-  alpha - d/2
                  kappa0 <-  log(8*nu)/2 -log(prior.range[1])
                  tau0 <-  0.5*(lgamma(nu) - lgamma(nu + d/2) - d/2*log(4*pi)) - nu*kappa0 - log(prior.sigma[1])
                  assign(paste0("spde1d_",variablesChosenUser[j]), inla.spde2.matern(mesh = mesh1d, B.tau = cbind(tau0, nu, -1), B.kappa = cbind(kappa0, -1, 0),
                                                                                     theta.prior.mean = c(prior.range[2],prior.sigma[2]), theta.prior.prec = c(prior.range[3],prior.sigma[3])))
                }
              } else {
                assign(paste0("spde1d_",variablesChosenUser[j]), inla.spde2.pcmatern(mesh=mesh1d, prior.range=c(abs(diff(range(Tot_cov)))/5, 0.5), prior.sigma=c(1,0.5)))
              }
              
              spde1d.index <- inla.spde.make.index(name=paste0(variablesChosenUser[j]), n.spde=eval(parse(text=paste0("spde1d_",variablesChosenUser[j])))$n.spde)
              formula_mod <- paste(formula_mod, paste0("f(",variablesChosenUser[j],", model=", paste0("spde1d_",variablesChosenUser[j]),  ")"), sep=" + ")

              Inf.lgcp.effects.list[[length(A_Inf.spde1)+3]] <- list()
              Pred.lgcp.effects.list[[length(A_Inf.spde1)+3]] <- list()
              Inf.lgcp.effects.list[[length(A_Inf.spde1)+3]][[variablesChosenUser[j]]] <- spde1d.index[[1]]
              Pred.lgcp.effects.list[[length(A_Inf.spde1)+3]][[variablesChosenUser[j]]] <- spde1d.index[[1]]
              
              A_Inf.spde1[[length(A_Inf.spde1)+1]] <- inla.spde.make.A(mesh=mesh1d, loc=Tot_cov[seq_len(nv+n)])
              A_Pred.spde1[[length(A_Pred.spde1)+1]] <- inla.spde.make.A(mesh=mesh1d, loc=Tot_cov[-seq_len(nv+n)])
              
            } else if(eval(parse(text=paste0("input$LgcpEffectCov",j)))=="linear"){
              if(eval(parse(text=paste0("input$LgcpEffectCustomPrior",j)))=="custom"){
                hyper <- as.numeric(unlist(strsplit(eval(parse(text=paste0("input$LgcpEffectCovPriorBaseValues",j))), ",")))
                formula_mod <- paste(formula_mod, paste0("f(",variablesChosenUser[j],", model='",eval(parse(text=paste0("input$LgcpEffectCov",j))),"', mean.linear=", hyper[1], ",prec.linear=", hyper[2],")"), sep=" + ")
              } else{
                formula_mod <- paste(formula_mod, paste0("f(",variablesChosenUser[j],", model='",eval(parse(text=paste0("input$LgcpEffectCov",j))),"')"), sep=" + ")
              }
              
              Inf.lgcp.effects.list[[2]][[variablesChosenUser[j]]] <- c(List.covariates.mesh[[variablesChosenUser[j]]], List.covariates.inf[[variablesChosenUser[j]]])
              Pred.lgcp.effects.list[[2]][[variablesChosenUser[j]]] <- c(List.covariates.pred[[variablesChosenUser[j]]])
              
            } else{
              showNotification(ui=paste("The effect of numerical covariates cannot possess an independent and identically distributed (iid) structure. If this is required, the variable values should be recorded as text, not as numerical input.."), duration = NULL)
            }
          } else{ # Section for factor variables
            if(eval(parse(text=paste0("input$LgcpEffectCov",j)))=="iid"){
              LgcpFactorModelDF <- data.frame(y=1, Lgcp=c(List.covariates.mesh[[variablesChosenUser[j]]], List.covariates.inf[[variablesChosenUser[j]]]))
              colnames(LgcpFactorModelDF) <- c("y", variablesChosenUser[j])
              idx.factor <- which(variablesChosenUser[j]==names(DFsample[!as.vector(unlist(lapply(X=DFsample, FUN=is.numeric)))])[names(DFsample[!as.vector(unlist(lapply(X=DFsample, FUN=is.numeric)))])%in%c(variablesChosenUser)])
              
              Inf.lgcp.effects.list[[2]][[variablesChosenUser[j]]] <- c(List.covariates.mesh[[variablesChosenUser[j]]], List.covariates.inf[[variablesChosenUser[j]]])
              Pred.lgcp.effects.list[[2]][[variablesChosenUser[j]]] <- if(eval(parse(text=paste0("input$LgcpKindPredictionFactorLevel",idx.factor)))=="reference"){
                rep( eval(parse(text=paste0("input$LgcpEffectCovFactorPred",idx.factor))), times=length(List.covariates.pred[[variablesChosenUser[j]]]))
              } else{List.covariates.pred[[variablesChosenUser[j]]]}
              
              #c(List.covariates.pred[[variablesChosenUser[j]]])
              
              if(eval(parse(text=paste0("input$LgcpEffectCustomPrior",j)))=="custom"){
                if(eval(parse(text=paste0("input$LgcpEffectCovKindPrior",j)))=="pc"){
                  hyper <- list(prec=list(prior="pc.prior",param=c( as.numeric(unlist(strsplit(eval(parse(text=paste0("input$LgcpEffectCovPriorPCValues",j))), ","))) )))
                  formula_mod <- paste(formula_mod, paste0("f(",variablesChosenUser[j],", model='",eval(parse(text=paste0("input$LgcpEffectCov",j))),"', hyper=",hyper, ", constr=",eval(parse(text=paste0("input$LgcpEffectCovConstr",j))),")"), sep=" + ")
                } else if(eval(parse(text=paste0("input$LgcpEffectCovKindPrior",j)))=="base"){
                  hyper <- list(prec=list(prior="loggamma",param=c( as.numeric(unlist(strsplit(eval(parse(text=paste0("input$LgcpEffectCovPriorBaseValues",j))), ","))) )))
                  formula_mod <- paste(formula_mod, paste0("f(",variablesChosenUser[j],", model='iid', hyper=",hyper,  ", constr=",eval(parse(text=paste0("input$LgcpEffectCovConstr",j))),")"), sep=" + ")
                } else if(eval(parse(text=paste0("input$LgcpEffectCovKindPrior",j)))=="unif"){
                  lim <- as.numeric(unlist(strsplit(eval(parse(text=paste0("input$LgcpEffectCovPriorUnif",j))), ",")))
                  sigma <- seq(lim[1]+1E-5, lim[2]*3, length.out=1E5)
                  theta <- -2*log(sigma)
                  logdens <- sapply(X=sigma, FUN=logdunif, lim=lim)
                  unif.prior <- list(theta=list(prior=paste0("table: ", paste(c(theta, logdens), collapse=" "))))
                  assign(paste0("hyper", variablesChosenUser[j]), unif.prior )
                  formula_mod <- paste(formula_mod, paste0("f(",variablesChosenUser[j],", model='iid', hyper=",paste0("hyper",variablesChosenUser[j]),  ", constr=",eval(parse(text=paste0("input$LgcpEffectCovConstr",j))),", scale.model=TRUE)"), sep=" + ")
                } else{ # flatunif
                  unifflat.prior= "expression:
                  log_dens = 0 - log(2) - theta/2
                  "
                  formula_mod <- paste(formula_mod, paste0("f(",variablesChosenUser[j],", model='iid', hyper=list(prec=list(prior=",unifflat.prior,")), constr=",eval(parse(text=paste0("input$LgcpEffectCovConstr",j))),")"), sep=" + ")
                }
              } else{
                formula_mod <- paste(formula_mod, paste0("f(", variablesChosenUser[j], ", model='iid'",  ", constr=",eval(parse(text=paste0("input$LgcpEffectCovConstr",j))),")"), sep=" + ")
              }
              
            } else{
              LgcpFactorModelDF <- data.frame(y=1, Lgcp=c(List.covariates.mesh[[variablesChosenUser[j]]], List.covariates.inf[[variablesChosenUser[j]]]))
              colnames(LgcpFactorModelDF) <- c("y", variablesChosenUser[j])
              FactorVariables <- data.frame( model.matrix(object=as.formula(paste0("y~-1+", variablesChosenUser[j])), LgcpFactorModelDF ) )
              idx.factor <- which(variablesChosenUser[j]==names(DFsample[!as.vector(unlist(lapply(X=DFsample, FUN=is.numeric)))])[names(DFsample[!as.vector(unlist(lapply(X=DFsample, FUN=is.numeric)))])%in%c(variablesChosenUser)])
              for(l in setdiff(colnames(FactorVariables),paste0(variablesChosenUser[j], eval(parse(text=paste0("input$LgcpEffectCovFactorPred",idx.factor))) ))){
                ll <- l
                Inf.lgcp.effects.list[[2]][[l]] <- FactorVariables[,l]
                Pred.lgcp.effects.list[[2]][[l]] <- rep(0, length(List.covariates.pred[[variablesChosenUser[j]]]))
                if(eval(parse(text=paste0("input$LgcpEffectCustomPrior",j)))=="custom"|eval(parse(text=paste0("input$LgcpEffectCov",j)))=="linear"){
                  hyper <- as.numeric(unlist(strsplit(eval(parse(text=paste0("input$LgcpEffectCovPriorBaseValues",j))), ",")))
                  formula_mod <- paste(formula_mod, paste0("f(", l, ", model='linear', mean.linear=", hyper[1], ",prec.linear=", hyper[2],")"), sep=" + ")
                } else{
                  formula_mod <- paste(formula_mod, paste0("f(", l, ", model='linear')"), sep=" + ")
                }
              }
            }
          }
        }
      }
      
      ### Stacks of the geostatistical, Lgcp and prediction layers ====
      
      ResponseVariable <- DFsample[,3]
      
      A_inf_tot <- c(A.inf,1)
      if(length(A_Inf.spde1)>0){
        for(i in seq_along(A_Inf.spde1)){
          A_inf_tot[[2+i]] <- A_Inf.spde1[[i]]
        }
      }
      
      A_pred_tot <- c(A.geo.pred,1)
      if(length(A_Inf.spde1)>0){
        for(i in seq_along(A_Inf.spde1)){
          A_pred_tot[[2+i]] <- A_Pred.spde1[[i]]
        }
      }
      
      Inf.LGCP.stack <- inla.stack(data=list(y=c(rep(0,times=nv), rep(1,times=n)), e=c(w, rep(0,n))),
                                   A=A_inf_tot,
                                   effects=Inf.lgcp.effects.list,
                                   tag="Inference_Lgcp"
      ) 
        
      Pred.LGCP.stack <- inla.stack(data=list(y=matrix(NA, nrow=nrow(A.geo.pred), ncol=1)),
                                   A=A_pred_tot,
                                   effects=Pred.lgcp.effects.list,
                                   tag="Prediction_Lgcp")
      
      Total.stack <- inla.stack(Inf.LGCP.stack, Pred.LGCP.stack)
      
      ### INLA model ====
      
      formula_inla <- as.formula(formula_mod)
      
      if(input$autocustomLgcpFamily=='custom'){
        if(input$LgcpFamilyPriorKind=="pc"){
          family.pcprec <- as.numeric(unlist(strsplit(input$LgcpFamilyHyper,",")))
          controlFamily <- list(hyper = list(prec = list(prior="pc.prec", param=family.pcprec)))
        } else if(input$LgcpFamilyPriorKind=="unif"){
          lim <- as.numeric(unlist(strsplit(input$LgcpFamilyHyper,",")))
          sigma <- seq(0, lim[2]*3, length.out=1E5)
          theta <- -2*log(sigma)
          logdens <- sapply(X=sigma, FUN=logdunif, lim=lim)
          family.unif.prior <- list(theta=list(prior=paste0("table: ", paste(c(theta, logdens), collapse=" "))))
          controlFamily <- list(hyper = list(theta = list(prior=family.unif.prior)))
        } else if(input$LgcpFamilyPriorKind=="unifflat"){
          unifflat.prior= "expression:
                  log_dens = 0 - log(2) - theta/2
                  "
          controlFamily <- list(hyper = list(prec = list(prior=unifflat.prior)))
        } else{
          controlFamily <- list(hyper = list(prec = list(prior="loggamma", param=as.numeric(unlist(strsplit(input$LgcpFamilyHyper,","))))), 
                                control.link=list(model="log"))
        }
      } else{
        controlFamily <- list(control.link=list(model="log"))
      }
      
      if(input$autocustomLgcpMode=='custom'){
        controlModeTheta <- list(theta=as.numeric(unlist(strsplit(input$LgcpModeHyper,","))), restart=TRUE)
      } else{controlModeTheta <- inla.set.control.mode.default()}
      
      if(input$INLAModeLgcp=="classic"){
        controlINLA <- list(strategy=input$strategyapproxINLALgcp,
                            int.strategy=input$strategyintINLALgcp)
      } else{
        controlINLA <- list()
      }
      
      Lgcp.model <- inla(formula=formula_inla, family = c(input$SelectLgcpFamily),
                         data = inla.stack.data(Total.stack),
                         E = inla.stack.data(Total.stack)$e,
                         control.inla = controlINLA,
                         control.predictor = list(A = inla.stack.A(Total.stack), compute = TRUE, link = 1),
                         control.family = controlFamily,
                         control.mode = controlModeTheta,
                         control.compute = list(cpo = TRUE, dic = TRUE, waic = TRUE, config = TRUE),
                         inla.mode=input$INLAModeLgcp,
                         verbose=FALSE)
      
      index.pred <- inla.stack.index(Total.stack, "Prediction_Lgcp")$data
      DFpred <- data.frame(Latitude=xy.pred[,1], Longitude=xy.pred[,2])
      DFpred$Abundance.mean <- Lgcp.model$summary.fitted.values[index.pred, "mean"]
      DFpred$Abundance.median <- Lgcp.model$summary.fitted.values[index.pred, "0.5quant"]
      DFpred$Abundance.sd <- Lgcp.model$summary.fitted.values[index.pred, "sd"]
      
      colnames(DFpred)[1:2] <- colnames(DFsample)[1:2]
      
      DFpredictorMeanMedianStdev <- data.frame(Latitude=xy.pred[,1], Longitude=xy.pred[,2],
                                               Predictor.mean=Lgcp.model$summary.linear.predictor[index.pred, "mean"],
                                               Predictor.median=Lgcp.model$summary.linear.predictor[index.pred, "0.5quant"],
                                               Predictor.stdev=Lgcp.model$summary.linear.predictor[index.pred, "sd"])
      
      colnames(DFpredictorMeanMedianStdev)[1:2] <- colnames(DFsample)[1:2]
      
      gridSpatial <- expand.grid(x=seq(min(DFpred[,1]), max(DFpred[,1]),length.out=input$dimLgcpmap),
                                 y=seq(min(DFpred[,2]), max(DFpred[,2]), length.out=input$dimLgcpmap))
      gridSpatial <- gridSpatial[which(!is.na(over(SpatialPoints(coords=gridSpatial),SpatialPolygons(Srl=list(Polygons(srl=list(Polygon(coords=mesh$loc[mesh$segm$int$idx[,1], 1:2])), ID="int")))))),1:2]
      
      A.spatial <- inla.spde.make.A(mesh=mesh, loc=as.matrix(gridSpatial))
      
      DFspatialMeanMedianStdev <- data.frame(Latitude=as.vector(gridSpatial[,1]), Longitude=as.vector(gridSpatial[,2]),
                                             Spatial.mean=as.vector(A.spatial%*%Lgcp.model$summary.random$Spatial$mean),
                                             Spatial.median=as.vector(A.spatial%*%Lgcp.model$summary.random$Spatial$`0.5quant`),
                                             Spatial.stdev=as.vector(A.spatial%*%Lgcp.model$summary.random$Spatial$sd))
      
      colnames(DFspatialMeanMedianStdev)[1:2] <- colnames(DFsample)[1:2]
      
      result <- list(DFpredAbunMeanMedianStdev=list(DFpredAbunMeanMedianStdev=DFpred),
                     DFpredictorMeanMedianStdev=list(DFpredictorMeanMedianStdev=DFpredictorMeanMedianStdev),
                     DFspatialMeanMedianStdev=list(DFspatialMeanMedianStdev=DFspatialMeanMedianStdev),
                     LgcpModel=list(LgcpModel=Lgcp.model),
                     DFPostFixed=list(DFPostFixed=Lgcp.model$marginals.fixed),
                     DFPostHyperpar=list(DFPostHyperpar=Lgcp.model$marginals.hyperpar),
                     Summary.fixed=list(Summary.fixed=Lgcp.model$summary.fixed),
                     Summary.hyperpar=list(Summary.hyperpar=Lgcp.model$summary.hyperpar),
                     SummaryInternalHyper=list(SummaryInternalHyper=Lgcp.model$internal.summary.hyperpar),
                     SummaryCPO=list(SummaryCPO=na.omit(Lgcp.model$cpo$cpo)),
                     DICmodel=list(DICmodel=data.frame(DIC=Lgcp.model$dic$family.dic, row.names=c("Point process"))),
                     WAICmodel=list(WAICmodel=data.frame(WAIC=cbind(unlist(lapply(na.omit(unique(Lgcp.model$dic$family)), function(i){sum(Lgcp.model$waic$local.waic[which(Lgcp.model$dic$family==i)])}))), row.names=c("Point process")))
                     )
      
      t2 <- Sys.time()
      difftime(t2,t1, units="secs")
      showNotification(ui=paste("The model has been fitted:", as.numeric(round(Lgcp.model$cpu.used[4])),
                                "(LGCP model) and", as.numeric(round(difftime(t2,t1, units="secs"))),
                                "(overall process) secs." ), duration = NULL)
      # showNotification(ui=paste("The model's DIC is", Lgcperential.model$dic$dic), duration = NULL)
      return(result)
    })
    
    
    dataggplotLgcpAbundanceMeanMedianStdevFit <- function(){
      DF <- LgcpModelFit()$DFpredAbunMeanMedianStdev$DFpredAbunMeanMedianStdev
      return(DF)
    }
    
    DFLgcpAbundance <- reactive({
      DF <- LgcpModelFit()$DFpredAbunMeanMedianStdev$DFpredAbunMeanMedianStdev
      condition <- !(length(unique(DF[,1]))*length(unique(DF[,2]))==nrow(DF))&val()
      if(condition){
        grid <- expand.grid(seq(min(DF[,1]),max(DF[,1]), length.out=200),seq(min(DF[,2]),max(DF[,2]), length.out=200))
        DFInter <- data.frame(Latitude=grid[,1],Longitude=grid[,2])
        DFInter$Abundance.mean <- InterpolateIrrGrid(z=DF$Abundance.mean,loc=DF[,1:2], gridInter=grid)$DataInter$z
        DFInter$Abundance.median <- InterpolateIrrGrid(z=DF$Abundance.median,loc=DF[,1:2], gridInter=grid)$DataInter$z
        DFInter$Abundance.sd <- InterpolateIrrGrid(z=DF$Abundance.sd,loc=DF[,1:2], gridInter=grid)$DataInter$z
        colnames(DFInter)[1:2] <- colnames(DF)[1:2]
        DF <- DFInter
      }
      
      return(DF)
    })
    
    ggplotLgcpAbundanceMeanMedianStdevFit <- function(){
      DF <- DFLgcpAbundance()
      g1 <- ggplot(DF) + geom_tile(aes(x=DF[,1], y=DF[,2], fill=Abundance.mean)) +
        scale_fill_viridis_c(option = "turbo") + theme_bw() + coord_fixed() + labs(fill="Values") +
        xlab(colnames(DF)[1]) + ylab(colnames(DF)[1]) + ggtitle("Response predicted (mean)") +
        theme(plot.title=element_text(color = "black", size = 16, face = "bold", hjust = 0.5))
      
      g2 <- ggplot(DF) + geom_tile(aes(x=DF[,1], y=DF[,2], fill=Abundance.median)) +
        scale_fill_viridis_c(option = "turbo") + theme_bw() + coord_fixed() + labs(fill="Values") +
        xlab(colnames(DF)[1]) + ylab(colnames(DF)[1]) + ggtitle("Response predicted (median)") +
        theme(plot.title=element_text(color = "black", size = 16, face = "bold", hjust = 0.5))
      
      g3 <- ggplot(DF) + geom_tile(aes(x=DF[,1], y=DF[,2], fill=Abundance.sd))+
        scale_fill_viridis_c(option = "turbo") + theme_bw() + coord_fixed() + labs(fill="Values") +
        xlab(colnames(DF)[1]) + ylab(colnames(DF)[1]) + ggtitle("Response predicted (stdev.)") +
        theme(plot.title=element_text(color = "black", size = 16, face = "bold", hjust = 0.5))
      
      gt <- g1+g2+g3
      
      return(gt)
    }
    
    downloadFile(
      id = "ggplotLgcpAbundanceMeanMedianStdevFit",
      logger = ss_userAction.Log,
      filenameroot = "ggplotLgcpAbundanceMeanMedianStdevFit",
      aspectratio  = 1,
      downloadfxns = list(png  = ggplotLgcpAbundanceMeanMedianStdevFit,
                          csv = dataggplotLgcpAbundanceMeanMedianStdevFit,
                          txt = dataggplotLgcpAbundanceMeanMedianStdevFit)
    )
    
    downloadablePlot(id = "ggplotLgcpAbundanceMeanMedianStdevFit",
                     logger = ss_userAction.Log,
                     filenameroot = "ggplotLgcpAbundanceMeanMedianStdevFit",
                     aspectratio  = 1,
                     downloadfxns = list(png = ggplotLgcpAbundanceMeanMedianStdevFit,
                                         csv = dataggplotLgcpAbundanceMeanMedianStdevFit,
                                         txt = dataggplotLgcpAbundanceMeanMedianStdevFit),
                     visibleplot  = ggplotLgcpAbundanceMeanMedianStdevFit)
    
    
    dataggplotLgcpSpatialMeanMedianStdev <- function(){
      DF <- LgcpModelFit()$DFspatialMeanMedianStdev$DFspatialMeanMedianStdev
      return(DF)
    }
    
    ggplotLgcpSpatialMeanMedianStdev <- function(){
      DF <- LgcpModelFit()$DFspatialMeanMedianStdev$DFspatialMeanMedianStdev
      g1 <- ggplot(DF) + geom_tile(aes(x=DF[,1], y=DF[,2], fill=Spatial.mean)) +
        scale_fill_viridis_c(option = "turbo") + theme_bw() + coord_fixed() + labs(fill="Values") +
        xlab(colnames(DF)[1]) + ylab(colnames(DF)[1]) + ggtitle("Posterior spatial (mean)") +
        theme(plot.title=element_text(color = "black", size = 16, face = "bold", hjust = 0.5))
      
      g2 <- ggplot(DF) + geom_tile(aes(x=DF[,1], y=DF[,2], fill=Spatial.median)) +
        scale_fill_viridis_c(option = "turbo") + theme_bw() + coord_fixed() + labs(fill="Values") +
        xlab(colnames(DF)[1]) + ylab(colnames(DF)[1]) + ggtitle("Posterior spatial (median)") +
        theme(plot.title=element_text(color = "black", size = 16, face = "bold", hjust = 0.5))
      
      g3 <- ggplot(DF) + geom_tile(aes(x=DF[,1], y=DF[,2], fill=Spatial.stdev))+
        scale_fill_viridis_c(option = "turbo") + theme_bw() + coord_fixed() + labs(fill="Values") +
        xlab(colnames(DF)[1]) + ylab(colnames(DF)[1]) + ggtitle("Posterior spatial (stdev.)") +
        theme(plot.title=element_text(color = "black", size = 16, face = "bold", hjust = 0.5))
      
      gt <- g1+g2+g3
      
      return(gt)
    }
    
    downloadFile(
      id = "ggplotLgcpSpatialMeanMedianStdev",
      logger = ss_userAction.Log,
      filenameroot = "ggplotLgcpSpatialMeanMedianStdev",
      aspectratio  = 1,
      downloadfxns = list(png  = ggplotLgcpSpatialMeanMedianStdev,
                          csv = dataggplotLgcpSpatialMeanMedianStdev,
                          txt = dataggplotLgcpSpatialMeanMedianStdev)
    )
    
    downloadablePlot(id = "ggplotLgcpSpatialMeanMedianStdev",
                     logger = ss_userAction.Log,
                     filenameroot = "ggplotLgcpSpatialMeanMedianStdev",
                     aspectratio  = 1,
                     downloadfxns = list(png=ggplotLgcpSpatialMeanMedianStdev,
                                         csv = dataggplotLgcpSpatialMeanMedianStdev,
                                         txt = dataggplotLgcpSpatialMeanMedianStdev),
                     visibleplot  = ggplotLgcpSpatialMeanMedianStdev)
    
    dataggplotLgcpPredictorMeanMedianStdevFit <- function(){
      DF <- LgcpModelFit()$DFpredictorMeanMedianStdev$DFpredictorMeanMedianStdev
      return(DF)
    }
    
    DFLgcpPredictor <- reactive({
      DF <- LgcpModelFit()$DFpredictorMeanMedianStdev$DFpredictorMeanMedianStdev
      condition <- !(length(unique(DF[,1]))*length(unique(DF[,2]))==nrow(DF))&val()
      if(condition){
        grid <- expand.grid(seq(min(DF[,1]),max(DF[,1]), length.out=200),seq(min(DF[,2]),max(DF[,2]), length.out=200))
        DFInter <- data.frame(Latitude=grid[,1],Longitude=grid[,2])
        DFInter$Predictor.mean <- InterpolateIrrGrid(z=DF$Predictor.mean,loc=DF[,1:2], gridInter=grid)$DataInter$z
        DFInter$Predictor.median <- InterpolateIrrGrid(z=DF$Predictor.median,loc=DF[,1:2], gridInter=grid)$DataInter$z
        DFInter$Predictor.stdev <- InterpolateIrrGrid(z=DF$Predictor.stdev,loc=DF[,1:2], gridInter=grid)$DataInter$z
        colnames(DFInter)[1:2] <- colnames(DF)[1:2]
        DF <- DFInter
      }
      
      return(DF)
    })
    
    ggplotLgcpPredictorMeanMedianStdevFit <- function(){
      DF <- DFLgcpPredictor()
      g1 <- ggplot(DF) + geom_tile(aes(x=DF[,1], y=DF[,2], fill=Predictor.mean)) +
        scale_fill_viridis_c(option = "turbo") + theme_bw() + coord_fixed() + labs(fill="Values") +
        xlab(colnames(DF)[1]) + ylab(colnames(DF)[1]) + ggtitle("Linear predictor (mean)") +
        theme(plot.title=element_text(color = "black", size = 16, face = "bold", hjust = 0.5))
      
      g2 <- ggplot(DF) + geom_tile(aes(x=DF[,1], y=DF[,2], fill=Predictor.median)) +
        scale_fill_viridis_c(option = "turbo") + theme_bw() + coord_fixed() + labs(fill="Values") +
        xlab(colnames(DF)[1]) + ylab(colnames(DF)[1]) + ggtitle("Linear predictor (median)") +
        theme(plot.title=element_text(color = "black", size = 16, face = "bold", hjust = 0.5))
      
      g3 <- ggplot(DF) + geom_tile(aes(x=DF[,1], y=DF[,2], fill=Predictor.stdev))+
        scale_fill_viridis_c(option = "turbo") + theme_bw() + coord_fixed() + labs(fill="Values") +
        xlab(colnames(DF)[1]) + ylab(colnames(DF)[1]) + ggtitle("Linear predictor (stdev.)") +
        theme(plot.title=element_text(color = "black", size = 16, face = "bold", hjust = 0.5))
      
      gt <- g1+g2+g3
      
      return(gt)
    }
    
    downloadFile(
      id = "ggplotLgcpPredictorMeanMedianStdevFit",
      logger = ss_userAction.Log,
      filenameroot = "ggplotLgcpPredictorMeanMedianStdevFit",
      aspectratio  = 1,
      downloadfxns = list(png  = ggplotLgcpPredictorMeanMedianStdevFit,
                          csv = dataggplotLgcpPredictorMeanMedianStdevFit,
                          txt = dataggplotLgcpPredictorMeanMedianStdevFit)
    )
    
    downloadablePlot(id = "ggplotLgcpPredictorMeanMedianStdevFit",
                     logger = ss_userAction.Log,
                     filenameroot = "ggplotLgcpPredictorMeanMedianStdevFit",
                     aspectratio  = 1,
                     downloadfxns = list(png = ggplotLgcpPredictorMeanMedianStdevFit,
                                         csv = dataggplotLgcpPredictorMeanMedianStdevFit,
                                         txt = dataggplotLgcpPredictorMeanMedianStdevFit),
                     visibleplot  = ggplotLgcpPredictorMeanMedianStdevFit)
    
    dataggplotLgcpFixParamFit <- function(){
      DF <- as.data.frame(LgcpModelFit()$DFPostFixed$DFPostFixed)
      return(DF)
    }
    
    ggplotLgcpFixParamFit <- function(){
      DF <- LgcpModelFit()$DFPostFixed$DFPostFixed
      gl <- c()
      for(i in 1:length(DF)){
        assign(paste0("g",i),
               ggplot(data=data.frame(x=DF[[i]][,1], y=DF[[i]][,2]), aes(x=x,y=y)) + geom_line() +
                 theme_bw() + xlab(names(DF)[i]) + ylab(HTML(paste("Density f(",names(DF)[i],")"))) +
                 ggtitle(names(DF)[i]) + theme(plot.title=element_text(color = "black", size = 16, face = "bold", hjust = 0.5))
        )
        gl[i] <- paste0("g",i)
      }
      gt <- eval(parse(text=paste(gl, collapse="+")))
      return(gt)
    }
    
    downloadFile(
      id = "ggplotLgcpFixParamFit",
      logger = ss_userAction.Log,
      filenameroot = "ggplotLgcpFixParamFit",
      aspectratio  = 1,
      downloadfxns = list(png  = ggplotLgcpFixParamFit,
                          csv = dataggplotLgcpFixParamFit,
                          txt = dataggplotLgcpFixParamFit)
    )
    
    downloadablePlot(id = "ggplotLgcpFixParamFit",
                     logger = ss_userAction.Log,
                     filenameroot = "ggplotLgcpFixParamFit",
                     aspectratio  = 1,
                     downloadfxns = list(png = ggplotLgcpFixParamFit,
                                         csv = dataggplotLgcpFixParamFit,
                                         txt = dataggplotLgcpFixParamFit),
                     visibleplot  = ggplotLgcpFixParamFit)
    
    
    dataggplotLgcpHyperParamFit <- function(){
      DF <- as.data.frame(LgcpModelFit()$DFPostHyperpar$DFPostHyperpar)
      return(DF)
    }
    
    ggplotLgcpHyperParamFit <- function(){
      DF <- LgcpModelFit()$DFPostHyperpar$DFPostHyperpar
      gl <- c()
      for(i in 1:length(DF)){
        nm <- strsplit(names(DF)[i], " ")[[1]]
        if(nm[1]=="Precision"){
          namesDFold <- names(DF)[i]
          names(DF)[i] <- paste("Stdev.", paste(nm[-1], collapse=" "), collapse=" ")
          assign(paste0("g",i), try(
            ggplot(data=data.frame(x=inla.tmarginal(function(x) 1/x**0.5, DF[[i]])[,1],
                                   y=inla.tmarginal(function(x) 1/x**0.5, DF[[i]])[,2]), aes(x=x,y=y)) + geom_line() +
              theme_bw()+ xlab(names(DF)[i]) +
              ylab(HTML(paste("Density f(",names(DF)[i],")"))) + ggtitle(names(DF)[i]) +
              theme(plot.title=element_text(color = "black", size = 16, face = "bold", hjust = 0.5)), silent=TRUE)
          )
          if(sum(class(eval(parse(text=paste0("g",i))))=="try-error")==1){
            assign(paste0("g",i),
                   ggplot(data=data.frame(x=inla.smarginal(marginal=DF[[i]])[[1]],
                                          y=inla.smarginal(marginal=DF[[i]])[[2]]), aes(x=x,y=y)) + geom_line() +
                     theme_bw()+ xlab(names(DF)[i]) +
                     ylab(HTML(paste("Density f(",namesDFold,")"))) + ggtitle(namesDFold) +
                     theme(plot.title=element_text(color = "black", size = 16, face = "bold", hjust = 0.5))
            )
          }
        } else{
          assign(paste0("g",i),
                 ggplot(data=data.frame(x=DF[[i]][,1], y=DF[[i]][,2]), aes(x=x,y=y)) + geom_line() +
                   theme_bw()+ xlab(names(DF)[i]) + xlim(quantile(DF[[i]][,1], probs = c(0.025,0.975))) +
                   ylab(HTML(paste("Density f(",names(DF)[i],")"))) + ggtitle(names(DF)[i]) +
                   theme(plot.title=element_text(color = "black", size = 16, face = "bold", hjust = 0.5))
          )
        }
        gl[i] <- paste0("g",i)
      }
      gt <- eval(parse(text=paste(gl, collapse="+")))
      return(gt)
    }
    
    downloadFile(
      id = "ggplotLgcpHyperParamFit",
      logger = ss_userAction.Log,
      filenameroot = "ggplotLgcpHyperParamFit",
      aspectratio  = 1,
      downloadfxns = list(png  = ggplotLgcpHyperParamFit,
                          csv = dataggplotLgcpHyperParamFit,
                          txt = dataggplotLgcpHyperParamFit)
    )
    
    downloadablePlot(id = "ggplotLgcpHyperParamFit",
                     logger = ss_userAction.Log,
                     filenameroot = "ggplotLgcpHyperParamFit",
                     aspectratio  = 1,
                     downloadfxns = list(png=ggplotLgcpHyperParamFit,
                                         csv = dataggplotLgcpHyperParamFit,
                                         txt = dataggplotLgcpHyperParamFit),
                     visibleplot  = ggplotLgcpHyperParamFit)
    
    
    tableLgcpModelFixedPar <- function(){
      DF <- LgcpModelFit()$Summary.fixed$Summary.fixed %>%
        mutate(across(where(is.numeric), round, digits = 2))
    }
    
    downloadableTable("tableLgcpModelFixedPar",
                      logger=ss_userAction.Log,
                      filenameroot="tableLgcpModelFixedPar",
                      downloaddatafxns=list(csv=tableLgcpModelFixedPar,
                                            tsv=tableLgcpModelFixedPar),
                      tabledata=tableLgcpModelFixedPar, rownames = TRUE,
                      caption="Summary fixed parameters")
    
    tableLgcpModelHyperPar <- function(){
      DF <- LgcpModelFit()$Summary.hyperpar$Summary.hyperpar %>%
        mutate(across(where(is.numeric), round, digits = 2))
    }
    
    downloadableTable("tableLgcpModelHyperPar",
                      logger=ss_userAction.Log,
                      filenameroot="tableLgcpModelHyperPar",
                      downloaddatafxns=list(csv=tableLgcpModelHyperPar,
                                            tsv=tableLgcpModelHyperPar),
                      tabledata=tableLgcpModelHyperPar, rownames = TRUE,
                      caption="Summary hyperparameters")
    
    tableLgcpModelInternalHyperPar <- function(){
      DF <- LgcpModelFit()$SummaryInternalHyper$SummaryInternalHyper %>%
        mutate(across(where(is.numeric), round, digits = 2))
    }
    
    downloadableTable("tableLgcpModelInternalHyperPar",
                      logger=ss_userAction.Log,
                      filenameroot="tableLgcpModelInternalHyperPar",
                      downloaddatafxns=list(csv=tableLgcpModelInternalHyperPar,
                                            tsv=tableLgcpModelInternalHyperPar),
                      tabledata=tableLgcpModelInternalHyperPar, rownames = TRUE,
                      caption="Summary internal hyperparameters")
    
    dataLgcpDICtable <- function(){
      DF <- LgcpModelFit()$DICmodel$DICmodel %>%
        mutate(across(where(is.numeric), round, digits = 2))
      return(DF)
    }
    
    downloadableTable("dataLgcpDICtable",
                      logger=ss_userAction.Log,
                      filenameroot="dataLgcpDICtable",
                      downloaddatafxns=list(csv=dataLgcpDICtable,
                                            tsv=dataLgcpDICtable),
                      tabledata=dataLgcpDICtable, rownames = TRUE,
                      caption="Model DIC")
    
    dataLgcpWAICtable <- function(){
      DF <- LgcpModelFit()$WAICmodel$WAICmodel %>%
        mutate(across(where(is.numeric), round, digits = 2))
      return(DF)
    }
    
    downloadableTable("dataLgcpWAICtable",
                      logger=ss_userAction.Log,
                      filenameroot="dataLgcpWAICtable",
                      downloaddatafxns=list(csv=dataLgcpWAICtable,
                                            tsv=dataLgcpWAICtable),
                      tabledata=dataLgcpWAICtable, rownames = TRUE,
                      caption="Model WAIC")
    
    dataLgcpCPOtable <- function(){
      CPO <- LgcpModelFit()$SummaryCPO$SummaryCPO
      DF <- data.frame(n=length(CPO), mean=mean(CPO), median=median(CPO), stdev=sd(CPO),
                       quantile2.5=quantile(CPO,probs=c(0.025,0.975))[1],
                       quantile97.5=quantile(CPO,probs=c(0.025,0.975))[2]) %>%
        mutate(across(where(is.numeric), round, digits = 2))
      return(DF)
    }
    
    downloadableTable("dataLgcpCPOtable",
                      logger=ss_userAction.Log,
                      filenameroot="dataLgcpCPOtable",
                      downloaddatafxns=list(csv=dataLgcpCPOtable,
                                            tsv=dataLgcpCPOtable),
                      tabledata=dataLgcpCPOtable, rownames = FALSE,
                      caption="Summary CPO")
    
    # Server_PreferentialModelling_Section ----
    
    PrefCheckBoxNames <- function(){
      if(input$PrefDataSimulatedLoaded=="load"){
        DF <- as.data.frame(datareadSample())
        if(input$SelectPrefFamily=="binomial"){DFnames <- names(DF)[c(5:ncol(DF))]}
        else{DFnames <- names(DF)[c(4:ncol(DF))]}
      } else if(input$PrefDataSimulatedLoaded=="sim"){
        DF <- Pref.sampling()
        DFnames <- names(DF)[c(4)]
      }
      return(DFnames)
    }
    
    observe({
      output$checkBoxPrefDataFrame <- renderUI({
        tagList(
          checkboxGroupInput(inputId="UserComponentsPref",
                             label="User Defined Components",
                             choices=PrefCheckBoxNames(),
                             selected=c()
          )
        )
      })
    })
    
    observe({
      if(input$PrefDataSimulatedLoaded=="sim"){DF <- Pref.sampling()}
      else{DF <- as.data.frame(datareadSample())}
      output$checkBoxPrefSharing <- renderUI({
          box(id="PrefCovariateSharingTerms", width=12, title="Sharing Components",
                checkboxGroupInput(inputId=paste0("UserComponentsPrefSharing"),
                                   label=paste("User Defined Sharing Components:"),
                                   choices=c(input$DefaultComponentsPref,input$UserComponentsPref),
                                   selected=c(input$DefaultComponentsPref,input$UserComponentsPref)
                                   )
              )
      })
    })
    
    observe({
      if(input$PrefDataSimulatedLoaded=="load"){
        output$SelectPrefEffectCov <- renderUI({

          if(length(input$UserComponentsPref)>0){
            box(id="PrefCovariateEffects", width=12, title="Covariate Effects",
                lapply(seq_along(input$UserComponentsPref), function(i){
                  list(
                    selectInput(inputId=paste0("PrefEffectCov",i), label=tags$span(style="color: blue; font-weight: bold;", paste(input$UserComponentsPref[i],"Effect")) ,
                                choices = list("Linear (or reference level)" = "linear", "Random Walk 1" = "rw1", "Random Walk 2" = "rw2", "SPDE 1" = "spde1", "IID"="iid"),
                                selected = "linear"),
                    conditionalPanel(condition=paste0("input.PrefEffectCov",i,"=='rw1'||","input.PrefEffectCov",i,"=='rw2'||","input.PrefEffectCov",i,"=='spde1'"),
                                     numericInput(inputId=paste0("PrefEffectCovNodes",i),
                                                  label="Number of nodes",
                                                  value=10, min=1, step=1
                                     )
                    ),
                    radioGroupButtons(
                      inputId = paste0("PrefEffectCustomPrior",i),
                      label = "Custom Prior",
                      choices = list("Auto" = "auto", "Custom" = "custom"),
                      status = "primary"
                    ),
                    conditionalPanel(condition=paste0("input.PrefEffectCustomPrior",i,"=='custom'"),
                                     list(
                                       selectInput(inputId = paste0("PrefEffectCovKindPrior",i),
                                                   label = "Prior distribution",
                                                   choices = list("Base" = "base", "PC prior" = "pc", "Uniform" = "unif", "Flat Uniform" = "flatunif")),
                                       conditionalPanel(condition=paste0("input.PrefEffectCovKindPrior",i,"=='base'"),
                                                        textInput(inputId = paste0("PrefEffectCovPriorBaseValues",i),
                                                                  label = "Base Prior Values",
                                                                  value = "0, 1e-3")),
                                       conditionalPanel(condition=paste0("input.PrefEffectCovKindPrior",i,"=='pc'"),
                                                        textInput(inputId = paste0("PrefEffectCovPriorPCValues",i),
                                                                  label = "PC-Prior Values",
                                                                  value = "0, 1e-3")),
                                       conditionalPanel(condition=paste0("input.PrefEffectCovKindPrior",i,"=='unif'"),
                                                        textInput(inputId = paste0("PrefEffectCovPriorUnif",i),
                                                                  label = "Uniform Lower and upper values",
                                                                  value = "0, 10",
                                                                  placeholder = "Lower and upper values: 'U(a,b)'"))
                                     )
                    ),
                    conditionalPanel(condition=paste0("input.PrefEffectCov",i,"=='rw1'||",
                                                      "input.PrefEffectCov",i,"=='rw2'||",
                                                      "input.PrefEffectCov",i,"=='iid'"),
                                     radioGroupButtons(
                                       inputId = paste0("PrefEffectCovConstr",i),
                                       label = "Sum to zero",
                                       choices = list("TRUE" = "TRUE", "FALSE" = "FALSE"),
                                       status = "primary"
                                     ))
                  )
                }))
          } else{}
        }) } else if(input$PrefDataSimulatedLoaded=="sim"){
          output$SelectPrefEffectCov <- renderUI({
            if(length(input$UserComponentsPref)>0){
              box(id="PrefCovariateEffects", width=12, title="Covariate Effects",
                  lapply(seq_along(input$UserComponentsPref), function(i){
                    list(
                      selectInput(inputId=paste0("PrefEffectCov",i), label=tags$span(style="color: blue; font-weight: bold;", paste(input$UserComponentsPref[i],"Effect")) ,
                                  choices = list("Linear (or reference level)" = "linear", "Random Walk 1" = "rw1", "Random Walk 2" = "rw2", "SPDE 1" = "spde1", "IID"="iid"),
                                  selected = "linear"),
                      conditionalPanel(condition=paste0("input.PrefEffectCov",i,"=='rw1'||","input.PrefEffectCov",i,"=='rw2'||","input.PrefEffectCov",i,"=='spde1'"),
                                       numericInput(inputId=paste0("PrefEffectCovNodes",i),
                                                    label="Number of nodes",
                                                    value=10, min=1, step=1
                                       )
                      ),
                      radioGroupButtons(
                        inputId = paste0("PrefEffectCustomPrior",i),
                        label = "Custom Prior",
                        choices = list("Auto" = "auto", "Custom" = "custom"),
                        status = "primary"
                      ),
                      conditionalPanel(condition=paste0("input.PrefEffectCustomPrior",i,"=='custom'"),
                                       list(
                                         selectInput(inputId = paste0("PrefEffectCovKindPrior",i),
                                                     label = "Prior distribution",
                                                     choices = list("Base" = "base", "PC prior" = "pc", "Uniform" = "unif", "Flat Uniform" = "unifflat")),
                                         conditionalPanel(condition=paste0("input.PrefEffectCovKindPrior",i,"=='base'"),
                                                          textInput(inputId = paste0("PrefEffectCovPriorBaseValues",i),
                                                                    label = "Base Prior Values",
                                                                    value = "0, 1e-3")),
                                         conditionalPanel(condition=paste0("input.PrefEffectCovKindPrior",i,"=='pc'"),
                                                          textInput(inputId = paste0("PrefEffectCovPriorPCValues",i),
                                                                    label = "PC-Prior Values",
                                                                    value = "0, 1e-3")),
                                         conditionalPanel(condition=paste0("input.PrefEffectCovKindPrior",i,"=='unif'"),
                                                          textInput(inputId = paste0("PrefEffectCovPriorUnif",i),
                                                                    label = "Uniform Lower and upper values",
                                                                    value = "0, 10",
                                                                    placeholder = "Lower and upper values: 'U(a,b)'"))
                                       )
                      ),
                      conditionalPanel(condition=paste0("input.PrefEffectCov",i,"=='rw1'||",
                                                        "input.PrefEffectCov",i,"=='rw2'||",
                                                        "input.PrefEffectCov",i,"=='iid'"),
                                       radioGroupButtons(
                                         inputId = paste0("PrefEffectCovConstr",i),
                                         label = "Sum to zero",
                                         choices = list("TRUE" = "TRUE", "FALSE" = "FALSE"),
                                         status = "primary"
                                       ))
                    )
                  }))
            } else{}
          }) }
    })
    
    observe({
      if(input$PrefDataSimulatedLoaded=="load"){
        DF <- as.data.frame(datareadSample())
        if(length(input$UserComponentsPref)>0){
          PrefUserComponent <- input$UserComponentsPref
          DF2 <- select(DF, PrefUserComponent[!as.vector(unlist(lapply(X=select(DF,PrefUserComponent), FUN=is.numeric)))])
          output$SelectLoadPrefEffectCovFactorPred <- renderUI({
            if(ncol(DF2)>0){
              box(id="PrefPredFactorLevel", width=12, title="Prediction Factor Level",
                  lapply(seq_along(names(DF2)), function(i){
                    choices <- unique(DF2[[names(DF2)[i]]])
                    list(
                      radioGroupButtons(inputId = paste0("PrefKindPredictionFactorLevel",i), label = tags$span(style="color: blue; font-weight: bold;", paste(names(DF2)[i], "(prediction protocol)")),
                                        choices = c("Reference level" = "reference", "Nearest level" = "nearest"), status = "success", justified = TRUE),
                      conditionalPanel(
                        condition=paste0("input.PrefKindPredictionFactorLevel",i,"=='reference'"),
                        selectInput(inputId=paste0("PrefEffectCovFactorPred",i), label=paste(names(DF2)[i],"Reference Factor"),
                                    choices = choices, selected = choices[1])
                      )
                    )
                  }))
            }
          })
        } else{
          output$SelectLoadPrefEffectCovFactorPred <- renderUI({})
        }
      }
    })
    
    JointPriorPreviewPref <- eventReactive(input$Prefpreviewpriordistributions,{
      d <- 2; nu <- 1
      rhoinputPC <- as.numeric(unlist(strsplit(input$Prefrangepcpriorprev, ",")))
      sigmainputPC <- as.numeric(unlist(strsplit(input$Prefsigmapcpriorprev, ",")))
      rhoPC <- seq(rhoinputPC[1], rhoinputPC[2], length.out=rhoinputPC[3])
      sigmaPC <- seq(sigmainputPC[1], sigmainputPC[2], length.out=sigmainputPC[3])
      rhosigmaPC <- expand.grid(rho=rhoPC, sigma=sigmaPC)
      rho0PC <- rhoinputPC[4]; alpha1 <- rhoinputPC[5]
      sigma0PC <- sigmainputPC[4]; alpha2 <- sigmainputPC[5]

      lambdarhoPC <- -log(alpha1)*rho0PC**(d/2) #-(rho0/sqrt(8*nu))**(d/2)*log(alpha1)
      lambdasigmaPC <- -log(alpha2)/sigma0PC #-(sqrt(8*nu)/rhosigma$rho)**(-nu)*sqrt(gamma(nu)/(gamma(nu+d/2)*(4*pi)**(d/2)))*log(alpha2)/sigma0
      pirhosigmaPCM <- d/2*lambdarhoPC*rhosigmaPC$rho**(-1-d/2)*exp(-lambdarhoPC*rhosigmaPC$rho**(-d/2))*lambdasigmaPC*exp(-lambdasigmaPC*rhosigmaPC$sigma)

      probIntPC <- diff(range(rhoPC))*diff(range(sigmaPC))/length(pirhosigmaPCM)*sum(pirhosigmaPCM)

      rhoinputBase <- as.numeric(unlist(strsplit(input$Prefrangebasepriorprev, ",")))
      sigmainputBase <- as.numeric(unlist(strsplit(input$Prefsigmabasepriorprev, ",")))
      rho <- seq(rhoinputBase[1], rhoinputBase[2], length.out=rhoinputBase[3])
      sigma <- seq(sigmainputBase[1], sigmainputBase[2], length.out=sigmainputBase[3])
      rhosigma <- expand.grid(rho=rho, sigma=sigma)
      rho0 <- rhoinputBase[4]; sigma0 <- sigmainputBase[4]
      meanrho <- rhoinputBase[5]; meansigma <- sigmainputBase[5]
      sdrho <- rhoinputBase[6]; sdsigma <- sigmainputBase[6]
      pirho <- (sqrt(2*pi*sdrho**2)*rhosigma$rho)**(-1)*exp(-(log(rhosigma$rho/rho0)**2-2*log(rhosigma$rho/rho0)*meanrho+meanrho**2)/(2*sdrho**2))
      pisigma <- (sqrt(2*pi*sdsigma**2)*rhosigma$sigma)**(-1)*exp(-(log(rhosigma$sigma/sigma0)**2-2*log(rhosigma$sigma/sigma0)*meansigma+meansigma**2)/(2*sdsigma**2))
      pirhosigmaM <- pirho*pisigma

      probInt <- sum(pirhosigmaM)*diff(range(rhosigma$rho))*diff(range(rhosigma$sigma))/length(pirhosigmaM)

      ListPreview <- list(rhosigmaPC=rhosigmaPC, pirhosigmaPCM=pirhosigmaPCM, probIntPC=probIntPC,
                          rhosigma=rhosigma, pirhosigmaM=pirhosigmaM, probInt=probInt)

      return(ListPreview)
    })

    PrefPreviewJointPriorPlot <- function(){
      ggplotPCprior <- ggplot() +
        geom_tile(data=data.frame(rho=JointPriorPreviewPref()$rhosigmaPC$rho, sigma=JointPriorPreviewPref()$rhosigmaPC$sigma, pi=JointPriorPreviewPref()$pirhosigmaPCM),
                  mapping=aes(x=rho, y=sigma, fill=pi)) +
        labs(title=paste("Joint PC-Prior. Cumulative Prob.=", round(JointPriorPreviewPref()$probIntPC, digits=2)),
             x=expression(rho), y=expression(sigma)) +
        scale_fill_viridis_c(option="turbo") + theme_bw() + theme(plot.title=element_text(face='bold'))

      ggplotENpriorAnalytic <- ggplot() +
        geom_tile(data=data.frame(rho=JointPriorPreviewPref()$rhosigma$rho, sigma=JointPriorPreviewPref()$rhosigma$sigma, pi=JointPriorPreviewPref()$pirhosigmaM),
                  mapping=aes(x=rho, y=sigma, fill=pi)) +
        scale_fill_viridis_c(option="turbo") +
        labs(title=paste("Joint Base Prior. Cumulative Prob.=", round(JointPriorPreviewPref()$probInt, digits=2)),
             x=expression(rho), y=expression(sigma)) +
        theme_bw() + theme(plot.title=element_text(face='bold'))

      ggplotJointPriorPreview <- ggplotPCprior + ggplotENpriorAnalytic
      return(ggplotJointPriorPreview)
    }

    downloadFile(
      id = "PrefPreviewJointPriorPlot",
      logger = ss_userAction.Log,
      filenameroot = "PrefPreviewJointPriorPlot",
      aspectratio  = 1,
      downloadfxns = list(png  = PrefPreviewJointPriorPlot)
    )

    downloadablePlot(id = "PrefPreviewJointPriorPlot",
                     logger = ss_userAction.Log,
                     filenameroot = "PrefPreviewJointPriorPlot",
                     aspectratio  = 1,
                     downloadfxns = list(png=PrefPreviewJointPriorPlot),
                     visibleplot  = PrefPreviewJointPriorPlot)
    
    ## Mesh construction for PreferentialModel ====
    
    PrefMeshBase <- reactive({
      if(input$PrefDataSimulatedLoaded=="sim"){
        DFsample <- as.data.frame(Pref.sampling())
        convexhull <- chull(DFsample[,1:2])
        convexhull <- c(convexhull, convexhull[1])
        qloc <- quantile(as.vector(dist(DFsample[sample(1:nrow(DFsample),size=min(c(50,nrow(DFsample)))),1:2])),probs=c(0.03,0.3))
        mesh <- fm_mesh_2d_inla(loc.domain=DFsample[convexhull,1:2], cutoff = qloc[1]/2, offset=c(-0.1, -0.2),
                                max.edge=c(qloc[1], qloc[2]))
        sample <- DFsample
      } else if(input$PrefDataSimulatedLoaded=="load"){
        DFsample <- datareadSample()
        if(input$PrefRasterSPDE=="raster"){
          rasterSample <- datareadRaster()[sample(1:nrow(datareadRaster()), min(c(50,nrow(datareadRaster())))),1:2]
          qloc <- quantile(as.vector(dist(rasterSample)),probs=c(0.03,0.3))
          convexhull <- chull(datareadRaster()[,1:2])
          convexhull <- c(convexhull, convexhull[1])
          mesh <- fm_mesh_2d_inla(loc.domain=datareadRaster()[convexhull,1:2], cutoff = qloc[1]/2, offset=c(-0.1, -0.2),
                                  max.edge=c(qloc[1], qloc[2]))
          sample <- rasterSample
        } else if(input$PrefRasterSPDE=="solvecov"){
          convexhull <- chull(rbind(DFsample[,1:2],datareadRaster()[,1:2]))
          convexhull <- c(convexhull, convexhull[1])
          qloc <- quantile(as.vector(dist(DFsample[sample(1:nrow(DFsample),size=min(c(50,nrow(DFsample)))),1:2])),probs=c(0.03,0.3))
          mesh <- fm_mesh_2d_inla(loc.domain=rbind(DFsample[,1:2],datareadRaster()[,1:2])[convexhull,], cutoff = qloc[1]/2, offset=c(-0.1, -0.2),
                                  max.edge=c(qloc[1], qloc[2]))
          sample <- DFsample
        }
      }
      result <- list(mesh=mesh, qloc=qloc, Sample=sample)
      return(result)
    })
    
    PrefMesh <- eventReactive(input$buildPrefMesh, {
      if(input$PrefDataSimulatedLoaded=="sim"){
        DFsample <- as.data.frame(Pref.sampling())
      } else if(input$PrefDataSimulatedLoaded=="load"){
        DFsample <- as.data.frame(datareadSample())
        DFraster <- try(datareadRaster(), silent=TRUE)
        if(class(DFraster)!="try-error"){
          DFsample <- rbind(DFsample[,1:2], DFraster[,1:2])
        }
      }
      
      interiorNonConvex <- function(x, condition, convex, resolution, file.read){
        if(condition=="interiorPrefMeshnonconvex"){
          interior <- inla.nonconvex.hull(points=x, convex=convex, resolution=resolution)
        } else if(condition=="interiorPrefMeshcustomboundary"){
          ext <- tools::file_ext(file.read$datapath)
          validate(need(ext == c("csv","rds"), "Please upload a csv or rds file"))
          if(ext=="csv"){coords <- read.csv(file.read$datapath)}
          else coords <- readRDS(file.read$datapath)
          innerBorderMesh <- SpatialPolygons(Srl=list(Polygons(srl=list(Polygon(coords=coords)), ID="interiorBoundary")))
          interior <- inla.sp2segment(sp=innerBorderMesh)
        } else{interior <- NULL}
        return(interior)
      }
      
      boundaryNonConvex <- function(x, condition, convex, resolution, file.read){
        if(condition=="PrefMeshnonconvex"){
          boundary <- inla.nonconvex.hull(points=x, convex=convex, resolution=resolution)
        } else if(condition=="PrefMeshcustomboundary"){
          ext <- tools::file_ext(file.read$datapath)
          validate(need(ext == c("csv","rds"), "Please upload a csv or rds file"))
          if(ext=="csv"){coords <- read.csv(file.read$datapath)}
          else coords <- readRDS(file.read$datapath)
          outerBorderMesh <- SpatialPolygons(Srl=list(Polygons(srl=list(Polygon(coords=coords)), ID="externalBoundary")))
          boundary <- inla.sp2segment(sp=outerBorderMesh)
        } else{boundary <- NULL}
        return(boundary)
      }
      
      if(input$selectionPrefMesh=="qlocation"){
        if(input$PrefRasterSPDE=="raster"){
          rasterSample <- datareadRaster()[sample(1:nrow(datareadRaster()), min(c(50,nrow(datareadRaster())))),1:2]
          qloc <- quantile(as.vector(dist(rasterSample)),probs=c(ifelse(input$PrefCustomMesh,input$PrefMeshQloc,0.03),0.3))
          convexhull <- chull(DFsample[,1:2])
          convexhull <- c(convexhull, convexhull[1])
          mesh <- fm_mesh_2d_inla(loc.domain=DFsample[convexhull,1:2], cutoff = qloc[1]/2, offset=c(-0.1, -0.2), max.edge=c(qloc[1], qloc[2]),
                                  boundary=list(interiorNonConvex(x=as.matrix(rasterSample[,1:2]), condition=input$interiorPrefMesh, convex=input$interiorcurvaturePrefMesh, resolution=input$interiorresolutionPrefMesh, file.read=input$interiorshapefilePrefMesh),
                                                boundaryNonConvex(x=as.matrix(rasterSample[,1:2]), condition=input$boundaryPrefMesh, convex=input$curvaturePrefMesh, resolution=input$resolutionPrefMesh, file.read=input$shapefilePrefMesh)))
          sample <- rasterSample
        } else if(input$PrefRasterSPDE=="solvecov"){
          qloc <- quantile(as.vector(dist(DFsample[sample(1:nrow(DFsample),size=min(c(50,nrow(DFsample)))),1:2])),probs=c(ifelse(input$PrefCustomMesh,input$PrefMeshQloc,0.03),0.3))
          convexhull <- chull(DFsample[,1:2])
          convexhull <- c(convexhull, convexhull[1])
          mesh <- fm_mesh_2d_inla(loc.domain=DFsample[convexhull,1:2], cutoff = qloc[1]/2, offset=c(-0.1, -0.2), max.edge=c(qloc[1], qloc[2]),
                                  boundary=list(interiorNonConvex(x=as.matrix(DFsample[,1:2]), condition=input$interiorPrefMesh, convex=input$interiorcurvaturePrefMesh, resolution=input$interiorresolutionPrefMesh, file.read=input$interiorshapefilePrefMesh),
                                                boundaryNonConvex(x=as.matrix(DFsample[,1:2]), condition=input$boundaryPrefMesh, convex=input$curvaturePrefMesh, resolution=input$resolutionPrefMesh, file.read=input$shapefilePrefMesh)))
          sample <- DFsample
        }
      } else if(input$selectionPrefMesh=="edgelength"){
        if(input$PrefRasterSPDE=="raster"){
          rasterSample <- datareadRaster()
          qloc <- input$EdgeLengthPrefMesh
          convexhull <- chull(DFsample[,1:2])
          convexhull <- c(convexhull, convexhull[1])
          mesh <- fm_mesh_2d_inla(loc.domain=DFsample[convexhull,1:2], max.edge=c(1,1.5)*input$EdgeLengthPrefMesh, cutoff=input$EdgeLengthPrefMesh/5, offset=c(-0.1, -0.2), max.edge=c(qloc[1], qloc[2]),
                                  boundary=list(interiorNonConvex(x=as.matrix(rasterSample[,1:2]), condition=input$interiorPrefMesh, convex=input$interiorcurvaturePrefMesh, resolution=input$interiorresolutionPrefMesh, file.read=input$interiorshapefilePrefMesh),
                                                boundaryNonConvex(x=as.matrix(rasterSample[,1:2]), condition=input$boundaryPrefMesh, convex=input$curvaturePrefMesh, resolution=input$resolutionPrefMesh, file.read=input$shapefilePrefMesh)))
          sample <- rasterSample
        } else if(input$PrefRasterSPDE=="solvecov"){
          qloc <- input$EdgeLengthPrefMesh
          convexhull <- chull(DFsample[,1:2])
          convexhull <- c(convexhull, convexhull[1])
          mesh <- fm_mesh_2d_inla(loc.domain=DFsample[convexhull,1:2], max.edge=c(1,1.5)*input$EdgeLengthPrefMesh,  cutoff=input$EdgeLengthPrefMesh/5, offset=c(-0.1, -0.2),
                                  boundary=list(interiorNonConvex(x=as.matrix(DFsample[,1:2]), condition=input$interiorPrefMesh, convex=input$interiorcurvaturePrefMesh, resolution=input$interiorresolutionPrefMesh, file.read=input$interiorshapefilePrefMesh),
                                                boundaryNonConvex(x=as.matrix(DFsample[,1:2]), condition=input$boundaryPrefMesh, convex=input$curvaturePrefMesh, resolution=input$resolutionPrefMesh, file.read=input$shapefilePrefMesh)))
          sample <- DFsample
        }
      }
      
      result <- list(mesh=mesh, qloc=qloc, Sample=sample)
      return(result)
    })

    ggplotPrefMesh <- function(){
      if(input$buildPrefMesh==0){
        ggplot()+ gg(PrefMeshBase()$mesh)+ theme_bw() + xlab("Latitude") + ylab("Longitude") +
          ggtitle("Mesh over the study region") +
          geom_point(data=PrefMeshBase()$Sample,
                     aes(x=PrefMeshBase()$Sample[,1],y=PrefMeshBase()$Sample[,2]), size=1) +
          theme(plot.title=element_text(color = "black", size = 16, face = "bold", hjust = 0.5))
      } else{
        ggplot()+ gg(PrefMesh()$mesh)+ theme_bw() + xlab("Latitude") + ylab("Longitude") +
          ggtitle("Mesh over the study region") +
          geom_point(data=PrefMesh()$Sample,
                     aes(x=PrefMesh()$Sample[,1],y=PrefMesh()$Sample[,2]), size=1) +
          theme(plot.title=element_text(color = "black", size = 16, face = "bold", hjust = 0.5))
      }

    }

    downloadFile(
      id = "ggplotPrefMesh",
      logger = ss_userAction.Log,
      filenameroot = "ggplotPrefMesh",
      aspectratio  = 1,
      downloadfxns = list(png  = ggplotPrefMesh)
    )

    downloadablePlot(
      id = "ggplotPrefMesh",
      logger = ss_userAction.Log,
      filenameroot = "ggplotPrefMesh",
      aspectratio  = 1,
      downloadfxns = list(png  = ggplotPrefMesh),
      visibleplot  = ggplotPrefMesh)
    
    ## Preferential modelling code ====
    
    PrefModelFit <- eventReactive(input$fitPref, {
      showNotification(ui=paste("Fitting the data."), duration = NULL)
      t1 <- Sys.time()
      #taking the data from simulation or from the loading tab
      if(input$PrefDataSimulatedLoaded=="sim"){
        DFsample <- as.data.frame(Pref.sampling())
      } else if(input$PrefDataSimulatedLoaded=="load"){
        DFsample <- as.data.frame(datareadSample())
        if(input$PrefRasterSPDE=="raster"){
          DFraster <- as.data.frame(datareadRaster())
        }
      }

      variablesChosenDefault <- input$DefaultComponentsPref
      variablesChosenUser <- input$UserComponentsPref
      variablesChosen <- c(variablesChosenDefault, variablesChosenUser)

      UserComponentsPrefSharing <- list()
      UserComponentsPrefSharing[[paste0("Pref")]] <- eval(parse(text=paste0("input$UserComponentsPrefSharing")))

      if(input$buildPrefMesh==0){
        server_mesh <- PrefMeshBase()
        mesh <- server_mesh$mesh
      } else{
        server_mesh <- PrefMesh()
        mesh <- server_mesh$mesh
      }

      if(input$optionPrefS=="auto"){
        if(input$KindPriorSpatialEffectPref=="PC.prior"){
          prior.range <- c(mean(c(diff(range(mesh$loc[mesh$segm$int$idx[,1],1])),diff(range(mesh$loc[mesh$segm$int$idx[,1],2]))))/5, 0.5)
          prior.sigma <- c(1,0.5)
          spde <- inla.spde2.pcmatern(mesh, prior.range = prior.range, prior.sigma = prior.sigma, alpha=2)
        } else{
          prior.range <- mean(c(diff(range(mesh$loc[mesh$segm$int$idx[,1],1])),diff(range(mesh$loc[mesh$segm$int$idx[,1],2]))))/5
          prior.sigma <- 1
          alpha <- 2; d <- 2
          nu <-  alpha - d/2
          kappa0 <-  log(8*nu)/2 -log(prior.range)
          tau0 <-  0.5*(lgamma(nu) - lgamma(nu + d/2) - d/2*log(4*pi)) - nu*kappa0 - log(prior.sigma)
          spde <-  inla.spde2.matern(mesh = mesh, B.tau = cbind(tau0, nu, -1), B.kappa = cbind(kappa0, -1, 0),
                                     theta.prior.mean = c(0.5,0.5), theta.prior.prec = c(1,1))
        }
      } else if(input$optionPrefS=="custom"){
        if(input$KindPriorSpatialEffectPref=="PC.prior"){
          prior.range <- as.numeric(unlist(strsplit(input$PrefPriorRangePC, split=",")))
          prior.sigma <- as.numeric(unlist(strsplit(input$PrefPriorStdevPC, split=",")))
          spde <- inla.spde2.pcmatern(mesh, prior.range = prior.range, prior.sigma = prior.sigma, alpha=2)
        } else{
          prior.range <- as.numeric(unlist(strsplit(input$PrefPriorRangeBase, split=",")))
          prior.sigma <- as.numeric(unlist(strsplit(input$PrefPriorStdevBase, split=",")))
          alpha <- 2; d <- 2
          nu <-  alpha - d/2
          kappa0 <-  log(8*nu)/2 -log(prior.range[1])
          tau0 <-  0.5*(lgamma(nu) - lgamma(nu + d/2) - d/2*log(4*pi)) - nu*kappa0 - log(prior.sigma[1])
          spde <-  inla.spde2.matern(mesh = mesh, B.tau = cbind(tau0, nu, -1), B.kappa = cbind(kappa0, -1, 0),
                                     theta.prior.mean = c(prior.range[2],prior.sigma[2]), theta.prior.prec = c(prior.range[3],prior.sigma[3]))
        }
      }

      spde.index <- inla.spde.make.index(name="Spatial", n.spde = spde$n.spde)

      # LGCP mesh operations
      ldomain <- unique(mesh$loc[mesh$segm$int$idx,1:2])
      dmesh <- mesh.dual(mesh = mesh)
      domain.polys <- Polygons(list(Polygon(ldomain)), '0')
      domainSP <- SpatialPolygons(list(domain.polys))
      w <- sapply(1:length(dmesh), function(i) {
        if (gIntersects(dmesh[i, ], domainSP))
          return(gArea(gIntersection(dmesh[i, ], domainSP)))
        else return(0)
      })

      n <- nrow(DFsample)
      nv <- mesh$n
      imat <- Diagonal(nv, rep(1, nv))
      lmat <- inla.spde.make.A(mesh, as.matrix(DFsample[,1:2]))
      A.inf <- rbind(imat, lmat)

      ### Prediction of covariates ====

      List.covariates.inf <- list()
      List.covariates.mesh <- list()
      List.covariates.pred <- list()


      if((input$PrefDataSimulatedLoaded=="load"&input$PrefRasterSPDE=="solvecov")|input$PrefDataSimulatedLoaded=="sim"|(input$PrefDataSimulatedLoaded=="load"&input$PrefRasterSPDE=="raster"|input$PrefRasterPred=="SPDEraster")){
        x.pred <- seq(range(mesh$loc[mesh$segm$int$idx[,1], 1])[1], range(mesh$loc[mesh$segm$int$idx[,1], 1])[2], length.out=input$PrefSPDErasterInput)
        y.pred <- seq(range(mesh$loc[mesh$segm$int$idx[,1], 2])[1], range(mesh$loc[mesh$segm$int$idx[,1], 2])[2], length.out=input$PrefSPDErasterInput)
        xy.pred <- expand.grid(x=x.pred,y=y.pred)
        xy.pred <- xy.pred[which(!is.na(over(SpatialPoints(coords=xy.pred),SpatialPolygons(Srl=list(Polygons(srl=list(Polygon(coords=mesh$loc[mesh$segm$int$idx[,1], 1:2])), ID="int")))))),1:2]
        A.geo.pred <- inla.spde.make.A(mesh=mesh, loc=as.matrix(xy.pred))
      } else{ #if(input$PrefDataSimulatedLoaded=="load"&input$PrefRasterSPDE=="raster"|input$PrefRasterPred=="rasterpred"){
        xy.pred <- DFraster
        A.geo.pred <- inla.spde.make.A(mesh=mesh, loc=as.matrix(xy.pred))
      }

      for(i in seq_along(variablesChosenUser)){
        if(!is.character(DFsample[,variablesChosenUser[i]])){
          prior.range.cov <- c(mean(c(diff(range(mesh$loc[mesh$segm$int$idx[,1],1])),diff(range(mesh$loc[mesh$segm$int$idx[,1],2]))))/5, 0.5)
          prior.sigma.cov <- c(1,0.5)
          spde.cov <- inla.spde2.pcmatern(mesh, prior.range = prior.range.cov, prior.sigma = prior.sigma.cov, alpha=2)
          spde.cov.index <- inla.spde.make.index(name="spatial.cov", n.spde = spde.cov$n.spde)
          formula.cov <- y ~ -1 + Intercept + f(spatial.cov, model=spde.cov)

          if((input$PrefDataSimulatedLoaded=="load"&input$PrefRasterSPDE=="solvecov")|input$PrefDataSimulatedLoaded=="sim"){
            x.pred <- seq(range(mesh$loc[mesh$segm$int$idx[,1], 1])[1], range(mesh$loc[mesh$segm$int$idx[,1], 1])[2], length.out=input$PrefSPDErasterInput)
            y.pred <- seq(range(mesh$loc[mesh$segm$int$idx[,1], 2])[1], range(mesh$loc[mesh$segm$int$idx[,1], 2])[2], length.out=input$PrefSPDErasterInput)
            xy.pred <- expand.grid(x=x.pred,y=y.pred)
            xy.pred <- xy.pred[which(!is.na(over(SpatialPoints(coords=xy.pred),SpatialPolygons(Srl=list(Polygons(srl=list(Polygon(coords=mesh$loc[mesh$segm$int$idx[,1], 1:2])), ID="int")))))),1:2]
            A.geo.pred <- inla.spde.make.A(mesh=mesh, loc=as.matrix(xy.pred))
            Inf.stack.cov <- inla.stack(data=list(y=DFsample[,variablesChosenUser[i]]),
                                        A=list(lmat, 1),
                                        effects=list(list(spatial.cov=spde.cov.index$spatial.cov),
                                                     list(Intercept=rep(1,nrow(DFsample)))
                                        ),
                                        tag="Inference.cov")
            Pred.mesh.stack.cov <- inla.stack(data=list(y=rep(NA, mesh$n)),
                                              A=list(diag(1,nrow=mesh$n), 1),
                                              effects=list(list(spatial.cov=spde.cov.index$spatial.cov),
                                                           list(Intercept=rep(1,mesh$n))
                                              ),
                                              tag="Mesh.cov")
            Pred.stack.cov <- inla.stack(data=list(y=rep(NA,nrow(xy.pred))),
                                         A=list(A.geo.pred, 1),
                                         effects=list(list(spatial.cov=spde.cov.index$spatial.cov),
                                                      list(Intercept=rep(1,nrow(xy.pred)))
                                         ),
                                         tag="Prediction.cov")
            Total.stack.cov <- inla.stack(Inf.stack.cov, Pred.mesh.stack.cov, Pred.stack.cov)
            mod.cov <- inla(formula=formula.cov, data=inla.stack.data(Total.stack.cov), family="gaussian",
                            control.predictor = list(compute=FALSE, A=inla.stack.A(Total.stack.cov)))
            indx.inf <- inla.stack.index(Total.stack.cov, tag="Inference.cov")$data
            indx.mesh <- inla.stack.index(Total.stack.cov, tag="Mesh.cov")$data
            indx.pred <- inla.stack.index(Total.stack.cov, tag="Prediction.cov")$data
            List.covariates.inf[[variablesChosenUser[i]]] <- as.numeric(!is.na(DFsample[,variablesChosenUser[i]]))*DFsample[,variablesChosenUser[i]] + as.numeric(is.na(DFsample[,variablesChosenUser[i]]))*mod.cov$summary.fitted.values[indx.inf,"mean"]
            List.covariates.mesh[[variablesChosenUser[i]]] <- mod.cov$summary.fitted.values[indx.mesh,"mean"]
            List.covariates.pred[[variablesChosenUser[i]]] <- mod.cov$summary.fitted.values[indx.pred,"mean"]
          } else if(input$PrefDataSimulatedLoaded=="load"&input$PrefRasterSPDE=="raster"|input$PrefRasterPred=="SPDEraster"){
            x.pred <- seq(range(mesh$loc[mesh$segm$int$idx[,1], 1])[1], range(mesh$loc[mesh$segm$int$idx[,1], 1])[2], length.out=input$PrefSPDErasterInput)
            y.pred <- seq(range(mesh$loc[mesh$segm$int$idx[,1], 2])[1], range(mesh$loc[mesh$segm$int$idx[,1], 2])[2], length.out=input$PrefSPDErasterInput)
            xy.pred <- expand.grid(x=x.pred,y=y.pred)
            xy.pred <- xy.pred[which(!is.na(over(SpatialPoints(coords=xy.pred),SpatialPolygons(Srl=list(Polygons(srl=list(Polygon(coords=mesh$loc[mesh$segm$int$idx[,1], 1:2])), ID="int")))))),1:2]
            A.geo.pred <- inla.spde.make.A(mesh=mesh, loc=as.matrix(xy.pred))
            if(variablesChosenUser[i] %in% colnames(DFraster)){
              Inf.stack.cov <- inla.stack(data=list(y=c(DFsample[,variablesChosenUser[i]], DFraster[!is.na(DFraster[,variablesChosenUser[i]]),variablesChosenUser[i]])),
                                          A=list(inla.spde.make.A(mesh=mesh, loc=as.matrix(rbind(DFsample[,1:2],DFraster[!is.na(DFraster[,variablesChosenUser[i]]),1:2]))), 1),
                                          effects=list(list(spatial.cov=spde.cov.index$spatial.cov),
                                                       list(Intercept=rep(1,nrow(DFsample)+nrow(DFraster)))
                                          ),
                                          tag="Inference.cov")
            } else{
              Inf.stack.cov <- inla.stack(data=list(y=DFsample[,variablesChosenUser[i]]),
                                          A=list(inla.spde.make.A(mesh=mesh, loc=as.matrix(DFsample[,1:2])), 1),
                                          effects=list(list(spatial.cov=spde.cov.index$spatial.cov),
                                                       list(Intercept=rep(1,nrow(DFsample)))
                                          ),
                                          tag="Inference.cov")
            }
            Pred.mesh.stack.cov <- inla.stack(data=list(y=rep(NA, mesh$n)),
                                              A=list(diag(1,nrow=mesh$n), 1),
                                              effects=list(list(spatial.cov=spde.cov.index$spatial.cov),
                                                           list(Intercept=rep(1,mesh$n))
                                              ),
                                              tag="Mesh.cov")
            Pred.stack.cov <- inla.stack(data=list(y=rep(NA,nrow(xy.pred))),
                                         A=list(A.geo.pred, 1),
                                         effects=list(list(spatial.cov=spde.cov.index$spatial.cov),
                                                      list(Intercept=rep(1,nrow(xy.pred)))
                                         ),
                                         tag="Prediction.cov")
            Total.stack.cov <- inla.stack(Inf.stack.cov, Pred.mesh.stack.cov, Pred.stack.cov)
            mod.cov <- inla(formula=formula.cov, data=inla.stack.data(Total.stack.cov), family="gaussian",
                            control.predictor = list(compute=FALSE, A=inla.stack.A(Total.stack.cov)))
            indx.inf <- inla.stack.index(Total.stack.cov, tag="Inference.cov")$data
            indx.mesh <- inla.stack.index(Total.stack.cov, tag="Mesh.cov")$data
            indx.pred <- inla.stack.index(Total.stack.cov, tag="Prediction.cov")$data
            List.covariates.inf[[variablesChosenUser[i]]] <- as.numeric(!is.na(DFsample[,variablesChosenUser[i]]))*DFsample[,variablesChosenUser[i]] + (as.numeric(is.na(DFsample[,variablesChosenUser[i]]))*(mod.cov$summary.fitted.values[indx.inf,"mean"])[1:nrow(DFsample)])
            List.covariates.mesh[[variablesChosenUser[i]]] <- mod.cov$summary.fitted.values[indx.mesh,"mean"]
            List.covariates.pred[[variablesChosenUser[i]]] <- mod.cov$summary.fitted.values[indx.pred,"mean"]
          } else if(input$PrefDataSimulatedLoaded=="load"&input$PrefRasterSPDE=="raster"|input$PrefRasterPred=="rasterpred"){
            xy.pred <- DFraste[,1:2]
            A.geo.pred <- inla.spde.make.A(mesh=mesh, loc=as.matrix(xy.pred))

            Inf.stack.cov <- inla.stack(data=list(y=c(DFsample[,variablesChosenUser[i]], DFraster[,variablesChosenUser[i]])),
                                        A=list(inla.spde.make.A(mesh=mesh, loc=as.matrix(rbind(DFsample[,1:2],DFraster[,1:2]))), 1),
                                        effects=list(list(spatial.cov=spde.cov.index$spatial.cov),
                                                     list(Intercept=rep(1,nrow(DFsample)+nrow(DFraster)))
                                        ),
                                        tag="Inference.cov")
            Pred.mesh.stack.cov <- inla.stack(data=list(y=rep(NA, mesh$n)),
                                              A=list(diag(1,nrow=mesh$n), 1),
                                              effects=list(list(spatial.cov=spde.cov.index$spatial.cov),
                                                           list(Intercept=rep(1,mesh$n))
                                              ),
                                              tag="Mesh.cov")
            Total.stack.cov <- inla.stack(Inf.stack.cov, Pred.mesh.stack.cov)
            mod.cov <- inla(formula=formula.cov, data=inla.stack.data(Total.stack.cov), family="gaussian",
                            control.predictor = list(compute=FALSE, A=inla.stack.A(Total.stack.cov)))
            indx.mesh <- inla.stack.index(Total.stack.cov, tag="Mesh.cov")$data
            List.covariates.mesh[[variablesChosenUser[i]]] <- mod.cov$summary.fitted.values[indx.mesh,"mean"]

            List.covariates.inf[[variablesChosenUser[i]]] <- DFsample[,variablesChosenUser[i]]
            List.covariates.pred[[variablesChosenUser[i]]] <- DFraster[,variablesChosenUser[i]]
          }
        } else{
          if((input$PrefDataSimulatedLoaded=="load"&input$PrefRasterSPDE=="solvecov")|input$PrefDataSimulatedLoaded=="sim"){
            cov.inf.indx <- as.vector(unlist(lapply(X=1:nrow(DFsample), FUN=function(Y){which.min(apply(X=(DFsample[!is.na(DFsample[,variablesChosenUser[i]]),1:2]-matrix(DFsample[Y,1:2],ncol=2))**2, MARGIN=1, FUN=sum))})))
            List.covariates.inf[[variablesChosenUser[i]]] <- DFsample[cov.inf.indx, variablesChosenUser[i]]
            mesh.indx <- as.vector(unlist(lapply(X=1:nrow(mesh$loc), FUN=function(Y){which.min(apply(X=(DFsample[!is.na(DFsample[,variablesChosenUser[i]]),1:2]-matrix(mesh$loc[Y,1:2],ncol=2))**2, MARGIN=1, FUN=sum))})))
            List.covariates.mesh[[variablesChosenUser[i]]] <- DFsample[mesh.indx, variablesChosenUser[i]]

            PrefFactorModelDF <- data.frame(y=1, Pref=c(List.covariates.mesh[[variablesChosenUser[i]]], List.covariates.inf[[variablesChosenUser[i]]]))
            colnames(PrefFactorModelDF) <- c("y", variablesChosenUser[i])
            idx.factor <- which(variablesChosenUser[i]==names(DFsample[!as.vector(unlist(lapply(X=DFsample, FUN=is.numeric)))])[names(DFsample[!as.vector(unlist(lapply(X=DFsample, FUN=is.numeric)))])%in%c(variablesChosenUser)])

            if(eval(parse(text=paste0("input$PrefKindPredictionFactorLevel",idx.factor)))=="nearest"){
              cov.pred.ind <- as.vector(unlist(lapply(X=1:nrow(xy.pred), FUN=function(Y){which.min(apply(X=(DFsample[!is.na(DFsample[,variablesChosenUser[i]]),1:2]-matrix(xy.pred[Y,1:2],ncol=2))**2, MARGIN=1, FUN=sum))})))
              List.covariates.pred[[variablesChosenUser[i]]] <- DFsample[cov.pred.ind, variablesChosenUser[i]]
            } else{
              List.covariates.pred[[variablesChosenUser[i]]] <- rep(NA, times=nrow(xy.pred))
            }

          } else {
            DFsampleraster <- rbind(DFsample[,c(1:2,which(variablesChosenUser[i]==colnames(DFsample)))], DFraster[,c(1:2, which(variablesChosenUser[i]==colnames(DFraster)))])
            cov.inf.indx <- as.vector(unlist(lapply(X=1:nrow(DFsample), FUN=function(Y){which.min(apply(X=(DFsampleraster[!is.na(DFsampleraster[,variablesChosenUser[i]]),1:2]-matrix(DFsample[Y,1:2],ncol=2))**2, MARGIN=1, FUN=sum))})))
            List.covariates.inf[[variablesChosenUser[i]]] <- DFsample[cov.inf.indx, variablesChosenUser[i]]
            mesh.indx <- as.vector(unlist(lapply(X=1:nrow(mesh$loc), FUN=function(Y){which.min(apply(X=(DFsampleraster[!is.na(DFsampleraster[,variablesChosenUser[i]]),1:2]-matrix(mesh$loc[Y,1:2],ncol=2))**2, MARGIN=1, FUN=sum))})))
            List.covariates.mesh[[variablesChosenUser[i]]] <- DFsample[mesh.indx, variablesChosenUser[i]]

            PrefFactorModelDF <- data.frame(y=1, Pref=c(List.covariates.mesh[[variablesChosenUser[i]]], List.covariates.inf[[variablesChosenUser[i]]]))
            colnames(PrefFactorModelDF) <- c("y", variablesChosenUser[i])
            idx.factor <- which(variablesChosenUser[i]==names(DFsample[!as.vector(unlist(lapply(X=DFsample, FUN=is.numeric)))])[names(DFsample[!as.vector(unlist(lapply(X=DFsample, FUN=is.numeric)))])%in%c(variablesChosenUser)])

            if(eval(parse(text=paste0("input$PrefKindPredictionFactorLevel",idx.factor)))=="nearest"){
              cov.pred.ind <- as.vector(unlist(lapply(X=1:nrow(xy.pred), FUN=function(Y){which.min(apply(X=(DFsampleraster[!is.na(DFsampleraster[,variablesChosenUser[i]]),1:2]-matrix(xy.pred[Y,1:2],ncol=2))**2, MARGIN=1, FUN=sum))})))
              List.covariates.pred[[variablesChosenUser[i]]] <- DFsample[cov.pred.ind, variablesChosenUser[i]]
            } else{
              List.covariates.pred[[variablesChosenUser[i]]] <- rep(NA, times=nrow(xy.pred))
            }

          }
        }
      }

      ### Building the main model and stacks structure ====

      Inf.geo.effects.list <- list(
        list(),
        list()
      )

      Pred.geo.effects.list <- list(
        list(),
        list()
      )


      Inf.Prefs.effects.list <- list()
      
      Inf.Prefs.effects.list[[paste0("Inf.", names(UserComponentsPrefSharing)[i],".effects.list")]] <- list(list(),list())

      formula_mod <- c("y ~ -1")

      A_Inf.spde1 <- list()
      A_Pred.spde1 <- list()

      for(i in seq_along(variablesChosen)){
        if(variablesChosen[i]=="Intercept"){
          formula_mod <- paste(formula_mod, "f(Intercept, model='linear')", sep=" + ")
          Inf.geo.effects.list[[2]][["Intercept"]] <- c(rep(1, times=nv), rep(1, times=n))
          Pred.geo.effects.list[[2]][["Intercept"]] <- rep(1, times=nrow(A.geo.pred))

          if(length(UserComponentsPrefSharing)>0){
            for(k in seq_along(UserComponentsPrefSharing)){
              if("Intercept" %in% UserComponentsPrefSharing[[k]]){
                showNotification(ui=paste0("From a theoretical perspective, sharing linear effects is feasible but doesn't make logical sense. Consequently, we will approach them as non-shared effects."), duration=10, closeButton=TRUE, type="warning")
                ComponentPrefGroup <- "Pref"
                formula_mod <- paste(formula_mod, paste0("f(Intercept_",ComponentPrefGroup[k],", model='linear')"), sep=" + ")
                Inf.Prefs.effects.list[[paste0("Inf.", names(UserComponentsPrefSharing)[k],".effects.list")]][[2]][[paste0("Intercept_", ComponentPrefGroup[k])]] <- c(rep(1, times=nv), rep(1, times=n))
              } else{
                ComponentPrefGroup <- "Pref"
                formula_mod <- paste(formula_mod, paste0("f(Intercept_",ComponentPrefGroup[k],", model='linear')"), sep=" + ")
                Inf.Prefs.effects.list[[paste0("Inf.", names(UserComponentsPrefSharing)[k],".effects.list")]][[2]][[paste0("Intercept_", ComponentPrefGroup[k])]] <- c(rep(1, times=nv), rep(1, times=n))
              }
            }
          }
        } else if(variablesChosen[i]=="Spatial Effect"){
          formula_mod <- paste(formula_mod, "f(Spatial, model=spde)", sep=" + ")
          Inf.geo.effects.list[[1]][["Spatial"]] <- spde.index$Spatial
          Pred.geo.effects.list[[1]][["Spatial"]] <- spde.index$Spatial

          if(length(UserComponentsPrefSharing)>0){
            for(k in seq_along(UserComponentsPrefSharing)){
              if("Spatial Effect" %in% UserComponentsPrefSharing[[k]]){
                ComponentPrefGroup <- "Pref"
                formula_mod <- paste(formula_mod, paste0("f(Spatial_",ComponentPrefGroup[k],"_copy, copy='Spatial', fixed=FALSE)"), sep=" + ")
                Inf.Prefs.effects.list[[paste0("Inf.", names(UserComponentsPrefSharing)[k],".effects.list")]][[1]][[paste0("Spatial_", ComponentPrefGroup[k], "_copy")]] <- spde.index$Spatial
              } else{
                ComponentPrefGroup <- "Pref"
                formula_mod <- paste(formula_mod, paste0("f(Spatial_",ComponentPrefGroup[k],")"), sep=" + ")
                Inf.Prefs.effects.list[[paste0("Inf.", names(UserComponentsPrefSharing)[k],".effects.list")]][[1]][[paste0("Spatial_", ComponentPrefGroup[k])]] <- spde.index$Spatial
              }
            }
          }
        } else{
          j <- which(variablesChosen[i]==variablesChosenUser)
          if(!is.character(DFsample[,variablesChosenUser[j]])){
            if(eval(parse(text=paste0("input$PrefEffectCov",j)))=="rw1"|eval(parse(text=paste0("input$PrefEffectCov",j)))=="rw2"){
              if(eval(parse(text=paste0("input$PrefEffectCustomPrior",j)))=="custom"){
                if(eval(parse(text=paste0("input$PrefEffectCovKindPrior",j)))=="pc"){
                  assign(paste0("pc.values", variablesChosenUser[j]), c( as.numeric(unlist(strsplit(eval(parse(text=paste0("input$PrefEffectCovPriorPCValues",j))), ","))) ))
                  formula_mod <- paste(formula_mod, paste0("f(",variablesChosenUser[j],", model='",eval(parse(text=paste0("input$PrefEffectCov",j))),"', hyper=list(prec=list(prior='pc.prec', param=", eval(parse(text=paste0("pc.values", variablesChosenUser[j]))), ")), constr=",eval(parse(text=paste0("input$PrefEffectCovConstr",j))),", scale.model=TRUE)"), sep=" + ")
                } else if(eval(parse(text=paste0("input$PrefEffectCovKindPrior",j)))=="unif"){
                  lim <- as.numeric(unlist(strsplit(eval(parse(text=paste0("input$PrefEffectCovPriorUnif",j))), ",")))
                  sigma <- seq(0, lim[2]*3, length.out=1E5)
                  theta <- -2*log(sigma)
                  logdens <- sapply(X=sigma, FUN=logdunif, lim=lim)
                  unif.prior <- list(theta=list(prior=paste0("table: ", paste(c(theta, logdens), collapse=" "))))
                  assign(paste0("hyper", variablesChosenUser[j]), unif.prior )
                  formula_mod <- paste(formula_mod, paste0("f(",variablesChosenUser[j],", model='",eval(parse(text=paste0("input$PrefEffectCov",j))),"', hyper=",paste0("hyper",variablesChosenUser[j]),  ", constr=",eval(parse(text=paste0("input$PrefEffectCovConstr",j))),", scale.model=TRUE)"), sep=" + ")
                } else if(eval(parse(text=paste0("input$PrefEffectCovKindPrior",j)))=="unifflat"){
                  assign(paste0("unifflat.prior",variablesChosenUser[j]), "expression:
                  log_dens = 0 - log(2) - theta/2
                  ")
                  formula_mod <- paste(formula_mod, paste0("f(",variablesChosenUser[j],", model='",eval(parse(text=paste0("input$PrefEffectCov",j))),"', hyper=list(prec=list(prior=",paste0("unifflat.prior",variablesChosenUser[j]),")), constr=",eval(parse(text=paste0("input$PrefEffectCovConstr",j))),", scale.model=TRUE)"), sep=" + ")
                } else{
                  assign(paste0("hyper",variablesChosenUser[j]), list(prec=list(prior="loggamma",param=c( as.numeric(unlist(strsplit(eval(parse(text=paste0("input$PrefEffectCovPriorBaseValues",j))), ","))) ))) )
                  formula_mod <- paste(formula_mod, paste0("f(",variablesChosenUser[j],", model='",eval(parse(text=paste0("input$PrefEffectCov",j))),"', hyper=",paste0("hyper",variablesChosenUser[j]),  ", constr=",eval(parse(text=paste0("input$PrefEffectCovConstr",j))),", scale.model=TRUE)"), sep=" + ")
                }
              } else{
                formula_mod <- paste(formula_mod, paste0("f(",variablesChosenUser[j],", model='",eval(parse(text=paste0("input$PrefEffectCov",j))),  "', constr=",eval(parse(text=paste0("input$PrefEffectCovConstr",j))),", scale.model=TRUE)"), sep=" + ")
              }
              group.cov <- inla.group(x=c(List.covariates.mesh[[variablesChosenUser[j]]], List.covariates.inf[[variablesChosenUser[j]]], List.covariates.pred[[variablesChosenUser[j]]]), n=eval(parse(text=paste0("input$PrefEffectCovNodes",j))), method="cut")

              Inf.geo.effects.list[[2]][[variablesChosenUser[j]]] <- group.cov[seq_len(n + nv)]
              Pred.geo.effects.list[[2]][[variablesChosenUser[j]]] <- group.cov[-seq_len(n + nv)]

              if(length(UserComponentsPrefSharing)>0){
                for(k in seq_along(UserComponentsPrefSharing)){
                  if(variablesChosenUser[j] %in% UserComponentsPrefSharing[[k]]){
                    ComponentPrefGroup <- "Pref"
                    formula_mod <- paste(formula_mod, paste0("f(",variablesChosenUser[j],"_",ComponentPrefGroup[k],"_copy, copy='",variablesChosenUser[j],"', fixed=FALSE)"), sep=" + ")
                    Inf.Prefs.effects.list[[paste0("Inf.", names(UserComponentsPrefSharing)[k],".effects.list")]][[2]][[paste0(variablesChosenUser[j], "_", ComponentPrefGroup[k], "_copy")]] <- group.cov[seq_len(n + nv)]
                  } else{
                    ComponentPrefGroup <- "Pref"
                    formula_mod <- paste(formula_mod, paste0("f(",variablesChosenUser[j],"_",ComponentPrefGroup[k],")"), sep=" + ")
                    Inf.Prefs.effects.list[[paste0("Inf.", names(UserComponentsPrefSharing)[k],".effects.list")]][[2]][[paste0(variablesChosenUser[j], "_", ComponentPrefGroup[k])]] <- group.cov[seq_len(n + nv)]
                  }
                }
              }

            } else if(eval(parse(text=paste0("input$PrefEffectCov",j)))=="spde1"){
              Tot_cov <- c(List.covariates.mesh[[variablesChosenUser[j]]], List.covariates.inf[[variablesChosenUser[j]]], List.covariates.pred[[variablesChosenUser[j]]])
              spde1_nodes <- seq(min(Tot_cov), max(Tot_cov), length.out=eval(parse(text=paste0("input$PrefEffectCovNodes",j))))
              mesh1d <- fm_mesh_1d(loc=spde1_nodes)

              if(eval(parse(text=paste0("input$PrefEffectCustomPrior",j)))=="custom"){
                if(eval(parse(text=paste0("input$PrefEffectCovKindPrior",j)))=="pc"){
                  hyper <- as.numeric(unlist(strsplit(eval(parse(text=paste0("input$PrefEffectCovPriorPCValues",i))), ",")))
                  assign(paste0("spde1d_",variablesChosenUser[j]) ,inla.spde2.pcmatern(mesh=mesh1d, prior.range=c(abs(diff(range(Tot_cov)))/5, 0.5), prior.sigma=c(1,0.5)))
                } else{
                  prior.range <- as.numeric(unlist(strsplit(eval(parse(text=paste0("input$PrefEffectCovPriorBaseValues",i))), ",")))[1:3]
                  prior.sigma <- as.numeric(unlist(strsplit(eval(parse(text=paste0("input$PrefEffectCovPriorBaseValues",i))), ",")))[4:6]
                  alpha <- 2; d <- 2
                  nu <-  alpha - d/2
                  kappa0 <-  log(8*nu)/2 -log(prior.range[1])
                  tau0 <-  0.5*(lgamma(nu) - lgamma(nu + d/2) - d/2*log(4*pi)) - nu*kappa0 - log(prior.sigma[1])
                  assign(paste0("spde1d_",variablesChosenUser[j]), inla.spde2.matern(mesh = mesh1d, B.tau = cbind(tau0, nu, -1), B.kappa = cbind(kappa0, -1, 0),
                                                                                     theta.prior.mean = c(prior.range[2],prior.sigma[2]), theta.prior.prec = c(prior.range[3],prior.sigma[3])))
                }
              } else {
                assign(paste0("spde1d_",variablesChosenUser[j]), inla.spde2.pcmatern(mesh=mesh1d, prior.range=c(abs(diff(range(Tot_cov)))/5, 0.5), prior.sigma=c(1,0.5)))
              }

              spde1d.index <- inla.spde.make.index(name=paste0(variablesChosenUser[j]), n.spde=eval(parse(text=paste0("spde1d_",variablesChosenUser[j])))$n.spde)
              formula_mod <- paste(formula_mod, paste0("f(",variablesChosenUser[j],", model=", paste0("spde1d_",variablesChosenUser[j]),  ")"), sep=" + ")
              Inf.geo.effects.list[[length(A_Inf.spde1)+3]] <- list()
              Pred.geo.effects.list[[length(A_Inf.spde1)+3]] <- list()
              Inf.geo.effects.list[[length(A_Inf.spde1)+3]][[variablesChosenUser[j]]] <- spde1d.index[[1]]
              Pred.geo.effects.list[[length(A_Inf.spde1)+3]][[variablesChosenUser[j]]] <- spde1d.index[[1]]

              if(length(UserComponentsPrefSharing)>0){
                for(k in seq_along(UserComponentsPrefSharing)){
                  if(variablesChosenUser[j] %in% UserComponentsPrefSharing[[k]]){
                    ComponentPrefGroup <- "Pref"
                    formula_mod <- paste(formula_mod, paste0("f(",variablesChosenUser[j],"_",ComponentPrefGroup[k],"_copy, copy='",variablesChosenUser[j],"', fixed=FALSE)"), sep=" + ")
                    Inf.Prefs.effects.list[[paste0("Inf.", names(UserComponentsPrefSharing)[k],".effects.list")]][[length(A_Inf.spde1)+3]] <- list()
                    Inf.Prefs.effects.list[[paste0("Inf.", names(UserComponentsPrefSharing)[k],".effects.list")]][[length(A_Inf.spde1)+3]][[paste0(variablesChosenUser[j], "_", ComponentPrefGroup[k], "_copy")]] <- spde1d.index[[1]]
                  } else{
                    ComponentPrefGroup <- "Pref"
                    formula_mod <- paste(formula_mod, paste0("f(",variablesChosenUser[j],"_",ComponentPrefGroup[k],", model=", paste0("spde1d_",variablesChosenUser[j]),")"), sep=" + ")
                    Inf.Prefs.effects.list[[paste0("Inf.", names(UserComponentsPrefSharing)[k],".effects.list")]][[length(A_Inf.spde1)+3]] <- list()
                    Inf.Prefs.effects.list[[paste0("Inf.", names(UserComponentsPrefSharing)[k],".effects.list")]][[length(A_Inf.spde1)+3]][[paste0(variablesChosenUser[j], "_", ComponentPrefGroup[k])]] <- spde1d.index[[1]]
                  }
                }
              }

              A_Inf.spde1[[length(A_Inf.spde1)+1]] <- inla.spde.make.A(mesh=mesh1d, loc=Tot_cov[seq_len(nv+n)])
              A_Pred.spde1[[length(A_Pred.spde1)+1]] <- inla.spde.make.A(mesh=mesh1d, loc=Tot_cov[-seq_len(nv+n)])

            } else if(eval(parse(text=paste0("input$PrefEffectCov",j)))=="linear"){
              if(eval(parse(text=paste0("input$PrefEffectCustomPrior",j)))=="custom"){
                hyper <- as.numeric(unlist(strsplit(eval(parse(text=paste0("input$PrefEffectCovPriorBaseValues",j))), ",")))
                formula_mod <- paste(formula_mod, paste0("f(",variablesChosenUser[j],", model='",eval(parse(text=paste0("input$PrefEffectCov",j))),"', mean.linear=", hyper[1], ",prec.linear=", hyper[2],")"), sep=" + ")
              } else{
                formula_mod <- paste(formula_mod, paste0("f(",variablesChosenUser[j],", model='",eval(parse(text=paste0("input$PrefEffectCov",j))),"')"), sep=" + ")
              }

              Inf.geo.effects.list[[2]][[variablesChosenUser[j]]] <- c(List.covariates.mesh[[variablesChosenUser[j]]], List.covariates.inf[[variablesChosenUser[j]]])
              Pred.geo.effects.list[[2]][[variablesChosenUser[j]]] <- c(List.covariates.pred[[variablesChosenUser[j]]])

              if(length(UserComponentsPrefSharing)>0){
                for(k in seq_along(UserComponentsPrefSharing)){
                  if(variablesChosenUser[j] %in% UserComponentsPrefSharing[[k]]){
                    showNotification(ui=paste0("From a theoretical perspective, sharing linear effects is feasible but doesn't make logical sense. Consequently, we will approach them as non-shared effects."), duration=10, closeButton=TRUE, type="warning")
                    ComponentPrefGroup <- "Pref"
                    formula_mod <- paste(formula_mod, paste0("f(",variablesChosenUser[j],"_",ComponentPrefGroup[k],", model='", eval(parse(text=paste0("input$PrefEffectCov",j))),"')"), sep=" + ")
                    Inf.Prefs.effects.list[[paste0("Inf.", names(UserComponentsPrefSharing)[k],".effects.list")]][[2]][[paste0(variablesChosenUser[j], "_", ComponentPrefGroup[k])]] <- c(List.covariates.mesh[[variablesChosenUser[j]]], List.covariates.inf[[variablesChosenUser[j]]])
                  } else{
                    ComponentPrefGroup <- "Pref"
                    formula_mod <- paste(formula_mod, paste0("f(",variablesChosenUser[j],"_",ComponentPrefGroup[k],", model='", eval(parse(text=paste0("input$PrefEffectCov",j))),"')"), sep=" + ")
                    Inf.Prefs.effects.list[[paste0("Inf.", names(UserComponentsPrefSharing)[k],".effects.list")]][[2]][[paste0(variablesChosenUser[j], "_", ComponentPrefGroup[k])]] <- c(List.covariates.mesh[[variablesChosenUser[j]]], List.covariates.inf[[variablesChosenUser[j]]])
                  }
                }
              }
            } else{
              showNotification(ui=paste("The effect of numerical covariates cannot possess an independent and identically distributed (iid) structure. If this is required, the variable values should be recorded as text, not as numerical input.."), duration = NULL)
            }
          } else{ # Section for factor variables
            if(eval(parse(text=paste0("input$PrefEffectCov",j)))=="iid"){
              PrefFactorModelDF <- data.frame(y=1, Pref=c(List.covariates.mesh[[variablesChosenUser[j]]], List.covariates.inf[[variablesChosenUser[j]]]))
              colnames(PrefFactorModelDF) <- c("y", variablesChosenUser[j])
              idx.factor <- which(variablesChosenUser[j]==names(DFsample[!as.vector(unlist(lapply(X=DFsample, FUN=is.numeric)))])[names(DFsample[!as.vector(unlist(lapply(X=DFsample, FUN=is.numeric)))])%in%c(variablesChosenUser)])

              Inf.geo.effects.list[[2]][[variablesChosenUser[j]]] <- c(List.covariates.mesh[[variablesChosenUser[j]]], List.covariates.inf[[variablesChosenUser[j]]])

              Pred.geo.effects.list[[2]][[variablesChosenUser[j]]] <- if(eval(parse(text=paste0("input$PrefKindPredictionFactorLevel",idx.factor)))=="reference"){
                rep( eval(parse(text=paste0("input$PrefEffectCovFactorPred",idx.factor))), times=length(List.covariates.pred[[variablesChosenUser[j]]]))
              } else{List.covariates.pred[[variablesChosenUser[j]]]}

              #c(List.covariates.pred[[variablesChosenUser[j]]])

              if(eval(parse(text=paste0("input$PrefEffectCustomPrior",j)))=="custom"){
                if(eval(parse(text=paste0("input$PrefEffectCovKindPrior",j)))=="pc"){
                  hyper <- list(prec=list(prior="pc.prior",param=c( as.numeric(unlist(strsplit(eval(parse(text=paste0("input$PrefEffectCovPriorPCValues",j))), ","))) )))
                  formula_mod <- paste(formula_mod, paste0("f(",variablesChosenUser[j],", model='",eval(parse(text=paste0("input$PrefEffectCov",j))),"', hyper=",hyper, ", constr=",eval(parse(text=paste0("input$PrefEffectCovConstr",j))),")"), sep=" + ")
                } else if(eval(parse(text=paste0("input$PrefEffectCovKindPrior",j)))=="base"){
                  hyper <- list(prec=list(prior="loggamma",param=c( as.numeric(unlist(strsplit(eval(parse(text=paste0("input$PrefEffectCovPriorBaseValues",j))), ","))) )))
                  formula_mod <- paste(formula_mod, paste0("f(",variablesChosenUser[j],", model='iid', hyper=",hyper,  ", constr=",eval(parse(text=paste0("input$PrefEffectCovConstr",j))),")"), sep=" + ")
                } else if(eval(parse(text=paste0("input$PrefEffectCovKindPrior",j)))=="unif"){
                  lim <- as.numeric(unlist(strsplit(eval(parse(text=paste0("input$PrefEffectCovPriorUnif",j))), ",")))
                  sigma <- seq(lim[1]+1E-5, lim[2]*3, length.out=1E5)
                  theta <- -2*log(sigma)
                  logdens <- sapply(X=sigma, FUN=logdunif, lim=lim)
                  unif.prior <- list(theta=list(prior=paste0("table: ", paste(c(theta, logdens), collapse=" "))))
                  assign(paste0("hyper", variablesChosenUser[j]), unif.prior )
                  formula_mod <- paste(formula_mod, paste0("f(",variablesChosenUser[j],", model='iid', hyper=",paste0("hyper",variablesChosenUser[j]),  ", constr=",eval(parse(text=paste0("input$PrefEffectCovConstr",j))),", scale.model=TRUE)"), sep=" + ")
                } else{ # flatunif
                  unifflat.prior= "expression:
                  log_dens = 0 - log(2) - theta/2
                  "
                  formula_mod <- paste(formula_mod, paste0("f(",variablesChosenUser[j],", model='iid', hyper=list(prec=list(prior=",unifflat.prior,")), constr=",eval(parse(text=paste0("input$PrefEffectCovConstr",j))),")"), sep=" + ")
                }
              } else{
                formula_mod <- paste(formula_mod, paste0("f(", variablesChosenUser[j], ", model='iid'",  ", constr=",eval(parse(text=paste0("input$PrefEffectCovConstr",j))),")"), sep=" + ")
              }

              if(length(UserComponentsPrefSharing)>0){
                for(k in seq_along(UserComponentsPrefSharing)){
                  if(variablesChosenUser[j] %in% UserComponentsPrefSharing[[k]]){
                    ComponentPrefGroup <- "Pref"
                    formula_mod <- paste(formula_mod, paste0("f(",variablesChosenUser[j],"_",ComponentPrefGroup[k],"_copy, copy='",variablesChosenUser[j],"')"), sep=" + ")
                    Inf.Prefs.effects.list[[paste0("Inf.", names(UserComponentsPrefSharing)[k],".effects.list")]][[2]][[paste0(variablesChosenUser[j], "_", ComponentPrefGroup[k], "_copy")]] <- c(List.covariates.mesh[[variablesChosenUser[j]]], List.covariates.inf[[variablesChosenUser[j]]])
                  } else{
                    ComponentPrefGroup <- "Pref"
                    formula_mod <- paste(formula_mod, paste0("f(",variablesChosenUser[j],"_",ComponentPrefGroup[k],")"), sep=" + ")
                    Inf.Prefs.effects.list[[paste0("Inf.", names(UserComponentsPrefSharing)[k],".effects.list")]][[2]][[paste0(variablesChosenUser[j], "_", ComponentPrefGroup[k])]] <- c(List.covariates.mesh[[variablesChosenUser[j]]], List.covariates.inf[[variablesChosenUser[j]]])
                  }
                }
              }

            } else{
              PrefFactorModelDF <- data.frame(y=1, Pref=c(List.covariates.mesh[[variablesChosenUser[j]]], List.covariates.inf[[variablesChosenUser[j]]]))
              colnames(PrefFactorModelDF) <- c("y", variablesChosenUser[j])
              FactorVariables <- data.frame( model.matrix(object=as.formula(paste0("y~-1+", variablesChosenUser[j])), PrefFactorModelDF ) )
              idx.factor <- which(variablesChosenUser[j]==names(DFsample[!as.vector(unlist(lapply(X=DFsample, FUN=is.numeric)))])[names(DFsample[!as.vector(unlist(lapply(X=DFsample, FUN=is.numeric)))])%in%c(variablesChosenUser)])
              for(l in setdiff(colnames(FactorVariables),paste0(variablesChosenUser[j], eval(parse(text=paste0("input$PrefEffectCovFactorPred",idx.factor))) ))){
                ll <- l
                Inf.geo.effects.list[[2]][[l]] <- FactorVariables[,l]
                Pred.geo.effects.list[[2]][[l]] <- rep(0, length(List.covariates.pred[[variablesChosenUser[j]]]))
                if(eval(parse(text=paste0("input$PrefEffectCustomPrior",j)))=="custom"|eval(parse(text=paste0("input$PrefEffectCov",j)))=="linear"){
                  hyper <- as.numeric(unlist(strsplit(eval(parse(text=paste0("input$PrefEffectCovPriorBaseValues",j))), ",")))
                  formula_mod <- paste(formula_mod, paste0("f(", l, ", model='linear', mean.linear=", hyper[1], ",prec.linear=", hyper[2],")"), sep=" + ")
                } else{
                  formula_mod <- paste(formula_mod, paste0("f(", l, ", model='linear')"), sep=" + ")
                }
              }

              if(length(UserComponentsPrefSharing)>0){
                for(k in seq_along(UserComponentsPrefSharing)){
                  for(l in setdiff(colnames(FactorVariables),paste0(variablesChosenUser[j], eval(parse(text=paste0("input$PrefEffectCovFactorPred",idx.factor))) ))){
                    if(variablesChosenUser[j] %in% UserComponentsPrefSharing[[k]]){
                      showNotification(ui=paste0("From a theoretical perspective, sharing linear effects is feasible but doesn't make logical sense. Consequently, we will approach them as non-shared effects."), duration=10, closeButton=TRUE, type="warning")
                      ComponentPrefGroup <- "Pref"
                      formula_mod <- paste(formula_mod, paste0("f(",l,"_",ComponentPrefGroup[k],")"), sep=" + ")
                      Inf.Prefs.effects.list[[paste0("Inf.", names(UserComponentsPrefSharing)[k],".effects.list")]][[2]][[paste0(l, "_", ComponentPrefGroup[k])]] <- c(List.covariates.mesh[[variablesChosenUser[j]]], List.covariates.inf[[variablesChosenUser[j]]])
                    } else{
                      ComponentPrefGroup <- "Pref"
                      formula_mod <- paste(formula_mod, paste0("f(",l,"_",ComponentPrefGroup[k],")"), sep=" + ")
                      Inf.Prefs.effects.list[[paste0("Inf.", names(UserComponentsPrefSharing)[k],".effects.list")]][[2]][[paste0(l, "_", ComponentPrefGroup[k])]] <- c(List.covariates.mesh[[variablesChosenUser[j]]], List.covariates.inf[[variablesChosenUser[j]]])
                    }
                  }
                }
              }


            }
          }
        }
      }

      ### Stacks of the geostatistical, Pref and prediction layers ====

      ResponseVariable <- DFsample[,3]

      A_inf_tot <- c(A.inf,1)
      if(length(A_Inf.spde1)>0){
        for(i in seq_along(A_Inf.spde1)){
          A_inf_tot[[2+i]] <- A_Inf.spde1[[i]]
        }
      }

      A_pred_tot <- c(A.geo.pred,1)
      if(length(A_Inf.spde1)>0){
        for(i in seq_along(A_Inf.spde1)){
          A_pred_tot[[2+i]] <- A_Pred.spde1[[i]]
        }
      }

      Inf.geo.stack <- inla.stack(data=list(y=cbind(c(rep(NA, nv),ResponseVariable), NA), e=rep(0,times=nv+n)),
                                   A=A_inf_tot,
                                   effects=Inf.geo.effects.list,
                                   tag="Inference_geo"
      )

      Pred.geo.stack <- inla.stack(data=list(y=matrix(NA, nrow=nrow(A.geo.pred), ncol=2)),
                                    A=A_pred_tot,
                                    effects=Pred.geo.effects.list,
                                    tag="Prediction_geo")


      if(length(UserComponentsPrefSharing)>0){
        Total.stack <- inla.stack(Inf.geo.stack)
        ComponentPrefGroup <- "Pref"
        for(k in seq_along(UserComponentsPrefSharing)){
          y.na.vec <- rep(NA, times=n)
          y.na.vec[ComponentPrefGroup[k]==DFsample[,input$UserComponentsPrefDependentGroup]] <- 1
          y.pp <- c(rep(0,times=nv), y.na.vec)
          assign(paste0("Inf.Pref",k,".stack"),
                 inla.stack(data=list(y=cbind(NA, y.pp), e=c(w, rep(0,n))),
                            A=A_inf_tot,
                            effects=Inf.Prefs.effects.list[[paste0("Inf.", names(UserComponentsPrefSharing)[k],".effects.list")]],
                            tag=paste0("Inf_Pref",k)
                 )
          )
          Total.stack <- inla.stack(Total.stack, eval(parse(text=paste0("Inf.Pref",k,".stack"))))
        }
        Total.stack <- inla.stack(Total.stack, Pred.geo.stack)
      } else{
        Total.stack <- inla.stack(Inf.geo.stack, Pred.geo.stack)
      }

      ### INLA model ====

      formula_inla <- as.formula(formula_mod)

      if(input$autocustomPrefFamily=='custom'){
        if(input$PrefFamilyPriorKind=="pc"){
          family.pcprec <- as.numeric(unlist(strsplit(input$PrefFamilyHyper,",")))
          controlFamily <- list(list(hyper = list(prec = list(prior="pc.prec", param=family.pcprec))), list())
        } else if(input$PrefFamilyPriorKind=="unif"){
          lim <- as.numeric(unlist(strsplit(input$PrefFamilyHyper,",")))
          sigma <- seq(0, lim[2]*3, length.out=1E5)
          theta <- -2*log(sigma)
          logdens <- sapply(X=sigma, FUN=logdunif, lim=lim)
          family.unif.prior <- list(theta=list(prior=paste0("table: ", paste(c(theta, logdens), collapse=" "))))
          controlFamily <- list(list(hyper = list(theta = list(prior=family.unif.prior))), list())
        } else if(input$PrefFamilyPriorKind=="unifflat"){
          unifflat.prior= "expression:
                  log_dens = 0 - log(2) - theta/2
                  "
          controlFamily <- list(list(hyper = list(prec = list(prior=unifflat.prior))), list())
        } else{
          controlFamily <- list(list(hyper = list(prec = list(prior="loggamma", param=as.numeric(unlist(strsplit(input$PrefFamilyHyper,",")))))), list())
        }
      } else{
        controlFamily <- list(list(), list())
      }

      if(input$autocustomPrefMode=='custom'){
        controlModeTheta <- list(theta=as.numeric(unlist(strsplit(input$PrefModeHyper,","))), restart=TRUE)
      } else{controlModeTheta <- inla.set.control.mode.default()}

      if(input$INLAModePref=="classic"){
        controlINLA <- list(strategy=input$strategyapproxINLAPref,
                            int.strategy=input$strategyintINLAPref)
      } else{
        controlINLA <- list()
      }

      Pref.model <- inla(formula=formula_inla, family = c(input$SelectPrefFamily,'poisson'),
                             data = inla.stack.data(Total.stack),
                             E = inla.stack.data(Total.stack)$e,
                             control.inla = controlINLA,
                             control.predictor = list(A = inla.stack.A(Total.stack), compute = TRUE, link = 1),
                             control.family = controlFamily,
                             control.mode = controlModeTheta,
                             control.compute = list(cpo = TRUE, dic = TRUE, waic = TRUE, config = TRUE),
                             inla.mode=input$INLAModePref,
                             verbose=FALSE)

      index.pred <- inla.stack.index(Total.stack, "Prediction_geo")$data
      DFpred <- data.frame(Latitude=xy.pred[,1], Longitude=xy.pred[,2])
      DFpred$Abundance.mean <- Pref.model$summary.fitted.values[index.pred, "mean"]
      DFpred$Abundance.median <- Pref.model$summary.fitted.values[index.pred, "0.5quant"]
      DFpred$Abundance.sd <- Pref.model$summary.fitted.values[index.pred, "sd"]

      colnames(DFpred)[1:2] <- colnames(DFsample)[1:2]

      DFpredictorMeanMedianStdev <- data.frame(Latitude=xy.pred[,1], Longitude=xy.pred[,2],
                                               Predictor.mean=Pref.model$summary.linear.predictor[index.pred, "mean"],
                                               Predictor.median=Pref.model$summary.linear.predictor[index.pred, "0.5quant"],
                                               Predictor.stdev=Pref.model$summary.linear.predictor[index.pred, "sd"])

      colnames(DFpredictorMeanMedianStdev)[1:2] <- colnames(DFsample)[1:2]

      gridSpatial <- expand.grid(x=seq(min(DFpred[,1]), max(DFpred[,1]),length.out=input$dimPrefmap),
                                 y=seq(min(DFpred[,2]), max(DFpred[,2]), length.out=input$dimPrefmap))
      gridSpatial <- gridSpatial[which(!is.na(over(SpatialPoints(coords=gridSpatial),SpatialPolygons(Srl=list(Polygons(srl=list(Polygon(coords=mesh$loc[mesh$segm$int$idx[,1], 1:2])), ID="int")))))),1:2]

      A.spatial <- inla.spde.make.A(mesh=mesh, loc=as.matrix(gridSpatial))

      DFspatialMeanMedianStdev <- data.frame(Latitude=as.vector(gridSpatial[,1]), Longitude=as.vector(gridSpatial[,2]),
                                              Spatial.mean=as.vector(A.spatial%*%Pref.model$summary.random$Spatial$mean),
                                              Spatial.median=as.vector(A.spatial%*%Pref.model$summary.random$Spatial$`0.5quant`),
                                              Spatial.stdev=as.vector(A.spatial%*%Pref.model$summary.random$Spatial$sd))

      colnames(DFspatialMeanMedianStdev)[1:2] <- colnames(DFsample)[1:2]

      result <- list(DFpredAbunMeanMedianStdev=list(DFpredAbunMeanMedianStdev=DFpred),
                     DFpredictorMeanMedianStdev=list(DFpredictorMeanMedianStdev=DFpredictorMeanMedianStdev),
                     DFspatialMeanMedianStdev=list(DFspatialMeanMedianStdev=DFspatialMeanMedianStdev),
                     PrefModel=list(PrefModel=Pref.model),
                     DFPostFixed=list(DFPostFixed=Pref.model$marginals.fixed),
                     DFPostHyperpar=list(DFPostHyperpar=Pref.model$marginals.hyperpar),
                     Summary.fixed=list(Summary.fixed=Pref.model$summary.fixed),
                     Summary.hyperpar=list(Summary.hyperpar=Pref.model$summary.hyperpar),
                     SummaryInternalHyper=list(SummaryInternalHyper=Pref.model$internal.summary.hyperpar),
                     SummaryCPO=list(SummaryCPO=na.omit(Pref.model$cpo$cpo)),
                     DICmodel=list(DICmodel=data.frame(DIC=Pref.model$dic$family.dic, row.names=if(length(UserComponentsPrefSharing)>0){c("Geostatistical", "Point process")}else{c("Geostatistical")})),
                     WAICmodel=list(WAICmodel=data.frame(WAIC=cbind(unlist(lapply(na.omit(unique(Pref.model$dic$family)), function(i){sum(Pref.model$waic$local.waic[which(Pref.model$dic$family==i)])}))), row.names=if(length(UserComponentsPrefSharing)>0){c("Geostatistical", "Point process")}else{c("Geostatistical")}))
                     )

      t2 <- Sys.time()
      difftime(t2,t1, units="secs")
      showNotification(ui=paste("The model has been fitted:", as.numeric(round(Pref.model$cpu.used[4])),
                                "(Preferential model) and", as.numeric(round(difftime(t2,t1, units="secs"))),
                                "(overall process) secs." ), duration = NULL)
      # showNotification(ui=paste("The model's DIC is", Preferential.model$dic$dic), duration = NULL)
      return(result)
    })


    dataggplotPrefAbundanceMeanMedianStdevFit <- function(){
      DF <- PrefModelFit()$DFpredAbunMeanMedianStdev$DFpredAbunMeanMedianStdev
      return(DF)
    }

    DFPrefAbundance <- reactive({
      DF <- PrefModelFit()$DFpredAbunMeanMedianStdev$DFpredAbunMeanMedianStdev
      condition <- !(length(unique(DF[,1]))*length(unique(DF[,2]))==nrow(DF))&val()
      if(condition){
        grid <- expand.grid(seq(min(DF[,1]),max(DF[,1]), length.out=200),seq(min(DF[,2]),max(DF[,2]), length.out=200))
        DFInter <- data.frame(Latitude=grid[,1],Longitude=grid[,2])
        DFInter$Abundance.mean <- InterpolateIrrGrid(z=DF$Abundance.mean,loc=DF[,1:2], gridInter=grid)$DataInter$z
        DFInter$Abundance.median <- InterpolateIrrGrid(z=DF$Abundance.median,loc=DF[,1:2], gridInter=grid)$DataInter$z
        DFInter$Abundance.sd <- InterpolateIrrGrid(z=DF$Abundance.sd,loc=DF[,1:2], gridInter=grid)$DataInter$z
        colnames(DFInter)[1:2] <- colnames(DF)[1:2]
        DF <- DFInter
      }

      return(DF)
    })

    ggplotPrefAbundanceMeanMedianStdevFit <- function(){
      DF <- DFPrefAbundance()
      g1 <- ggplot(DF) + geom_tile(aes(x=DF[,1], y=DF[,2], fill=Abundance.mean)) +
        scale_fill_viridis_c(option = "turbo") + theme_bw() + coord_fixed() + labs(fill="Values") +
        xlab(colnames(DF)[1]) + ylab(colnames(DF)[1]) + ggtitle("Response predicted (mean)") +
        theme(plot.title=element_text(color = "black", size = 16, face = "bold", hjust = 0.5))

      g2 <- ggplot(DF) + geom_tile(aes(x=DF[,1], y=DF[,2], fill=Abundance.median)) +
        scale_fill_viridis_c(option = "turbo") + theme_bw() + coord_fixed() + labs(fill="Values") +
        xlab(colnames(DF)[1]) + ylab(colnames(DF)[1]) + ggtitle("Response predicted (median)") +
        theme(plot.title=element_text(color = "black", size = 16, face = "bold", hjust = 0.5))

      g3 <- ggplot(DF) + geom_tile(aes(x=DF[,1], y=DF[,2], fill=Abundance.sd))+
        scale_fill_viridis_c(option = "turbo") + theme_bw() + coord_fixed() + labs(fill="Values") +
        xlab(colnames(DF)[1]) + ylab(colnames(DF)[1]) + ggtitle("Response predicted (stdev.)") +
        theme(plot.title=element_text(color = "black", size = 16, face = "bold", hjust = 0.5))

      gt <- g1+g2+g3

      return(gt)
    }

    downloadFile(
      id = "ggplotPrefAbundanceMeanMedianStdevFit",
      logger = ss_userAction.Log,
      filenameroot = "ggplotPrefAbundanceMeanMedianStdevFit",
      aspectratio  = 1,
      downloadfxns = list(png  = ggplotPrefAbundanceMeanMedianStdevFit,
                          csv = dataggplotPrefAbundanceMeanMedianStdevFit,
                          txt = dataggplotPrefAbundanceMeanMedianStdevFit)
    )

    downloadablePlot(id = "ggplotPrefAbundanceMeanMedianStdevFit",
                     logger = ss_userAction.Log,
                     filenameroot = "ggplotPrefAbundanceMeanMedianStdevFit",
                     aspectratio  = 1,
                     downloadfxns = list(png = ggplotPrefAbundanceMeanMedianStdevFit,
                                         csv = dataggplotPrefAbundanceMeanMedianStdevFit,
                                         txt = dataggplotPrefAbundanceMeanMedianStdevFit),
                     visibleplot  = ggplotPrefAbundanceMeanMedianStdevFit)


    dataggplotPrefSpatialMeanMedianStdev <- function(){
      DF <- PrefModelFit()$DFspatialMeanMedianStdev$DFspatialMeanMedianStdev
      return(DF)
    }

    ggplotPrefSpatialMeanMedianStdev <- function(){
      DF <- PrefModelFit()$DFspatialMeanMedianStdev$DFspatialMeanMedianStdev
      g1 <- ggplot(DF) + geom_tile(aes(x=DF[,1], y=DF[,2], fill=Spatial.mean)) +
        scale_fill_viridis_c(option = "turbo") + theme_bw() + coord_fixed() + labs(fill="Values") +
        xlab(colnames(DF)[1]) + ylab(colnames(DF)[1]) + ggtitle("Posterior spatial (mean)") +
        theme(plot.title=element_text(color = "black", size = 16, face = "bold", hjust = 0.5))

      g2 <- ggplot(DF) + geom_tile(aes(x=DF[,1], y=DF[,2], fill=Spatial.median)) +
        scale_fill_viridis_c(option = "turbo") + theme_bw() + coord_fixed() + labs(fill="Values") +
        xlab(colnames(DF)[1]) + ylab(colnames(DF)[1]) + ggtitle("Posterior spatial (median)") +
        theme(plot.title=element_text(color = "black", size = 16, face = "bold", hjust = 0.5))

      g3 <- ggplot(DF) + geom_tile(aes(x=DF[,1], y=DF[,2], fill=Spatial.stdev))+
        scale_fill_viridis_c(option = "turbo") + theme_bw() + coord_fixed() + labs(fill="Values") +
        xlab(colnames(DF)[1]) + ylab(colnames(DF)[1]) + ggtitle("Posterior spatial (stdev.)") +
        theme(plot.title=element_text(color = "black", size = 16, face = "bold", hjust = 0.5))

      gt <- g1+g2+g3

      return(gt)
    }

    downloadFile(
      id = "ggplotPrefSpatialMeanMedianStdev",
      logger = ss_userAction.Log,
      filenameroot = "ggplotPrefSpatialMeanMedianStdev",
      aspectratio  = 1,
      downloadfxns = list(png  = ggplotPrefSpatialMeanMedianStdev,
                          csv = dataggplotPrefSpatialMeanMedianStdev,
                          txt = dataggplotPrefSpatialMeanMedianStdev)
    )

    downloadablePlot(id = "ggplotPrefSpatialMeanMedianStdev",
                     logger = ss_userAction.Log,
                     filenameroot = "ggplotPrefSpatialMeanMedianStdev",
                     aspectratio  = 1,
                     downloadfxns = list(png=ggplotPrefSpatialMeanMedianStdev,
                                         csv = dataggplotPrefSpatialMeanMedianStdev,
                                         txt = dataggplotPrefSpatialMeanMedianStdev),
                     visibleplot  = ggplotPrefSpatialMeanMedianStdev)

    dataggplotPrefPredictorMeanMedianStdevFit <- function(){
      DF <- PrefModelFit()$DFpredictorMeanMedianStdev$DFpredictorMeanMedianStdev
      return(DF)
    }

    DFPrefPredictor <- reactive({
      DF <- PrefModelFit()$DFpredictorMeanMedianStdev$DFpredictorMeanMedianStdev
      condition <- !(length(unique(DF[,1]))*length(unique(DF[,2]))==nrow(DF))&val()
      if(condition){
        grid <- expand.grid(seq(min(DF[,1]),max(DF[,1]), length.out=200),seq(min(DF[,2]),max(DF[,2]), length.out=200))
        DFInter <- data.frame(Latitude=grid[,1],Longitude=grid[,2])
        DFInter$Predictor.mean <- InterpolateIrrGrid(z=DF$Predictor.mean,loc=DF[,1:2], gridInter=grid)$DataInter$z
        DFInter$Predictor.median <- InterpolateIrrGrid(z=DF$Predictor.median,loc=DF[,1:2], gridInter=grid)$DataInter$z
        DFInter$Predictor.stdev <- InterpolateIrrGrid(z=DF$Predictor.stdev,loc=DF[,1:2], gridInter=grid)$DataInter$z
        colnames(DFInter)[1:2] <- colnames(DF)[1:2]
        DF <- DFInter
      }

      return(DF)
    })

    ggplotPrefPredictorMeanMedianStdevFit <- function(){
      DF <- DFPrefPredictor()
      g1 <- ggplot(DF) + geom_tile(aes(x=DF[,1], y=DF[,2], fill=Predictor.mean)) +
        scale_fill_viridis_c(option = "turbo") + theme_bw() + coord_fixed() + labs(fill="Values") +
        xlab(colnames(DF)[1]) + ylab(colnames(DF)[1]) + ggtitle("Linear predictor (mean)") +
        theme(plot.title=element_text(color = "black", size = 16, face = "bold", hjust = 0.5))

      g2 <- ggplot(DF) + geom_tile(aes(x=DF[,1], y=DF[,2], fill=Predictor.median)) +
        scale_fill_viridis_c(option = "turbo") + theme_bw() + coord_fixed() + labs(fill="Values") +
        xlab(colnames(DF)[1]) + ylab(colnames(DF)[1]) + ggtitle("Linear predictor (median)") +
        theme(plot.title=element_text(color = "black", size = 16, face = "bold", hjust = 0.5))

      g3 <- ggplot(DF) + geom_tile(aes(x=DF[,1], y=DF[,2], fill=Predictor.stdev))+
        scale_fill_viridis_c(option = "turbo") + theme_bw() + coord_fixed() + labs(fill="Values") +
        xlab(colnames(DF)[1]) + ylab(colnames(DF)[1]) + ggtitle("Linear predictor (stdev.)") +
        theme(plot.title=element_text(color = "black", size = 16, face = "bold", hjust = 0.5))

      gt <- g1+g2+g3

      return(gt)
    }

    downloadFile(
      id = "ggplotPrefPredictorMeanMedianStdevFit",
      logger = ss_userAction.Log,
      filenameroot = "ggplotPrefPredictorMeanMedianStdevFit",
      aspectratio  = 1,
      downloadfxns = list(png  = ggplotPrefPredictorMeanMedianStdevFit,
                          csv = dataggplotPrefPredictorMeanMedianStdevFit,
                          txt = dataggplotPrefPredictorMeanMedianStdevFit)
    )

    downloadablePlot(id = "ggplotPrefPredictorMeanMedianStdevFit",
                     logger = ss_userAction.Log,
                     filenameroot = "ggplotPrefPredictorMeanMedianStdevFit",
                     aspectratio  = 1,
                     downloadfxns = list(png = ggplotPrefPredictorMeanMedianStdevFit,
                                         csv = dataggplotPrefPredictorMeanMedianStdevFit,
                                         txt = dataggplotPrefPredictorMeanMedianStdevFit),
                     visibleplot  = ggplotPrefPredictorMeanMedianStdevFit)

    dataggplotPrefFixParamFit <- function(){
      DF <- as.data.frame(PrefModelFit()$DFPostFixed$DFPostFixed)
      return(DF)
    }

    ggplotPrefFixParamFit <- function(){
      DF <- PrefModelFit()$DFPostFixed$DFPostFixed
      gl <- c()
      for(i in 1:length(DF)){
        assign(paste0("g",i),
               ggplot(data=data.frame(x=DF[[i]][,1], y=DF[[i]][,2]), aes(x=x,y=y)) + geom_line() +
                 theme_bw() + xlab(names(DF)[i]) + ylab(HTML(paste("Density f(",names(DF)[i],")"))) +
                 ggtitle(names(DF)[i]) + theme(plot.title=element_text(color = "black", size = 16, face = "bold", hjust = 0.5))
        )
        gl[i] <- paste0("g",i)
      }
      gt <- eval(parse(text=paste(gl, collapse="+")))
      return(gt)
    }

    downloadFile(
      id = "ggplotPrefFixParamFit",
      logger = ss_userAction.Log,
      filenameroot = "ggplotPrefFixParamFit",
      aspectratio  = 1,
      downloadfxns = list(png  = ggplotPrefFixParamFit,
                          csv = dataggplotPrefFixParamFit,
                          txt = dataggplotPrefFixParamFit)
    )

    downloadablePlot(id = "ggplotPrefFixParamFit",
                     logger = ss_userAction.Log,
                     filenameroot = "ggplotPrefFixParamFit",
                     aspectratio  = 1,
                     downloadfxns = list(png = ggplotPrefFixParamFit,
                                         csv = dataggplotPrefFixParamFit,
                                         txt = dataggplotPrefFixParamFit),
                     visibleplot  = ggplotPrefFixParamFit)


    dataggplotPrefHyperParamFit <- function(){
      DF <- as.data.frame(PrefModelFit()$DFPostHyperpar$DFPostHyperpar)
      return(DF)
    }

    ggplotPrefHyperParamFit <- function(){
      DF <- PrefModelFit()$DFPostHyperpar$DFPostHyperpar
      gl <- c()
      for(i in 1:length(DF)){
        nm <- strsplit(names(DF)[i], " ")[[1]]
        if(nm[1]=="Precision"){
          namesDFold <- names(DF)[i]
          names(DF)[i] <- paste("Stdev.", paste(nm[-1], collapse=" "), collapse=" ")
          assign(paste0("g",i), try(
            ggplot(data=data.frame(x=inla.tmarginal(function(x) 1/x**0.5, DF[[i]])[,1],
                                   y=inla.tmarginal(function(x) 1/x**0.5, DF[[i]])[,2]), aes(x=x,y=y)) + geom_line() +
              theme_bw()+ xlab(names(DF)[i]) +
              ylab(HTML(paste("Density f(",names(DF)[i],")"))) + ggtitle(names(DF)[i]) +
              theme(plot.title=element_text(color = "black", size = 16, face = "bold", hjust = 0.5)), silent=TRUE)
          )
          if(sum(class(eval(parse(text=paste0("g",i))))=="try-error")==1){
            assign(paste0("g",i),
                   ggplot(data=data.frame(x=inla.smarginal(marginal=DF[[i]])[[1]],
                                          y=inla.smarginal(marginal=DF[[i]])[[2]]), aes(x=x,y=y)) + geom_line() +
                     theme_bw()+ xlab(names(DF)[i]) +
                     ylab(HTML(paste("Density f(",namesDFold,")"))) + ggtitle(namesDFold) +
                     theme(plot.title=element_text(color = "black", size = 16, face = "bold", hjust = 0.5))
            )
          }
        } else{
          assign(paste0("g",i),
                 ggplot(data=data.frame(x=DF[[i]][,1], y=DF[[i]][,2]), aes(x=x,y=y)) + geom_line() +
                   theme_bw()+ xlab(names(DF)[i]) + xlim(quantile(DF[[i]][,1], probs = c(0.025,0.975))) +
                   ylab(HTML(paste("Density f(",names(DF)[i],")"))) + ggtitle(names(DF)[i]) +
                   theme(plot.title=element_text(color = "black", size = 16, face = "bold", hjust = 0.5))
          )
        }
        gl[i] <- paste0("g",i)
      }
      gt <- eval(parse(text=paste(gl, collapse="+")))
      return(gt)
    }

    downloadFile(
      id = "ggplotPrefHyperParamFit",
      logger = ss_userAction.Log,
      filenameroot = "ggplotPrefHyperParamFit",
      aspectratio  = 1,
      downloadfxns = list(png  = ggplotPrefHyperParamFit,
                          csv = dataggplotPrefHyperParamFit,
                          txt = dataggplotPrefHyperParamFit)
    )

    downloadablePlot(id = "ggplotPrefHyperParamFit",
                     logger = ss_userAction.Log,
                     filenameroot = "ggplotPrefHyperParamFit",
                     aspectratio  = 1,
                     downloadfxns = list(png=ggplotPrefHyperParamFit,
                                         csv = dataggplotPrefHyperParamFit,
                                         txt = dataggplotPrefHyperParamFit),
                     visibleplot  = ggplotPrefHyperParamFit)


    tablePrefModelFixedPar <- function(){
      DF <- PrefModelFit()$Summary.fixed$Summary.fixed %>%
        mutate(across(where(is.numeric), round, digits = 2))
    }

    downloadableTable("tablePrefModelFixedPar",
                      logger=ss_userAction.Log,
                      filenameroot="tablePrefModelFixedPar",
                      downloaddatafxns=list(csv=tablePrefModelFixedPar,
                                            tsv=tablePrefModelFixedPar),
                      tabledata=tablePrefModelFixedPar, rownames = TRUE,
                      caption="Summary fixed parameters")

    tablePrefModelHyperPar <- function(){
      DF <- PrefModelFit()$Summary.hyperpar$Summary.hyperpar %>%
        mutate(across(where(is.numeric), round, digits = 2))
    }

    downloadableTable("tablePrefModelHyperPar",
                      logger=ss_userAction.Log,
                      filenameroot="tablePrefModelHyperPar",
                      downloaddatafxns=list(csv=tablePrefModelHyperPar,
                                            tsv=tablePrefModelHyperPar),
                      tabledata=tablePrefModelHyperPar, rownames = TRUE,
                      caption="Summary hyperparameters")

    tablePrefModelInternalHyperPar <- function(){
      DF <- PrefModelFit()$SummaryInternalHyper$SummaryInternalHyper %>%
        mutate(across(where(is.numeric), round, digits = 2))
    }

    downloadableTable("tablePrefModelInternalHyperPar",
                      logger=ss_userAction.Log,
                      filenameroot="tablePrefModelInternalHyperPar",
                      downloaddatafxns=list(csv=tablePrefModelInternalHyperPar,
                                            tsv=tablePrefModelInternalHyperPar),
                      tabledata=tablePrefModelInternalHyperPar, rownames = TRUE,
                      caption="Summary internal hyperparameters")

    dataPrefDICtable <- function(){
      DF <- PrefModelFit()$DICmodel$DICmodel %>% # data.frame(DIC=PreferentialModelFit()$DICmodel$DICmodel) %>%
        mutate(across(where(is.numeric), round, digits = 2))
      return(DF)
    }

    downloadableTable("dataPrefDICtable",
                      logger=ss_userAction.Log,
                      filenameroot="dataPrefDICtable",
                      downloaddatafxns=list(csv=dataPrefDICtable,
                                            tsv=dataPrefDICtable),
                      tabledata=dataPrefDICtable, rownames = TRUE,
                      caption="Model DIC")
    
    dataPrefWAICtable <- function(){
      DF <- PrefModelFit()$WAICmodel$WAICmodel %>% # data.frame(DIC=PreferentialModelFit()$DICmodel$DICmodel) %>%
        mutate(across(where(is.numeric), round, digits = 2))
      return(DF)
    }
    
    downloadableTable("dataPrefWAICtable",
                      logger=ss_userAction.Log,
                      filenameroot="dataPrefWAICtable",
                      downloaddatafxns=list(csv=dataPrefWAICtable,
                                            tsv=dataPrefWAICtable),
                      tabledata=dataPrefWAICtable, rownames = TRUE,
                      caption="Model WAIC")

    dataPrefCPOtable <- function(){
      CPO <- PrefModelFit()$SummaryCPO$SummaryCPO
      DF <- data.frame(n=length(CPO), mean=mean(CPO), median=median(CPO), stdev=sd(CPO),
                       quantile2.5=quantile(CPO,probs=c(0.025,0.975))[1],
                       quantile97.5=quantile(CPO,probs=c(0.025,0.975))[2]) %>%
        mutate(across(where(is.numeric), round, digits = 2))
      return(DF)
    }

    downloadableTable("dataPrefCPOtable",
                      logger=ss_userAction.Log,
                      filenameroot="dataPrefCPOtable",
                      downloaddatafxns=list(csv=dataPrefCPOtable,
                                            tsv=dataPrefCPOtable),
                      tabledata=dataPrefCPOtable, rownames = FALSE,
                      caption="Summary CPO")
    
    # Server_Mixture_Section ----
    
    MixtureCheckBoxNames <- function(){
      if(input$MixtureDataSimulatedLoaded=="load"){
        DF <- as.data.frame(datareadSample())
        if(input$SelectMixtureFamily=="binomial"){DFnames <- names(DF)[c(5:ncol(DF))]}
        else{DFnames <- names(DF)[c(4:ncol(DF))]}
      } else if(input$MixtureDataSimulatedLoaded=="sim"){
        DF <- Mixture.sampling()$MixtureSample
        DFnames <- names(DF)[c(4,5)]
      }
      return(DFnames)
    }
    
    observe({
      output$checkBoxMixtureDataFrame <- renderUI({
        tagList(
          checkboxGroupInput(inputId="UserComponentsMixture",
                             label="User Defined Components",
                             choices=MixtureCheckBoxNames(),
                             selected=c()
          )
        )
      })
    })
    
    observe({
      output$checkBoxSelectMixtureDependents <- renderUI({
        tagList(
          checkboxGroupInput(inputId="UserComponentsMixtureDependentGroup",
                             label="Selection of Mixture Group ID",
                             choices=MixtureCheckBoxNames(),
                             selected=c()
          )
        )
      })
    })
    
    observe({
      if(input$MixtureDataSimulatedLoaded=="sim"){DF <- Mixture.sampling()$MixtureSample[,input$UserComponentsMixtureDependentGroup]}
      else{DF <- as.data.frame(datareadSample())[,input$UserComponentsMixtureDependentGroup]}
      output$checkBoxMixtureSharing <- renderUI({
        if(length(which(unique(DF)!="Ind"))>0){
          box(id="MixtureCovariateSharingTerms", width=12, title="Sharing Components Mixtures",
              lapply(seq_along(which(unique(DF)!="Ind")), function(i){
                checkboxGroupInput(inputId=paste0("UserComponentsMixtureSharing_",i),
                                   label=paste("User Defined Sharing Components: Mixture",i, "(",unique(DF)[unique(DF)!="Ind"][i],")"),
                                   choices=c(input$DefaultComponentsMixture,input$UserComponentsMixture),
                                   selected=c(input$DefaultComponentsMixture,input$UserComponentsMixture)
                )
              }))
        } else{}
        })
    })

    observe({
      if(input$MixtureDataSimulatedLoaded=="load"){
        output$SelectMixtureEffectCov <- renderUI({
          
            if(length(input$UserComponentsMixture)>0){
              box(id="MixtureCovariateEffects", width=12, title="Covariate Effects",
                  lapply(seq_along(input$UserComponentsMixture), function(i){
                    list(
                      selectInput(inputId=paste0("MixtureEffectCov",i), label=tags$span(style="color: blue; font-weight: bold;", paste(input$UserComponentsMixture[i],"Effect")) ,
                                  choices = list("Linear (or reference level)" = "linear", "Random Walk 1" = "rw1", "Random Walk 2" = "rw2", "SPDE 1" = "spde1", "IID"="iid"),
                                  selected = "linear"),
                      conditionalPanel(condition=paste0("input.MixtureEffectCov",i,"=='rw1'||","input.MixtureEffectCov",i,"=='rw2'||","input.MixtureEffectCov",i,"=='spde1'"),
                                       numericInput(inputId=paste0("MixtureEffectCovNodes",i),
                                                    label="Number of nodes",
                                                    value=10, min=1, step=1
                                                    )
                      ),
                      radioGroupButtons(
                        inputId = paste0("MixtureEffectCustomPrior",i),
                        label = "Custom Prior",
                        choices = list("Auto" = "auto", "Custom" = "custom"),
                        status = "primary"
                      ),
                      conditionalPanel(condition=paste0("input.MixtureEffectCustomPrior",i,"=='custom'"),
                                       list(
                                         selectInput(inputId = paste0("MixtureEffectCovKindPrior",i),
                                                 label = "Prior distribution",
                                                 choices = list("Base" = "base", "PC prior" = "pc", "Uniform" = "unif", "Flat Uniform" = "flatunif")),
                                         conditionalPanel(condition=paste0("input.MixtureEffectCovKindPrior",i,"=='base'"),
                                                          textInput(inputId = paste0("MixtureEffectCovPriorBaseValues",i),
                                                                    label = "Base Prior Values",
                                                                    value = "0, 1e-3")),
                                         conditionalPanel(condition=paste0("input.MixtureEffectCovKindPrior",i,"=='pc'"),
                                         textInput(inputId = paste0("MixtureEffectCovPriorPCValues",i),
                                                 label = "PC-Prior Values",
                                                 value = "0, 1e-3")),
                                         conditionalPanel(condition=paste0("input.MixtureEffectCovKindPrior",i,"=='unif'"),
                                                          textInput(inputId = paste0("MixtureEffectCovPriorUnif",i),
                                                                    label = "Uniform Lower and upper values",
                                                                    value = "0, 10",
                                                                    placeholder = "Lower and upper values: 'U(a,b)'"))
                                         )
                      ),
                      conditionalPanel(condition=paste0("input.MixtureEffectCov",i,"=='rw1'||",
                                                        "input.MixtureEffectCov",i,"=='rw2'||",
                                                        "input.MixtureEffectCov",i,"=='iid'"),
                                       radioGroupButtons(
                                         inputId = paste0("MixtureEffectCovConstr",i),
                                         label = "Sum to zero",
                                         choices = list("TRUE" = "TRUE", "FALSE" = "FALSE"),
                                         status = "primary"
                                       ))
                    )
                  }))
            } else{}
        }) } else if(input$MixtureDataSimulatedLoaded=="sim"){
          output$SelectMixtureEffectCov <- renderUI({
            if(length(input$UserComponentsMixture)>0){
              box(id="MixtureCovariateEffects", width=12, title="Covariate Effects",
                  lapply(seq_along(input$UserComponentsMixture), function(i){
                    list(
                      selectInput(inputId=paste0("MixtureEffectCov",i), label=tags$span(style="color: blue; font-weight: bold;", paste(input$UserComponentsMixture[i],"Effect")) ,
                                  choices = list("Linear (or reference level)" = "linear", "Random Walk 1" = "rw1", "Random Walk 2" = "rw2", "SPDE 1" = "spde1", "IID"="iid"),
                                  selected = "linear"),
                      conditionalPanel(condition=paste0("input.MixtureEffectCov",i,"=='rw1'||","input.MixtureEffectCov",i,"=='rw2'||","input.MixtureEffectCov",i,"=='spde1'"),
                                       numericInput(inputId=paste0("MixtureEffectCovNodes",i),
                                                    label="Number of nodes",
                                                    value=10, min=1, step=1
                                                    )
                      ),
                      radioGroupButtons(
                        inputId = paste0("MixtureEffectCustomPrior",i),
                        label = "Custom Prior",
                        choices = list("Auto" = "auto", "Custom" = "custom"),
                        status = "primary"
                      ),
                      conditionalPanel(condition=paste0("input.MixtureEffectCustomPrior",i,"=='custom'"),
                                       list(
                                         selectInput(inputId = paste0("MixtureEffectCovKindPrior",i),
                                                 label = "Prior distribution",
                                                 choices = list("Base" = "base", "PC prior" = "pc", "Uniform" = "unif", "Flat Uniform" = "unifflat")),
                                         conditionalPanel(condition=paste0("input.MixtureEffectCovKindPrior",i,"=='base'"),
                                                          textInput(inputId = paste0("MixtureEffectCovPriorBaseValues",i),
                                                                    label = "Base Prior Values",
                                                                    value = "0, 1e-3")),
                                         conditionalPanel(condition=paste0("input.MixtureEffectCovKindPrior",i,"=='pc'"),
                                         textInput(inputId = paste0("MixtureEffectCovPriorPCValues",i),
                                                 label = "PC-Prior Values",
                                                 value = "0, 1e-3")),
                                         conditionalPanel(condition=paste0("input.MixtureEffectCovKindPrior",i,"=='unif'"),
                                                          textInput(inputId = paste0("MixtureEffectCovPriorUnif",i),
                                                                    label = "Uniform Lower and upper values",
                                                                    value = "0, 10",
                                                                    placeholder = "Lower and upper values: 'U(a,b)'"))
                                         )
                      ),
                      conditionalPanel(condition=paste0("input.MixtureEffectCov",i,"=='rw1'||",
                                                        "input.MixtureEffectCov",i,"=='rw2'||",
                                                        "input.MixtureEffectCov",i,"=='iid'"),
                                       radioGroupButtons(
                                         inputId = paste0("MixtureEffectCovConstr",i),
                                         label = "Sum to zero",
                                         choices = list("TRUE" = "TRUE", "FALSE" = "FALSE"),
                                         status = "primary"
                                       ))
                    )
                  }))
            } else{}
          }) }
    })
    
    observe({
      if(input$MixtureDataSimulatedLoaded=="load"){
        DF <- as.data.frame(datareadSample())
        if(length(input$UserComponentsMixture)>0){
          MixtureUserComponent <- input$UserComponentsMixture
          DF2 <- select(DF, MixtureUserComponent[!as.vector(unlist(lapply(X=select(DF,MixtureUserComponent), FUN=is.numeric)))])
          output$SelectLoadMixtureEffectCovFactorPred <- renderUI({
            if(ncol(DF2)>0){
              box(id="MixturePredFactorLevel", width=12, title="Prediction Factor Level",
                  lapply(seq_along(names(DF2)), function(i){
                    choices <- unique(DF2[[names(DF2)[i]]])
                    list(
                      radioGroupButtons(inputId = paste0("MixtureKindPredictionFactorLevel",i), label = tags$span(style="color: blue; font-weight: bold;", paste(names(DF2)[i], "(prediction protocol)")), 
                                        choices = c("Reference level" = "reference", "Nearest level" = "nearest"), status = "success", justified = TRUE),
                      conditionalPanel(
                        condition=paste0("input.MixtureKindPredictionFactorLevel",i,"=='reference'"),
                        selectInput(inputId=paste0("MixtureEffectCovFactorPred",i), label=paste(names(DF2)[i],"Reference Factor"),
                                    choices = choices, selected = choices[1])
                      )
                    )
                  }))
            }
          })
        } else{
          output$SelectLoadMixtureEffectCovFactorPred <- renderUI({})
        }
      }
    })
    
    JointPriorPreviewMixture <- eventReactive(input$Mixturepreviewpriordistributions,{
      d <- 2; nu <- 1
      rhoinputPC <- as.numeric(unlist(strsplit(input$Mixturerangepcpriorprev, ",")))
      sigmainputPC <- as.numeric(unlist(strsplit(input$Mixturesigmapcpriorprev, ",")))
      rhoPC <- seq(rhoinputPC[1], rhoinputPC[2], length.out=rhoinputPC[3])
      sigmaPC <- seq(sigmainputPC[1], sigmainputPC[2], length.out=sigmainputPC[3])
      rhosigmaPC <- expand.grid(rho=rhoPC, sigma=sigmaPC)
      rho0PC <- rhoinputPC[4]; alpha1 <- rhoinputPC[5]
      sigma0PC <- sigmainputPC[4]; alpha2 <- sigmainputPC[5]
      
      lambdarhoPC <- -log(alpha1)*rho0PC**(d/2) #-(rho0/sqrt(8*nu))**(d/2)*log(alpha1)
      lambdasigmaPC <- -log(alpha2)/sigma0PC #-(sqrt(8*nu)/rhosigma$rho)**(-nu)*sqrt(gamma(nu)/(gamma(nu+d/2)*(4*pi)**(d/2)))*log(alpha2)/sigma0
      pirhosigmaPCM <- d/2*lambdarhoPC*rhosigmaPC$rho**(-1-d/2)*exp(-lambdarhoPC*rhosigmaPC$rho**(-d/2))*lambdasigmaPC*exp(-lambdasigmaPC*rhosigmaPC$sigma)
      
      probIntPC <- diff(range(rhoPC))*diff(range(sigmaPC))/length(pirhosigmaPCM)*sum(pirhosigmaPCM)
      
      rhoinputBase <- as.numeric(unlist(strsplit(input$Mixturerangebasepriorprev, ",")))
      sigmainputBase <- as.numeric(unlist(strsplit(input$Mixturesigmabasepriorprev, ",")))
      rho <- seq(rhoinputBase[1], rhoinputBase[2], length.out=rhoinputBase[3])
      sigma <- seq(sigmainputBase[1], sigmainputBase[2], length.out=sigmainputBase[3])
      rhosigma <- expand.grid(rho=rho, sigma=sigma)
      rho0 <- rhoinputBase[4]; sigma0 <- sigmainputBase[4]
      meanrho <- rhoinputBase[5]; meansigma <- sigmainputBase[5]
      sdrho <- rhoinputBase[6]; sdsigma <- sigmainputBase[6]
      pirho <- (sqrt(2*pi*sdrho**2)*rhosigma$rho)**(-1)*exp(-(log(rhosigma$rho/rho0)**2-2*log(rhosigma$rho/rho0)*meanrho+meanrho**2)/(2*sdrho**2))
      pisigma <- (sqrt(2*pi*sdsigma**2)*rhosigma$sigma)**(-1)*exp(-(log(rhosigma$sigma/sigma0)**2-2*log(rhosigma$sigma/sigma0)*meansigma+meansigma**2)/(2*sdsigma**2))
      pirhosigmaM <- pirho*pisigma
      
      probInt <- sum(pirhosigmaM)*diff(range(rhosigma$rho))*diff(range(rhosigma$sigma))/length(pirhosigmaM)
      
      ListPreview <- list(rhosigmaPC=rhosigmaPC, pirhosigmaPCM=pirhosigmaPCM, probIntPC=probIntPC,
                          rhosigma=rhosigma, pirhosigmaM=pirhosigmaM, probInt=probInt)
      
      return(ListPreview)
    })
    
    MixturePreviewJointPriorPlot <- function(){
      ggplotPCprior <- ggplot() +
        geom_tile(data=data.frame(rho=JointPriorPreviewMixture()$rhosigmaPC$rho, sigma=JointPriorPreviewMixture()$rhosigmaPC$sigma, pi=JointPriorPreviewMixture()$pirhosigmaPCM),
                  mapping=aes(x=rho, y=sigma, fill=pi)) + 
        labs(title=paste("Joint PC-Prior. Cumulative Prob.=", round(JointPriorPreviewMixture()$probIntPC, digits=2)),
             x=expression(rho), y=expression(sigma)) +
        scale_fill_viridis_c(option="turbo") + theme_bw() + theme(plot.title=element_text(face='bold'))
      
      ggplotENpriorAnalytic <- ggplot() +
        geom_tile(data=data.frame(rho=JointPriorPreviewMixture()$rhosigma$rho, sigma=JointPriorPreviewMixture()$rhosigma$sigma, pi=JointPriorPreviewMixture()$pirhosigmaM),
                  mapping=aes(x=rho, y=sigma, fill=pi)) +
        scale_fill_viridis_c(option="turbo") + 
        labs(title=paste("Joint Base Prior. Cumulative Prob.=", round(JointPriorPreviewMixture()$probInt, digits=2)), 
             x=expression(rho), y=expression(sigma)) +
        theme_bw() + theme(plot.title=element_text(face='bold'))
      
      ggplotJointPriorPreview <- ggplotPCprior + ggplotENpriorAnalytic
      return(ggplotJointPriorPreview)
    }
    
    downloadFile(
      id = "MixturePreviewJointPriorPlot",
      logger = ss_userAction.Log,
      filenameroot = "MixturePreviewJointPriorPlot",
      aspectratio  = 1,
      downloadfxns = list(png  = MixturePreviewJointPriorPlot)
    )
    
    downloadablePlot(id = "MixturePreviewJointPriorPlot",
                     logger = ss_userAction.Log,
                     filenameroot = "MixturePreviewJointPriorPlot",
                     aspectratio  = 1,
                     downloadfxns = list(png=MixturePreviewJointPriorPlot),
                     visibleplot  = MixturePreviewJointPriorPlot)
    
    ## Mesh construction for MixtureModel ====
    
    MixtureMeshBase <- reactive({
      if(input$MixtureDataSimulatedLoaded=="sim"){
        DFsample <- as.data.frame(Mixture.sampling())
        convexhull <- chull(DFsample[,1:2])
        convexhull <- c(convexhull, convexhull[1])
        qloc <- quantile(as.vector(dist(DFsample[sample(1:nrow(DFsample),size=min(c(50,nrow(DFsample)))),1:2])),probs=c(0.03,0.3))
        mesh <- fm_mesh_2d_inla(loc.domain=DFsample[convexhull,1:2], cutoff = qloc[1]/2, offset=c(-0.1, -0.2),
                                max.edge=c(qloc[1], qloc[2]))
        sample <- DFsample
      } else if(input$MixtureDataSimulatedLoaded=="load"){
        DFsample <- datareadSample()
        if(input$MixtureRasterSPDE=="raster"){
          rasterSample <- datareadRaster()[sample(1:nrow(datareadRaster()), min(c(50,nrow(datareadRaster())))),1:2]
          qloc <- quantile(as.vector(dist(rasterSample)),probs=c(0.03,0.3))
          convexhull <- chull(datareadRaster()[,1:2])
          convexhull <- c(convexhull, convexhull[1])
          mesh <- fm_mesh_2d_inla(loc.domain=datareadRaster()[convexhull,1:2], cutoff = qloc[1]/2, offset=c(-0.1, -0.2),
                                  max.edge=c(qloc[1], qloc[2]))
          sample <- rasterSample
        } else if(input$MixtureRasterSPDE=="solvecov"){
          convexhull <- chull(rbind(DFsample[,1:2],datareadRaster()[,1:2]))
          convexhull <- c(convexhull, convexhull[1])
          qloc <- quantile(as.vector(dist(DFsample[sample(1:nrow(DFsample),size=min(c(50,nrow(DFsample)))),1:2])),probs=c(0.03,0.3))
          mesh <- fm_mesh_2d_inla(loc.domain=rbind(DFsample[,1:2],datareadRaster()[,1:2])[convexhull,], cutoff = qloc[1]/2, offset=c(-0.1, -0.2),
                                  max.edge=c(qloc[1], qloc[2]))
          sample <- DFsample
        }
      }
      result <- list(mesh=mesh, qloc=qloc, Sample=sample)
      return(result)
    })
    
    MixtureMesh <- eventReactive(input$buildMixtureMesh, {
      if(input$MixtureDataSimulatedLoaded=="sim"){
        DFsample <- as.data.frame(Mixture.sampling())
      } else if(input$MixtureDataSimulatedLoaded=="load"){
        DFsample <- as.data.frame(datareadSample())
        DFraster <- try(datareadRaster(), silent=TRUE)
        if(class(DFraster)!="try-error"){
          DFsample <- rbind(DFsample[,1:2], DFraster[,1:2])
        }
      }
      
      interiorNonConvex <- function(x, condition, convex, resolution, file.read){
        if(condition=="interiorMixtureMeshnonconvex"){
          interior <- inla.nonconvex.hull(points=x, convex=convex, resolution=resolution)
        } else if(condition=="interiorMixtureMeshcustomboundary"){
          ext <- tools::file_ext(file.read$datapath)
          validate(need(ext == c("csv","rds"), "Please upload a csv or rds file"))
          if(ext=="csv"){coords <- read.csv(file.read$datapath)}
          else coords <- readRDS(file.read$datapath)
          innerBorderMesh <- SpatialPolygons(Srl=list(Polygons(srl=list(Polygon(coords=coords)), ID="interiorBoundary")))
          interior <- inla.sp2segment(sp=innerBorderMesh)
        } else{interior <- NULL}
        return(interior)
      }
      
      boundaryNonConvex <- function(x, condition, convex, resolution, file.read){
        if(condition=="MixtureMeshnonconvex"){
          boundary <- inla.nonconvex.hull(points=x, convex=convex, resolution=resolution)
        } else if(condition=="MixtureMeshcustomboundary"){
          ext <- tools::file_ext(file.read$datapath)
          validate(need(ext == c("csv","rds"), "Please upload a csv or rds file"))
          if(ext=="csv"){coords <- read.csv(file.read$datapath)}
          else coords <- readRDS(file.read$datapath)
          outerBorderMesh <- SpatialPolygons(Srl=list(Polygons(srl=list(Polygon(coords=coords)), ID="externalBoundary")))
          boundary <- inla.sp2segment(sp=outerBorderMesh)
        } else{boundary <- NULL}
        return(boundary)
      }
      
      if(input$selectionMixtureMesh=="qlocation"){
        if(input$MixtureRasterSPDE=="raster"){
          rasterSample <- datareadRaster()[sample(1:nrow(datareadRaster()), min(c(50,nrow(datareadRaster())))),1:2]
          qloc <- quantile(as.vector(dist(rasterSample)),probs=c(ifelse(input$MixtureCustomMesh,input$MixtureMeshQloc,0.03),0.3))
          convexhull <- chull(DFsample[,1:2])
          convexhull <- c(convexhull, convexhull[1])
          mesh <- fm_mesh_2d_inla(loc.domain=DFsample[convexhull,1:2], cutoff = qloc[1]/2, offset=c(-0.1, -0.2), max.edge=c(qloc[1], qloc[2]),
                                  boundary=list(interiorNonConvex(x=as.matrix(rasterSample[,1:2]), condition=input$interiorMixtureMesh, convex=input$interiorcurvatureMixtureMesh, resolution=input$interiorresolutionMixtureMesh, file.read=input$interiorshapefileMixtureMesh),
                                                boundaryNonConvex(x=as.matrix(rasterSample[,1:2]), condition=input$boundaryMixtureMesh, convex=input$curvatureMixtureMesh, resolution=input$resolutionMixtureMesh, file.read=input$shapefileMixtureMesh)))
          sample <- rasterSample
        } else if(input$MixtureRasterSPDE=="solvecov"){
          qloc <- quantile(as.vector(dist(DFsample[sample(1:nrow(DFsample),size=min(c(50,nrow(DFsample)))),1:2])),probs=c(ifelse(input$MixtureCustomMesh,input$MixtureMeshQloc,0.03),0.3))
          convexhull <- chull(DFsample[,1:2])
          convexhull <- c(convexhull, convexhull[1])
          mesh <- fm_mesh_2d_inla(loc.domain=DFsample[convexhull,1:2], cutoff = qloc[1]/2, offset=c(-0.1, -0.2), max.edge=c(qloc[1], qloc[2]),
                                  boundary=list(interiorNonConvex(x=as.matrix(DFsample[,1:2]), condition=input$interiorMixtureMesh, convex=input$interiorcurvatureMixtureMesh, resolution=input$interiorresolutionMixtureMesh, file.read=input$interiorshapefileMixtureMesh),
                                                boundaryNonConvex(x=as.matrix(DFsample[,1:2]), condition=input$boundaryMixtureMesh, convex=input$curvatureMixtureMesh, resolution=input$resolutionMixtureMesh, file.read=input$shapefileMixtureMesh)))
          sample <- DFsample
        }
      } else if(input$selectionMixtureMesh=="edgelength"){
        if(input$MixtureRasterSPDE=="raster"){
          rasterSample <- datareadRaster()
          qloc <- input$EdgeLengthMixtureMesh
          convexhull <- chull(DFsample[,1:2])
          convexhull <- c(convexhull, convexhull[1])
          mesh <- fm_mesh_2d_inla(loc.domain=DFsample[convexhull,1:2], max.edge=c(1,1.5)*input$EdgeLengthMixtureMesh, cutoff=input$EdgeLengthMixtureMesh/5, offset=c(-0.1, -0.2), max.edge=c(qloc[1], qloc[2]),
                                  boundary=list(interiorNonConvex(x=as.matrix(rasterSample[,1:2]), condition=input$interiorMixtureMesh, convex=input$interiorcurvatureMixtureMesh, resolution=input$interiorresolutionMixtureMesh, file.read=input$interiorshapefileMixtureMesh),
                                                boundaryNonConvex(x=as.matrix(rasterSample[,1:2]), condition=input$boundaryMixtureMesh, convex=input$curvatureMixtureMesh, resolution=input$resolutionMixtureMesh, file.read=input$shapefileMixtureMesh)))
          sample <- rasterSample
        } else if(input$MixtureRasterSPDE=="solvecov"){
          qloc <- input$EdgeLengthMixtureMesh
          convexhull <- chull(DFsample[,1:2])
          convexhull <- c(convexhull, convexhull[1])
          mesh <- fm_mesh_2d_inla(loc.domain=DFsample[convexhull,1:2], max.edge=c(1,1.5)*input$EdgeLengthMixtureMesh,  cutoff=input$EdgeLengthMixtureMesh/5, offset=c(-0.1, -0.2),
                                  boundary=list(interiorNonConvex(x=as.matrix(DFsample[,1:2]), condition=input$interiorMixtureMesh, convex=input$interiorcurvatureMixtureMesh, resolution=input$interiorresolutionMixtureMesh, file.read=input$interiorshapefileMixtureMesh),
                                                boundaryNonConvex(x=as.matrix(DFsample[,1:2]), condition=input$boundaryMixtureMesh, convex=input$curvatureMixtureMesh, resolution=input$resolutionMixtureMesh, file.read=input$shapefileMixtureMesh)))
          sample <- DFsample
        }
      }
      
      result <- list(mesh=mesh, qloc=qloc, Sample=sample)
      return(result)
    })
    
    ggplotMixtureMesh <- function(){
      if(input$buildMixtureMesh==0){
        ggplot()+ gg(MixtureMeshBase()$mesh)+ theme_bw() + xlab("Latitude") + ylab("Longitude") +
          ggtitle("Mesh over the study region") + 
          geom_point(data=MixtureMeshBase()$Sample, 
                     aes(x=MixtureMeshBase()$Sample[,1],y=MixtureMeshBase()$Sample[,2]), size=1) + 
          theme(plot.title=element_text(color = "black", size = 16, face = "bold", hjust = 0.5))
      } else{
        ggplot()+ gg(MixtureMesh()$mesh)+ theme_bw() + xlab("Latitude") + ylab("Longitude") +
          ggtitle("Mesh over the study region") + 
          geom_point(data=MixtureMesh()$Sample, 
                     aes(x=MixtureMesh()$Sample[,1],y=MixtureMesh()$Sample[,2]), size=1) + 
          theme(plot.title=element_text(color = "black", size = 16, face = "bold", hjust = 0.5))
      }
      
    }
    
    downloadFile(
      id = "ggplotMixtureMesh",
      logger = ss_userAction.Log,
      filenameroot = "ggplotMixtureMesh",
      aspectratio  = 1,
      downloadfxns = list(png  = ggplotMixtureMesh)
    )
    
    downloadablePlot(
      id = "ggplotMixtureMesh",
      logger = ss_userAction.Log,
      filenameroot = "ggplotMixtureMesh",
      aspectratio  = 1,
      downloadfxns = list(png  = ggplotMixtureMesh),
      visibleplot  = ggplotMixtureMesh)
    
    ## Mixture modelling code ====
    
    MixtureModelFit <- eventReactive(input$fitMixture, {
      showNotification(ui=paste("Fitting the data."), duration = NULL)
      t1 <- Sys.time()
      #taking the data from simulation or from the loading tab
      if(input$MixtureDataSimulatedLoaded=="sim"){
        DFsample <- Mixture.sampling()$MixtureSample
      } else if(input$MixtureDataSimulatedLoaded=="load"){
        DFsample <- as.data.frame(datareadSample())
        if(input$MixtureRasterSPDE=="raster"){
          DFraster <- as.data.frame(datareadRaster())
        }
      }

      variablesChosenDefault <- input$DefaultComponentsMixture
      variablesChosenUser <- input$UserComponentsMixture
      variablesChosen <- c(variablesChosenDefault, variablesChosenUser)

      idMixtures <- unique(DFsample[,input$UserComponentsMixtureDependentGroup])[unique(DFsample[,input$UserComponentsMixtureDependentGroup])!="Ind"]
      nMixtures <- length(idMixtures)
      UserComponentsMixtureSharing <- list()
      for(i in seq_len(nMixtures)){
        UserComponentsMixtureSharing[[paste0("Mixture_",i)]] <- eval(parse(text=paste0("input$UserComponentsMixtureSharing_",i)))
      }
      
      if(input$buildMixtureMesh==0){
        server_mesh <- MixtureMeshBase()
        mesh <- server_mesh$mesh
      } else{
        server_mesh <- MixtureMesh()
        mesh <- server_mesh$mesh
      }
      
      if(input$optionMixtureS=="auto"){
        if(input$KindPriorSpatialEffectMixture=="PC.prior"){
          prior.range <- c(mean(c(diff(range(mesh$loc[mesh$segm$int$idx[,1],1])),diff(range(mesh$loc[mesh$segm$int$idx[,1],2]))))/5, 0.5)
          prior.sigma <- c(1,0.5)
          spde <- inla.spde2.pcmatern(mesh, prior.range = prior.range, prior.sigma = prior.sigma, alpha=2)
        } else{
          prior.range <- mean(c(diff(range(mesh$loc[mesh$segm$int$idx[,1],1])),diff(range(mesh$loc[mesh$segm$int$idx[,1],2]))))/5
          prior.sigma <- 1
          alpha <- 2; d <- 2
          nu <-  alpha - d/2
          kappa0 <-  log(8*nu)/2 -log(prior.range)
          tau0 <-  0.5*(lgamma(nu) - lgamma(nu + d/2) - d/2*log(4*pi)) - nu*kappa0 - log(prior.sigma)
          spde <-  inla.spde2.matern(mesh = mesh, B.tau = cbind(tau0, nu, -1), B.kappa = cbind(kappa0, -1, 0),
                                     theta.prior.mean = c(0.5,0.5), theta.prior.prec = c(1,1))
        }
      } else if(input$optionMixtureS=="custom"){
        if(input$KindPriorSpatialEffectMixture=="PC.prior"){
          prior.range <- as.numeric(unlist(strsplit(input$MixturePriorRangePC, split=",")))
          prior.sigma <- as.numeric(unlist(strsplit(input$MixturePriorStdevPC, split=",")))
          spde <- inla.spde2.pcmatern(mesh, prior.range = prior.range, prior.sigma = prior.sigma, alpha=2)
        } else{
          prior.range <- as.numeric(unlist(strsplit(input$MixturePriorRangeBase, split=",")))
          prior.sigma <- as.numeric(unlist(strsplit(input$MixturePriorStdevBase, split=",")))
          alpha <- 2; d <- 2
          nu <-  alpha - d/2
          kappa0 <-  log(8*nu)/2 -log(prior.range[1])
          tau0 <-  0.5*(lgamma(nu) - lgamma(nu + d/2) - d/2*log(4*pi)) - nu*kappa0 - log(prior.sigma[1])
          spde <-  inla.spde2.matern(mesh = mesh, B.tau = cbind(tau0, nu, -1), B.kappa = cbind(kappa0, -1, 0),
                                     theta.prior.mean = c(prior.range[2],prior.sigma[2]), theta.prior.prec = c(prior.range[3],prior.sigma[3]))
        }
      }
      
      spde.index <- inla.spde.make.index(name="Spatial", n.spde = spde$n.spde)
      
      # LGCP mesh operations
      ldomain <- unique(mesh$loc[mesh$segm$int$idx,1:2])
      dmesh <- mesh.dual(mesh = mesh)
      domain.polys <- Polygons(list(Polygon(ldomain)), '0')
      domainSP <- SpatialPolygons(list(domain.polys))
      w <- sapply(1:length(dmesh), function(i) {
        if (gIntersects(dmesh[i, ], domainSP))
          return(gArea(gIntersection(dmesh[i, ], domainSP)))
        else return(0)
      })
      
      n <- nrow(DFsample)
      nv <- mesh$n
      imat <- Diagonal(nv, rep(1, nv))
      lmat <- inla.spde.make.A(mesh, as.matrix(DFsample[,1:2]))
      A.inf <- rbind(imat, lmat) 
      
      ### Prediction of covariates ====
      
      List.covariates.inf <- list()
      List.covariates.mesh <- list()
      List.covariates.pred <- list()
      
      
      if((input$MixtureDataSimulatedLoaded=="load"&input$MixtureRasterSPDE=="solvecov")|input$MixtureDataSimulatedLoaded=="sim"|(input$MixtureDataSimulatedLoaded=="load"&input$MixtureRasterSPDE=="raster"|input$MixtureRasterPred=="SPDEraster")){
        x.pred <- seq(range(mesh$loc[mesh$segm$int$idx[,1], 1])[1], range(mesh$loc[mesh$segm$int$idx[,1], 1])[2], length.out=input$MixtureSPDErasterInput)
        y.pred <- seq(range(mesh$loc[mesh$segm$int$idx[,1], 2])[1], range(mesh$loc[mesh$segm$int$idx[,1], 2])[2], length.out=input$MixtureSPDErasterInput)
        xy.pred <- expand.grid(x=x.pred,y=y.pred)
        xy.pred <- xy.pred[which(!is.na(over(SpatialPoints(coords=xy.pred),SpatialPolygons(Srl=list(Polygons(srl=list(Polygon(coords=mesh$loc[mesh$segm$int$idx[,1], 1:2])), ID="int")))))),1:2]          
        A.geo.pred <- inla.spde.make.A(mesh=mesh, loc=as.matrix(xy.pred))
      } else{ #if(input$MixtureDataSimulatedLoaded=="load"&input$MixtureRasterSPDE=="raster"|input$MixtureRasterPred=="rasterpred"){
        xy.pred <- DFraster        
        A.geo.pred <- inla.spde.make.A(mesh=mesh, loc=as.matrix(xy.pred))
      }
      
      for(i in seq_along(variablesChosenUser)){
        if(!is.character(DFsample[,variablesChosenUser[i]])){
          prior.range.cov <- c(mean(c(diff(range(mesh$loc[mesh$segm$int$idx[,1],1])),diff(range(mesh$loc[mesh$segm$int$idx[,1],2]))))/5, 0.5)
          prior.sigma.cov <- c(1,0.5)
          spde.cov <- inla.spde2.pcmatern(mesh, prior.range = prior.range.cov, prior.sigma = prior.sigma.cov, alpha=2)
          spde.cov.index <- inla.spde.make.index(name="spatial.cov", n.spde = spde.cov$n.spde)
          formula.cov <- y ~ -1 + Intercept + f(spatial.cov, model=spde.cov)
          
          if((input$MixtureDataSimulatedLoaded=="load"&input$MixtureRasterSPDE=="solvecov")|input$MixtureDataSimulatedLoaded=="sim"){
            x.pred <- seq(range(mesh$loc[mesh$segm$int$idx[,1], 1])[1], range(mesh$loc[mesh$segm$int$idx[,1], 1])[2], length.out=input$MixtureSPDErasterInput)
            y.pred <- seq(range(mesh$loc[mesh$segm$int$idx[,1], 2])[1], range(mesh$loc[mesh$segm$int$idx[,1], 2])[2], length.out=input$MixtureSPDErasterInput)
            xy.pred <- expand.grid(x=x.pred,y=y.pred)
            xy.pred <- xy.pred[which(!is.na(over(SpatialPoints(coords=xy.pred),SpatialPolygons(Srl=list(Polygons(srl=list(Polygon(coords=mesh$loc[mesh$segm$int$idx[,1], 1:2])), ID="int")))))),1:2]          
            A.geo.pred <- inla.spde.make.A(mesh=mesh, loc=as.matrix(xy.pred))
            Inf.stack.cov <- inla.stack(data=list(y=DFsample[,variablesChosenUser[i]]),
                                        A=list(lmat, 1),
                                        effects=list(list(spatial.cov=spde.cov.index$spatial.cov),
                                                     list(Intercept=rep(1,nrow(DFsample)))
                                                     ),
                                        tag="Inference.cov")
            Pred.mesh.stack.cov <- inla.stack(data=list(y=rep(NA, mesh$n)),
                                              A=list(diag(1,nrow=mesh$n), 1),
                                              effects=list(list(spatial.cov=spde.cov.index$spatial.cov),
                                                           list(Intercept=rep(1,mesh$n))
                                              ),
                                              tag="Mesh.cov")
            Pred.stack.cov <- inla.stack(data=list(y=rep(NA,nrow(xy.pred))),
                                         A=list(A.geo.pred, 1),
                                         effects=list(list(spatial.cov=spde.cov.index$spatial.cov),
                                                      list(Intercept=rep(1,nrow(xy.pred)))
                                                      ),
                                         tag="Prediction.cov")
            Total.stack.cov <- inla.stack(Inf.stack.cov, Pred.mesh.stack.cov, Pred.stack.cov)
            mod.cov <- inla(formula=formula.cov, data=inla.stack.data(Total.stack.cov), family="gaussian",
                            control.predictor = list(compute=FALSE, A=inla.stack.A(Total.stack.cov)))
            indx.inf <- inla.stack.index(Total.stack.cov, tag="Inference.cov")$data
            indx.mesh <- inla.stack.index(Total.stack.cov, tag="Mesh.cov")$data
            indx.pred <- inla.stack.index(Total.stack.cov, tag="Prediction.cov")$data
            List.covariates.inf[[variablesChosenUser[i]]] <- as.numeric(!is.na(DFsample[,variablesChosenUser[i]]))*DFsample[,variablesChosenUser[i]] + as.numeric(is.na(DFsample[,variablesChosenUser[i]]))*mod.cov$summary.fitted.values[indx.inf,"mean"]
            List.covariates.mesh[[variablesChosenUser[i]]] <- mod.cov$summary.fitted.values[indx.mesh,"mean"]
            List.covariates.pred[[variablesChosenUser[i]]] <- mod.cov$summary.fitted.values[indx.pred,"mean"]
          } else if(input$MixtureDataSimulatedLoaded=="load"&input$MixtureRasterSPDE=="raster"|input$MixtureRasterPred=="SPDEraster"){
            x.pred <- seq(range(mesh$loc[mesh$segm$int$idx[,1], 1])[1], range(mesh$loc[mesh$segm$int$idx[,1], 1])[2], length.out=input$MixtureSPDErasterInput)
            y.pred <- seq(range(mesh$loc[mesh$segm$int$idx[,1], 2])[1], range(mesh$loc[mesh$segm$int$idx[,1], 2])[2], length.out=input$MixtureSPDErasterInput)
            xy.pred <- expand.grid(x=x.pred,y=y.pred)
            xy.pred <- xy.pred[which(!is.na(over(SpatialPoints(coords=xy.pred),SpatialPolygons(Srl=list(Polygons(srl=list(Polygon(coords=mesh$loc[mesh$segm$int$idx[,1], 1:2])), ID="int")))))),1:2]          
            A.geo.pred <- inla.spde.make.A(mesh=mesh, loc=as.matrix(xy.pred))
            if(variablesChosenUser[i] %in% colnames(DFraster)){
              Inf.stack.cov <- inla.stack(data=list(y=c(DFsample[,variablesChosenUser[i]], DFraster[!is.na(DFraster[,variablesChosenUser[i]]),variablesChosenUser[i]])),
                                          A=list(inla.spde.make.A(mesh=mesh, loc=as.matrix(rbind(DFsample[,1:2],DFraster[!is.na(DFraster[,variablesChosenUser[i]]),1:2]))), 1),
                                          effects=list(list(spatial.cov=spde.cov.index$spatial.cov),
                                                       list(Intercept=rep(1,nrow(DFsample)+nrow(DFraster)))
                                          ),
                                          tag="Inference.cov")
            } else{
              Inf.stack.cov <- inla.stack(data=list(y=DFsample[,variablesChosenUser[i]]),
                                          A=list(inla.spde.make.A(mesh=mesh, loc=as.matrix(DFsample[,1:2])), 1),
                                          effects=list(list(spatial.cov=spde.cov.index$spatial.cov),
                                                       list(Intercept=rep(1,nrow(DFsample)))
                                          ),
                                          tag="Inference.cov")
            }
            Pred.mesh.stack.cov <- inla.stack(data=list(y=rep(NA, mesh$n)),
                                              A=list(diag(1,nrow=mesh$n), 1),
                                              effects=list(list(spatial.cov=spde.cov.index$spatial.cov),
                                                           list(Intercept=rep(1,mesh$n))
                                              ),
                                              tag="Mesh.cov")
            Pred.stack.cov <- inla.stack(data=list(y=rep(NA,nrow(xy.pred))),
                                         A=list(A.geo.pred, 1),
                                         effects=list(list(spatial.cov=spde.cov.index$spatial.cov),
                                                      list(Intercept=rep(1,nrow(xy.pred)))
                                         ),
                                         tag="Prediction.cov")
            Total.stack.cov <- inla.stack(Inf.stack.cov, Pred.mesh.stack.cov, Pred.stack.cov)
            mod.cov <- inla(formula=formula.cov, data=inla.stack.data(Total.stack.cov), family="gaussian",
                            control.predictor = list(compute=FALSE, A=inla.stack.A(Total.stack.cov)))
            indx.inf <- inla.stack.index(Total.stack.cov, tag="Inference.cov")$data
            indx.mesh <- inla.stack.index(Total.stack.cov, tag="Mesh.cov")$data
            indx.pred <- inla.stack.index(Total.stack.cov, tag="Prediction.cov")$data
            List.covariates.inf[[variablesChosenUser[i]]] <- as.numeric(!is.na(DFsample[,variablesChosenUser[i]]))*DFsample[,variablesChosenUser[i]] + (as.numeric(is.na(DFsample[,variablesChosenUser[i]]))*(mod.cov$summary.fitted.values[indx.inf,"mean"])[1:nrow(DFsample)])
            List.covariates.mesh[[variablesChosenUser[i]]] <- mod.cov$summary.fitted.values[indx.mesh,"mean"]
            List.covariates.pred[[variablesChosenUser[i]]] <- mod.cov$summary.fitted.values[indx.pred,"mean"]
          } else if(input$MixtureDataSimulatedLoaded=="load"&input$MixtureRasterSPDE=="raster"|input$MixtureRasterPred=="rasterpred"){
            xy.pred <- DFraste[,1:2]          
            A.geo.pred <- inla.spde.make.A(mesh=mesh, loc=as.matrix(xy.pred))
            
            Inf.stack.cov <- inla.stack(data=list(y=c(DFsample[,variablesChosenUser[i]], DFraster[,variablesChosenUser[i]])),
                                        A=list(inla.spde.make.A(mesh=mesh, loc=as.matrix(rbind(DFsample[,1:2],DFraster[,1:2]))), 1),
                                        effects=list(list(spatial.cov=spde.cov.index$spatial.cov),
                                                     list(Intercept=rep(1,nrow(DFsample)+nrow(DFraster)))
                                        ),
                                        tag="Inference.cov")
            Pred.mesh.stack.cov <- inla.stack(data=list(y=rep(NA, mesh$n)),
                                              A=list(diag(1,nrow=mesh$n), 1),
                                              effects=list(list(spatial.cov=spde.cov.index$spatial.cov),
                                                           list(Intercept=rep(1,mesh$n))
                                              ),
                                              tag="Mesh.cov")
            Total.stack.cov <- inla.stack(Inf.stack.cov, Pred.mesh.stack.cov)
            mod.cov <- inla(formula=formula.cov, data=inla.stack.data(Total.stack.cov), family="gaussian",
                            control.predictor = list(compute=FALSE, A=inla.stack.A(Total.stack.cov)))
            indx.mesh <- inla.stack.index(Total.stack.cov, tag="Mesh.cov")$data
            List.covariates.mesh[[variablesChosenUser[i]]] <- mod.cov$summary.fitted.values[indx.mesh,"mean"]
            
            List.covariates.inf[[variablesChosenUser[i]]] <- DFsample[,variablesChosenUser[i]]
            List.covariates.pred[[variablesChosenUser[i]]] <- DFraster[,variablesChosenUser[i]]
            }
        } else{
          if((input$MixtureDataSimulatedLoaded=="load"&input$MixtureRasterSPDE=="solvecov")|input$MixtureDataSimulatedLoaded=="sim"){
            cov.inf.indx <- as.vector(unlist(lapply(X=1:nrow(DFsample), FUN=function(Y){which.min(apply(X=(DFsample[!is.na(DFsample[,variablesChosenUser[i]]),1:2]-matrix(DFsample[Y,1:2],ncol=2))**2, MARGIN=1, FUN=sum))})))
            List.covariates.inf[[variablesChosenUser[i]]] <- DFsample[cov.inf.indx, variablesChosenUser[i]]
            mesh.indx <- as.vector(unlist(lapply(X=1:nrow(mesh$loc), FUN=function(Y){which.min(apply(X=(DFsample[!is.na(DFsample[,variablesChosenUser[i]]),1:2]-matrix(mesh$loc[Y,1:2],ncol=2))**2, MARGIN=1, FUN=sum))})))
            List.covariates.mesh[[variablesChosenUser[i]]] <- DFsample[mesh.indx, variablesChosenUser[i]]
            
            MixtureFactorModelDF <- data.frame(y=1, Mixture=c(List.covariates.mesh[[variablesChosenUser[i]]], List.covariates.inf[[variablesChosenUser[i]]]))
            colnames(MixtureFactorModelDF) <- c("y", variablesChosenUser[i])
            idx.factor <- which(variablesChosenUser[i]==names(DFsample[!as.vector(unlist(lapply(X=DFsample, FUN=is.numeric)))])[names(DFsample[!as.vector(unlist(lapply(X=DFsample, FUN=is.numeric)))])%in%c(variablesChosenUser)])
            
            if(eval(parse(text=paste0("input$MixtureKindPredictionFactorLevel",idx.factor)))=="nearest"){
              cov.pred.ind <- as.vector(unlist(lapply(X=1:nrow(xy.pred), FUN=function(Y){which.min(apply(X=(DFsample[!is.na(DFsample[,variablesChosenUser[i]]),1:2]-matrix(xy.pred[Y,1:2],ncol=2))**2, MARGIN=1, FUN=sum))})))
              List.covariates.pred[[variablesChosenUser[i]]] <- DFsample[cov.pred.ind, variablesChosenUser[i]]
            } else{
              List.covariates.pred[[variablesChosenUser[i]]] <- rep(NA, times=nrow(xy.pred))
            }
            
          } else {
          DFsampleraster <- rbind(DFsample[,c(1:2,which(variablesChosenUser[i]==colnames(DFsample)))], DFraster[,c(1:2, which(variablesChosenUser[i]==colnames(DFraster)))])
          cov.inf.indx <- as.vector(unlist(lapply(X=1:nrow(DFsample), FUN=function(Y){which.min(apply(X=(DFsampleraster[!is.na(DFsampleraster[,variablesChosenUser[i]]),1:2]-matrix(DFsample[Y,1:2],ncol=2))**2, MARGIN=1, FUN=sum))})))
          List.covariates.inf[[variablesChosenUser[i]]] <- DFsample[cov.inf.indx, variablesChosenUser[i]]
          mesh.indx <- as.vector(unlist(lapply(X=1:nrow(mesh$loc), FUN=function(Y){which.min(apply(X=(DFsampleraster[!is.na(DFsampleraster[,variablesChosenUser[i]]),1:2]-matrix(mesh$loc[Y,1:2],ncol=2))**2, MARGIN=1, FUN=sum))})))
          List.covariates.mesh[[variablesChosenUser[i]]] <- DFsample[mesh.indx, variablesChosenUser[i]]
          
          MixtureFactorModelDF <- data.frame(y=1, Mixture=c(List.covariates.mesh[[variablesChosenUser[i]]], List.covariates.inf[[variablesChosenUser[i]]]))
          colnames(MixtureFactorModelDF) <- c("y", variablesChosenUser[i])
          idx.factor <- which(variablesChosenUser[i]==names(DFsample[!as.vector(unlist(lapply(X=DFsample, FUN=is.numeric)))])[names(DFsample[!as.vector(unlist(lapply(X=DFsample, FUN=is.numeric)))])%in%c(variablesChosenUser)])
          
          if(eval(parse(text=paste0("input$MixtureKindPredictionFactorLevel",idx.factor)))=="nearest"){
            cov.pred.ind <- as.vector(unlist(lapply(X=1:nrow(xy.pred), FUN=function(Y){which.min(apply(X=(DFsampleraster[!is.na(DFsampleraster[,variablesChosenUser[i]]),1:2]-matrix(xy.pred[Y,1:2],ncol=2))**2, MARGIN=1, FUN=sum))})))
            List.covariates.pred[[variablesChosenUser[i]]] <- DFsample[cov.pred.ind, variablesChosenUser[i]]
          } else{
            List.covariates.pred[[variablesChosenUser[i]]] <- rep(NA, times=nrow(xy.pred))
          }
          
        }
      }
      }
      
      ### Building the main model and stacks structure ====
      
      Inf.geo.effects.list <- list(
        list(),
        list()
      )
      
      Pred.geo.effects.list <- list(
        list(),
        list()
      )
      
      
      Inf.Mixtures.effects.list <- list()
      
      if(length(UserComponentsMixtureSharing)>0){
        for(i in seq_along(UserComponentsMixtureSharing)){
          Inf.Mixtures.effects.list[[paste0("Inf.", names(UserComponentsMixtureSharing)[i],".effects.list")]] <- list(list(),list()) 
        }
      } else{
        showNotification(ui=paste("No mixture group selected." ), duration = NULL)
      }
      
      
      formula_mod <- c("y ~ -1")
      
      A_Inf.spde1 <- list()
      A_Pred.spde1 <- list()
      
      for(i in seq_along(variablesChosen)){
        if(variablesChosen[i]=="Intercept"){
          formula_mod <- paste(formula_mod, "f(Intercept, model='linear')", sep=" + ")
          Inf.geo.effects.list[[2]][["Intercept"]] <- c(rep(1, times=nv), rep(1, times=n))
          Pred.geo.effects.list[[2]][["Intercept"]] <- rep(1, times=nrow(A.geo.pred))
          
          if(length(UserComponentsMixtureSharing)>0){
            for(k in seq_along(UserComponentsMixtureSharing)){
              if("Intercept" %in% UserComponentsMixtureSharing[[k]]){
                showNotification(ui=paste0("From a theoretical perspective, sharing linear effects is feasible but doesn't make logical sense. Consequently, we will approach them as non-shared effects."), duration=10, closeButton=TRUE, type="warning")
                ComponentMixtureGroup <- unique(DFsample[,input$UserComponentsMixtureDependentGroup])[unique(DFsample[,input$UserComponentsMixtureDependentGroup])!="Ind"]
                formula_mod <- paste(formula_mod, paste0("f(Intercept_",ComponentMixtureGroup[k],", model='linear')"), sep=" + ")
                Inf.Mixtures.effects.list[[paste0("Inf.", names(UserComponentsMixtureSharing)[k],".effects.list")]][[2]][[paste0("Intercept_", ComponentMixtureGroup[k])]] <- c(rep(1, times=nv), rep(1, times=n))
              } else{
                ComponentMixtureGroup <- unique(DFsample[,input$UserComponentsMixtureDependentGroup])[unique(DFsample[,input$UserComponentsMixtureDependentGroup])!="Ind"]
                formula_mod <- paste(formula_mod, paste0("f(Intercept_",ComponentMixtureGroup[k],", model='linear')"), sep=" + ")
                Inf.Mixtures.effects.list[[paste0("Inf.", names(UserComponentsMixtureSharing)[k],".effects.list")]][[2]][[paste0("Intercept_", ComponentMixtureGroup[k])]] <- c(rep(1, times=nv), rep(1, times=n))
              }
            }
          }
        } else if(variablesChosen[i]=="Spatial Effect"){
          formula_mod <- paste(formula_mod, "f(Spatial, model=spde)", sep=" + ")
          Inf.geo.effects.list[[1]][["Spatial"]] <- spde.index$Spatial
          Pred.geo.effects.list[[1]][["Spatial"]] <- spde.index$Spatial
          
          if(length(UserComponentsMixtureSharing)>0){
            for(k in seq_along(UserComponentsMixtureSharing)){
              if("Spatial Effect" %in% UserComponentsMixtureSharing[[k]]){
                ComponentMixtureGroup <- unique(DFsample[,input$UserComponentsMixtureDependentGroup])[unique(DFsample[,input$UserComponentsMixtureDependentGroup])!="Ind"]
                formula_mod <- paste(formula_mod, paste0("f(Spatial_",ComponentMixtureGroup[k],"_copy, copy='Spatial', fixed=FALSE)"), sep=" + ")
                Inf.Mixtures.effects.list[[paste0("Inf.", names(UserComponentsMixtureSharing)[k],".effects.list")]][[1]][[paste0("Spatial_", ComponentMixtureGroup[k], "_copy")]] <- spde.index$Spatial
              } else{
                ComponentMixtureGroup <- unique(DFsample[,input$UserComponentsMixtureDependentGroup])[unique(DFsample[,input$UserComponentsMixtureDependentGroup])!="Ind"]
                formula_mod <- paste(formula_mod, paste0("f(Spatial_",ComponentMixtureGroup[k],")"), sep=" + ")
                Inf.Mixtures.effects.list[[paste0("Inf.", names(UserComponentsMixtureSharing)[k],".effects.list")]][[1]][[paste0("Spatial_", ComponentMixtureGroup[k])]] <- spde.index$Spatial
              }
            }
          }
        } else{
          j <- which(variablesChosen[i]==variablesChosenUser)
          if(!is.character(DFsample[,variablesChosenUser[j]])){
            if(eval(parse(text=paste0("input$MixtureEffectCov",j)))=="rw1"|eval(parse(text=paste0("input$MixtureEffectCov",j)))=="rw2"){
              if(eval(parse(text=paste0("input$MixtureEffectCustomPrior",j)))=="custom"){
                if(eval(parse(text=paste0("input$MixtureEffectCovKindPrior",j)))=="pc"){
                  assign(paste0("pc.values", variablesChosenUser[j]), c( as.numeric(unlist(strsplit(eval(parse(text=paste0("input$MixtureEffectCovPriorPCValues",j))), ","))) ))
                  formula_mod <- paste(formula_mod, paste0("f(",variablesChosenUser[j],", model='",eval(parse(text=paste0("input$MixtureEffectCov",j))),"', hyper=list(prec=list(prior='pc.prec', param=", eval(parse(text=paste0("pc.values", variablesChosenUser[j]))), ")), constr=",eval(parse(text=paste0("input$MixtureEffectCovConstr",j))),", scale.model=TRUE)"), sep=" + ")
                } else if(eval(parse(text=paste0("input$MixtureEffectCovKindPrior",j)))=="unif"){
                  lim <- as.numeric(unlist(strsplit(eval(parse(text=paste0("input$MixtureEffectCovPriorUnif",j))), ",")))
                  sigma <- seq(0, lim[2]*3, length.out=1E5)
                  theta <- -2*log(sigma)
                  logdens <- sapply(X=sigma, FUN=logdunif, lim=lim)
                  unif.prior <- list(theta=list(prior=paste0("table: ", paste(c(theta, logdens), collapse=" "))))
                  assign(paste0("hyper", variablesChosenUser[j]), unif.prior )
                  formula_mod <- paste(formula_mod, paste0("f(",variablesChosenUser[j],", model='",eval(parse(text=paste0("input$MixtureEffectCov",j))),"', hyper=",paste0("hyper",variablesChosenUser[j]),  ", constr=",eval(parse(text=paste0("input$MixtureEffectCovConstr",j))),", scale.model=TRUE)"), sep=" + ")
                } else if(eval(parse(text=paste0("input$MixtureEffectCovKindPrior",j)))=="unifflat"){
                  assign(paste0("unifflat.prior",variablesChosenUser[j]), "expression:
                  log_dens = 0 - log(2) - theta/2
                  ")
                  formula_mod <- paste(formula_mod, paste0("f(",variablesChosenUser[j],", model='",eval(parse(text=paste0("input$MixtureEffectCov",j))),"', hyper=list(prec=list(prior=",paste0("unifflat.prior",variablesChosenUser[j]),")), constr=",eval(parse(text=paste0("input$MixtureEffectCovConstr",j))),", scale.model=TRUE)"), sep=" + ")
                } else{
                  assign(paste0("hyper",variablesChosenUser[j]), list(prec=list(prior="loggamma",param=c( as.numeric(unlist(strsplit(eval(parse(text=paste0("input$MixtureEffectCovPriorBaseValues",j))), ","))) ))) )
                  formula_mod <- paste(formula_mod, paste0("f(",variablesChosenUser[j],", model='",eval(parse(text=paste0("input$MixtureEffectCov",j))),"', hyper=",paste0("hyper",variablesChosenUser[j]),  ", constr=",eval(parse(text=paste0("input$MixtureEffectCovConstr",j))),", scale.model=TRUE)"), sep=" + ")
                }
              } else{
                formula_mod <- paste(formula_mod, paste0("f(",variablesChosenUser[j],", model='",eval(parse(text=paste0("input$MixtureEffectCov",j))),  "', constr=",eval(parse(text=paste0("input$MixtureEffectCovConstr",j))),", scale.model=TRUE)"), sep=" + ")
              }
              group.cov <- inla.group(x=c(List.covariates.mesh[[variablesChosenUser[j]]], List.covariates.inf[[variablesChosenUser[j]]], List.covariates.pred[[variablesChosenUser[j]]]), n=eval(parse(text=paste0("input$MixtureEffectCovNodes",j))), method="cut")
              
              Inf.geo.effects.list[[2]][[variablesChosenUser[j]]] <- group.cov[seq_len(n + nv)]
              Pred.geo.effects.list[[2]][[variablesChosenUser[j]]] <- group.cov[-seq_len(n + nv)]
              
              if(length(UserComponentsMixtureSharing)>0){
                for(k in seq_along(UserComponentsMixtureSharing)){
                  if(variablesChosenUser[j] %in% UserComponentsMixtureSharing[[k]]){
                    ComponentMixtureGroup <- unique(DFsample[,input$UserComponentsMixtureDependentGroup])[unique(DFsample[,input$UserComponentsMixtureDependentGroup])!="Ind"]
                    formula_mod <- paste(formula_mod, paste0("f(",variablesChosenUser[j],"_",ComponentMixtureGroup[k],"_copy, copy='",variablesChosenUser[j],"', fixed=FALSE)"), sep=" + ")
                    Inf.Mixtures.effects.list[[paste0("Inf.", names(UserComponentsMixtureSharing)[k],".effects.list")]][[2]][[paste0(variablesChosenUser[j], "_", ComponentMixtureGroup[k], "_copy")]] <- group.cov[seq_len(n + nv)]
                  } else{
                    ComponentMixtureGroup <- unique(DFsample[,input$UserComponentsMixtureDependentGroup])[unique(DFsample[,input$UserComponentsMixtureDependentGroup])!="Ind"]
                    formula_mod <- paste(formula_mod, paste0("f(",variablesChosenUser[j],"_",ComponentMixtureGroup[k],")"), sep=" + ")
                    Inf.Mixtures.effects.list[[paste0("Inf.", names(UserComponentsMixtureSharing)[k],".effects.list")]][[2]][[paste0(variablesChosenUser[j], "_", ComponentMixtureGroup[k])]] <- group.cov[seq_len(n + nv)]
                  }
                }
              }
              
            } else if(eval(parse(text=paste0("input$MixtureEffectCov",j)))=="spde1"){
              Tot_cov <- c(List.covariates.mesh[[variablesChosenUser[j]]], List.covariates.inf[[variablesChosenUser[j]]], List.covariates.pred[[variablesChosenUser[j]]])
              spde1_nodes <- seq(min(Tot_cov), max(Tot_cov), length.out=eval(parse(text=paste0("input$MixtureEffectCovNodes",j))))
              mesh1d <- fm_mesh_1d(loc=spde1_nodes)
              
              if(eval(parse(text=paste0("input$MixtureEffectCustomPrior",j)))=="custom"){
                if(eval(parse(text=paste0("input$MixtureEffectCovKindPrior",j)))=="pc"){
                  hyper <- as.numeric(unlist(strsplit(eval(parse(text=paste0("input$MixtureEffectCovPriorPCValues",i))), ",")))
                  assign(paste0("spde1d_",variablesChosenUser[j]) ,inla.spde2.pcmatern(mesh=mesh1d, prior.range=c(abs(diff(range(Tot_cov)))/5, 0.5), prior.sigma=c(1,0.5)))
                } else{
                  prior.range <- as.numeric(unlist(strsplit(eval(parse(text=paste0("input$MixtureEffectCovPriorBaseValues",i))), ",")))[1:3]
                  prior.sigma <- as.numeric(unlist(strsplit(eval(parse(text=paste0("input$MixtureEffectCovPriorBaseValues",i))), ",")))[4:6]
                  alpha <- 2; d <- 2
                  nu <-  alpha - d/2
                  kappa0 <-  log(8*nu)/2 -log(prior.range[1])
                  tau0 <-  0.5*(lgamma(nu) - lgamma(nu + d/2) - d/2*log(4*pi)) - nu*kappa0 - log(prior.sigma[1])
                  assign(paste0("spde1d_",variablesChosenUser[j]), inla.spde2.matern(mesh = mesh1d, B.tau = cbind(tau0, nu, -1), B.kappa = cbind(kappa0, -1, 0),
                                              theta.prior.mean = c(prior.range[2],prior.sigma[2]), theta.prior.prec = c(prior.range[3],prior.sigma[3])))
                }
              } else { 
                assign(paste0("spde1d_",variablesChosenUser[j]), inla.spde2.pcmatern(mesh=mesh1d, prior.range=c(abs(diff(range(Tot_cov)))/5, 0.5), prior.sigma=c(1,0.5)))
              }
              
              spde1d.index <- inla.spde.make.index(name=paste0(variablesChosenUser[j]), n.spde=eval(parse(text=paste0("spde1d_",variablesChosenUser[j])))$n.spde)
              formula_mod <- paste(formula_mod, paste0("f(",variablesChosenUser[j],", model=", paste0("spde1d_",variablesChosenUser[j]),  ")"), sep=" + ")

              Inf.geo.effects.list[[length(A_Inf.spde1)+3]] <- list()
              Pred.geo.effects.list[[length(A_Inf.spde1)+3]] <- list()
              Inf.geo.effects.list[[length(A_Inf.spde1)+3]][[variablesChosenUser[j]]] <- spde1d.index[[1]]
              Pred.geo.effects.list[[length(A_Inf.spde1)+3]][[variablesChosenUser[j]]] <- spde1d.index[[1]]
              
              if(length(UserComponentsMixtureSharing)>0){
                for(k in seq_along(UserComponentsMixtureSharing)){
                  if(variablesChosenUser[j] %in% UserComponentsMixtureSharing[[k]]){
                    ComponentMixtureGroup <- unique(DFsample[,input$UserComponentsMixtureDependentGroup])[unique(DFsample[,input$UserComponentsMixtureDependentGroup])!="Ind"]
                    formula_mod <- paste(formula_mod, paste0("f(",variablesChosenUser[j],"_",ComponentMixtureGroup[k],"_copy, copy='",variablesChosenUser[j],"', fixed=FALSE)"), sep=" + ")
                    Inf.Mixtures.effects.list[[paste0("Inf.", names(UserComponentsMixtureSharing)[k],".effects.list")]][[length(A_Inf.spde1)+3]] <- list()
                    Inf.Mixtures.effects.list[[paste0("Inf.", names(UserComponentsMixtureSharing)[k],".effects.list")]][[length(A_Inf.spde1)+3]][[paste0(variablesChosenUser[j], "_", ComponentMixtureGroup[k], "_copy")]] <- spde1d.index[[1]]
                  } else{
                    ComponentMixtureGroup <- unique(DFsample[,input$UserComponentsMixtureDependentGroup])[unique(DFsample[,input$UserComponentsMixtureDependentGroup])!="Ind"]
                    formula_mod <- paste(formula_mod, paste0("f(",variablesChosenUser[j],"_",ComponentMixtureGroup[k],", model=", paste0("spde1d_",variablesChosenUser[j]),")"), sep=" + ")
                    Inf.Mixtures.effects.list[[paste0("Inf.", names(UserComponentsMixtureSharing)[k],".effects.list")]][[length(A_Inf.spde1)+3]] <- list()
                    Inf.Mixtures.effects.list[[paste0("Inf.", names(UserComponentsMixtureSharing)[k],".effects.list")]][[length(A_Inf.spde1)+3]][[paste0(variablesChosenUser[j], "_", ComponentMixtureGroup[k])]] <- spde1d.index[[1]]
                  }
                }
              }

              A_Inf.spde1[[length(A_Inf.spde1)+1]] <- inla.spde.make.A(mesh=mesh1d, loc=Tot_cov[seq_len(nv+n)])
              A_Pred.spde1[[length(A_Pred.spde1)+1]] <- inla.spde.make.A(mesh=mesh1d, loc=Tot_cov[-seq_len(nv+n)])
              
            } else if(eval(parse(text=paste0("input$MixtureEffectCov",j)))=="linear"){
              if(eval(parse(text=paste0("input$MixtureEffectCustomPrior",j)))=="custom"){
                hyper <- as.numeric(unlist(strsplit(eval(parse(text=paste0("input$MixtureEffectCovPriorBaseValues",j))), ",")))
                formula_mod <- paste(formula_mod, paste0("f(",variablesChosenUser[j],", model='",eval(parse(text=paste0("input$MixtureEffectCov",j))),"', mean.linear=", hyper[1], ",prec.linear=", hyper[2],")"), sep=" + ")
              } else{
                formula_mod <- paste(formula_mod, paste0("f(",variablesChosenUser[j],", model='",eval(parse(text=paste0("input$MixtureEffectCov",j))),"')"), sep=" + ")
              }
              
              Inf.geo.effects.list[[2]][[variablesChosenUser[j]]] <- c(List.covariates.mesh[[variablesChosenUser[j]]], List.covariates.inf[[variablesChosenUser[j]]])
              Pred.geo.effects.list[[2]][[variablesChosenUser[j]]] <- c(List.covariates.pred[[variablesChosenUser[j]]])
              
              if(length(UserComponentsMixtureSharing)>0){
                for(k in seq_along(UserComponentsMixtureSharing)){
                  if(variablesChosenUser[j] %in% UserComponentsMixtureSharing[[k]]){
                    showNotification(ui=paste0("From a theoretical perspective, sharing linear effects is feasible but doesn't make logical sense. Consequently, we will approach them as non-shared effects."), duration=10, closeButton=TRUE, type="warning")
                    ComponentMixtureGroup <- unique(DFsample[,input$UserComponentsMixtureDependentGroup])[unique(DFsample[,input$UserComponentsMixtureDependentGroup])!="Ind"]
                    formula_mod <- paste(formula_mod, paste0("f(",variablesChosenUser[j],"_",ComponentMixtureGroup[k],", model='", eval(parse(text=paste0("input$MixtureEffectCov",j))),"')"), sep=" + ")
                    Inf.Mixtures.effects.list[[paste0("Inf.", names(UserComponentsMixtureSharing)[k],".effects.list")]][[2]][[paste0(variablesChosenUser[j], "_", ComponentMixtureGroup[k])]] <- c(List.covariates.mesh[[variablesChosenUser[j]]], List.covariates.inf[[variablesChosenUser[j]]])
                  } else{
                    ComponentMixtureGroup <- unique(DFsample[,input$UserComponentsMixtureDependentGroup])[unique(DFsample[,input$UserComponentsMixtureDependentGroup])!="Ind"]
                    formula_mod <- paste(formula_mod, paste0("f(",variablesChosenUser[j],"_",ComponentMixtureGroup[k],", model='", eval(parse(text=paste0("input$MixtureEffectCov",j))),"')"), sep=" + ")
                    Inf.Mixtures.effects.list[[paste0("Inf.", names(UserComponentsMixtureSharing)[k],".effects.list")]][[2]][[paste0(variablesChosenUser[j], "_", ComponentMixtureGroup[k])]] <- c(List.covariates.mesh[[variablesChosenUser[j]]], List.covariates.inf[[variablesChosenUser[j]]])
                  }
                }
              }
            } else{
              showNotification(ui=paste("The effect of numerical covariates cannot possess an independent and identically distributed (iid) structure. If this is required, the variable values should be recorded as text, not as numerical input.."), duration = NULL)
            }
          } else{ # Section for factor variables
            if(eval(parse(text=paste0("input$MixtureEffectCov",j)))=="iid"){
              MixtureFactorModelDF <- data.frame(y=1, Mixture=c(List.covariates.mesh[[variablesChosenUser[j]]], List.covariates.inf[[variablesChosenUser[j]]]))
              colnames(MixtureFactorModelDF) <- c("y", variablesChosenUser[j])
              idx.factor <- which(variablesChosenUser[j]==names(DFsample[!as.vector(unlist(lapply(X=DFsample, FUN=is.numeric)))])[names(DFsample[!as.vector(unlist(lapply(X=DFsample, FUN=is.numeric)))])%in%c(variablesChosenUser)])
              
              Inf.geo.effects.list[[2]][[variablesChosenUser[j]]] <- c(List.covariates.mesh[[variablesChosenUser[j]]], List.covariates.inf[[variablesChosenUser[j]]])
              
              Pred.geo.effects.list[[2]][[variablesChosenUser[j]]] <- if(eval(parse(text=paste0("input$MixtureKindPredictionFactorLevel",idx.factor)))=="reference"){
                rep( eval(parse(text=paste0("input$MixtureEffectCovFactorPred",idx.factor))), times=length(List.covariates.pred[[variablesChosenUser[j]]]))
              } else{List.covariates.pred[[variablesChosenUser[j]]]}
              
              #c(List.covariates.pred[[variablesChosenUser[j]]])
              
              if(eval(parse(text=paste0("input$MixtureEffectCustomPrior",j)))=="custom"){
                if(eval(parse(text=paste0("input$MixtureEffectCovKindPrior",j)))=="pc"){
                  hyper <- list(prec=list(prior="pc.prior",param=c( as.numeric(unlist(strsplit(eval(parse(text=paste0("input$MixtureEffectCovPriorPCValues",j))), ","))) )))
                  formula_mod <- paste(formula_mod, paste0("f(",variablesChosenUser[j],", model='",eval(parse(text=paste0("input$MixtureEffectCov",j))),"', hyper=",hyper, ", constr=",eval(parse(text=paste0("input$MixtureEffectCovConstr",j))),")"), sep=" + ")
                } else if(eval(parse(text=paste0("input$MixtureEffectCovKindPrior",j)))=="base"){
                  hyper <- list(prec=list(prior="loggamma",param=c( as.numeric(unlist(strsplit(eval(parse(text=paste0("input$MixtureEffectCovPriorBaseValues",j))), ","))) )))
                  formula_mod <- paste(formula_mod, paste0("f(",variablesChosenUser[j],", model='iid', hyper=",hyper,  ", constr=",eval(parse(text=paste0("input$MixtureEffectCovConstr",j))),")"), sep=" + ")
                } else if(eval(parse(text=paste0("input$MixtureEffectCovKindPrior",j)))=="unif"){
                  lim <- as.numeric(unlist(strsplit(eval(parse(text=paste0("input$MixtureEffectCovPriorUnif",j))), ",")))
                  sigma <- seq(lim[1]+1E-5, lim[2]*3, length.out=1E5)
                  theta <- -2*log(sigma)
                  logdens <- sapply(X=sigma, FUN=logdunif, lim=lim)
                  unif.prior <- list(theta=list(prior=paste0("table: ", paste(c(theta, logdens), collapse=" "))))
                  assign(paste0("hyper", variablesChosenUser[j]), unif.prior )
                  formula_mod <- paste(formula_mod, paste0("f(",variablesChosenUser[j],", model='iid', hyper=",paste0("hyper",variablesChosenUser[j]),  ", constr=",eval(parse(text=paste0("input$MixtureEffectCovConstr",j))),", scale.model=TRUE)"), sep=" + ")
                } else{ # flatunif
                  unifflat.prior= "expression:
                  log_dens = 0 - log(2) - theta/2
                  "
                  formula_mod <- paste(formula_mod, paste0("f(",variablesChosenUser[j],", model='iid', hyper=list(prec=list(prior=",unifflat.prior,")), constr=",eval(parse(text=paste0("input$MixtureEffectCovConstr",j))),")"), sep=" + ")
                }
              } else{
                formula_mod <- paste(formula_mod, paste0("f(", variablesChosenUser[j], ", model='iid'",  ", constr=",eval(parse(text=paste0("input$MixtureEffectCovConstr",j))),")"), sep=" + ") 
              }
              
              if(length(UserComponentsMixtureSharing)>0){
                for(k in seq_along(UserComponentsMixtureSharing)){
                  if(variablesChosenUser[j] %in% UserComponentsMixtureSharing[[k]]){
                    ComponentMixtureGroup <- unique(DFsample[,input$UserComponentsMixtureDependentGroup])[unique(DFsample[,input$UserComponentsMixtureDependentGroup])!="Ind"]
                    formula_mod <- paste(formula_mod, paste0("f(",variablesChosenUser[j],"_",ComponentMixtureGroup[k],"_copy, copy='",variablesChosenUser[j],"')"), sep=" + ")
                    Inf.Mixtures.effects.list[[paste0("Inf.", names(UserComponentsMixtureSharing)[k],".effects.list")]][[2]][[paste0(variablesChosenUser[j], "_", ComponentMixtureGroup[k], "_copy")]] <- c(List.covariates.mesh[[variablesChosenUser[j]]], List.covariates.inf[[variablesChosenUser[j]]])
                  } else{
                    ComponentMixtureGroup <- unique(DFsample[,input$UserComponentsMixtureDependentGroup])[unique(DFsample[,input$UserComponentsMixtureDependentGroup])!="Ind"]
                    formula_mod <- paste(formula_mod, paste0("f(",variablesChosenUser[j],"_",ComponentMixtureGroup[k],")"), sep=" + ")
                    Inf.Mixtures.effects.list[[paste0("Inf.", names(UserComponentsMixtureSharing)[k],".effects.list")]][[2]][[paste0(variablesChosenUser[j], "_", ComponentMixtureGroup[k])]] <- c(List.covariates.mesh[[variablesChosenUser[j]]], List.covariates.inf[[variablesChosenUser[j]]])
                  }
                }
              }
              
            } else{
              MixtureFactorModelDF <- data.frame(y=1, Mixture=c(List.covariates.mesh[[variablesChosenUser[j]]], List.covariates.inf[[variablesChosenUser[j]]]))
              colnames(MixtureFactorModelDF) <- c("y", variablesChosenUser[j])
              FactorVariables <- data.frame( model.matrix(object=as.formula(paste0("y~-1+", variablesChosenUser[j])), MixtureFactorModelDF ) )
              idx.factor <- which(variablesChosenUser[j]==names(DFsample[!as.vector(unlist(lapply(X=DFsample, FUN=is.numeric)))])[names(DFsample[!as.vector(unlist(lapply(X=DFsample, FUN=is.numeric)))])%in%c(variablesChosenUser)])
              for(l in setdiff(colnames(FactorVariables),paste0(variablesChosenUser[j], eval(parse(text=paste0("input$MixtureEffectCovFactorPred",idx.factor))) ))){
                ll <- l
                Inf.geo.effects.list[[2]][[l]] <- FactorVariables[,l]
                Pred.geo.effects.list[[2]][[l]] <- rep(0, length(List.covariates.pred[[variablesChosenUser[j]]]))
                if(eval(parse(text=paste0("input$MixtureEffectCustomPrior",j)))=="custom"|eval(parse(text=paste0("input$MixtureEffectCov",j)))=="linear"){
                  hyper <- as.numeric(unlist(strsplit(eval(parse(text=paste0("input$MixtureEffectCovPriorBaseValues",j))), ",")))
                  formula_mod <- paste(formula_mod, paste0("f(", l, ", model='linear', mean.linear=", hyper[1], ",prec.linear=", hyper[2],")"), sep=" + ")
                } else{
                  formula_mod <- paste(formula_mod, paste0("f(", l, ", model='linear')"), sep=" + ") 
                }
              }
              
              if(length(UserComponentsMixtureSharing)>0){
                for(k in seq_along(UserComponentsMixtureSharing)){
                  for(l in setdiff(colnames(FactorVariables),paste0(variablesChosenUser[j], eval(parse(text=paste0("input$MixtureEffectCovFactorPred",idx.factor))) ))){
                    if(variablesChosenUser[j] %in% UserComponentsMixtureSharing[[k]]){
                      showNotification(ui=paste0("From a theoretical perspective, sharing linear effects is feasible but doesn't make logical sense. Consequently, we will approach them as non-shared effects."), duration=10, closeButton=TRUE, type="warning")
                      ComponentMixtureGroup <- unique(DFsample[,input$UserComponentsMixtureDependentGroup])[unique(DFsample[,input$UserComponentsMixtureDependentGroup])!="Ind"]
                      formula_mod <- paste(formula_mod, paste0("f(",l,"_",ComponentMixtureGroup[k],")"), sep=" + ")
                      Inf.Mixtures.effects.list[[paste0("Inf.", names(UserComponentsMixtureSharing)[k],".effects.list")]][[2]][[paste0(l, "_", ComponentMixtureGroup[k])]] <- c(List.covariates.mesh[[variablesChosenUser[j]]], List.covariates.inf[[variablesChosenUser[j]]])
                    } else{
                      ComponentMixtureGroup <- unique(DFsample[,input$UserComponentsMixtureDependentGroup])[unique(DFsample[,input$UserComponentsMixtureDependentGroup])!="Ind"]
                      formula_mod <- paste(formula_mod, paste0("f(",l,"_",ComponentMixtureGroup[k],")"), sep=" + ")
                      Inf.Mixtures.effects.list[[paste0("Inf.", names(UserComponentsMixtureSharing)[k],".effects.list")]][[2]][[paste0(l, "_", ComponentMixtureGroup[k])]] <- c(List.covariates.mesh[[variablesChosenUser[j]]], List.covariates.inf[[variablesChosenUser[j]]])
                    }
                  }
                }
              }
              
              
            }
          }
        }
      }
      
      ### Stacks of the geostatistical, mixtures and prediction layers ====
      
      ResponseVariable <- DFsample[,3]
      
      A_inf_tot <- c(A.inf,1)
      if(length(A_Inf.spde1)>0){
        for(i in seq_along(A_Inf.spde1)){
          A_inf_tot[[2+i]] <- A_Inf.spde1[[i]] 
        }
      } 
      
      A_pred_tot <- c(A.geo.pred,1)
      if(length(A_Inf.spde1)>0){
        for(i in seq_along(A_Inf.spde1)){
          A_pred_tot[[2+i]] <- A_Pred.spde1[[i]] 
        }
      } 
      
      Inf.geo.stack <- inla.stack(data=list(y=cbind(c(rep(NA, nv),ResponseVariable), NA), e=rep(0,times=nv+n)),
                                  A=A_inf_tot,
                                  effects=Inf.geo.effects.list,
                                  tag="Inference_geo"
                                  )
      
      Pred.geo.stack <- inla.stack(data=list(y=matrix(NA, nrow=nrow(A.geo.pred), ncol=2)),
                                   A=A_pred_tot,
                                   effects=Pred.geo.effects.list,
                                   tag="Prediction_geo")
      
      
      if(length(UserComponentsMixtureSharing)>0){
        Total.stack <- inla.stack(Inf.geo.stack)
        ComponentMixtureGroup <- unique(DFsample[,input$UserComponentsMixtureDependentGroup])[unique(DFsample[,input$UserComponentsMixtureDependentGroup])!="Ind"]
        for(k in seq_along(UserComponentsMixtureSharing)){
          y.na.vec <- rep(NA, times=n)
          y.na.vec[ComponentMixtureGroup[k]==DFsample[,input$UserComponentsMixtureDependentGroup]] <- 1
          y.pp <- c(rep(0,times=nv), y.na.vec)
          assign(paste0("Inf.Mixture",k,".stack"), 
                 inla.stack(data=list(y=cbind(NA, y.pp), e=c(w, rep(0,n))),
                            A=A_inf_tot,
                            effects=Inf.Mixtures.effects.list[[paste0("Inf.", names(UserComponentsMixtureSharing)[k],".effects.list")]],
                            tag=paste0("Inf_Mixture",k)
                            )
                 )
          Total.stack <- inla.stack(Total.stack, eval(parse(text=paste0("Inf.Mixture",k,".stack"))))
        }
        Total.stack <- inla.stack(Total.stack, Pred.geo.stack)
      } else{
        Total.stack <- inla.stack(Inf.geo.stack, Pred.geo.stack)
      }
      
      ### INLA model ====
      
      formula_inla <- as.formula(formula_mod)

      if(input$autocustomMixtureFamily=='custom'){
        if(input$MixtureFamilyPriorKind=="pc"){
          family.pcprec <- as.numeric(unlist(strsplit(input$MixtureFamilyHyper,",")))
          controlFamily <- list(list(hyper = list(prec = list(prior="pc.prec", param=family.pcprec))), list())
        } else if(input$MixtureFamilyPriorKind=="unif"){
          lim <- as.numeric(unlist(strsplit(input$MixtureFamilyHyper,",")))
          sigma <- seq(0, lim[2]*3, length.out=1E5)
          theta <- -2*log(sigma)
          logdens <- sapply(X=sigma, FUN=logdunif, lim=lim)
          family.unif.prior <- list(theta=list(prior=paste0("table: ", paste(c(theta, logdens), collapse=" "))))
          controlFamily <- list(list(hyper = list(theta = list(prior=family.unif.prior))), list())
        } else if(input$MixtureFamilyPriorKind=="unifflat"){
          unifflat.prior= "expression:
                  log_dens = 0 - log(2) - theta/2
                  "
          controlFamily <- list(list(hyper = list(prec = list(prior=unifflat.prior))), list())
        } else{
          controlFamily <- list(list(hyper = list(prec = list(prior="loggamma", param=as.numeric(unlist(strsplit(input$MixtureFamilyHyper,",")))))), list())
        }
      } else{
        controlFamily <- list(list(), list())
      }

      if(input$autocustomMixtureMode=='custom'){
        controlModeTheta <- list(theta=as.numeric(unlist(strsplit(input$PrefModeHyper,","))), restart=TRUE)
      } else{controlModeTheta <- inla.set.control.mode.default()}
      
      if(input$INLAModeMixture=="classic"){
        controlINLA <- list(strategy=input$strategyapproxINLAMixture,
                            int.strategy=input$strategyintINLAMixture)
      } else{
        controlINLA <- list()
      }

      Mixture.model <- inla(formula=formula_inla, family = c(input$SelectMixtureFamily,'poisson'),
                            data = inla.stack.data(Total.stack),
                            E = inla.stack.data(Total.stack)$e,
                            control.inla = controlINLA,
                            control.predictor = list(A = inla.stack.A(Total.stack), compute = TRUE, link = 1),
                            control.family = controlFamily,
                            control.mode = controlModeTheta,
                            control.compute = list(cpo = TRUE, dic = TRUE, waic = TRUE, config = TRUE),
                            inla.mode=input$INLAModeMixture,
                            verbose=FALSE)

      index.pred <- inla.stack.index(Total.stack, "Prediction_geo")$data
      DFpred <- data.frame(Latitude=xy.pred[,1], Longitude=xy.pred[,2])
      DFpred$Abundance.mean <- Mixture.model$summary.fitted.values[index.pred, "mean"]
      DFpred$Abundance.median <- Mixture.model$summary.fitted.values[index.pred, "0.5quant"]
      DFpred$Abundance.sd <- Mixture.model$summary.fitted.values[index.pred, "sd"]

      colnames(DFpred)[1:2] <- colnames(DFsample)[1:2]
      
      DFpredictorMeanMedianStdev <- data.frame(Latitude=xy.pred[,1], Longitude=xy.pred[,2],
                                               Predictor.mean=Mixture.model$summary.linear.predictor[index.pred, "mean"],
                                               Predictor.median=Mixture.model$summary.linear.predictor[index.pred, "0.5quant"],
                                               Predictor.stdev=Mixture.model$summary.linear.predictor[index.pred, "sd"])
      
      colnames(DFpredictorMeanMedianStdev)[1:2] <- colnames(DFsample)[1:2]
      
      gridSpatial <- expand.grid(x=seq(min(DFpred[,1]), max(DFpred[,1]),length.out=input$dimMixturemap),
                                 y=seq(min(DFpred[,2]), max(DFpred[,2]), length.out=input$dimMixturemap))
      gridSpatial <- gridSpatial[which(!is.na(over(SpatialPoints(coords=gridSpatial),SpatialPolygons(Srl=list(Polygons(srl=list(Polygon(coords=mesh$loc[mesh$segm$int$idx[,1], 1:2])), ID="int")))))),1:2]          
      
      A.spatial <- inla.spde.make.A(mesh=mesh, loc=as.matrix(gridSpatial))

      DFspatialMeanMedianStdev <- data.frame(Latitude=as.vector(gridSpatial[,1]), Longitude=as.vector(gridSpatial[,2]),
                                             Spatial.mean=as.vector(A.spatial%*%Mixture.model$summary.random$Spatial$mean),
                                             Spatial.median=as.vector(A.spatial%*%Mixture.model$summary.random$Spatial$`0.5quant`),
                                             Spatial.stdev=as.vector(A.spatial%*%Mixture.model$summary.random$Spatial$sd))

      colnames(DFspatialMeanMedianStdev)[1:2] <- colnames(DFsample)[1:2]
      
      result <- list(DFpredAbunMeanMedianStdev=list(DFpredAbunMeanMedianStdev=DFpred),
                     DFpredictorMeanMedianStdev=list(DFpredictorMeanMedianStdev=DFpredictorMeanMedianStdev),
                     DFspatialMeanMedianStdev=list(DFspatialMeanMedianStdev=DFspatialMeanMedianStdev),
                     MixtureModel=list(MixtureModel=Mixture.model),
                     DFPostFixed=list(DFPostFixed=Mixture.model$marginals.fixed),
                     DFPostHyperpar=list(DFPostHyperpar=Mixture.model$marginals.hyperpar),
                     Summary.fixed=list(Summary.fixed=Mixture.model$summary.fixed),
                     Summary.hyperpar=list(Summary.hyperpar=Mixture.model$summary.hyperpar),
                     SummaryInternalHyper=list(SummaryInternalHyper=Mixture.model$internal.summary.hyperpar),
                     SummaryCPO=list(SummaryCPO=na.omit(Mixture.model$cpo$cpo)),
                     DICmodel=list(DICmodel=data.frame(DIC=Mixture.model$dic$family.dic, row.names=if(length(UserComponentsMixtureSharing)>0){c("Geostatistical", "Point process")}else{c("Geostatistical")})),
                     WAICmodel=list(WAICmodel=data.frame(WAIC=cbind(unlist(lapply(na.omit(unique(Mixture.model$dic$family)), function(i){sum(Mixture.model$waic$local.waic[which(Mixture.model$dic$family==i)])}))), row.names=if(length(UserComponentsMixtureSharing)>0){c("Geostatistical", "Point process")}else{c("Geostatistical")}))
                     )

      t2 <- Sys.time()
      difftime(t2,t1, units="secs")
      showNotification(ui=paste("The model has been fitted:", as.numeric(round(Mixture.model$cpu.used[4])),
                                "(Mixture model) and", as.numeric(round(difftime(t2,t1, units="secs"))),
                                "(overall process) secs." ), duration = NULL)
      # showNotification(ui=paste("The model's DIC is", Preferential.model$dic$dic), duration = NULL)
      return(result)
    })

    
    dataggplotMixtureAbundanceMeanMedianStdevFit <- function(){
      DF <- MixtureModelFit()$DFpredAbunMeanMedianStdev$DFpredAbunMeanMedianStdev
      return(DF)
    }
    
    DFMixtureAbundance <- reactive({
      DF <- MixtureModelFit()$DFpredAbunMeanMedianStdev$DFpredAbunMeanMedianStdev
      condition <- !(length(unique(DF[,1]))*length(unique(DF[,2]))==nrow(DF))&val()
      if(condition){
        grid <- expand.grid(seq(min(DF[,1]),max(DF[,1]), length.out=200),seq(min(DF[,2]),max(DF[,2]), length.out=200))
        DFInter <- data.frame(Latitude=grid[,1],Longitude=grid[,2])
        DFInter$Abundance.mean <- InterpolateIrrGrid(z=DF$Abundance.mean,loc=DF[,1:2], gridInter=grid)$DataInter$z
        DFInter$Abundance.median <- InterpolateIrrGrid(z=DF$Abundance.median,loc=DF[,1:2], gridInter=grid)$DataInter$z
        DFInter$Abundance.sd <- InterpolateIrrGrid(z=DF$Abundance.sd,loc=DF[,1:2], gridInter=grid)$DataInter$z
        colnames(DFInter)[1:2] <- colnames(DF)[1:2]
        DF <- DFInter
      }
      
      return(DF)
    })
    
    ggplotMixtureAbundanceMeanMedianStdevFit <- function(){
      DF <- DFMixtureAbundance()
      g1 <- ggplot(DF) + geom_tile(aes(x=DF[,1], y=DF[,2], fill=Abundance.mean)) +
        scale_fill_viridis_c(option = "turbo") + theme_bw() + coord_fixed() + labs(fill="Values") +
        xlab(colnames(DF)[1]) + ylab(colnames(DF)[1]) + ggtitle("Response predicted (mean)") +
        theme(plot.title=element_text(color = "black", size = 16, face = "bold", hjust = 0.5))
      
      g2 <- ggplot(DF) + geom_tile(aes(x=DF[,1], y=DF[,2], fill=Abundance.median)) +
        scale_fill_viridis_c(option = "turbo") + theme_bw() + coord_fixed() + labs(fill="Values") +
        xlab(colnames(DF)[1]) + ylab(colnames(DF)[1]) + ggtitle("Response predicted (median)") +
        theme(plot.title=element_text(color = "black", size = 16, face = "bold", hjust = 0.5))
      
      g3 <- ggplot(DF) + geom_tile(aes(x=DF[,1], y=DF[,2], fill=Abundance.sd))+
        scale_fill_viridis_c(option = "turbo") + theme_bw() + coord_fixed() + labs(fill="Values") +
        xlab(colnames(DF)[1]) + ylab(colnames(DF)[1]) + ggtitle("Response predicted (stdev.)") +
        theme(plot.title=element_text(color = "black", size = 16, face = "bold", hjust = 0.5))
      
      gt <- g1+g2+g3
      
      return(gt)
    }
    
    downloadFile(
      id = "ggplotMixtureAbundanceMeanMedianStdevFit",
      logger = ss_userAction.Log,
      filenameroot = "ggplotMixtureAbundanceMeanMedianStdevFit",
      aspectratio  = 1,
      downloadfxns = list(png  = ggplotMixtureAbundanceMeanMedianStdevFit,
                          csv = dataggplotMixtureAbundanceMeanMedianStdevFit,
                          txt = dataggplotMixtureAbundanceMeanMedianStdevFit)
    )
    
    downloadablePlot(id = "ggplotMixtureAbundanceMeanMedianStdevFit",
                     logger = ss_userAction.Log,
                     filenameroot = "ggplotMixtureAbundanceMeanMedianStdevFit",
                     aspectratio  = 1,
                     downloadfxns = list(png = ggplotMixtureAbundanceMeanMedianStdevFit,
                                         csv = dataggplotMixtureAbundanceMeanMedianStdevFit,
                                         txt = dataggplotMixtureAbundanceMeanMedianStdevFit),
                     visibleplot  = ggplotMixtureAbundanceMeanMedianStdevFit)
    
    
    dataggplotMixtureSpatialMeanMedianStdev <- function(){
      DF <- MixtureModelFit()$DFspatialMeanMedianStdev$DFspatialMeanMedianStdev
      return(DF)
    }
    
    ggplotMixtureSpatialMeanMedianStdev <- function(){
      DF <- MixtureModelFit()$DFspatialMeanMedianStdev$DFspatialMeanMedianStdev
      g1 <- ggplot(DF) + geom_tile(aes(x=DF[,1], y=DF[,2], fill=Spatial.mean)) +
        scale_fill_viridis_c(option = "turbo") + theme_bw() + coord_fixed() + labs(fill="Values") +
        xlab(colnames(DF)[1]) + ylab(colnames(DF)[1]) + ggtitle("Posterior spatial (mean)") +
        theme(plot.title=element_text(color = "black", size = 16, face = "bold", hjust = 0.5))
      
      g2 <- ggplot(DF) + geom_tile(aes(x=DF[,1], y=DF[,2], fill=Spatial.median)) +
        scale_fill_viridis_c(option = "turbo") + theme_bw() + coord_fixed() + labs(fill="Values") +
        xlab(colnames(DF)[1]) + ylab(colnames(DF)[1]) + ggtitle("Posterior spatial (median)") +
        theme(plot.title=element_text(color = "black", size = 16, face = "bold", hjust = 0.5))
      
      g3 <- ggplot(DF) + geom_tile(aes(x=DF[,1], y=DF[,2], fill=Spatial.stdev))+
        scale_fill_viridis_c(option = "turbo") + theme_bw() + coord_fixed() + labs(fill="Values") +
        xlab(colnames(DF)[1]) + ylab(colnames(DF)[1]) + ggtitle("Posterior spatial (stdev.)") +
        theme(plot.title=element_text(color = "black", size = 16, face = "bold", hjust = 0.5))
      
      gt <- g1+g2+g3
      
      return(gt)
    }
    
    downloadFile(
      id = "ggplotMixtureSpatialMeanMedianStdev",
      logger = ss_userAction.Log,
      filenameroot = "ggplotMixtureSpatialMeanMedianStdev",
      aspectratio  = 1,
      downloadfxns = list(png  = ggplotMixtureSpatialMeanMedianStdev,
                          csv = dataggplotMixtureSpatialMeanMedianStdev,
                          txt = dataggplotMixtureSpatialMeanMedianStdev)
    )
    
    downloadablePlot(id = "ggplotMixtureSpatialMeanMedianStdev",
                     logger = ss_userAction.Log,
                     filenameroot = "ggplotMixtureSpatialMeanMedianStdev",
                     aspectratio  = 1,
                     downloadfxns = list(png=ggplotMixtureSpatialMeanMedianStdev,
                                         csv = dataggplotMixtureSpatialMeanMedianStdev,
                                         txt = dataggplotMixtureSpatialMeanMedianStdev),
                     visibleplot  = ggplotMixtureSpatialMeanMedianStdev)
    
    dataggplotMixturePredictorMeanMedianStdevFit <- function(){
      DF <- MixtureModelFit()$DFpredictorMeanMedianStdev$DFpredictorMeanMedianStdev
      return(DF)
    }
    
    DFMixturePredictor <- reactive({
      DF <- MixtureModelFit()$DFpredictorMeanMedianStdev$DFpredictorMeanMedianStdev
      condition <- !(length(unique(DF[,1]))*length(unique(DF[,2]))==nrow(DF))&val()
      if(condition){
        grid <- expand.grid(seq(min(DF[,1]),max(DF[,1]), length.out=200),seq(min(DF[,2]),max(DF[,2]), length.out=200))
        DFInter <- data.frame(Latitude=grid[,1],Longitude=grid[,2])
        DFInter$Predictor.mean <- InterpolateIrrGrid(z=DF$Predictor.mean,loc=DF[,1:2], gridInter=grid)$DataInter$z
        DFInter$Predictor.median <- InterpolateIrrGrid(z=DF$Predictor.median,loc=DF[,1:2], gridInter=grid)$DataInter$z
        DFInter$Predictor.stdev <- InterpolateIrrGrid(z=DF$Predictor.stdev,loc=DF[,1:2], gridInter=grid)$DataInter$z
        colnames(DFInter)[1:2] <- colnames(DF)[1:2]
        DF <- DFInter
      }
      
      return(DF)
    })
    
    ggplotMixturePredictorMeanMedianStdevFit <- function(){
      DF <- DFMixturePredictor()
      g1 <- ggplot(DF) + geom_tile(aes(x=DF[,1], y=DF[,2], fill=Predictor.mean)) +
        scale_fill_viridis_c(option = "turbo") + theme_bw() + coord_fixed() + labs(fill="Values") +
        xlab(colnames(DF)[1]) + ylab(colnames(DF)[1]) + ggtitle("Linear predictor (mean)") +
        theme(plot.title=element_text(color = "black", size = 16, face = "bold", hjust = 0.5))
      
      g2 <- ggplot(DF) + geom_tile(aes(x=DF[,1], y=DF[,2], fill=Predictor.median)) +
        scale_fill_viridis_c(option = "turbo") + theme_bw() + coord_fixed() + labs(fill="Values") +
        xlab(colnames(DF)[1]) + ylab(colnames(DF)[1]) + ggtitle("Linear predictor (median)") +
        theme(plot.title=element_text(color = "black", size = 16, face = "bold", hjust = 0.5))
      
      g3 <- ggplot(DF) + geom_tile(aes(x=DF[,1], y=DF[,2], fill=Predictor.stdev))+
        scale_fill_viridis_c(option = "turbo") + theme_bw() + coord_fixed() + labs(fill="Values") +
        xlab(colnames(DF)[1]) + ylab(colnames(DF)[1]) + ggtitle("Linear predictor (stdev.)") +
        theme(plot.title=element_text(color = "black", size = 16, face = "bold", hjust = 0.5))
      
      gt <- g1+g2+g3
      
      return(gt)
    }
    
    downloadFile(
      id = "ggplotMixturePredictorMeanMedianStdevFit",
      logger = ss_userAction.Log,
      filenameroot = "ggplotMixturePredictorMeanMedianStdevFit",
      aspectratio  = 1,
      downloadfxns = list(png  = ggplotMixturePredictorMeanMedianStdevFit,
                          csv = dataggplotMixturePredictorMeanMedianStdevFit,
                          txt = dataggplotMixturePredictorMeanMedianStdevFit)
    )
    
    downloadablePlot(id = "ggplotMixturePredictorMeanMedianStdevFit",
                     logger = ss_userAction.Log,
                     filenameroot = "ggplotMixturePredictorMeanMedianStdevFit",
                     aspectratio  = 1,
                     downloadfxns = list(png = ggplotMixturePredictorMeanMedianStdevFit,
                                         csv = dataggplotMixturePredictorMeanMedianStdevFit,
                                         txt = dataggplotMixturePredictorMeanMedianStdevFit),
                     visibleplot  = ggplotMixturePredictorMeanMedianStdevFit)
    
    dataggplotMixtureFixParamFit <- function(){
      DF <- as.data.frame(MixtureModelFit()$DFPostFixed$DFPostFixed)
      return(DF)
    }
    
    ggplotMixtureFixParamFit <- function(){
      DF <- MixtureModelFit()$DFPostFixed$DFPostFixed
      gl <- c()
      for(i in 1:length(DF)){
        assign(paste0("g",i),
               ggplot(data=data.frame(x=DF[[i]][,1], y=DF[[i]][,2]), aes(x=x,y=y)) + geom_line() +
                 theme_bw() + xlab(names(DF)[i]) + ylab(HTML(paste("Density f(",names(DF)[i],")"))) + 
                 ggtitle(names(DF)[i]) + theme(plot.title=element_text(color = "black", size = 16, face = "bold", hjust = 0.5))
        )
        gl[i] <- paste0("g",i)
      }
      gt <- eval(parse(text=paste(gl, collapse="+")))
      return(gt)
    }
    
    downloadFile(
      id = "ggplotMixtureFixParamFit",
      logger = ss_userAction.Log,
      filenameroot = "ggplotMixtureFixParamFit",
      aspectratio  = 1,
      downloadfxns = list(png  = ggplotMixtureFixParamFit,
                          csv = dataggplotMixtureFixParamFit,
                          txt = dataggplotMixtureFixParamFit)
    )
    
    downloadablePlot(id = "ggplotMixtureFixParamFit",
                     logger = ss_userAction.Log,
                     filenameroot = "ggplotMixtureFixParamFit",
                     aspectratio  = 1,
                     downloadfxns = list(png = ggplotMixtureFixParamFit,
                                         csv = dataggplotMixtureFixParamFit,
                                         txt = dataggplotMixtureFixParamFit),
                     visibleplot  = ggplotMixtureFixParamFit)
    
    
    dataggplotMixtureHyperParamFit <- function(){
      DF <- as.data.frame(MixtureModelFit()$DFPostHyperpar$DFPostHyperpar)
      return(DF)
    }
    
    ggplotMixtureHyperParamFit <- function(){
      DF <- MixtureModelFit()$DFPostHyperpar$DFPostHyperpar
      gl <- c()
      for(i in 1:length(DF)){
        nm <- strsplit(names(DF)[i], " ")[[1]]
        if(nm[1]=="Precision"){
          namesDFold <- names(DF)[i]
          names(DF)[i] <- paste("Stdev.", paste(nm[-1], collapse=" "), collapse=" ")
          assign(paste0("g",i), try(
                 ggplot(data=data.frame(x=inla.tmarginal(function(x) 1/x**0.5, DF[[i]])[,1], 
                                        y=inla.tmarginal(function(x) 1/x**0.5, DF[[i]])[,2]), aes(x=x,y=y)) + geom_line() +
                   theme_bw()+ xlab(names(DF)[i]) +
                   ylab(HTML(paste("Density f(",names(DF)[i],")"))) + ggtitle(names(DF)[i]) +
                   theme(plot.title=element_text(color = "black", size = 16, face = "bold", hjust = 0.5)), silent=TRUE)
          )
          if(sum(class(eval(parse(text=paste0("g",i))))=="try-error")==1){
            assign(paste0("g",i),
              ggplot(data=data.frame(x=inla.smarginal(marginal=DF[[i]])[[1]], 
                                     y=inla.smarginal(marginal=DF[[i]])[[2]]), aes(x=x,y=y)) + geom_line() +
                theme_bw()+ xlab(names(DF)[i]) +
                ylab(HTML(paste("Density f(",namesDFold,")"))) + ggtitle(namesDFold) +
                theme(plot.title=element_text(color = "black", size = 16, face = "bold", hjust = 0.5))
            )
          }
        } else{
          assign(paste0("g",i),
                 ggplot(data=data.frame(x=DF[[i]][,1], y=DF[[i]][,2]), aes(x=x,y=y)) + geom_line() +
                   theme_bw()+ xlab(names(DF)[i]) + xlim(quantile(DF[[i]][,1], probs = c(0.025,0.975))) +
                   ylab(HTML(paste("Density f(",names(DF)[i],")"))) + ggtitle(names(DF)[i]) +
                   theme(plot.title=element_text(color = "black", size = 16, face = "bold", hjust = 0.5))
          )
        }
        gl[i] <- paste0("g",i)
      }
      gt <- eval(parse(text=paste(gl, collapse="+")))
      return(gt)
    }
    
    downloadFile(
      id = "ggplotMixtureHyperParamFit",
      logger = ss_userAction.Log,
      filenameroot = "ggplotMixtureHyperParamFit",
      aspectratio  = 1,
      downloadfxns = list(png  = ggplotMixtureHyperParamFit,
                          csv = dataggplotMixtureHyperParamFit,
                          txt = dataggplotMixtureHyperParamFit)
    )
    
    downloadablePlot(id = "ggplotMixtureHyperParamFit",
                     logger = ss_userAction.Log,
                     filenameroot = "ggplotMixtureHyperParamFit",
                     aspectratio  = 1,
                     downloadfxns = list(png=ggplotMixtureHyperParamFit,
                                         csv = dataggplotMixtureHyperParamFit,
                                         txt = dataggplotMixtureHyperParamFit),
                     visibleplot  = ggplotMixtureHyperParamFit)
    
    
    tableMixtureModelFixedPar <- function(){
      DF <- MixtureModelFit()$Summary.fixed$Summary.fixed %>%
        mutate(across(where(is.numeric), round, digits = 2))
    }
    
    downloadableTable("tableMixtureModelFixedPar",
                      logger=ss_userAction.Log,
                      filenameroot="tableMixtureModelFixedPar",
                      downloaddatafxns=list(csv=tableMixtureModelFixedPar,
                                            tsv=tableMixtureModelFixedPar),
                      tabledata=tableMixtureModelFixedPar, rownames = TRUE,
                      caption="Summary fixed parameters")
    
    tableMixtureModelHyperPar <- function(){
      DF <- MixtureModelFit()$Summary.hyperpar$Summary.hyperpar %>%
        mutate(across(where(is.numeric), round, digits = 2))
    }
    
    downloadableTable("tableMixtureModelHyperPar",
                      logger=ss_userAction.Log,
                      filenameroot="tableMixtureModelHyperPar",
                      downloaddatafxns=list(csv=tableMixtureModelHyperPar,
                                            tsv=tableMixtureModelHyperPar),
                      tabledata=tableMixtureModelHyperPar, rownames = TRUE,
                      caption="Summary hyperparameters")
    
    tableMixtureModelInternalHyperPar <- function(){
      DF <- MixtureModelFit()$SummaryInternalHyper$SummaryInternalHyper %>%
        mutate(across(where(is.numeric), round, digits = 2))
    }
    
    downloadableTable("tableMixtureModelInternalHyperPar",
                      logger=ss_userAction.Log,
                      filenameroot="tableMixtureModelInternalHyperPar",
                      downloaddatafxns=list(csv=tableMixtureModelInternalHyperPar,
                                            tsv=tableMixtureModelInternalHyperPar),
                      tabledata=tableMixtureModelInternalHyperPar, rownames = TRUE,
                      caption="Summary internal hyperparameters")
    
    dataMixtureDICtable <- function(){
      DF <- MixtureModelFit()$DICmodel$DICmodel %>%
        mutate(across(where(is.numeric), round, digits = 2))
      return(DF)
    }
    
    downloadableTable("dataMixtureDICtable",
                      logger=ss_userAction.Log,
                      filenameroot="dataMixtureDICtable",
                      downloaddatafxns=list(csv=dataMixtureDICtable,
                                            tsv=dataMixtureDICtable),
                      tabledata=dataMixtureDICtable, rownames = TRUE,
                      caption="Model DIC")
    
    dataMixtureWAICtable <- function(){
      DF <- MixtureModelFit()$WAICmodel$WAICmodel %>%
        mutate(across(where(is.numeric), round, digits = 2))
      return(DF)
    }
    
    downloadableTable("dataMixtureWAICtable",
                      logger=ss_userAction.Log,
                      filenameroot="dataMixtureWAICtable",
                      downloaddatafxns=list(csv=dataMixtureWAICtable,
                                            tsv=dataMixtureWAICtable),
                      tabledata=dataMixtureWAICtable, rownames = TRUE,
                      caption="Model WAIC")
    
    dataMixtureCPOtable <- function(){
      CPO <- MixtureModelFit()$SummaryCPO$SummaryCPO
      DF <- data.frame(n=length(CPO), mean=mean(CPO), median=median(CPO), stdev=sd(CPO), 
                       quantile2.5=quantile(CPO,probs=c(0.025,0.975))[1], 
                       quantile97.5=quantile(CPO,probs=c(0.025,0.975))[2]) %>%
        mutate(across(where(is.numeric), round, digits = 2))
      return(DF)
    }
    
    downloadableTable("dataMixtureCPOtable",
                      logger=ss_userAction.Log,
                      filenameroot="dataMixtureCPOtable",
                      downloaddatafxns=list(csv=dataMixtureCPOtable,
                                            tsv=dataMixtureCPOtable),
                      tabledata=dataMixtureCPOtable, rownames = FALSE,
                      caption="Summary CPO")
    
    # ----
    
    periscope:::fw_server_setup(input, output, session, ss_userAction.Log)
    
    loginfo("%s started with log level <%s>",
            periscope:::fw_get_title(), periscope:::fw_get_loglevel(),
            logger = ss_userAction.Log)
})
