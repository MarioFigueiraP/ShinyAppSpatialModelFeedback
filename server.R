ss_userAction.Log <- periscope:::fw_get_user_log()

# Define server logic required to draw a histogram
shinyServer(function(input, output, session) {
    
    ### Update buttons 
    # observeEvent(input$parametersMapSim, ({
    #     updateButton(session, "parametersMapSim", style = ifelse(input$parametersMapSim, "primary", "default"))
    # }))
    # observeEvent(input$parametersSampleSim, ({
    #     updateButton(session, "parametersSampleSim", style = ifelse(input$parametersSampleSim, "primary", "default"))
    # }))
    # observeEvent(input$infoSim, ({
    #     updateButton(session, "infoSim", style = ifelse(input$infoSim, "warning", "default"))
    # }))
    
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
    
    
    # observeEvent(input$makeSim, {
    #     updateBox("resultbox",
    #               action = "toggle")
    #               #options = list(collapsible=ifelse(input$makeSim>0,FALSE,TRUE)))
    # })
    
    # observe({
    #     if (input$close > 0) stopApp() # close shiny
    # })
    
    ### Some functions
    
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
        lattice <- inla.mesh.lattice(seq(lim1[1],lim1[2],length.out=2),seq(lim2[1],lim2[2],length.out=2))
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
    
    ######################################################################
    ##### ################ ###### Simulation ###### ################ #####  
    ######################################################################
    
    ### Simulation of the mesh, covariates and abundance

    meshSim <- eventReactive(input$makeSim ,{
        xlattice0 <- seq(input$limlattice[1], input$limlattice[2], length.out=input$lengthlattice)
        ylattice0 <- seq(input$limlattice[1], input$limlattice[2], length.out=input$lengthlattice)
        lattice0 <-  inla.mesh.lattice(xlattice0, ylattice0)
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
        if(!is.na(input$seedGlobal)){set.seed(input$seedGlobal)}
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
                    #print("Está llevando mucho tiempo")
                    k <- k+1
                    warningMessage <- paste0("Simulate the random walk (rw1) is taking a while (more than ", 
                                   round(as.numeric(difftime(tinwhile,tprewhile, units="secs")), digits=0),
                           " secs). Maybe would be better restart the process with a wider interval between 
                           extreme values or with a higher precision.")
                    showNotification(ui=warningMessage, duration=10, closeButton=TRUE, type="warning")
                    #break
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
                bat.eff.seq <- RW1(N=input$nknots.rw2, x0=input$init.rw2, mu=0, stdev=1/(input$prec.rw2)**0.5)
                min.eff <- min(bat.eff.seq);max.eff <- max(bat.eff.seq)
                tinwhile <- Sys.time()
                if (as.numeric(difftime(tinwhile,tprewhile, units="secs"))>20*k){
                    #print("Está llevando mucho tiempo")
                    k <- k+1
                    warningMessage <- paste0("Simulate the random walk (rw1) is taking a while (more than ", 
                                             round(as.numeric(difftime(tinwhile,tprewhile, units="secs")), digits=0),
                                             " secs). Maybe would be better restart the process with a wider interval between 
                           extreme values or with a higher precision.")
                    showNotification(ui=warningMessage, duration=10, closeButton=TRUE, type="warning")
                    #break
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
    
    # Data and plot of the spatial effect
    
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
    
    #Data and relation between the bathymetry and its effect
    
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
    
    # Bathymetric effect chart
    
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
    
    # Abundance map
    
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
    
    Random.sampling <- eventReactive(input$makeSample, {
        if(!is.na(input$seedSampleR)){set.seed(input$seedSampleR)}
        indx <- sample(1:length(SimMap()$DataSim$DataSim$gridx), size=input$niid.samples)
        RandomSample <- list(Latitude=SimMap()$DataSim$DataSim$gridx[indx], 
                                   Longitude=SimMap()$DataSim$DataSim$gridy[indx], 
                                   Bathymetry=SimMap()$DataSim$DataSim$bathymetry.field[indx], 
                                   Abundance=SimMap()$DataSim$DataSim$abundance.field[indx])
        return(RandomSample)
    })
    
    dataggplotRandomSampleSim <- function(){
        data.frame(Latitude=Random.sampling()$Latitude, Longitude=Random.sampling()$Longitude,
                   Bathymetry=Random.sampling()$Bathymetry, Abundance=Random.sampling()$Abundance)
    }
    
    ggplotRandomSampleSim <- function(){
        DFM <- data.frame(Latitude=SimMap()$DataSim$DataSim$gridx, Longitude=SimMap()$DataSim$DataSim$gridy,
                          Abundance.field=SimMap()$DataSim$DataSim$abundance.field)
        DFS <- data.frame(Latitude=Random.sampling()$Latitude, Longitude=Random.sampling()$Longitude)
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
        if(!is.na(input$seedSampleP)){set.seed(input$seedSampleP)}
        indx <- sample(1:length(SimMap()$DataSim$DataSim$gridx), size=input$nps.samples,
                       prob=exp(input$r.scale*SimMap()$DataSim$DataSim$abundance.field/max(SimMap()$DataSim$DataSim$abundance.field)))
        PreferentialSample <- list(Latitude=SimMap()$DataSim$DataSim$gridx[indx],
                             Longitude=SimMap()$DataSim$DataSim$gridy[indx],
                             Bathymetry=SimMap()$DataSim$DataSim$bathymetry.field[indx],
                             Abundance=SimMap()$DataSim$DataSim$abundance.field[indx])
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

    
    #######################################################################
    ##### ################ ###### Upload Data ###### ################ #####
    #######################################################################

    # output$content <- renderTable({input$file.uploadData})

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
    
    quiltplotSampleReadRasterBath <- function(){
        DF <- datareadSample()
        DFggplot <- ggplot(DF) + geom_point(aes(x=DF[,1],y=DF[,2],colour=DF[,3])) +
            scale_colour_viridis_c(option = "turbo") + labs(colour=paste(names(DF)[3],"\nValues")) + theme_bw() +
            xlab(names(DF)[1]) + ylab(names(DF)[2]) + ggtitle(names(DF)[3]) +
            theme(plot.title=element_text(color = "black", size = 16, face = "bold", hjust = 0.5))
        return(DFggplot)
    }
    
    quiltplotSampleReadRasterAbundance <- function(){
        DF <- datareadSample()
        DFggplot <- ggplot(DF) + geom_point(aes(x=DF[,1],y=DF[,2],colour=DF[,4])) +
            scale_colour_viridis_c(option = "turbo") + labs(colour=paste(names(DF)[4],"\nValues")) + theme_bw() +
            xlab(names(DF)[1]) + ylab(names(DF)[2]) + ggtitle(names(DF)[4]) +
            theme(plot.title=element_text(color = "black", size = 16, face = "bold", hjust = 0.5))
        return(DFggplot)
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
        id = "quilt.plot.sample.bathymetry",
        logger = ss_userAction.Log,
        filenameroot = "quilt.plot.sample.bathymetry",
        aspectratio  = 1,
        downloadfxns = list(png  = quiltplotSampleReadRasterBath)
    )
    
    downloadablePlot(id = "quilt.plot.sample.bathymetry",
                     logger = ss_userAction.Log,
                     filenameroot = "quilt.plot.sample.bathymetry",
                     aspectratio  = 1,
                     downloadfxns = list(png=quiltplotSampleReadRasterBath),
                     visibleplot  = quiltplotSampleReadRasterBath)
    
    downloadFile(
        id = "quilt.plot.sample.abundance",
        logger = ss_userAction.Log,
        filenameroot = "quilt.plot.sample.abundance",
        aspectratio  = 1,
        downloadfxns = list(png  = quiltplotSampleReadRasterAbundance)
    )
    
    downloadablePlot(id = "quilt.plot.sample.abundance",
                     logger = ss_userAction.Log,
                     filenameroot = "quilt.plot.sample.abundance",
                     aspectratio  = 1,
                     downloadfxns = list(png=quiltplotSampleReadRasterAbundance),
                     visibleplot  = quiltplotSampleReadRasterAbundance,
                     caption="Observational Data")
    
    file.BathymetryRaster.read <- reactive({
        if(is.null(input$file.uploadDataBathymetryRaster)) return(NULL)
        input$file.uploadDataBathymetryRaster
    })
    
    output$fileUploadedRaster <- reactive({
        return(!is.null(file.BathymetryRaster.read()))
    })
    outputOptions(output, 'fileUploadedRaster', suspendWhenHidden=FALSE)
    
    datareadBathymetricRaster <- function(){
        file <- file.BathymetryRaster.read()
        ext <- tools::file_ext(file$datapath)
        req(file)
        validate(need(ext == c("csv","rds"), "Please upload a csv or rds file"))
        if(ext=="csv"){file.read <- read.csv(file$datapath)}
        else file.read <- readRDS(file$datapath)
        return(file.read)
        # return(as.data.frame(file.read))
    }
    
    downloadableTable("table.read.bathymetry.raster",
                      logger=ss_userAction.Log,
                      filenameroot="table.read.bathymetry.raster",
                      downloaddatafxns=list(csv=datareadBathymetricRaster,
                                            tsv=datareadBathymetricRaster),
                      tabledata=datareadBathymetricRaster,
                      caption="Bathymetry Data Raster"
    )
    
    val <- reactiveVal(FALSE)
    
    observeEvent(input$EnhancementIrrGrid, {
        val(!val())
    })
    
    
    DFbathyRaster <- reactive({
        DF <- datareadBathymetricRaster()
        condition <- !(length(unique(DF[,1]))*length(unique(DF[,2]))==nrow(DF))&val()
        if(condition){
            if(input$RefineInterpolationIrrGrid){
                names <- names(DF)
                DF <- RefineInterMesh(z=DF[,3],loc=DF[,1:2])$DataInter; names(DF) <- names[1:3]
            } else{
                names <- names(DF)
                grid <- expand.grid(seq(min(DF[,1]),max(DF[,1]), length.out=200),seq(min(DF[,2]),max(DF[,2]), length.out=200))
                DF <- InterpolateIrrGrid(z=DF[,3],loc=DF[,1:2], gridInter=grid)$DataInter; names(DF) <- names[1:3]
            }
            
        }

        return(DF)
    })
    
    quiltplotBathymetryRaster <- function(){
        # DF <- datareadBathymetricRaster()
        # grid <- expand.grid(seq(min(DF[,1]),max(DF[,1]), length.out=200),seq(min(DF[,2]),max(DF[,2]), length.out=200))
        # if(!(length(unique(DF[,1]))*length(unique(DF[,2]))==nrow(DF))&val()){
        #     DF <- InterpolateIrrGrid(z=DF[,3],loc=DF[,1:2], gridInter=grid)$DataInter
        # }
        DF_base <- datareadBathymetricRaster()
        DF <- DFbathyRaster()
        
        DFggplot <- ggplot(DF) + geom_tile(aes(x=DF[,1],y=DF[,2], fill=DF[,3])) + 
            geom_point(data=DF_base, aes(x=DF_base[,1],y=DF_base[,2]), size=1) +
            scale_fill_viridis_c(option = "turbo") + labs(fill=paste(names(DF)[3],"\nValues")) + theme_bw() +
            xlab(names(DF)[1]) + ylab(names(DF)[2]) + ggtitle(paste(names(DF)[3],"Raster")) +
            theme(plot.title=element_text(color = "black", size = 16, face = "bold", hjust = 0.5))
        return(DFggplot)
    }
    
    downloadFile(
        id = "quilt.plot.bathymetry",
        logger = ss_userAction.Log,
        filenameroot = "quilt.plot.bathymetry",
        aspectratio  = 1,
        downloadfxns = list(png  = quiltplotBathymetryRaster)
    )
    
    downloadablePlot(id = "quilt.plot.bathymetry",
                     logger = ss_userAction.Log,
                     filenameroot = "quilt.plot.bathymetry",
                     aspectratio  = 1,
                     downloadfxns = list(png=quiltplotBathymetryRaster),
                     visibleplot  = quiltplotBathymetryRaster)
    
    #############################################################################################
    ##### ################ ###### Independent and Random Model Data ###### ################ #####
    #############################################################################################
    
    RandomCMesh <- reactive({
      #taking the data from simulation or from the loading tab
      if(input$DataSimulatedLoaded=="sim"){
        DFsample <- as.data.frame(Random.sampling())
      } else if(input$DataSimulatedLoaded=="load"){
        DFsample <- datareadSample()
      }
      
      if(input$BathymetryRasterSPDE=="raster"){
        rasterSample <- datareadBathymetricRaster()[sample(1:nrow(datareadBathymetricRaster()), min(c(60,nrow(datareadBathymetricRaster())))),1:2]
        qloc <- quantile(as.vector(dist(rasterSample)),probs=c(ifelse(input$RandomCustomMesh,input$RandomMeshQloc,0.03),0.3))
          # quantile(as.vector(dist(rasterSample)),probs=c(0.01,0.3))
        mesh <- inla.mesh.2d(loc=cbind(rasterSample[,1],rasterSample[,2]), cutoff = qloc[1]/2, offset=c(-0.1, -0.4), 
                             max.edge=c(qloc[1], qloc[2]))
        sample <- rasterSample
      } else if(input$BathymetryRasterSPDE=="solvebathy"){
        qloc <- quantile(as.vector(dist(DFsample[sample(1:nrow(DFsample),size=min(c(50,nrow(DFsample)))),1:2])),probs=c(ifelse(input$RandomCustomMesh,input$RandomMeshQloc,0.03),0.3))
          # quantile(as.vector(dist(DFsample[sample(1:nrow(DFsample),size=60),1:2])),probs=c(0.01,0.3))
        mesh <- inla.mesh.2d(loc=cbind(DFsample[,1],DFsample[,2]), cutoff = qloc[1]/2, offset=c(-0.1, -0.4), 
                             max.edge=c(qloc[1], qloc[2]))
        sample <- DFsample
      }
      result <- list(mesh=mesh, qloc=qloc, Sample=sample)
      return(result)
    })
    
    
    ggplotRandomMesh <- function(){
      ggplot()+ gg(RandomCMesh()$mesh)+ theme_bw() + xlab("Latitude") + ylab("Longitude") +
        ggtitle("Mesh over the study region") + 
        geom_point(data=RandomCMesh()$Sample, 
                   aes(x=RandomCMesh()$Sample[,1],y=RandomCMesh()$Sample[,2]), size=1) + 
        theme(plot.title=element_text(color = "black", size = 16, face = "bold", hjust = 0.5))
    }
    
    downloadFile(
      id = "ggplotRandomMesh",
      logger = ss_userAction.Log,
      filenameroot = "ggplotRandomMesh",
      aspectratio  = 1,
      downloadfxns = list(png  = ggplotRandomMesh)
    )
    
    downloadablePlot(
      id = "ggplotRandomMesh",
      logger = ss_userAction.Log,
      filenameroot = "ggplotRandomMesh",
      aspectratio  = 1,
      downloadfxns = list(png  = ggplotRandomMesh),
      visibleplot  = ggplotRandomMesh)
    
    RandomModelFit <- eventReactive(input$fitIndRandom, {
        showNotification(ui=paste("Fitting the data."), duration = NULL)
        t1 <- Sys.time()
        
        #taking the data from simulation or from the loading tab
        if(input$DataSimulatedLoaded=="sim"){
          DFsample <- as.data.frame(Random.sampling())
        } else if(input$DataSimulatedLoaded=="load"){
          DFsample <- datareadSample()
        }
        
        mesh <- RandomCMesh()$mesh
        qloc <- RandomCMesh()$qloc


        #now we'll make the spde structure
        if(input$optionRPS=="auto"){
            prior.range <- c(qloc[2],0.5)
            prior.sigma <- c(1,0.5)
        } else if(input$optionRPS=="custom"){
            prior.range <- as.numeric(unlist(strsplit(input$RandomPriorRange, split=",")))
            prior.sigma <- as.numeric(unlist(strsplit(input$RandomPriorStdev, split=",")))
        }

        spde <- inla.spde2.pcmatern(mesh, prior.range = prior.range, prior.sigma = prior.sigma, alpha=2)
        A.sample <- inla.spde.make.A(mesh = mesh, loc = cbind(DFsample[,1], DFsample[,2]))
        s.index <- inla.spde.make.index(name = "spatial", n.spde = spde$n.spde)
        
        
        
        # Now the bathymetric stage to the prediction
        if(input$BathymetryRasterSPDE=="raster"&input$rasterBatymetryPred=='SPDEraster'){
            DFraster <- datareadBathymetricRaster()

            gridPred <- expand.grid(x=seq(min(DFraster[,1]), max(DFraster[,2]), length.out=input$SPDErasterInput),
                                    y=seq(min(DFraster[,1]), max(DFraster[,2]), length.out=input$SPDErasterInput))

            A.raster.inf <- inla.spde.make.A(mesh=mesh, loc=as.matrix(DFraster[,1:2]))
            Bat.raster.inf <- inla.stack(data=list(y=DFraster[,3]),
                                  A = list(A.raster.inf,1),
                                  effects = list(list(spatial=s.index$spatial),
                                                 list(b0 = rep(1, nrow(DFraster)))),
                                  tag="raster.inf")

            A.raster.pred <- inla.spde.make.A(mesh = mesh, loc = as.matrix(gridPred))
            Bat.raster.pred <- inla.stack(data = list(y = NA),
                                   A = list(A.raster.pred, 1),
                                   effects = list(list(spatial=s.index$spatial),
                                                  list (b0 = rep(1, dim(gridPred)[1]))),
                                   tag = "raster.pred")
            Bat.total.stack <- inla.stack(Bat.raster.inf, Bat.raster.pred)

            batrandom.model <- inla(y ~ -1 + b0 + f(spatial, model = spde),
                                    data = inla.stack.data(Bat.total.stack),
                                    family = "gaussian",
                                    control.inla = list(strategy="simplified.laplace", int.strategy="eb"),
                                    control.predictor = list(A = inla.stack.A(Bat.total.stack), compute = TRUE, link = 1),
                                    control.compute = list(cpo = FALSE, dic = FALSE),
                                    verbose=FALSE)
            index.pred <- inla.stack.index(Bat.total.stack, "raster.pred")$data
            bathymetry.pred <- batrandom.model$summary.fitted.values[index.pred, "mean"]
            bathymetry.pred.sd <- batrandom.model$summary.fitted.values[index.pred, "sd"]
        } else if(input$BathymetryRasterSPDE=="raster"&input$rasterBatymetryPred=='smoothraster'){
            DFraster <- datareadBathymetricRaster()
            gridPred <- expand.grid(x=seq(min(DFraster[,1]), max(DFraster[,2]), length.out=input$SPDErasterInput),
                                    y=seq(min(DFraster[,1]), max(DFraster[,2]), length.out=input$SPDErasterInput))
            bathymetry.pred <- ProjectionFunctionRegularGrid(lim1=range(DFraster[,1]), lim2=range(DFraster[,2]),
                                                             loc=as.matrix(DFraster[,1:2]), z=DFraster[,3], proj.grid=gridPred)
        } else if(input$BathymetryRasterSPDE=="raster"&input$rasterBatymetryPred=='rasterpred'){
            DFraster <- datareadBathymetricRaster()
            gridPred <- DFraster[,1:2]
            bathymetry.pred <- DFraster[,3]
        } else if(input$BathymetryRasterSPDE=="solvebathy"){
            gridPred <- expand.grid(x=seq(min(DFsample[,1]), max(DFsample[,2]), length.out=input$dimrandompred),
                                    y=seq(min(DFsample[,1]), max(DFsample[,2]), length.out=input$dimrandompred))

            A.raster.inf <- inla.spde.make.A(mesh=mesh, loc=as.matrix(DFsample[,1:2]))
            Bat.raster.inf <- inla.stack(data=list(y=DFsample[,3]),
                                         A = list(A.sample ,1),
                                         effects = list(list(spatial=s.index$spatial),
                                                        list(b0 = rep(1, nrow(DFsample)))),
                                         tag="raster.inf")

            A.raster.pred <- inla.spde.make.A(mesh = mesh, loc = as.matrix(gridPred))
            Bat.raster.pred <- inla.stack(data = list(y = NA),
                                          A = list(A.raster.pred, 1),
                                          effects = list(list(spatial=s.index$spatial),
                                                         list (b0 = rep(1, dim(gridPred)[1]))),
                                          tag = "raster.pred")
            Bat.total.stack <- inla.stack(Bat.raster.inf, Bat.raster.pred)

            batrandom.model <- inla(y ~ -1 + b0 + f(spatial, model = spde),
                                    data = inla.stack.data(Bat.total.stack),
                                    #family = "gaussian",
                                    control.inla = list(strategy="auto", int.strategy="auto"),
                                    control.predictor = list(A = inla.stack.A(Bat.total.stack), compute = TRUE, link = 1),
                                    control.compute = list(cpo = FALSE, dic = FALSE),
                                    verbose=FALSE)
            index.pred <- inla.stack.index(Bat.total.stack, "raster.pred")$data
            bathymetry.pred <- batrandom.model$summary.fitted.values[index.pred, "mean"]
            bathymetry.pred.sd <- batrandom.model$summary.fitted.values[index.pred, "sd"]
        }

        
        if(input$optionRPB=='custom'){
            CovariatesPriorParam <- list(mean=list(Intercept=as.numeric(unlist(strsplit(input$RandomMeanPriorBeta,",")))[1]),
                                         prec=list(Intercept=as.numeric(unlist(strsplit(input$RandomMeanPriorBeta,",")))[2]))
        } else if(input$optionRPB=='default'&input$RPModelBathy!="lin"){
            CovariatesPriorParam <- inla.set.control.fixed.default()
        } else {
            CovariatesPriorParam <- list(mean=list(Intercept=NA),
                                         prec=list(Intercept=NA))
        }
        
        if(input$RPModelBathy=="bs"){
            
            knots <- seq(min(bathymetry.pred), max(bathymetry.pred), length.out=input$bsRandomnKnots)
            mod.bathymetry <- "bs(Bathymetry, knots=knots)"
            
            effects.list <- list(list(spatial=s.index$spatial),
                                 list(Intercept = rep(1, nrow(DFsample)[1]),
                                      Bathymetry = DFsample[,3]))
            A.inf <- list(A.sample ,1)
            
            Pred.effects.list <- list(list(spatial=s.index$spatial),
                                      list(Intercept = rep(1, nrow(gridPred)[1]),
                                           Bathymetry = bathymetry.pred))
            DFpred <- data.frame(Latitude=gridPred[,1], Longitude=gridPred[,2],
                                 Bathymetry=bathymetry.pred)
            
            A.pred <- inla.spde.make.A(mesh = mesh, loc = as.matrix(gridPred))
            AA.pred <- list(A.pred,1)
            
            rn.for <- paste("y ~ -1", "Intercept", mod.bathymetry, "f(spatial, model=spde)", sep = " + ")
            
        } else if(input$RPModelBathy=="rw1"){
            
            group <- inla.group(c(DFsample[,3],bathymetry.pred), n = input$rw1RandomnKnots, method = "quantile")
            group.inf <- group[1:length(DFsample[,3])]
            group.pred <- group[-(1:length(DFsample[,3]))]
            
            if(input$autocustomRandomRw1=='custom'){
                loggammapar <- as.numeric(unlist(strsplit(input$RandomPrecRw1, split=",")))
                mod.bathymetry <- paste0("f(Bathymetry, model='rw1', hyper = list(prec = list(prior='loggamma',param=c(",
                                         loggammapar[1],",",loggammapar[2],"))))")
            } else{mod.bathymetry <- "f(Bathymetry, model='rw1')"}
                
            
            
            effects.list <- list(list(spatial=s.index$spatial),
                                 list(Intercept = rep(1, nrow(DFsample)[1]),
                                      Bathymetry = group.inf))
            A.inf <- list(A.sample ,1)
            
            Pred.effects.list <- list(list(spatial=s.index$spatial),
                                      list(Intercept = rep(1, nrow(gridPred)[1]),
                                           Bathymetry = group.pred))
            DFpred <- data.frame(Latitude=gridPred[,1], Longitude=gridPred[,2],
                                 Bathymetry=bathymetry.pred)
            
            A.pred <- inla.spde.make.A(mesh = mesh, loc = as.matrix(gridPred))
            AA.pred <- list(A.pred,1)
            rn.for <- paste("y ~ -1", "Intercept", mod.bathymetry, "f(spatial, model=spde)", sep = " + ")
            
        } else if(input$RPModelBathy=="rw2"){
            
            group <- inla.group(c(DFsample[,3],bathymetry.pred), n = input$rw2RandomnKnots, method = "quantile")
            group.inf <- group[1:length(DFsample[,3])]
            group.pred <- group[-(1:length(DFsample[,3]))]
            
            if(input$autocustomRandomRw2=='custom'){
                loggammapar <- as.numeric(unlist(strsplit(input$RandomPrecRw2, split=",")))
                mod.bathymetry <- paste0("f(Bathymetry, model='rw2', hyper = list(prec = list(prior='loggamma',param=c(",
                                         loggammapar[1],",",loggammapar[2],"))))")
            } else {mod.bathymetry <- "f(Bathymetry, model='rw2')"}
            A.inf <- list(A.sample ,1)
            
            effects.list <- list(list(spatial=s.index$spatial),
                                 list(Intercept = rep(1, nrow(DFsample)[1]),
                                      Bathymetry = group.inf))
            
            Pred.effects.list <- list(list(spatial=s.index$spatial),
                                      list(Intercept = rep(1, nrow(gridPred)[1]),
                                           Bathymetry = group.pred))
            
            DFpred <- data.frame(Latitude=gridPred[,1], Longitude=gridPred[,2],
                                 Bathymetry=group.pred)
            
            A.pred <- inla.spde.make.A(mesh = mesh, loc = as.matrix(gridPred))
            AA.pred <- list(A.pred,1)
            rn.for <- paste("y ~ -1", "Intercept", mod.bathymetry, "f(spatial, model=spde)", sep = " + ")
            
        } else if(input$RPModelBathy=="spde"){
            
            knots <- seq(min(bathymetry.pred), max(bathymetry.pred), length.out=input$spdeRandomnKnots)
            mesh1d <- inla.mesh.1d(knots)
            
            A.infmesh1d <- inla.spde.make.A(mesh1d, DFsample[,3])
            if(input$autocustomRandomSpde=='custom'){
                theta.prior1 <- as.numericunlist(strsplit(input$RandomTheta1Spde, ","))
                theta.prior2 <- as.numericunlist(strsplit(input$RandomTheta2Spde, ","))
                spde1 <- inla.spde2.matern(mesh1d, theta.prior.mean=c(theta.prior1[1],theta.prior2[1]), 
                                           theta.prior.prec=c(theta.prior1[2],theta.prior2[2]), constr = FALSE)
            } else {spde1 <- inla.spde2.matern(mesh1d, constr = FALSE)}
            
            spde1.index <- inla.spde.make.index(name="sp1", n.spde = spde1$n.spde)
            mod.bathymetry <- "f(sp1, model=spde1)"
            
            A.inf <- list(A.sample, 1, A.infmesh1d)
            
            A.pred <- inla.spde.make.A(mesh = mesh, loc = as.matrix(gridPred))
            A.predmesh1d <- inla.spde.make.A(mesh1d, bathymetry.pred)
            AA.pred <- list(A.pred, 1, A.predmesh1d)
            effects.list <- list(list(spatial=s.index$spatial),
                                 list(Intercept = rep(1, nrow(DFsample)[1])),
                                      list(sp1=spde1.index$sp1))
            Pred.effects.list <- list(list(spatial=s.index$spatial),
                                      list(Intercept = rep(1, nrow(gridPred)[1])),
                                           list(sp1=spde1.index$sp1))
            DFpred <- data.frame(Latitude=gridPred[,1], Longitude=gridPred[,2],
                                 Bathymetry=bathymetry.pred)
            rn.for <- paste("y ~ -1", "Intercept", "f(spatial, model=spde)", mod.bathymetry, sep = " + ")
            
        } else if(input$RPModelBathy=="lin"){
            
            if(input$autocustomLinBathy=='custom'){
                CovariatesPriorParam$mean$Bathymetry <- input$RandomMeanPriorLinBathymetry
                CovariatesPriorParam$prec$Bathymetry <- input$RandomPrecPrioLinrBathymetry
            } else{
                CovariatesPriorParam$mean$Bathymetry <- NA
                CovariatesPriorParam$prec$Bathymetry <- NA
                }
            
            mod.bathymetry <- "Bathymetry"
            effects.list <- list(list(spatial=s.index$spatial),
                                 list(Intercept = rep(1, nrow(DFsample)[1]),
                                      Bathymetry = DFsample[,3]))
            A.inf <- list(A.sample,1)
            
            Pred.effects.list <- list(list(spatial=s.index$spatial),
                                      list(Intercept = rep(1, nrow(gridPred)[1]),
                                           Bathymetry = bathymetry.pred))
            DFpred <- data.frame(Latitude=gridPred[,1], Longitude=gridPred[,2],
                                 Bathymetry=bathymetry.pred)
            
            A.pred <- inla.spde.make.A(mesh = mesh, loc = as.matrix(gridPred))
            AA.pred <- list(A.pred, 1)
            rn.for <- paste("y ~ -1", "Intercept", mod.bathymetry, "f(spatial, model=spde)", sep = " + ")
            
        } else {
            showNotification(ui=paste("Something is wrong with the bathymetric model specification." ), duration = NULL)
            break
        }
        
        
        if(input$autocustomRandomFamily=='custom'){
            controlFamily <- list(hyper = list(prec = list(param = as.numeric(unlist(strsplit(input$RandomFamilyHyper,","))))))
        } else{controlFamily <- inla.set.control.family.default()}
        
        if(input$autocustomRandomMode=='custom'){
            controlModeTheta <- list(theta=as.numeric(unlist(strsplit(input$RandomModeHyper,","))), restart=TRUE)
        } else{controlModeTheta <- inla.set.control.mode.default()}
        
        # if(ncol(DFsample)==4){
        #     effects.list <- list(list(spatial=s.index$spatial),
        #                          list(Intercept = rep(1, nrow(DFsample)[1]),
        #                               Bathymetry = DFsample[,3]))
        # } else if(ncol(DFsample)>4){
        #     effects.list <- list(list(spatial=s.index$spatial),
        #                          list(Intercept = rep(1, nrow(DFsample)[1]),
        #                               Bathymetry = DFsample[,3]))
        #     for(i in 5:ncol(DFsample)){
        #         effects.list[[2]][[names(DFsample)[i]]] <- DFsample[,i]
        #     }
        # }
        # 
        # inference.stack <- inla.stack(data  = list(y = DFsample[,4]),
        #                               A = list(A.sample ,1),
        #                               effects = effects.list,
        #                               tag = "inference")
        # 
        # 
        # 
        # if(ncol(DFsample)==4){
        #     Pred.effects.list <- list(list(spatial=s.index$spatial),
        #                                list(Intercept = rep(1, nrow(gridPred)[1]),
        #                                     Bathymetry = bathymetry.pred))
        #     DFpred <- data.frame(Latitude=gridPred[,1], Longitude=gridPred[,2],
        #                            Bathymetry=bathymetry.pred)
        # 
        #     A.pred <- inla.spde.make.A(mesh = mesh, loc = as.matrix(gridPred))
        #     prediction.stack <- inla.stack(data = list(y = NA),
        #                              A = list(A.pred, 1),
        #                              effects = Pred.effects.list,
        #                              tag = "prediction")
        # } else if(ncol(DFsample)>4){
        #     effects.list <- list(list(spatial=s.index$spatial),
        #                          list(Intercept = rep(1, nrow(DFsample)[1]),
        #                               Batimetria = DFsample[,3]))
        #     for(i in 5:ncol(DFsample)){
        #         effects.list[[2]][[names(DFsample)[i]]] <- DFsample[,i]
        #     }
        # }
        
        inference.stack <- inla.stack(data  = list(y = DFsample[,4]),
                                      A = A.inf,
                                      effects = effects.list,
                                      tag = "inference")
        
        prediction.stack <- inla.stack(data = list(y = NA),
                                       A = AA.pred,
                                       effects = Pred.effects.list,
                                       tag = "prediction")
        
        total.stack <- inla.stack(inference.stack, prediction.stack)

        # random.formula <- as.formula(paste("y","~", "-1 + Intercept", "+ Bathymetry", "+f(spatial, model=spde)"))
        random.formula <- as.formula(rn.for)

        Random.model <- inla(formula= random.formula,#y~-1+Intercept+Bathymetry+f(spatial,model=spde),
                             data = inla.stack.data(total.stack),
                             family = input$SelectRandomFamily,
                             control.inla = list(strategy=input$strategyapproxINLARandom,
                                                 int.strategy=input$strategyintINLARandom),
                             control.predictor = list(A = inla.stack.A(total.stack), compute = TRUE, link = 1),
                             control.fixed = CovariatesPriorParam,
                             control.family = controlFamily,
                             control.mode = controlModeTheta,
                             control.compute = list(cpo = TRUE, dic = TRUE),
                             verbose=FALSE)

        t2 <- Sys.time()

        index.pred <- inla.stack.index(total.stack, "prediction")$data
        DFpred$Abundance.mean <- Random.model$summary.fitted.values[index.pred, "mean"]
        DFpred$Abundance.median <- Random.model$summary.fitted.values[index.pred, "0.5quant"]
        DFpred$Abundance.sd <- Random.model$summary.fitted.values[index.pred, "sd"]

        DFpostFixed <- Random.model$marginals.fixed
        DFpostHyper <- Random.model$marginals.hyperpar
        
        rownames(Random.model$summary.hyperpar)[1] <- ifelse(rownames(Random.model$summary.hyperpar)[1]=="Precision parameter for the Gamma observations", "Precision Gamma", "Precision Gaussian")
        
        DFsummaryFixed <- Random.model$summary.fixed
        DFsummaryHyper <- Random.model$summary.hyperpar
        DFsummaryInternalHyper <- Random.model$internal.summary.hyperpar
        
        gridSpatial <- expand.grid(x=seq(min(DFpred[,1]), max(DFpred[,1]),length.out=input$dimrandommap),
                                   y=seq(min(DFpred[,2]), max(DFpred[,2]), length.out=input$dimrandommap))
        A.spatial <- inla.spde.make.A(mesh=mesh, loc=as.matrix(gridSpatial))
        
        DFspatialMeanMedianStdev <- data.frame(Latitude=as.vector(gridSpatial[,1]), Longitude=as.vector(gridSpatial[,2]), 
                                               Spatial.mean=as.vector(A.spatial%*%Random.model$summary.random$spatial$mean),
                                               Spatial.median=as.vector(A.spatial%*%Random.model$summary.random$spatial$`0.5quant`),
                                               Spatial.stdev=as.vector(A.spatial%*%Random.model$summary.random$spatial$sd))
        
        # DFAbundanceSPGrid <- data.frame(Latitude=as.vector(gridSpatial[,1]), Longitude=as.vector(gridSpatial[,2]), 
        #                                 Abundance.mean=as.vector(A.spatial%*%Random.model$summary.random$spatial$mean),
        #                                 Abundance.median=as.vector(A.spatial%*%Random.model$summary.random$spatial$`0.5quant`),
        #                                 Abundance.stdev=as.vector(A.spatial%*%Random.model$summary.random$spatial$sd))
        
        showNotification(ui=paste("The model has been fitted:", as.numeric(round(Random.model$cpu.used[4])), 
                                  "(abundance model) and", as.numeric(round(difftime(t2,t1, units="secs"))), 
                                  "(overall process) secs." ), duration = NULL)
        
        result <- list(DFpred=list(DFpred=DFpred),
                       DFspatialMeanMedianStdev=list(DFspatialMeanMedianStdev=DFspatialMeanMedianStdev),
                       DFpostFixed=list(DFpostFixed=DFpostFixed),
                       DFpostHyper=list(DFpostHyper=DFpostHyper),
                       DFsummaryFixed=list(DFsummaryFixed=DFsummaryFixed),
                       DFsummaryHyper=list(DFsummaryHyper=DFsummaryHyper),
                       DFsummaryInternalHyper=list(DFsummaryInternalHyper=DFsummaryInternalHyper),
                       SummaryCPO=list(SummaryCPO=na.omit(Random.model$cpo$cpo)),
                       DICmodel=list(DICmodel=Random.model$dic$dic))
        return(result)
    })
    
    
    dataDICtable <- function(){
      DF <- data.frame(DIC=RandomModelFit()$DICmodel$DICmodel) %>%
        mutate(across(where(is.numeric), round, digits = 2))
      return(DF)
    }
    
    downloadableTable("dataDICtable",
                      logger=ss_userAction.Log,
                      filenameroot="dataDICtable",
                      downloaddatafxns=list(csv=dataDICtable,
                                            tsv=dataDICtable),
                      tabledata=dataDICtable, rownames = FALSE,
                      caption="Model DIC")
    
    dataCPOtable <- function(){
      CPO <- RandomModelFit()$SummaryCPO$SummaryCPO
      DF <- data.frame(n=length(CPO), mean=mean(CPO), median=median(CPO), stdev.=sd(CPO), 
                       quantile2.5=quantile(CPO,probs=c(0.025,0.975))[1], 
                       quantile97.5=quantile(CPO,probs=c(0.025,0.975))[2]) %>%
        mutate(across(where(is.numeric), round, digits = 2))
      return(DF)
    }
    
    downloadableTable("dataCPOtable",
                      logger=ss_userAction.Log,
                      filenameroot="dataCPOtable",
                      downloaddatafxns=list(csv=dataCPOtable,
                                            tsv=dataCPOtable),
                      tabledata=dataCPOtable, rownames = FALSE,
                      caption="Summary CPO")
    
    
    dataggplotAbundanceMeanMedianStedevFit <- function(){
        DF <- RandomModelFit()$DFpred$DFpred
        return(DF)
    }
    
    DFRandomAbundance <- reactive({
        DF <- RandomModelFit()$DFpred$DFpred
        condition <- !(length(unique(DF[,1]))*length(unique(DF[,2]))==nrow(DF))&val()
        if(condition){
            grid <- expand.grid(seq(min(DF[,1]),max(DF[,1]), length.out=200),seq(min(DF[,2]),max(DF[,2]), length.out=200))
            DFInter <- data.frame(Latitude=grid[,1],Longitude=grid[,2])
            DFInter$Abundance.mean <- InterpolateIrrGrid(z=DF$Abundance.mean,loc=DF[,1:2], gridInter=grid)$DataInter$z
            DFInter$Abundance.median <- InterpolateIrrGrid(z=DF$Abundance.median,loc=DF[,1:2], gridInter=grid)$DataInter$z
            DFInter$Abundance.sd <- InterpolateIrrGrid(z=DF$Abundance.sd,loc=DF[,1:2], gridInter=grid)$DataInter$z
            DF <- DFInter
        }

        return(DF)
    })

    ggplotAbundanceMeanMedianStedevFit <- function(){
        # DF <- RandomModelFit()$DFpred$DFpred
        DF <- DFRandomAbundance()
        g1 <- ggplot(DF) + geom_tile(aes(x=Latitude, y=Longitude, fill=Abundance.mean)) +
            scale_fill_viridis_c(option = "turbo") + theme_bw() + coord_fixed() + labs(fill="Values") +
            xlab("Latitude") + ylab("Longitude") + ggtitle("Abundance predicted (mean)") +
            theme(plot.title=element_text(color = "black", size = 16, face = "bold", hjust = 0.5))
        
        g2 <- ggplot(DF) + geom_tile(aes(x=Latitude, y=Longitude, fill=Abundance.median)) +
            scale_fill_viridis_c(option = "turbo") + theme_bw() + coord_fixed() + labs(fill="Values") +
            xlab("Latitude") + ylab("Longitude") + ggtitle("Abundance predicted (median)") +
            theme(plot.title=element_text(color = "black", size = 16, face = "bold", hjust = 0.5))
        
        g3 <- ggplot(DF) + geom_tile(aes(x=Latitude, y=Longitude, fill=Abundance.sd))+
            scale_fill_viridis_c(option = "turbo") + theme_bw() + coord_fixed() + labs(fill="Values") +
            xlab("Latitude") + ylab("Longitude") + ggtitle("Abundance predicted (stdev.)") +
            theme(plot.title=element_text(color = "black", size = 16, face = "bold", hjust = 0.5))
        
        gt <- g1+g2+g3
        
        return(gt)
    }

    downloadFile(
        id = "ggplotAbundanceMeanMedianStedevFit",
        logger = ss_userAction.Log,
        filenameroot = "ggplotAbundanceMeanMedianStedevFit",
        aspectratio  = 1,
        downloadfxns = list(png  = ggplotAbundanceMeanMedianStedevFit,
                            csv = dataggplotAbundanceMeanMedianStedevFit,
                            txt = dataggplotAbundanceMeanMedianStedevFit)
    )
    
    downloadablePlot(id = "ggplotAbundanceMeanMedianStedevFit",
                     logger = ss_userAction.Log,
                     filenameroot = "ggplotAbundanceMeanMedianStedevFit",
                     aspectratio  = 1,
                     downloadfxns = list(png=ggplotAbundanceMeanMedianStedevFit,
                                         csv = dataggplotAbundanceMeanMedianStedevFit,
                                         txt = dataggplotAbundanceMeanMedianStedevFit),
                     visibleplot  = ggplotAbundanceMeanMedianStedevFit)
    

    dataggplotSpatialMeanMedianStdev <- function(){
        DF <- RandomModelFit()$DFspatialMeanMedianStdev$DFspatialMeanMedianStdev
        return(DF)
    }
    
    ggplotSpatialMeanMedianStdev <- function(){
        DF <- RandomModelFit()$DFspatialMeanMedianStdev$DFspatialMeanMedianStdev
        g1 <- ggplot(DF) + geom_tile(aes(x=Latitude, y=Longitude, fill=Spatial.mean)) +
            scale_fill_viridis_c(option = "turbo") + theme_bw() + coord_fixed() + labs(fill="Values") +
            xlab("Latitude") + ylab("Longitude") + ggtitle("Posterior spatial (mean)") +
            theme(plot.title=element_text(color = "black", size = 16, face = "bold", hjust = 0.5))
        
        g2 <- ggplot(DF) + geom_tile(aes(x=Latitude, y=Longitude, fill=Spatial.median)) +
            scale_fill_viridis_c(option = "turbo") + theme_bw() + coord_fixed() + labs(fill="Values") +
            xlab("Latitude") + ylab("Longitude") + ggtitle("Posterior spatial (median)") +
            theme(plot.title=element_text(color = "black", size = 16, face = "bold", hjust = 0.5))
        
        g3 <- ggplot(DF) + geom_tile(aes(x=Latitude, y=Longitude, fill=Spatial.stdev))+
            scale_fill_viridis_c(option = "turbo") + theme_bw() + coord_fixed() + labs(fill="Values") +
            xlab("Latitude") + ylab("Longitude") + ggtitle("Posterior spatial (stdev.)") +
            theme(plot.title=element_text(color = "black", size = 16, face = "bold", hjust = 0.5))
        
        gt <- g1+g2+g3
        
        return(gt)
    }
    
    downloadFile(
        id = "ggplotSpatialMeanMedianStdev",
        logger = ss_userAction.Log,
        filenameroot = "ggplotSpatialMeanMedianStdev",
        aspectratio  = 1,
        downloadfxns = list(png  = ggplotSpatialMeanMedianStdev,
                            csv = dataggplotSpatialMeanMedianStdev,
                            txt = dataggplotSpatialMeanMedianStdev)
    )
    
    downloadablePlot(id = "ggplotSpatialMeanMedianStdev",
                     logger = ss_userAction.Log,
                     filenameroot = "ggplotSpatialMeanMedianStdev",
                     aspectratio  = 1,
                     downloadfxns = list(png=ggplotSpatialMeanMedianStdev,
                                         csv = dataggplotSpatialMeanMedianStdev,
                                         txt = dataggplotSpatialMeanMedianStdev),
                     visibleplot  = ggplotSpatialMeanMedianStdev)
    
    dataggplotFixParamFit <- function(){
        DF <- as.data.frame(RandomModelFit()$DFpostFixed$DFpostFixed)
        return(DF)
    }
    
    ggplotFixParamFit <- function(){
        DF <- RandomModelFit()$DFpostFixed$DFpostFixed
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
        id = "ggplotFixParamFit",
        logger = ss_userAction.Log,
        filenameroot = "ggplotFixParamFit",
        aspectratio  = 1,
        downloadfxns = list(png  = ggplotFixParamFit,
                            csv = dataggplotFixParamFit,
                            txt = dataggplotFixParamFit)
    )
    
    downloadablePlot(id = "ggplotFixParamFit",
                     logger = ss_userAction.Log,
                     filenameroot = "ggplotFixParamFit",
                     aspectratio  = 1,
                     downloadfxns = list(png=ggplotFixParamFit,
                                         csv = dataggplotFixParamFit,
                                         txt = dataggplotFixParamFit),
                     visibleplot  = ggplotFixParamFit)
    
    
    dataggplotHyperParamFit <- function(){
        DF <- as.data.frame(RandomModelFit()$DFpostHyper$DFpostHyper)
        return(DF)
    }
    
    ggplotHyperParamFit <- function(){
        DF <- RandomModelFit()$DFpostHyper$DFpostHyper
        gl <- c("g1")
        title <- ifelse(names(DF)[1]=="Precision parameter for the Gamma observations",
                        "Stdev. Gamma", "Stdev. Gaussian")
        g1 <- ggplot(data=data.frame(x=inla.tmarginal(function(x) 1/x**0.5, DF[[1]])[,1], 
                                     y=inla.tmarginal(function(x) 1/x**0.5, DF[[1]])[,2]), aes(x=x,y=y)) + 
            geom_line() + theme_bw() + xlab(title) +
            ylab(HTML(paste("Density f(", title ,")"))) + ggtitle(title) +
            theme(plot.title=element_text(color = "black", size = 16, face = "bold", hjust = 0.5))
        for(i in 2:length(DF)){
          nm <- strsplit(names(DF)[i], " ")[[1]]
          if(nm[1]=="Precision"){
            names(DF)[i] <- paste("Stdev.", paste(nm[-1], collapse=" "), collapse=" ")
            assign(paste0("g",i),
                   ggplot(data=data.frame(x=inla.tmarginal(function(x) 1/x**0.5, DF[[i]])[,1], 
                                          y=inla.tmarginal(function(x) 1/x**0.5, DF[[i]])[,2]), aes(x=x,y=y)) + geom_line() +
                     theme_bw()+ xlab(names(DF)[i]) +
                     ylab(HTML(paste("Density f(",names(DF)[i],")"))) + ggtitle(names(DF)[i]) +
                     theme(plot.title=element_text(color = "black", size = 16, face = "bold", hjust = 0.5))
            )
          }
          else{
            assign(paste0("g",i),
                   ggplot(data=data.frame(x=DF[[i]][,1], y=DF[[i]][,2]), aes(x=x,y=y)) + geom_line() +
                       theme_bw()+ xlab(names(DF)[i]) +
                       ylab(HTML(paste("Density f(",names(DF)[i],")"))) + ggtitle(names(DF)[i]) +
                       theme(plot.title=element_text(color = "black", size = 16, face = "bold", hjust = 0.5))
            )}
            gl[i] <- paste0("g",i)
        }
        gt <- eval(parse(text=paste(gl, collapse="+")))
        return(gt)
    }
    
    downloadFile(
        id = "ggplotHyperParamFit",
        logger = ss_userAction.Log,
        filenameroot = "ggplotHyperParamFit",
        aspectratio  = 1,
        downloadfxns = list(png  = ggplotHyperParamFit,
                            csv = dataggplotHyperParamFit,
                            txt = dataggplotHyperParamFit)
    )
    
    downloadablePlot(id = "ggplotHyperParamFit",
                     logger = ss_userAction.Log,
                     filenameroot = "ggplotHyperParamFit",
                     aspectratio  = 1,
                     downloadfxns = list(png=ggplotHyperParamFit,
                                         csv = dataggplotHyperParamFit,
                                         txt = dataggplotHyperParamFit),
                     visibleplot  = ggplotHyperParamFit)
    
    tableRandomFixedPar <- function(){
        DF <- RandomModelFit()$DFsummaryFixed$DFsummaryFixed %>%
            mutate(across(where(is.numeric), round, digits = 2))
    }
    
    
    downloadableTable("tableRandomFixedPar",
                      logger=ss_userAction.Log,
                      filenameroot="tableRandomFixedPar",
                      downloaddatafxns=list(csv=tableRandomFixedPar,
                                            tsv=tableRandomFixedPar),
                      tabledata=tableRandomFixedPar, rownames = TRUE,
                      caption="Summary fixed parameters")
    
    tableRandomHyperPar <- function(){
        DF <- RandomModelFit()$DFsummaryHyper$DFsummaryHyper %>%
            mutate(across(where(is.numeric), round, digits = 2))
    }
    
    
    downloadableTable("tableRandomHyperPar",
                      logger=ss_userAction.Log,
                      filenameroot="tableRandomHyperPar",
                      downloaddatafxns=list(csv=tableRandomHyperPar,
                                            tsv=tableRandomHyperPar),
                      tabledata=tableRandomHyperPar, rownames = TRUE,
                      caption="Summary hyperparameters")
    
    tableRandomInternalHyperPar <- function(){
        DF <- RandomModelFit()$DFsummaryInternalHyper$DFsummaryInternalHyper %>%
            mutate(across(where(is.numeric), round, digits = 2))
    }
    
    downloadableTable("tableRandomInternalHyperPar",
                      logger=ss_userAction.Log,
                      filenameroot="tableRandomInternalHyperPar",
                      downloaddatafxns=list(csv=tableRandomInternalHyperPar,
                                            tsv=tableRandomInternalHyperPar),
                      tabledata=tableRandomInternalHyperPar, rownames = TRUE,
                      caption="Summary internal hyperparameters")
    

    #############################################################################################
    ##### ##################### ###### Point Process Model Data ###### #################### #####
    #############################################################################################
    
    PrefPointProcessMesh <- reactive({
        if(input$PrefDataSimulatedLoaded=="sim"){
            DFsample <- as.data.frame(Pref.sampling())
        } else if(input$PrefDataSimulatedLoaded=="load"){
            DFsample <- datareadSample()
        }
        if(input$PrefBathymetryRasterSPDE=="raster"){
            rasterSample <- datareadBathymetricRaster()[sample(1:nrow(datareadBathymetricRaster()), min(c(50,nrow(datareadBathymetricRaster())))),1:2]
            qloc <- quantile(as.vector(dist(rasterSample)),probs=c(ifelse(input$PPPCustomMesh,input$PrefPPMeshQloc,0.03),0.3))
            mesh <- inla.mesh.2d(loc=cbind(rasterSample[,1],rasterSample[,2]), cutoff = qloc[1]/2, offset=c(-0.1, -0.4), 
                                  max.edge=c(qloc[1], qloc[2]))
            sample <- rasterSample
        } else if(input$PrefBathymetryRasterSPDE=="solvebathy"){
            qloc <- quantile(as.vector(dist(DFsample[sample(1:nrow(DFsample),size=min(c(50,nrow(DFsample)))),1:2])),probs=c(ifelse(input$PPPCustomMesh,input$PrefPPMeshQloc,0.03),0.3))
            mesh <- inla.mesh.2d(loc=cbind(DFsample[,1],DFsample[,2]), cutoff = qloc[1]/2, offset=c(-0.1, -0.4), 
                                  max.edge=c(qloc[1], qloc[2]))
            sample <- DFsample
        }
        result <- list(mesh=mesh, qloc=qloc, Sample=sample)
        return(result)
    })
    
    ggplotPPPMesh <- function(){
        ggplot()+ gg(PrefPointProcessMesh()$mesh)+ theme_bw() + xlab("Latitude") + ylab("Longitude") +
            ggtitle("Mesh over the study region") + 
            geom_point(data=PrefPointProcessMesh()$Sample, 
                       aes(x=PrefPointProcessMesh()$Sample[,1],y=PrefPointProcessMesh()$Sample[,2]), size=1) + 
            theme(plot.title=element_text(color = "black", size = 16, face = "bold", hjust = 0.5))
    }
    
    downloadFile(
        id = "ggplotPPPMesh",
        logger = ss_userAction.Log,
        filenameroot = "ggplotPPPMesh",
        aspectratio  = 1,
        downloadfxns = list(png  = ggplotPPPMesh)
    )
    
    downloadablePlot(
        id = "ggplotPPPMesh",
        logger = ss_userAction.Log,
        filenameroot = "ggplotPPPMesh",
        aspectratio  = 1,
        downloadfxns = list(png  = ggplotPPPMesh),
        visibleplot  = ggplotPPPMesh)
    
    
    PointProcessModelFit <- eventReactive(input$FitPointProcess, {

        showNotification(ui=paste("Fitting the data."), duration = NULL)
        t1 <- Sys.time()
        #taking the data from simulation or from the loading tab
        if(input$PrefDataSimulatedLoaded=="sim"){
            DFsample <- as.data.frame(Pref.sampling())
        } else if(input$PrefDataSimulatedLoaded=="load"){
            DFsample <- datareadSample()
        }
        
        mesh <- PrefPointProcessMesh()$mesh
        qloc <- PrefPointProcessMesh()$qloc
        
        if(input$optionPPB=='custom'){
            CovariatesPriorParam <- list(mean=list(Intercept=as.numeric(unlist(strsplit(input$RandomMeanPriorBeta,",")))[1]),
                                         prec=list(Intercept=as.numeric(unlist(strsplit(input$RandomMeanPriorBeta,",")))[2]))
        } else if(input$optionPPB=='default'&input$RPModelBathy!="lin"){
            CovariatesPriorParam <- inla.set.control.fixed.default()
        } else {
            CovariatesPriorParam <- list(mean=list(Intercept=NA),
                                         prec=list(Intercept=NA))
        }
        
        prior.range <- c(qloc[2],0.5)
        prior.sigma <- c(1,0.5)
        spde <- inla.spde2.pcmatern(mesh, prior.range = prior.range, prior.sigma = prior.sigma, alpha=2)
        
        ###############################################################
        
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
        y.pp <- rep(0:1, c(nv, n))
        e.pp <- c(w, rep(0, n))
        imat <- Diagonal(nv, rep(1, nv))
        xy <- as.matrix(DFsample[,1:2])
        lmat <- inla.spde.make.A(mesh, xy)
        A.pp <- rbind(imat, lmat)
        
        s.index <- inla.spde.make.index(name="spatial", n.spde = spde$n.spde)
        gridPred <- as.matrix(mesh$loc[,1:2])
        
        ###############################################################
        
        # Now the bathymetric stage to the prediction
        if(input$PrefBathymetryRasterSPDE=="raster"&(input$PrefrasterBatymetryPred=='SPDEraster'|input$PrefrasterBatymetryPred=='rasterpred')){
            DFraster <- datareadBathymetricRaster()
            
            A.raster.inf <- inla.spde.make.A(mesh=mesh, loc=as.matrix(DFraster[,1:2]))
            Bat.raster.inf <- inla.stack(data=list(y=DFraster[,3]),
                                         A = list(A.raster.inf,1),
                                         effects = list(list(spatial=s.index$spatial),
                                                        list(b0 = rep(1, nrow(DFraster)))),
                                         tag="raster.inf")
            
            A.raster.pred <- inla.spde.make.A(mesh = mesh, loc = as.matrix(gridPred))
            Bat.raster.pred <- inla.stack(data = list(y = NA),
                                          A = list(A.raster.pred, 1),
                                          effects = list(list(spatial=s.index$spatial),
                                                         list (b0 = rep(1, dim(gridPred)[1]))),
                                          tag = "raster.pred")
            Bat.total.stack <- inla.stack(Bat.raster.inf, Bat.raster.pred)
            
            batrandom.model <- inla(y ~ -1 + b0 + f(spatial, model = spde),
                                    data = inla.stack.data(Bat.total.stack),
                                    family = "gaussian",
                                    control.inla = list(strategy="simplified.laplace", int.strategy="eb"),
                                    control.predictor = list(A = inla.stack.A(Bat.total.stack), compute = TRUE, link = 1),
                                    control.compute = list(cpo = FALSE, dic = FALSE),
                                    verbose=FALSE)
            index.pred <- inla.stack.index(Bat.total.stack, "raster.pred")$data
            bathymetry.pred <- batrandom.model$summary.fitted.values[index.pred, "mean"]
            bathymetry.pred.sd <- batrandom.model$summary.fitted.values[index.pred, "sd"]
        } else if(input$PrefBathymetryRasterSPDE=="raster"&input$PrefrasterBatymetryPred=='smoothraster'){
            DFraster <- datareadBathymetricRaster()
            bathymetry.pred <- ProjectionFunctionRegularGrid(lim1=range(DFraster[,1]), lim2=range(DFraster[,2]),
                                                             loc=as.matrix(DFraster[,1:2]), z=DFraster[,3], proj.grid=gridPred)
        } else if(input$PrefBathymetryRasterSPDE=="solvebathy"){
            A.raster.inf <- inla.spde.make.A(mesh=mesh, loc=as.matrix(DFsample[,1:2]))
            Bat.raster.inf <- inla.stack(data=list(y=DFsample[,3]),
                                         A = list(lmat ,1),
                                         effects = list(list(spatial=s.index$spatial),
                                                        list(b0 = rep(1, nrow(DFsample)))),
                                         tag="raster.inf")
            
            A.raster.pred <- inla.spde.make.A(mesh = mesh, loc = as.matrix(gridPred))
            Bat.raster.pred <- inla.stack(data = list(y = NA),
                                          A = list(A.raster.pred, 1),
                                          effects = list(list(spatial=s.index$spatial),
                                                         list (b0 = rep(1, dim(gridPred)[1]))),
                                          tag = "raster.pred")
            Bat.total.stack <- inla.stack(Bat.raster.inf, Bat.raster.pred)
            
            batrandom.model <- inla(y ~ -1 + b0 + f(spatial, model = spde),
                                    data = inla.stack.data(Bat.total.stack),
                                    #family = "gaussian",
                                    control.inla = list(strategy="auto", int.strategy="auto"),
                                    control.predictor = list(A = inla.stack.A(Bat.total.stack), compute = TRUE, link = 1),
                                    control.compute = list(cpo = FALSE, dic = FALSE),
                                    verbose=FALSE)
            index.pred <- inla.stack.index(Bat.total.stack, "raster.pred")$data
            bathymetry.pred <- batrandom.model$summary.fitted.values[index.pred, "mean"]
            bathymetry.pred.sd <- batrandom.model$summary.fitted.values[index.pred, "sd"]
        }
        
        if(input$PPModelBathy=="bs"){
            
            bathymetry.pp <- c(bathymetry.pred, DFsample[,3])
            knots <- seq(min(bathymetry.pp), max(bathymetry.pp), length.out=input$bsRandomnKnots)
            mod.bathymetry <- "bs(Bathymetry, knots=knots)"
            
            effects.list <- list(list(spatial=s.index$spatial),
                                 list(Intercept = rep(1, nv + n),
                                      Bathymetry = bathymetry.pp))
            A.inf <- list(A.pp,1)
            rn.for <- paste("y ~ -1", "Intercept", mod.bathymetry, "f(spatial, model=spde)", sep = " + ")
            
        } else if(input$PPModelBathy=="rw1"){
            
            bathymetry.pp <- c(bathymetry.pred, DFsample[,3])
            group.inf <- inla.group(bathymetry.pp, n = input$rw1PrefnKnots, method = "quantile")

            if(input$autocustomPrefRw1=='custom'){
                loggammapar <- as.numeric(unlist(strsplit(input$PrefPrecRw1, split=",")))
                mod.bathymetry <- paste0("f(Bathymetry, model='rw1', hyper = list(prec = list(prior='loggamma',param=c(",
                                         loggammapar[1],",",loggammapar[2],"))))")
            } else{mod.bathymetry <- "f(Bathymetry, model='rw1')"}
            
            effects.list <- list(list(spatial=s.index$spatial),
                                 list(Intercept = rep(1, nv + n),
                                      Bathymetry = group.inf))
            A.inf <- list(A.pp,1)
            rn.for <- paste("y ~ -1", "Intercept", mod.bathymetry, "f(spatial, model=spde)", sep = " + ")
            
        } else if(input$PPModelBathy=="rw2"){
            bathymetry.pp <- c(bathymetry.pred, DFsample[,3])
            group.inf <- inla.group(bathymetry.pp, n = input$rw2PrefnKnots, method = "quantile")

            if(input$autocustomPrefRw2=='custom'){
                loggammapar <- as.numeric(unlist(strsplit(input$RandomPrecRw2, split=",")))
                mod.bathymetry <- paste0("f(Bathymetry, model='rw2', hyper = list(prec = list(prior='loggamma',param=c(",
                                         loggammapar[1],",",loggammapar[2],"))))")
            } else {mod.bathymetry <- "f(Bathymetry, model='rw2')"}

            A.inf <- list(A.pp,1)
            effects.list <- list(list(spatial=s.index$spatial),
                                 list(Intercept = rep(1, nv + n),
                                      Bathymetry = group.inf))
            rn.for <- paste("y ~ -1", "Intercept", mod.bathymetry, "f(spatial, model=spde)", sep = " + ")
            
        } else if(input$PPModelBathy=="spde"){
            bathymetry.pp <- c(bathymetry.pred, DFsample[,3])
            knots <- seq(min(bathymetry.pp), max(bathymetry.pp), length.out=input$spdePrefnKnots)
            mesh1d <- inla.mesh.1d(knots)
            
            A.infmesh1d <- inla.spde.make.A(mesh1d, bathymetry.pp)
            if(input$autocustomRandomSpde=='custom'){
                theta.prior1 <- as.numericunlist(strsplit(input$PrefTheta1Spde, ","))
                theta.prior2 <- as.numericunlist(strsplit(input$PrefTheta2Spde, ","))
                spde1 <- inla.spde2.matern(mesh1d, theta.prior.mean=c(theta.prior1[1],theta.prior2[1]), 
                                           theta.prior.prec=c(theta.prior1[2],theta.prior2[2]), constr = FALSE)
            } else {spde1 <- inla.spde2.matern(mesh1d, constr = FALSE)}
            
            spde1.index <- inla.spde.make.index(name="sp1", n.spde = spde1$n.spde)
            mod.bathymetry <- "f(sp1, model=spde1)"
            
            A.inf <- list(A.pp,1, A.infmesh1d)
            
            effects.list <- list(list(spatial=s.index$spatial),
                                 list(Intercept = rep(1, nv + n)),
                                 list(sp1=spde1.index$sp1))
            
            rn.for <- paste("y ~ -1", "Intercept", "f(spatial, model=spde)", mod.bathymetry, sep = " + ")
            
        } else if(input$PPModelBathy=="lin"){
            bathymetry.pp <- c(bathymetry.pred, DFsample[,3])
            
            if(input$autocustomPrefLinBathy=='custom'){
                CovariatesPriorParam$mean$Bathymetry <- input$PrefMeanPriorLinBathymetry
                CovariatesPriorParam$prec$Bathymetry <- input$PrefPrecPrioLinrBathymetry
            } else{
                CovariatesPriorParam$mean$Bathymetry <- NA
                CovariatesPriorParam$prec$Bathymetry <- NA
            }
            
            mod.bathymetry <- "Bathymetry"
            A.inf <- list(A.pp,1)
            effects.list <- list(list(spatial=s.index$spatial),
                                 list(Intercept = rep(1, nv + n),
                                      Bathymetry = bathymetry.pp))
            rn.for <- paste("y ~ -1", "Intercept", mod.bathymetry, "f(spatial, model=spde)", sep = " + ")
            
        } else {
            showNotification(ui=paste("Something is wrong with the bathymetric model specification." ), duration = NULL)
            break
        }
        
        stk2.pp <- inla.stack(data = list(y = y.pp, e = e.pp),
                              A = A.inf,
                              effects = effects.list,
                              tag = 'pp2')
        
        pp.formula <- as.formula(rn.for)
        
        #Fitting the point process model
        pp.inf <- inla(formula=pp.formula, family = c('poisson'),
                          data = inla.stack.data(stk2.pp),
                          E = inla.stack.data(stk2.pp)$e,
                          control.predictor = list(A = inla.stack.A(stk2.pp))
        )
        result <- list(Summary.fixed=pp.inf$summary.fixed,
                    Summary.hyperpar=pp.inf$summary.hyperpar)
        
        
        t2 <- Sys.time()
        difftime(t2,t1, units="secs")
        showNotification(ui=paste("The model has been fitted:", "it took", as.numeric(round(difftime(t2,t1, units="secs"))), 
                                  "secs." ), duration = NULL)
        
        return(result)
        
    })
    
    tablePointProcessFixedPar <- function(){
        DF <- PointProcessModelFit()$Summary.fixed %>%
            mutate(across(where(is.numeric), round, digits = 2))
    }
    
    downloadableTable("tablePointProcessFixedPar",
                      logger=ss_userAction.Log,
                      filenameroot="tablePointProcessFixedPar",
                      downloaddatafxns=list(csv=tablePointProcessFixedPar,
                                            tsv=tablePointProcessFixedPar),
                      tabledata=tablePointProcessFixedPar, rownames = TRUE,
                      caption="Summary fixed parameters")
    
    tablePointProcessHyperPar <- function(){
        DF <- PointProcessModelFit()$Summary.hyperpar %>%
            mutate(across(where(is.numeric), round, digits = 2))
    }
    
    downloadableTable("tablePointProcessHyperPar",
                      logger=ss_userAction.Log,
                      filenameroot="tablePointProcessHyperPar",
                      downloaddatafxns=list(csv=tablePointProcessHyperPar,
                                            tsv=tablePointProcessHyperPar),
                      tabledata=tablePointProcessHyperPar, rownames = TRUE,
                      caption="Summary hyperparameters")
    

    
    #############################################################################################
    ##### ##################### ###### Preferential Model Data ###### ##################### #####
    #############################################################################################
    
    PreferentialModelFit <- eventReactive(input$fitPref, {
        showNotification(ui=paste("Fitting the data."), duration = NULL)
        t1 <- Sys.time()
        #taking the data from simulation or from the loading tab
        if(input$PrefDataSimulatedLoaded=="sim"){
            DFsample <- as.data.frame(Pref.sampling())
        } else if(input$PrefDataSimulatedLoaded=="load"){
            DFsample <- datareadSample()
        }
        
        mesh <- PrefPointProcessMesh()$mesh
        qloc <- PrefPointProcessMesh()$qloc
        
        if(input$optionPPB=='custom'){
            CovariatesPriorParam <- list(mean=list(Intercept=as.numeric(unlist(strsplit(input$RandomMeanPriorBeta,",")))[1],
                                                   Intercept.pp=as.numeric(unlist(strsplit(input$RandomMeanPriorBeta,",")))[3]),
                                         prec=list(Intercept=as.numeric(unlist(strsplit(input$RandomMeanPriorBeta,",")))[2],
                                                   Intercept.pp=as.numeric(unlist(strsplit(input$RandomMeanPriorBeta,",")))[4]))
        } else if(input$optionPPB=='default'&input$RPModelBathy!="lin"){
            CovariatesPriorParam <- inla.set.control.fixed.default()
        } else {
            CovariatesPriorParam <- list(mean=list(Intercept=NA, Intercept.pp=NA),
                                         prec=list(Intercept=NA, Intercept.pp=NA))
        }
        
        #now we'll make the spde structure
        if(input$optionPPS=="auto"){
            prior.range <- c(qloc[2],0.5)
            prior.sigma <- c(1,0.5)
        } else if(input$optionPPS=="custom"){
            prior.range <- as.numeric(unlist(strsplit(input$PrefPriorRange, split=",")))
            prior.sigma <- as.numeric(unlist(strsplit(input$PrefPriorStdev, split=",")))
        }
        
        prior.range <- c(qloc[2],0.5)
        prior.sigma <- c(1,0.5)
        spde <- inla.spde2.pcmatern(mesh, prior.range = prior.range, prior.sigma = prior.sigma, alpha=2)
        
        ###############################################################
        
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
        y.pp <- rep(0:1, c(nv, n))
        e.pp <- c(w, rep(0, n))
        imat <- Diagonal(nv, rep(1, nv))
        xy <- as.matrix(DFsample[,1:2])
        lmat <- inla.spde.make.A(mesh, xy)
        A.pp <- rbind(imat, lmat)
        
        s.index <- inla.spde.make.index(name="spatial", n.spde = spde$n.spde)
        gridPred <- as.matrix(mesh$loc[,1:2])
        
        ###############################################################
        
        # Now the bathymetric stage to the prediction
        if(input$PrefBathymetryRasterSPDE=="raster"&input$PrefrasterBatymetryPred=='SPDEraster'){
            DFraster <- datareadBathymetricRaster()
            
            gridPredPred <- expand.grid(x=seq(min(DFraster[,1]), max(DFraster[,2]), length.out=input$SPDErasterInput),
                                     y=seq(min(DFraster[,1]), max(DFraster[,2]), length.out=input$SPDErasterInput))
            
            A.raster.inf <- inla.spde.make.A(mesh=mesh, loc=as.matrix(DFraster[,1:2]))
            Bat.raster.inf <- inla.stack(data=list(y=DFraster[,3]),
                                         A = list(A.raster.inf,1),
                                         effects = list(list(spatial=s.index$spatial),
                                                        list(b0 = rep(1, nrow(DFraster)))),
                                         tag="raster.inf")
            
            A.raster.pred <- inla.spde.make.A(mesh = mesh, loc = as.matrix(gridPred))
            Bat.raster.pred <- inla.stack(data = list(y = NA),
                                          A = list(A.raster.pred, 1),
                                          effects = list(list(spatial=s.index$spatial),
                                                         list (b0 = rep(1, dim(gridPred)[1]))),
                                          tag = "raster.pred")
            
            A.raster.predPred <- inla.spde.make.A(mesh = mesh, loc = as.matrix(gridPredPred))
            Bat.raster.predPred <- inla.stack(data = list(y = NA),
                                          A = list(A.raster.predPred, 1),
                                          effects = list(list(spatial=s.index$spatial),
                                                         list (b0 = rep(1, dim(gridPredPred)[1]))),
                                          tag = "raster.predPred")
            
            Bat.total.stack <- inla.stack(Bat.raster.inf, Bat.raster.pred, Bat.raster.predPred)

            batrandom.model <- inla(y ~ -1 + b0 + f(spatial, model = spde),
                                    data = inla.stack.data(Bat.total.stack),
                                    family = "gaussian",
                                    control.inla = list(strategy="simplified.laplace", int.strategy="eb"),
                                    control.predictor = list(A = inla.stack.A(Bat.total.stack), compute = TRUE, link = 1),
                                    control.compute = list(cpo = FALSE, dic = FALSE),
                                    verbose=FALSE)

            index.pred <- inla.stack.index(Bat.total.stack, "raster.pred")$data
            bathymetry.pred <- batrandom.model$summary.fitted.values[index.pred, "mean"]
            bathymetry.pred.sd <- batrandom.model$summary.fitted.values[index.pred, "sd"]
            
            index.predPred <- inla.stack.index(Bat.total.stack, "raster.predPred")$data
            bathymetry.predPred <- batrandom.model$summary.fitted.values[index.predPred, "mean"]
            bathymetry.predPred.sd <- batrandom.model$summary.fitted.values[index.predPred, "sd"]
            
        } else if(input$PrefBathymetryRasterSPDE=="raster"&input$PrefrasterBatymetryPred=='smoothraster'){
            DFraster <- datareadBathymetricRaster()
            bathymetry.pred <- ProjectionFunctionRegularGrid(lim1=range(DFraster[,1]), lim2=range(DFraster[,2]),
                                                             loc=as.matrix(DFraster[,1:2]), z=DFraster[,3], proj.grid=gridPred)
        } else if(input$PrefBathymetryRasterSPDE=="raster"&input$PrefrasterBatymetryPred=='rasterpred'){
            DFraster <- datareadBathymetricRaster()
            bathymetry.predPred <- DFraster[,3]
        } else if(input$PrefBathymetryRasterSPDE=="solvebathy"){
            gridPredPred <- expand.grid(x=seq(min(DFsample[,1]), max(DFsample[,1]), length.out=input$dimprefmap),
                                      y=seq(min(DFsample[,2]), max(DFsample[,2]), length.out=input$dimprefmap))
          
            

            A.raster.inf <- inla.spde.make.A(mesh=mesh, loc=as.matrix(DFsample[,1:2]))
            Bat.raster.inf <- inla.stack(data=list(y=DFsample[,3]),
                                         A = list(lmat ,1),
                                         effects = list(list(spatial=s.index$spatial),
                                                        list(b0 = rep(1, nrow(DFsample)))),
                                         tag="raster.inf")
            
            A.raster.pred <- inla.spde.make.A(mesh = mesh, loc = as.matrix(gridPred))
            Bat.raster.pred <- inla.stack(data = list(y = NA),
                                          A = list(A.raster.pred, 1),
                                          effects = list(list(spatial=s.index$spatial),
                                                         list (b0 = rep(1, dim(gridPred)[1]))),
                                          tag = "raster.pred")
            
            A.raster.predPred <- inla.spde.make.A(mesh = mesh, loc = as.matrix(gridPredPred))
            Bat.raster.predPred <- inla.stack(data = list(y = NA),
                                              A = list(A.raster.predPred, 1),
                                              effects = list(list(spatial=s.index$spatial),
                                                             list (b0 = rep(1, dim(gridPredPred)[1]))),
                                              tag = "raster.predPred")
            
            Bat.total.stack <- inla.stack(Bat.raster.inf, Bat.raster.pred, Bat.raster.predPred)
            

            batrandom.model <- inla(y ~ -1 + b0 + f(spatial, model = spde),
                                    data = inla.stack.data(Bat.total.stack),
                                    #family = "gaussian",
                                    control.inla = list(strategy="auto", int.strategy="auto"),
                                    control.predictor = list(A = inla.stack.A(Bat.total.stack), compute = TRUE, link = 1),
                                    control.compute = list(cpo = FALSE, dic = FALSE),
                                    verbose=FALSE)
            index.pred <- inla.stack.index(Bat.total.stack, "raster.pred")$data
            bathymetry.pred <- batrandom.model$summary.fitted.values[index.pred, "mean"]
            bathymetry.pred.sd <- batrandom.model$summary.fitted.values[index.pred, "sd"]
            
            index.predPred <- inla.stack.index(Bat.total.stack, "raster.predPred")$data
            bathymetry.predPred <- batrandom.model$summary.fitted.values[index.predPred, "mean"]
            bathymetry.predPred.sd <- batrandom.model$summary.fitted.values[index.predPred, "sd"]
        }
        
        if(input$PPModelBathy=="bs"){
          
            bathymetry.pp <- c(bathymetry.pred, DFsample[,3])
            knots <- seq(min(c(bathymetry.pp, bathymetry.predPred)), max(c(bathymetry.pp, bathymetry.predPred)), 
                         length.out=input$bsRandomnKnots)
            
            mod.bathymetry <- "bs(Bathymetry, knots=knots)"
            mod.bathymetry.pp <- "bs(Bathymetry.pp, knots=knots)"
            
            
            A.pp.inf <- list(A.pp,1)
            effects.list.pp <- list(list(spatial.pp=s.index$spatial),
                                    list(Intercept.pp = rep(1, nv + n),
                                         Bathymetry.pp = bathymetry.pp))
            
            A.y.inf <- list(lmat,1)
            Inf.effects.list.y <- list(list(spatial=s.index$spatial),
                                       list(Intercept = rep(1, n),
                                            Bathymetry = DFsample[,3]))
            
            A.pred <- inla.spde.make.A(mesh = mesh, loc = as.matrix(gridPredPred))
            AA.y.pred <- list(A.pred,1)
            Pred.effects.list.y <- list(list(spatial=s.index$spatial),
                                        list(Intercept = rep(1, length(bathymetry.predPred)),
                                             Bathymetry = bathymetry.predPred))
            
            
            gaus.prior <- list(prior = 'gaussian', param = c(0, 0.001))
            rn.for.y <- paste("y ~ -1","(Intercept", mod.bathymetry,  "f(spatial, model=spde))", sep = " + ")
            rn.for.pp <- paste("(Intercept.pp", mod.bathymetry.pp, 
                               "f(spatial.pp, copy = 'spatial', fixed = FALSE, hyper = list(beta = gaus.prior)))", sep = " + ")
            rn.for <- paste(rn.for.y, rn.for.pp, sep = " + ")
            
        } else if(input$PPModelBathy=="rw1"){
            
            bathymetry.pp <- c(bathymetry.pred, DFsample[,3])
            group <- inla.group(c(bathymetry.pp, bathymetry.predPred), n = input$rw1PrefnKnots, method = "quantile")
            group.pp <- group[1:length(bathymetry.pp)]
            group.inf <- group[(length(bathymetry.pred)+1):length(bathymetry.pp)]
            group.pred <- group[-(1:length(bathymetry.pp))]
            
            if(input$autocustomPrefRw1=='custom'){
                loggammapar <- as.numeric(unlist(strsplit(input$PrefPrecRw1, split=",")))
                mod.bathymetry <- paste0("f(Bathymetry, model='rw1', hyper = list(prec = list(prior='loggamma',param=c(",
                                         loggammapar[1],",",loggammapar[2],"))))")
                mod.bathymetry.pp <- paste0("f(Bathymetry.pp, model='rw1')")
            } else{
              mod.bathymetry <- "f(Bathymetry, model='rw1')"
              mod.bathymetry.pp <- "f(Bathymetry.pp, model='rw1')"
            }
            
            A.pp.inf <- list(A.pp,1)
            effects.list.pp <- list(list(spatial.pp=s.index$spatial),
                                    list(Intercept.pp = rep(1, nv + n),
                                         Bathymetry.pp = group.pp))
            
            A.y.inf <- list(lmat,1)
            Inf.effects.list.y <- list(list(spatial=s.index$spatial),
                                       list(Intercept = rep(1, n),
                                            Bathymetry = group.inf))
            
            A.pred <- inla.spde.make.A(mesh = mesh, loc = as.matrix(gridPredPred))
            AA.y.pred <- list(A.pred,1)
            Pred.effects.list.y <- list(list(spatial=s.index$spatial),
                                        list(Intercept = rep(1, length(bathymetry.predPred)),
                                             Bathymetry = group.pred))
            
            DFpred <- data.frame(Latitude=gridPredPred[,1], Longitude=gridPredPred[,2],
                                 Bathymetry=bathymetry.predPred)
            
            gaus.prior <- list(prior = 'gaussian', param = c(0, 0.001))
            rn.for.y <- paste("y ~ -1","(Intercept", mod.bathymetry,  "f(spatial, model=spde))", sep = " + ")
            rn.for.pp <- paste("(Intercept.pp", mod.bathymetry.pp, 
                               "f(spatial.pp, copy = 'spatial', fixed = FALSE, hyper = list(beta = gaus.prior)))", sep = " + ")
            rn.for <- paste(rn.for.y, rn.for.pp, sep = " + ")
            
        } else if(input$PPModelBathy=="rw2"){
          
          bathymetry.pp <- c(bathymetry.pred, DFsample[,3])
          group <- inla.group(c(bathymetry.pp, bathymetry.predPred), n = input$rw2PrefnKnots, method = "quantile")
          group.pp <- group[1:length(bathymetry.pp)]
          group.inf <- group[(length(bathymetry.pred)+1):length(bathymetry.pp)]
          group.pred <- group[-(1:length(bathymetry.pp))]
          
          if(input$autocustomPrefRw2=='custom'){
            loggammapar <- as.numeric(unlist(strsplit(input$PrefPrecRw2, split=",")))
            mod.bathymetry <- paste0("f(Bathymetry, model='rw2', hyper = list(prec = list(prior='loggamma',param=c(",
                                     loggammapar[1],",",loggammapar[2],"))))")
            mod.bathymetry.pp <- paste0("f(Bathymetry.pp, model='rw2')")
          } else {
            mod.bathymetry <- "f(Bathymetry, model='rw2')"
            mod.bathymetry.pp <- "f(Bathymetry.pp, model='rw2')"
          }
          
          A.pp.inf <- list(A.pp,1)
          effects.list.pp <- list(list(spatial.pp=s.index$spatial),
                                  list(Intercept.pp = rep(1, nv + n),
                                       Bathymetry.pp = group.pp))
          
          A.y.inf <- list(lmat,1)
          Inf.effects.list.y <- list(list(spatial=s.index$spatial),
                                     list(Intercept = rep(1, n),
                                          Bathymetry = group.inf))
          
          A.pred <- inla.spde.make.A(mesh = mesh, loc = as.matrix(gridPredPred))
          AA.y.pred <- list(A.pred,1)
          Pred.effects.list.y <- list(list(spatial=s.index$spatial),
                                      list(Intercept = rep(1, length(bathymetry.predPred)),
                                           Bathymetry = group.pred))
          
          DFpred <- data.frame(Latitude=gridPredPred[,1], Longitude=gridPredPred[,2],
                               Bathymetry=bathymetry.predPred)
          
          gaus.prior <- list(prior = 'gaussian', param = c(0, 0.001))
          rn.for.y <- paste("y ~ -1","(Intercept", mod.bathymetry,  "f(spatial, model=spde))", sep = " + ")
          rn.for.pp <- paste("(Intercept.pp", mod.bathymetry.pp, 
                             "f(spatial.pp, copy = 'spatial', fixed = FALSE, hyper = list(beta = gaus.prior)))", sep = " + ")
          rn.for <- paste(rn.for.y, rn.for.pp, sep = " + ")
          
        } else if(input$PPModelBathy=="spde"){
            bathymetry.pp <- c(bathymetry.pred, DFsample[,3])
            knots <- seq(min(c(bathymetry.pp, bathymetry.predPred)), max(c(bathymetry.pp, bathymetry.predPred)), 
                         length.out=input$spdePrefnKnots)
            
            mesh1d <- inla.mesh.1d(knots)
            
            if(input$autocustomRandomSpde=='custom'){
                theta.prior1 <- as.numericunlist(strsplit(input$PrefTheta1Spde, ","))
                theta.prior2 <- as.numericunlist(strsplit(input$PrefTheta2Spde, ","))
                spde1 <- inla.spde2.matern(mesh1d, theta.prior.mean=c(theta.prior1[1],theta.prior2[1]), 
                                           theta.prior.prec=c(theta.prior1[2],theta.prior2[2]), constr = FALSE)
            } else {spde1 <- inla.spde2.matern(mesh1d, constr = FALSE)}
            
            spde1.index <- inla.spde.make.index(name="sp1", n.spde = spde1$n.spde)
            
            mod.bathymetry <- "f(sp1, model=spde1)"
            mod.bathymetry.pp <- "f(sp1.pp, model=spde1)"
            # mod.bathymetry.pp <- "f(sp1.pp, copy='sp1', fixed=FALSE)"
            
            A.infmesh1d.pp <- inla.spde.make.A(mesh1d, bathymetry.pp)
            A.pp.inf <- list(A.pp,1,A.infmesh1d.pp)
            effects.list.pp <- list(list(spatial.pp=s.index$spatial),
                                    list(Intercept.pp = rep(1, nv + n)),
                                    list(sp1.pp=spde1.index$sp1))
            
            A.infmesh1d.y <- inla.spde.make.A(mesh1d, DFsample[,3])
            A.y.inf <- list(lmat,1, A.infmesh1d.y)
            Inf.effects.list.y <- list(list(spatial=s.index$spatial),
                                       list(Intercept = rep(1, n)),
                                       list(sp1=spde1.index$sp1))
            
            A.predmesh1d.y <- inla.spde.make.A(mesh1d, bathymetry.predPred)
            A.y.pred <- inla.spde.make.A(mesh = mesh, loc = as.matrix(gridPredPred))
            AA.y.pred <- list(A.y.pred,1,A.predmesh1d.y)
            Pred.effects.list.y <- list(list(spatial=s.index$spatial),
                                        list(Intercept = rep(1, nrow(gridPredPred))),
                                        list(sp1 = spde1.index$sp1))
            
            DFpred <- data.frame(Latitude=gridPredPred[,1], Longitude=gridPredPred[,2],
                         Bathymetry=bathymetry.predPred)
            
            gaus.prior <- list(prior = 'gaussian', param = input$PrefpriorBetacopy)
            # gaus.prior <- list(prior = 'gaussian', param = ifelse(input$PrefBetaCopy, input$PrefpriorBetacopy, c(0, 0.001)))
            rn.for.y <- paste("y ~ -1", "(Intercept", mod.bathymetry,  "f(spatial, model=spde))", sep = " + ")
            rn.for.pp <- paste("(Intercept.pp", mod.bathymetry.pp, "f(spatial.pp, copy = 'spatial', fixed = FALSE, hyper = list(beta = gaus.prior)))", sep = " + ")
            rn.for <- paste(rn.for.y, rn.for.pp, sep = " + ")
            
            
            ####################################################################################################
            
        } else if(input$PPModelBathy=="lin"){
            

            bathymetry.pp <- c(bathymetry.pred, DFsample[,3])
            
            if(input$autocustomPrefLinBathy=='custom'){
                CovariatesPriorParam$mean$Bathymetry <- input$PrefMeanPriorLinBathymetry
                CovariatesPriorParam$prec$Bathymetry <- input$PrefPrecPrioLinrBathymetry
            } else{
                CovariatesPriorParam$mean$Bathymetry <- NA
                CovariatesPriorParam$prec$Bathymetry <- NA
                CovariatesPriorParam$mean$Bathymetry.pp <- NA
                CovariatesPriorParam$prec$Bathymetry.pp <- NA
            }
            
            mod.bathymetry <- "Bathymetry"
            mod.bathymetry.pp <- "Bathymetry.pp"
            
            A.pp.inf <- list(A.pp,1)
            effects.list.pp <- list(list(spatial.pp=s.index$spatial),
                                 list(Intercept.pp = rep(1, nv + n),
                                      Bathymetry.pp = bathymetry.pp))
            
            A.y.inf <- list(lmat,1)
            
            Inf.effects.list.y <- list(list(spatial=s.index$spatial),
                                    list(Intercept = rep(1, n),
                                         Bathymetry = DFsample[,3]))
            
            Pred.effects.list.y <- list(list(spatial=s.index$spatial),
                                      list(Intercept = rep(1, nrow(gridPredPred)),
                                           Bathymetry = bathymetry.predPred))
            DFpred <- data.frame(Latitude=gridPredPred[,1], Longitude=gridPredPred[,2],
                                 Bathymetry=bathymetry.predPred)
            
            A.y.pred <- inla.spde.make.A(mesh = mesh, loc = as.matrix(gridPredPred))
            AA.y.pred <- list(A.y.pred,1)
            
            gaus.prior <- list(prior = 'gaussian', param = input$PrefpriorBetacopy)
            # gaus.prior <- list(prior = 'gaussian', param = ifelse(input$PrefBetaCopy, input$PrefpriorBetacopy, c(0, 0.001)))
            rn.for.y <- paste("y ~ -1", "(Intercept", mod.bathymetry, "f(spatial, model=spde))", sep = " + ")
            rn.for.pp <- paste("(Intercept.pp", mod.bathymetry.pp, "f(spatial.pp, copy = 'spatial', fixed = FALSE, hyper = list(beta = gaus.prior)))", sep = " + ")
            rn.for <- paste(rn.for.y, rn.for.pp, sep = " + ")
            
        } else {
            showNotification(ui=paste("Something is wrong with the bathymetric model specification." ), duration = NULL)
            break
        }
        

        stk2.inf.y <- inla.stack(data = list(y = cbind(DFsample[,4], NA), e = rep(0, n)),
                                 A=A.y.inf, effects=Inf.effects.list.y, tag='inf.y')
        
        
        stk2.pred.y <- inla.stack(data=list(y=matrix(NA, nrow=nrow(gridPredPred), ncol=2)),
                                  A=AA.y.pred, effects=Pred.effects.list.y, tag='pred.y')
        
        stk2.pp <- inla.stack(data = list(y = cbind(NA, y.pp), e = e.pp),
                              A = A.pp.inf, effects = effects.list.pp, tag = 'inf.pp')
        
        stk2.total <- inla.stack(stk2.inf.y, stk2.pred.y, stk2.pp)
        
        formula.pref <- as.formula(rn.for)
        
        if(input$autocustomRandomFamily=='custom'){
          controlFamily <- list(list(hyper = list(prec = list(param = as.numeric(unlist(strsplit(input$PrefFamilyHyper,",")))))), list())
        } else{controlFamily <- list(list(), list())
        #inla.set.control.family.default()
        }
        
        if(input$autocustomRandomMode=='custom'){
          controlModeTheta <- list(theta=as.numeric(unlist(strsplit(input$PrefModeHyper,","))), restart=TRUE)
        } else{controlModeTheta <- inla.set.control.mode.default()}
        
        #Fitting the point process model
        Preferential.model <- inla(formula=formula.pref, family = c(input$SelectPrefFamily,'poisson'),
                                   data = inla.stack.data(stk2.total),
                                   E = inla.stack.data(stk2.total)$e,
                                   control.inla = list(strategy=input$strategyapproxINLAPref,
                                                       int.strategy=input$strategyintINLAPref),
                                   control.predictor = list(A = inla.stack.A(stk2.total), compute = TRUE, link = 1),
                                   control.fixed = CovariatesPriorParam,
                                   control.family = controlFamily,
                                   control.mode = controlModeTheta,
                                   control.compute = list(cpo = TRUE, dic = TRUE),
                                   verbose=FALSE)
                       
        rownames(Preferential.model$summary.hyperpar)[1] <- ifelse(rownames(Preferential.model$summary.hyperpar)[1]=="Precision parameter for the Gamma observations", "Precision Gamma", "Precision Gaussian")
        
        index.pred <- inla.stack.index(stk2.total, "pred.y")$data
        DFpred$Abundance.mean <- Preferential.model$summary.fitted.values[index.pred, "mean"]
        DFpred$Abundance.median <- Preferential.model$summary.fitted.values[index.pred, "0.5quant"]
        DFpred$Abundance.sd <- Preferential.model$summary.fitted.values[index.pred, "sd"]
        
        gridSpatial <- expand.grid(x=seq(min(DFpred[,1]), max(DFpred[,1]),length.out=input$dimprefmap),
                                    y=seq(min(DFpred[,2]), max(DFpred[,2]), length.out=input$dimprefmap))
        A.spatial <- inla.spde.make.A(mesh=mesh, loc=as.matrix(gridSpatial))
        
        DFspatialMeanMedianStdev <- data.frame(Latitude=as.vector(gridSpatial[,1]), Longitude=as.vector(gridSpatial[,2]), 
                                               Spatial.mean=as.vector(A.spatial%*%Preferential.model$summary.random$spatial$mean),
                                               Spatial.median=as.vector(A.spatial%*%Preferential.model$summary.random$spatial$`0.5quant`),
                                               Spatial.stdev=as.vector(A.spatial%*%Preferential.model$summary.random$spatial$sd))
        
        result <- list(DFpredAbunMeanMedianStdev=list(DFpredAbunMeanMedianStdev=DFpred),
                       DFspatialMeanMedianStdev=list(DFspatialMeanMedianStdev=DFspatialMeanMedianStdev),
                       DFPostFixed=list(DFPostFixed=Preferential.model$marginals.fixed),
                       DFPostHyperpar=list(DFPostHyperpar=Preferential.model$marginals.hyperpar),
                       Summary.fixed=list(Summary.fixed=Preferential.model$summary.fixed),
                       Summary.hyperpar=list(Summary.hyperpar=Preferential.model$summary.hyperpar),
                       SummaryInternalHyper=list(SummaryInternalHyper=Preferential.model$internal.summary.hyperpar),
                       SummaryCPO=list(SummaryCPO=na.omit(Preferential.model$cpo$cpo)),
                       DICmodel=list(DICmodel=data.frame(DIC=Preferential.model$dic$family.dic, row.names=c("Geostatistical", "Point process"))))
        
        
        t2 <- Sys.time()
        difftime(t2,t1, units="secs")
        showNotification(ui=paste("The model has been fitted:", as.numeric(round(Preferential.model$cpu.used[4])), 
                                  "(abundance model) and", as.numeric(round(difftime(t2,t1, units="secs"))), 
                                  "(overall process) secs." ), duration = NULL)
        # showNotification(ui=paste("The model's DIC is", Preferential.model$dic$dic), duration = NULL)
        return(result)
    })
    
    
    dataggplotPrefAbundanceMeanMedianStedevFit <- function(){
      DF <- PreferentialModelFit()$DFpredAbunMeanMedianStdev$DFpredAbunMeanMedianStdev
      return(DF)
    }
    
    DFPrefAbundance <- reactive({
      DF <- PreferentialModelFit()$DFpredAbunMeanMedianStdev$DFpredAbunMeanMedianStdev
      condition <- !(length(unique(DF[,1]))*length(unique(DF[,2]))==nrow(DF))&val()
      if(condition){
        grid <- expand.grid(seq(min(DF[,1]),max(DF[,1]), length.out=200),seq(min(DF[,2]),max(DF[,2]), length.out=200))
        DFInter <- data.frame(Latitude=grid[,1],Longitude=grid[,2])
        DFInter$Abundance.mean <- InterpolateIrrGrid(z=DF$Abundance.mean,loc=DF[,1:2], gridInter=grid)$DataInter$z
        DFInter$Abundance.median <- InterpolateIrrGrid(z=DF$Abundance.median,loc=DF[,1:2], gridInter=grid)$DataInter$z
        DFInter$Abundance.sd <- InterpolateIrrGrid(z=DF$Abundance.sd,loc=DF[,1:2], gridInter=grid)$DataInter$z
        DF <- DFInter
      }
      
      return(DF)
    })
    
    ggplotPrefAbundanceMeanMedianStedevFit <- function(){
      # DF <- PreferentialModelFit()$DFpredAbunMeanMedianStdev$DFpredAbunMeanMedianStdev
      DF <- DFPrefAbundance()
      g1 <- ggplot(DF) + geom_tile(aes(x=Latitude, y=Longitude, fill=Abundance.mean)) +
        scale_fill_viridis_c(option = "turbo") + theme_bw() + coord_fixed() + labs(fill="Values") +
        xlab("Latitude") + ylab("Longitude") + ggtitle("Abundance predicted (mean)") +
        theme(plot.title=element_text(color = "black", size = 16, face = "bold", hjust = 0.5))
      
      g2 <- ggplot(DF) + geom_tile(aes(x=Latitude, y=Longitude, fill=Abundance.median)) +
        scale_fill_viridis_c(option = "turbo") + theme_bw() + coord_fixed() + labs(fill="Values") +
        xlab("Latitude") + ylab("Longitude") + ggtitle("Abundance predicted (median)") +
        theme(plot.title=element_text(color = "black", size = 16, face = "bold", hjust = 0.5))
      
      g3 <- ggplot(DF) + geom_tile(aes(x=Latitude, y=Longitude, fill=Abundance.sd))+
        scale_fill_viridis_c(option = "turbo") + theme_bw() + coord_fixed() + labs(fill="Values") +
        xlab("Latitude") + ylab("Longitude") + ggtitle("Abundance predicted (stdev.)") +
        theme(plot.title=element_text(color = "black", size = 16, face = "bold", hjust = 0.5))
      
      gt <- g1+g2+g3
      
      return(gt)
    }
    
    downloadFile(
      id = "ggplotPrefAbundanceMeanMedianStedevFit",
      logger = ss_userAction.Log,
      filenameroot = "ggplotPrefAbundanceMeanMedianStedevFit",
      aspectratio  = 1,
      downloadfxns = list(png  = ggplotPrefAbundanceMeanMedianStedevFit,
                          csv = dataggplotPrefAbundanceMeanMedianStedevFit,
                          txt = dataggplotPrefAbundanceMeanMedianStedevFit)
    )
    
    downloadablePlot(id = "ggplotPrefAbundanceMeanMedianStedevFit",
                     logger = ss_userAction.Log,
                     filenameroot = "ggplotPrefAbundanceMeanMedianStedevFit",
                     aspectratio  = 1,
                     downloadfxns = list(png=ggplotPrefAbundanceMeanMedianStedevFit,
                                         csv = dataggplotPrefAbundanceMeanMedianStedevFit,
                                         txt = dataggplotPrefAbundanceMeanMedianStedevFit),
                     visibleplot  = ggplotPrefAbundanceMeanMedianStedevFit)
    
    
    dataggplotPrefSpatialMeanMedianStdev <- function(){
      DF <- PreferentialModelFit()$DFspatialMeanMedianStdev$DFspatialMeanMedianStdev
      return(DF)
    }
    
    ggplotPrefSpatialMeanMedianStdev <- function(){
      DF <- PreferentialModelFit()$DFspatialMeanMedianStdev$DFspatialMeanMedianStdev
      g1 <- ggplot(DF) + geom_tile(aes(x=Latitude, y=Longitude, fill=Spatial.mean)) +
        scale_fill_viridis_c(option = "turbo") + theme_bw() + coord_fixed() + labs(fill="Values") +
        xlab("Latitude") + ylab("Longitude") + ggtitle("Posterior spatial (mean)") +
        theme(plot.title=element_text(color = "black", size = 16, face = "bold", hjust = 0.5))
      
      g2 <- ggplot(DF) + geom_tile(aes(x=Latitude, y=Longitude, fill=Spatial.median)) +
        scale_fill_viridis_c(option = "turbo") + theme_bw() + coord_fixed() + labs(fill="Values") +
        xlab("Latitude") + ylab("Longitude") + ggtitle("Posterior spatial (median)") +
        theme(plot.title=element_text(color = "black", size = 16, face = "bold", hjust = 0.5))
      
      g3 <- ggplot(DF) + geom_tile(aes(x=Latitude, y=Longitude, fill=Spatial.stdev))+
        scale_fill_viridis_c(option = "turbo") + theme_bw() + coord_fixed() + labs(fill="Values") +
        xlab("Latitude") + ylab("Longitude") + ggtitle("Posterior spatial (stdev.)") +
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
    
    dataggplotPrefFixParamFit <- function(){
      DF <- as.data.frame(PreferentialModelFit()$DFPostFixed$DFPostFixed)
      return(DF)
    }
    
    ggplotPrefFixParamFit <- function(){
      DF <- PreferentialModelFit()$DFPostFixed$DFPostFixed
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
                     downloadfxns = list(png=ggplotPrefFixParamFit,
                                         csv = dataggplotPrefFixParamFit,
                                         txt = dataggplotPrefFixParamFit),
                     visibleplot  = ggplotPrefFixParamFit)
    
    
    dataggplotPrefHyperParamFit <- function(){
      DF <- as.data.frame(PreferentialModelFit()$DFPostHyperpar$DFPostHyperpar)
      return(DF)
    }
    
    ggplotPrefHyperParamFit <- function(){
      DF <- PreferentialModelFit()$DFPostHyperpar$DFPostHyperpar
      gl <- c("g1")
      title <- ifelse(names(DF)[1]=="Precision parameter for the Gamma observations",
                      "Stdev. Gamma", "Stdev. Gaussian")
      g1 <- ggplot(data=data.frame(x=inla.tmarginal(function(x) 1/x**0.5, DF[[1]])[,1], 
                                   y=inla.tmarginal(function(x) 1/x**0.5, DF[[1]])[,2]), aes(x=x,y=y)) + 
        geom_line() + theme_bw() + xlab(title) +
        ylab(HTML(paste("Density f(", title ,")"))) + ggtitle(title) +
        theme(plot.title=element_text(color = "black", size = 16, face = "bold", hjust = 0.5))
      for(i in 2:length(DF)){
        nm <- strsplit(names(DF)[i], " ")[[1]]
        if(nm[1]=="Precision"){
          names(DF)[i] <- paste("Stdev.", paste(nm[-1], collapse=" "), collapse=" ")
          assign(paste0("g",i),
                 ggplot(data=data.frame(x=inla.tmarginal(function(x) 1/x**0.5, DF[[i]])[,1], 
                                        y=inla.tmarginal(function(x) 1/x**0.5, DF[[i]])[,2]), aes(x=x,y=y)) + geom_line() +
                   theme_bw()+ xlab(names(DF)[i]) +
                   ylab(HTML(paste("Density f(",names(DF)[i],")"))) + ggtitle(names(DF)[i]) +
                   theme(plot.title=element_text(color = "black", size = 16, face = "bold", hjust = 0.5))
          )
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
      DF <- PreferentialModelFit()$Summary.fixed$Summary.fixed %>%
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
      DF <- PreferentialModelFit()$Summary.hyperpar$Summary.hyperpar %>%
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
      DF <- PreferentialModelFit()$SummaryInternalHyper$SummaryInternalHyper %>%
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
      DF <- PreferentialModelFit()$DICmodel$DICmodel %>% # data.frame(DIC=PreferentialModelFit()$DICmodel$DICmodel) %>%
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
    
    dataPrefCPOtable <- function(){
      CPO <- PreferentialModelFit()$SummaryCPO$SummaryCPO
      DF <- data.frame(n=length(CPO), mean=mean(CPO), median=median(CPO), stdev.=sd(CPO), 
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
    
    
    
    
    
    
    periscope:::fw_server_setup(input, output, session, ss_userAction.Log)
    
    loginfo("%s started with log level <%s>",
            periscope:::fw_get_title(), periscope:::fw_get_loglevel(),
            logger = ss_userAction.Log)
})
