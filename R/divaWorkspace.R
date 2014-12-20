setClass("divaWorkspace", contains = "flowJoWorkspace")

setMethod("getSamples","divaWorkspace",function(x){
      selectMethod("getSampleGroups","divaWorkspace")(x)
    })

setMethod("getSampleGroups","divaWorkspace",function(x){
      .getSampleGroupsDiva(x)
    })
.getSampleGroupsDiva<-function(x){
    ldply(
        xpathApply(x@doc, "/bdfacs/experiment/specimen",function(specimen){
              samples <- xpathApply(specimen, "tube",function(tube){
                                            c(tube = xmlGetAttr(tube,"name")
                                                , sampleName = xmlValue(xmlElementsByTagName(tube,"data_filename")[[1]])
                                            )
                                  })
                              
              samples <- ldply(samples)
              samples[["specimen"]] <- xmlGetAttr(specimen, "name")
              samples
            })
      )
  
}
    
setMethod("show",c("divaWorkspace"),function(object){
      cat("Diva Workspace Version ",object@version,"\n");
      cat("File location: ",object@path,"\n");
      cat("File name: ",object@file,"\n");
      if(object@.cache$flag){
        cat("Workspace is open.","\n");
        cat("\nGroups in Workspace\n");
        
        sg <- getSampleGroups(object)
        
        tbl <- data.frame(table(sg$specimen))
        colnames(tbl) <- c("specimen", "samples")
        print(tbl)
      }else{	
        cat("Workspace is closed.","\n")
      }
    })

setMethod("parseWorkspace",signature("divaWorkspace"),function(obj, ...){
      .preprocessorDiva(obj, ...)
    })

.preprocessorDiva<- function(obj, name = NULL
                                    , subset = NULL
                                    , path = obj@path
                                    , ...)
{	
  	
  #sample info  
  sg <- getSamples(obj)
  
  # filter by group name
  sg[["specimen"]] <- factor(sg[["specimen"]])
  groups <- levels(sg[["specimen"]])
  
  if(is.null(name)){
    message("Choose which group of samples to import:\n");
    groupInd <- menu(groups,graphics=FALSE);
  }else if(is.numeric(name)){
    if(length(groups)<name)
      stop("Invalid sample group index.")
    groupInd <- name
  }else if(is.character(name)){
    if(is.na(match(name,groups)))
      stop("Invalid sample group name.")
    groupInd <- match(name,groups)
  }
  group.name <- groups[groupInd]
  
  sg <- subset(sg, specimen == group.name)
#    browser()    
  #filter by subset (sample name or numeric index)
  if(is.factor(subset)){
    subset<-as.character(subset)
  }
  if(is.character(subset)){
    sg <- subset(sg, name %in% subset)
  }else if(is.numeric(subset))
    sg <- sg[subset, ] 
  
  #check if there are samples to parse
  sn <- sg[["sampleName"]]
  nSample <- length(sn)
  if(nSample == 0)
    stop("No samples in this workspace to parse!")
  
  
  #check duplicated sample names
  
  isDup <- duplicated(sn)
  if(any(isDup))
    stop("Duplicated sample names detected within group: ", paste(sampleSelected[isDup], collapse = " "), "\n Please check if the appropriate group is selected.")
  
  
  message("Parsing ", nSample," samples");
  .parseDivaWorkspace(xmlFileName=file.path(obj@path,obj@file)
                      ,samples = sn
                      , groupName = group.name
                      ,path=path
                      ,xmlParserOption = obj@options
                      ,ws = obj
                      ,...)
  
}

.parseDivaWorkspace <- function(xmlFileName,samples,path,xmlParserOption, ws, groupName,...){
  
  if(!file.exists(xmlFileName))
    stop(xmlFileName," not found!")
#  gs <- new("GatingSet", guid = .uuid_gen(), flag = FALSE)
  
  
  dataPaths <- vector("character")
  excludefiles<-vector("logical")
  for(file in samples){
    
    #########################################################
    #get full path for each fcs
    #########################################################
    ##escape "illegal" characters
    file<-gsub("\\?","\\\\?",gsub("\\]","\\\\]",gsub("\\[","\\\\[",gsub("\\-","\\\\-",gsub("\\+","\\\\+",gsub("\\)","\\\\)",gsub("\\(","\\\\(",file)))))))
    absPath <- list.files(pattern=paste("^",file,"",sep=""),path=path,recursive=TRUE,full.names=TRUE)
    nFound <- length(absPath)
    if(nFound == 0){
      warning("Can't find ",file," in directory: ",path,"\n");
      excludefiles<-c(excludefiles,TRUE);
      
    }else if(nFound > 1){
      stop('Multiple files found for:', file) 
    }else{
      
      dataPaths<-c(dataPaths,dirname(absPath))
      excludefiles<-c(excludefiles,FALSE);
    }
  }
  #Remove samples where files don't exist.
  if(length(which(excludefiles))>0){
    message("Removing ",length(which(excludefiles))," samples from the analysis since we can't find their FCS files.");
    samples<-samples[!excludefiles];
  }
  
  
  files<-file.path(dataPaths,samples)

  if(length(files)==0)
    stop("not sample to be added to GatingSet!")

  #load the raw data from FCS
  fs <- read.ncdfFlowSet(files,isWriteSlice=FALSE,...)

  rootDoc <- ws@doc
  # get comp (assuming identical for entire experiment)
  browser()       
  comp <- xpathApply(groupNode, "/")
  xpathLogParam <- paste0("/bdfacs/experiment/instrument_settings/parameter[is_log='true']")
  logParamNodes <- xpathApply(rootDoc, xpathLogParam)
  
  axis <- lapply(files,function(file){
        
        
        sampleName <- basename(file)
        browser()       
        #get tube node
        xpathGroup <- paste0("/bdfacs/experiment/specimen[@name='", groupName, "']")
        xpathSample <- paste0(xpathGroup, "/tube[data_filename='", sampleName, "']")
        sampleNode <- xpathApply(rootDoc, xpathSample)
         
        
        
        ##################################
        #Compensating the data
        ##################################
        if(execute)
        {
          
          cnd <- colnames(fs)
          message("loading data: ",file);
          if(isNcdf)
            data <- read.FCS(file)[, cnd]
          else
            data <- fs[[sampleName]]
          
          #alter colnames(replace "/" with "_") for flowJo X
          if(wsType == "vX"){
            new_cnd <- gsub("/","_",cnd)
            if(!all(new_cnd == cnd)){ #check if needs to update colnames to avoid unneccessary expensive colnames<- call
              cnd <- new_cnd
              colnames(data) <- cnd
              
            }
            
          }
          
          
          
          if(cid=="")
            cid=-2
          
          if(cid!="-1" && cid!="-2"){
            message("Compensating");
            
            marker <- comp$parameters
            
            if(is.null(compensation)){
              ## try to match marker from comp with flow data in case flowJo is not consistent with data
              markerInd <- sapply(marker, function(thisMarker)grep(thisMarker, cnd, ignore.case = ignore.case))
              matchedMarker <- cnd[markerInd]
              if(length(matchedMarker) != length(marker))
                stop("channels mismatched between compensation and flow data!")
              marker <- matchedMarker
              
              compobj <- compensation(matrix(comp$spillOver,nrow=length(marker),ncol=length(marker),byrow=TRUE,dimnames=list(marker,marker)))
            }else
              compobj <- compensation#TODO: to update compensation information in C part
            #TODO this compensation will fail if the parameters have <> braces (meaning the data is stored compensated).
            #I need to handle this case properly.
            
            res <- try(compensate(data,compobj),silent=TRUE)
            if(inherits(res,"try-error")){
              message("Data is probably stored already compensated");
            }else{
              data <- res
              rm(res);
            }
            
          }
          else if(cid=="-2"){
            #TODO the matrix may be acquisition defined.
            message("No compensation");
          }
          else if(cid=="-1")
          {
            ##Acquisition defined compensation.
            nm <- comp$comment
            
            
            if(grepl("Acquisition-defined",nm)){
              ###Code to compensate the sample using the acquisition defined compensation matrices.
              message("Compensating with Acquisition defined compensation matrix");
              #browser()
              if(is.null(compensation))
              {
                compobj <- compensation(spillover(data)$SPILL)
                
              }else
              {
                compobj <- compensation
                
              }
              
              res <- try(compensate(data,compobj),silent=TRUE)
              
              if(inherits(res,"try-error")){
                message("Data is probably stored already compensated");
              }else{
                data<-res
                rm(res);
                
              }
              
            }
            
          }
        }else{
          # get kw from ws (We are not sure yet if this R API will always 
          # return keywords from workspace successfully, thus it is currently
          # only used when execute = FALSE
          if(!is.null(ws))
            kw <- getKeywords(ws, sampleName, sampNloc = sampNloc)
          #use $PnB to determine the number of parameters since {N,R,S) could be
          #redundant in some workspaces
          key_names <- unique(names(kw[grep("\\$P[0-9]{1,}B", names(kw))]))
          key_names <- gsub("B", "N", key_names, fixed = TRUE)
          cnd <- as.vector(unlist(kw[key_names]))   
          
        }
        
        ##################################
        #alter the colnames
        ##################################
        if(cid!="-2")
        {
          
          #get prefix if it is not set yet
          if(is.null(prefixColNames)&&prefix){
            
            if(is.null(cnd)){
              cnd <- as.vector(parameters(data)@data$name)
            }
            prefixColNames <- cnd
            if(execute)
              comp_param <- parameters(compobj)
            else
              comp_param <- comp$parameters
            
            wh <- match(comp_param, prefixColNames)
            
            prefixColNames[wh] <- paste(comp$prefix,comp_param,comp$suffix,sep="")
            
            
            
          }
        }else{
          prefixColNames <- cnd
        }
        ##################################
        #transforming and gating
        ##################################
        if(execute)
        {
          
          message(paste("gating ..."))
          #stop using gating API of cdf-version because c++ doesn't store the view of ncdfFlowSet anymore
          mat <- data@exprs #using @ is faster than exprs()
          #get gains from keywords
          # for now we still parse it from data
          # once confirmed that workspace is a reliable source for this info
          # we can parse it from ws as well
          this_pd <- pData(parameters(data))
          paramIDs <- rownames(pData(parameters(data)))
          key_names <- paste(paramIDs,"G",sep="")
          kw <- keyword(data)
          if(as.numeric(kw[["FCSversion"]])>=3){
            kw_gains <- kw[key_names]
            
            # For keywords where the gain is not set, the gain is NULL.
            # We replace these instances with the default of 1.
            kw_gains[sapply(kw_gains, is.null)] <- 1
            
            gains <- as.numeric(kw_gains)
          }else{
            gains <- rep(1,length(paramIDs))
          }
          
          names(gains) <- this_pd$name
          
          #update colnames in order for the gating to find right dims
          if(!is.null(prefixColNames)){
            dimnames(mat) <- list(NULL, prefixColNames)
          }
          recompute <- FALSE
          nodeInd <- 0
          
          .Call("R_gating",G@pointer, mat, sampleName, gains, nodeInd, recompute, extend_val, ignore.case, leaf.bool)
#            browser()
          #restore the non-prefixed colnames for updating data in fs with [[<-
          #since colnames(fs) is not udpated yet.
          if(!is.null(prefixColNames)){
            #restore the orig colnames(replace "_" with "/") for flowJo X
            if(wsType == "vX"){
              old_cnd <- gsub("_","/",cnd)
              if(!all(old_cnd == cnd)){ #check if needs to update colnames to avoid unneccessary expensive colnames<- call
                cnd <- old_cnd
                colnames(data) <- cnd #restore colnames for flowFrame as well for flowJo vX  
              }
              
            }
            dimnames(mat) <- list(NULL, cnd)
          }
          
          data@exprs <- mat #circumvent the validity check of exprs<- to speed up
          if(isNcdf){
            fs[[sampleName]] <- data
            
          }else{
            assign(sampleName,data,fs@frames)
          }
          #range info within parameter object is not always the same as the real data range
          #it is used to display the data.
          #so we need update this range info by transforming it
          tInd <- grepl("[Tt]ime",cnd)
          if(any(tInd))
            tRg  <- range(mat[,tInd])
          else
            tRg <- NULL
          axis.labels <- .transformRange(G,sampleName,wsType,fs@frames,timeRange = tRg)
          
        }else{
          #extract gains from keyword of ws
          #currently it is only used for extracting gates without gating
          #In future we want to use it for gating as well
          #once we have confirmed that ws is a reliable source of keyword
#          browser()
          
          #get gains from keywords
          kw_gains <- grep("P[0-9]{1,}G", names(kw))
          
          if(length(kw_gains) > 0){
            key_names <- unique(names(kw[kw_gains]))
            kw_gains <- kw[key_names]
            
            # For keywords where the gain is not set, the gain is NULL.
            # We replace these instances with the default of 1.
            kw_gains[sapply(kw_gains, is.null)] <- 1
            
            gains <- as.numeric(kw_gains)                      
          }else{
            gains <- rep(1,length(cnd))
          }
          
          
          names(gains) <- prefixColNames
          #transform and adjust the gates without gating
          .Call("R_computeGates",G@pointer, sampleName, gains, extend_val, extend_to)
          axis.labels <- list()
        }
        
        #set global variable
        tempenv$prefixColNames <- prefixColNames
        
        #return axis.labels
        axis.labels
      })
  
  names(axis) <- basename(files)
  G@axis <- axis
  G@flag <- execute #assume the excution would succeed if the entire G gets returned finally
  
  if(execute)
  {
#		browser()
    #update colnames
    #can't do it before fs fully compensated since
    #compensate function check the consistency colnames between input flowFrame and fs
    if(!is.null(tempenv$prefixColNames))
      colnames(fs) <- tempenv$prefixColNames
    
    #attach filename and colnames to internal stucture for gating
#		browser()
  }
  
  flowData(G) <- fs
  G
  
  
  
  message("done!")
  
  
  
  #we don't want to return the splitted gs since they share the same cdf and externalptr
  #thus should be handled differently(more efficiently) from the regular gslist
  
#    # try to post process the GatingSet to split the GatingSets(based on different the gating trees) if needed                
  gslist <- suppressMessages(.groupByTree(G))
  if(length(gslist) > 1)
    warning("GatingSet contains different gating tree structures and must be cleaned before using it! ")
#    if(length(gslist) == 1){
#      return (gslist[[1]])      
#    }else
  {
#      warning("Due to the different gating tree structures, a list of GatingSets is returned instead!")
#      return (gslist)
  }
  G
  
} 

