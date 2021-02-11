library(FactoMineR)

# screening for variables with 0 variance
varZero <- function(dat) {
  out <- apply(dat,2, function(x) length(unique(x)))
  want <- which(!out > 1)
  as.numeric(unlist(want))
}


# Elastic net
en <- function(dat, method = "cv", kfold = 5, metric = "ROC"){
  set.seed(202101)
  
  nr_x <- ncol(dat)-2
  names <- colnames(dat)[1:nr_x]
  coeff <- setNames(rep(0, nr_x), names) # set initial coefficients to 0
  
  # drop columns with 0 variance to avoid error
  drop <- varZero(dat)
  if(length(drop) > 0) {
    dat <- dat[,-drop]
  }
  
  # Standardize including covariates
  dat[,-ncol(dat)] %<>% scale(scale = T)

  # update number of X variables
  nr_x <- ncol(dat)-2
  names <- colnames(dat)[1:nr_x]
  
  # 0 indicates no penalization. Age is not penalized
  penalty <- c(rep(1,nr_x),0)
  
  # Grid search for hyperparameters
  glmnet_grid <- expand.grid(alpha = (1:5)*0.2,
                             lambda = seq(.0001, 1, length = 200))
  glmnet_ctrl <- trainControl(method = method, number = kfold, classProbs = T, 
                              summaryFunction=twoClassSummary)
  
  glmnet_fit <- suppressWarnings(train(x = dat[,-ncol(dat)],
                                       y = factor(dat[,ncol(dat)], levels = c(1,2), labels = c("yes","no")),
                                       method = "glmnet",
                                       penalty.factor = penalty, #--> This is the penalty
                                       family = "binomial",
                                       preProcess = c("center", "scale"),
                                       metric = metric,
                                       tuneGrid = glmnet_grid,
                                       trControl = glmnet_ctrl))

  # Fit a model
  summary(glmnet_fit)
  glmnet_fit$bestTune
  
  fit<- glmnet(x = dat[,-ncol(dat)],
                      y = factor(dat[,ncol(dat)]),
                      family = "binomial", 
                      lambda = glmnet_fit$bestTune$lambda, 
                      alpha = glmnet_fit$bestTune$alpha)
  
  nm <- c(names,'age')
  
  sel_names <- nm[which((coef(fit) %>% as.vector)[-1] !=0)] # -1 is the intercept
  sel_coef <- (coef(fit) %>% as.vector)[(coef(fit) %>% as.vector)!=0][-1]
  
  coeff[names(coeff) %in% sel_names] <- sel_coef
  
  return(coeff)
}

# Random forest. The main hyperparameter "mtry" is searched
rf <- function(dat, method = "oob", search="grid"){
  set.seed(202101)
  
  # Search for "mtry"
  rf_ctrl <- trainControl(method=method,search=search)
  tunegrid <- expand.grid(.mtry=c(1:30)) # mtry cannot be larger than number of independent variables
  rf_random <- train(x = dat[,-ncol(dat)],
                     y = factor(dat[,ncol(dat)], levels = c(1,2), labels = c("yes","no")),
                     method="rf", tuneGrid=tunegrid, trControl=rf_ctrl,
                     ntree = 500, replace = T)
  
  # matrix to dataframe to use ranger
  dat <- data.frame(dat)
  dat$y %<>% as.factor
  
  fit <- ranger::ranger(formula = y~., data = dat,
                 importance='impurity', 
                 mtry = as.numeric(rf_random$bestTune),
                 always.split.variables = 'age')
  
  rank_rf <- -fit$variable.importance %>% rank(ties.method = "min")
  
  return(rank_rf)
}

# Single point, firth is by default used to avoid saturation, last two column of dat should be covariate and outcome
sp <- function(dat, firth=F){
  set.seed(202101)
  
  nr_x <- ncol(dat)-2
  names <- colnames(dat)[1:nr_x]
  p_val <- setNames(rep(1, nr_x), names) # set initial p-value to 1

  # drop columns with 0 variance to avoid error
  drop <- varZero(dat)
  if(length(drop) > 0) {
    dat <- dat[,-drop]
  }
  
  # update number of X variables
  nr_x <- ncol(dat)-2
  names <- colnames(dat)[1:nr_x]

  # Firth to prevent saturation
  if(firth){
    for(i in 1:nr_x) {
      fit <- tryCatch({logistf::logistf(as.numeric(as.vector(factor(dat[,ncol(dat)], levels = c(1,2), labels = c('1','0'))))~ 
                                          dat[ ,(ncol(dat)-1)] + dat[,i], firth=T, family = "binomial")},
                      error = function(e) print(dat[,i]))
      p_val[names(p_val) == names[i]] <- fit$prob[3] # save p
    }
  }else{
    dat %<>% data.frame()
    for(i in 1:nr_x) {
      fit <-suppressWarnings({glm(factor(y, levels = c(1,2), labels = c('yes','no'))~age+dat[,i], data=dat,
                          family = "binomial")})
      p_val[names(p_val) == names[i]] <- summary(fit)$coef[3,4] # 3rd row (1.intercept, 2.age, 3.bacteria); 4th column p-value
    }
  }

  return(p_val)
}

# Correlation heatmap plotting
myHeatmap <- function(data, type='spearman', name) {
  library(reshape2)
  library(ggcorrplot)
  cormat <- round(cor(data,method=type),2)

  cormat <- cormat[order(colSums(cormat)),order(colSums(cormat))]
  med <- as.character(median(as.vector(cormat)))
  iqr <- as.character(IQR(as.vector(cormat)))
  
  p <- ggcorrplot(cormat,hc.order=FALSE,type='upper',ggtheme=theme_void(),tl.cex = 5) + theme(axis.text = element_blank()) +
    annotate(geom="text",label=paste0("Median: ",med),x=ncol(data)*0.7,y=ncol(data)*0.3, size = 5) +
    annotate(geom="text",label=paste0("IQR: ",iqr),x=ncol(data)*0.7,y=ncol(data)*0.22, size = 5)
  return(p)
}

# Correlation heatmap plotting when p is small (position of notation adjusted)
myHeatmap_sub <- function(data, type='spearman', name) {
  library(reshape2)
  library(ggcorrplot)
  cormat <- cor(data,method=type)
  #  cormat <- round(cor(data,method="pearson"),2)

  # cormat <- cormat[order(colSums(cormat)),order(colSums(cormat))]
  med <- as.character(round(median(as.vector(cormat)), digits = 2))
  iqr <- as.character(round(IQR(as.vector(cormat)), digits = 2))
  
  p <- ggcorrplot(cormat,hc.order=FALSE,type='upper', lab = T)
  return(p)
}

# PCA plots with specified outliers highlighted
myPCA <- function(data,groups,name,nPC=3,robust=FALSE) {
  # Please change outliers list here
  # outlls<- c('S00140','S00114','S00019','S00024','S00110','S00183','S00055','S00061')
  outlls <- NA
  if (robust==TRUE) {
    library(pcaPP)
    fit <- pcaPP::PCAgrid(data,k=3)
    loadings <- fit$loadings
    var <- NULL
    plot_data <- data.frame("Condition"=as.factor(groups),fit$scores)
    names(plot_data)[2:4] <- c("Dim.1","Dim.2","Dim.3")
    
  } else {
    fit <- FactoMineR::PCA(data,scale.unit = TRUE, graph = F)
    loadings <- loadings<-sweep(fit$var$coord,2,sqrt(fit$eig[1:5,1]),FUN="/")
    var <- round(fit$eig[,2],2)
    plot_data <- data.frame("Condition"=as.factor(groups),fit$ind$coord[,1:nPC])
    plot_data %<>% mutate(sample = rownames(plot_data))
  }
  
  p1 <- ggplot(data=plot_data, aes(x=Dim.1,y=Dim.2,color=Condition,label=sample)) +
    geom_point(aes(shape=Condition)) +
    geom_text(data=subset(plot_data, sample %in% outlls), aes(label=sample)) +
    stat_ellipse() +
    xlab(paste0("PC1 (",var[1],"%)")) +
    ylab(paste0("PC2 (",var[2],"%)")) + 
    theme(axis.text.y = element_text(size = 15),
          axis.title = element_text(size = 15),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 15)) +
    guides(color=guide_legend("CVID"),
           shape=F)

  # p2 <- ggplot(data=plot_data, aes(x=Dim.1,y=Dim.3,color=Condition)) +
  #   geom_point(aes(shape=Condition)) +
  #   stat_ellipse() +
  #   xlab(paste0("PC1 - ",var[1],"% Explained Variance")) +
  #   ylab(paste0("PC3 - ",var[3],"% Explained Variance"))
  
  p3 <- ggplot(data=plot_data, aes(x=Dim.2,y=Dim.3,color=Condition,label=sample)) +
    geom_point(aes(shape=Condition)) +
    geom_text(data=subset(plot_data, sample %in% outlls), aes(label=sample)) +
    stat_ellipse() +
    xlab(paste0("PC2 - ",var[2],"% Explained Variance")) +
    ylab(paste0("PC3 - ",var[3],"% Explained Variance")) + 
    ggtitle(paste("Score plot of", name, "data"))
  
  
  out <- list("PCA"=fit,"Loadings"=loadings,"P1"=p1, "P3"=p3, "plot_data"=plot_data)
  return(out)
}