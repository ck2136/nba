#' Net Benefit Analysis of Prediction Models and Predictive Markers
#'
#' \code{nba} implements Vicker's et al's net benefit calculation
#' (based on weighing the true and false positive cases) for prediction models
#' and predictive markers. It can be used to compare the clinical utility
#' of prediction models.
#'
#' @param formula A formula containing the outcome and predictor variables.
#' @param data A data frame containing the variables in the model.
#' @param xstart The starting probability threshold value.
#' @param xstop The ending probability threshold value.
#' @param xby The interval between \code{xstart} and \code{xstop}.
#' @param ymin The minimum value of net benefit (for displaying/graphing).
#' @param harm A vector of length equal to the vector of predictors. Specify the value of harm into model assessment.
#' @param graph If set to \code{TRUE}, net benefit graphic of the prediction model will be outputted on the device graphics output.
#' @param intervention If set to \code{TRUE}, net reduction in interventions per 100 patients are plotted.
#' @param pred.model Prediction models to use. This should be a vector of models such as \code{"glm"} or \code{"svm"}.
#' @param interventionper Number used to calculate the number of interventions avoided by. Default is \code{100}.
#' @param smooth If set to \code{TRUE}, a smoothed net benefit values are calculated.
#' @param loess.span The value of \code{span} for the loess smoothing.
#' @return If the nba() function is stored in an object, the resulting output will contain an object of type list
#'   containing values of the net benefit of the specified prediction models in addition to the reduction in
#'   interventions.
#' @author Chong Kim \email{chong.kim@ucdenver.edu} and Andrew Vickers \email{vickersa@mskcc.org} based on the \href{http://journals.sagepub.com/doi/abs/10.1177/0272989x06295361}{2006 Article}
#' @references Vickers, A. (2008), Decision Curve Analysis: A Novel Method for Evaluating Prediction Models Vol 26,
#'   Issue 6, 2006.
#' @seealso \code{\link[DecisionCurve]{DecisionCurve}}
#' @examples
#' ## standard net benefit analysis
#' nba(cancer ~ famhistory + age + marker, data.set, pred.mode= c("glm","rpart","svm","rf"))
#'
#' ## loess smooth
#' modified_dca(cancer ~ famhistory + age + marker, data.set, pred.mode= c("glm","rpart","svm","rf"), smooth = TRUE, loess.span = 0.5)
#'
#' ## plot intervention reduced
#' modified_dca(cancer ~ famhistory + age + marker, data.set, pred.mode= c("glm","rpart","svm","rf"), smooth = TRUE, loess.span = 0.5, intervention = TRUE)




nba <- function(formula, data, xstart=0.01, xstop=0.99, xby=0.01,
                         ymin=-0.05, harm=NULL,graph=TRUE, intervention=FALSE, pred.model,
                         interventionper=100, smooth=FALSE,loess.span=0.10) {

  # LOADING REQUIRED LIBRARIES
  require(stats)
  require(xgboost)
  require(e1071)
  require(rpart)
  require(randomForest)
  require(dplyr)

  # data MUST BE A DATA FRAME
  if (class(data)!="data.frame") {
    stop("Input data must be class data.frame")
  }

  # formula MUST BE A FORMULA
  if (class(formula) != "formula") {
    stop("Formula is not a valid formula class!")
  }

  #ONLY KEEPING COMPLETE CASES
  data=data[complete.cases(data[, all.vars(formula)]), all.vars(formula)]

  # outcome MUST BE CODED AS 0 AND 1
  if (max(data[[all.vars(formula)[1]]])>1 | min(data[[all.vars(formula)[1]]])<0) {
    stop("outcome cannot be less than 0 or greater than 1")
  }

  # xstart IS BETWEEN 0 AND 1
  if (xstart<0 | xstart>1) {
    stop("xstart must lie between 0 and 1")
  }

  # xstop IS BETWEEN 0 AND 1
  if (xstop<0 | xstop>1) {
    stop("xstop must lie between 0 and 1")
  }

  # xby IS BETWEEN 0 AND 1
  if (xby<=0 | xby>=1) {
    stop("xby must lie between 0 and 1")
  }

  # xstart IS BEFORE xstop
  if (xstart>=xstop) {
    stop("xstop must be larger than xstart")
  }

  # STORING THE NUMBER OF PREDICTORS SPECIFIED.. here it is based on the number of algorithms
  pred.n=length(pred.model)


  # Based on the algorithm selected attach the predicted probabilities onto the original dataset
  if("glm" %in% pred.model) {
    fit <- glm(formula, data, family=binomial(link = "logit"))
    data$pred_glm <- fit$fitted.values
  }
  if("rpart" %in% pred.model) {
    fit <- rpart(formula, data)
    data$pred_rpart <- predict(fit, data)
  }
  if("rf" %in% pred.model) {
    fit <- randomForest(formula, data)
    data$pred_rf <- predict(fit, data)
  }
  if("svm" %in% pred.model) {
    fit <- svm(formula, data)
    data$pred_svm <- predict(fit, data)
  }
  if("gbm" %in% pred.model) {
    fit <- xgboost(formula, data)
    data$pred_gbm <- predict(fit, data)
  }

  #IF harm SPECIFIED ENSURING THAT EACH PREDICTOR HAS A SPECIFIED HARM
  if (length(harm)>0 & pred.n!=length(harm)) {
    stop("Number of harms specified must be the same as the number of predictors being checked.")
  }

  #INITIALIZING DEFAULT VALUES FOR PROBABILITES AND HARMS IF NOT SPECIFIED
  if (length(harm)==0) {
    harm=rep(0,pred.n)
  }


  #########  CALCULATING NET BENEFIT   #########
  N=dim(data)[1]
  event.rate=colMeans(data[all.vars(formula)[1]])



  # CREATING DATAFRAME THAT IS ONE LINE PER THRESHOLD PER all AND none STRATEGY
  nb=data.frame(seq(from=xstart, to=xstop, by=xby))
  names(nb)="threshold"
  interv=nb

  nb["all"]=event.rate - (1-event.rate)*nb$threshold/(1-nb$threshold)
  nb["none"]=0

  # CYCLING THROUGH EACH PREDICTION MODEL AND CALCULATING NET BENEFIT
  for(m in 1:pred.n){
    for(t in 1:length(nb$threshold)){
      # COUNTING TRUE POSITIVES AT EACH THRESHOLD
      tp=mean(data[data[[paste0("pred_",pred.model[m])]]>=nb$threshold[t],all.vars(formula)[1]])*sum(data[[paste0("pred_",pred.model[m])]]>=nb$threshold[t])
      # COUNTING FALSE POSITIVES AT EACH THRESHOLD
      fp=(1-mean(data[data[[paste0("pred_",pred.model[m])]]>=nb$threshold[t],all.vars(formula)[1]]))*sum(data[[paste0("pred_",pred.model[m])]]>=nb$threshold[t])
      #setting TP and FP to 0 if no observations meet threshold prob.
      if (sum(data[[paste0("pred_",pred.model[m])]]>=nb$threshold[t])==0) {
        tp=0
        fp=0
      }

      # CALCULATING NET BENEFIT
      nb[t,paste0("pred_",pred.model[m])]=tp/N - fp/N*(nb$threshold[t]/(1-nb$threshold[t])) - harm[m]
    }
    interv[paste0("pred_",pred.model[m])]=(nb[paste0("pred_",pred.model[m])] - nb["all"])*interventionper/(interv$threshold/(1-interv$threshold))
  }

  # CYCLING THROUGH EACH PREDICTOR AND SMOOTH NET BENEFIT AND INTERVENTIONS AVOIDED
  for(m in 1:pred.n) {
    if (smooth==TRUE){
      lws=loess(data.matrix(nb[!is.na(nb[[paste0("pred_",pred.model[m])]]),paste0("pred_",pred.model[m])]) ~ data.matrix(nb[!is.na(nb[[paste0("pred_",pred.model[m])]]),"threshold"]),span=loess.span)
      nb[!is.na(nb[[paste0("pred_",pred.model[m])]]),paste0("pred_",pred.model[m],"_sm")]=lws$fitted

      lws=loess(data.matrix(interv[!is.na(nb[[paste0("pred_",pred.model[m])]]),paste0("pred_",pred.model[m])]) ~ data.matrix(interv[!is.na(nb[[paste0("pred_",pred.model[m])]]),"threshold"]),span=loess.span)
      interv[!is.na(nb[[paste0("pred_",pred.model[m])]]),paste0("pred_",pred.model[m],"_sm")]=lws$fitted
    }
  }

  # PLOTTING GRAPH IF REQUESTED
  if (graph==TRUE) {
    require(graphics)

    # PLOTTING INTERVENTIONS AVOIDED IF REQUESTED
    if(intervention==TRUE) {
      # initialize the legend label, color, and width using the standard specs of the none and all lines
      legendlabel <- NULL
      legendcolor <- NULL
      legendwidth <- NULL
      legendpattern <- NULL

      #getting maximum number of avoided interventions
      ymax=max(interv[paste0("pred_",pred.model)],na.rm = TRUE)

      #INITIALIZING EMPTY PLOT WITH LABELS
      plot(x=nb$threshold, y=nb$all, type="n" ,xlim=c(xstart, xstop), ylim=c(ymin, ymax), xlab="Threshold probability", ylab=paste("Net reduction in interventions per",interventionper,"patients"))

      #PLOTTING INTERVENTIONS AVOIDED FOR EACH PREDICTOR
      for(m in 1:pred.n) {
        if (smooth==TRUE){
          lines(interv$threshold,data.matrix(interv[paste0("pred_",pred.model[m],"_sm")]),col=m,lty=2)
        } else {
          lines(interv$threshold,data.matrix(interv[paste0("pred_",pred.model[m])]),col=m,lty=2)
        }

        # adding each model to the legend
        legendlabel <- c(legendlabel, paste0("pred_",pred.model[m]))
        legendcolor <- c(legendcolor, m)
        legendwidth <- c(legendwidth, 1)
        legendpattern <- c(legendpattern, 2)
      }
    } else {
      # PLOTTING NET BENEFIT IF REQUESTED

      # initialize the legend label, color, and width using the standard specs of the none and all lines
      legendlabel <- c("None", "All")
      legendcolor <- c(17, 8)
      legendwidth <- c(2, 2)
      legendpattern <- c(1, 1)

      #getting maximum net benefit
      ymax=max(nb[names(nb)!="threshold"],na.rm = TRUE)

      # inializing new benfit plot with treat all option
      plot(x=nb$threshold, y=nb$all, type="l", col=8, lwd=2 ,xlim=c(xstart, xstop), ylim=c(ymin, ymax), xlab="Threshold probability", ylab="Net benefit")
      # adding treat none option
      lines(x=nb$threshold, y=nb$none,lwd=2)
      #PLOTTING net benefit FOR EACH PREDICTOR
      for(m in 1:pred.n) {
        if (smooth==TRUE){
          lines(nb$threshold,data.matrix(nb[paste0("pred_",pred.model[m],"_sm")]),col=m,lty=2)
        } else {
          lines(nb$threshold,data.matrix(nb[paste0("pred_",pred.model[m])]),col=m,lty=2)
        }
        # adding each model to the legend
        legendlabel <- c(legendlabel, paste0("pred_",pred.model[m]))
        legendcolor <- c(legendcolor, m)
        legendwidth <- c(legendwidth, 1)
        legendpattern <- c(legendpattern, 2)
      }
    }
    # then add the legend
    legend("topright", legendlabel, cex=0.8, col=legendcolor, lwd=legendwidth, lty=legendpattern)

  }

  #RETURNING RESULTS
  results=list()
  results$N=N
  results$pred.model=data.frame(cbind(pred.model,harm,probability))
  names(results$pred.model)=c("models","harm.applied","probability")
  results$interventions.avoided.per=interventionper
  results$net.benefit=nb
  results$interventions.avoided=interv

  return(results)

}
