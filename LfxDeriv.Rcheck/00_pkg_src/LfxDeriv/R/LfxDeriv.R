# NOTE: Modified Deriv_ starting line 520 to wrap the M and I in calls rather than do the chain rule or ignore them
# NOTE: changed nd2expr line 1569 to NOT sort alphabetically, this way only the variable that is derived with respect to is moved to the front or the back only
# NOTE: changed nd2expr line 1569 to NOT sort by function call, this way only the variable that is derived with respect to is moved to the front or the back only
# NOTE: Commented out the sub expression sub in in Deriv at line 340 ish, not sure this is useful, just looks neater to have {.e1 <- log(age); .e1*exp(age)}
# Just makes it difficult to parse

#' reArrance
#'
#' @param dsSp Deriviates
#' @param oriDexs Original index of the variable we are deriving with respect to
#' @param newDexs New index of the variable we are deriving with respect to
#'
#' @return Deriviate with all variables back in their original positions
#'
#'
#' @examples
#' varList <- ret
#' oriDexs <- lapply(varList, function(x) grep(vName,unlist(x)))
#' varList <- lapply(varList, paste, collapse = " * ")
#' varArr <- array(varList,c(1,length(varList)))
#' ds <- lapply(varArr,function(x) Deriv(x,vName))
#' dsSp <- lapply(ds, function (x) gsub("\\) \\* M","\\)  \\*  M", x))
#' dsSp <- lapply(dsSp, function (x) unlist(strsplit(x, " \\* |(?>\\(.*?\\).*?\\K(, |$))", perl = TRUE)))
#' dsSp <- lapply(dsSp, function(x) unlist(strsplit(x, "  \\*  ")))
#' newDexs <- lapply(dsSp, function(x) grep(vName,x))
#' dsFormatted <- reArrange(dsSp,oriDexs,newDexs)
#'
reArrange <- function(dsSp, oriDexs = oriDexs, newDexs = newDexs)
{
  # Need to shuffle
  # Output always seems first or last
  # Take the value beyond the new spot for the variable up to the end of the previous spot (-1 as the oriDexs have the 1 * as the first element of the list)
  # If the derivate doesn't involve the variable we should still keep the deriviate results
  shuf <- sapply(seq_along(1:length(dsSp)),function(i) {
    if((length(newDexs[[i]]) > 0))
    {
      if(newDexs[[i]] != (oriDexs[[i]]))
      {
        if(newDexs[[i]] > oriDexs[[i]])
        {
          c(dsSp[[i]][newDexs[[i]]],dsSp[[i]][1:max(1,(oriDexs[[i]]-1))]) #dsSp[[i]][oriDexs[[i]]:max(oriDexs[[i]],(length(dsSp[[i]])-1))])
        }
        else
        {
          c(dsSp[[i]][(newDexs[[i]]+1):max(oriDexs[[i]],newDexs[[i]]+1)],dsSp[[i]][newDexs[[i]]]) #,dsSp[[i]][oriDexs[[i]]:max(length(dsSp[[i]]),oriDexs[[i]])])
        }
      }
      else
      {
        dsSp[[i]] <- dsSp[[i]]
      }
    }
    else
    {
      dsSp[[i]] <- dsSp[[i]]
    }
  } )

  rmNA <- lapply(shuf, function(x) x[!is.na(x)])

  strV <- lapply(rmNA, paste, collapse=" * ")
  res <- lapply(seq_along(1:length(strV)), function(i) {
    if(strV[[i]] == "0")
    {
      strV[[i]] = strV[[i]]
    }
    else
    {
      paste("1 * ", strV[[i]], collapse="")
    }
  }
  )

  res
}

#' LfxDerivative
#' @title Creates an L matrix for wald testing via symbollic differentiation of a lm, glm, or lme R model
#' @description Creates an L matrix for wald testing via symbollic differentiation of a lm, glm, or lme R model
#'
#' @param fit model
#' @param vName variable to take the derivative with respect to
#' @param yName the LHS of the model formula if you wish to perform implicit differentiation
#' @param factors factors involved in model, will determine this via getFactoNames(fit) by default
#'
#' @return a matrix of strings, the first column is the original formula broken up into
#' the additive poritions of the model, the second is the derivative of each additive porition of the model
#' this can be used with evalExpr <- do.call(cbind,lapply(dsFormatted, function (x) eval(parse(text=x),data))) to evaluate the derivative across a data
#' frame
#'
#' @examples
LfxDerivative <- function(fit, vName, yName = NULL, factors = getFactorNames(fit))
{
  ts <- colnames(attr(terms(fit), "factors"))
  ts <- strsplit(ts, ":")
  factors <- getFactorNames(fit)
  facl <- lapply(ts, function(x) x %in% factors | grepl(")",
                                                        x))
  ret <- lapply(seq_along(facl), function(i) {
    lapply(seq_along(facl[[i]]), function(j) if (facl[[i]][j])
      paste("M(", ts[[i]][j], ")", sep = "")
      else ts[[i]][j])
  })

  # Implicit Differentation Required
  if(!is.null(yName))
  {
    lhs <- colnames(getData(fit))[1]
    lhsD <- Deriv(lhs,yName)
    lhsD <- paste("(",lhsD,")^(-1)",sep="")
  }

  varList <- ret
  oriDexs <- lapply(varList, function(x) grep(vName,unlist(x)))
  varList <- lapply(varList, paste, collapse = " * ")
  varArr <- array(varList,c(1,length(varList)))

  ds <- lapply(varArr,function(x) Deriv(x,vName))
  # https://stackoverflow.com/questions/35347537/using-strsplit-in-r-ignoring-anything-in-parentheses/35347645
  dsSp <- lapply(ds, function (x) gsub("\\) \\* M","\\)  \\*  M", x))
  dsSp <- lapply(dsSp, function (x) unlist(strsplit(x, " \\* |(?>\\(.*?\\).*?\\K(, |$))", perl = TRUE)))
  dsSp <- lapply(dsSp, function(x) unlist(strsplit(x, "  \\*  ")))
  newDexs <- lapply(dsSp, function(x) grep(vName,x))

  dsFormatted <- reArrange(dsSp,oriDexs,newDexs)

  # Implicit Differentiation so multiple by inverse of LHS
  if(!is.null(yName))
  {
    dsFormatted <- lapply(dsFormatted,function (x) if( x != "0") paste(x," * ", lhsD, sep="") else x)
  }

  cbind(c(varList),c(dsFormatted))
}

######################################################################
# Compare to actual usage case from notes
######################################################################

LfxDerivData <- function(fit, vName, data, yName = NULL, factors = getFactorNames(fit), wrap = FALSE, debug = debug)
{
  ts <- colnames(attr(terms(fit), "factors"))
  ts <- strsplit(ts, ":")
  factors <- getFactorNames(fit)
  facl <- lapply(ts, function(x) x %in% factors | grepl(")",
                                                        x))
  ret <- lapply(seq_along(facl), function(i) {
    lapply(seq_along(facl[[i]]), function(j) if (facl[[i]][j])
      paste("M(", ts[[i]][j], ")", sep = "")
      else ts[[i]][j])
  })

  # Implicit Differentation Required
  if(!is.null(yName))
  {
    lhs <- colnames(getData(fit))[1]
    lhsD <- Deriv(lhs,yName)
    lhsD <- paste("(",lhsD,")^(-1)",sep="")
  }

  varList <- ret
  oriDexs <- lapply(varList, function(x) grep(vName,unlist(x)))
  varList <- lapply(varList, paste, collapse = " * ")
  varArr <- array(varList,c(1,length(varList)))

  ds <- lapply(varArr,function(x) Deriv(x,vName))
  # https://stackoverflow.com/questions/35347537/using-strsplit-in-r-ignoring-anything-in-parentheses/35347645
  dsSp <- lapply(ds, function (x) gsub("\\) \\* M","\\)  \\*  M", x))
  dsSp <- lapply(dsSp, function (x) unlist(strsplit(x, " \\* |(?>\\(.*?\\).*?\\K(, |$))", perl = TRUE)))
  dsSp <- lapply(dsSp, function(x) unlist(strsplit(x, "  \\*  ")))
  newDexs <- lapply(dsSp, function(x) grep(vName,x))

  dsFormatted <- reArrange(dsSp,oriDexs,newDexs)

  # Implicit Differentiation so multiple by inverse of LHS
  if(!is.null(yName))
  {
    dsFormatted <- lapply(dsFormatted,function (x) if( x != "0") paste(x," * ", lhsD, sep="") else x)
  }

  evalExpr <- do.call(cbind,lapply(dsFormatted, function (x) eval(parse(text=x),data)))
  # Add 0 column for the intercept
  L <- cbind(0,evalExpr)
  gg <- getFix(fit)
  rownames(L) <- rownames(data)
  colnames(L) <- names(gg$fixed)
  attr(L, "data") <- data
  L
}

# library(spida2)
# library(effects)
# ds <- Arrests
# ds$numVar <- rnorm(dim(ds[1]))
# #asc <- glm(released ~ I(2*exp(age)*log(age)*age), data=ds, family=binomial)
# #asc <- glm(released ~ I(exp(age) + log(age) + age) + I(2*exp(age)*log(age)):numVar + I(2*exp(age)*log(age)*age):numVar, data=ds, family=binomial)
# asc <- glm(released ~ (log(age) + sex + colour)^3 + I(age^2) + I(age^2)*sex*colour + I(1.0/age) + I(exp(age) + log(age) + age) + I(exp(age)*log(age)*age):numVar, data=ds, family=binomial)
# summary(asc)
#
# vName <- "age"
# dsForm <- LfxDerivative(asc, vName)
# dsForm
#
# vName <- "sex"
# dsForm <- LfxDerivative(asc, vName)
# dsForm
#
# library(car)
# mod <- lm(log(income) ~  education * type * women * prestige, data = Prestige)
#
# vName <- "education"
# dsForm <- LfxDerivate(mod, vName, yName = "income")
# dsForm
#
# server <- 'blackwell.math.yorku.ca'
# server <- '3.83.113.57'
#
# dall <- read.csv(paste0("http://",server,"/data/Smoking3.csv"))
# dd <- subset( dall, sex == 'BTSX')   # subset of a data frame (combined sexex)
# dd$LifeExp <- dd$lifeexp.Birth # Life expectancy at birth
# dd$LE <- dd$LifeExp
# dd$smoke <- dd$consumption.cigPC # cigarette consumption per adult per year
# dd$HE <- dd$HealthExpPC.Tot.ppp  # health expenditures per capita in US$ PPP
# dd$hiv <- dd$hiv_prev15_49  # prevalence of HIV in population 15 to 49
# dd$special <- ifelse(
#   dd$country %in% c('Angola','Sierra Leone','Equatorial Guinea'),
#   1,
#   0)  # indicator variable for 3 outlying countries
# head(dd)
#
# fit.hiv2 <- lm( LifeExp ~ log(HE) * (smoke + I(smoke^2)) + hiv+special, dd ,
#                 na.action = na.exclude)
#
# pred <- expand.grid(
#   HE = c(50,150,500, 1000, 1500, 5000),
#   smoke = seq(10,2000,20),
#   hiv = 0,
#   special = 0)
#
# vName <- "smoke"
# autoL <- LfxDerivData(fit.hiv2,vName,pred)
# wwA <- wald(fit.hiv2,autoL)
# wwA <- as.data.frame(wwA)
# head(wwA)
#
# L <- Lfx(fit.hiv2,
#          list( 0,
#                0 * M(log(HE)),
#                1 ,
#                1 * M(I(2*smoke)),
#                0 * hiv,
#                0 * special,
#                1 * M(log(HE)) * 1,
#                1 * M(log(HE)) * M(I(2*smoke))
#          ), pred)
#
# wwO <- wald(fit.hiv2, L)
# wwO <- as.data.frame(wwO)
# head(wwO)
#
# sum(wwA == wwO) == sum(dim(wwO)[1] * dim(wwO)[2])
#
