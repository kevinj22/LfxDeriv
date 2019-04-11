# NOTE: Modified Deriv_ starting line 520 to wrap the M and I in calls rather than do the chain rule or ignore them
# NOTE: changed nd2expr line 1569 to NOT sort alphabetically, this way only the variable that is derived with respect to is moved to the front or the back only
# NOTE: changed nd2expr line 1569 to NOT sort by function call, this way only the variable that is derived with respect to is moved to the front or the back only
# NOTE: Commented out the sub expression sub in in Deriv at line 340 ish, not sure this is useful, just looks neater to have {.e1 <- log(age); .e1*exp(age)}
# Just makes it difficult to parse

#' reArrange
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
reArrange <- function(dsSp, oriVars, oriDexs = oriDexs, newDexs = newDexs)
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
      # Need the M(categorical) to create all proper columns for the L matrices
      strV[[i]] = paste(strV[[i]], "*", as.character(oriVars[[i]]))
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
#' @return a matrix of strings that is the derivative of each additive porition of the model
#' this can be used with evalExpr <- do.call(cbind,lapply(dsFormatted, function (x) eval(parse(text=x),data))) to evaluate the derivative across a data
#' frame. For an end to end function see LfxDerivData.
#'
#' @examples
#' ds <- Arrests
#' ds$numVar <- rnorm(dim(ds[1]))
#' asc <- glm(released ~ (log(age) + sex + colour)^3 + I(age^2) + I(age^2)*sex*colour + I(1.0/age) +
#'             I(exp(age) + log(age) + age) + I(exp(age)*log(age)*age):numVar, data=ds, family=binomial)
#' vName <- "sex"
#' dsForm <- LfxDerivative(asc, vName)
#'
#'
#' Implicit Differentiation:
#' mod <- lm(log(income) ~  education * type * women * prestige, data = Prestige)
#' vName <- "education"
#' dsForm <- LfxDerivative(mod, vName, yName = "income")

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

  dsFormatted <- reArrange(dsSp,varArr,oriDexs,newDexs)
  dsFormatted <- lapply(dsFormatted, function (x) gsub(" {2,}"," ",x))

  # Implicit Differentiation so multiple by inverse of LHS
  if(!is.null(yName))
  {
    dsFormatted <- lapply(dsFormatted,function (x) if( x != "0") paste(x," * ", lhsD, sep="") else x)
  }

  #cbind(c(varList),c(dsFormatted))
  dsFormatted
}

#' LfxDerivData
#' @title Creates an L matrix for wald testing via symbollic differentiation of a lm, glm, or lme R model on the given data.
#' @description Creates an L matrix for wald testing via symbollic differentiation of a lm, glm, or lme R model
#'
#' @param fit model
#' @param vName variable to take the derivative with respect to
#' @param data data set used to fit the model by default, or one you wish to make predictions on
#' @param yName the LHS of the model formula if you wish to perform implicit differentiation
#' @param factors factors involved in model, will determine this via getFactoNames(fit) by default
#'
#' @return L hypothesis matrix
#'
#' @examples
#' mod <- lm(log(income) ~  education*type, data = Prestige)
#' vName <- "education"
#' autoL <- LfxDerivData(mod, vName, data = getData(mod))
#' wwA <- wald(mod,autoL)
LfxDerivData <- function(fit, vName, data = getData(fit), yName = NULL, factors = getFactorNames(fit), wrap = FALSE, debug = debug)
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

  dsFormatted <- reArrange(dsSp,varArr, oriDexs,newDexs)
  dsFormatted <- lapply(dsFormatted, function (x) gsub(" {2,}"," ",x))

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
