library(spida2)
library(effects)
library(car)

context("LFX Deriv Tests")

test_that("Loaded Properly", {
  expect_equal(1,1)
})

test_that("Arrests Age", {

  ans <- list("1 * M(1/age)","0 * M(sex)","0 * M(colour)","1 * M(I(2 * age))","1 * M(I(-(1/age^2)))",
    "1 * M(I(1 + 1/age + exp(age)))","1 * M(1/age) * M(sex)","1 * M(1/age) * M(colour)",
    "0 * M(sex) * M(colour)", "1 * M(sex) * M(I(2 * age))", "1 * M(colour) * M(I(2 * age))",
    "1 * M(I(exp(age) * ((1/age + log(age)) * age + log(age)))) * numVar",
    "1 * M(1/age) * M(sex) * M(colour)", "1 * M(sex) * M(colour) * M(I(2 * age))")

  ds <- Arrests
  ds$numVar <- rnorm(dim(ds[1]))
  asc <- glm(released ~ (log(age) + sex + colour)^3 + I(age^2) + I(age^2)*sex*colour + I(1.0/age) +
               I(exp(age) + log(age) + age) + I(exp(age)*log(age)*age):numVar, data=ds, family=binomial)
  vName <- "age"
  dsForm <- LfxDerivative(asc, vName)
  expect_equal(dsForm, ans)
})

test_that("Arrests Sex", {

  ans <- list("0 * M(log(age))","1 * M(1)","0 * M(colour)","0 * M(I(age^2))","0 * M(I(1/age))","0 * M(I(exp(age) + log(age) + age))",
              "1 * M(log(age))","0 * M(log(age)) * M(colour)","1 * M(colour)",
              "1 * M(I(age^2))", "0 * M(colour) * M(I(age^2))", "0 * M(I(exp(age) * log(age) * age)) * numVar",
              "1 * M(log(age)) * M(colour)", "1 * M(colour) * M(I(age^2))")

  ds <- Arrests
  ds$numVar <- rnorm(dim(ds[1]))
  asc <- glm(released ~ (log(age) + sex + colour)^3 + I(age^2) + I(age^2)*sex*colour + I(1.0/age) +
               I(exp(age) + log(age) + age) + I(exp(age)*log(age)*age):numVar, data=ds, family=binomial)
  vName <- "sex"
  dsForm <- LfxDerivative(asc, vName)
  expect_equal(dsForm, ans)
})

test_that("Implicit Diff", {

  mod <- lm(log(income) ~  education * type * women * prestige, data = Prestige)

  ans <- list("1 * 1 * (1/income)^(-1)", "0 * M(type) * (1/income)^(-1)", "0 * women * (1/income)^(-1)",
              "0 * prestige * (1/income)^(-1)", "1 * M(type) * (1/income)^(-1)",
              "1 * women * (1/income)^(-1)", "0 * M(type) * women * (1/income)^(-1)",
              "1 * prestige * (1/income)^(-1)", "0 * M(type) * prestige * (1/income)^(-1)",
              "0 * women * prestige * (1/income)^(-1)", "1 * M(type) * women * (1/income)^(-1)", "1 * M(type) * prestige * (1/income)^(-1)",
              "1 * women * prestige * (1/income)^(-1)", "0 * M(type) * women * prestige * (1/income)^(-1)",
              "1 * M(type) * women * prestige * (1/income)^(-1)")
  vName <- "education"
  dsForm <- LfxDerivative(mod, vName, yName = "income")
  expect_equal(dsForm, ans)
})


test_that("Equivalent Wald", {

  server <- 'blackwell.math.yorku.ca'
  server <- '3.83.113.57'

  dall <- read.csv(paste0("http://",server,"/data/Smoking3.csv"))
  dd <- subset( dall, sex == 'BTSX')   # subset of a data frame (combined sexex)
  dd$LifeExp <- dd$lifeexp.Birth # Life expectancy at birth
  dd$LE <- dd$LifeExp
  dd$smoke <- dd$consumption.cigPC # cigarette consumption per adult per year
  dd$HE <- dd$HealthExpPC.Tot.ppp  # health expenditures per capita in US$ PPP
  dd$hiv <- dd$hiv_prev15_49  # prevalence of HIV in population 15 to 49
  dd$special <- ifelse(
    dd$country %in% c('Angola','Sierra Leone','Equatorial Guinea'),
    1,
    0)  # indicator variable for 3 outlying countries

  fit.hiv2 <- lm( LifeExp ~ log(HE) * (smoke + I(smoke^2)) + hiv+special, dd ,
                  na.action = na.exclude)

  pred <- expand.grid(
    HE = c(50,150,500, 1000, 1500, 5000),
    smoke = seq(10,2000,20),
    hiv = 0,
    special = 0)

  vName <- "smoke"
  autoL <- LfxDerivData(fit.hiv2,vName,pred)
  wwA <- wald(fit.hiv2,autoL)
  wwA <- as.data.frame(wwA)

  L <- Lfx(fit.hiv2,
           list( 0,
                 0 * M(log(HE)),
                 1 ,
                 1 * M(I(2*smoke)),
                 0 * hiv,
                 0 * special,
                 1 * M(log(HE)) * 1,
                 1 * M(log(HE)) * M(I(2*smoke))
           ), pred)

  wwO <- wald(fit.hiv2, L)
  wwO <- as.data.frame(wwO)

  expect_equal(sum(wwA == wwO) == sum(dim(wwO)[1] * dim(wwO)[2]),TRUE)
})

test_that("Three Factor",{
  mod <- lm(log(income) ~  education*type, data = Prestige)
  vName <- "education"
  autoL <- LfxDerivData(mod, vName, data = getData(mod))
  wwA <- wald(mod,autoL)
  wwA <- as.data.frame(wwA)

  L <- Lfx(mod, list(0, 1, 0 * M(type), 1 * M(type)))

  wwO <- wald(mod, L)
  wwO <- as.data.frame(wwO)

  expect_equal(sum(wwA == wwO) == sum(dim(wwO)[1] * dim(wwO)[2]),TRUE)
})
