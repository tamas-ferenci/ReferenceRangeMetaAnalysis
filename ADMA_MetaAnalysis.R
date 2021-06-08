library(metafor)

###

eta <- function(Q) sapply(Q, function(Q) if(Q <= 115)
  2*prod((3*Q+1):(4*Q+1))/prod(1:Q)*integrate(function(z, Q) z*pnorm(z)^(3*Q)*(1-pnorm(z))^(Q)*dnorm(z),
                                              -Inf, Inf, Q = Q)$value else 1.35)

plotstudy <- function(y, m, ri.low, ri.up) {
  lines(c(ri.low[1], ri.low[3]), c(y, y), lwd = 10, col = "gray48")
  lines(c(ri.up[1], ri.up[3]), c(y, y), lwd = 10, col = "gray48")
  lines(c(ri.low[2], m[1]), c(y, y))
  lines(c(m[3], ri.up[2]), c(y, y))
  polygon(c(m, m[2]), c(y, y-0.3, y, y+0.3))
  points(ri.low[2], y, pch = "(") 
  points(ri.up[2], y, pch = ")") 
}

plotstudy2 <- function(y, data) {
  m <- data$m
  ri.low <- data$ri.low
  ri.up <- data$ri.up
  lines(c(ri.low[1], ri.low[3]), c(y, y), lwd = 10, col = "gray48")
  lines(c(ri.up[1], ri.up[3]), c(y, y), lwd = 10, col = "gray48")
  lines(c(ri.low[2], m[1]), c(y, y))
  lines(c(m[3], ri.up[2]), c(y, y))
  polygon(c( m, m[2]), c(y, y-0.3, y, y+0.3))
  points(ri.low[2], y, pch = "(") 
  points(ri.up[2], y, pch = ")")
}

plotsummary <- function(fit, fit.low, fit.up) {
  return(data.frame(m = c(max(0, fit$ci.lb), fit$b, fit$ci.ub),
                    ri.low = c(max(0, fit.low$ci.lb), fit.low$b, fit.low$ci.ub),
                    ri.up = c(fit.up$ci.lb, fit.up$b, fit.up$ci.ub)))
}

###

RawData <- rbind(rio::import("full_data_nemeth_ADMA_revision_NB_proc.xlsx", sheet = "HPLC",
                             .name_repair = "universal"),
                 rio::import("full_data_nemeth_ADMA_revision_NB_proc.xlsx", sheet = "ELISA",
                             .name_repair = "universal"))

RawData$no.weight <- RawData$Male+RawData$Female
RawData$no.weight[61] <- 28
RawData$MaleRatio <- RawData$Male/RawData$no.weight

RawData$Q1 <- as.numeric(gsub(",", ".", sapply(strsplit(RawData$ADMA.interquartile, "-"), `[`, 1)))
RawData$Q3 <- as.numeric(gsub(",", ".", sapply(strsplit(RawData$ADMA.interquartile, "-"), `[`, 2)))

RawData$Min <- as.numeric(gsub(",", ".", sapply(strsplit(RawData$ADMA.minmax, "-"), `[`, 1)))
RawData$Max <- as.numeric(gsub(",", ".", sapply(strsplit(RawData$ADMA.minmax, "-"), `[`, 2)))

RawData$ADMA.SD <- ifelse(is.na(RawData$ADMA.SD), RawData$ADMA.SEM*sqrt(RawData$no.weight), RawData$ADMA.SD)

# DOI: 10.1186/1471-2288-14-135, Eq. 14, 15
RawData$mean.ADMA.umol.l <- ifelse(!is.na(RawData$ADMA.median)&!is.na(RawData$Q1)&!is.na(RawData$Q3)&is.na(RawData$mean.ADMA.umol.l),
                                   (RawData$ADMA.median+RawData$Q1+RawData$Q3)/3, RawData$mean.ADMA.umol.l)
RawData$ADMA.SD <- ifelse(!is.na(RawData$ADMA.median)&!is.na(RawData$Q1)&!is.na(RawData$Q3)&is.na(RawData$ADMA.SD),
                          (RawData$Q3-RawData$Q1)/eta(round((RawData$no.weight-1)/4)), RawData$ADMA.SD)

# DOI: 10.1186/1471-2288-14-135, Eq. 3, 9
RawData$mean.ADMA.umol.l <- ifelse(!is.na(RawData$ADMA.median)&!is.na(RawData$Min)&!is.na(RawData$Max)&is.na(RawData$mean.ADMA.umol.l),
                                   (2*RawData$ADMA.median+RawData$Min+RawData$Max)/4, RawData$mean.ADMA.umol.l)
RawData$ADMA.SD <- ifelse(!is.na(RawData$ADMA.median)&!is.na(RawData$Min)&!is.na(RawData$Max)&is.na(RawData$ADMA.SD),
                          (RawData$Max-RawData$Min)/(2*qnorm((RawData$no.weight-0.375)/(RawData$no.weight+0.25))), RawData$ADMA.SD)

RawData$AgeQ1 <- as.numeric(gsub(",", ".", sapply(strsplit(RawData$AGE.IQR, "-"), `[`, 1)))
RawData$AgeQ3 <- as.numeric(gsub(",", ".", sapply(strsplit(RawData$AGE.IQR, "-"), `[`, 2)))
RawData$AgeMin <- as.numeric(gsub(",", ".", sapply(strsplit(RawData$AGE.minmax, "-"), `[`, 1)))
RawData$AgeMax <- as.numeric(gsub(",", ".", sapply(strsplit(RawData$AGE.minmax, "-"), `[`, 2)))

RawData$AGE.SD <- ifelse(is.na(RawData$AGE.SEM), RawData$AGE.SEM*sqrt(RawData$no.weight), RawData$AGE.SD)

# DOI: 10.1186/1471-2288-14-135, Eq. 14, 15
RawData$AGE <- ifelse(!is.na(RawData$AGE.median)&!is.na(RawData$AgeQ1)&!is.na(RawData$AgeQ3)&is.na(RawData$AGE),
                      (RawData$AGE.median+RawData$AgeQ1+RawData$AgeQ3)/3, RawData$AGE)
RawData$AGE.SD <- ifelse(!is.na(RawData$AGE.median)&!is.na(RawData$AgeQ1)&!is.na(RawData$AgeQ3)&is.na(RawData$AGE.SD),
                         (RawData$AgeQ3-RawData$AgeQ1)/eta(round((RawData$no.weight-1)/4)), RawData$AGE.SD)

# DOI: 10.1186/1471-2288-14-135, Eq. 3, 9
RawData$AGE <- ifelse(!is.na(RawData$AGE.median)&!is.na(RawData$AgeMin)&!is.na(RawData$AgeMax)&is.na(RawData$AGE),
                      (2*RawData$AGE.median+RawData$AgeMin+RawData$AgeMax)/4, RawData$AGE)
RawData$AGE.SD <- ifelse(!is.na(RawData$AGE.median)&!is.na(RawData$AgeMin)&!is.na(RawData$AgeMax)&is.na(RawData$AGE.SD),
                         (RawData$AgeMax-RawData$AgeMin)/(2*qnorm((RawData$no.weight-0.375)/(RawData$no.weight+0.25))), RawData$AGE.SD)

RawData <- RawData[!is.na(RawData$mean.ADMA.umol.l)&!is.na(RawData$ADMA.SD), ]

nrow(RawData)
length(unique(interaction(RawData$Author, RawData$DATE)))
table(RawData$fotipus[!duplicated(interaction(RawData$Author, RawData$DATE))])
sum(RawData$no.weight)
sum(RawData$Male, na.rm = TRUE)
sum(RawData$Female, na.rm = TRUE)
mean(RawData$AGE, na.rm = TRUE)
sd(RawData$AGE, na.rm = TRUE)
tapply(RawData$no.weight, RawData$fotipus, sum)
table(RawData$tipus)
tapply(RawData$no.weight, RawData$tipus, sum)
tapply(RawData$no.weight, RawData$fotipus, sum)

quantile((RawData$mean.ADMA.umol.l/RawData$ADMA.median-1)*100, na.rm = TRUE)
quantile((RawData$ADMA.median/((RawData$Q1+RawData$Q3)/2)-1)*100, na.rm = TRUE)
plot((RawData$Q1+RawData$Q3)/2,RawData$ADMA.median)
min(RawData$no.weight)

fit <- rma.uni(mi = mean.ADMA.umol.l, sdi = ADMA.SD, ni = no.weight, data = RawData, measure = "MN",
               slab = paste0(Author, " (", DATE, ")"), method = "REML")
fit.hplc <- update(fit, subset = fotipus=="HPLC")
fit.hplc.fluor <- update(fit.hplc, subset = tipus=="Fluorescence")
fit.hplc.lcms <- update(fit.hplc, subset = tipus=="LC/MS")
fit.elisa <- update(fit, subset = fotipus=="ELISA")
fit.elisa.dld <- update(fit.elisa, subset = tipus=="DLD")
fit.elisa.id <- update(fit.elisa, subset = tipus=="Immundiagnostik")

summary(fit)

ri_low <- pmax(0, RawData$mean.ADMA.umol.l-1.96 * RawData$ADMA.SD)
ri_up <- RawData$mean.ADMA.umol.l + 1.96*RawData$ADMA.SD
se <- RawData$ADMA.SD/sqrt(RawData$no.weight)

ki_low <- RawData$mean.ADMA.umol.l - 1.96 * RawData$ADMA.SD/sqrt(RawData$no.weight)
ki_up <- RawData$mean.ADMA.umol.l + 1.96 * RawData$ADMA.SD/sqrt(RawData$no.weight)

varchidist <- function(k) return(sapply(k, function(i) if(i<=300) i-2*(gamma((i+1)/2)/gamma(i/2))^2 else 0.5))

se_ri <- RawData$ADMA.SD*sqrt(1/RawData$no.weight+1.96^2*varchidist(RawData$no.weight-1)/(RawData$no.weight-1))
df <- RawData$no.weight-1
t <- qt(0.975, df)

ri_low_ki_low <- pmax(0, ri_low - t*se_ri)
ri_low_ki_up <- ri_low + t * se_ri
ri_up_ki_low <- ri_up - t*se_ri
ri_up_ki_up <- ri_up + t * se_ri

low.data <- data.frame(mi = ri_low, sdi = se_ri*sqrt(RawData$no.weight), ni = RawData$no.weight,
                       fotipus = RawData$fotipus, tipus = RawData$tipus)
up.data <- data.frame(mi = ri_up, sdi = se_ri*sqrt(RawData$no.weight), ni = RawData$no.weight,
                      fotipus = RawData$fotipus, tipus = RawData$tipus)

fit.low <- rma.uni(mi = mi, sdi = sdi, ni = ni, data = low.data, measure = "MN", method = "REML")
fit.up <- rma.uni(mi = mi, sdi = sdi, ni = ni, data = up.data, measure = "MN", method = "REML")

fit.low.hplc <- update(fit.low, subset = fotipus=="HPLC")
fit.low.hplc.fluor <- update(fit.low.hplc, subset = tipus=="Fluorescence")
fit.low.hplc.lcms <- update(fit.low.hplc, subset = tipus=="LC/MS")
fit.low.elisa <- update(fit.low, subset = fotipus=="ELISA")
fit.low.elisa.dld <- update(fit.low.elisa, subset = tipus=="DLD")
fit.low.elisa.id <- update(fit.low.elisa, subset = tipus=="Immundiagnostik")

fit.up.hplc <- update(fit.up, subset = fotipus=="HPLC")
fit.up.hplc.fluor <- update(fit.up.hplc, subset = tipus=="Fluorescence")
fit.up.hplc.lcms <- update(fit.up.hplc, subset = tipus=="LC/MS")
fit.up.elisa <- update(fit.up, subset = fotipus=="ELISA")
fit.up.elisa.dld <- update(fit.up.elisa, subset = tipus=="DLD")
fit.up.elisa.id <- update(fit.up.elisa, subset = tipus=="Immundiagnostik")

RawData$slab <- paste0(RawData$Author, " (", RawData$DATE, ")")

maxrow <- nrow(RawData)+10

plot(RawData$mean.ADMA.umol.l, 1:nrow(RawData), xlim = c(-3, 4.5), ylim = c(-5, maxrow), type = "n",
     xaxt = "n", yaxt = "n", axes = FALSE, xlab = "ADMA concentration [umol/l]", ylab = "")
sapply(1:nrow(RawData[RawData$fotipus=="HPLC",]),
       function(i) plotstudy(maxrow-2-i, c(ki_low[i], RawData$mean.ADMA.umol.l[i], ki_up[i]),
                             c(ri_low_ki_low[i], ri_low[i], ri_low_ki_up[i]),
                             c(ri_up_ki_low[i], ri_up[i], ri_up_ki_up[i])))
sapply((nrow(RawData[RawData$fotipus=="HPLC",])+1):nrow(RawData),
       function(i) plotstudy(maxrow-9-i, c( ki_low[i], RawData$mean.ADMA.umol.l[i], ki_up[i]),
                             c(ri_low_ki_low[i], ri_low[i], ri_low_ki_up[i]),
                             c(ri_up_ki_low[i], ri_up[i], ri_up_ki_up[i])))
text(-3, (maxrow-3):((maxrow-3)-length(RawData$slab[ RawData$fotipus=="HPLC" ])+1),
     RawData$slab[ RawData$fotipus=="HPLC" ], pos = 4)
text(-3, ((maxrow-3)-length(RawData$slab[RawData$fotipus=="HPLC"])-1):((maxrow-3)-length(RawData$slab[RawData$fotipus=="HPLC"])-3),
     c("HPLC", "   Fluorescence", "   LC/MS"), pos = 4, font = 3)
plotstudy2(((maxrow-3)-length(RawData$slab[ RawData$fotipus=="HPLC" ])-1),
           plotsummary(fit.hplc, fit.low.hplc, fit.up.hplc))
plotstudy2(((maxrow-3)-length(RawData$slab[ RawData$fotipus=="HPLC" ])-2),
           plotsummary(fit.hplc.fluor, fit.low.hplc.fluor, fit.up.hplc.fluor))
plotstudy2(((maxrow-3)-length(RawData$slab[ RawData$fotipus=="HPLC" ])-3),
           plotsummary( fit.hplc.lcms, fit.low.hplc.lcms, fit.up.hplc.lcms))
text(-3, ((maxrow-3)-length(RawData$slab[ RawData$fotipus=="HPLC" ])-7):1,
     RawData$slab[ RawData$fotipus=="ELISA" ], pos = 4)
text(-3.1, maxrow-1, "HPLC", font = 2, pos = 4)
text(-3.1, ((maxrow-3)-length(RawData$slab[ RawData$fotipus=="HPLC" ])-5), "ELISA", font = 2, pos = 4)
text(-3, -(1:3), c( "ELISA", "   DLD", "   Immundiagnostik" ), pos = 4, font = 3)
plotstudy2(-1, plotsummary(fit.elisa, fit.low.elisa, fit.up.elisa))
plotstudy2(-2, plotsummary(fit.elisa.dld, fit.low.elisa.dld, fit.up.elisa.dld))
plotstudy2(-3, plotsummary(fit.elisa.id, fit.low.elisa.id, fit.up.elisa.id))
text(-3, -5, "Overall", font = 2, pos = 4)
plotstudy2(-5, plotsummary( fit, fit.low, fit.up))
axis(1, seq(0, 4.5, 0.5))
dev.off()

rbind(data.frame(Type = "Overall", plotsummary(fit, fit.low, fit.up)),
      data.frame(Type = "HPLC", plotsummary(fit.hplc, fit.low.hplc, fit.up.hplc)),
      data.frame(Type = "HPLC - Fluorescence",
                 plotsummary(fit.hplc.fluor, fit.low.hplc.fluor, fit.up.hplc.fluor)),
      data.frame(Type = "HPLC - LC/MS", plotsummary(fit.hplc.lcms, fit.low.hplc.lcms, fit.up.hplc.lcms)),
      data.frame(Type = "ELISA", plotsummary(fit.elisa, fit.low.elisa, fit.up.elisa)),
      data.frame(Type = "ELISA - DLD", plotsummary(fit.elisa.dld, fit.low.elisa.dld, fit.up.elisa.dld)),
      data.frame(Type = "ELISA - Immundiagnostik",
                 plotsummary(fit.elisa.id, fit.low.elisa.id, fit.up.elisa.id)))

influence(fit)
plot(influence(fit), slab.style = 1, las = 3)
dev.off()
infl <- influence(fit)
data.frame(id = infl$ids, infl$inf, dfbs = infl$dfbs$intrcpt, influential = ifelse(infl$is.infl, "*", ""))

qqnorm(fit)
dev.off()

fit <- rma.uni(mi = mean.ADMA.umol.l, sdi = ADMA.SD, ni = no.weight, data = RawData, measure = "MN",
               slab = paste0(Author, " (", DATE, ")"), method = "REML",
               mods = ~AGE+MaleRatio+fotipus+BMI.mean)
summary(fit)
plot(influence(fit))
qqnorm(fit)

RawData <- RawData[!infl$is.infl,]

# Re-run 109 to 207

infl2 <- influence(fit)

cbind(sapply(infl$inf, range), range(infl$dfbs))
cbind(sapply(infl2$inf, range), range(infl2$dfbs))