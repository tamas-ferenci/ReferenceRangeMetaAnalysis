# Meta-analysis of reference ranges

The method is inspired by that of 'Rhinomanometric reference intervals for normal total nasal airflow resistance' (Rhinology. 2014 Dec;52(4):292-9. doi: 10.4193/Rhino13.220.), and we wish to express our gratitude to Martin Hellmich (Institut für Medizinische Statistik, Informatik und Epidemiologie, Universitätsklinikum Köln) for selflessly discussing the methodology.

Our paper where the present approach is used: Balázs Németh, Zénó Ajtay, László Hejjel, Tamás Ferenci, Zoltán Ábrám, Edit Murányi, István Kiss. The issue of plasma asymmetric dimethylarginine reference range -- A systematic review and meta-analysis. (PLoS One. 2017 May 11;12(5):e0177493. doi: 10.1371/journal.pone.0177493.) Available at: [https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0177493](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0177493).

## Methodology

Normal approximation was used for the meta-analysis, both for the mean and for the reference interval.

For studies where only median and IQR or median and minimum/maximum were given, the mean and SD were approximated using the method of Wan et al (Xiang Wan, Wenqian Wang, Jiming Liu, Tiejun Tong. Estimating the sample mean and standard deviation from the sample size, median, range and/or interquartile range. BMC Med Res Methodol. 2014 Dec 19;14:135. doi: 10.1186/1471-2288-14-135.) which should be preferred over the method of Hozo et al.

Confidence interval for the mean was calculated as mean ± 1.96·SD$/\sqrt{n}$, and was visualized with diamonds.

Reference interval was calculated as mean ± 1.96·SD (under the assumption of normality, this has a coverage of 95%). Confidence interval for the endpoints of the reference interval was also calculated with normal approximation of Bland and Altmann (J M Bland, D G Altman. Measuring agreement in method comparison studies. Stat Methods Med Res. 1999 Jun;8(2):135-60.doi: 10.1177/096228029900800204.).

In short, the mean and the variance are independent for normal distribution, thus the variance of mean ± 1.96·SD is the variance of mean + 1.96$^2$ times the variance of SD. The former can be estimated as variance$/\sqrt{n}$ (as sample variance is a consistent estimator of the population variance). As far as the latter is concerned, SD follows a $\chi$-distribution (as variance follows a $\chi^2$-distribution) after appropriate scaling, namely $\frac{n-1}{\sigma^2}SD^2\sim \chi_{n-1}^2$, thus the variance of SD can be approximated as
\[
  \frac{SD^2}{n-1}\mathbb{D}^2\left(\chi_{n-1}\right)=\frac{SD^2}{n-1}\cdot\left[n-1-2\left(\frac{\Gamma\left(\frac{n}{2}\right)}{\Gamma\left(\frac{n-1}{2}\right)}\right)^2\right]=SD^2\cdot\left[1-\frac{2}{n-1}\left(\frac{\Gamma\left(\frac{n}{2}\right)}{\Gamma\left(\frac{n-1}{2}\right)}\right)^2\right].
\]

Confidence interval was then constructed using the quantiles of $t$-distribution with $n-1$ degrees of freedom (truncated at zero if it were negative); it will be visualized as shaded areas around the endpoints of the reference interval.

Normality can be very roughly assessed -- without detailed distribution data -- by comparing mean with median (where both was given) and median with the midpoint of IQR (when those were given). Also, studies with non-small samples allow for central limit theorem to become effective (on one level higher, for the means).

Extreme heterogeneity suggests that the results of the random effects models should be used instead of the fixed effects model. The models were estimated using restricted maximum likelihood. Effects of moderator variables were studied with standard meta-regression approach.

## Scripts

The script is available for download here, and complete reproduction is possible with the raw data available for download here.