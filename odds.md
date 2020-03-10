The Association Between Exposure and Disease Outcome
====================================================

Consider a standard epidemiological table that cross-classifies members
of a study population by their disease and exposure status. The columns
of the two-way table are the disease state (“disease”, “no disease”) and
the rows are the exposure state (“exposed”,“not exposed”).

    library(knitr)
    epitab <- matrix(c("a", "b",
                    "c", "d"), nr=2, nc=2, byrow=TRUE)
    colnames(epitab) <- c("Disease","No Disease")
    rownames(epitab) <- c("Exposed", "Not Exposed")
    kable(epitab, align="c")

<table>
<thead>
<tr class="header">
<th></th>
<th style="text-align: center;">Disease</th>
<th style="text-align: center;">No Disease</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td>Exposed</td>
<td style="text-align: center;">a</td>
<td style="text-align: center;">b</td>
</tr>
<tr class="even">
<td>Not Exposed</td>
<td style="text-align: center;">c</td>
<td style="text-align: center;">d</td>
</tr>
</tbody>
</table>

The values in the cells of this table are counts. For example, `a` is
the count of disease cases who experienced exposure. The row and column
sums are known collectively as the *marginal totals* of the table: `a+b`
and `c+d` are the row sums and `a+c` and `b+d` are the column sums. The
grand sum is `n=a+b+c+d`.

There are two broad ways that we could construct an epidemiological
table like this. First, we could select a population and follow it
through time, noting exposures and the occurance of disease. This is
known as a *cohort design*. Cohort designs are difficult, expensive, and
impractical for many infectious diseases, particularly emerging
infectious diseases. A common design begins with a group of cases of the
disease and then seeks a comparable sample of non-diseased individuals
as controls, ascertaining for both groups their exposure status. This is
known as a *case-control* design.

Prevalence, Absolute Risk, and Risk Difference
----------------------------------------------

The *prevalence* is simply the number of cases divided by the total
population: (*a* + *c*)/(*a* + *b* + *c* + *d*). The *absolute risk* or
*attack rate* is *a*/(*a* + *b*). The Risk difference is the difference
between the absolute risk of exposed and unexposed individuals:

<!--
\[ RD = \frac{a}{a+b} - \frac{c}{c+d} 
\]
-->
<img src="https://latex.codecogs.com/svg.latex?\large&space;RD&space;=&space;\frac{a}{a&plus;b}&space;-&space;\frac{c}{c&plus;d}" title="\large RD = \frac{a}{a+b} - \frac{c}{c+d}" />

In a cohort design, the risk difference is typically known as the
*attributable risk*.

Something not easily captured in our simple 2 × 2 table is *incidence*,
which is the number of new cases in some time interval divided by the
person-time units. For example, if there are 3 new cases in a week in a
population with 100 people at risk, then the incidence is 3/100. Suppose
that each case is fatal and that the average number of person-weeks
lived by those who contract the disease is 0.5 weeks. In this case, the
incidence would be slightly higher at 3/98.5 (i.e., 97 person-weeks for
those not afflicted and a total of 3 × 0.5 = 1.5 person-weeks for those
dying).

Risk Ratios and Odds Ratios
---------------------------

The *risk ratio* or *relative risk* is the ratio of the risk to those
exposed relative to the risk of those not exposed:

<!--\[ \widehat{RR} = \frac{a/(a+b)}{c/(c+d)} \]-->
<img src="https://latex.codecogs.com/svg.latex?\large&space;\widehat{RR}&space;=&space;\frac{a/(a&plus;b)}{c/(c&plus;d)}" title="\large \widehat{RR} = \frac{a/(a+b)}{c/(c+d)}" />

If the entries of our table are large enough a normal approximation for
the binomial distribution applies and we can use normal theory to
calculate standard errors and place confidence bounds around
$\\widehat{RR}$. It turns out that $\\log \\widehat{RR}$ has a sampling
distribution better approximated by a normal, so we work with that. The
standard error of $\\log \\widehat{RR}$ is

<!--\[ 
\mathrm{se}[\log \widehat{RR}] \approx \sqrt{\frac{b}{an_1} +
  \frac{d}{cn_2}}, 
\] -->
<img src="https://latex.codecogs.com/svg.latex?\large&space;\mathrm{se}[\log&space;\widehat{RR}]&space;\approx&space;\sqrt{\frac{b}{an_1}&space;&plus;&space;\frac{d}{cn_2}}" title="\large \mathrm{se}[\log \widehat{RR}] \approx \sqrt{\frac{b}{an_1} + \frac{d}{cn_2}}" />

where *n*<sub>1</sub> = *a* + *b* and *n*<sub>2</sub> = *c* + *d* are
the row sums.

Using the standard approach for calculating confidence intervals from
the normal approximation, the 95% confidence interval on the relative
risk are calculated from:

<!--\[ 
c_1 = \log \widehat{RR} - 1.96 \sqrt{\frac{b}{an_1} +
  \frac{d}{cn_2}}
\]

\[ 
c_2 = \log \widehat{RR} + 1.96 \sqrt{\frac{b}{an_1} +
  \frac{d}{cn_2}}
\]-->
<img src="https://latex.codecogs.com/svg.latex?\large&space;c_1&space;=&space;\log&space;\widehat{RR}&space;-&space;1.96&space;\sqrt{\frac{b}{an_1}&space;&plus;&space;\frac{d}{cn_2}}" title="\large c_1 = \log \widehat{RR} - 1.96 \sqrt{\frac{b}{an_1} + \frac{d}{cn_2}}" />

<img src="https://latex.codecogs.com/svg.latex?\large&space;c_2&space;=&space;\log&space;\widehat{RR}&space;&plus;&space;1.96&space;\sqrt{\frac{b}{an_1}&space;&plus;&space;\frac{d}{cn_2}}" title="\large c_2 = \log \widehat{RR} + 1.96 \sqrt{\frac{b}{an_1} + \frac{d}{cn_2}}" />

And the CI is thus,
\[*e*<sup>*c*<sub>1</sub></sup>, *e*<sup>*c*<sub>2</sub></sup>\]. Just
as a reminder, we multiply by 1.96 when calculating a 95% confidence
interval using a normal approximation because 95% of the probability
mass of a normal distribution lies within 1.96 standard deviations of
the mean.

    round(c(pnorm(1.96), pnorm(-1.96)),3)

    ## [1] 0.975 0.025

    round(c(qnorm(0.975), qnorm(0.025)),3)

    ## [1]  1.96 -1.96

The relative risk is a pretty intuitive idea, but a problem with it is
that it is constrained by the denominator. If *c*/(*c* + *d*) = 0.5,
then the biggest the relative risk could be is $\\widehat{RR}=2$. This
is because the biggest *a*/(*a* + *b*) can be is unity since this would
happen when *b* = 0, therefore $\\widehat{RR} = 2$. In addition, we
often (usually?) don’t have a cohort design for our data collection.
Without the prospective design of a cohort study, relative risks don’t
make much sense since we shouldn’t believe that *a*/(*a* + *b*) is a
good estimator of the risk associated with the exposure of interest. We
can get around this by using an odds ratio.

The *odds ratio* is the ratio of odds of exposure among cases to those
among controls:

<!--\[ OR = \frac{a/c}{b/d} = \frac{ad}{bc} \] -->
<img src="https://latex.codecogs.com/svg.latex?\large&space;OR&space;=&space;\frac{a/c}{b/d}&space;=&space;\frac{ad}{bc}" title="\large OR = \frac{a/c}{b/d} = \frac{ad}{bc}" />

It is the cross-product of the cells of our epidemiological table. For
rare diseases, the OR is a very good approximation of the RR. This is
because, if the disease is rare, *a* ≪ *b* and *c* ≪ *d*. This means:

<!--\[ \frac{a/(a+b)}{c/(c+d)} \approx \frac{a/b}{c/d} = \frac{ad}{bc}.\]-->
<img src="https://latex.codecogs.com/svg.latex?\large&space;\frac{a/(a&plus;b)}{c/(c&plus;d)}&space;\approx&space;\frac{a/b}{c/d}&space;=&space;\frac{ad}{bc}" title="\large \frac{a/(a+b)}{c/(c+d)} \approx \frac{a/b}{c/d} = \frac{ad}{bc}" />

We can calculate confidence intervals on an odds ratio using the normal
approximation (assuming the cells are large enough). Again, we work with
the logarithm of the measure of association. The standard error of the
logarithm of the odds ratio is,

<!--\[ \mathrm{se}[\log OR] = \sqrt{\frac{1}{a} + \frac{1}{b} +
  \frac{1}{c} + \frac{1}{d}}, 
\]-->
<img src="https://latex.codecogs.com/svg.latex?\large&space;\mathrm{se}[\log&space;OR]&space;=&space;\sqrt{\frac{1}{a}&space;&plus;&space;\frac{1}{b}&space;&plus;&space;\frac{1}{c}&space;&plus;&space;\frac{1}{d}}" title="\large \mathrm{se}[\log OR] = \sqrt{\frac{1}{a} + \frac{1}{b} + \frac{1}{c} + \frac{1}{d}}" />

and the 95% confidence intervals for the log-odds ratio are thus:

<!--\[ c_1 = \log OR - 1.96 \sqrt{\frac{1}{a} + \frac{1}{b} +
  \frac{1}{c} + \frac{1}{d}} 
\]

\[ c_2 = \log OR + 1.96 \sqrt{\frac{1}{a} + \frac{1}{b} +
  \frac{1}{c} + \frac{1}{d}} 
\]-->
<img src="https://latex.codecogs.com/svg.latex?\large&space;c_1&space;=&space;\log&space;OR&space;-&space;1.96&space;\sqrt{\frac{1}{a}&space;&plus;&space;\frac{1}{b}&space;&plus;&space;\frac{1}{c}&space;&plus;&space;\frac{1}{d}}" title="\large c_1 = \log OR - 1.96 \sqrt{\frac{1}{a} + \frac{1}{b} + \frac{1}{c} + \frac{1}{d}}" />

<img src="https://latex.codecogs.com/svg.latex?\large&space;c_2&space;=&space;\log&space;OR&space;&plus;&space;1.96&space;\sqrt{\frac{1}{a}&space;&plus;&space;\frac{1}{b}&space;&plus;&space;\frac{1}{c}&space;&plus;&space;\frac{1}{d}}" title="\large c_2 = \log OR + 1.96 \sqrt{\frac{1}{a} + \frac{1}{b} + \frac{1}{c} + \frac{1}{d}}" />

Back-transforming, we get the confidence interval on the unit scale of
\[*e*<sup>*c*<sub>1</sub></sup>, *e*<sup>*c*<sub>2</sub></sup>\].

SARS
----

As China increased its cooperation with WHO and other health
organizations, the details of the early epidemic in Guangdong Province
became known. The first case in Guangdong to have signs of the disease
consistent with the WHO case definition occurred in Foshan city.
Unfortunately, the patient died before samples could be collected for
virological investigation. The second case presented itself on 17
December 2002. A chef from Heyuan presented with severe atypical
pneumonia. It seems that this chef worked in a restaurant that
specialized in exotic game and had regular contact with several species
of live caged animals.

Sampling effort thus turned to the exotic animal markets of Guangdong. A
team from the University of Hong Kong and the Guangdong CDC sampled
animals from a market in Shenzhen. A total of seven common game animals
to the market were subjected to both nasal and fecal sampling. Virus was
repeatedly detected in palm civets. These researchers also took nasal
samples from the animal traders, in addition to hospital workers,
Guandong CDC employees, and health visitors to a clinic as a controls. A
total of 70 people in the market were SCoV+, 66 of which were animal
traders. 579 were SCoV-, of which 442 were animal traders ([Yu et
al. 2003](https://www.cdc.gov/mmwr/preview/mmwrhtml/mm5241a2.htm)).

    cc <- matrix(c(66L, 442L,
                    4L, 133L), nr=2, nc=2, byrow=TRUE)
    colnames(cc) <- c("SCoV+","SCoV-")
    rownames(cc) <- c("Animal Trader", "Hospital Worker")
    cc

    ##                 SCoV+ SCoV-
    ## Animal Trader      66   442
    ## Hospital Worker     4   133

    kable(cc, align="c")

<table>
<thead>
<tr class="header">
<th></th>
<th style="text-align: center;">SCoV+</th>
<th style="text-align: center;">SCoV-</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td>Animal Trader</td>
<td style="text-align: center;">66</td>
<td style="text-align: center;">442</td>
</tr>
<tr class="even">
<td>Hospital Worker</td>
<td style="text-align: center;">4</td>
<td style="text-align: center;">133</td>
</tr>
</tbody>
</table>

Calculate the odds ratio and confidence intervals.

    cc <- matrix(c(66, 442,
                    4, 133), nr=2, nc=2, byrow=TRUE)
    colnames(cc) <- c("SCoV+","SCoV-")
    rownames(cc) <- c("Animal", "Hospital")
    cc

    ##          SCoV+ SCoV-
    ## Animal      66   442
    ## Hospital     4   133

    (row.totals <- apply(cc,1,sum))

    ##   Animal Hospital 
    ##      508      137

    (col.totals <- apply(cc,2,sum))

    ## SCoV+ SCoV- 
    ##    70   575

    ## absolute risk
    (ar <- cc[1,1]/row.totals[1])

    ##    Animal 
    ## 0.1299213

    # risk difference
    (rd <- cc[1,1]/row.totals[1] - cc[2,1]/row.totals[2])

    ##    Animal 
    ## 0.1007242

    ## rate ratio (watch the parentheses!)
    (rr <- (cc[1,1]/row.totals[1])/(cc[2,1]/row.totals[2]))

    ##   Animal 
    ## 4.449803

    #odds ratio
    (or <- (cc[1,1]*cc[2,2])/(cc[1,2]*cc[2,1]))

    ## [1] 4.964932

    se <- sqrt(1/cc[1,1]+1/cc[2,2]+1/cc[1,2]+1/cc[2,1])
    #lower
    exp(log(or)-(1.96*se))

    ## [1] 1.776584

    #upper
    exp(log(or)+(1.96*se))

    ## [1] 13.87525

Fisher Exact Test
-----------------

When the cells are small and the normal approximation won’t work. Note
that the odds ratio will be undefined if *c* = 0 but you could imagine a
situation where there is a really strong association between exposure
and disease, and where the sample size is relatively small, where there
are no cases in unexposed people (i.e., *c* = 0). We can use Fisher’s
exact test in situations like this. Fisher’s exact test calculates the
probabilities of the table directly from a hypergeometric distribution.

[Keele et al. (2009)](http://dx.doi.org/10.1038/nature08200) examined
the demographic effects of SIV infection among wild chimpanzees. They
report that, in the nine years of their study, four of nine SIV-infected
females give birth, whereas 22 of 30 uninfected females gave birth in
that same period. We can set up an epidemiological table from this
information and test for an association between infection status and
fertility.

    gombe.fert <- matrix(c(4,22,5,8),2,2,byrow=TRUE)
    colnames(gombe.fert) <- c("SIV+","SIV-")
    rownames(gombe.fert) <- c("Birth", "No Birth")
    gombe.fert

    ##          SIV+ SIV-
    ## Birth       4   22
    ## No Birth    5    8

    fisher.test(gombe.fert, alternative="less")

    ## 
    ##  Fisher's Exact Test for Count Data
    ## 
    ## data:  gombe.fert
    ## p-value = 0.1146
    ## alternative hypothesis: true odds ratio is less than 1
    ## 95 percent confidence interval:
    ##  0.000000 1.401445
    ## sample estimates:
    ## odds ratio 
    ##  0.3014952

We can see that while the odds ratio is substantially less than one
(meaning that SIV infection appears to reduce fertility), the confidence
intervals on the odds ratio cross one, suggesting that this result can
not be distinguished from a null association. Keele et al. (2009)
actually used a different method – that accounted for chimp-years of
exposure – to calculate the association. This example is just for
illustration of the method.

Logistic Regression
-------------------

Logistic regression gives you the same answer as cross-multiplying a
two-way table. We can analyze the SARS table using a logistic regression
model (i.e., a binomial-family `glm`). In `R`, the left-hand side of a
formula for a binomial-family `glm` can be specified in one of several
ways. One of those is to pass a two-column matrix in which the first
column is the number of “successes” and the second column is the number
of “failures”. This happens to be exactly the form of our 2 × 2 table,
which is convenient. All we need is to construct a factor for exposure
and call the model:

    expose <- factor(c("yes","no"))
    out <- glm(cc ~ expose, family=binomial)
    summary(out)

    ## 
    ## Call:
    ## glm(formula = cc ~ expose, family = binomial)
    ## 
    ## Deviance Residuals: 
    ## [1]  0  0
    ## 
    ## Coefficients:
    ##             Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)  -3.5041     0.5075  -6.905 5.02e-12 ***
    ## exposeyes     1.6024     0.5243   3.056  0.00224 ** 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for binomial family taken to be 1)
    ## 
    ##     Null deviance:  1.4450e+01  on 1  degrees of freedom
    ## Residual deviance: -3.1086e-14  on 0  degrees of freedom
    ## AIC: 13.127
    ## 
    ## Number of Fisher Scoring iterations: 4

    ## odds
    exp(out$coef[2])

    ## exposeyes 
    ##  4.964932

    ## now compare to cross-multiplication
    (or <- (cc[1,1]*cc[2,2])/(cc[1,2]*cc[2,1]))

    ## [1] 4.964932

We can see that the two ways of calculating the odds-ratio give
essentially identical answers.
