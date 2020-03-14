The Next-Generation Matrix
==========================

The basic reproduction number, *R*<sub>0</sub>, is the expected number
of secondary cases caused by a single typical case at the outset of an
epidemic. When there are multiple types of hosts, we need to think hard
about what “typical” means. Fortunately, as we discuss in the [notes on
*R*<sub>0</sub>](https://eehh-stanford.github.io/gceid/Jones-R0-notes2020.pdf),
there is a simple way to extend the calculation of *R*<sub>0</sub> to
structured populations first introduced by Diekmann, Heesterbeek, and
Metz (1990). The *next-generation matrix* is a square matrix **G** in
which the *ij*th element of **G**, *g*<sub>ij</sub>, is the expected
number of secondary infections of type *i* caused by a single infected
individual of type *j*, again assuming that the population of type *i*
is entirely susceptible. **G** is a positive matrix and this means that,
from the Perron-Frobenius Theorem, there will be one eigenvalue that is
real, positive, and strictly greater than all others. We call this
eigenvalue “the dominant eigenvalue” of matrix **G**. *R*<sub>0</sub> is
simply the dominant eigenvalue of **G**.

Working with 2 × 2 Matrices
---------------------------

When you’re going to construct a lot of next-generation matrices, it can
be handy to write a quick function to facilitate this. The `R` function
`matrix()` requires an input vector (subject to standard recycling
rules), the number of rows, and number of columns. Most people will also
want to add the argument `byrow=TRUE` to override the default columnwise
filling of matrices (a holdover from FORTRAN).

    make2x2 <- function(a,b,c,d) matrix( c(a,b,c,d), nr=2, nc=2, byrow=TRUE)
    (g <- make2x2(15,2,4,8))

    ##      [,1] [,2]
    ## [1,]   15    2
    ## [2,]    4    8

    eigen(g)

    ## eigen() decomposition
    ## $values
    ## [1] 16  7
    ## 
    ## $vectors
    ##           [,1]       [,2]
    ## [1,] 0.8944272 -0.2425356
    ## [2,] 0.4472136  0.9701425

The Fastest Way to Reduce *R*<sub>0</sub> in a Multi-Host Epidemic
==================================================================

We want to bring the effective reproduction number *R*<sub>0</sub> below
one in a multi-host epidemic. Clearly, reducing the number of infections
in any cell in the next-generation matrix will reduce *R*<sub>0</sub>,
but which will provide the greatest reduction for a given effort? To
figure this out, we calculate the *sensitivities* of *R*<sub>0</sub> to
a change in the various entries of the next-generation matrix. These are
the partial derivatives (i.e., the rate of change in *R*<sub>0</sub>
with respect to *g*<sub>*i**j*</sub> holding everything else constant).

<img src="https://latex.codecogs.com/svg.latex?s_{ij}&space;=&space;\frac{\partial&space;R_0}{\partial&space;g_{ij}}" title="s_{ij} = \frac{\partial R_0}{\partial g_{ij}}" />
<!-- \[ s_{ij} = \frac{\partial R_0}{\partial g_{ij}} \] -->

    calc.sens <- function(G){
      ev <- eigen(G)
      lmax <- which(Re(ev$values)==max(Re(ev$values)))
      U <- ev$vectors
      u <- abs(Re(U[,lmax]))
      V <- solve(Conj(U))
      v <- abs(Re(V[lmax,]))
      
      s <- v%o%u
      return(s)
    }

Chagas Disease Example
----------------------

Consider a three-species model introduced by Cohen and Gürtler (2001).
The three species included in this household-based model are humans,
triatomine bugs, and dogs. Define **G**:

<img src="https://latex.codecogs.com/svg.latex?\mathbf{G}&space;=&space;\left[&space;\begin{array}{ccc}&space;0&space;&&space;g_{B&space;\leftarrow&space;D}&space;&&space;g_{B&space;\leftarrow&space;H}&space;\\&space;g_{D&space;\leftarrow&space;B}&space;&&space;0&space;&&space;0&space;\\&space;g_{H&space;\leftarrow&space;B}&space;&&space;0&space;&&space;0&space;\end{array}&space;\right]," title="\mathbf{G} = \left[ \begin{array}{ccc} 0 & g_{B \leftarrow D} & g_{B \leftarrow H} \\ g_{D \leftarrow B} & 0 & 0 \\ g_{H \leftarrow B} & 0 & 0 \end{array} \right]," />

<!-- \[ \mathbf{G} = \left[ \begin{array}{ccc}
            0 & g_{B \leftarrow D} & g_{B \leftarrow H} \\
            g_{D \leftarrow B} & 0 & 0 \\
            g_{H \leftarrow B}  & 0 & 0
            \end{array} \right],
\]
-->
where (e.g.)
<img src="https://latex.codecogs.com/svg.latex?g_{B&space;\leftarrow&space;D}" title="g_{B \leftarrow D}" />
is the dog-infecting-bug reproduction number.

Cohen and Gürtler (2001) give some data that allow us to estimate the
next generation matrix. We fill in the missing parameters with educated
guesses (more on this later) and calculate the sensitivities.

    G <- matrix( c(0,122.01,7.47, 0.1992,0,0, 0.1992,0,0), nr=3, nc=3, byrow=TRUE)
    # first, what is R_0?
    eigen(G)$values

    ## [1] -5.078623  5.078623  0.000000

    (S <- calc.sens(G))

    ##            [,1]       [,2]       [,3]
    ## [1,]  0.5000000 0.01961161 0.01961161
    ## [2,] 12.0121133 0.47115385 0.47115385
    ## [3,]  0.7354355 0.02884615 0.02884615

We can see that *R*<sub>0</sub> = 5.08 and that, by far, the element of
**G** with the highest sensitivity is the bug-to-dog entry. Note that
the ordering of the eigenvalues in the output actually lists the
negative one first. This is why can’t just take the first eigenvalue but
need the line in the code to find the dominant eigenvalue (i.e., the one
that is positive, real, and strictly greater than all the others).

Elasticities
============

Sensitivities tell us about the absolute change in *R*<sub>0</sub> given
a small change in one of the entries in **G**. Sometimes we want to know
about the proportional sensitivities though. Use a trick invented by
Carl Jacobi in the mid-nineteenth century. Define a proportional
sensitivity (or *elasticity*) as

<!-- \[ e_{ij} = \frac{g_{ij}}{R_0} \cdot \frac{\partial R_0}{\partial
  g_{ij}} = \frac{\partial \log(R_0)}{\partial \log(g_{ij})}.\] -->
<img src="https://latex.codecogs.com/svg.latex?e_{ij}&space;=&space;\frac{g_{ij}}{R_0}&space;\cdot&space;\frac{\partial&space;R_0}{\partial&space;g_{ij}}&space;=&space;\frac{\partial&space;\log(R_0)}{\partial&space;\log(g_{ij})}" title="e_{ij} = \frac{g_{ij}}{R_0} \cdot \frac{\partial R_0}{\partial g_{ij}} = \frac{\partial \log(R_0)}{\partial \log(g_{ij})}" />

Elasticities answer the following question: if *g*<sub>*i**j*</sub>
increases by 1%, by what percentage will *R*<sub>0</sub> increase?
Elasticities also have the convenient property of summing to one across
all elements of the next generation matrix. They thus represent the
fraction of the total potential reduction embodied in each cell of the
next generation matrix (assuming everything else stays the same).

    calc.elas <- function(G){
      ev <- eigen(G)
      lmax <- which(Re(ev$values)==max(Re(ev$values)))
      R0 <- ev$values[lmax]
      U <- ev$vectors
      u <- abs(Re(U[,lmax]))
      V <- solve(Conj(U))
      v <- abs(Re(V[lmax,]))
      
      s <- v%o%u
      e <- s*G/R0
      return(e)
    }

    (E <- calc.elas(G))

    ##            [,1]      [,2]       [,3]
    ## [1,] 0.00000000 0.4711538 0.02884615
    ## [2,] 0.47115385 0.0000000 0.00000000
    ## [3,] 0.02884615 0.0000000 0.00000000

So our elasticity matrix looks like this:

<!-- \[ \mathbf{E} = \left[ \begin{array}{ccc}
0.00 & 0.47 & 0.03 \\
0.47 & 0.00 & 0.00 \\
0.03 & 0.00 & 0.00 
 \end{array} \right].
\] -->
<img src="https://latex.codecogs.com/svg.latex?\mathbf{E}&space;=&space;\left[&space;\begin{array}{ccc}&space;0.00&space;&&space;0.47&space;&&space;0.03&space;\\&space;0.47&space;&&space;0.00&space;&&space;0.00&space;\\&space;0.03&space;&&space;0.00&space;&&space;0.00&space;\end{array}&space;\right]" title="\mathbf{E} = \left[ \begin{array}{ccc} 0.00 & 0.47 & 0.03 \\ 0.47 & 0.00 & 0.00 \\ 0.03 & 0.00 & 0.00 \end{array} \right]" />

Transmission from bugs to dogs and dogs to bugs accounts for more than
94% of the elasticity. Clearly, we need to stop the canine transmission
cycle! The elasticities of all pairs (e.g., bugs to dogs/dogs to bugs)
are identical because of a feature of elasticities identified by
Groenendael et al. (1994). In brief, the sum of elasticities of the
incoming transitions to a vertex in the transmission graph will equal
the sum of the elasticities of the outgoing transitions. Because humans
and dogs don’t infect each other, this means that the elasticities
involving bugs and the two mammal species (i.e., humans, dogs) will be
equal to each other. More later.

Other Structured Epidemic Models
================================

All of the examples so far have involved multi-host infections. There
are many important ways that disease transmission is structured *within*
human populations as well. For example, the transmission of respiratory
infections like influenza is quite strongly age-structured (Hill and
Longini 2003). Indeed, social categories can have a profound impact on
the transmission dynamics of epidemics. Kucharski et al. (2014) show
that during the 2009 H1N1 Influenza pandemic, infection risk among
people in Hong Kong was more influenced by the average reported social
mixing behavior of an individual’s age group, than it was by their own
reported contacts. For easily-transmitted respiratory pathogens like
influenza A, the aggregate patterns can overwhelm the particular
behavior of individuals.

We can use sensitivity/elasticity analysis of the next-generation matrix
to help us determine the categories of people on whom we should
concentrate our epidemic-control efforts if we want to have maximum
effect for a given effort.

Hill and Longini (2003) present a next-generation matrix for influenza-A
transmission calculated from data presented in Longini, Ackerman, and
Elveback (1978):

<!-- \[ \mathbf{G} = \left[ \begin{array}{ccccc}
    0.60 & 0.10 & 0.10 & 0.10 & 0.10 \\
    0.20 & 1.70 & 0.30 & 0.20 & 0.20 \\
    0.40 & 0.30 & 0.50 & 0.40 & 0.30 \\
    0.20 & 0.10 & 0.30 & 0.20 & 0.10 \\
    0.10 & 0.10 & 0.10 & 0.10 & 0.10 
               \end{array} \right] 
\]
-->
<img src="https://latex.codecogs.com/svg.latex?\mathbf{G}&space;=&space;\left[&space;\begin{array}{ccccc}&space;0.60&space;&&space;0.10&space;&&space;0.10&space;&&space;0.10&space;&&space;0.10&space;\\&space;0.20&space;&&space;1.70&space;&&space;0.30&space;&&space;0.20&space;&&space;0.20&space;\\&space;0.40&space;&&space;0.30&space;&&space;0.50&space;&&space;0.40&space;&&space;0.30&space;\\&space;0.20&space;&&space;0.10&space;&&space;0.30&space;&&space;0.20&space;&&space;0.10&space;\\&space;0.10&space;&&space;0.10&space;&&space;0.10&space;&&space;0.10&space;&&space;0.10&space;\end{array}&space;\right]" title="\mathbf{G} = \left[ \begin{array}{ccccc} 0.60 & 0.10 & 0.10 & 0.10 & 0.10 \\ 0.20 & 1.70 & 0.30 & 0.20 & 0.20 \\ 0.40 & 0.30 & 0.50 & 0.40 & 0.30 \\ 0.20 & 0.10 & 0.30 & 0.20 & 0.10 \\ 0.10 & 0.10 & 0.10 & 0.10 & 0.10 \end{array} \right]" />

where the classes are: (1) pre-school age, (2) school age, (3) young
adult, (4) mid-adult, (5) old adult.

    G <- matrix( 
      c(0.60,0.10,0.10,0.10,0.10,
        0.20,1.70,0.30,0.20,0.20,
        0.40,0.30,0.50,0.40,0.30,
        0.20,0.10,0.30,0.20,0.10,
        0.10,0.10,0.10,0.10,0.10), nr=5, nc=5,
      byrow=TRUE)

    (E <- calc.elas(G))

    ##             [,1]        [,2]        [,3]        [,4]         [,5]
    ## [1,] 0.010155405 0.013880374 0.004378325 0.001901094 0.0012386153
    ## [2,] 0.011068898 0.771575605 0.042949451 0.012432602 0.0081001829
    ## [3,] 0.006770904 0.041645019 0.021893672 0.007605088 0.0037161936
    ## [4,] 0.002477231 0.010157612 0.009612131 0.002782429 0.0009064146
    ## [5,] 0.001081376 0.008868129 0.002797298 0.001214603 0.0007913475

    round(E,2)

    ##      [,1] [,2] [,3] [,4] [,5]
    ## [1,] 0.01 0.01 0.00 0.00 0.00
    ## [2,] 0.01 0.77 0.04 0.01 0.01
    ## [3,] 0.01 0.04 0.02 0.01 0.00
    ## [4,] 0.00 0.01 0.01 0.00 0.00
    ## [5,] 0.00 0.01 0.00 0.00 0.00

We can see that, overwhelmingly, the most important task should be
controlling transmission between school-age children. This should not be
that surprising, given the fact that the entry for the **G**-matrix for
within-school-age transmission is more that three times higher than for
any other pairwise value. Note that the elasticity for this transmission
is over 19 times higher than any other dyad in the population.

Longini, Ackerman, and Elveback (1978) note that different infections
will have different optimal patterns of control and that these depend
critically on the age-structure of transmission for the different
pathogens.

Simple Model for Disease with High- and Low-Risk Subpopulations
---------------------------------------------------------------

As we saw in the [notes on structured epidemic
models](https://eehh-stanford.github.io/gceid/struct.html), reducing
contacts with the low-risk/high-activity subpopulation is the most
effective way to flatten the epidemic curve.

    # beta matrix: (0.5,, 0.15, 0.05, 0.5)
    # population structure is 2/3 low-risk, 1/3 high-risk
    # everyone recovers at rate 1/5
    (G <- make2x2(0.5*.67, 0.15*0.33, 0.05*0.67, 0.5*0.33)*5)

    ##        [,1]   [,2]
    ## [1,] 1.6750 0.2475
    ## [2,] 0.1675 0.8250

    (E <- calc.elas(G))

    ##           [,1]       [,2]
    ## [1,] 0.9253696 0.02555400
    ## [2,] 0.0255540 0.02352238

We see that, overwhelmingly, the best epidemic-reduction bang for our
public-health buck comes from reducing the contacts within the low-risk
population. This means social distancing.

A generalization that emerges from applying this sort of analysis to a
range of systems is that the most efficient method for epidemic control
is to reduce the number of infections in the the subpopulation that has
the most infections. This makes sense: you want to reduce the largest
source of potential infections in the population. Examples include:

-   Reducing transmission among the high-risk segment (e.g., sex
    workers, injection-drug users), often known as “the core,” for
    sexually-transmitted infections like gonorrhea (Hethcote and
    Yorke 1984).
-   Reducing transmission among school-age children for influenza
-   Reducing transmission among domestic fowl in bird-flu
-   Reducing transmission from triatomine bugs to peri-domestic dogs in
    Chagas disease - Reducing transmission within the
    low-risk/high-activity subpopulation in a [simple model of COVID-19
    dynamics](https://eehh-stanford.github.io/gceid/struct.html)

Going Further
=============

Caswell (2019) provides a comprehensive account of the use of
sensitivity analysis in matrix population models. The book is
open-access.

References
==========

Caswell, H. 2019. *Sensitivity Analysis: Matrix Methods in Demography
and Ecology*. Demographic Research Monographs (a Series of the Max
Planck Institute for Demographic Research). Cham, Switzerland: Springer.
doi:[10.1007/978-3-030-10534-1\_1](https://doi.org/10.1007/978-3-030-10534-1_1).

Cohen, J.E., and R.E. Gürtler. 2001. “Modeling Household Transmission of
American Trypanosomiasis.” *Science* 293: 694–98.
doi:[10.1126/science.1060638](https://doi.org/10.1126/science.1060638).

Diekmann, O., J. A. P. Heesterbeek, and J. A. J. Metz. 1990. “On the
Definition and the Computation of the Basic Reproduction Ratio
*R*<sub>0</sub> in Models for Infectious Diseases in Heterogeneous
Populations.” *Journal of Mathematical Biology* 28 (4): 365–82.
doi:[10.1007/BF00178324](https://doi.org/10.1007/BF00178324).

Groenendael, Jan van, Hans de Kroon, Susan Kalisz, and Shripad
Tuljapurkar. 1994. “Loop Analysis: Evaluating Life History Pathways in
Population Projection Matrices.” *Ecology* 75 (8): 2410–5.
doi:[10.2307/1940894](https://doi.org/10.2307/1940894).

Hethcote, H. W., and J. A. Yorke. 1984. *Gonorrhea: Transmission
Dynamics and Control*. Vol. 56. Lecture Notes in Biomathematics.

Hill, A. N., and I. M. Longini. 2003. “The Critical Vaccination Fraction
for Heterogeneous Epidemic Models.” *Mathematical Biosciences* 181 (1):
85–106.
doi:[10.1016/S0025-5564(02)00129-3](https://doi.org/10.1016/S0025-5564(02)00129-3).

Kucharski, Adam J., Kin O. Kwok, Vivian W. I. Wei, Benjamin J. Cowling,
Jonathan M. Read, Justin Lessler, Derek A. Cummings, and Steven Riley.
2014. “The Contribution of Social Behaviour to the Transmission of
Influenza A in a Human Population.” *PLOS Pathogens* 10 (6): e1004206.
doi:[10.1371/journal.ppat.1004206](https://doi.org/10.1371/journal.ppat.1004206).

Longini, Ira M., Eugene Ackerman, and Lila R. Elveback. 1978. “An
Optimization Model for Influenza a Epidemics.” *Mathematical
Biosciences* 38 (1): 141–57.
doi:[10.1016/0025-5564(78)90023-8](https://doi.org/10.1016/0025-5564(78)90023-8).
