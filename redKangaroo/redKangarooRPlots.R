library(ggplot2)

# the redKangarooInt4+5*3+4.psmc run
kanga1 = read.delim("/media/alex/My Passport/bioinformaticsScholarship/Data/redKangarooInt4+5*3+4.txt")
# the redKangarooInt40.psmc run
# THIS ONE JUST GAVE CONSTANT POP
# kanga2 = read.delim("/media/alex/My Passport/bioinformaticsScholarship/Data/redKangarooInt40.txt")


# recall that time (x) must be divided by 2 to give in terms on N_0 generations having past
# taking log of time for clarity
ggplot(data=kanga1, mapping=aes(x = log(t_k/2), y = lambda_k)) + geom_step(colour="red") + labs(xlab("log(t_k/2) (log N_0 generations)"), ylab("lambda_k")) + ggtitle("History of the red kangaroo's population dynamics") # + geom_step(data = kanga2, mapping=aes(x=log(t_k/2), y=lambda_k), colour="black")
ggplot(data=kanga1, mapping=aes(x = t_k/2, y = lambda_k)) + geom_step(colour="red") + labs(xlab("t_k/2 (N_0 generations)"), ylab("lambda_k")) + ggtitle("History of the red kangaroo's population dynamics")

### Can we work out when back in time this was? NOTE: STRONG POSSIBILITY THIS IS COMPLETELY WRONG

# According to the Department of the Environment https://www.environment.gov.au/biodiversity/wildlife-trade/natives/wild-harvest/kangaroo-wallaby-statistics/kangaroo-population (accessed 3/2/15), the 2011 population estimates for kangaroos within the commercial harvest areas was 11514298.
# Would it be reasonable to divide by 2 to get effective population size? I have no idea.
N_0 = ceiling(11514298/2) # I haven't found an upper and lower estimate at the moment.
# Generation time: 7-10 years (Dawson 2012) http://www.environment.nsw.gov.au/resources/threatenedspecies/determinations/PDRedKangReject.pdf
gen_low = 7
gen_hi = 10

# let's unscale time. At the moment, 1 unit of scaled time is equal to N_0 generations having passed.
kanga1$t_unscaled_low = kanga1$t_k/2 * N_0 * gen_low
kanga1$t_unscaled_hi = kanga1$t_k/2 * N_0 * gen_hi
kanga1$pop_hist = kanga1$lambda_k * N_0

# when was the LGM... according to wikipedia lol,
LGM_start = 26500
LGM_end = 19000

# let's see if the population dynamics can be matched with the LGM?
ggplot(data=kanga1, mapping=aes(x = log(t_unscaled_low), y = pop_hist)) + geom_step(colour="red") + labs(xlab("log(years in past)"), ylab("population size")) + ggtitle("The red kangaroo's population dynamics, and the LGM") + geom_step(data = kanga1, mapping=aes(x=log(t_unscaled_hi), y=pop_hist), colour="black") + geom_vline(xintercept = log(LGM_end)) + geom_vline(xintercept = log(LGM_start))
