################################################################################
# Part of the R/EpiILM package
#
# AUTHORS:
#         Vineetha Warriyar. K. V. <vineethawarriyar.kod@ucalgary.ca> and
#         Rob Deardon <robert.deardon@ucalgary.ca>
#
# HISTORY:
#          Version 1.0: 2017-04-14
#          Version 1.1: 2017-04-17
#
# Free software under the terms of the GNU General Public License, version 2,
# a copy of which is available at http://www.r-project.org/Licenses/.
################################################################################

epidic <- function(burnin, niter, LLchain, LLpostmean) {
 
 # loglikelihood values after burn-in
  loglike <- LLchain[burnin:niter]
  D.bar   <- -2 * mean(loglike)
  D.hat   <- -2 * LLpostmean
  pD      <- D.bar - D.hat
  dic     <- pD + D.bar
  return(dic = dic)
}

