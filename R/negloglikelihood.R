#' Calculates the negative log-likelihood of a chromosome or ploidy model
#'
#' @param log.par A vector of parameter values in ln scale. Length of the vector depends on model defined by Q.FUN
#' @param phy Binary rooted phylogenetic tree
#' @param tip.values A data frame with at least two columns first column has tip labels in the same order than phylogenetic tree, second column ploidy or chromosome numbers depending on Q.FUN
#' @param pi.0 Probability distribution to be used at the root of phylogeny
#' @param Q.FUN A function that calculates the Q-matrix of the parameters of interest
#' @param ... Extra parameters called in Q.FUN
#' @return Result is the negative log-likelihood of parameters log.par given the model in Q.FUN
#' @author Rosana Zenil-Ferguson
#' @export
negloglikelihood <- function(log.par, phy, tip.values, pi.0,Q.FUN, ...) {
	Q.mat<-Q.FUN(log.par, ...) #Q.mat<-Q_bichrom(log.theta=log.params, size=max.chromosomes)
	if (all(abs(rowSums(Q.mat)) < 1e-8)==TRUE){
		print("Q matrix has all row sums equal to zero")
	}else{
		print("Q matrix not well defined. Rows don't add to zero")
		break
	}
	charnum=1
	nb.tip <- length(phy$tip.label)
	nb.node <- phy$Nnode
	nl <- nrow(Q.mat)
	#Now we need to build the matrix of likelihoods to pass to dev.raydisc:
	liks <- matrix(0, nb.tip + nb.node, nl)
	#Now loop through the tips.
	for(i in 1:nb.tip){
		#The codon at a site for a species is not NA, then just put a 1 in the appropriate column.
		#Note: We add charnum+1, because the first column in the data is the species labels:
		if(!is.na(tip.values[i,charnum+1])){
			liks[i,tip.values[i,charnum+1]] <- 1
		}else{
			#If here, then the site has no data, so we treat it as ambiguous for all possible codons. Likely things might be more complicated, but this can modified later:
			liks[i,] <- 1
		}
	}
	#The result here is just the likelihood:
	result <- pruning_like(p=NULL, phy=phy, liks=liks, Q=Q.mat, rate=NULL, root.p=pi.0)
	return(result)
}
