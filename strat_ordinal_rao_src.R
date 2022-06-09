## Author: Yiling Huang, yilingh@umich.edu
## Note: This file provides the implementation of the function
## that calculates the Rao score test statistic for the adjacent-category
## logit model, under a stratified design of observations

library(tidyr)
library(Matrix)

#######################################################################
## The following are implementations for the program that calculates ##
## the first derivative (score) vector of the parameters	     ##
#######################################################################

# Requires: the labels in n_js are consecutive numbers 
#						such as 1,2,3; 1,2; 1,2,3,4; ...;
#           J is the total number of categories.
y_to_indicator = function(y_labels, J) {
	ts = matrix(data = 0, nrow = length(y_labels), ncol = J)
	# Convert ordinal labels into a (n x J) matrix of indicator variables
	for (i in 1:length(y_labels)) {
		ts[i, y_labels[i]] = 1
	}
	return(ts)
}

# Requires: the labels in strata_labels are consecutive numbers 
#  					such as 1,2,3; 1,2; 1,2,3,4; ...;
#           s is the total number of strata.
strata_to_indicator = function(strata_labels, s) {
	I = matrix(data = 0, nrow = length(strata_labels), ncol = s)
	# Convert strata into a (n x s) matrix of indicator variables
	for (i in 1:length(strata_labels)) {
		I[i, strata_labels[i]] = 1
	}
	return(I)
}

# Requires: y_indicators is the output of y_to_indicator()
#           j is the category of interest
# 					J is the total number of categories
z_j = function(y_indicators, j, J) {
	if (j == 1) {
		return(y_indicators[,1])
	} else {
		return(rowSums(y_indicators[,1:j]))
	}
}

# Requires: strata_indicators is the output of strata_to_indicator()
# 					s is the total number of starta
# Effect:   Calculates the matrix of strata-specific mean vectors, 
#           of dimension s x p
stratified_means = function(X, strata_indicators, s) {
	means = matrix(0, nrow = s, ncol = ncol(X))
	# Calculate strata-specific mean vectors
	for (i in 1:s) {
		means[i,] = colMeans(X[as.logical(strata_indicators[,i]),])
	}
	return(means)
}

# Requires: strata_labels is the vector indicating strata membership
#						y_labels is the label indicating category membership
# Effect: Calculate the entire score vector, of dimension p(J-1) x 1
stratified_score = function(X, strata_labels, y_labels) {
	s = length(unique(strata_labels))
	J = length(unique(y_labels))
	p = ncol(X)
	
	y_ind = y_to_indicator(y_labels = y_labels, J = J)
	strata_ind = strata_to_indicator(strata_labels = strata_labels, s = s)
	means = stratified_means(X = X, strata_indicators = strata_ind, s = s)
	
	score = matrix(0, nrow = p*(J-1))
	
	# Assign values by block
	for (j in 1:(J-1)) {
		# The j-th cumulative category count
		zj = z_j(y_indicators = y_ind, j = j)
		# The j-th score vector
		score[((j-1)*p + 1):(j*p),] = t(X - strata_ind %*% means) %*% zj
	}
	return(score)
}

#######################################################################
## The following are implementations for the program that calculates ##
## the second derivative (Hessian) matrix of the parameters	     ##
#######################################################################
# Requires: y_indicators is the output of y_to_indicator()
# Effect: 	Horizontally stack z_j's into a N x J matrix,
#           where the (i,j)-th entry indicates whether y_i <= j.
cumulative_counts = function(y_indicators) {
	return(t(apply(y_indicators,1,cumsum)))
}

# Requires: strata_labels is the vector indicating strata membership
#						y_labels is the label indicating category membership
# Effect: Calculate the entire Hessian, of dimension p(J-1) x p(J-1)
stratified_hessian = function(X, strata_labels, y_labels) {
	s = length(unique(strata_labels))
	J = length(unique(y_labels))
	p = ncol(X)
	
	y_ind = y_to_indicator(y_labels = y_labels, J = J)
	strata_ind = strata_to_indicator(strata_labels = strata_labels, s = s)
	cum_counts = cumulative_counts(y_ind)
	
	hessian = matrix(0, nrow = p*(J-1), ncol = p*(J-1))
	
	# Assign values by block
	for (b in 1:s) {
		# c_b: collection of c_jb, c_lb's.
		c_b = colSums(cum_counts[as.logical(strata_ind[,b]),])
		n_b = sum(strata_ind[,b])
		# A_b: matrix of multipliers
		A_b = matrix(0, nrow = J-1, ncol = J-1)
		for (j in 1:(J-1)) {
			for (l in j:(J-1)) {
				A_b[j,l] = c_b[j] * (c_b[l] - n_b) / n_b
				A_b[l,j] = A_b[j,l]
			}
		}
		# Covariance matrix at the b-th strata
		cov_b = cov(X[as.logical(strata_ind[,b]),])
		# Update
		hessian = hessian + kronecker(A_b, cov_b)
	}
	return(hessian)
}

calculateStatistic = function(score, Hessian, tol=1e-10) {
	Hessian_eigen = eigen(Hessian)
	V = Hessian_eigen$vectors
	lambda = Hessian_eigen$values
	lambda_inv = numeric(length = length(lambda))
	rank = 0
	for (i in 1:length(lambda)) {
		if (abs(lambda[i]) < tol) {
			lambda_inv[i] = 0
		} else {
			lambda_inv[i] = 1/lambda[i]
			rank = rank + 1
		}
	}
	
	T2 = - t(score) %*% V %*% diag(lambda_inv) %*% t(V) %*% score
	
	return(list(T_sq = T2, rank = rank))
}
