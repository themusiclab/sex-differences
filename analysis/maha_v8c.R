#
# maha_v8c.R
#
# Marco Del Giudice (2018). Version 8c. Contact: marcodg@unm.edu. 
#
# Two simple R functions to compute multivariate standardized differences between two groups (Mahalanobis' D), either from raw data with function maha() or from Cohen's d values and correlations with function maha.summary(). The functions return uncorrected and bias-corrected values of D, exact confidence intervals, bootstrapped confidence intervals (only from raw data), heterogeneity coefficients, and a number of additional diagnostics and effect sizes. Disattenuated estimates can also be computed if desired. For more information see Del Giudice (2009, 2013, 2019).
#
# Confidence intervals: exact confidence intervals for D are computed with Reiser's (2001) method. Especially for small D values, the equations may not be solvable; in those cases, one or both CI bounds are set to NA. For more information see Reiser (2001). The bootstrapped CIs are bias-corrected and accelerated; for details see Kelley (2005).
#
# Additional effect sizes: additional effect sizes are calculated from D under the assumptions of multivariate normality and equality of covariance matrices. These are: the overlapping coefficient OVL; Cohen’s alternative coefficient of overlap OVL2 calculated on the joint distribution (equivalent to 1–U1); the common language effect size CL (McGraw & Wong, 1992); and the probability of correct classification (PCC) based on linear discriminant analysis (assuming equal group sizes). For details and discussion see Del Giudice (2019).
#
# Heterogeneity coefficients: the heterogeneity coefficient H2 quantifies heterogeneity in the variables' contribution to the multivariate effect size (0 = max homogeneity; 1 = max heterogeneity). Coefficient EPV2 represents heterogeneity as the proportion of contributing variables that would result in the same heterogeneity (equivalent proportion of variables), in a hypothetical scenario where a certain proportion of variables contribute equally to D while the remaining ones make no contribution (range: 0-1). See Del Giudice (2017, 2018).
#
# Similarity between correlation matrices: before computing D and other effect sizes, the correlation matrices of the two groups are pooled by taking weighted averages, unless a common correlation matrix is directly provided to function maha.summary(). If two matrices are provided to function maha.summary() without specifying the size of the groups, the unweighted average is taken. Tucker’s Congruence Coefficient (CC) can be used as an index of similarity between the two sample correlation matrices (–1 = max dissimilarity; 1 = max similarity). Values above .85 indicate fair similarity; values above .95 indicate high similarity (Lorenzo-Seva & ten Berge, 2006). The significance of the difference between correlation matrices can be tested with Steiger's test of equality (Steiger, 1980); however, significance tests can be overly sensitive to small differences when sample size is large. See Del Giudice (2019) for more details.
#
# Bias correction: the bias-corrected distance Du removes the upward bias in D, which is especially large when sample size is small relative to the number of variables (see Hess et al., 2007; Del Giudice, 2019). The bias-corrected version of Cohen’s d is du (also known as Hedges’ g); the approximate formula is used here (see Kelley, 2005).
#
# Disattenuation: disattenuated estimates of D and other effect sizes can be obtained by supplying a vector of reliability coefficients (e.g., Cronbach’s α). For details see Del Giudice (2019). 
#
#
# Del Giudice, M. (2009). On the real magnitude of psychological sex differences. Evolutionary Psychology, 7, 264-279. doi:10.1177/147470490900700209
# Del Giudice, M. (2013). Multivariate misgivings: Is D a valid measure of group and sex differences? Evolutionary Psychology, 11, 1067-1076. doi:10.1177/147470491301100511
# Del Giudice, M. (2017). Heterogeneity coefficients for Mahalanobis' D as a multivariate effect size. Multivariate Behavioral Research, 52, 216-221. doi:10.1080/00273171.2016.1262237
# Del Giudice, M. (2018). Addendum to: Heterogeneity coefficients for Mahalanobis' D as a multivariate effect size. Multivariate Behavioral Research, 53, 571-573. doi: 10.1080/00273171.2018.1462138
# Del Giudice, M. (2019). Measuring sex differences and similarities. In D. P. VanderLaan & W. I. Wong (Eds.), Gender and sexuality development: Contemporary theory and research. New York, NY: Springer.
# Hess, M. R., Hogarty, K. Y., Ferron, J. M., & Kromrey, J. D. (2007). Interval estimates of multivariate effect sizes: Coverage and interval width estimates under variance heterogeneity and nonnormality. Educational and Psychological Measurement, 67, 21-40. doi:10.1177/0013164406288159
# Kelley, K. (2005). The effects of nonnormal distributions on confidence intervals around the standardized mean difference: Bootstrap and parametric confidence intervals. Educational and Psychological Measurement, 65, 51-69. doi:10.1177/0013164404264850
# Lorenzo-Seva, U., & Ten Berge, J. M. (2006). Tucker's congruence coefficient as a meaningful index of factor similarity. Methodology, 2, 57-64. doi:10.1027/1614-2241.2.2.57
# McGraw, K. O., & Wong, S. P. (1992). A common language effect size statistic. Psychological Bulletin, 111, 361-365.
# Reiser, B. (2001). Confidence intervals for the Mahalanobis distance. Communications in Statistics: Simulation and Computation, 30, 37–45. doi:10.1081/SAC-100001856
# Steiger, J. H. (1980). Testing pattern hypotheses on correlation matrices: Alternative statistics and some empirical results. Multivariate Behavioral Research, 15(3), 335-352. doi:10.1207/s15327906mbr1503_7



# function maha(dataA, dataB, rel_values=NULL, conf.level=.95, boot.n=NULL, round.digits=3)
#
# Computes the Mahalanobis distance D, confidence intervals, and various relevant indices from raw data. Can compute disattenuated estimates if desired. 
#
# Arguments
#
# dataA, dataB 		raw matrices/data frames for the two groups
# rel_values		vector of reliability coefficients (optional: required for disattenuation)
# conf.level 		CI width (default is 95%)
# boot.n 			number of bootstrap samples (optional; recommended: 5,000-10,000)
# round.digits		rounding digits for the output (default is 3)
#
# Value
# returns a list object containing some or all of the following:
#
# D 				Mahalanobis' D
# CI_exact_D		exact CI for D (lower bound, upper bound); NA if not solvable
# CI_boot_D 		bootstrapped CI for D (lower bound, upper bound)
# OVL 				overlapping coefficient OVL (single distribution)
# OVL2 				coefficient of overlap OVL2 (joint distribution; equals 1-U1) 
# CL				Common language effect size CL
# PCC 				Probability of correct classification (PCC; linear discriminant analysis)
# H2 				heterogeneity coefficient H2
# EPV2 				EPV2 coefficient (equivalent proportion of variables)
# CC_cor 			Tucker's Congruence Coefficient CC (similarity between correlation matrices)
# steiger.p			p value for Steiger's test of equality of correlation matrices (null hypothesis: equality)
# d_values 			vector of Cohen's d values
#
# Du				bias-corrected Du
# CI_boot_Du 		bootstrapped CI for Du (lower bound, upper bound)
# OVLu				overlapping coefficient OVL based on Du
# OVL2u 			coefficient of overlap OVL2 based on Du
# CLu 				Common language effect size CL based on Du
# PCCu 				Probability of correct classification based on Du
# du_values 		vector of bias-corrected du values (equivalent to Hedges’ g)
#
# Dc 				disattenuated Dc
# OVLc 				disattenuated overlapping coefficient (single distribution)
# OVL2c 			disattenuated Cohen's coefficient of overlap (joint distribution)
# CLc 				disattenuated common language effect size 
# PCCc 				disattenuated probability of correct classification
# H2c 				disattenuated heterogeneity coefficient H2
# EPV2c 			disattenuated EPV2 coefficient
# dc_values 		vector of disattenuated Cohen's d values
#
# Dcu 				disattenuated, bias-corrected Dcu
# OVLcu 			disattenuated overlapping coefficient (single distribution) based on Dcu
# OVL2cu 			disattenuated Cohen’s coefficient of overlap (joint distribution) based on Dcu
# CLcu 				disattenuated common language effect size based on Dcu
# PCCcu 			disattenuated probability of correct classification based on Dcu
# dcu_values 		vector of disattenuated, bias-corrected dcu values (disattenuated Hedges’ g)




# function maha.summary(d_values, corA, corB=NULL, nA=NULL, nB=NULL, rel_values=NULL, conf.level=.95, round.digits=3)
#
# Returns the Mahalanobis distance D, the bias-corrected distance Du, confidence intervals, heterogeneity coefficients, and coefficients of overlap, computed from summary statistics. Can compute disattenuated estimates (Dc, Dcu) if desired.
#
# Arguments
#
# d_values 			(row) vector of standardized differences (Cohen's d)
# corA 				either the common correlation matrix, or the correlation matrix for group A 
# corB 				correlation matrix for group B (optional) 
# alpha 			vector of reliability coefficients (optional: only for disattenuation)
# nA, nB 			sample size of the two groups (optional: required for bias correction, exact CIs, and Steiger's test)
# conf.level 		CI width (optional)
# round.digits		rounding digits for the output (default is 3)
#
# Value
# returns a list object containing some or all of the following:
#
# D 				Mahalanobis' D
# CI_exact_D		exact CI for D (lower bound, upper bound); NA if not solvable
# OVL 				overlapping coefficient OVL (single distribution)
# OVL2 				coefficient of overlap OVL2 (joint distribution; equals 1-U1) 
# CL				Common language effect size CL
# PCC 				Probability of correct classification (PCC; linear discriminant analysis)
# H2 				heterogeneity coefficient H2
# EPV2 				EPV2 coefficient (equivalent proportion of variables)
# CC_cor 			Tucker's Congruence Coefficient CC (similarity between correlation matrices)
# steiger.p			p value for Steiger's test of equality of correlation matrices (null hypothesis: equality)
# d_values 			vector of Cohen's d values
#
# Du				bias-corrected Du
# OVLu				overlapping coefficient OVL based on Du
# OVL2u 			coefficient of overlap OVL2 based on Du
# CLu 				Common language effect size CL based on Du
# PCCu 				Probability of correct classification based on Du
# du_values 		vector of bias-corrected du values (equivalent to Hedges’ g)
#
# Dc 				disattenuated Dc
# OVLc 				disattenuated overlapping coefficient (single distribution)
# OVL2c 			disattenuated Cohen's coefficient of overlap (joint distribution)
# CLc 				disattenuated common language effect size 
# PCCc 				disattenuated probability of correct classification
# H2c 				disattenuated heterogeneity coefficient H2
# EPV2c 			disattenuated EPV2 coefficient
# dc_values 		vector of disattenuated Cohen's d values
#
# Dcu 				disattenuated, bias-corrected Dcu
# OVLcu 			disattenuated overlapping coefficient (single distribution) based on Dcu
# OVL2cu 			disattenuated Cohen’s coefficient of overlap (joint distribution) based on Dcu
# CLcu 				disattenuated common language effect size based on Dcu
# PCCcu 			disattenuated probability of correct classification based on Dcu
# dcu_values 		vector of disattenuated, bias-corrected dcu values (disattenuated Hedges’ g)







############ START Function maha()

maha <- function(dataA, dataB, rel_values=NULL, conf.level=.95, boot.n=NULL, round.digits=3) {

############ preliminary computations: Ns, pooled variances, pooled correlation matrix

nA = as.numeric(length(complete.cases(dataA)))
nB = as.numeric(length(complete.cases(dataB)))

corA = cor(dataA, use="pairwise.complete.obs")
corB = cor(dataB, use="pairwise.complete.obs")
pooled_cor = (corA*nA+corB*nB)/(nA+nB)

pooled_variances = ((nA-1)*diag(var(dataA, na.rm=TRUE))+(nB-1)*diag(var(dataB, na.rm=TRUE)))/(nA+nB-2)

############ compute Mahalanobis D (uncorrected)

d_values = (colMeans(dataA, na.rm=TRUE)-colMeans(dataB, na.rm=TRUE))/sqrt(pooled_variances)

p = length(d_values)

D2 = mahalanobis(d_values, pooled_cor, center=FALSE)
output = list()
output$D = sqrt(D2)

############ compute exact CI (uncorrected)

Fcal = (D2)*(nA*nB*(nA+nB-p-1))/(p*(nA+nB)*(nA+nB-2))

lower_prob = conf.level+(1-conf.level)/2
upper_prob = (1-conf.level)/2
critical_F_upper = qf(upper_prob, p, nA+nB-p-1)
critical_F_lower = qf(lower_prob, p, nA+nB-p-1)

# lower
if (D2>critical_F_lower) { 
	ncp_est_max=10000/(1/nA+1/nB)
	ncp_est_min = 0
	ncp_est = (ncp_est_max-ncp_est_min)/2

	est_p = pf(Fcal, p, (nA+nB-p-1), lower.tail=TRUE, ncp=ncp_est)

	while(abs(est_p-lower_prob) > .00001) {
		if ((est_p-lower_prob) < 0) {ncp_est_max=ncp_est
			ncp_est = ncp_est_min+(ncp_est_max-ncp_est_min)/2
			}
		else {ncp_est_min = ncp_est
			ncp_est = ncp_est_min+(ncp_est_max-ncp_est_min)/2
			}
		est_p = pf(Fcal, p, (nA+nB-p-1), lower.tail=TRUE, ncp=ncp_est);
		}
	lower = sqrt(ncp_est*(1/nA+1/nB))
}
else lower=NA

# upper
if (D2>critical_F_upper) { 
	ncp_est_max=10000/(1/nA+1/nB)
	ncp_est_min = 0
	ncp_est = (ncp_est_max-ncp_est_min)/2
	
	est_p = pf(Fcal, p, (nA+nB-p-1), lower.tail=TRUE, ncp=ncp_est)
	
	while(abs(est_p-upper_prob) > .00001) {
		if ((est_p-upper_prob) < 0) {ncp_est_max = ncp_est
			ncp_est=ncp_est_min+(ncp_est_max-ncp_est_min)/2
			}
		else {ncp_est_min = ncp_est
			ncp_est = ncp_est_min+(ncp_est_max-ncp_est_min)/2
			}
		est_p = pf(Fcal, p, (nA+nB-p-1), lower.tail=TRUE, ncp=ncp_est);
		}
	output$CI_exact_D = c(lower, sqrt(ncp_est*(1/nA+1/nB)) )
	}
else output$CI_exact_D = c(lower, NA)

############ compute bootstrapped CI (uncorrected)

if (is.null(boot.n)==FALSE) {
	boot_D = numeric(boot.n)
	alpha_level = 1-conf.level
	
	for (sample in 1:boot.n) {
		data_A.sampled = dataA[sample(seq(1:nrow(dataA)), replace = TRUE),]
		data_B.sampled = dataB[sample(seq(1:nrow(dataB)), replace = TRUE),]
		
		nA.sampled = length(complete.cases(data_A.sampled))
		nB.sampled = length(complete.cases(data_B.sampled))
		cor_A.sampled = cor(data_A.sampled, use="pairwise.complete.obs")
		cor_B.sampled = cor(data_B.sampled, use="pairwise.complete.obs")
		pooled_cor.sampled = (cor_A.sampled*nA.sampled+cor_B.sampled*nB.sampled)/(nA.sampled+nB.sampled)
		pooled_variances.sampled = ((nA.sampled-1)*diag(var(data_A.sampled, na.rm=TRUE))+(nB.sampled-1)*diag(var(data_B.sampled, na.rm=TRUE)))/(nA.sampled+nB.sampled-2)
		d_values.sampled = (colMeans(data_A.sampled, na.rm=TRUE)-colMeans(data_B.sampled, na.rm=TRUE))/sqrt(pooled_variances.sampled)
		boot_D[sample] = sqrt(mahalanobis(d_values.sampled, pooled_cor.sampled, center=FALSE))
	}
	
	jackknife_results = numeric(nrow(dataA)+nrow(dataB))
	n.1 = nrow(dataA)
	n.2 = nrow(dataB)
	
	Marker.1 = seq(1, n.1, 1)
	
	for (sample in 1:n.1) {
		data_A.jack = dataA[Marker.1[-sample],]
		data_B.jack = dataB
		
		nA.jack = length(complete.cases(data_A.jack))
		nB.jack = length(complete.cases(data_B.jack))
		cor_A.jack = cor(data_A.jack, use="pairwise.complete.obs")
		cor_B.jack = cor(data_B.jack, use="pairwise.complete.obs")
		pooled_cor.jack = (cor_A.jack*nA.jack+cor_B.jack*nB.jack)/(nA.jack+nB.jack)
		pooled_variances.jack = ((nA.jack-1)*diag(var(data_A.jack, na.rm=TRUE))+(nB.jack-1)*diag(var(data_B.jack, na.rm=TRUE)))/(nA.jack+nB.jack-2)
		d_values.jack = (colMeans(data_A.jack, na.rm=TRUE)-colMeans(data_B.jack, na.rm=TRUE))/sqrt(pooled_variances.jack)	
		jackknife_results[sample] = sqrt(mahalanobis(d_values.jack, pooled_cor.jack, center=FALSE))
	}
	
	Marker.2 = seq(1, n.2, 1)
	
	for (sample in 1:n.2) {
		data_A.jack = dataA
		data_B.jack = dataB[Marker.2[-sample],]
		
		nA.jack = length(complete.cases(data_A.jack))
		nB.jack = length(complete.cases(data_B.jack))
		cor_A.jack = cor(data_A.jack, use="pairwise.complete.obs")
		cor_B.jack = cor(data_B.jack, use="pairwise.complete.obs")
		pooled_cor.jack = (cor_A.jack*nA.jack+cor_B.jack*nB.jack)/(nA.jack+nB.jack)
		pooled_variances.jack = ((nA.jack-1)*diag(var(data_A.jack, na.rm=TRUE))+(nB.jack-1)*diag(var(data_B.jack, na.rm=TRUE)))/(nA.jack+nB.jack-2)
		d_values.jack = (colMeans(data_A.jack, na.rm=T)-colMeans(data_B.jack, na.rm=TRUE))/sqrt(pooled_variances.jack)
		jackknife_results[n.1+sample] = sqrt(mahalanobis(d_values.jack, pooled_cor.jack, center=FALSE))
	}
	
	Mean.Jackknife = mean(jackknife_results)
	
	a = (sum((Mean.Jackknife-jackknife_results)^3))/(6*sum((Mean.Jackknife-jackknife_results)^2)^(3/2))
	
	z0 = qnorm(sum(boot_D < output$D)/boot.n)
	
	CI.Low.BCa = pnorm(z0 + (z0+qnorm(alpha_level/2))/(1-a*(z0+qnorm(alpha_level/2))))
	CI.Up.BCa = pnorm(z0 + (z0+qnorm(1-alpha_level/2))/(1-a*(z0+qnorm(1-alpha_level/2))))
	
	output$CI_boot_D = c(quantile(boot_D, CI.Low.BCa, names=FALSE), quantile(boot_D, CI.Up.BCa, names=FALSE))
}

############ compute additional effect sizes (uncorrected)

output$OVL = 2*pnorm(-output$D/2)
output$OVL2 = output$OVL/(2-output$OVL)
output$CL = pnorm(output$D/sqrt(2))
output$PCC = pnorm(output$D/2)

############ compute heterogeneity coefficients H2, EPV2

C_values = (d_values%*%solve(pooled_cor))*d_values
C_values = sort(abs(C_values))
N = length(C_values)
sum_Ci = sum(C_values)
sum_iCi = sum(C_values*seq(1:N))

output$H2 = ( (2/N)*sum_iCi -((N+1)/N)*sum_Ci) / ((N-1)*mean(C_values))
output$EPV2 = 1-(output$H2*(N-1)/N)

############ compute Tucker's CC, Steiger test, and Cohen's d values

cor_ab = corA*corB
cor_a2 = corA^2
cor_b2 = corB^2

output$CC_cor = sum(cor_ab[lower.tri(cor_ab)])/sqrt(sum(cor_a2[lower.tri(cor_a2)])*sum(cor_b2[lower.tri(cor_b2)]))

z_diff_matrix = 0.5*log((1+corA)/(1-corA)) - 0.5*log((1+corB)/(1-corB))   # z = Fisher's transformation
df = (p^2-p)/2
chi2 = ((nA*nB)/(nA+nB)-3)*sum(z_diff_matrix[lower.tri(z_diff_matrix)]^2)
p_value = pchisq(chi2, df, lower.tail=FALSE)

output$steiger.p = p_value

output$d_values=d_values

############ compute bias-corrected Du

output$Du = sqrt(max(0,(D2*(nA+nB-p-3)/(nA+nB-2)-p*(nA+nB)/(nA*nB))))
		
############ compute bootstrapped CI based on Du (bias-corrected)

if (is.null(boot.n)==FALSE) {
	boot_D = numeric(boot.n)
	alpha_level = 1-conf.level
	
	for (sample in 1:boot.n) {
		data_A.sampled = dataA[sample(seq(1:nrow(dataA)), replace = TRUE),]
		data_B.sampled = dataB[sample(seq(1:nrow(dataB)), replace = TRUE),]
		
		nA.sampled = length(complete.cases(data_A.sampled))
		nB.sampled = length(complete.cases(data_B.sampled))
		cor_A.sampled = cor(data_A.sampled, use="pairwise.complete.obs")
		cor_B.sampled = cor(data_B.sampled, use="pairwise.complete.obs")
		pooled_cor.sampled = (cor_A.sampled*nA.sampled+cor_B.sampled*nB.sampled)/(nA.sampled+nB.sampled)
		pooled_variances.sampled = ((nA.sampled-1)*diag(var(data_A.sampled, na.rm=TRUE))+(nB.sampled-1)*diag(var(data_B.sampled, na.rm=TRUE)))/(nA.sampled+nB.sampled-2)
		d_values.sampled = (colMeans(data_A.sampled, na.rm=TRUE)-colMeans(data_B.sampled, na.rm=TRUE))/sqrt(pooled_variances.sampled)
		unbiased_D2 = mahalanobis(d_values.sampled, pooled_cor.sampled, center=FALSE)*(nA+nB-p-3)/(nA+nB-2)-p*(nA+nB)/(nA*nB)
		boot_D[sample] = sqrt(max(0, unbiased_D2))
	}
	
	jackknife_results = numeric(nrow(dataA)+nrow(dataB))
	n.1 = nrow(dataA)
	n.2 = nrow(dataB)
	
	Marker.1 = seq(1, n.1, 1)
	
	for (sample in 1:n.1) {
		data_A.jack = dataA[Marker.1[-sample],]
		data_B.jack = dataB
		
		nA.jack = length(complete.cases(data_A.jack))
		nB.jack = length(complete.cases(data_B.jack))
		cor_A.jack = cor(data_A.jack, use="pairwise.complete.obs")
		cor_B.jack = cor(data_B.jack, use="pairwise.complete.obs")
		pooled_cor.jack = (cor_A.jack*nA.jack+cor_B.jack*nB.jack)/(nA.jack+nB.jack)
		pooled_variances.jack = ((nA.jack-1)*diag(var(data_A.jack, na.rm=TRUE))+(nB.jack-1)*diag(var(data_B.jack, na.rm=TRUE)))/(nA.jack+nB.jack-2)
		d_values.jack = (colMeans(data_A.jack, na.rm=TRUE)-colMeans(data_B.jack, na.rm=TRUE))/sqrt(pooled_variances.jack)	
		unbiased_D2 = mahalanobis(d_values.jack, pooled_cor.jack, center=FALSE)*(nA+nB-p-3)/(nA+nB-2)-p*(nA+nB)/(nA*nB)
		jackknife_results[sample] = sqrt(max(0, unbiased_D2))
	}
	
	Marker.2 = seq(1, n.2, 1)
	
	for (sample in 1:n.2) {
		data_A.jack = dataA
		data_B.jack = dataB[Marker.2[-sample],]
		
		nA.jack = length(complete.cases(data_A.jack))
		nB.jack = length(complete.cases(data_B.jack))
		cor_A.jack = cor(data_A.jack, use="pairwise.complete.obs")
		cor_B.jack = cor(data_B.jack, use="pairwise.complete.obs")
		pooled_cor.jack = (cor_A.jack*nA.jack+cor_B.jack*nB.jack)/(nA.jack+nB.jack)
		pooled_variances.jack = ((nA.jack-1)*diag(var(data_A.jack, na.rm=TRUE))+(nB.jack-1)*diag(var(data_B.jack, na.rm=TRUE)))/(nA.jack+nB.jack-2)
		d_values.jack = (colMeans(data_A.jack, na.rm=T)-colMeans(data_B.jack, na.rm=TRUE))/sqrt(pooled_variances.jack)
		unbiased_D2 = mahalanobis(d_values.jack, pooled_cor.jack, center=FALSE)*(nA+nB-p-3)/(nA+nB-2)-p*(nA+nB)/(nA*nB)
		jackknife_results[n.1+sample] = sqrt(max(0, unbiased_D2))
	}
	
	Mean.Jackknife = mean(jackknife_results)
	
	a = (sum((Mean.Jackknife-jackknife_results)^3))/(6*sum((Mean.Jackknife-jackknife_results)^2)^(3/2))
	
	z0 = qnorm(sum(boot_D < output$Du)/boot.n)
	
	CI.Low.BCa = pnorm(z0 + (z0+qnorm(alpha_level/2))/(1-a*(z0+qnorm(alpha_level/2))))
	CI.Up.BCa = pnorm(z0 + (z0+qnorm(1-alpha_level/2))/(1-a*(z0+qnorm(1-alpha_level/2))))
	
	output$CI_boot_Du = c(quantile(boot_D, CI.Low.BCa, names=FALSE), quantile(boot_D, CI.Up.BCa, names=FALSE))
}

############ compute additional effect sizes based on Du (bias-corrected)

output$OVLu = 2*pnorm(-output$Du/2)
output$OVL2u = output$OVLu/(2-output$OVLu)
output$CLu = pnorm(output$Du/sqrt(2))
output$PCCu = pnorm(output$Du/2)

############ compute du / Hedges' g values

output$du_values=d_values*(1-3/(4*(nA+nB-2)-1))

############ compute disattenuated Mahalanobis D (uncorrected)

if (is.null(rel_values)==FALSE) {
	d_values.d = d_values/sqrt(rel_values)
	rel_matrix = sqrt(crossprod(t(rel_values), t(rel_values)))
	diag(rel_matrix) = 1
	pooled_cor.d = pooled_cor/rel_matrix
	D2.d = mahalanobis(d_values.d, pooled_cor.d, center=FALSE)
	output$Dc = sqrt(D2.d)
}

############ compute disattenuated effect sizes

if (is.null(rel_values)==FALSE) {
	output$OVLc = 2*pnorm(-output$Dc/2)
	output$OVL2c = output$OVLc/(2-output$OVLc)
	output$CLc = pnorm(output$Dc/sqrt(2))
	output$PCCc = pnorm(output$Dc/2)
}

############ compute disattenuated heterogeneity coefficients H2, EPV2

if (is.null(rel_values)==FALSE) {
	C_values = (d_values.d%*%solve(pooled_cor.d))*d_values.d
	C_values = sort(abs(C_values))
	N = length(C_values)
	sum_Ci = sum(C_values)
	sum_iCi = sum(C_values*seq(1:N))
	
	output$H2c = ( (2/N)*sum_iCi -((N+1)/N)*sum_Ci) / ((N-1)*mean(C_values))
	output$EPV2c = 1-(output$H2c*(N-1)/N)	
}

############ compute disattenuated Cohen's d values

if (is.null(rel_values)==FALSE) {
	output$dc_values = d_values.d
}

############ compute disattenuated, bias-corrected Dcu

if (is.null(rel_values)==FALSE) {
	output$Dcu = sqrt(max(0,(D2.d*(nA+nB-p-3)/(nA+nB-2)-p*(nA+nB)/(nA*nB))))
}

############ compute disattenuated effect sizes based on Dcu and dcu / disattenuated Hedges' g values

if (is.null(rel_values)==FALSE) {
	output$OVLcu = 2*pnorm(-output$Dcu/2)
	output$OVL2cu = output$OVLcu/(2-output$OVLcu)
	output$CLcu = pnorm(output$Dcu/sqrt(2))
	output$PCCcu = pnorm(output$Dcu/2)
	output$dcu_values = d_values.d*(1-3/(4*(nA+nB-2)-1))
}

############ output

return(lapply(output, round, round.digits))
}

############ END Function maha()







############ START Function maha.summary()

maha.summary <- function(d_values, corA, corB=NULL, nA=NULL, nB=NULL, rel_values=NULL, conf.level=.95, round.digits=3) {

############ preliminary computations: pooled correlation matrix

if (is.null(nA)==FALSE & is.null(corB)==FALSE)	{
	pooled_cor = (corA*nA+corB*nB)/(nA+nB)
}

if (is.null(nA)==TRUE & is.null(corB)==FALSE) {
	pooled_cor = (corA+corB)/(2)
}

if (is.null(corB)==TRUE) {
	pooled_cor = corA
}

############ compute Mahalanobis D

D2 = mahalanobis(d_values, pooled_cor, center=FALSE)
output = list()
output$D = sqrt(D2)
p = length(d_values)		

############ compute exact CI

if (is.null(nA)==FALSE)	{
	Fcal = (D2)*(nA*nB*(nA+nB-p-1))/(p*(nA+nB)*(nA+nB-2))
	lower_prob = conf.level+(1-conf.level)/2
	upper_prob = (1-conf.level)/2
	critical_F_upper = qf(upper_prob, p, nA+nB-p-1)
	critical_F_lower = qf(lower_prob, p, nA+nB-p-1)

	# lower
	if (D2>critical_F_lower) { 
		ncp_est_max = 10000/(1/nA+1/nB)
		ncp_est_min = 0
		ncp_est = (ncp_est_max-ncp_est_min)/2
		
		est_p = pf(Fcal, p, (nA+nB-p-1), lower.tail=TRUE, ncp=ncp_est)
		
		while(abs(est_p-lower_prob) > .00001) {
			if ((est_p-lower_prob) < 0) {ncp_est_max = ncp_est
				ncp_est = ncp_est_min+(ncp_est_max-ncp_est_min)/2
				}
			else {ncp_est_min = ncp_est
				ncp_est = ncp_est_min+(ncp_est_max-ncp_est_min)/2
				}
			est_p = pf(Fcal, p, (nA+nB-p-1), lower.tail=TRUE, ncp=ncp_est);
			}
		lower = sqrt(ncp_est*(1/nA+1/nB))
	}
	else lower = NA
	
	# upper
	if (D2>critical_F_upper) { 
		ncp_est_max = 10000/(1/nA+1/nB)
		ncp_est_min = 0
		ncp_est = (ncp_est_max-ncp_est_min)/2
		
		est_p = pf(Fcal, p, (nA+nB-p-1), lower.tail=TRUE, ncp=ncp_est)
		
		while(abs(est_p-upper_prob) > .00001 ) {
			if ((est_p-upper_prob) < 0) {ncp_est_max = ncp_est
				ncp_est = ncp_est_min+(ncp_est_max-ncp_est_min)/2
				}
			else {ncp_est_min = ncp_est
				ncp_est = ncp_est_min+(ncp_est_max-ncp_est_min)/2
				}
			est_p = pf(Fcal, p, (nA+nB-p-1), lower.tail=TRUE, ncp=ncp_est);
			}
		output$CI_exact_D = c(lower, sqrt(ncp_est*(1/nA+1/nB)) )
	}
	else output$CI_exact_D = c(lower, NA)
}

############ compute additional effect sizes (uncorrected)

output$OVL = 2*pnorm(-output$D/2)
output$OVL2 = output$OVL/(2-output$OVL)
output$CL = pnorm(output$D/sqrt(2))
output$PCC = pnorm(output$D/2)

############ compute heterogeneity coefficients H2, EPV2 and Cohen's d values

C_values = (d_values%*%solve(pooled_cor))*d_values
C_values = sort(abs(C_values))
N = length(C_values)
sum_Ci = sum(C_values)
sum_iCi = sum(C_values*seq(1:N))

output$H2 = ( (2/N)*sum_iCi -((N+1)/N)*sum_Ci) / ((N-1)*mean(C_values))
output$EPV2 = 1-(output$H2*(N-1)/N)

############ compute Tucker's CC, Steiger test and Cohen's d values

if (is.null(corB)==FALSE) {
	cor_ab = corA*corB
	cor_a2 = corA^2
	cor_b2 = corB^2

	output$CC_cor = sum(cor_ab[lower.tri(cor_ab)])/sqrt(sum(cor_a2[lower.tri(cor_a2)])*sum(cor_b2[lower.tri(cor_b2)]))
}

if (is.null(nA)==FALSE & is.null(corB)==FALSE)	{
	z_diff_matrix = 0.5*log((1+corA)/(1-corA)) - 0.5*log((1+corB)/(1-corB))   # z = Fisher's transformation
	df = (p^2-p)/2
	chi2 = ((nA*nB)/(nA+nB)-3)*sum(z_diff_matrix[lower.tri(z_diff_matrix)]^2)
	p_value = pchisq(chi2, df, lower.tail=FALSE)

	output$steiger.p = p_value
}

output$d_values=d_values

############ compute bias-corrected Du

if (is.null(nA)==FALSE)	{
	output$Du = sqrt(max(0,(D2*(nA+nB-p-3)/(nA+nB-2)-p*(nA+nB)/(nA*nB))))
}	

############ compute additional effect sizes based on Du (bias-corrected)

if (is.null(nA)==FALSE)	{
	output$OVLu = 2*pnorm(-output$Du/2)
	output$OVL2u = output$OVLu/(2-output$OVLu)
	output$CLu = pnorm(output$Du/sqrt(2))
	output$PCCu = pnorm(output$Du/2)
}

############ compute du / Hedges' g values

if (is.null(nA)==FALSE)	{
	output$du_values=d_values*(1-3/(4*(nA+nB-2)-1))
}

############ compute disattenuated Mahalanobis D (uncorrected)

if (is.null(rel_values)==FALSE) {
	d_values.d = d_values/sqrt(rel_values)
	rel_matrix = sqrt(crossprod(t(rel_values), t(rel_values)))
	diag(rel_matrix) = 1
	pooled_cor.d = pooled_cor/rel_matrix
	D2.d = mahalanobis(d_values.d, pooled_cor.d, center=FALSE)
	output$Dc = sqrt(D2.d)
}

############ compute disattenuated effect sizes

if (is.null(rel_values)==FALSE) {
	output$OVLc = 2*pnorm(-output$Dc/2)
	output$OVL2c = output$OVLc/(2-output$OVLc)
	output$CLc = pnorm(output$Dc/sqrt(2))
	output$PCCc = pnorm(output$Dc/2)
}

############ compute disattenuated heterogeneity coefficients H2, EPV2

if (is.null(rel_values)==FALSE) {
	C_values = (d_values.d%*%solve(pooled_cor.d))*d_values.d
	C_values = sort(abs(C_values))
	N = length(C_values)
	sum_Ci = sum(C_values)
	sum_iCi = sum(C_values*seq(1:N))
	
	output$H2c = ( (2/N)*sum_iCi -((N+1)/N)*sum_Ci) / ((N-1)*mean(C_values))
	output$EPV2c = 1-(output$H2c*(N-1)/N)	
}

############ compute disattenuated Cohen's d values

if (is.null(rel_values)==FALSE) {
	output$dc_values = d_values.d
}

############ compute disattenuated, bias-corrected Dcu

if (is.null(rel_values)==FALSE & is.null(nA)==FALSE) {
	output$Dcu = sqrt(max(0,(D2.d*(nA+nB-p-3)/(nA+nB-2)-p*(nA+nB)/(nA*nB))))
}

############ compute disattenuated effect sizes based on Dcu and dcu / disattenuated Hedges' g values

if (is.null(rel_values)==FALSE & is.null(nA)==FALSE) {
	output$OVLcu = 2*pnorm(-output$Dcu/2)
	output$OVL2cu = output$OVLcu/(2-output$OVLcu)
	output$CLcu = pnorm(output$Dcu/sqrt(2))
	output$PCCcu = pnorm(output$Dcu/2)
	output$dcu_values = d_values.d*(1-3/(4*(nA+nB-2)-1))
}
	
############ output

return(lapply(output, round, round.digits))
}

############ END Function maha.summary()
