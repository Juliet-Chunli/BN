#!/usr/bin/python3.6
import math as ma

def varToMedian(mean, var):
	median = mean / ma.sqrt(var / mean**2 + 1)
	return(median)

def medianToVar(mean, median):
	var = ((mean / median)**2 - 1) * mean ** 2
	return var

def meanMedianToMuStd(mean, median):
	mu = ma.log(median)
	std = (2*(ma.log(mean) - ma.log(median)))**0.5
	return (mu, std)

def meanVarToMuStd(mean, var):
	std = ma.sqrt(ma.log(1 + var / mean**2))
	mu = ma.log(mean) - ma.log(1 + var / mean**2) / 2
	return (mu, std)



