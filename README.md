# Spatial_Scales_2017
Code for spatial population dynamics used in Hart et al. 2017
################################################################################
#
# This is a spatially explicit simulation of N species competing in a 
# heterogeneous environment, numerical analysis of coexistence, and a simulated 
# sampling scheme to reproduce coexistence-area and species-area curves. Phase 1
# produces an equilibrium community, Phase 2 tests the coexistence of each 
# species in this community more formally by allowing it to invade against the 
# remaining residents, and Phase 3 repeatedly subsamples the equilibrium 
# community over larger extents and tests coexistence of the subsampled community.
#
# The model is a discrete space lottery model, where multiple 
# populations may occupy the same site in varying proportions. The model allows 
# dispersal of seeds based on a dispersal kernel, and competition from 
# neighboring sites according to a competition kernel. Sedentary adults do not 
# compete, but experience a  constant, density-independent rate of survival.   
#
# Currently set for short-range dispersal (and competition), and long range
# spatial correlation. 
# See section "Variables of model that can be tuned" to change these and other
# parameter values. 
