#=======================================================
#-------------------------------------------------------
#               Random Forests to prioritize/reduce
#               Climate variables
#               Leander DL Anderegg
#                   12/07/15
#-------------------------------------------------------
#=======================================================

## the goal of this script is to use Random Forest regression algorithms to identify the most
# important climate variables, either on a per mountain or per species basis (possibly both in the same analysis?)
# From some online reading (e.g. http://www.r-bloggers.com/a-brief-tour-of-the-trees-and-forests/)
# it seems like it might even be possible to throw down random forest regression with random effects 
# and different error structures from {nlme}

# last updated: 12/07/15


##############################

