
## Piecewise regression around a guessed break point "psi" 
# from bill 1/23/14

dir="C:/Users/Bill Anderegg/Desktop/Projects/Hindcasting SAD/Data/0---Model runs---0/"
x7=read.csv(file=paste(dir, "Ks-CWD-reg.csv", sep=""),header=TRUE)

#------------  PieceWise regression

install.packages("segmented")
library(segmented)

#Combined 2010-2011 measurements
CWD <- x7$CWD[c(1:4,6:16)]        # X variable
Ks <- x7$Ks[c(1:4,6:16)]*100            # Y Variable
modelw <- lm(Ks ~ CWD)                #Initial linear model that you will then "break" into pieces
segmented.mod <- segmented(modelw, seg.Z = ~CWD, psi=410)    #Piece-wise model with an initial guess of 410, R2 = 0.56, AIC=38.5
AIC(segmented.mod)

This will find the best break-point around an initial guess (called "psi", which I have here of a Climatic Water Deficit of 410 mm). To compare between models (e.g. normal linear or polynomial models versus a piece-wise model), I guess you could use AIC? Anyway, hope this helps!
  