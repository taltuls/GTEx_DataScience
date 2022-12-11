# Regression Colon adjusted by SMRIN and AGE
qn.colon.t=data.frame(t(qn.colon))
qn.colon.t$SMRIN=filtered.pheno.colon$SMRIN
qn.colon.t$AGE=filtered.pheno.colon$AGE


for (s in colnames(qn.colon.t)[1:(length(colnames(qn.colon.t))-1)]){
  lelm=lm(as.numeric(unlist(qn.colon.t[s]))~qn.colon.t$SMRIN+qn.colon.t$AGE)
  qn.colon.t[s]=lelm$residuals
}

drop <- c("SMRIN","AGE")
qn.colon.t = qn.colon.t[,!(names(qn.colon.t) %in% drop)]
qn.colon.res=t(qn.colon.t)

# Regression Brain SMRIN and AGE
qn.brain.t=data.frame(t(qn.brain))
qn.brain.t$SMRIN=filtered.pheno.brain$SMRIN
qn.brain.t$AGE=filtered.pheno.brain$AGE


for (s in colnames(qn.brain.t)[1:(length(colnames(qn.brain.t))-1)]){
  lelm=lm(as.numeric(unlist(qn.brain.t[s]))~qn.brain.t$SMRIN+qn.brain.t$AGE)
  qn.brain.t[s]=lelm$residuals
}

drop <- c("SMRIN","AGE")
qn.brain.t = qn.brain.t[,!(names(qn.brain.t) %in% drop)]
qn.brain.res=t(qn.brain.t)
