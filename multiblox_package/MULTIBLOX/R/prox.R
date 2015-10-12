prox <-
function(beta,t,alpha) ifelse(abs(beta)<t*alpha,0,beta-t*alpha*sign(beta))
