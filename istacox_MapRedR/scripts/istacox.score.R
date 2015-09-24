istacox.score <- function(y.test, y.hat, type=c("spll", "deviance", "logrank")){
  
  if(type=="spll"){
    return(perf=y.hat)
  } else if(type=="deviance"){
    
  } else if(type=="logrank"){
    
  } else {
    print("Sorry, this type of score is not yet implemented !")
  }
}