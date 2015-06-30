istacox.score <- function(y.test, y.hat, type="spll"){
  
  if(type=="spll"){
    return(perf=y.hat)
  } else {
    print("Sorry, this type of score is not yet implemented !")
  }
}