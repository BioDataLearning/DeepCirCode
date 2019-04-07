#' @title <convStringToMatrix>
#'
#' @param encodedSeqs a vector of string, each string is the one-hot encoded sequence
#' @return x_array_new a 3D array
#' @description This function converts the one-hot encoded sequence into a matrix/3D array as the input for DeepCirCode.
#' @export

convStringToMatrix <- function(encodedSeqs){
  # ensure the character type of encodedSeqs
  encodedSeqs <- as.character(encodedSeqs)
  # create the feature matrix:
  x_array <- array(data = 0, dim = c(4,200, length(encodedSeqs)))
  s <- 1 # sequence/instance index
  r <- 1 # row of the matrix, each row represents A,T/U, G, C
  c <- 1 # column of the matrix, each column represents each nucleotide in the 200nt sequence
  j <- 1 # index of character in the one-hot encoded string
  # store each character into the right place of 3D matrix
  while (s <= length(encodedSeqs)) {
    c <- 1
    while (c <= 200) {
      r <- 1
      while (r <= 4) {
        x_array[r,c,s] <- as.integer(substr(encodedSeqs[s], j,j))
        r <- r + 1
        j <- j + 1
      }
      c <- c + 1
    }
    s <- s + 1
    j <- 1
  }

  #change the index order of x_array to the one keras package required:
  x_array_new <- aperm(x_array,c(3,2,1))
  return(x_array_new)
}
