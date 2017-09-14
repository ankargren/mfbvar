#' Create a new object of the S3 class MFBVAR (the class name can be changed)
#'
#' An object of the S3 class MFBVAR will be created inside the function.
#'
#' The members (attributes) of the object give the settings for the mfbvar model.
#' The functions in the package take the object as the input.
#' The members of the object can be initialized by the function, and can also be changed by the following "set..." functions.
#'
#' 
new_MFBVAR <- function(...)
{
  ret = list()
  class(ret) = "MFBVAR"
  
  return(ret)
}

#'
#'
#'
set_Member <- function(obj, ...)
{
  if(class(obj)!="MFBVAR") stop("message")
  
  #
  
  return(obj)
}
