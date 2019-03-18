

############################################# put app on server #############################################

library(rsconnect)
# http://docs.rstudio.com/shinyapps.io/getting-started.html#CreateAccount
rsconnect::setAccountInfo(name='heidelberg', token='A8E1B9CEDA6D829860FD5357F3C7241B', secret='r5aX2qshRXaSP5j/g2fjczLhs8CyBbWcAmAVSA0+')
# Tools -> Global Options -> Publishing.

# Rstudio-> Prefeences-> publishing
rsconnect::deployApp("/Users/miyang/Documents/RWTH_Aachen/DREAM_CPTAC/ProteoExplorer/ProteoExplorer" )

rsconnect::appDependencies("/Users/miyang/Documents/RWTH_Aachen/DREAM_CPTAC/ProteoExplorer/ProteoExplorer") # will report the packages, versions, and the repository.

## remove rsconnect folder if it does not work, or in Rstudio delete other account.

## NEW address:  https://heidelberg.shinyapps.io/proteoexplorer




