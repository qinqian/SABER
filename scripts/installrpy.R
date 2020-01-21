verify_installation <- function(package_name){
    if (!require(package_name,character.only=TRUE))
    {
        install.packages(package_name,repos='http://cran.us.r-project.org')
        if(!require(package_name,character.only=TRUE)) stop(paste0("Package ",package_name," not found"))
    }
}

verify_installation("rPython")
verify_installation("stringi")
verify_installation("FField")
verify_installation("Matrix")
