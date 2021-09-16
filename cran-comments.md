## Test environments 

### devtools::check()
* **local installation on x86_64-pc-linux-gnu Ubuntu 20.04.3 LTS  R 4.1.1**  
**NOTE:** checking CRAN incoming feasibility ...  
      Maintainer: ‘Laura Fancello <laura.fancello@cea.fr>’  


### devtools::check_rhub()
* **Windows Server 2008 R2 SP1, R-devel, 32/64 bit**   
**ERROR:** packages 'graph', 'BiocStyle' are not available for this version of R  
       Bioconductor does not yet build and check packages for R version 4.2  
**NOTE:** checking CRAN incoming feasibility ...   
       Maintainer: 'Laura Fancello <laura.fancello@cea.fr>'  
       possibly mispelled words in DESCRIPTION  

* **Windows Server 2008 R2 SP1, R-release, 32/64 bit**  
**NOTE:** checking CRAN incoming feasibility ...   
       Maintainer: 'Laura Fancello <laura.fancello@cea.fr>'  
       possibly mispelled words in DESCRIPTION  
       
* **Fedora Linux, R-devel, clang, gfortran**  
**ERROR:** packages ‘graph’, ‘BiocStyle’ are not available for this version of R  
       Bioconductor does not yet build and check packages for R version 4.2  

* **Ubuntu Linux 20.04.1 LTS, R-release, GCC**  
**NOTE:** checking CRAN incoming feasibility ...   
       Maintainer: 'Laura Fancello <laura.fancello@cea.fr>'  
       possibly mispelled words in DESCRIPTION  
       
* **Ubuntu Linux 20.04.1 LTS, R-devel, GCC**  
**ERROR:** packages ‘graph’, ‘BiocStyle’ are not available for this version of R  
       Bioconductor does not yet build and check packages for R version 4.2  

* **Debian Linux, R-release, GCC**  
**NOTE:** checking CRAN incoming feasibility ...   
       Maintainer: 'Laura Fancello <laura.fancello@cea.fr>'  
       possibly mispelled words in DESCRIPTION  
       
* **Debian Linux, R-devel, GCC**  
**ERROR:** packages ‘graph’, ‘BiocStyle’ are not available for this version of R  
       Bioconductor does not yet build and check packages for R version 4.2  

  
### devtools::check_win_xxx()
* **devtools::check_win_devel()**  
**NOTE:** checking CRAN incoming feasibility ...   
       Maintainer: 'Laura Fancello <laura.fancello@cea.fr>'  
       possibly mispelled words in DESCRIPTION  

* **devtools::check_win_release()**  
**NOTE:** checking CRAN incoming feasibility ...   
       Maintainer: 'Laura Fancello <laura.fancello@cea.fr>'  
       possibly mispelled words in DESCRIPTION  

* **devtools::check_win_oldrelease()**  
**NOTE:** checking CRAN incoming feasibility ...   
       Maintainer: 'Laura Fancello <laura.fancello@cea.fr>'  
       possibly mispelled words in DESCRIPTION  


## General comments on R CMD check results
**Concerning NOTE: "checking CRAN incoming feasibility ... Maintainer: 'Laura Fancello <laura.fancello@cea.fr>'"**  
According to the CRAN maintainer Uwe Ligges: "This is just a note that reminds CRAN maintainers to check that the submission comes actually from his maintainer and not anybody else".  
  
**Concerning NOTE "Possibly mis-spelled words in DESCRIPTION".**   
The words were not mispelled.  

**Concerning ERRORS "packages ‘graph’, ‘BiocStyle’ are not available for this version of R. Bioconductor does not yet build and check packages for R version 4.2"**   
They only occur on R-devel. I don't intend to deploy on bioconductor.  

## Downstream dependencies
There are currently no downstream dependencies for this package. (new package, first submission)
