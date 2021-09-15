## Test environments 

### devtools::check()
* local installation on x86_64-pc-linux-gnu Ubuntu 18.04.5 LTS  R 3.6.1  
* local installation on x86_64-pc-linux-gnu Ubuntu 20.04.3 LTS  R 4.1.1  
0 ERRORS, 0 WARNINGS, 0 NOTES

### devtools::check_rhub()
* Windows Server 2008 R2 SP1, R-release, 32/64 bit      1 NOTE
    checking CRAN incoming feasibility ...          
    Maintainer: 'Laura Fancello <laura.fancello@cea.fr>'
    
* Debian Linux, R-release, GCC                          1 NOTE
    checking CRAN incoming feasibility ...          
    Maintainer: 'Laura Fancello <laura.fancello@cea.fr>'
    
* CentOS 8, stock R from EPEL                           1 ERROR    
    checking package dependencies ... ERROR
    Package suggested but not available: ‘roxygen2’
  
* Ubuntu Linux 20.04.1 LTS, R-release, GCC              1 NOTE
    checking CRAN incoming feasibility ...          
    Maintainer: 'Laura Fancello <laura.fancello@cea.fr>'
    
* Oracle Solaris 10, x86, 32 bit, R-release             1 NOTE
    checking CRAN incoming feasibility ... 
    Maintainer: ‘Laura Fancello <laura.fancello@cea.fr>’
    
* Windows Server 2008 R2 SP1, R-devel, 32/64 bit        1 NOTE, 1 ERROR
  N checking CRAN incoming feasibility
    Maintainer: 'Laura Fancello <laura.fancello@cea.fr>'
  E checking package dependencies
   Package required but not available: 'graph'
   Package suggested but not available: 'BiocStyle'
   
* Fedora Linux, R-devel, clang, gfortran                1 ERROR
  Bioconductor does not yet build and check packages for R version 4.2; 
  
<<<<<<< HEAD
### devtools::check_win_xxx()
devtools::check_win_devel()
devtools::check_win_release()
devtools::check_win_oldrelease()

1 NOTE:
 checking CRAN incoming feasibility ... NOTE
     Maintainer: 'Laura Fancello <laura.fancello@cea.fr>'

### rhub::check_for_cran()
* Windows Server 2008 R2 SP1, R-devel, 32/64 bit  1 ERROR
    checking package dependencies ...
  Package suggested but not available: 'BiocStyle'
  Package required but not available: 'graph'

* Fedora Linux, R-devel, clang, gfortran          1 ERROR
    Bioconductor does not yet build and check packages for R version 4.2;
    packages ‘graph’, ‘BiocStyle’ are not available for this version of R
    
* Ubuntu Linux 20.04.1 LTS, R-release, GCC        1 NOTE
    checking CRAN incoming feasibility ... 
    Maintainer: ‘Laura Fancello <laura.fancello@cea.fr>’

## Gneeral comments on R CMD check results
Always same NOTE:
checking CRAN incoming feasibility ...
     Maintainer: 'Laura Fancello <laura.fancello@cea.fr>'
According to the CRAN maintainer Uwe Ligges: "This is just a note that reminds CRAN maintainers to check that the submission comes actually from his maintainer and not anybody else"

All ERRORS but one concern Bioconductor package availability (graph and biocStyle) on R-devel.
(see:
"Error: Bioconductor does not yet build and check packages for R version 4.2; see https://bioconductor.org/install" 
"Package suggested but not available: 'BiocStyle'.  Package required but not available: 'graph'")
I don't intend to deploy on bioconductor.

I believe that the error from concerning the roxygen2 package is most likely a win-builder issue since 
checking CRAN incoming feasibility ... NOTE
Maintainer: 'Laura Fancello <laura.fancello@cea.fr>'
=======

### Check using win-builder:
* x86_64-w64-mingw32 (64-bit)   R-devel                 1 NOTE
     checking CRAN incoming feasibility ... NOTE
     Maintainer: 'Laura Fancello <laura.fancello@cea.fr>'
* x86_64-w64-mingw32 (64-bit)   R current release       1 ERROR
      checking whether package 'CCs4prot' can be installed ... ERROR
      Installation failed. ('magrittr' not found: maybe not installed for this architecture?)
* x86_64-w64-mingw32 (64-bit)   R old release           1 NOTE
      checking CRAN incoming feasibility ... NOTE
     Maintainer: 'Laura Fancello <laura.fancello@cea.fr>'

  
### rhub::check_for_cran()
* Windows Server 2008 R2 SP1, R-devel, 32/64 bit  1 ERROR
    checking package dependencies ...
  Package suggested but not available: 'BiocStyle'
  Package required but not available: 'graph'

* Fedora Linux, R-devel, clang, gfortran          1 ERROR
    Bioconductor does not yet build and check packages for R version 4.2;
    packages ‘graph’, ‘BiocStyle’ are not available for this version of R
    
* Ubuntu Linux 20.04.1 LTS, R-release, GCC        1 NOTE
    checking CRAN incoming feasibility ... 
    Maintainer: ‘Laura Fancello <laura.fancello@cea.fr>’

## Gneeral comments on R CMD check results
Always same NOTE:
checking CRAN incoming feasibility ...
     Maintainer: 'Laura Fancello <laura.fancello@cea.fr>'
According to the CRAN maintainer Uwe Ligges: "This is just a note that reminds CRAN maintainers to check that the submission comes actually from his maintainer and not anybody else"

All ERRORS but two concern Bioconductor package availability (graph and biocStyle) on R-devel.
(see:
"Error: Bioconductor does not yet build and check packages for R version 4.2; see https://bioconductor.org/install" 
"Package suggested but not available: 'BiocStyle'.  Package required but not available: 'graph'")
I don't intend to deploy on bioconductor.

>>>>>>> bdbb68381ff76d2073e983c6b3f19cb452cf15bb

## Downstream dependencies
There are currently no downstream dependencies for this package. (new package, first submission)
