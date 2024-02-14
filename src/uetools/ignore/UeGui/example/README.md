# UeGui example

The UeGui module requires that you have already installed Uedge.

To run this example just start Python and import the run.py script.

<pre>
Python 3.9.13 (main, Aug 25 2022, 18:29:29) 
[Clang 12.0.0 ] :: Anaconda, Inc. on darwin
Type "help", "copyright", "credits" or "license" for more information.
>>> import run

Uedge version 7.0.9.3.0 is newer than available with pip (7.0.9.2.3)

File attributes:
No file attributes, trying to restore
 UEDGE $Name:  $                                                                       
 Wrote file "gridue" with runid:    EFITD    09/07/90      # 66832 ,2384ms                    

 ***** Grid generation has been completed
  Updating Jacobian, npe =                      1
 iter=    0 fnrm=     0.1892981621012411     nfe=      1
  Updating Jacobian, npe =                      2
 iter=    1 fnrm=     0.2825350362855697E-01 nfe=      7
 iter=    2 fnrm=     0.1093103397540877E-02 nfe=     15
 iter=    3 fnrm=     0.1064213093873644E-05 nfe=     25
 iter=    4 fnrm=     0.1317816471728206E-11 nfe=     40


 nksol ---  iterm = 1.
            maxnorm(sf*f(u)) .le. ftol, where maxnorm() is
            the maximum norm function.  u is probably an
            approximate root of f.
 Interpolants created; mype =                   -1
UEDGE>>> 
</pre>


You should receive a popup tabbed notebook on your desktop with buttons you can click to produce various plots.

# Something to try

After doing the initial run above try:

<pre>

UEDGE>>> from uedge.double import uedouble
UEDGE>>> uedouble()
UEDGE>>> run.uedge.bbb.exmain()
 Wrote file "gridue" with runid:    EFITD    09/07/90      # 66832 ,2384ms                    

 ***** Grid generation has been completed
  Updating Jacobian, npe =                      1
 iter=    0 fnrm=      3.494290481079158     nfe=      1
  Updating Jacobian, npe =                      2
 iter=    1 fnrm=      2.327450993943378     nfe=      6
 iter=    2 fnrm=      1.258219584051095     nfe=     13
 iter=    3 fnrm=     0.5666468377340491     nfe=     23
 iter=    4 fnrm=     0.1187548449190309     nfe=     37
 iter=    5 fnrm=     0.1232586397166010E-02 nfe=     52
  Updating Jacobian, npe =                      3
 iter=    6 fnrm=     0.1144783395627690E-05 nfe=     59
 iter=    7 fnrm=     0.1549360109569676E-12 nfe=     69


 nksol ---  iterm = 1.
            maxnorm(sf*f(u)) .le. ftol, where maxnorm() is
            the maximum norm function.  u is probably an
            approximate root of f.
 Interpolants created; mype =                   -1
UEDGE>>> 

</pre>

Then redo some of the plot.
