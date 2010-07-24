##  this creates the documentation, needs: GAPDoc package, latex, pdflatex,
##  mkindex, dvips
##  
##  Call this with GAP.
##

RequirePackage("GAPDoc");

MakeGAPDocDoc("doc", "specht", [], "specht");

GAPDocManualLab("specht");

quit;

