#######################################################################
##  Hecke - output.gd : printing functions                           ##
##                                                                   ##
##     A GAP package for calculating the decomposition numbers of    ##
##     Hecke algebras of type A (over fields of characteristic       ##
##     zero). The functions provided are primarily combinatorial in  ##
##     nature. Many of the combinatorial tools from the (modular)    ##
##     representation theory of the symmetric groups appear in the   ##
##     package.                                                      ##
##                                                                   ##
##     These programs, and the enclosed libraries, are distributed   ##
##     under the usual licensing agreements and conditions of GAP.   ##
##                                                                   ##
##     Dmitriy Traytel                                               ##
##     (heavily using the GAP3-package SPECHT 2.4 by Andrew Mathas)  ##
##                                                                   ##
#######################################################################

## Hecke 1.0: June 2010:
##   - initial

DeclareOperation("AlgebraString",[IsAlgebraObj]);
DeclareOperation("ModuleString",[IsAlgebraObjModule,IsBool]);
DeclareOperation("DecompositionMatrixString",[IsDecompositionMatrix]);
DeclareOperation("LabelPartition",[IsList]);
DeclareOperation("StringPartition",[IsList]);
DeclareOperation("TightStringList",[IsList]);
MakeDispatcherFunc("BUG",[[IsString,IsInt],[IsString,IsInt]],[3,0],[3,3]);

