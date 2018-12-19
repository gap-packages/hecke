#######################################################################
##  Hecke - specht.gd : the kernel of Hecke                          ##
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

## Hecke 1.0: August 2010:
##   - initial

BindGlobal("AlgebraObjFamily", NewFamily("AlgebraObjFamily"));
DeclareCategory("IsAlgebraObj", IsComponentObjectRep and IsAttributeStoringRep);
DeclareCategory("IsAlgebraObjModule", IsComponentObjectRep and IsAttributeStoringRep);
DeclareCategory("IsDecompositionMatrix", IsComponentObjectRep and IsAttributeStoringRep);
DeclareCategory("IsCrystalDecompositionMatrix",
  IsDecompositionMatrix and IsComponentObjectRep and IsAttributeStoringRep);

DeclareAttribute("Characteristic", IsAlgebraObj);
DeclareAttribute("OrderOfQ",IsAlgebraObj);
DeclareAttribute("OrderOfQ",IsAlgebraObjModule);
DeclareOperation("CopyDecompositionMatrix",[IsDecompositionMatrix]);

DeclareProperty("IsZeroCharacteristic", IsAlgebraObj);

DeclareCategory("IsHecke", IsAlgebraObj);
BindGlobal("HeckeType", NewType(AlgebraObjFamily, IsHecke));
##
DeclareCategory("IsSchur", IsHecke);
BindGlobal("SchurType", NewType(AlgebraObjFamily, IsSchur));

BindGlobal("DecompositionMatrixType", NewType(AlgebraObjFamily, IsDecompositionMatrix));
BindGlobal("CrystalDecompositionMatrixType", NewType(AlgebraObjFamily, IsCrystalDecompositionMatrix));

DeclareCategory("IsHeckeModule", IsAlgebraObjModule);
DeclareCategory("IsHeckeSpecht", IsHeckeModule);
DeclareCategory("IsHeckePIM", IsHeckeModule);
DeclareCategory("IsHeckeSimple", IsHeckeModule);

DeclareCategory("IsFockModule", IsHeckeModule);
DeclareSynonym("IsFockSpecht", IsHeckeSpecht and IsFockModule);
DeclareSynonym("IsFockPIM", IsHeckePIM and IsFockModule);
DeclareSynonym("IsFockSimple", IsHeckeSimple and IsFockModule);
##
DeclareCategory("IsSchurModule", IsHeckeModule);
DeclareCategory("IsSchurWeyl", IsHeckeSpecht and IsSchurModule);
DeclareSynonym("IsSchurPIM", IsHeckePIM and IsSchurModule);
DeclareSynonym("IsSchurSimple", IsHeckeSimple and IsSchurModule);

DeclareSynonym("IsFockSchurModule", IsFockModule and IsSchurModule);
DeclareSynonym("IsFockSchurWeyl", IsFockSpecht and IsSchurWeyl);
DeclareSynonym("IsFockSchurPIM", IsFockPIM and IsSchurPIM);
DeclareSynonym("IsFockSchurSimple", IsFockSimple and IsSchurSimple);

##
BindGlobal("HeckeSpechtType", NewType(AlgebraObjFamily, IsHeckeSpecht));
BindGlobal("HeckePIMType", NewType(AlgebraObjFamily, IsHeckePIM));
BindGlobal("HeckeSimpleType", NewType(AlgebraObjFamily, IsHeckeSimple));
##
BindGlobal("HeckeSpechtFockType", NewType(AlgebraObjFamily, IsFockSpecht));
BindGlobal("HeckePIMFockType", NewType(AlgebraObjFamily, IsFockPIM));
BindGlobal("HeckeSimpleFockType", NewType(AlgebraObjFamily, IsFockSimple));
##
BindGlobal("SchurWeylType", NewType(AlgebraObjFamily, IsSchurWeyl));
BindGlobal("SchurPIMType", NewType(AlgebraObjFamily, IsSchurPIM));
BindGlobal("SchurSimpleType", NewType(AlgebraObjFamily, IsSchurSimple));
##
BindGlobal("SchurWeylFockType", NewType(AlgebraObjFamily, IsFockSchurWeyl));
BindGlobal("SchurPIMFockType", NewType(AlgebraObjFamily, IsFockSchurPIM));
BindGlobal("SchurSimpleFockType", NewType(AlgebraObjFamily, IsFockSchurSimple));

DeclareOperation("Specht", [IsInt]);
DeclareOperation("Specht", [IsInt,IsInt]);
DeclareOperation("Specht", [IsInt,IsInt,IsFunction]);
DeclareOperation("Specht", [IsInt,IsInt,IsFunction,IsString]);
DeclareOperation("Schur", [IsInt]);
DeclareOperation("Schur", [IsInt,IsInt]);
DeclareOperation("Schur", [IsInt,IsInt,IsFunction]);
DeclareOperation("Schur", [IsInt,IsInt,IsFunction,IsString]);
DeclareOperation("Specht_GenAlgebra", [IsString,IsInt]);
DeclareOperation("Specht_GenAlgebra", [IsString,IsInt,IsInt]);
DeclareOperation("Specht_GenAlgebra", [IsString,IsInt,IsInt,IsFunction]);
DeclareOperation("Specht_GenAlgebra", [IsString,IsInt,IsInt,IsFunction,IsString]);

MakeDispatcherFunc("Hook",[[IsInt]],[2],[2]);
DeclareOperation("DoubleHook",[IsInt,IsInt,IsInt,IsInt]);
DeclareOperation("HeckeOmega",[IsAlgebraObj,IsString,IsInt]);
DeclareOperation("Module",[IsAlgebraObj,IsString,IsInt,IsList]);
DeclareOperation("Module",[IsAlgebraObj,IsString,IsLaurentPolynomial,IsList]);
DeclareOperation("Module",[IsAlgebraObj,IsString,IsList,IsList]);
DeclareOperation("Collect",[IsAlgebraObj,IsString,IsList,IsList]);

MakeDispatcherFunc("MakeSpecht",
[[IsAlgebraObj],
 [IsDecompositionMatrix],
 [IsDecompositionMatrix,IsAlgebraObjModule],
 [IsAlgebraObjModule],
 [IsAlgebraObjModule,IsBool]],
 [2,2,0,0,0],[2,2,2,1,2]);
MakeDispatcherFunc("MakePIM",
[[IsAlgebraObj],
 [IsDecompositionMatrix],
 [IsDecompositionMatrix,IsAlgebraObjModule],
 [IsAlgebraObjModule],
 [IsAlgebraObjModule,IsBool]],
 [2,2,0,0,0],[2,2,2,1,2]);
MakeDispatcherFunc("MakeSimple",
[[IsAlgebraObj],
 [IsDecompositionMatrix],
 [IsDecompositionMatrix,IsAlgebraObjModule],
 [IsAlgebraObjModule],
 [IsAlgebraObjModule,IsBool]],
 [2,2,0,0,0],[2,2,2,1,2]);
MakeDispatcherFunc("MakePIMSpecht",[[IsDecompositionMatrix]],[2],[2]);
MakeDispatcherFunc("MakeFockSpecht",[[IsAlgebraObj]],[2],[2]);
MakeDispatcherFunc("MakeFockPIM",[[IsAlgebraObj]],[2],[2]);

DeclareOperation("InnerProduct",[IsAlgebraObjModule,IsAlgebraObjModule]);
MakeDispatcherFunc("Coefficient",[[IsAlgebraObjModule]],[2],[2]);
DeclareOperation("PositiveCoefficients",[IsAlgebraObjModule]);
DeclareOperation("IntegralCoefficients",[IsAlgebraObjModule]);

#DeclareOperation("\=",[IsAlgebraObjModule,IsAlgebraObjModule]);
DeclareOperation("\+",[IsAlgebraObjModule,IsAlgebraObjModule]);
DeclareOperation("\*",[IsAlgebraObjModule,IsAlgebraObjModule]);
DeclareOperation("\*",[IsScalar,IsAlgebraObjModule]);
DeclareOperation("\*",[IsAlgebraObjModule,IsScalar]);
DeclareOperation("\-",[IsAlgebraObjModule,IsAlgebraObjModule]);
DeclareOperation("\/",[IsAlgebraObjModule,IsScalar]);

DeclareOperation("SetOrdering",[IsAlgebraObj,IsFunction]);

DeclareOperation("SpechtPartitions",[IsHeckeSpecht]);
DeclareOperation("SpechtCoefficients",[IsHeckeSpecht]);

MakeDispatcherFunc("SimpleDimension",
  [[IsDecompositionMatrix],[IsAlgebraObj],[IsAlgebraObj,IsInt],[IsDecompositionMatrix]],
  [2,                       2,             0,                   0],
  [2,                       2,             2,                   1]
);
DeclareOperation("ListERegulars",[IsAlgebraObjModule]);
DeclareOperation("ERegulars",[IsAlgebraObjModule]);
DeclareOperation("ERegulars",[IsDecompositionMatrix]);
MakeDispatcherFunc("SplitECores",
	[[IsAlgebraObjModule],[IsAlgebraObjModule],[IsAlgebraObjModule,IsAlgebraObjModule]],
	[ 2									 , 0									,	0],
	[ 2									 , 1									,	2]);
MakeDispatcherFunc("IsSimpleModule", [[IsAlgebraObj]],[2],[2]);
MakeDispatcherFunc("MullineuxMap",
	[[IsAlgebraObj],[IsInt],[IsDecompositionMatrix],[IsAlgebraObjModule],],
	[ 2						 , 2		 , 2                     , 0								 ],
	[ 2						 , 2		 , 2                     , 1								 ]);
MakeDispatcherFunc("Schaper", [[IsAlgebraObj]],[2],[2]);
DeclareOperation("SchaperMatrix",[IsDecompositionMatrix]);

DeclareOperation("DecompositionNumber",[IsDecompositionMatrix,IsList,IsList]);
DeclareOperation("DecompositionNumber",[IsAlgebraObj,IsList,IsList]);
DeclareOperation("Specht_DecompositionNumber",[IsAlgebraObj,IsList,IsList]);
DeclareOperation("Obstructions",[IsDecompositionMatrix,IsAlgebraObjModule]);
MakeDispatcherFunc("IsNewIndecomposable",
  [[IsAlgebraObj,IsDecompositionMatrix,IsAlgebraObjModule,IsDecompositionMatrix],
   [IsDecompositionMatrix,IsAlgebraObjModule],
   [IsDecompositionMatrix,IsAlgebraObjModule]],
   [5,3,0],[5,3,2]);
MakeDispatcherFunc("RemoveIndecomposable",[[IsDecompositionMatrix]],[2],[2]);
DeclareOperation("MissingIndecomposables",[IsDecompositionMatrix]);
DeclareOperation("CalculateDecompositionMatrix",[IsAlgebraObj,IsInt]);
DeclareOperation("CrystalDecompositionMatrix",[IsAlgebraObj,IsInt]);
DeclareOperation("CrystalDecompositionMatrix",[IsAlgebraObj,IsInt,IsFunction]);
DeclareOperation("InducedDecompositionMatrix",[IsDecompositionMatrix]);
DeclareOperation("InvertDecompositionMatrix",[IsDecompositionMatrix]);
DeclareOperation("SaveDecompositionMatrix",[IsDecompositionMatrix,IsString]);
DeclareOperation("SaveDecompositionMatrix",[IsDecompositionMatrix]);
DeclareOperation("AdjustmentMatrix",[IsDecompositionMatrix,IsDecompositionMatrix]);
DeclareOperation("MatrixDecompositionMatrix",[IsDecompositionMatrix]);
DeclareOperation("DecompositionMatrixMatrix",[IsAlgebraObj,IsMatrix,IsInt]);

MakeDispatcherFunc("RInducedModule",
[[IsAlgebraObjModule],
 [IsAlgebraObj,IsAlgebraObjModule],
 [IsAlgebraObj,IsHeckeSpecht,IsInt,IsInt]],[2,3,0],[2,3,4]);
MakeDispatcherFunc("RRestrictedModule",
[[IsAlgebraObjModule],
 [IsAlgebraObj,IsAlgebraObjModule],
 [IsAlgebraObj,IsHeckeSpecht,IsInt,IsInt]],[2,3,0],[2,3,4]);
MakeDispatcherFunc("SInducedModule",
[[IsAlgebraObjModule],
 [IsAlgebraObj,IsAlgebraObjModule],
 [IsAlgebraObj,IsHeckeSpecht,IsInt,IsInt,IsInt]],[2,3,0],[2,3,5]);
MakeDispatcherFunc("SRestrictedModule",
[[IsAlgebraObjModule],
 [IsAlgebraObj,IsAlgebraObjModule],
 [IsAlgebraObj,IsHeckeSpecht,IsInt,IsInt,IsInt]],[2,3,0],[2,3,5]);
DeclareOperation("qSInducedModule",[IsAlgebraObj,IsAlgebraObjModule,IsInt,IsInt]);
DeclareOperation("qSRestrictedModule",[IsAlgebraObj,IsAlgebraObjModule,IsInt,IsInt]);

#DeclareOperation("\=",[IsDecompositionMatrix,IsDecompositionMatrix]);
DeclareOperation("Store",[IsDecompositionMatrix,IsInt]);
DeclareOperation("Specialized",[IsCrystalDecompositionMatrix,IsInt]);
DeclareOperation("Specialized",[IsCrystalDecompositionMatrix]);
DeclareOperation("Specialized",[IsFockModule,IsInt]);
DeclareOperation("Specialized",[IsFockModule]);
DeclareOperation("Invert",[IsDecompositionMatrix,IsList]);
DeclareOperation("AddIndecomposable",[IsDecompositionMatrix,IsAlgebraObjModule,IsBool]);
DeclareOperation("AddIndecomposable",[IsDecompositionMatrix,IsAlgebraObjModule]);

DeclareOperation("ReadDecompositionMatrix",[IsAlgebraObj,IsString,IsBool]);
DeclareOperation("ReadDecompositionMatrix",[IsAlgebraObj,IsInt,IsBool]);
DeclareOperation("KnownDecompositionMatrix",[IsAlgebraObj,IsInt]);
DeclareOperation("FindDecompositionMatrix",[IsAlgebraObj,IsInt]);

DeclareOperation("FindPq",[IsAlgebraObj,IsList]);
DeclareOperation("FindSq",[IsAlgebraObj,IsList]);
DeclareOperation("FindDq",[IsAlgebraObj,IsList]);

