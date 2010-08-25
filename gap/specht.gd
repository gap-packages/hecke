#######################################################################
##  SPECHT - specht.g : the kernel of SPECHT                         ##
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
##     (heavily using the GAP3-version by Andrew Mathas)             ##
##                                                                   ##
#######################################################################

BindGlobal("AlgebraObjFamily", NewFamily("AlgebraObjFamily"));
DeclareCategory("IsAlgebraObj", IsComponentObjectRep and IsAttributeStoringRep);
DeclareCategory("IsAlgebraObjModule", IsComponentObjectRep and IsAttributeStoringRep);
DeclareCategory("IsDecompositionMatrix", IsComponentObjectRep and IsAttributeStoringRep);
DeclareCategory("IsCrystalDecompositionMatrix", IsComponentObjectRep and IsAttributeStoringRep);

DeclareAttribute("Characteristic", IsAlgebraObj);
DeclareAttribute("OrderOfQ",IsAlgebraObj);
DeclareAttribute("OrderOfQ",IsAlgebraObjModule);

DeclareProperty("IsZeroCharacteristic", IsAlgebraObj);

DeclareCategory("IsHecke", IsAlgebraObj);
BindGlobal("HeckeType", NewType(AlgebraObjFamily, IsHecke));
##
DeclareCategory("IsSchur", IsAlgebraObj);
BindGlobal("SchurType", NewType(AlgebraObjFamily, IsSchur));

BindGlobal("DecompositionMatrixType", NewType(AlgebraObjFamily, IsDecompositionMatrix));
BindGlobal("CrystalDecompositionMatrixType", NewType(AlgebraObjFamily, IsCrystalDecompositionMatrix));

DeclareCategory("IsHeckeModule", IsAlgebraObjModule);
DeclareCategory("IsHeckeSpecht", IsHeckeModule);
DeclareCategory("IsHeckePIM", IsHeckeModule);
DeclareCategory("IsHeckeSimple", IsHeckeModule);
##
DeclareCategory("IsFockModule", IsAlgebraObjModule);
DeclareCategory("IsFockSpecht", IsFockModule);
DeclareCategory("IsFockPIM", IsFockModule);
DeclareCategory("IsFockSimple", IsFockModule);
##
DeclareCategory("IsSchurModule", IsAlgebraObjModule);
DeclareCategory("IsSchurWeyl", IsSchurModule);
DeclareCategory("IsSchurPIM", IsSchurModule);
DeclareCategory("IsSchurSimple", IsSchurModule);

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

DeclareOperation("Specht", [IsInt]);
DeclareOperation("Specht", [IsInt,IsInt]);
DeclareOperation("Specht", [IsInt,IsInt,IsFunction]);
DeclareOperation("Specht", [IsInt,IsInt,IsFunction,IsString]);
##
DeclareOperation("Schur", [IsInt]);
DeclareOperation("Schur", [IsInt,IsInt]);
DeclareOperation("Schur", [IsInt,IsInt,IsFunction]);
DeclareOperation("Schur", [IsInt,IsInt,IsFunction,IsString]);

DeclareOperation("NewModule",[IsAlgebraObj,IsString,IsList]);
DeclareOperation("NewModule",[IsAlgebraObj,IsString,IsAlgebraObjModule]);
DeclareOperation("NewModule",[IsAlgebraObj,IsString,IsDecompositionMatrix,IsList]);
DeclareOperation("NewModule",[IsAlgebraObj,IsString,IsDecompositionMatrix,IsAlgebraObjModule]);

MakeDispatcherFunc("Hook",[[IsInt]],[2],[2]);
DeclareOperation("DoubleHook",[IsInt,IsInt,IsInt,IsInt]);
DeclareOperation("HeckeOmega",[IsAlgebraObj,IsString,IsInt]);
DeclareOperation("Module",[IsAlgebraObj,IsString,IsInt,IsList]);
DeclareOperation("Module",[IsAlgebraObj,IsString,IsUnivariatePolynomial,IsList]);
DeclareOperation("Module",[IsAlgebraObj,IsString,IsList,IsList]);
DeclareOperation("Collect",[IsAlgebraObj,IsString,IsList,IsList]);

DeclareOperation("MakeSpecht",[IsAlgebraObjModule,IsBool]);
DeclareOperation("MakePIM",[IsAlgebraObjModule,IsBool]);
DeclareOperation("MakeSimple",[IsAlgebraObjModule,IsBool]);

DeclareOperation("InnerProduct",[IsAlgebraObjModule,IsAlgebraObjModule]);
DeclareOperation("Coefficient",[IsAlgebraObjModule,IsList]);
DeclareOperation("PositiveCoefficients",[IsAlgebraObjModule]);
DeclareOperation("IntegralCoefficients",[IsAlgebraObjModule]);

DeclareOperation("\=",[IsAlgebraObjModule,IsAlgebraObjModule]);
DeclareOperation("\+",[IsAlgebraObjModule,IsAlgebraObjModule]);
DeclareOperation("\*",[IsAlgebraObjModule,IsAlgebraObjModule]);
DeclareOperation("\*",[IsScalar,IsAlgebraObjModule]);
DeclareOperation("\*",[IsAlgebraObjModule,IsScalar]);
DeclareOperation("\-",[IsAlgebraObjModule,IsAlgebraObjModule]);
DeclareOperation("\/",[IsAlgebraObjModule,IsScalar]);

DeclareOperation("SetOrdering",[IsAlgebraObj,IsFunction]);

DeclareOperation("SpechtPartitions",[IsHeckeSpecht]);
DeclareOperation("SpechtCoefficients",[IsHeckeSpecht]);

DeclareOperation("ListERegulars",[IsAlgebraObjModule]);
## DeclareProperty("IsERegular",[IsDecompositionMatrix]); ## TODO
MakeDispatcherFunc("SplitECores",
	[[IsAlgebraObjModule],[IsAlgebraObjModule],[IsAlgebraObjModule,IsHeckeSpecht]],
	[ 2									 , 0									,	0],
	[ 2									 , 1									,	2]);
MakeDispatcherFunc("IsSimpleModule", [[IsAlgebraObj]],[2],[2]);
MakeDispatcherFunc("MullineuxMap", 
	[[IsAlgebraObj],[IsInt],[IsAlgebraObjModule],[IsDecompositionMatrix]],
	[ 2						 , 2		 , 0									, 2											],
	[ 2						 , 2		 , 1									, 2											]);
MakeDispatcherFunc("Schaper", [[IsAlgebraObj]],[2],[2]);

DeclareOperation("RInducedModule",[IsAlgebraObj,IsHeckeSpecht,IsInt,IsInt]);
DeclareOperation("SInducedModule",[IsAlgebraObj,IsHeckeSpecht,IsInt,IsInt,IsInt]);
DeclareOperation("RRestrictedModule",[IsAlgebraObj,IsHeckeSpecht,IsInt,IsInt]);
DeclareOperation("SRestrictedModule",[IsAlgebraObj,IsHeckeSpecht,IsInt,IsInt,IsInt]);

DeclareOperation("RInducedModule",[IsAlgebraObj,IsAlgebraObjModule,IsList]);
DeclareOperation("SInducedModule",[IsAlgebraObj,IsAlgebraObjModule,IsList]);
DeclareOperation("RRestrictedModule",[IsAlgebraObj,IsAlgebraObjModule,IsList]);
DeclareOperation("SRestrictedModule",[IsAlgebraObj,IsAlgebraObjModule,IsList]);

DeclareOperation("ReadDecompositionMatrix",[IsAlgebraObj,IsString,IsBool]);
DeclareOperation("ReadDecompositionMatrix",[IsAlgebraObj,IsInt,IsBool]);

