BindGlobal("AlgebraObjFamily", NewFamily("AlgebraObjFamily"));
DeclareCategory("IsAlgebraObj", IsComponentObjectRep and IsAttributeStoringRep);
DeclareCategory("IsAlgebraObjModule", IsComponentObjectRep);
DeclareCategory("IsDecompositionMatrix", IsComponentObjectRep and IsAttributeStoringRep);

DeclareCategory("IsHecke", IsAlgebraObj);
DeclareCategory("IsHeckeDecompositionMatrix", IsDecompositionMatrix);
BindGlobal("HeckeType", NewType(AlgebraObjFamily, IsHecke));
BindGlobal("HeckeDecompositionMatrixType", NewType(AlgebraObjFamily, IsHeckeDecompositionMatrix));
##
DeclareCategory("IsSchur", IsAlgebraObj);
DeclareCategory("IsSchurDecompositionMatrix", IsDecompositionMatrix);
BindGlobal("SchurType", NewType(AlgebraObjFamily, IsSchur));
BindGlobal("SchurDecompositionMatrixType", NewType(AlgebraObjFamily, IsSchurDecompositionMatrix));

DeclareCategory("IsHeckeModule", IsAlgebraObjModule);
DeclareCategory("IsHeckeSpecht", IsHeckeModule); ##TODO chack naming, maybe replace this to IsSpecht
DeclareCategory("IsHeckePIM", IsHeckeModule);
DeclareCategory("IsHeckeSimple", IsHeckeModule); ## Immediate Method -> simple?
##
DeclareCategory("IsSchurModule", IsAlgebraObjModule);
DeclareCategory("IsSchurWeyl", IsSchurModule);
DeclareCategory("IsSchurPIM", IsSchurModule);
DeclareCategory("IsSchurSimple", IsSchurModule); ## Immediate Method -> simple?

BindGlobal("HeckeSpechtType", NewType(AlgebraObjFamily, IsHeckeSpecht));
BindGlobal("HeckePIMType", NewType(AlgebraObjFamily, IsHeckePIM));
BindGlobal("HeckeSimpleType", NewType(AlgebraObjFamily, IsHeckeSimple));
##
BindGlobal("SchurWeylType", NewType(AlgebraObjFamily, IsSchurWeyl));
BindGlobal("SchurPIMType", NewType(AlgebraObjFamily, IsSchurPIM));
BindGlobal("SchurSimpleType", NewType(AlgebraObjFamily, IsSchurSimple));

DeclareOperation("Specht", [IsInt]);
DeclareOperation("Specht", [IsInt,IsPrime]);
#  DeclareOperation("Specht", [IsInt], , );
#  DeclareOperation("Specht", [IsInt], , );
##
DeclareOperation("Schur", [IsInt]);
DeclareOperation("Schur", [IsInt,IsPrime]);
#  DeclareOperation("Schur", [IsInt], , );
#  DeclareOperation("Schur", [IsInt], , );

DeclareOperation("NewModule",[IsAlgebraObj,IsString,IsList]);
DeclareOperation("NewModule",[IsAlgebraObj,IsString,IsAlgebraObjModule]);
DeclareOperation("NewModule",[IsAlgebraObj,IsString,IsDecompositionMatrix,IsList]);
DeclareOperation("NewModule",[IsAlgebraObj,IsString,IsDecompositionMatrix,IsAlgebraObjModule]);

DeclareOperation("OrderOfQ",[IsAlgebraObj]);
DeclareOperation("OrderOfQ",[IsAlgebraObjModule]);
DeclareOperation("SetOrdering",[IsAlgebraObj,IsFunction]);

DeclareOperation("SpechtPartitions",[IsHeckeSpecht]);
DeclareOperation("SpechtCoefficients",[IsHeckeSpecht]);

DeclareOperation("ListERegulars",[IsAlgebraObjModule]);
