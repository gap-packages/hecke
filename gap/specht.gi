#######################################################################
##  Hecke - specht.gi : the kernel of Hecke                          ##
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
##   - Translated to GAP4

## SPECHT Change log
## 2.4:
##  - fixed more bugs in H.valuation; returned incorrect answers before
##    when e=0 or e=p (symmetric group case).
##  - fixed bug in Dq(), via Sasha Kleshchev.

## 2.3:
##  - fixed bug in H.valuation  reported by Johannes Lipp
##  - fixed bug in Sq() reported by Johannes Lipp.
##  - corrected FindDecompositionMatrix() so that it updates the matrices
##    CrystalMatrices[] and DecompositionMatrices[] after calculating a
##    crystalized decomposition matrix.

## 2.2: June 1996: various changes requested by the referee.
##  - mainly changing function names.
##  - DecompositionMatrix changed so that it no longer attempts to
##    calculate decomposition matrices in the finite field case; added
##    CalculateDecompositionMatrix() to do this.
##  - replaced Matrix() function with MatrixDecompositionMatrix() and
##    added the function DecompositionMatrixMatix().

## 2.1: April 1996:
##  - Added a filename argument to SaveDecompositionMatrix and made it
##    save the non-decomposition matrices under (more) sensible names; the
##    later is done using the existence of a record component d.matname.
##    Need to do something here about reading such matrices in (this can be
##    done using DecompositionMatrix)..
##  - Changed ReadDecompositionMatrix so that it automatically reads the
##    format of version 1.0 files (and deleted ReadOldDecompositionMatrix).
##    Also fixed a bug here which was causing confusion between Hecke algebra
##    and Schur algebra matrices.
##  - Renamed FindDecompositionMatrix as KnownDecompositionMatrix because it
##    doesn't try to calculate a crytalized decomposition matrix; the new
##    FindDecompositionMatrix will return the decomposition matrix if at all
##    possible. In particular, this fixes a 'bug' in SimpleDimension.
##  - Rewrote AdjustmentMatrix so that it actually works.
##  - Changed P()->S() module conversions so that it no longer requires
##    the projective module to have positive coefficients (and fixed bug in
##    matrix ops).

## 2.0: March 1996:
##   - LLT algorithm implemented for calculating crystal basis of
##     the Fock space and hence by specialization the decomposition
##     matrices for all Hecke algebras over fields of characteristic
##     zero; most of the work is done by the internal function Pq().
##     This required a new set of 'module' types ("Pq", "Sq", and "Dq"),
##     with corresponding operation sets H.operations.X. In particular,
##     these include q-induction and q-restriction, effectively allowing
##     computations with $U_q(\hat sl_e)$ on the Fock space.
##   - crystallized decomposition matrices added and decomposition
##     matrix type overhauled. d.d[c] is now a list of records with
##     partitions replaced by references to d.rows etc. Changed format
##     of decomposition matrix library files so as to handle polynomial
##     entries in the matrices..
##   - printing of Specht records changed allowing more compact and
##     readable notation. eg. S(1,1,1,1)->S(1^4). (see SpechtPrintFn).
##   - reversed order of parts and coeffs records in H.S(), d.d etc;
##     these lists are now sets which improves use of Position().
##   - reorganised the Specht function and the record H() which it returns;
##     in particular added H.info.
##   - extended SInducedModule() and SRestrict to allow multiple inducing and
##     restricting (all residues).

## 1.0: December 1995: initial release.

######################################################################

## Here is a description of the structure of the main records used
## in SPECHT
#### In GAP4 all "operations"-fields do not exist. They are replaced
#### by top-level operations with some filter restrictions

## 1. Specht()
## Specht() is the main function in specht.g, it returns a record 'H'
## which represents the family of Hecke algebras, or Schur algebras,
## corresponding to some fixed field <R> and parameter <q>. 'H' has
## the components:
##   IsSpecht      : is either 'true' or 'false' depending on whether
##                   'H' is a Hecke algebra or Schur algebra resp.
#### This function is obsolete in the GAP4 version - it is replaced
#### by the filter IsHecke
##   S(), P(), D() : these three functions return records which
##                   represent Specht modules, PIMs, and simple
##                   'H' modules respectively. These functions exist
##                   only when 'H' is a Hecke algebra record.
#### Use MakeSpecht(H,...) instead of H.S(...)
##   W(), P(), F() : these are the corresponding functions for Weyl
##                   modules, PIMs, and simple modules of Schur algbras.
#### Use MakeSpecht(S,...) instead of S.W(...)
##   info          : this is a record with components
##                     version: SPECHT version number,
##                     Library: path to SPECHT library files
##                     SpechtDirectory: addition directory searched by
##                            SPECHT (defaults to current directory)
##                     Indeterminate: the indeterminate used by SPECHT
##                            in the LLT algorithm when H.p=0.
##   operations    : apart from the obvious things like the Print()
##                   function for 'H' this record also contains the
##                   operation records for the modules S(), P() etc.
##                   as well as functions for manipulating these records
##                   and accessing decomposition matrix files. The most
##                   most important of these are:
##                     S, P, D, Pq, Sq, Dq : operations records for modules
#### Most functions are now on toplevel
##                     New : creation function for modules. Internally
##                       modules are created via
##                         H.operations.new(module,coeffs,parts)
##                       where module is one of "S", "P", "D", "Sq", "Pq",
##                       or "Dq" ("S" and "D" are used even for Schur
##                       algebras), coeffs is an integer or a *set* of
##                       integers and parts is a partition or a *set* of
##                       partitions. In any programs the use of New is
##                       better than H.S(mu), for example, because the
##                       function names are different for Hecke and Schur
##                       algebras. Note that coeffs and parts must be
##                       ordered reverse lexicographically (ie. they are
##                       *sets*).
#### Use Module(H,...) instead of H.operations.New(...)
##                     Collect: like New() except that  coeffs and parts
##                       need not be sets (and may contain repeats).
#### Use Collect(H,...) instead of H.operations.Collect(...)
##                     NewDecompositionMatrix : creates a decomposition
##                       matrix.
#### Use ReadDecompositionMatrix(H,...) instead of H.operations.ReadDecompositionMatrix(...)
##                     ReadDecompositionMatrix : reads, and returns, a
##                       decomposition matrix file.
#### Use KnownDecompositionMatrix(H,...) instead of H.operations.KnownDecompositionMatrix(...)
##                     KnownDecompositionMatrix : returns a decomposition
##                       matrix of a given size; will either extract this
##                       matrix from Specht's internal lists or call
##                       ReadDecompositionMatrix(), or try to calculate
##                       the decomposition matrix (without using the
##                       crystalized decomposition matrices).
#### Use FindDecompositionMatrix(H,...) instead of H.operations.FindDecompositionMatrix(...)
##                     FindDecompositionMatrix : like KnownDM except that
##                       it will calculate the crystalized DM if needed.
##   Ordering      : a function for ordering partitions; controls how
##                   decomposition matrices for H are printed.
#### Use SetOrdering(...) to control the output
##   e             : order of <q> in <R>
#### Use OrderOfQ(...) to extract the e from an algebra or a module
#### corresponding to an algebra
##   p             : characteristic of <R>
##   valuation     : the valuation map of [JM2]; used primarily by
##                   the q-Schaper theorem.D
##   HeckeRing     : bookkeeping string used primarily in searching for
##                   library files.
##   Pq(), Sq()    : Functions for computing elements of the Fock space
##                   when H.p=0 (used in LLT algorithm). Note that there is
##                   no Dq; also unlike their counter parts S(), P(), and
##                   D() they accept only partitions as arguments.
#### Use MakeFockPIM(H,...) instead of H.Pq(...)
#### Use MakeFockSpecht(H,...) instead of H.Sq(...)
##
## 2. The module functions S(), P() and D() (and Schur equivalents)
## These functions return record 'x' which represents some 'H'--module.
## 'x' is a record with the following components:
##   H      : a pointer back to the corresponding algebra
##   module : one of "S", "P", "D", "Sq", "Pq", or "Dq", (not "W", or "F").
##   coeffs : a *set* of coefficients
##   parts  : the corresponding *set* of partitions
##   operations :
##       + - * / : for algebraic manipulations
##       Print : calls PrintModule
##       Coefficient : returns the coefficient of a given partition
#### Use Coefficient(x,...) instead of x.operations.Coefficient(...)
##       PositiveCoefficients : true if all coefficients are non-negative
##       IntegralCoefficients : true if all coefficients are integral
##       InnerProduct : computes the 'Kronecker' inner product
##       Induce, Restrict, SInduce, SRestrict : induction and restriction
##                 functions taylored to 'x'. These functions convert 'x'
##                 to a linear combination of Specht modules, induce and
##                 then convert back to the type of 'x' (if possible).
##                 Quantized versions are applied as appropriate.
##       S, P, D : functions for rewriting 'x' into the specified type
##                 (so, for example, S('x') rewrites 'x' as a linear
##                 combination of Specht modules).
## 'x'.operations is a pointer to 'H'.operations.('x'.module).
##
## 3. DecompositionMatrices
## Decomposition matrices 'd' in Specht are represented as records with the
## following components:
##
##   d : a list, indexed by d.cols, where each entry is a record
##       corresponding to a column of 'd'; this record has components
##       two * sets* coeffs and parts, where parts is the index of the
##       corresponding partition in d.rows.
##   rows : the *set* of the partitions which make up the rows or 'd'.
##   cols : the *set* of the partitions which make up the rows or 'd'.
##   inverse : a list of records containing the inverse of 'd'. These
##          records are computed only as needed.
##   dimensions : a list of the dimensions of the simple modules; again
##          computed only as needed.
##   IsDecompositionMatrix : false if 'd' is a crystallized decomposition
##          matrix, and true otherwise.
##   H :a pointer back to the corresponding algebra
##   operations :
##     = : equality.
##     Print, TeX, Matrix : printing, TeX, and a GAP Matrix.
##     AddIndecomposable : for adding a PIM to 'd'.
##     Store : for updating Specht's internal record of 'd'.
##     S, P, D: for accessing the entries of 'd' and using 'd' to
##              convert between the various types of 'H'--modules.
##              These are actually records, each containing three
##              functions S(), P(), and D(); so X.Y() tells 'd' how
##              to write an X-module as a linear combination of Y-modules.
##     Invert : calculates D(mu) using 'd'.
##     IsNewIndecomposable : the heart of the 'IsNewIndecomposable'
##              function.
##     Induce : for inducing decomposition matrices (non--crystallized).
##   P : a short-hand for d.H.P('d',<mu>).

###########################################################################

## Specht() is the main function in the package, although in truth it is
## little more than a wrapper for the functions S(), P(), and D().
## Originally, I had these as external functions, but decided that it
## was better to tie these functions to e=H.e as strongly as possible.
InstallMethod(Specht_GenAlgebra,"generate a type-Algebra object",
  [IsString,IsInt,IsInt,IsFunction,IsString],
  function(type,e,p,valuation,HeckeRing)
    local H;
    if not IsPrime(p) and p<>0
    then Error("Specht(<e>,<p>,<val>), <p> must be a prime number");
    fi;

    H := rec(
      e:=e,
      p:=p,
      valuation:=valuation,
      HeckeRing:=HeckeRing,
      ## bits and pieces about H
      info:=rec(version:=PackageInfo("hecke")[1].Version,
                Library:=Directory(
                  Concatenation(DirectoriesPackageLibrary("hecke")[1]![1],"e",
                  String(e),"/")),
      ## We keep a copy of SpechtDirectory in H so that we have a
      ## chance of finding new decomposition matrices when it changes.
                SpechtDirectory:=Directory(".")),

      ## This record will hold any decomposition matrices which Specht()
      ## (or rather its derivatives) read in. This used to be a public
      ## record; it is now private because q-Schur algebra matrices and
      ## Hecke algebra matrices might need to coexist.
      DecompositionMatrices:=[],

      ## for ordering the rows of the decomposition matrices
      ## (as it is common to all decomposition matrices it lives here)
      Ordering:=Lexicographic,
    );

    if p = 0
    then
      ## This list will hold the crystallized decomposition matrices (p=0)
      H.CrystalMatrices:=[];
      H.Indeterminate:=Indeterminate(Integers);
      SetName(H.Indeterminate,"v");
    else H.Indeterminate:=1;
    fi;

    if type = "Schur" then
      Objectify(SchurType,H);
    else
      Objectify(HeckeType,H);
    fi;

    return H;
  end
);

InstallMethod(Specht_GenAlgebra,"generate a type-Algebra object",
  [IsString,IsInt,IsInt,IsFunction],
  function(type,e,p,val) local ring;
    if not IsPrime(p)
    then Error("Specht(<e>,<p>,<val>), <p> must be a prime number");
    fi;
    if e=p then
      ring:=Concatenation("p",String(p),"sym");
      return Specht_GenAlgebra(type,e,p,val,ring);
    else
      return Specht_GenAlgebra(type,e,p,val,"unknown");
    fi;
  end
);

InstallMethod(Specht_GenAlgebra,"generate a type-Algebra object",
  [IsString,IsInt,IsInt],
  function(type,e,p) local val, ring;
    if not IsPrime(p)
    then Error("Specht(<e>,<p>,<val>), <p> must be a prime number");
    fi;
    if e=p then
      ring:=Concatenation("p",String(p),"sym");
      ## return the exponent of the maximum power of p dividing x
      val:=function(x) local i;
        i:=0;
        while x mod p=0 do
          i:=i+1;
          x:=x/p;
        od;
        return i;
      end;
    else
      ring:=Concatenation("e",String(e), "p",String(p));
      ## return the exponent of the maximum power of p that
      ## divides e^x-1.
      val:=function(x) local i;
        if x mod e=0 then return 0;
        else
          i:=0;
          while x mod p=0 do
            i:=i+1;
            x:=x/p;
          od;
          return p^i;
        fi;
      end;
    fi;
    return Specht_GenAlgebra(type,e,p,val,ring);
  end
);

InstallMethod(Specht_GenAlgebra,"generate a type-Algebra object",
  [IsString,IsInt],
  function(type,e) local val;
      if e=0 then val:=x->x;
      else
        val:=function(x)
          if x mod e = 0 then return 1;
          else return 0;
          fi;
        end;
      fi;
      return Specht_GenAlgebra(type,e,0,val,Concatenation("e",String(e), "p0"));
  end
);

InstallMethod(Specht,"generate a Hecke-Algebra object",[IsInt],
  function(e) return Specht_GenAlgebra("Specht",e); end
);
InstallMethod(Specht,"generate a Hecke-Algebra object",[IsInt,IsInt],
  function(e,p) return Specht_GenAlgebra("Specht",e,p); end
);
InstallMethod(Specht,"generate a Hecke-Algebra object",[IsInt,IsInt,IsFunction],
  function(e,p,val) return Specht_GenAlgebra("Specht",e,p,val); end
);
InstallMethod(Specht,"generate a Hecke-Algebra object",[IsInt,IsInt,IsFunction,IsString],
  function(e,p,val,ring) return Specht_GenAlgebra("Specht",e,p,val,ring); end
);
InstallMethod(Schur,"generate a Schur-Algebra object",[IsInt],
  function(e) return Specht_GenAlgebra("Schur",e); end
);
InstallMethod(Schur,"generate a Schur-Algebra object",[IsInt,IsInt],
  function(e,p) return Specht_GenAlgebra("Schur",e,p); end
);
InstallMethod(Schur,"generate a Schur-Algebra object",[IsInt,IsInt,IsFunction],
  function(e,p,val) return Specht_GenAlgebra("Schur",e,p,val); end
);
InstallMethod(Schur,"generate a Schur-Algebra object",[IsInt,IsInt,IsFunction,IsString],
  function(e,p,val,ring) return Specht_GenAlgebra("Schur",e,p,val,ring); end
);

InstallImmediateMethod(Characteristic, IsAlgebraObj, 0,
  function(H) return H!.p; end
);

InstallImmediateMethod(IsZeroCharacteristic, IsAlgebraObj, 0,
  function(H) return H!.p = 0; end
);

InstallImmediateMethod(OrderOfQ, IsAlgebraObj, 0,
  function(H) return H!.e; end
);

InstallImmediateMethod(OrderOfQ, IsAlgebraObjModule, 0,
  function(x) return x!.H!.e; end
);

InstallMethod(SetOrdering,"writing access to H.Ordering",
  [IsAlgebraObj,IsFunction],
  function(H,ord) H!.Ordering := ord; end
);

InstallMethod(SpechtCoefficients,"reading access to S.coeffs",[IsHeckeSpecht],
  function(S) return S!.coeffs; end
);

InstallMethod(SpechtPartitions,"reading access to S.parts",[IsHeckeSpecht],
  function(S) return S!.parts; end
);

## FORMER TOPLEVEL FUNCTIONS ###################################################
#F Calculates the dimensions of the simple modules in d
## Usage:  SimpleDimension(d)   -> prints all simple dimensions
##         SimpleDimension(H,n) -> prints all again
##         SimpleDimension(H,mu) or SimpleDimension(d,mu) -> dim D(mu)
InstallMethod(SimpleDimensionOp,
  "all simple dimensions from decomposition matrix",[IsDecompositionMatrix],
  function(d) local cols, collabel, M, c, x;
    if IsSchur(d!.H) then
      Print("# SimpleDimension() not implemented for Schur algebras\n");
      return fail;
    fi;
    cols:=StructuralCopy(d!.cols);
    if d!.H!.Ordering=Lexicographic then
      cols:=cols{[Length(cols),Length(cols)-1..1]};
    else Sort(cols, d!.H!.Ordering);
    fi;
    cols:=List(cols, c->Position(d!.cols,c));
    collabel:=List([1..Length(cols)], c->LabelPartition(d!.cols[cols[c]]));
    M:=Maximum(List(collabel, Length))+1;

    for c in [1..Length(cols)] do
      Print(String(collabel[c],-M),": ");
      if IsBound(d!.dimensions[cols[c]]) then
        Print(d!.dimensions[cols[c]],"\n");
      else
        x:=MakeSimpleOp(d,d!.cols[cols[c]]);
        if x=fail then Print("not known\n");
        else
          d!.dimensions[cols[c]]:=Sum([1..Length(x!.parts)],
                             r->x!.coeffs[r]*SpechtDimension(x!.parts[r]));
          Print(d!.dimensions[cols[c]],"\n");
        fi;
      fi;
    od;
    return true;
  end
);

InstallMethod(SimpleDimensionOp,
  "simple dimensions of a partition from decomposition matrix",
  [IsDecompositionMatrix,IsList],
  function(d,mu) local c, x;
    c:=Position(d!.cols,mu);
    if c=fail then
      Print("# SimpleDimension(<d>,<mu>), <mu> is not in <d>.cols\n");
      return fail;
    else
      if not IsBound(d!.dimensions[c]) then
        x:=MakeSimpleOp(d,d!.cols[c]);
        if x=fail then return fail;
        else d!.dimensions[c]:=Sum([1..Length(x!.parts)],
                            r->x!.coeffs[r]*SpechtDimension(x!.parts[r]));
        fi;
      fi;
      return d!.dimensions[c];
    fi;
  end
);

InstallMethod(SimpleDimensionOp,
  "all simple dimensions from algebra",[IsAlgebraObj,IsInt],
  function(H,n) local d;
    d:=FindDecompositionMatrix(H,n);
    if d=fail then
      Print("# SimpleDimension(H,n), the decomposition matrix of H_n is ",
            "not known.\n");
      return fail;
    fi;
    return SimpleDimensionOp(d);
  end
);

InstallMethod(SimpleDimensionOp,
  "simple dimensions of a partition from algebra",[IsAlgebraObj,IsList],
  function(H,mu) local d;
    d:=FindDecompositionMatrix(H,Sum(mu));
    if d=fail then
      Print("# SimpleDimension(H,mu), the decomposition matrix of H_Sum(mu) is",
            " not known.\n");
      return fail;
    fi;
    return SimpleDimensionOp(d,mu);
  end
); # SimpleDimension

#P returns a list of the e-regular partitions occurring in x
InstallMethod(ListERegulars,"e-regular partitions of a module",
  [IsAlgebraObjModule],
  function(x) local e,parts,coeffs,p;
    e:=x!.H!.e;
    parts:=x!.parts;
    coeffs:=x!.coeffs;
    if e=0 then return parts;
    elif x=0*x then return [];
    else return List(Filtered([Length(parts),Length(parts)-1..1],
           p->IsERegular(e,parts[p])),p->[coeffs[p], parts[p]]);
    fi;
  end
); # ListERegulars


##P Print the e-regular partitions in x if IsSpecht(x); on the other hand,
### if IsDecompositionMatrix(x) then return the e-regular part of the
### decomposition matrix.
InstallMethod(ERegulars,"e-regular part of the given decomposition matrix",
  [IsDecompositionMatrix],
  function(d) local regs, y, r, len;
    regs:=DecompositionMatrix(d!.H,d!.rows,d!.cols,not IsCrystalDecompositionMatrix(d));
    regs!.d:=[]; #P returns a list of the e-regular partitions occurring in x
    for y in [1..Length(d!.cols)] do
      if IsBound(d!.d[y]) then
        regs!.d[y]:=rec(parts:=[], coeffs:=[]);
        for r in [1..Length(d!.d[y].parts)] do
          len:=Position(d!.cols,d!.rows[d!.d[y].parts[r]]);
          if len<>fail then
            Add(regs!.d[y].parts,len);
            Add(regs!.d[y].coeffs,d!.d[y].coeffs[r]);
          fi;
        od;
      fi;
    od;
    regs!.rows:=regs!.cols;
    return regs;
  end
);

InstallMethod(ERegulars, "print e-regular partitions of a module",
  [IsAlgebraObjModule],
  function(x) local len, regs, y;
    len:=0;
    regs:=ListERegulars(x);
    if regs=[] or IsInt(regs[1]) then Print(regs, "\n");
    else
      for y in regs do
        if (len + 5 + 4*Length(y[2])) > 75 then len:=0; Print("\n"); fi;
        if y[1]<>1 then Print(y[1], "*"); len:=len + 3; fi;
        Print(y[2], "  ");
        len:=len + 5 + 4*Length(y[2]);
      od;
      Print("\n");
    fi;
  end
); # ERegulars

#F Returns true if S(mu)=D(mu) - note that this implies that mu is e-regular
## (if mu is not e-regular, fail is returned).     -- see [JM2]
## IsSimple(H,mu)
##   ** uses H.valuation
InstallMethod(IsSimpleModuleOp,
  "test whether the given partition defines a simple module",
  [IsAlgebraObj,IsList],
  function(H,mu) local mud, simple, r, c, v;
    if not IsERegular(H!.e,mu) then return false;
    elif mu=[] then return true; fi;

    mud:=ConjugatePartition(mu);
    simple:=true; c:=1;
    while simple and c <=mu[1] do
      v:=H!.valuation(mu[1]+mud[c]-c);
      simple:=ForAll([2..mud[c]], r->v=H!.valuation(mu[r]+mud[c]-c-r+1));
      c:=c+1;
    od;
    return simple;
  end
); #IsSimpleModule

#F Split an element up into components which have the same core.
## Usage: SplitECores(x) - returns as list of all block components
##        SplitECores(x,lambda) - returns a list with (i) core lambda,
## (ii) the same core as lambda, or (iii) the same core as the first
## element in lambda if IsSpecht(lambda).
InstallMethod(SplitECoresOp,"for a single module",[IsAlgebraObjModule],
  function(x) local cores, c, cpos, y, cmp;
    if x=fail or x=0*x then return []; fi;

    cores:=[]; cmp:=[];
    for y in [1..Length(x!.parts)] do
      c:=ECore(x!.H!.e, x!.parts[y]);
      cpos:=Position(cores, c);
      if cpos=fail then
        Add(cores, c);
        cpos:=Length(cores);
        cmp[cpos]:=[[],[]];
      fi;
      Add(cmp[cpos][1], x!.coeffs[y]);
      Add(cmp[cpos][2], x!.parts[y]);
    od;
    for y in [1..Length(cmp)] do
      cmp[y]:=Module(x!.H,x!.module,cmp[y][1],cmp[y][2]);
    od;
    return cmp;
  end
);

InstallMethod(SplitECoresOp,"for a module and a partition",
  [IsAlgebraObjModule,IsList],
  function(x,mu) local c, cpos, y, cmp;
    c:=ECore(x!.H!.e, mu);
    cmp:=[ [],[] ];
    for y in [1..Length(x!.parts)] do
      if ECore(x!.H!.e, x!.parts[y])=c then
        Add(cmp[1], x!.coeffs[y]);
        Add(cmp[2], x!.parts[y]);
      fi;
    od;
    cmp:=Module(x!.H,x!.module, cmp[1], cmp[2]);
    return cmp;
  end
);

InstallMethod(SplitECoresOp,"for a module and a specht module",
  [IsAlgebraObjModule,IsAlgebraObjModule],
  function(x,s) local c, cpos, y, cmp;
    c:=ECore(s!.H!.e, s!.parts[Length(x!.parts)]);
    cmp:=[ [],[] ];
    for y in [1..Length(x!.parts)] do
      if ECore(x!.H!.e, x!.parts[y])=c then
        Add(cmp[1], x!.coeffs[y]);
        Add(cmp[2], x!.parts[y]);
      fi;
    od;
    cmp:=Module(x!.H,x!.module, cmp[1], cmp[2]);
    return cmp;
  end
); #SplitECores

#F This function returns the image of <mu> under the Mullineux map using
## the Kleshcehev(-James) algorithm, or the supplied decomposition matrix.
## Alternatively, given a "module" x it works out the image of x under
## Mullineux.
## Usage:  MullineuxMap(e|H|d, mu) or MullineuxMap(x)
InstallMethod(MullineuxMapOp,"image of x under Mullineux",[IsAlgebraObjModule],
  function(x) local e, v;
    e := x!.H!.e;
    if x=fail or not IsERegular(e,x!.parts[Length(x!.parts)]) then
      Print("# The Mullineux map is defined only for e-regular partitions\n");
      return fail;
    fi;
    if x=fail or x=0*x then return fail; fi;
    if IsHeckeSpecht(x) then
      if Length(x!.module)=1 then
        return Collect(x!.H,x!.module,x!.coeffs,
                 List(x!.parts, ConjugatePartition));
      else
        v:=x!.H!.info.Indeterminate;
        return Collect(x!.H,x!.module,
             List([1..Length(x!.coeffs)],
                mu->Value(v^-EWeight(e,x!.parts[mu])*x!.coeffs[mu],v^-1)),
             List(x!.parts,ConjugatePartition) );
      fi;
    elif Length(x!.module)=1 then
      return Sum([1..Length(x!.coeffs)],
               mu->Module(x!.H,x!.module,x!.coeffs[mu],
                     MullineuxMap(e,x!.parts[mu])));
    else
      v:=x!.H!.info.Indeterminate;
      return Sum([1..Length(x!.coeffs)],
               mu->Module(x!.H,x!.module,
                     Value(v^-EWeight(e,x!.parts[mu])*x!.coeffs[mu]),
                     MullineuxMap(e,x!.parts[mu])));
    fi;
  end
);

InstallMethod(MullineuxMapOp,"for ints: image of <mu> under the Mullineux map",
  [IsInt,IsList],
  function(e,mu)
    if not IsERegular(e,mu) then                     ## q-Schur algebra
      Error("# The Mullineux map is defined only for e-regular ",
            "partitions\n");
    fi;
    return PartitionGoodNodeSequence(e,
                  List(GoodNodeSequence(e,mu),x->-x mod e));
  end
);

InstallMethod(MullineuxMapOp,
  "for algebras: image of <mu> under the Mullineux map",
  [IsAlgebraObj,IsList],
  function(H,mu)
    return MullineuxMapOp(H!.e,mu);
  end
);

InstallMethod(MullineuxMapOp,
  "for decomposition matrices: image of <mu> under the Mullineux map",
  [IsDecompositionMatrix,IsList],
  function(d,mu) local e, x;
    e := d!.H!.e;
    if not IsERegular(e,mu) then                     ## q-Schur algebra
      Error("# The Mullineux map is defined only for e-regular ",
            "partitions\n");
    fi;
    x:=d!.H!.P(d,mu);
    if x=fail or x=0*x then
      Print("MullineuxMap(<d>,<mu>), P(<d>,<mu>) not known\n");
      return fail;
    else return ConjugatePartition(x!.parts[1]);
    fi;
    return true;
  end
); # MullineuxMap

#F Calculates the Specht modules in sum_{i>0}S^lambda(i) using the
## q-analogue of Schaper's theorem.
## Uses H.valuation.
##   Usage:  Schaper(H,mu);
InstallMethod(SchaperOp,"calculates Specht modules",[IsAlgebraObj,IsList],
  function(H,mu)
    local mud, schaper, hooklen, c, row, r, s, v;

    Sort(mu); mu:=mu{[Length(mu),Length(mu)-1..1]};
    mud:=ConjugatePartition(mu);
    hooklen:=[];
    for r in [1..Length(mu)] do
      hooklen[r]:=[];
      for c in [1..mu[r]] do
        hooklen[r][c]:=mu[r] + mud[c] - r - c + 1;
      od;
    od;

    schaper:=Module(H,"S",0,[]);
    for c in [1..mu[1]] do
      for row in [1..mud[1]] do
        for r in [row+1..mud[1]] do
          if mu[row] >=c and mu[r] >=c then
            v:=H!.valuation(hooklen[row][c])
                  - H!.valuation(hooklen[r][c]);
            if v<>0 then
              s:=AddRimHook(RemoveRimHook(mu,r,c,mud),row,hooklen[r][c]);
              if s<>fail then
                schaper:=schaper+Module(H,"S",
                                    (-1)^(s[2]+mud[c]-r)*v,s[1]);
              fi;
            fi;
          fi;
        od;
      od;
    od;
    return schaper;
  end
);  #Schaper

#F returns the matrix of upper bounds on the entries in the decomposition
## matrix <d> given by the q-Schaper theorem
## *** undocumented
InstallMethod(SchaperMatrix,"upper bounds of entries of a decomposition matrix",
  [IsDecompositionMatrix],
  function(d) local r, C, c, coeff, sh, shmat;
    shmat:=DecompositionMatrix(d!.H,d!.rows,d!.cols,true);
    shmat!.d:=List(shmat!.cols, c->rec(parts:=[],coeffs:=[]));
    C:=Length(d!.cols)+1; ## this keeps track of which column we're up to
    for r in [Length(d!.rows),Length(d!.rows)-1..1] do
      if d!.rows[r] in d!.cols then C:=C-1; fi;
      sh:=Schaper(d!.H,d!.rows[r]);
      for c in [C..Length(d!.cols)] do
        coeff:=InnerProduct(sh,MakePIMSpechtOp(d,d!.cols[c]));
        if coeff<>fail and coeff<>0*coeff then
          Add(shmat!.d[c].parts,r);
          Add(shmat!.d[c].coeffs,coeff);
        fi;
      od;
    od;
    sh:=[];
    for c in [1..Length(d!.d)] do
      Add(shmat!.d[c].parts, Position(shmat!.rows,shmat!.cols[c]));
      Add(shmat!.d[c].coeffs,1);
    od;
    shmat!.matname:="Schaper matrix";
    return shmat;
  end
);

## FORMER DECOMPOSITION MATRICES TOPLEVEL FUNCTIONS ############################
##############################################################
## Next some functions for accessing decomposition matrices ##
##############################################################

#F Returns the decomposition number d_{mu,nu}; row and column removal
## are used if the projective P(nu) is not already known.
## Usage: DecompositionNumber(H,mu,nu), or
##        DecompositionNumber(d,mu,nu);
## If unable to calculate the decomposition number we return false.
## Note that H.IsSpecht is false if we are looking at decomposition matrices
## of a q-Schur algebra and true for a Hecke algebra.
InstallMethod(DecompositionNumber,
  "for a decomposition matrix and two partitions",
  [IsDecompositionMatrix,IsList,IsList],
  function(d,mu,nu) local Pnu;
    if mu=nu then return 1;
    elif not Dominates(nu,mu) then return 0;
    else
      Pnu:=MakePIMSpechtOp(d,nu);
      if Pnu<>fail then return Coefficient(Pnu,mu); fi;
      return Specht_DecompositionNumber(d!.H,mu,nu);
    fi;
  end
);

InstallMethod(DecompositionNumber,"for an algebra and two partitions",
  [IsAlgebraObj,IsList,IsList],
  function(H,mu,nu) local Pnu;
    Pnu:=MakeSpechtOp(Module(H,"P",1,nu),true);
    if Pnu<>fail then return Coefficient(Pnu,mu); fi;
    if not IsSchur(H) and not IsERegular(H!.e, nu) then
      Error("DecompositionNumber(H,mu,nu), <nu> is not ",H!.e,"-regular");
    fi;
    return Specht_DecompositionNumber(H,mu,nu);
  end
);

InstallMethod(Specht_DecompositionNumber,
  "internal: for an algebra and two partitions",
  [IsAlgebraObj,IsList,IsList],
  function(H,mu,nu) local Pnu, RowAndColumnRemoval;

    ## Next we try row and column removal (James, Theorem 6.18)
    ## (here fn is either the identity or conjugation).
    RowAndColumnRemoval:=function(fn) local m,n,i,d1,d2;
      ## x, mu, and nu as above

      mu:=fn(mu); nu:=fn(nu);

      m:=0; n:=0; i:=1;
      while i<Length(nu) and i<Length(mu) do
        m:=m+mu[i]; n:=n+nu[i];
        if m=n then
          d2:=DecompositionNumber(H, fn(mu{[i+1..Length(mu)]}),
                     fn(nu{[i+1..Length(nu)]}));
          if d2=0 then return d2;
          elif IsInt(d2) then
            d1:=DecompositionNumber(H, fn(mu{[1..i]}),fn(nu{[1..i]}));
            if IsInt(d1) then return d1*d2; fi;
          fi;
        fi;
        i:=i+1;
      od;
      return fail;
    end;

    Pnu:=RowAndColumnRemoval(a->a);
    if Pnu=fail then Pnu:=RowAndColumnRemoval(ConjugatePartition); fi;
    return Pnu;
  end
);

#F Returns a list of those e-regular partitions mu such that Px-P(mu)
## has positive coefficients (ie. those partitions mu such that P(mu)
## could potentially split off Px). Simple minded, but useful.
InstallMethod(Obstructions,"for a decomposition matrix and a module",
  [IsDecompositionMatrix,IsAlgebraObjModule],
  function(d,Px) local obs, mu, Pmu, possibles;
    obs:=[];
    if not IsSchur(d!.H) then
      possibles:=Filtered(Px!.parts, mu->IsERegular(Px!.H!.e, mu));
    else possibles:=Px!.parts;
    fi;
    for mu in possibles do
      if mu<>Px!.parts[Length(Px!.parts)] then
        Pmu:=MakePIMSpechtOp(d,mu);
        if Pmu=fail or PositiveCoefficients(Px-Pmu) then Add(obs,mu); fi;
      fi;
    od;
    return obs{[Length(obs),Length(obs)-1..1]};
  end
);

## Interface to d.operations.IsNewDecompositionMatrix. Returns true
## if <Px> contains an indecomposable not listed in <d> and false
## otherwise. Note that the value of <Px> may well be changed by
## this function. If the argument <mu> is used then we assume
## that all of the decomposition numbers down given by <Px> down to
## <mu> are correct. Note also that if d is the decomposition matrix
## for H(Sym_{r+1}) then the decomposition matrix for H(Sym_r) is passed
## to IsNewDecompositionMatrix.
##   Usage: IsNewIndecomposable(<d>,<Px> [,<mu>]);
## If <mu> is not supplied then we set mu:=true; this
## turns on the message printing in IsNewIndecomposable().
InstallMethod(IsNewIndecomposableOp,
  "for a decomposition matrix and a module",
  [IsDecompositionMatrix,IsAlgebraObjModule],
  function(d,Px) local oldd;
    oldd:=FindDecompositionMatrix(d!.H,Sum(d!.rows[1])-1);
    return IsNewIndecomposableOp(d!.H,d,Px,oldd,[]);
  end
);

InstallMethod(IsNewIndecomposableOp,
  "for a decomposition matrix, a module and a partition",
  [IsDecompositionMatrix,IsAlgebraObjModule,IsList],
  function(d,Px,mu) local oldd;
    oldd:=FindDecompositionMatrix(d!.H,Sum(d!.rows[1])-1);
    return IsNewIndecomposableOp(d!.H,d,Px,oldd,mu);
  end
); # IsNewIndecomposable (toplevel)

##P Removes the columns for <Px> in <d>
InstallMethod(RemoveIndecomposableOp,
  "for a decomposition matrix and a partition",
  [IsDecompositionMatrix,IsList],
  function(d,mu) local r, c;
    c:=Position(d!.cols, mu);
    if c=fail then
      Print("RemoveIndecomposable(<d>,<mu>), <mu> is not listed in <d>\n");
    else Unbind(d!.d[c]);
    fi;
  end
); # RemoveIndecomposable

### Prints a list of the indecomposable missing from d
InstallMethod(MissingIndecomposables,
  "missing entries of a decomposition matrix",
  [IsDecompositionMatrix],
  function(d) local c, missing;
    missing:=List([1..Length(d!.cols)], c->not IsBound(d!.d[c]) );
    if true in missing then
      Print("The following projectives are missing from <d>:\n  ");
      for c in [Length(missing),Length(missing)-1..1] do
        if missing[c] then Print("  ", d!.cols[c]); fi;
      od;
      Print("\n");
    fi;
  end
); # MissingIndecomposables

## When no ordering is supplied then rows are ordered first by length and
## then lexicographically. The rows and columns may also be explicitly
## assigned.
## Usage:
##   DecompositionMatrix(H, n [,ordering]);
##   DecompositionMatrix(H, <file>) ** force Specht() to read <file>
InstallOtherMethod(DecompositionMatrix,"for an algebra and an integer",
  [IsAlgebraObj,IsInt],
  function(H,n) local Px, d, c;
    d:=FindDecompositionMatrix(H,n);

    if d=fail then
      if H!.p>0 and n>2*H!.e then  ## no point even trying
        Print("# This decomposition matrix is not known; use ",
              "CalculateDecompositionMatrix()\n# or ",
              "InducedDecompositionMatrix() to calculate with this matrix.",
              "\n");
        return d;
      fi;
      if not IsSchur(H) then c:=ERegularPartitions(H!.e,n);
      else c:=Partitions(n);
      fi;
      d:=DecompositionMatrix(H,Partitions(n),c,true);
    fi;
    if ForAny([1..Length(d!.cols)],c->not IsBound(d!.d[c])) then
      for c in [1..Length(d!.cols)] do
        if not IsBound(d!.d[c]) then
          Px:=MakeSpechtOp(Module(H,"P",1,d!.cols[c]),true);
          if Px<>fail then AddIndecomposable(d,Px,false);
          else Print("# Projective indecomposable P(",
                     TightStringList(d!.cols[c]),") not known.\n");
          fi;
        fi;
      od;
      Store(d,n);
    fi;
    if d<>fail then   ## can't risk corrupting the internal matrix lists
      d:=CopyDecompositionMatrix(d);
    fi;
    return d;
  end
);

InstallOtherMethod(DecompositionMatrix,
  "for an algebra, an integer and an ordering",
  [IsAlgebraObj,IsInt,IsFunction],
  function(H,n,ord)
    H!.Ordering := ord;
    return DecompositionMatrix(H,n);
  end
);

InstallOtherMethod(DecompositionMatrix,"for an algebra and a filename",
  [IsAlgebraObj,IsString],
  function(H,file) local d;
    d:=ReadDecompositionMatrix(H,file,false);
    if d<>fail and not IsBound(d!.matname) then ## override and copy
      Store(d,Sum(d!.cols[1]));
      MissingIndecomposables(d);
    fi;
    if d<>fail then   ## can't risk corrupting the internal matrix lists
      d:=CopyDecompositionMatrix(d);
    fi;
    return d;
  end
);

InstallOtherMethod(DecompositionMatrix,
  "for an algebra, a filename and an ordering",
  [IsAlgebraObj,IsString,IsFunction],
  function(H,file,ord)
      H!.Ordering := ord;
    return DecompositionMatrix(H,file);
  end
); # DecompositionMatrix

#F Tries to calculate the decomposition matrix d_{H,n} from scratch.
## At present will return only those column indexed by the partitions
## of e-weight less than 2.
InstallMethod(CalculateDecompositionMatrix,"for an algebra and an integer",
  [IsAlgebraObj,IsInt],
  function(H,n) local d, c, Px;
    if not IsSchur(H) then c:=ERegularPartitions(H!.e,n);
    else c:=Partitions(n);
    fi;
    d:=DecompositionMatrix(H,Partitions(n),c,true);
    for c in [1..Length(d!.cols)] do
      if not IsBound(d!.d[c]) then
        Px:=MakeSpechtOp(Module(H,"P",1,d!.cols[c]),true);
        if Px<>fail then AddIndecomposable(d,Px,false);
         else Print("# Projective indecomposable P(",
                    TightStringList(d!.cols[c]),") not known.\n");
        fi;
      fi;
    od;
    return d;
  end
); # CalculateDecompositionMatrix

#F Returns a crystallized decomposition matrix
InstallMethod(CrystalDecompositionMatrix,"for an algebra and an integer",
  [IsAlgebraObj,IsInt],
  function(H,n) local d, Px, c;
    if not IsZeroCharacteristic(H) or IsSchur(H) then
      Error("Crystal decomposition matrices are defined only ",
		         "for Hecke algebras\n         with H!.p=0\n");
    fi;

    d:=ReadDecompositionMatrix(H,n,true);
    if d<>fail then d:=CopyDecompositionMatrix(d);
    else d:=DecompositionMatrix(H,
                Partitions(n),ERegularPartitions(H!.e,n),false);
    fi;
    for c in [1..Length(d!.cols)] do
      if not IsBound(d!.d[c]) then
        AddIndecomposable(d,FindPq(H,d!.cols[c]),false);
      fi;
    od;
    return d;
  end
);

InstallMethod(CrystalDecompositionMatrix,
  "for an algebra, an integer and an ordering",
  [IsAlgebraObj,IsInt,IsFunction],
  function(H,n,ord)
    H!.Ordering := ord;
    return CrystalDecompositionMatrix(H,n);
  end
); # CrystalDecompositionMatrix

## Given a decomposition matrix induce it to find as many columns as
## possible of the next higher matrix using simple minded induction.
## Not as simple minded as it was originally, as it now tries to use
## Schaper's theorem [JM2] to break up troublesome projectives. The
## function looks deceptively simple because all of the work is now
## done by IsNewIndecomposable().
## Usage: InducedDecompositionMatrix(dn)
## in the second form new columns are added to d{n+1}.
InstallMethod(InducedDecompositionMatrix,"induce from decomposition matrix",
  [IsDecompositionMatrix],
  function(d)
    local newd, mu, nu, Px, Py, n,r, needNewline;

    if IsCrystalDecompositionMatrix(d)
    then Error("InducedDecompositionMatrix(d): ",
                 "<d> must be a decomposition matrix.");
    fi;

    needNewline := false;
    n:=Sum(d!.rows[1])+1;
    if n>8 then                            ## print dots to let the user
      PrintTo("*stdout*","# Inducing.");   ## know something is happening.
      needNewline := true;
    fi;

    nu:=Partitions(n);
    if IsSchur(d!.H) then
      newd:=DecompositionMatrix(d!.H, nu, nu, true);
    else newd:=DecompositionMatrix(d!.H, nu,
                ERegularPartitions(d!.H!.e,n),true);
    fi;

    ## add any P(mu)'s with EWeight(mu)<=1 or P(mu)=S(mu) <=> S(mu')=D(mu')
    for mu in newd!.cols do
      if EWeight(d!.H!.e,mu)<=1 then
       AddIndecomposable(newd,
           MakeSpechtOp(Module(d!.H,"P",1,mu),true),false);
      elif IsSimpleModule(d!.H,ConjugatePartition(mu)) then
        AddIndecomposable(newd, Module(d!.H,"S",1,mu),false);
      fi;
    od;

    ## next we r-induce all of the partitions in d so we can just add
    ## them up as we need them later.
    ## (note that this InducedModule() is Specht()'s and not the generic one)
    d!.ind:=List(d!.rows, mu->List([0..d!.H!.e-1],
              r->RInducedModule(d!.H,Module(d!.H,"S",1,mu),d!.H!.e,r)));

    if n<9 then n:=Length(d!.cols)+1; fi; ## fudge for user friendliness

    for mu in [1..Length(d!.cols)] do
      if IsBound(d!.d[mu]) then
        for r in [1..d!.H!.e] do   ## really the e-residues; see ind above
          ## Here we calculate InducedModule(P(mu),H.e,r).
          Px:=Sum([1..Length(d!.d[mu].parts)],
                     nu->d!.d[mu].coeffs[nu]*d!.ind[d!.d[mu].parts[nu]][r]);
          if IsNewIndecomposableOp(d!.H,newd,Px,d,[fail]) then
            if IsERegular(Px!.H!.e,Px!.parts[Length(Px!.parts)]) then
              # can apply MullineuxMap
              nu:=ConjugatePartition(Px!.parts[1]);
              if nu<>MullineuxMap(d!.H!.e,Px!.parts[Length(Px!.parts)]) then
                ## wrong image under the Mullineux map
                BUG("Induce", 7, "nu = ", nu, ", Px = ", Px);
              else   ## place the Mullineux image of Px as well
                AddIndecomposable(newd,MullineuxMap(Px),false);
              fi;
            fi;
            AddIndecomposable(newd,Px,false);
          fi;
        od;
        if mu mod n = 0 then PrintTo("*stdout*","."); needNewline := true; fi;
      fi;
    od;
    Unbind(d!.ind); Unbind(d!.simples); ## maybe we should leave these.

    if needNewline then Print("\n"); fi;
    MissingIndecomposables(newd);
    return newd;
  end
); # InducedDecompositionMatrix

#F Returns the inverse of (the e-regular part of) d. We invert the matrix
## 'by hand' because the matrix routines can't handle polynomial entries.
## This should be much faster than it is???
InstallMethod(InvertDecompositionMatrix,"for a decomposition matrix",
  [IsDecompositionMatrix],
  function(d) local inverse, c, r;
    inverse:=DecompositionMatrix(d!.H,d!.cols,d!.cols,
                                          not IsCrystalDecompositionMatrix(d));

    ## for some reason I can't put this inside the second loop (deleting
    ## the first because d.inverse is not updated this way around...).
    for c in [1..Length(inverse!.cols)] do
      Invert(d,d!.cols[c]);
    od;
    for c in [1..Length(inverse!.cols)] do
      if IsBound(d!.inverse[c]) then
        inverse!.d[c]:=rec(parts:=[], coeffs:=[]);
        for r in [1..c] do
          if IsBound(d!.inverse[r]) and c in d!.inverse[r].parts then
            Add(inverse!.d[c].parts,r);
            Add(inverse!.d[c].coeffs,
                d!.inverse[r].coeffs[Position(d!.inverse[r].parts,c)]);
          fi;
        od;
        if inverse!.d[c]=rec(parts:=[], coeffs:=[]) then Unbind(inverse!.d[c]); fi;
      fi;
    od;
    inverse!.matname:="Inverse matrix";
    return inverse;
  end
); # InvertDecompositionMatrix

#P Saves a full decomposition matrix; actually, only the d, rows, and cols
## records components are saved and the rest calculated when read back in.
## The decomposition matrices are saved in the following format:
##   A_Specht_Decomposition_Matrix:=rec(
##   d:=[[r1,...,rk,d1,...dk],[...],...[]],rows:=[..],cols:=[...]);
## where r1,...,rk are the rows in the first column with corresponding
## decomposition numbers d1,...,dk (if di is a polynomial then it is saved
## as a list [di.valuation,<sequence of di.coffcients]; in particular we
## don't save the polynomial name).
## Usage: SaveDecompositionMatrix(<d>)
##    or  SaveDecompositionMatrix(<d>,<filename>);
InstallMethod(SaveDecompositionMatrix,
  "for a decomposition matrix and a filename",
  [IsDecompositionMatrix,IsString],
  function(d,file)
    local TightList,n,SaveDm,size, r, c, str, tmp;

    n:=Sum(d!.rows[1]);

    size:=SizeScreen();    ## SizeScreen(0 shouldn't affect PrintTo()
    SizeScreen([80,40]);  ## but it does; this is our protection.

    TightList:=function(list) local l, str;
      str:="[";
      for l in list{[1..Length(list)-1]} do
        if IsList(l) then
          l:=TightList(l);
          str:=Concatenation(str,l);
        else str:=Concatenation(str,String(l));
        fi;
        str:=Concatenation(str,",");
      od;
      l:=list[Length(list)];
      if IsList(l) then
        l:=TightList(l);
        str:=Concatenation(str,l);
      else str:=Concatenation(str,String(l));
      fi;
      return Concatenation(str,"]");
    end;

    if d=fail then Error("SaveDecompositionMatrix(<d>), d=fail!!!\n"); fi;

    SaveDm:=function(file)
      AppendTo(file,"## This is a GAP library file generated by \n## Hecke ",
            d!.H!.info.version, "\n\n## This file contains ");
      if IsBound(d!.matname) then
        AppendTo(file,"a(n) ", d!.matname, " for n = ", Sum(d!.rows[1]),"\n");
      else
        if IsCrystalDecompositionMatrix(d) then AppendTo(file,"the crystallized "); fi;
        AppendTo(file,"the decomposition matrix\n## of the ");
        if not IsSchur(d!.H) then
          if d!.H!.e<>d!.H!.p then AppendTo(file,"Hecke algebra of ");
          else AppendTo(file,"symmetric group ");
          fi;
        else AppendTo(file,"q-Schur algebra of ");
        fi;
        AppendTo(file,"Sym(",n,") over a field\n## ");
        if d!.H!.p=0 then AppendTo(file,"of characteristic 0 with ");
        elif d!.H!.p=d!.H!.e then AppendTo(file,"of characteristic ",d!.H!.p,".\n\n");
        else AppendTo(file,"with HeckeRing = ", d!.H!.HeckeRing, ", and ");
        fi;
        if d!.H!.p<>d!.H!.e then AppendTo(file,"e=", d!.H!.e, ".\n\n");fi;
      fi;

      AppendTo(file,"A_Specht_Decomposition_Matrix:=rec(\nd:=[");
      str:="[";
      for c in [1..Length(d!.cols)] do
        if not IsBound(d!.d[c]) then AppendTo(file,str,"]");
        else
          for r in d!.d[c].coeffs do
            if IsLaurentPolynomial(r) then
            tmp:=ShallowCopy(CoefficientsOfLaurentPolynomial(r));
            AppendTo(file,str,"[",tmp[2],",",TightStringList(tmp[1]),"]");
            else AppendTo(file,str,r);
            fi;
            str:=",";
          od;
          for r in d!.d[c].parts do
            AppendTo(file,str,r);
          od;
          AppendTo(file,"]");
          str:=",[";
        fi;
      od;
      AppendTo(file,"],rows:=",TightList(d!.rows));
      AppendTo(file,",cols:=",TightList(d!.cols));
      if IsCrystalDecompositionMatrix(d) then
        AppendTo(file,",crystal:=true");
      fi;
      if IsBound(d!.matname) then AppendTo(file,",matname:=\"",d!.matname,"\""); fi;
      AppendTo(file,");\n");
    end;

    ## the actual saving of d
    SaveDm(file);

    ## now we put d into DecompositionMatrices
    if not IsBound(d!.matname) then Store(d,n); fi;

    SizeScreen(size); # restore screen.
  end
);

InstallMethod(SaveDecompositionMatrix,"for a decomposition matrix",
  [IsDecompositionMatrix],
  function(d) local n,file;
    if d!.H!.HeckeRing="unknown" then
      Print("SaveDecompositionMatrix(d): \n     the base ring of the Hecke ",
            "algebra is unknown.\n     You must set <d>.H.HeckeRing in ",
            "order to save <d>.\n");
      return;
    fi;

    n:=Sum(d!.rows[1]);

    if IsBound(d!.matname) then
      file:=Concatenation(d!.H!.HeckeRing,".",d!.matname{[1]},String(n));
    elif not IsCrystalDecompositionMatrix(d) then
      file:=Concatenation(d!.H!.HeckeRing,".",String(n));
    else  ## crystallized decomposition matrix
      file:=Concatenation("e", String(d!.H!.e), "crys.", String(n));
    fi;

    SaveDecompositionMatrix(d,file);
  end
);# SaveDecompositionMatrix()

#F Returns the 'adjustment matrix' [J] for <d> and <dp>. ie the
## matrix <a> such that <dp>=<a>*<d>.
InstallMethod(AdjustmentMatrix,"for two decomposition matrices",
  [IsDecompositionMatrix,IsDecompositionMatrix],
  function(dp,d) local ad, c, x;
    if d!.cols<>dp!.cols or d!.rows<>dp!.rows then return fail; fi;

    ad:=DecompositionMatrix(dp!.H,dp!.cols,dp!.cols,true);
    ad!.matname:="Adjustment matrix";
    c:=1;
    while ad<>fail and c<=Length(d!.cols) do
      if IsBound(dp!.cols[c]) then
        x:=MakePIMSpechtOp(dp, dp!.cols[c]);
        x!.H:=d!.H;
        x:=MakePIMOp(d,x);
        if x=fail then ad:=fail;
        else AddIndecomposable(ad,x,false);
        fi;
      fi;
      c:=c+1;
    od;
    return ad;
  end
); # AdjustmentMatrix

## Returns the a GAP matrix for the decomposition matrix <d>. Note that
## the rows and columns and <d> are ordered according to H.info.Ordering.
InstallMethod(MatrixDecompositionMatrix,"decomposition matrix -> matrix",
  [IsDecompositionMatrix],
  function(d) local r,c, rows, cols, m;
    rows:=StructuralCopy(d!.rows);
    if d!.H!.Ordering<>Lexicographic then
      Sort(rows,d!.H!.Ordering);
      rows:=List(rows,r->Position(d!.rows,r));
    else rows:=[Length(rows),Length(rows)-1..1];
    fi;
    cols:=StructuralCopy(d!.cols);
    if d!.H!.Ordering<>Lexicographic then
      Sort(cols,d!.H!.Ordering);
      cols:=List(cols,r->Position(d!.cols,r));
    else cols:=[Length(cols),Length(cols)-1..1];
    fi;
    m:=[];
    for r in [1..Length(rows)] do
      m[r]:=[];
      for c in [1..Length(cols)] do
        if IsBound(d!.d[cols[c]]) and rows[r] in d!.d[cols[c]].parts then
          m[r][c]:=d!.d[cols[c]].coeffs[Position(d!.d[cols[c]].parts,rows[r])];
        else m[r][c]:=0;
        fi;
      od;
    od;
    return m;
  end
);

## Given a GAP matrix this function returns a Specht decomposition matrix.
##   H = Specht() record
##   m = matrix: either #reg x #reg, #parts x #reg, or #parts x #parts
##   n = Sym(n)
InstallMethod(DecompositionMatrixMatrix,"matrix -> decomposition matrix",
  [IsAlgebraObj,IsMatrix,IsInt],
  function(H,m,n) local r, c, rows, cols, d;
    rows:=Partitions(n);
    cols:=ERegularPartitions(H,n);
    if Length(rows)<>Length(m) then rows:=cols; fi;
    if Length(cols)<>Length(m[1]) then cols:=rows; fi;
    if Length(rows)<>Length(m) or Length(cols)<>Length(m[1]) then
       Print("# usage: DecompositionMatrixMatrix(H, m, n)\n",
             "   where m is a matrix of an appropriate size.\n");
       return fail;
    fi;
    if ForAll(m, r->ForAll(r, c->IsInt(c) )) then
      d:=DecompositionMatrix(H,StructuralCopy(rows),StructuralCopy(cols),true);
    else  ## presumably crystalized
      d:=DecompositionMatrix(H,StructuralCopy(rows),StructuralCopy(cols),false);
    fi;
    ## now we order the rows and columns properly
    if H!.Ordering<>Lexicographic then Sort(rows, H!.Ordering);
    else rows:=rows{[Length(rows),Length(rows)-1..1]};
    fi;
    rows:=List(d!.rows, r->Position(rows, r) );
    if H!.Ordering<>Lexicographic then Sort(cols, H!.Ordering);
    else cols:=cols{[Length(cols),Length(cols)-1..1]};
    fi;
    cols:=List(d!.cols, c->Position(cols, c) );
    for c in [1..Length(cols)] do
       d!.d[c]:=rec(parts:=[], coeffs:=[]);
       for r in [1..Length(rows)] do
         if m[rows[r]][cols[c]]<>0*m[rows[r]][cols[c]] then ## maybe polynomial
           Add(d!.d[c].parts, r);
           Add(d!.d[c].coeffs, m[rows[r]][cols[c]]);
         fi;
       od;
       if d!.d[c].parts=[] then Unbind(d!.d[c]); fi;
    od;
    return d;
  end
); # DecompositionMatrixMatrix

## OPERATIONS OF FORMER SPECHT RECORD ##########################################

## The following two functions are used by P(), and elsewhere.
##   generate the hook (k,1^n-k)  - as a list - where k=arg
##   actually not quite a hook since if arg is a list (n,k1,k2,...)
##   this returns (k1,k2,...,1^(n-sum k_i))
InstallMethod(HookOp,"for an integer and a list of lists",[IsInt,IsList],
  function(n,K) local k, i;
    k:=Sum(K);
    if k < n then Append(K, List([1..(n-k)], i->1));
    elif k > n then Error("hook, partition ", k, " bigger than ",n, "\n");
    fi;
    return K;
  end
); #Hook

InstallMethod(DoubleHook,"for four integers",[IsInt,IsInt,IsInt,IsInt],
  function(n,x,y,a) local s, i;
    s:=[x];
    if y<>0 then Add(s,y); fi;
    if a<>0 then Append(s, List([1..a],i->2)); fi;
    i:=Sum(s);
    if i < n then
      Append(s, List([1..n-i], i->1));
      return s;
    elif i=n then return s;
    else return [];
    fi;
  end
); #DoubleHook

#### RENAMED Omega -> HeckeOmega
## Returns p(n) - p(n-1,1) + p(n-2,1^2) - ... + (-1)^(n-1)*p(1^n).
## So, S(mu)*Omega(n) is the linear combination of the S(nu)'s where
## nu is obtained by wrapping an n-hook onto mu and attaching the
## sign of the leg length.
InstallMethod(HeckeOmega,"for an algebra, a string and an integer",
  [IsAlgebraObj,IsString,IsInt],
  function(H,module,n)
    return Module(H,module,List([1..n],x->(-1)^(x)),
                     List([1..n],x->Hook(n,x)));
  end
); #HeckeOmega

## MODULES #####################################################################
InstallMethod(Module,"create new module",[IsAlgebraObj,IsString,IsList,IsList],
  function(H,m,c,p)
    local module;
    module := rec(H:=H,module:=m,coeffs:=c,parts:=p);

    if m = "S" and not IsSchur(H) then Objectify(HeckeSpechtType,module);
    elif m = "P" and not IsSchur(H) then Objectify(HeckePIMType,module);
    elif m = "D" and not IsSchur(H) then Objectify(HeckeSimpleType,module);
    elif m = "Sq" and not IsSchur(H) then Objectify(HeckeSpechtFockType,module);
    elif m = "Pq" and not IsSchur(H) then Objectify(HeckePIMFockType,module);
    elif m = "Dq" and not IsSchur(H) then Objectify(HeckeSimpleFockType,module);
    elif m = "W" or m = "S" then Objectify(SchurWeylType,module); module!.module:="W";
    elif m = "P" then Objectify(SchurPIMType,module);
    elif m = "F" or m = "D" then Objectify(SchurSimpleType,module); module!.module:="F";
    elif m = "Wq" or m = "Sq" then Objectify(SchurWeylFockType,module); module!.module:="Wq";
    elif m = "Pq" then Objectify(SchurPIMFockType,module);
    elif m = "Fq" or m = "Dq" then Objectify(SchurSimpleFockType,module); module!.module:="Fq";
    fi;

    return module;
  end
);

InstallMethod(Module,"create new module",[IsAlgebraObj,IsString,IsInt,IsList],
  function(H,m,c,p)
    return Module(H,m,[c],[p]);
  end
);

InstallMethod(Module,"create new module",[IsAlgebraObj,IsString,IsLaurentPolynomial,IsList],
  function(H,m,c,p)
    return Module(H,m,[c],[p]);
  end
);

## Takes two lists, one containing coefficients and the other the
## corresponding partitions, and orders them lexicographically collecting
## like terms on the way. We use a variation on quicksort which is
## induced by the lexicographic order (if parts contains partitions of
## different integers this can lead to an error - which we don't catch).
InstallMethod(Collect,"utility for module generation",
  [IsAlgebraObj,IsString,IsList,IsList],
  function(H, module, coeffs, parts)
    local newx, i, Place, Unplace, places;

    ## inserting parts[i] into places. if parts[i]=[p1,p2,...] then
    ## we insert it into places at places[p1][[p2][...] stopping
    ## at the fist empty position (say places[p1], or places[p1][p2]
    ## etc.). Here we are trying to position parts[i] at
    ## places[p1]...[pd]. Actually, we just insert i rather than
    ## parts[i] and if we find that parts[i]=parts[j] for some j then
    ## we set coeffs[i]:=coeffs[i]+coeffs[j] and don't place j.
    Place:=function(i, places, d) local t;
      if IsInt(places) then
        t:=places;
        if parts[t]=parts[i] then coeffs[t]:=coeffs[t]+coeffs[i];
        else
          places:=[];
          places[parts[t][d]]:=t;
          if parts[i][d]<>parts[t][d] then places[parts[i][d]]:=i;
          else places:=Place(i, places, d);
          fi;
        fi;
      elif places=[] or not IsBound(places[parts[i][d]]) then
        # must be a list
        places[parts[i][d]]:=i;
      else places[parts[i][d]]:=Place(i, places[parts[i][d]], d+1);
      fi;
      return places;
    end;

    Unplace:=function(places) local p, newp, np;
      newp:=[[],[]];
      for p in places do
        if IsInt(p) then if coeffs[p]<>0*coeffs[p] then
          Add(newp[1], coeffs[p]);
          Add(newp[2], StructuralCopy(parts[p])); fi;
        else
          np:=Unplace(p);
          Append(newp[1], np[1]);
          Append(newp[2], np[2]);
        fi;
      od;
      return newp;
    end;

   if parts=[] then return Module(H,module,0,[]);
   elif Length(parts)=1 then return Module(H,module,coeffs,parts);
   else places:=[];
      for i in [1..Length(parts)] do places:=Place(i, places, 1); od;
      newx:=Unplace(places);
      if newx=[[],[]] then return Module(H,module,0,[]);
      else return Module(H,module,newx[1],newx[2]);
      fi;
    fi;
  end  ## H.operations.Collect
);

## MODULE CONVERSION ###########################################################

## Finally the conversion functions S(), P() and D(). All take
## a linear combination of Specht modules and return corresponding
## linear combinations of Specht, indecomposables, and simples resp.
## If they have a problem they return false and print an error
## message unless silent=true.
InstallMethod(MakeSpechtOp,"S()->S()",[IsHeckeSpecht,IsBool],
  function(x,silent) return x; end
);

InstallMethod(MakeSpechtOp,"S[q]()->S[q]()",[IsDecompositionMatrix,IsHeckeSpecht],
  function(d,x) return x; end
);

## Here I only allow for linear combinations of projectives which
## have non-negative coefficients; the reason for this is that I
## can't see how to do it otherwise. The problem is that in the
## Grothendieck ring there are many ways to write a given linear
## combination of Specht modules (or PIMs).
InstallMethod(MakePIMOp,"S()->P()",[IsHeckeSpecht,IsBool],
  function(x,silent) local proj, tmp;
    if x=fail or x=0*x then return x;
    elif x!.parts=[[]] then return Module(x!.H,"P",x!.coeffs[1],[]);
    fi;

    proj:=Module(x!.H,"P",0,[]);
    while x<>fail and x<>0*x and
    ( IsSchur(x!.H) or IsERegular(x!.H!.e,x!.parts[Length(x!.parts)]) ) do
      proj:=proj+Module(x!.H,"P",x!.coeffs[Length(x!.parts)],
                                      x!.parts[Length(x!.parts)]);
      tmp:=MakeSpechtOp(
                Module(x!.H,"P",-x!.coeffs[Length(x!.parts)],
                                      x!.parts[Length(x!.parts)]),true);
      if tmp<>fail then x:=x+tmp; else x:=fail; fi;
    od;
    if x=fail or x<>0*x then
      if not silent then
        Print("# P(<x>), unable to rewrite <x> as a sum of projectives\n");
      fi;
    else return proj;
    fi;
    return fail;
  end
);

InstallMethod(MakePIMOp,"S[q]()->P[q]()",[IsDecompositionMatrix,IsHeckeSpecht],
  function(d,x) local nx, r, c, P, S;
    if x=fail or x=0*x then return x; fi;
    if IsCrystalDecompositionMatrix(d) then P:="Pq"; S:="Sq";
    else P:="P"; S:="S";
    fi;

    nx:=Module(x!.H,P,0,[]);
    while x<>fail and x<>0*x do
      c:=Position(d!.cols,x!.parts[Length(x!.parts)]);
      if c=fail or not IsBound(d!.d[c]) then return fail; fi;
      nx:=nx+Module(x!.H,P,x!.coeffs[Length(x!.parts)],d!.cols[c]);
      x:=x+Module(x!.H,S,-x!.coeffs[Length(x!.parts)]*d!.d[c].coeffs,
                            List(d!.d[c].parts,r->d!.rows[r]));
    od;
    return nx;
  end
);

InstallMethod(MakeSimpleOp,"S()->D()",[IsHeckeSpecht,IsBool],
  function(x,silent) local y, d, simples, r, c;
    if x=fail or x=0*x then return x;
    elif x!.parts=[[]] then return Module(x!.H,"D",x!.coeffs[1],[]);
    fi;

    d:=KnownDecompositionMatrix(x!.H,Sum(x!.parts[1]));
    if d<>fail then
      y:=MakeSimpleOp(d,x);
      if y<>fail then return y; fi;
    fi;

    ## since that didn't work, we use the LLT algorithm when IsBound(H.Pq)
    if IsZeroCharacteristic(x!.H) and not IsSchur(x!.H) then
      return Sum([1..Length(x!.parts)],
                 r->x!.coeffs[r]*Specialized(FindSq(x!.H,x!.parts[r])));
    fi;

    # next, see if we can calculate the answer.
    d:=Concatenation(x!.H!.HeckeRing,"D");
    # finally, we can hope that only partitions of e-weight<2 appear in x
    r:=1; simples:=Module(x!.H,"D",0,[]);
    while simples<>fail and r <= Length(x!.parts) do
      if IsSimpleModule(x!.H, x!.parts[r]) then
        simples:=simples+Module(x!.H,"D",x!.coeffs[r], x!.parts[r]);
      elif IsERegular(x!.H!.e,x!.parts[r]) and EWeight(x!.H!.e,x!.parts[r])=1
      then
        y:=Module(x!.H,"S",1,ECore(x!.H!.e,x!.parts[r]))
                           * HeckeOmega(x!.H,"S",x!.H!.e);
        c:=Position(y!.parts,x!.parts[r]); ## >1 since not IsSimpleModule
        simples:=simples
                  +Module(x!.H,"D",[1,1],[y!.parts[c],y!.parts[c-1]]);
      #### elif IsBound(x.operations.(d)) then ## FIXME not needed anymore?
      ####   simples:=simples+x.operations.(d)(x.parts[r]);
      else simples:=fail;
      fi;
      r:=r+1;
    od;
    if simples=fail and not silent then
      Print("# D(<x>), unable to rewrite <x> as a sum of simples\n");
      return fail;
    else return simples;
    fi;
  end
);

InstallMethod(MakeSimpleOp,"S[q]()->D[q]()",[IsDecompositionMatrix,IsHeckeSpecht],
  function(d,x) local nx, y, r, rr, c, D, core;
    if x=fail or x=0*x then return x; fi;
    if IsCrystalDecompositionMatrix(d) then D:="Dq"; else D:="D"; fi;

    nx:=Module(x!.H,D,0,[]);
    for y in [1..Length(x!.parts)] do
      r:=Position(d!.rows, x!.parts[y]);
      if r=fail then return fail; fi;
      core:=ECore(x!.H!.e,x!.parts[y]);
      c:=Length(d!.cols);
      while c>0 and d!.cols[c]>=x!.parts[y] do
        if IsBound(d!.d[c]) then
          rr:=Position(d!.d[c].parts,r);
          if rr<>fail then nx:=nx+Module(x!.H,D,
                         x!.coeffs[y]*d!.d[c].coeffs[rr],d!.cols[c]);
          fi;
        elif ECore(x!.H!.e,d!.cols[c])=core then return fail;
        fi;
        c:=c-1;
      od;
    od;
    return nx;
  end
);

## The P->S functions are quite involved.

#F Writes x, which is a sum of indecomposables, as a sum of S(nu)'s if
## possible. We first check to see if the decomposition matrix for x is
## stored somewhere, and if not we try to calculate what we need. If we
## can't do this we return false.
InstallMethod(MakeSpechtOp,"P()->S()",[IsHeckePIM,IsBool],
  function(x,silent) local y, c, d, mu, specht;
    if x=fail or x=0*x then return x;
    elif x!.parts=[[]] then return Module(x!.H,"S",x!.coeffs[1],[]);
    fi;
    d:=KnownDecompositionMatrix(x!.H,Sum(x!.parts[1]));
    if d<>fail then
      y:=MakeSpechtOp(d,x);
      if y<>fail then return y; fi;
    fi;

    ## since that didn't work, we use the LLT algorithm when
    ## IsBound(H.Pq)
    if IsZeroCharacteristic(x!.H) then
      if not IsSchur(x!.H) or ForAll(x!.parts, c->IsERegular(x!.H!.e,c)) then
        return Sum([1..Length(x!.parts)],c->
                 x!.coeffs[c]*Specialized(FindPq(x!.H,x!.parts[c])));
      fi;
    fi;

    d:=Concatenation(x!.H!.HeckeRing,"S");
    mu:=1; specht:=Module(x!.H,"S",0,[]);
    while specht<>fail and mu<=Length(x!.parts) do
      if IsSimpleModule(x!.H,ConjugatePartition(x!.parts[mu])) then
        specht:=specht+Module(x!.H,"S",1,x!.parts[mu]);
      elif EWeight(x!.H!.e,x!.parts[mu])=1 then ## wrap e-hooks onto c
        c:=Module(x!.H,"S",1,ECore(x!.H!.e, x!.parts[mu]))
                 * HeckeOmega(x!.H,"S",x!.H!.e);
        y:=Position(c!.parts, x!.parts[mu]);
        specht:=specht+Module(x!.H,"S",[1,1],
                                        [c!.parts[y-1],c!.parts[y]]);
      #### elif IsBound(H.operations.P.(d)) then ## FIXME not needed anymore?
      ####   specht:=specht+x.H.operations.P.(d)(x.parts[mu]);
      else specht:=fail;
      fi;
      mu:=mu+1;
    od;
    if specht<>fail then return specht;
    elif not silent then
      Print("# P(<x>), unable to rewrite <x> as a sum of specht modules\n");
    fi;
    return fail;
  end
);

InstallMethod(MakeSpechtOp,"P[q]()->S[q]()",[IsDecompositionMatrix,IsHeckePIM],
  function(d,x) local S, nx, y, r, c;
    if x=fail or x=0*x then return x; fi;
    if IsCrystalDecompositionMatrix(d) then S:="Sq"; else S:="S"; fi;

    nx:=Module(x!.H,S,0,[]);
    for y in [1..Length(x!.parts)] do
      c:=Position(d!.cols,x!.parts[y]);
      if c=fail or not IsBound(d!.d[c]) then return fail; fi;
      nx:=nx+Module(x!.H,S,x!.coeffs[y]*d!.d[c].coeffs,
                              List(d!.d[c].parts,r->d!.rows[r]));
    od;
    return nx;
  end
);

InstallMethod(MakePIMOp,"P()->P()",[IsHeckePIM,IsBool],
  function(x,silent) return x; end
);

InstallMethod(MakePIMOp,"P[q]()->P[q]()",[IsDecompositionMatrix,IsHeckePIM],
  function(d,x) return x; end
);

InstallMethod(MakeSimpleOp,"P()->D()",[IsHeckePIM,IsBool],
    function(x,silent)
      x:=MakeSpechtOp(x,silent);
      if x=fail then return x;
      else return MakeSimpleOp(x,silent);
      fi;
    end
);

InstallMethod(MakeSimpleOp,"P[q]()->D[q]()",[IsDecompositionMatrix,IsHeckePIM],
  function(d,x)
    x:=MakeSpechtOp(d,x);
    if x=fail then return x;
    else return MakeSimpleOp(d,x);
    fi;
  end
);

#F Writes D(mu) as a sum of S(nu)'s if possible. We first check to see
## if the decomposition matrix for Sum(mu) is stored in the library, and
## then try to calculate it directly. If we are unable to do this either
## we return fail.
InstallMethod(MakeSpechtOp,"D()->S()",[IsHeckeSimple,IsBool],
  function(x,silent) local c, d, y, a;
    if x=fail or x=0*x then return x;
    elif x!.parts=[[]] then return Module(x!.H,"S",x!.coeffs[1],[]);
    fi;

    ## look for the decomposition matrix
    d:=KnownDecompositionMatrix(x!.H,Sum(x!.parts[1]));
    if d<>fail then
      y:=MakeSpechtOp(d,x);
      if y<>fail then return y; fi;
    fi;

    ## since that didn't work, we use the LLT algorithm when IsBound(H.Pq)
    if IsZeroCharacteristic(x!.H) and not IsSchur(x!.H) then
      return Sum([1..Length(x!.parts)],
                   c->x!.coeffs[c]*Specialized(FindDq(x!.H,x!.parts[c])));
    fi;

    #### ## Next, see if we can calculate it. ## FIXME not needed anymore?
    #### d:=Concatenation(H.HeckeRing, "S");
    #### if IsBound(H.operations.D.(d)) then
    ####  return H.operations.D.(d)(x,silent);
    #### fi;

    ## Finally, hope only e-weights<2 are involved.
    c:=1; d:=true; y:=Module(x!.H,"S",0,[]);
    while d and c<=Length(x!.parts) do
      if IsSimpleModule(x!.H, x!.parts[c]) then
        y:=y+Module(x!.H,"S",x!.coeffs[c],x!.parts[c]);
      elif IsERegular(x!.H!.e, x!.parts[c]) and EWeight(x!.H!.e,x!.parts[c])=1
      then ## wrap e-hooks onto c
        a:=Module(x!.H,"S",1,ECore(x!.H!.e,x!.parts[c]))
             * HeckeOmega(x!.H,"S",x!.H!.e);
        a!.parts:=a!.parts{[1..Position(a!.parts, x!.parts[c])]};
        a!.coeffs:=a!.coeffs{[1..Length(a!.parts)]}*(-1)^(1+Length(a!.parts));
        y:=y+a;
      else d:=fail;
      fi;
      c:=c+1;
    od;
    if d<>fail then return y;
    elif not silent then
      Print("# Unable to calculate D(mu)\n");
    fi;
    return fail;
  end
);

InstallMethod(MakeSpechtOp,"D[q]()->S[q]()",[IsDecompositionMatrix,IsHeckeSimple],
  function(d,x) local S, nx, y, c, inv;
    if x=fail or x=0*x then return x; fi;
    if IsCrystalDecompositionMatrix(d) then S:="Sq"; else S:="S"; fi;

    nx:=Module(x!.H,S,0,[]);
    for y in [1..Length(x!.parts)] do
      c:=Position(d!.cols,x!.parts[y]);
      if c=fail then return fail; fi;
      inv:=Invert(d,x!.parts[y]);
      if inv=fail then return inv;
      else nx:=nx+x!.coeffs[y]*inv;
      fi;
    od;
    return nx;
  end
);

InstallMethod(MakePIMOp,"D()->P()",[IsHeckeSimple,IsBool],
  function(x,silent)
      x:=MakeSpechtOp(x,silent);
      if x=fail then return x;
      else return MakePIMOp(x,silent);
      fi;
    end
);

InstallMethod(MakePIMOp,"D[q]()->P[q]()",[IsDecompositionMatrix,IsHeckeSimple],
  function(d,x)
    x:=MakeSimpleOp(d,x);
    if x=fail then return x;
    else return MakePIMOp(d,x);
    fi;
  end
);

InstallMethod(MakeSimpleOp,"D()->D()",[IsHeckeSimple,IsBool],
  function(x,silent) return x; end
);

InstallMethod(MakeSimpleOp,"D[q]()->D[q]()",[IsDecompositionMatrix,IsHeckeSimple],
  function(d,x) return x; end
);

## Finally, change the various conversion functions X()->Y();
## in fact, we only have to change the four non-trivial ones:
##   P() <-> S() <-> D().

InstallMethod(MakePIMOp,"Sq()->Pq()",[IsFockSpecht,IsBool],
  function(x,silent) local proj;
    if x=fail or x=0*x then return x;
    elif x!.parts=[[]] then return Module(x!.H,"Pq",x!.coeffs[1],[]);
    fi;

    proj:=Module(x!.H,"Pq",0,[]);
    while x<>0*x and PositiveCoefficients(x) do
      proj:=proj+Module(x!.H,"Pq",x!.coeffs[Length(x!.parts)],
                                       x!.parts[Length(x!.parts)]);
      x:=x-x!.coeffs[Length(x!.parts)]*FindPq(x!.H,x!.parts[Length(x!.parts)]);
    od;
    if x=0*x then return proj;
    elif not silent then
      Print("# P(<x>), unable to rewrite <x> as a sum of projectives\n");
    fi;
    return fail;
  end
);

InstallMethod(MakeSimpleOp,"Sq()->Dq()",[IsFockSpecht,IsBool],
  function(x,silent) local mu;
    if x=fail or x=0*x then return x;
    elif x!.parts=[[]] then return Module(x!.H,"Dq",x!.coeffs[1],[]);
    fi;
    return Sum([1..Length(x!.parts)],
      mu->x!.coeffs[mu]*FindSq(x!.H,x!.parts[mu]) );
  end
);

InstallMethod(MakeSpechtOp,"Pq()->Sq()",[IsFockPIM,IsBool],
  function(x,silent) local mu;
    if x=fail or x=0*x then return x;
    elif x!.parts=[[]] then return Module(x!.H,"Sq",x!.coeffs[1],[]);
    fi;

    return Sum([1..Length(x!.parts)],
      mu->x!.coeffs[mu]*FindPq(x!.H,x!.parts[mu]) );
  end
);

InstallMethod(MakeSpechtOp,"Dq()->Sq()",[IsFockPIM,IsBool],
  function(x,silent) local mu;
    if x=fail or x=0*x then return x;
    elif x!.parts=[[]] then return Module(x!.H,"Sq",x!.coeffs[1],[]);
    fi;

    return Sum([1..Length(x!.coeffs)],
      mu->x!.coeffs[mu]*FindDq(x!.H,x!.parts[mu]) );
  end
);

## Make<module> now also plays the role of H.<module> functions
InstallMethod(MakeSpechtOp,"H.S(mu)",[IsAlgebraObj,IsList],
  function(H,mu) local z;
    if mu = [] then return Module(H,"S",1,[]);
    else
      if not ForAll(mu,z->IsInt(z)) then return fail; fi;
      z:=StructuralCopy(mu);
      Sort(mu, function(a,b) return a>b;end); # non-increasing
      if mu<>z then
        Print("## S(mu), warning <mu> is not a partition.\n");
      fi;
      if Length(mu)>0 and mu[Length(mu)]<0 then
        Error("## S(mu): <mu> contains negative parts.\n");
      fi;
      z:=Position(mu,0);
      if z<>fail then mu:=mu{[1..z-1]}; fi;  ## remove any zeros from mu
    fi;
    return Module(H,"S", 1, mu);
  end
);

InstallMethod(MakeSpechtOp,"H.S(d,mu)",[IsDecompositionMatrix,IsList],
  function(d,mu) local z;
    if mu = [] then return Module(d!.H,"S",1,[]);
    else
      if not ForAll(mu,z->IsInt(z)) then return fail; fi;
      z:=StructuralCopy(mu);
      Sort(mu, function(a,b) return a>b;end); # non-increasing
      if mu<>z then
        Print("## S(mu), warning <mu> is not a partition.\n");
      fi;
      if Length(mu)>0 and mu[Length(mu)]<0 then
        Error("## S(mu): <mu> contains negative parts.\n");
      fi;
      z:=Position(mu,0);
      if z<>fail then mu:=mu{[1..z-1]}; fi;  ## remove any zeros from mu
    fi;
    return MakeSimpleOp(d,Module(d!.H,"S", 1, mu));
  end
);

InstallMethod(MakePIMOp,"H.P(mu)",[IsAlgebraObj,IsList],
  function(H,mu) local z;
    if mu = [] then return Module(H,"P",1,[]);
    else
      if not ForAll(mu,z->IsInt(z)) then return fail; fi;
      z:=StructuralCopy(mu);
      Sort(mu, function(a,b) return a>b;end); # non-increasing
      if mu<>z then
        Print("## P(mu), warning <mu> is not a partition.\n");
      fi;
      if Length(mu)>0 and mu[Length(mu)]<0 then
        Error("## P(mu): <mu> contains negative parts.\n");
      fi;
      z:=Position(mu,0);
      if z<>fail then mu:=mu{[1..z-1]}; fi;  ## remove any zeros from mu
    fi;
    return Module(H,"P", 1, mu);
  end
);

InstallMethod(MakePIMOp,"H.P(d,mu)",[IsDecompositionMatrix,IsList],
  function(d,mu) local z;
    if mu = [] then return Module(d!.H,"P",1,[]);
    else
      if not ForAll(mu,z->IsInt(z)) then return fail; fi;
      z:=StructuralCopy(mu);
      Sort(mu, function(a,b) return a>b;end); # non-increasing
      if mu<>z then
        Print("## P(mu), warning <mu> is not a partition.\n");
      fi;
      if Length(mu)>0 and mu[Length(mu)]<0 then
        Error("## P(mu): <mu> contains negative parts.\n");
      fi;
      z:=Position(mu,0);
      if z<>fail then mu:=mu{[1..z-1]}; fi;  ## remove any zeros from mu
    fi;
    if not IsSchur(d!.H) and not IsERegular(d!.H!.e,mu) then
      Error("P(mu): <mu>=[",TightStringList(mu),
              "] must be ", d!.H!.e,"-regular\n\n");
    fi;
    return MakeSpechtOp(d,Module(d!.H,"P", 1, mu));
  end
);

InstallMethod(MakeSimpleOp,"H.D(mu)",[IsAlgebraObj,IsList],
  function(H,mu) local z;
    if mu = [] then return Module(H,"D",1,[]);
    else
      if not ForAll(mu,z->IsInt(z)) then return fail; fi;
      z:=StructuralCopy(mu);
      Sort(mu, function(a,b) return a>b;end); # non-increasing
      if mu<>z then
        Print("## D(mu), warning <mu> is not a partition.\n");
      fi;
      if Length(mu)>0 and mu[Length(mu)]<0 then
        Error("## D(mu): <mu> contains negative parts.\n");
      fi;
      z:=Position(mu,0);
      if z<>fail then mu:=mu{[1..z-1]}; fi;  ## remove any zeros from mu
    fi;
    if not IsSchur(H) and not IsERegular(H!.e,mu) then
      Error("D(mu): <mu>=[",TightStringList(mu),
              "] must be ", H!.e,"-regular\n\n");
    fi;
    return Module(H,"D", 1, mu);
  end
);

InstallMethod(MakeSimpleOp,"H.D(d,mu)",[IsDecompositionMatrix,IsList],
  function(d,mu) local z;
    if mu = [] then return Module(d!.H,"D",1,[]);
    else
      if not ForAll(mu,z->IsInt(z)) then return fail; fi;
      z:=StructuralCopy(mu);
      Sort(mu, function(a,b) return a>b;end); # non-increasing
      if mu<>z then
        Print("## D(mu), warning <mu> is not a partition.\n");
      fi;
      if Length(mu)>0 and mu[Length(mu)]<0 then
        Error("## D(mu): <mu> contains negative parts.\n");
      fi;
      z:=Position(mu,0);
      if z<>fail then mu:=mu{[1..z-1]}; fi;  ## remove any zeros from mu
    fi;
    if not IsSchur(d!.H) and not IsERegular(d!.H!.e,mu) then
      Error("D(mu): <mu>=[",TightStringList(mu),
              "] must be ", d!.H!.e,"-regular\n\n");
    fi;
    return MakeSpechtOp(d,Module(d!.H,"D", 1, mu));
  end
);

InstallMethod(MakeSpechtOp,"H.S(x)",[IsAlgebraObjModule],
  function(x) return MakeSpechtOp(x,false); end
);

InstallMethod(MakePIMOp,"H.P(x)",[IsAlgebraObjModule],
  function(x) return MakePIMOp(x,false); end
);

InstallMethod(MakeSimpleOp,"H.D(x)",[IsAlgebraObjModule],
  function(x) return MakeSimpleOp(x,false); end
);

#a lazy helper
InstallMethod(MakePIMSpechtOp,"mu->P(mu)->S(mu)",[IsDecompositionMatrix,IsList],
  function(d,mu)
    return MakeSpechtOp(d,Module(d!.H,"P",1,mu));
  end
);

InstallMethod(MakeFockSpechtOp,"H.Sq(mu)",[IsAlgebraObj,IsList],
  function(H,mu) local z;
    if not ForAll(mu,z->IsInt(z)) then
      Error("usage: H.Sq(<mu1,mu2,...>)\n");
    fi;
    z:=StructuralCopy(mu);
    Sort(mu, function(a,b) return a>b;end); # non-increasing
    if mu<>z then
      Print("## Sq(mu), warning <mu> is not a partition.\n");
    fi;
    if Length(mu)>0 and mu[Length(mu)]<0 then
      Error("## B Sq(mu): <mu> contains negative parts.\n");
    fi;
    z:=Position(mu,0);
    if z<>fail then mu:=mu{[1..z-1]}; fi;  ## remove any zeros from mu
    return Module(H,"Sq",1,mu);
  end
);

InstallMethod(MakeFockPIMOp,"H.Pq(mu)",[IsAlgebraObj,IsList],
  function(H,mu) local z;
    if not ForAll(mu,z->IsInt(z)) then
      Error("usage: H.Pq(<mu1,mu2,...>)\n");
    fi;
    z:=StructuralCopy(mu);
    Sort(mu, function(a,b) return a>b;end); # non-increasing
    if mu<>z then
      Print("## Pq(mu), warning <mu> is not a partition.\n");
    fi;
    if Length(mu)>0 and mu[Length(mu)]<0 then
      Error("## B Pq(mu): <mu> contains negative parts.\n");
    fi;
    z:=Position(mu,0);
    if z<>fail then mu:=mu{[1..z-1]}; fi;  ## remove any zeros from mu
    if not IsERegular(H!.e,mu) then
          Error("Pq(mu): <mu>=[",TightStringList(mu),
                "] must be ", H!.e,"-regular\n\n");
    else return FindPq(H,mu);
    fi;
  end
);

## ARITHMETICS #################################################################
InstallMethod(\=,"compare modules",[IsAlgebraObjModule,IsAlgebraObjModule],
  function(a,b) return a!.H=b!.H and a!.module=b!.module
    and Length(a!.parts)=Length(b!.parts) and Length(a!.coeffs)=Length(b!.coeffs)
    and ForAll(Zip(a!.coeffs,a!.parts),a->a in Zip(b!.coeffs,b!.parts)); end
);

InstallMethod(\+,"add modules",[IsAlgebraObjModule,IsAlgebraObjModule],
  function(a,b)
    local i, j, ab, x;

    if a=fail or b=fail then return fail;
    elif a=0*a then return b;
    elif b=0*b then return a;
    elif a!.H<>b!.H then
      Error("modules belong to different Grothendieck rings");
    fi;

    if a!.module<>b!.module then # only convert to Specht modules if different
      if Length(a!.module) <> Length(b!.module) then
        Error("AddModule(<a>,<b>): can only add modules of same type.");
      fi;
      a:=MakeSpechtOp(a,false);
      b:=MakeSpechtOp(b,false);
      if a=fail or b=fail then return fail; fi;
    fi;

    ## profiling shows _convincingly_ that the method used below to add
    ## a and b is faster than using SortParallel or H.operations.Collect.
    ab:=[[],[]];
    i:=1; j:=1;
    while i <=Length(a!.parts) and j <=Length(b!.parts) do
      if a!.parts[i]=b!.parts[j] then
        x:=a!.coeffs[i]+b!.coeffs[j];
        if x<>0*x then
          Add(ab[1],x);
          Add(ab[2], a!.parts[i]);
        fi;
        i:=i+1; j:=j+1;
      elif a!.parts[i] < b!.parts[j] then
        if a!.coeffs[i]<>0*a!.coeffs[i] then
          Add(ab[1], a!.coeffs[i]);
          Add(ab[2], a!.parts[i]);
        fi;
        i:=i+1;
      else
        if b!.coeffs[j]<>0*b!.coeffs[j] then
          Add(ab[1], b!.coeffs[j]);
          Add(ab[2], b!.parts[j]);
        fi;
        j:=j+1;
      fi;
    od;
    if i <=Length(a!.parts) then
      Append(ab[1], a!.coeffs{[i..Length(a!.coeffs)]});
      Append(ab[2], a!.parts{[i..Length(a!.parts)]});
    elif j <=Length(b!.parts) then
      Append(ab[1], b!.coeffs{[j..Length(b!.coeffs)]});
      Append(ab[2], b!.parts{[j..Length(b!.parts)]});
    fi;
    if ab=[[],[]] then ab:=[ [0],[[]] ]; fi;
    return Module(a!.H, a!.module, ab[1], ab[2]);
  end
); # AddModules

InstallMethod(\-,[IsAlgebraObjModule,IsAlgebraObjModule],
  function(a,b)
    if a=fail or b=fail then return fail;
    else
      b:=Module(b!.H, b!.module, -b!.coeffs, b!.parts);
      return a+b;
    fi;
  end
); # SubModules

InstallMethod(\*,"multiply module by scalar",[IsScalar,IsAlgebraObjModule],
  function(n,b)
    if n = 0
    then return Module(b!.H, b!.module, 0, []);
    else return Module(b!.H, b!.module, n*b!.coeffs, b!.parts);
    fi;
  end
);

InstallMethod(\*,"multiply module by scalar",[IsAlgebraObjModule,IsScalar],
  function(a,n)
    if n = 0
    then return Module(a!.H, a!.module, 0, []);
    else return Module(a!.H, a!.module, n*a!.coeffs, a!.parts);
    fi;
  end
);

InstallMethod(\*,"multiply specht modules",[IsHeckeSpecht,IsHeckeSpecht],
  function(a,b) local x, y, ab, abcoeff, xy, z;
    if a=fail or b=fail then return fail;
    elif a!.H<>b!.H then
      Error("modules belong to different Grothendieck rings");
    fi;
    #a:=MakeSpechtOp(a,false);
    ab:=[[],[]];
    for x in [1..Length(a!.parts)] do
      for y in [1..Length(b!.parts)] do
        abcoeff:=a!.coeffs[x]*b!.coeffs[y];
        if abcoeff<>0*abcoeff then
          z:=LittlewoodRichardsonRule(a!.parts[x], b!.parts[y]);
          Append(ab[1], List(z, xy->abcoeff));
          Append(ab[2], z);
        fi;
      od;
    od;
    if ab=[] then return Module(a!.H, a!.module, 0, []);
    else return Collect(b!.H, b!.module, ab[1], ab[2]);
    fi;
  end
);

InstallMethod(\*,"multiply projective indecomposable modules",[IsHeckePIM,IsHeckePIM],
  function(a,b) local x, nx;
    x:=MakeSpechtOp(a,false) * MakeSpechtOp(b,false);
    nx:=MakePIMOp(x,true);
    if nx<>fail then return nx; else return x; fi;
  end
);

InstallMethod(\*,"multiply simple modules",[IsHeckeSimple,IsHeckeSimple],
  function(a,b) local x, nx;
    x:=MakeSpechtOp(a,false) * MakeSpechtOp(b,false);
    nx:=MakeSimpleOp(x,true);
    if nx<>fail then return nx; else return x; fi;
  end
);

InstallMethod(\*,"multiply modules",[IsAlgebraObjModule,IsAlgebraObjModule],
  function(a,b)
    return MakeSpechtOp(a,false) * MakeSpechtOp(b,false);
  end
); # MulModules

InstallMethod(\/,"divide module by scalar",[IsAlgebraObjModule,IsScalar],
  function(b,n) local x;
    if n=0 then Error("can't divide by 0!\n");
    else return Module(b!.H, b!.module, b!.coeffs/n, b!.parts);
    fi;
  end
); # DivModules

################################################################################

InstallMethod(InnerProduct,"inner product of modules",
  [IsAlgebraObjModule,IsAlgebraObjModule],
  function(a,b) local pr, x, y;
    if a=0*a or b=0*b then return 0;
    elif a!.module<>b!.module then
      a:=MakeSpechtOp(a,true);
      b:=MakeSpechtOp(b,true);
    fi;

    pr:=0; x:=1; y:=1;  # use the fact that a.parts and b.parts are ordered
    while x <=Length(a!.parts) and y <=Length(b!.parts) do
      if a!.parts[x]=b!.parts[y] then
        pr:=pr + a!.coeffs[x]*b!.coeffs[y];
        x:=x + 1; y:=y + 1;
      elif a!.parts[x]<b!.parts[y] then x:=x + 1;
      else y:=y + 1;
      fi;
    od;
    return pr;
  end
);  # InnerProduct

#F Returns the Coefficient of p in x
InstallMethod(CoefficientOp, "extract coefficient of a partition from module",
  [IsAlgebraObjModule,IsList],
  function(x,p) local pos;
    pos:=Position(x!.parts, p);
    if pos=fail then return 0;
    else return x!.coeffs[pos];
    fi;
  end
);  # Coefficient

#F Returns true if all coefficients are non-negative
InstallMethod(PositiveCoefficients, "test if all coefficients are non-negative",
  [IsAlgebraObjModule],
  function(x) local c;
    return ForAll(x!.coeffs, c->c>=0);
  end
);

InstallMethod(PositiveCoefficients,
  "test if all coefficients of a Fock module are non-negative",
  [IsFockModule],
  function(x)
    return ForAll(x!.coeffs, p->( IsInt(p) and p>=0 ) or
                      ForAll(CoefficientsOfLaurentPolynomial(p)[1], c->c>=0) );
  end
); # PositiveCoefficients

#F Returns true if all coefficients are integral
InstallMethod(IntegralCoefficients, "test if all coefficients are integral",
  [IsAlgebraObjModule],
  function(x) local c;
    return ForAll(x!.coeffs, c->IsInt(c));
  end
);

InstallMethod(IntegralCoefficients,
  "test if all coefficients of a Fock module are non-negative",
  [IsFockModule],
  function(x)
    return ForAll(x!.coeffs, p-> ( IsInt(p) and p>=0 ) or
                  ForAll(CoefficientsOfLaurentPolynomial(p)[1], c->IsInt(c)) );
  end
); # IntegralCoefficients

## INDUCTION AND RESTRICTION ###################################################

## The next functions are for restricting and inducing Specht
## modules. They all assume that their arguments are indeed Specht
## modules; conversations are done in H.operations.X.Y() as necessary.

## r-induction: on Specht modules:
InstallMethod(RInducedModuleOp, "r-induction on specht modules",
  [IsAlgebraObj,IsHeckeSpecht,IsInt,IsInt],
  function(H, a, e, r) local ind, x, i, j, np;
    ind:=[[],[]];
    for x in [1..Length(a!.parts)] do
      for i in [1..Length(a!.parts[x])] do
        if (a!.parts[x][i] + 1 - i) mod e=r then
          if i=1 or a!.parts[x][i-1] > a!.parts[x][i] then
            np:=StructuralCopy(a!.parts[x]);
            np[i]:=np[i]+1;
            Add(ind[1], a!.coeffs[x]);
            Add(ind[2], np);
          fi;
        fi;
      od;
      if ( -Length(a!.parts[x]) mod e)=r then
        np:=StructuralCopy(a!.parts[x]); Add(np, 1);
        Add(ind[1],a!.coeffs[x]);
        Add(ind[2],np);
      fi;
    od;
    if ind=[ [],[] ] then return Module(H,"S",0,[]);
    else return Collect(H,"S", ind[1], ind[2]);
    fi;
  end
); # RInducedModule

## String-induction: add s r's from each partition in x (ignoring
## multiplicities). Does both standard and q-induction.

## We look at the size of x.module to decide whether we want to use
## ordinary induction or q-induction (in the Fock space). We could
## write H.operations.X.SInduce to as to make this choice for us, or
## do q-induction always, setting v=1 afterwards, but this seems the
## better choice.
InstallMethod(SInducedModuleOp,"string-induction on specht modules",
  [IsAlgebraObj,IsHeckeSpecht,IsInt,IsInt,IsInt],
  function(H, x, e, s, r) local coeffs, parts, y, z, sinduced;
    # add n nodes of residue r to the partition y from the i-th row down
    sinduced:=function(y, n, e, r, i) local ny, j, z;
      ny:=[];
      for j in [i..Length(y)-n+1] do
        if r=(y[j] - j + 1) mod e then
          if j=1 or y[j] < y[j-1] then
            z:=StructuralCopy(y);
            z[j]:=z[j] + 1; # only one node of residue r can be added
            if n=1 then Add(ny, z);   # no more nodes to add
            else Append(ny, sinduced(z, n-1, e, r, j+1));
            fi;
          fi;
        fi;
      od;
      return ny;
    end;

    if s=0 then return Module(x!.H,x!.module,1,[]); fi;
    coeffs:=[]; parts:=[];
    for y in [1..Length(x!.parts)] do
      Append(parts,sinduced(x!.parts[y], s, e, r, 1));
      Append(coeffs,List([1..Length(parts)-Length(coeffs)],r->x!.coeffs[y]));
      if r=( -Length(x!.parts[y]) mod e) then # add a node to the bottom
        z:=StructuralCopy(x!.parts[y]);
        Add(z,1);
        if s > 1 then                        # need to add some more nodes
          Append(parts,sinduced(z, s-1, e, r, 1));
          Append(coeffs,List([1..Length(parts)-Length(coeffs)],
                             r->x!.coeffs[y]));
        else Add(coeffs, x!.coeffs[y]); Add(parts, z);
        fi;
      fi;
    od;

    if coeffs=[] then return Module(H, x!.module,0,[]);
    else return Collect(H, x!.module, coeffs, parts);
    fi;
  end
);  # SInducedModule

## r-restriction
InstallMethod(RRestrictedModuleOp,"r-restriction on specht modules",
  [IsAlgebraObj,IsHeckeSpecht,IsInt,IsInt],
  function(H, a, e, r) local ind, x, i, j, np;
    ind:=[[],[]];
    for x in [1..Length(a!.parts)] do
      for i in [1..Length(a!.parts[x])] do
        if (a!.parts[x][i] - i) mod e=r then
          np:=StructuralCopy(a!.parts[x]);
          if i=Length(a!.parts[x]) or np[i] > np[i+1] then
            np[i]:=np[i] - 1;
            if np[i]=0 then Unbind(np[i]); fi;
            Add(ind[1], a!.coeffs[x]); Add(ind[2], np);
          fi;
        fi;
      od;
    od;
    if ind=[ [],[] ] then return Module(H,"S",0,[]);
    else return Collect(H,"S", ind[1], ind[2]);
    fi;
  end
); #RRestrictedModule

## string-restriction: remove m r's from each partition in x
InstallMethod(SRestrictedModuleOp,"string-restriction on specht modules",
  [IsAlgebraObj,IsHeckeSpecht,IsInt,IsInt,IsInt],
  function(H,x,e,s,r) local coeffs, parts, y, i, srestricted;
    ## remove n nodes from y from the ith row down
    srestricted:=function(y, n, e, r, i) local ny, j, z;
      ny:=[];
      for j in [i..Length(y)-n+1] do
        if r=(y[j] - j) mod e then
          if j=Length(y) or y[j] > y[j+1] then
            z:=StructuralCopy(y);
            z[j]:=z[j] - 1;
            if z[j]=0 then   # n must be 1
              Unbind(z[j]);
              Add(ny, z);
            elif n=1 then Add(ny, z); # no mode nodes to remove
            else Append(ny, srestricted(z, n-1, e, r, j+1));
            fi;
          fi;
        fi;
      od;
      return ny;
    end;

    coeffs:=[]; parts:=[];
    for y in [1..Length(x!.parts)] do
      Append(parts, srestricted(x!.parts[y], s, e, r, 1));
      Append(coeffs,List([1..Length(parts)-Length(coeffs)],i->x!.coeffs[y]));
    od;
    if parts=[] then return Module(H,"S",0,[]);
    else return Collect(H,"S", coeffs, parts);
    fi;
  end
);  # SRestrictedModule

## Induction and restriction; for S()
InstallMethod(RInducedModuleOp, "r-induction for specht modules",
  [IsAlgebraObj, IsHeckeSpecht, IsList],
  function(H, x, list) local r;
    if x=fail or x=0*x then return x;
    elif list=[] then return RInducedModule(H,x,1,0);
    elif H!.e=0 then
      Error("Induce, r-induction is not defined when e=0.");
    elif ForAny(list,r-> r>=H!.e or r<0) then
      Error("Induce, r-induction is defined only when 0<=r<e.\n");
    else
      for r in list do
        x:=RInducedModule(H,x,H!.e,r);
      od;
      return x;
    fi;
  end
);

InstallMethod(RInducedModuleOp, "r-induction for projective indecomposable modules",
  [IsAlgebraObj, IsHeckePIM, IsList],
  function(H, x, list) local nx;
    x:=RInducedModule(H,MakeSpechtOp(x,false),list);
    if x=fail or x=0*x then return x; fi;
    nx:=MakePIMOp(x,false);
    if nx<>fail then return nx; else return x; fi;
  end
);

InstallMethod(RInducedModuleOp, "r-induction for simple modules",
  [IsAlgebraObj, IsHeckeSimple, IsList],
  function(H, x, list) local nx;
    x:=RInducedModule(H,MakeSpechtOp(x,false),list);
    if x=fail or x=0*x then return x; fi;
    nx:=MakeSimpleOp(x,false);
    if nx<>fail then return nx; else return x; fi;
  end
);

InstallMethod(RInducedModuleOp, "r-induction for Fock space specht modules",
  [IsAlgebraObj,IsFockSpecht,IsList],
  function(H, x, list) local r;
    if list=[] then return Sum([0..H!.e-1],r->qSInducedModule(H,x,1,r));
    elif H!.e=0 then
      Error("Induce, r-induction is not defined when e=0.");
    elif ForAny(list,r-> r>=H!.e or r<0) then
      Error("Induce, r-induction is defined only when 0<=r<e.\n");
    else
      for r in list do   ## we could do slightly better here
        x:= qSInducedModule(H,x,1,r);
      od;
      return x;
    fi;
  end
);

InstallMethod(RInducedModuleOp,
  "r-induction for Fock space projective indecomposable modules",
  [IsAlgebraObj,IsFockPIM,IsList],
  function(H,x,list)
    return MakePIMOp(RInducedModule(H,MakeSpechtOp(x,false),list),false);
  end
);

InstallMethod(RInducedModuleOp,
  "r-induction for Fock space simple modules",
  [IsAlgebraObj,IsFockSimple,IsList],
  function(H,x,list)
    return MakeSimpleOp(RInducedModule(H,MakeSpechtOp(x,false),list),false);
  end
); # RInducedModule

InstallMethod(RRestrictedModuleOp, "r-restriction for specht modules",
  [IsAlgebraObj, IsHeckeSpecht, IsList],
  function(H, x, list) local r;
    if x=fail or x=0*x then return x;
    elif list=[] then return RRestrictedModule(H,x,1,0);
    elif H!.e=0 then
      Error("Restrict, r-restriction is not defined when e=0.");
   elif ForAny(list,r-> r>=H!.e or r<0) then
      Error("Restrict, r-restriction is defined only when 0<=r<e.\n");
    else
      for r in list do
        x:=RRestrictedModule(H,x,H!.e,r);
      od;
      return x;
    fi;
  end
);

InstallMethod(RRestrictedModuleOp, "r-restriction for projective indecomposable modules",
  [IsAlgebraObj, IsHeckePIM, IsList],
  function(H, x, list) local nx;
    x:=RRestrictedModule(H,MakeSpechtOp(x,false),list);
    if x=fail or x=0*x then return x; fi;
    nx:=MakePIMOp(x,false);
    if nx<>fail then return nx; else return x; fi;
  end
);

InstallMethod(RRestrictedModuleOp, "r-restriction for simple modules",
  [IsAlgebraObj, IsHeckeSimple, IsList],
  function(H, x, list) local nx;
    x:=RRestrictedModule(H,MakeSpechtOp(x,false),list);
    if x=fail or x=0*x then return x; fi;
    nx:=MakeSimpleOp(x,false);
    if nx<>fail then return nx; else return x; fi;
  end
);

InstallMethod(RRestrictedModuleOp, "r-restriction for Fock space specht modules",
  [IsAlgebraObj,IsFockSpecht,IsList],
  function(H, x, list) local r;
    if list=[] then return Sum([0..H!.e-1],r->qSRestrictedModule(H,x,1,r));
    elif H!.e=0 then
      Error("Restrict, r-restriction is not defined when e=0.");
    elif ForAny(list,r-> r>=H!.e or r<0) then
      Error("Restrict, r-restriction is defined only when 0<=r<e.\n");
    else
      for r in list do   ## we could do slightly better here
        x:= qSRestrictedModule(H,x,1,r);
      od;
      return x;
    fi;
  end
);

InstallMethod(RRestrictedModuleOp,
  "r-restriction for Fock space projective indecomposable modules",
  [IsAlgebraObj,IsFockPIM,IsList],
  function(H,x,list)
    return MakePIMOp(RRestrictedModule(H,MakeSpechtOp(x,false),list),false);
  end
);

InstallMethod(RRestrictedModuleOp,
  "r-restriction for Fock space simple modules",
  [IsAlgebraObj,IsFockSimple,IsList],
  function(H,x,list)
    return MakeSimpleOp(RRestrictedModule(H,MakeSpechtOp(x,false),list),false);
  end
); # RRestrictedModule

InstallMethod(SInducedModuleOp,"string induction for specht modules",
  [IsAlgebraObj, IsHeckeSpecht, IsList],
  function(H, x, list) local r;
    if x=fail or x=0*x then return x;
    elif Length(list)=1 then
      list:=list[1];
      if list=0 then return Module(H,"Sq",1,[]); fi;
      while list > 0 do
        x:=SInducedModule(H,x,1,1,0);
        list:=list-1;
      od;
      return x;
    elif H!.e=0 then
      Error("SInduce, string induction is not defined when e=0.");
    elif list[2]>H!.e or list[2]<0 then
      Error("SInduce, string induction is defined only when 0<=r<e.\n");
    else return SInducedModule(H, x, H!.e, list[1], list[2]);
    fi;
  end
);

InstallMethod(SInducedModuleOp, "string induction for projective indecomposable modules",
  [IsAlgebraObj, IsHeckePIM, IsList],
  function(H, x, list) local nx;
    x:=SInducedModule(H,MakeSpechtOp(x,false),list);
    if x=fail or x=0*x then return x; fi;
    nx:=MakePIMOp(x,false);
    if nx<>fail then return nx; else return x; fi;
  end
);

InstallMethod(SInducedModuleOp, "string induction for simple modules",
  [IsAlgebraObj, IsHeckeSimple, IsList],
  function(H, x, list) local nx;
    x:=SInducedModule(H,MakeSpechtOp(x,false),list);
    if x=fail or x=0*x then return x; fi;
    nx:=MakeSimpleOp(x,false);
    if nx<>fail then return nx; else return x; fi;
  end
);

InstallMethod(SInducedModuleOp, "string induction for Fock space specht modules",
  [IsAlgebraObj,IsFockSpecht,IsList],
  function(H, x, list) local r;
    if Length(list)=1 then
      list:=list[1];
      if list=0 then return Module(H,"Sq",1,[]); fi;
      while list > 0 do
        x:=Sum([0..H!.e-1],r->qSInducedModule(H,x,1,r));
        list:=list-1;
      od;
      return x;
    elif H!.e=0 then
      Error("SInduce, string induction is not defined when e=0.");
    elif list[2]>H!.e or list[2]<0 then
      Error("SInduce, string induction is defined only when 0<=r<e.\n");
    else return qSInducedModule(H, x, list[1], list[2]);
    fi;
  end
);

InstallMethod(SInducedModuleOp,
  "string induction for Fock space projective indecomposable modules",
  [IsAlgebraObj,IsFockPIM,IsList],
  function(H,x,list)
    return MakePIMOp(SInducedModule(H,MakeSpechtOp(x,false),list),false);
  end
);

InstallMethod(SInducedModuleOp,
  "string induction for Fock space simple modules",
  [IsAlgebraObj,IsFockSimple,IsList],
  function(H,x,list)
    return MakeSimpleOp(SInducedModule(H,MakeSpechtOp(x,false),list),false);
  end
); # SInducedModule

InstallMethod(SRestrictedModuleOp,"string restriction for specht modules",
  [IsAlgebraObj, IsHeckeSpecht, IsList],
  function(H, x, list) local r;
    if x=fail or x=0*x then return x;
    elif Length(list)=1 then
      list:=list[1];
      if list=0 then return Module(H,"Sq",1,[]); fi;
      while list > 0 do
        x:=SRestrictedModule(H,x,1,1,0);
        list:=list-1;
      od;
      return x;
    elif H!.e=0 then
      Error("SRestrict, r-restriction is not defined when e=0.");
    elif list[2]>H!.e or list[2]<0 then
      Error("SRestrict, r-restriction is defined only when 0<=r<e.\n");
    else return SRestrictedModule(H, x, H!.e, list[1], list[2]);
    fi;
  end
);

InstallMethod(SRestrictedModuleOp, "string restriction for projective indecomposable modules",
  [IsAlgebraObj, IsHeckePIM, IsList],
  function(H, x, list) local nx;
    x:=SRestrictedModule(H,MakeSpechtOp(x,false),list);
    if x=fail or x=0*x then return x; fi;
    nx:=MakePIMOp(x,false);
    if nx<>fail then return nx; else return x; fi;
  end
);

InstallMethod(SRestrictedModuleOp, "string restriction for simple modules",
  [IsAlgebraObj, IsHeckeSimple, IsList],
  function(H, x, list) local nx;
    x:=SRestrictedModule(H,MakeSpechtOp(x,false),list);
    if x=fail or x=0*x then return x; fi;
    nx:=MakeSimpleOp(x,false);
    if nx<>fail then return nx; else return x; fi;
  end
);

InstallMethod(SRestrictedModuleOp, "string restriction for Fock space specht modules",
  [IsAlgebraObj,IsFockSpecht,IsList],
  function(H, x, list) local r;
    if Length(list)=1 then
      list:=list[1];
      if list=0 then return Module(H,"Sq",1,[]); fi;
      while list > 0 do
        x:=Sum([0..H!.e-1],r->qSRestrictedModule(H,x,1,r));
        list:=list-1;
      od;
      return x;
    elif H!.e=0 then
      Error("SRestrict, string restriction is not defined when e=0.");
    elif list[2]>H!.e or list[2]<0 then
      Error("SRestrict, string restriction is defined only when 0<=r<e.\n");
    else return qSRestrictedModule(H, x, list[1], list[2]);
    fi;
  end
);

InstallMethod(SRestrictedModuleOp,
  "string restriction for Fock space projective indecomposable modules",
  [IsAlgebraObj,IsFockPIM,IsList],
  function(H,x,list)
    return MakePIMOp(SRestrictedModule(H,MakeSpechtOp(x,false),list),false);
  end
);

InstallMethod(SRestrictedModuleOp,
  "string restriction for Fock space simple modules",
  [IsAlgebraObj,IsFockSimple,IsList],
  function(H,x,list)
    return MakeSimpleOp(SRestrictedModule(H,MakeSpechtOp(x,false),list),false);
  end
); #SRestrictedModule

InstallMethod(RInducedModuleOp,
  "toplevel r-induction",
  [IsAlgebraObjModule,IsList],
  function(x,list)
    return RInducedModuleOp(x!.H,x,list);
  end
); #RInducedModule

InstallMethod(RRestrictedModuleOp,
  "toplevel r-restriction",
  [IsAlgebraObjModule,IsList],
  function(x,list)
    return RRestrictedModuleOp(x!.H,x,list);
  end
); #RRestrictedModule

InstallMethod(SInducedModuleOp,
  "toplevel string induction",
  [IsAlgebraObjModule,IsList],
  function(x,list)
    return SInducedModuleOp(x!.H,x,list);
  end
); #SInducedModule

InstallMethod(SRestrictedModuleOp,
  "toplevel string restriction",
  [IsAlgebraObjModule,IsList],
  function(x,list)
    return SRestrictedModuleOp(x!.H,x,list);
  end
); #SRestrictedModule

## Q-INDUCTION AND Q-RESTRICTION ###############################################

## notice that we can't pull the trick to induce all residues at once
## that we used in InducedModule() etc. as we have to keep track of the
## number of addable and removable nodes of each residue. Rather than
## do this we just call this function as many times as necessary.
InstallMethod(qSInducedModule,"q-induction for modules",
  [IsAlgebraObj,IsAlgebraObjModule,IsInt,IsInt],
  function(H,x,s,r) local coeffs, parts, y, z, qsinduced, v;
    v:=H!.Indeterminate;

    # add n nodes of residue r to the partition y from the i-th row down
    # here exp is the exponent of the indeterminate
    qsinduced:=function(y, n, r, i, exp) local ny, j, z;
     ny:=[];
     for j in [i..Length(y)-n+1] do
       if y[j]>0 and r=(y[j]-j) mod H!.e and (j=Length(y) or y[j]>y[j+1])
       then exp:=exp-n;               ## removable node of residue r
       elif r=(y[j] - j + 1) mod H!.e then
         if j=1 or y[j] < y[j-1] then ## addable node of residue r
           z:=StructuralCopy(y);
           z[j]:=z[j] + 1; # only one node of residue r can be added
           if n=1 then Add(ny, [exp,z] );   # no more nodes to add
           else Append(ny, qsinduced(z, n-1, r, j+1, exp));
           fi;
           exp:=exp+n;
         fi;
       fi;
     od;
     return ny;
    end;

    coeffs:=[]; parts:=[];
    for y in [1..Length(x!.parts)] do
      if r=( -Length(x!.parts[y]) mod H!.e) then # add a node to the bottom
        z:=StructuralCopy(x!.parts[y]);
        Add(z,0);
        for z in qsinduced(z,s,r,1,0)  do
          if z[1]=0 then Add(coeffs, x!.coeffs[y]);
          else Add(coeffs, x!.coeffs[y]*v^z[1]);
          fi;
          if z[2][Length(z[2])]=0 then Unbind(z[2][Length(z[2])]); fi;
          Add(parts, z[2]);
        od;
      else
        for z in qsinduced(x!.parts[y],s,r,1,0) do
          if z[1]=0 then Add(coeffs, x!.coeffs[y]);
          else Add(coeffs, x!.coeffs[y]*v^z[1]);
          fi;
          Add(parts, z[2]);
        od;
      fi;
    od;

    if coeffs=[] then return Module(H,x!.module,0,[]);
    else return Collect(H,x!.module, coeffs, parts);
    fi;
  end
);  # qSInducedModule

## string-restriction: remove m r's from each partition in x
## ** should allow restricting s-times also
InstallMethod(qSRestrictedModule,"q-restriction for modules",
  [IsAlgebraObj,IsAlgebraObjModule,IsInt,IsInt],
  function(H, x, s, r) local coeffs, parts, z, y, i, e, exp, v, qrestricted;
    v:=H!.Indeterminate;

    qrestricted:=function(y, n, r, i, exp) local ny, j, z;
      ny:=[];
      for j in [i,i-1..n] do
        if y[j]>0 and r=(y[j]+1-j) mod H!.e and (j=1 or y[j]<y[j-1]) then
           exp:=exp-n;                 ## an addable node of residue r
        elif r=(y[j] - j) mod H!.e then   ## removable node of residue r
          if j=Length(y) or y[j] > y[j+1] then
            z:=StructuralCopy(y);
            z[j]:=z[j]-1;
            if z[j]=0 then Unbind(z[j]); fi;
            if n=1 then Add(ny, [exp,z]); # no mode nodes to remove
            else Append(ny, qrestricted(z, n-1, r, j-1, exp));
            fi;
            exp:=exp+n;
          fi;
        fi;
      od;
      return ny;
    end;

    e:=x!.H!.e;
    coeffs:=[]; parts:=[];
    if s=0 then return Module(H,"Sq",1,[]); fi;
    for y in [1..Length(x!.parts)] do
      if -Length(x!.parts[y]) mod H!.e = r then exp:=-s; else exp:=0; fi;
      for z in qrestricted(x!.parts[y], s, r, Length(x!.parts[y]), exp) do
        if z[1]=0 then Add(coeffs, x!.coeffs[y]);
        else Add(coeffs, x!.coeffs[y]*v^z[1]);
        fi;
        Add(parts, z[2]);
      od;
    od;
    if parts=[] then return Module(H,"Sq",0,[]);
    else return Collect(H, "Sq", coeffs, parts);
    fi;
  end
);  # qSRestrict

## DECOMPOSITION MATRICES ######################################################

## The following directory is searched by ReadDecompositionMatrix()
## when it is looking for decomposition matrices. By default, it points
## to the current directory (if set, the current directory is not
## searched).
if not IsBound(SpechtDirectory) then SpechtDirectory:=Directory("."); fi;

## This variable is what is used in the decomposition matrices files saved
## by SaveDecompositionMatrix() (and also the variable which contains them
## when they are read back in).
A_Specht_Decomposition_Matrix:=fail;

## Finally, we can define the creation function for decomposition matrices
## (note that NewDM() does not add the partition labels to the decomp.
## matrix; this used to be done here but now happens in PrintDM() because
## crystallized matrices may never be printed and this operation is
## expensive).
## **NOTE: we assume when extracting entries from d that d.rows is
## ordered lexicographically. If this is not the case then addition
## will not work properly.
InstallOtherMethod(DecompositionMatrix,"creates a new decomposition matrix",
  [IsAlgebraObj,IsList,IsList,IsBool],
  function(H, rows, cols, decompmat) local d;
    d := rec(d:=[],      # matrix entries
       rows:=rows, # matrix rows
       cols:=cols, # matrix cols
       inverse:=[], dimensions:=[], ## inverse matrix and dimensions
       H:=H
    );

    if decompmat then
      return Objectify(DecompositionMatrixType,d);
    else
      return Objectify(CrystalDecompositionMatrixType,d);
    fi;
  end
); # DecompositonMatrix

InstallMethod(\=, "for decomposition matrices",
  [IsDecompositionMatrix,IsDecompositionMatrix],
  function(d1,d2)
    return d1!.d=d2!.d and d1!.cols = d2!.cols and d1!.rows = d2!.rows;
  end
); # Equal matrices

InstallMethod(CopyDecompositionMatrix, "for decomposition matrices",
  [IsDecompositionMatrix],
  function(d) local newd;
    newd:=DecompositionMatrix(
      d!.H, StructuralCopy(d!.rows), StructuralCopy(d!.cols),
      not IsCrystalDecompositionMatrix(d)
      );
    newd!.d:=StructuralCopy(d!.d);
    newd!.inverse:=StructuralCopy(d!.inverse);
    newd!.dimensions:=StructuralCopy(d!.dimensions);
    return newd;
  end
); # CopyDecompositionMatrix

## Used by SaveDecomposition matrix to update CrystalMatrices[]
InstallMethod(Store,"for crystal decomposition matrices",
  [IsCrystalDecompositionMatrix,IsInt],
  function(d,n) d!.H!.CrystalMatrices[n]:=d; end
);

## Used by SaveDecomposition matrix to update DecompositionMatrices[]
InstallMethod(Store,"for decomposition matrices",[IsDecompositionMatrix,IsInt],
  function(d,n) d!.H!.DecompositionMatrices[n]:=d; end
);

## This also needs to be done for crystal matrices (this wasn't
## put in above because normal decomposition matrices don't
## specialise.
InstallMethod(Specialized,"specialise CDM",
  [IsCrystalDecompositionMatrix,IsInt],
  function(d,a) local sd, c, p, coeffs, v;
    v:=d!.H!.Indeterminate;
    sd:=DecompositionMatrix(d!.H,d!.rows,d!.cols,true);
    for c in [1..Length(d!.cols)] do
      if IsBound(d!.d[c]) then
        sd!.d[c]:=rec();
        coeffs:=List(d!.d[c].coeffs*v^0,p->Value(p,a));
        p:=List(coeffs,p->p<>0);
        if true in p then
          sd!.d[c]:=rec(coeffs:=ListBlist(coeffs,p),
                        parts:=ListBlist(d!.d[c].parts,p) );
        else sd!.d[c]:=rec(coeffs:=[0], parts:=[ [] ] );
        fi;
      fi;
      if IsBound(d!.inverse[c]) then
        coeffs:=List(d!.inverse[c].coeffs*v^0,p->Value(p,a));
        p:=List(coeffs,p->p<>0);
        if true in p then
          sd!.inverse[c]:=rec(coeffs:=ListBlist(coeffs,p),
                          parts:=ListBlist(d!.inverse[c].parts,p) );
        else sd!.d[c]:=rec(coeffs:=[0], parts:=[ [] ] );
        fi;
      fi;
    od;
    return sd;
  end
);

InstallMethod(Specialized,"specialise CDM",
  [IsCrystalDecompositionMatrix], a->Specialized(a,1)
);

## Specialization taking Xq -> X
InstallMethod(Specialized,"specialise Fock module",
  [IsFockModule,IsInt],
  function(x,a) local coeffs, c;
    coeffs:=List(x!.coeffs*(x!.H!.Indeterminate)^0,c->Value(c,a));
    c:=List(coeffs,c->c<>0);
    if true in c then return Module(x!.H,x!.module{[1]},
                                ListBlist(coeffs,c),ListBlist(x!.parts,c));
    else return Module(x!.H,x!.module{[1]},0,[]);
    fi;
  end
);

InstallMethod(Specialized,"specialise Fock module",
  [IsFockModule], x->Specialized(x,1)
); # Specialized

## writes D(mu) as a linear combination of S(nu)'s if possible
InstallMethod(Invert, "for a decomposition matrix and a partition",
  [IsDecompositionMatrix,IsList],
  function(d,mu) local c, S, D, inv, smu, l, tmp;
    if IsCrystalDecompositionMatrix(d) then S:="Sq"; D:="Dq";
    else S:="S"; D:="D";
    fi;

    c:=Position(d!.cols,mu);
    if c=fail then return fail;
    elif IsBound(d!.inverse[c]) then
      return Module(d!.H,S,d!.inverse[c].coeffs,
                    List(d!.inverse[c].parts,l->d!.cols[l]));
    fi;

    inv:=Module(d!.H,S,1,mu);
    tmp:=MakeSimpleOp(d,inv);
    if tmp=fail then
      smu:=fail;
    else
      smu:=Module(d!.H,D,1,mu)-tmp;
    fi;

    while smu<>fail and smu<>0*smu do
      inv:=inv+Module(d!.H,S,smu!.coeffs[1],smu!.parts[1]);
      tmp:=MakeSimpleOp(d, Module(d!.H,S,-smu!.coeffs[1],smu!.parts[1]));
      if tmp=fail then
        smu:=fail;
      else
        smu:=smu+tmp;
      fi;
    od;
    if smu=fail then return fail; fi;

    d!.inverse[c]:=rec(coeffs:=inv!.coeffs,
                  parts:=List(inv!.parts,l->Position(d!.cols,l)));
    return inv;
  end
);

#P This function adds the column for Px to the decomposition matrix <d>.
## if <checking>=true then Px is checked against its image under the
## Mullineux map; if this image is not there then we also insert it
## into <d>.
InstallMethod(AddIndecomposable,"fill out entries of decomposition matrix",
  [IsDecompositionMatrix,IsAlgebraObjModule,IsBool],
  function(d,Px,checking) local mPx, r, c;
    c:=Position(d!.cols, Px!.parts[Length(Px!.parts)]);
    if checking then
      ## first look to see if <Px> already exists in <d>
      if IsBound(d!.d[c]) then ## Px already exists
        Print("# AddIndecomposable: overwriting old value of P(",
              TightStringList(Px!.parts[Length(Px!.parts)]),") in <d>\n");
        Unbind(d!.inverse);       # just in case these were bound
        Unbind(d!.dimensions);
      fi;
      ## now looks at the image of <Px> under Mullineux
      if (not IsSchur(Px!.H) and IsERegular(Px!.H!.e,Px!.parts[Length(Px!.parts)]))
      and Px!.parts[Length(Px!.parts)]<>ConjugatePartition(Px!.parts[1]) then
        mPx:=MullineuxMap(Px);
        if IsBound(d!.d[Position(d!.cols,mPx!.parts[Length(Px!.parts)])])
        and MakePIMSpechtOp(d,mPx!.parts[Length(Px!.parts)]) <> mPx then
          Print("# AddIndecomposable(<d>,<Px>), WARNING: P(",
                TightStringList(Px!.parts[Length(Px!.parts)]), ") and P(",
                TightStringList(mPx!.parts[Length(Px!.parts)]),
                ") in <d> are incompatible\n");
        else
          Print("# AddIndecomposable(<d>,<Px>): adding MullineuxMap(<Px>) ",
                "to <d>\n");
          AddIndecomposable(d,mPx,false);
        fi;
      fi;
    fi;      # end of check
    d!.d[c]:=rec(coeffs:=Px!.coeffs,
                parts:=List(Px!.parts,r->Position(d!.rows,r)) );
  end
);

#F This function checks to see whether Px is indecomposable using the
## decomposition matrix d. The basic idea is to loop through all of the
## (e-regular, unless IsSpecht=false) projectives Py in Px such that Px-Py
## has non-negative coefficients and to then apply the q-Schaper theorem
## and induce simples, together with a few other tricks to decide
## whether ir not Py slits off. If yes, then we move on; if we can't
## decide (or don't know the value of Py), we spit the dummy and return
## false, printing a "message from our sponsor" explaining why we failed
## if called from outside InducedModule(). If we can count for all
## projectives then we return true. Note that in this case the value
## of Px may have changed, but we update the original value of px=Px
## before leaving (whether we return false or true).
InstallMethod(IsNewIndecomposableOp, "checks whether the given module is indecomposable",
  [IsAlgebraObj, IsDecompositionMatrix,IsAlgebraObjModule,IsDecompositionMatrix,IsList],
  function(H,d,px,oldd,Mu)
    local Px,nu,regs,Py,y,z,a,b,n,mu,m,M,Message;

    if px=fail and px=0*px then return false; fi;

    if Mu=[fail] then Message:=Ignore;
    else Message:=Print;
    fi;

    Px:=px; Py:=true; # strip PIMs from the top of Px
    while Py<>fail and Px<>0*Px do
      Py:=MakePIMSpechtOp(d,Px!.parts[Length(Px!.parts)]);
      if Py<>fail then Px:=Px-Px!.coeffs[Length(Px!.parts)]*Py; fi;
    od;
    if Px=0*Px then
      Message("# This module is a sum of known indecomposables.\n");
      return false;
    fi;
    if IntegralCoefficients(Px / Px!.coeffs[Length(Px!.parts)]) then
      Px:=Px / Px!.coeffs[Length(Px!.parts)];
   fi;

    regs:=Obstructions(d,Px);

    if Mu<>[] then regs:=Filtered(regs,mu->mu<Mu); fi;
    regs:=Obstructions(d,Px);
    for y in regs do   ## loop through projectives that might split

      if not IsSchur(H) and MullineuxMap(H!.e,ConjugatePartition(Px!.parts[1]))
      <>Px!.parts[Length(Px!.parts)] then
        Py:=true;  ## strip any known indecomposables off the bottom of Px
        while Py<>fail and Px<>0*Px do
          Py:=MakePIMSpechtOp(d,ConjugatePartition(Px!.parts[1]));
          if Py<>fail and IsERegular(Py!.H!.e,Py!.parts[Length(Py!.parts)]) then
            Px:=Px-Px!.coeffs[1]*MullineuxMap(Py);
          else Py:=fail;
          fi;
        od;
      fi;

      m:=0;            ## lower and upper bounds on decomposition the number
      if Px=0*Px then M:=0;
      else M:=Coefficient(Px,y)/Px!.coeffs[Length(Px!.parts)];
      fi;
      if M<>0 then Py:=MakePIMSpechtOp(d,y); fi;
      if not ( m=M or Px!.parts[Length(Px!.parts)]>=y ) then
        if Py=fail then
          Message("# The multiplicity of S(", TightStringList(y),
            ") in <Px> is zero; however, S(", TightStringList(y),
            ") is not known\n");
          px!.coeffs:=Px!.coeffs; px!.parts:=Px!.parts; return false;
        else
          Px:=Px-Coefficient(Px,y)*Py;
          if not PositiveCoefficients(Px) then BUG("IsNewIndecomposable",1);fi;
          M:=0;
        fi;
      fi;

      if Px<>0*Px and IntegralCoefficients(Px/Px!.coeffs[Length(Px!.parts)]) then
        Px:=Px/Px!.coeffs[Length(Px!.parts)];
      fi;

      ## remember that Px.coeffs[Length(Px.parts)] could be greater than 1
      if m<>M and (Coefficient(Px,y) mod Px!.coeffs[Length(Px!.parts)])<>0 then
        if Py=fail then
          Message("# <Px> is not indecomposable, as at least ",
            Coefficient(Px,y) mod Px!.coeffs[Length(Px!.parts)], " copies of P(",
            TightStringList(y), ") split off.\n# However, P(",
            TightStringList(y), ") is not known\n");
          px!.coeffs:=Px!.coeffs; px!.parts:=Px!.parts; return false;
        else
          ## this is at least projective; perhaps more still come off though
          Px:=Px-(Coefficient(Px,y) mod Px!.coeffs[Length(Px!.parts)])*Py;
          if not PositiveCoefficients(Px) then BUG("IsNewIndecomposable",2);fi;
          if IntegralCoefficients(Px/Px!.coeffs[Length(Px!.parts)]) then
            Px:=Px/Px!.coeffs[Length(Px!.parts)];
          fi;
        fi;
      fi;

      ## At this point the coefficient of Sx in Px divides the coefficient of
      ## Sy in Px. If Py splits off then it does so in multiples of
      ## (Px:Sx)=Px.coeffs[Length(Px.parts)].
      if m<>M and (Py=fail
        or PositiveCoefficients(Px-Px!.coeffs[Length(Px!.parts)]*Py))
      then
        ## use the q-Schaper theorem to test whether Sy is contained in Px
        M:=Minimum(M,InnerProduct(Px/Px!.coeffs[Length(Px!.parts)],
                                  Schaper(H, y)));
        if M=0 then # NO!
          ## Px-(Px:Sy)Py is still projective so subtract Py if it is
          ## known. If Py=false (ie. not known), then at least we know
          ## that Px is not indecomposable, even though we couldn't
          ## calculate Py.
          if Py=fail then
            Message("# The multiplicity of S(", TightStringList(y),
              ") in P(", TightStringList(Px!.parts[Length(Px!.parts)]),
             ") is zero;\n#  however, P(", TightStringList(y),
             ") is not known.\n");
            px!.coeffs:=Px!.coeffs; px!.parts:=Px!.parts; return false;
          else
            Px:=Px-Coefficient(Px,y)*Py;
            if not PositiveCoefficients(Px) then BUG("IsNewIndecomposable",3);
            fi;
          fi;
        elif Px!.coeffs[Length(Px!.parts)]<>Coefficient(Px,y) then
          ## We know that (Px:Sy)>=m>0, but perhaps some Py's still split off

          m:=1;
          if m=M then
            if Coefficient(Px,y)<>m*Px!.coeffs[Length(Px!.parts)] then
              if Py=fail then
                Message("# The multiplicity of S(", TightStringList(y),
                    ") in P(",TightStringList(Px!.parts[Length(Px!.parts)]),
                    ") is ", m, "however, P(",TightStringList(y),
                    ") is unknown.\n");
                px!.coeffs:=Px!.coeffs; px!.parts:=Px!.parts;
                return false;
              else
                Px:=Px-(Coefficient(Px,y)-Px!.coeffs[Length(Px!.parts)]*m)*Py;
                if not PositiveCoefficients(Px) then
                  BUG("IsNewIndecomposable",4);
                fi;
              fi;
            fi;
          fi;

          if m<>M then
            ## see if we can calculate this decomposition number (this uses
            ## row and column removal)
            a:=DecompositionNumber(Px!.H, y, Px!.parts[Length(Px!.parts)]);
            if a<>fail then
              if Px!.coeffs[Length(Px!.parts)]*a=Coefficient(Px,y) then m:=a; M:=a;
              elif Py<>fail then
                ## precisely this many Py's come off
                Px:=Px-(Coefficient(Px,y)-Px!.coeffs[Length(Px!.parts)]*a)*Py;
                m:=a; M:=a; # upper and lower bounds are equal
                if not PositiveCoefficients(Px) then
                  BUG("IsNewIndecomposable",5);
                fi;
              fi;
            fi;
          fi;

          if m<>M and Py=fail then ## nothing else we can do
            Message("# The multiplicity of S(", TightStringList(y),
                    ") in P(",TightStringList(Px!.parts[Length(Px!.parts)]),
                    ") is at least ", m, " and at most ", M,
                    ";\n# however, P(",TightStringList(y),") is unknown.\n");
            px!.coeffs:=Px!.coeffs; px!.parts:=Px!.parts;
            return false;
          fi;

          if m<>M then
            ## Maybe the Mullineux map can lower M...
            M:=Minimum(M,Coefficient(Px,
                 Py!.parts[Length(Py!.parts)]/Px!.coeffs[Length(Px!.parts)]));
            if Coefficient(Px,y)>M*Px!.coeffs[Length(Px!.parts)] then
              Px:=Px-(Coefficient(Px,y)-M*Px!.coeffs[Length(Px!.parts)])*Py;
            fi;
            while m<M and not
            PositiveCoefficients(Px-(Px!.coeffs[Length(Px!.parts)]*(M-m))*Py) do
              m:=m+1;
            od;
          fi;

          ## finally, we take a look at inducing the simples in oldd
          if m<>M and oldd<>fail then
            if not IsBound(oldd!.simples) then ## we have to first induce them
              oldd!.simples:=[];
              if IsBound(oldd!.ind) then       ## defined in InducedModule()
               for mu in [1..Length(oldd!.cols)] do
                  a:=Invert(oldd,oldd!.cols[mu]);
                  if a<>fail then
                    for n in [1..H!.e] do  ## induce simples of oldd
                      z:=Sum([1..Length(a!.coeffs)],b->a!.coeffs[b]
                           *oldd!.ind[Position(oldd!.rows,a!.parts[b])][n]);
                      if z<>0*z then Add(oldd!.simples,z);fi;
                    od;
                  fi;
                od;
              else   ## do everything by hand
                for mu in [1..Length(oldd!.cols)] do
                  a:=Invert(oldd,oldd!.cols[mu]);
                  if a<>fail then
                    for n in [0..H!.e-1] do
                      z:=RInducedModule(H,a,H!.e,n);
                      if z<>0*z then Add(oldd!.simples,z); fi;
                    od;
                  fi;
                od;
              fi;
            fi;
            mu:=Length(oldd!.simples);
            while mu >0 and m<M do
              z:=oldd!.simples[mu];
              if y=regs[Length(regs)]
              or Lexicographic(z!.parts[1],Py!.parts[Length(Py!.parts)]) then
                a:=InnerProduct(z,Py);
                if a<>0 then
                  b:=InnerProduct(z,Px)/Px!.coeffs[Length(Px!.parts)];
                  m:=Maximum(m,M-Int(b/a));
                fi;
              fi;
              mu:=mu-1;
            od;
            if Coefficient(Px,y)>M*Px!.coeffs[Length(Px!.parts)] then
              Px:=Px-(Coefficient(Px,y)-M*Px!.coeffs[Length(Px!.parts)])*Py;
            fi;
          fi;
          if m<>M then ## nothing else we can do
            px!.coeffs:=Px!.coeffs; px!.parts:=Px!.parts;
            Message("# The multiplicity of S(", TightStringList(y),
                ") in P(",TightStringList(Px!.parts[Length(Px!.parts)]),
                ") is at least ", m, " and at most ", M, ".\n");
            return false;
          fi;
        fi;   ## q-Schaper test
      fi;
    od;
    if Px=0*Px then
      Message("# This module is a sum of known indecomposables.\n");
      return false;
    elif Px!.coeffs[Length(Px!.parts)]<>1 then BUG("IsNewIndecomposable",6);
    else px!.coeffs:=Px!.coeffs; px!.parts:=Px!.parts; return true;
   fi;
  end
);  # IsNewIndecomposable

#P A front end to d.operations.AddIndecomposable. This function adds <Px>
## into the decomposition matrix <d> and checks that it is compatible with
## its image under the Mullineux map, if this is already in <d>, and
## inserts it if it is not.
InstallMethod(AddIndecomposable,"fill out entries of decomposition matrix",
  [IsDecompositionMatrix,IsHeckeSpecht],
  function(d, Px)
    if Position(d!.cols, Px!.parts[Length(Px!.parts)])=fail then
      Print("# The projective P(",TightStringList(Px!.parts[Length(Px!.parts)]),
            ") is not listed in <D>\n");
    else AddIndecomposable(d,Px,true);
    fi;
  end
);# AddIndecomposable

## This function will read the file "n" if n is a string; otherwise
## it will look for the relevant file for H(Sym_n). If crystal=true
## this will be a crystal decomposition matrix, otherwise it will be
## a 'normal' decomposition matrix. When IsInt(n) the lists
## CrystalMatrices[] and DecompositionMatrices[] are also checked, and
## the crystal decomposition matrix is calculated if it is not found
## and crystal=true (and IsInt(n)).
## Question: if crystal=false but IsBound(CrystalMatrices[n]) should
## we still try and read e<H.e>p0.n or specialize CrystalMatrices[n]
## immediately. We try and read the file first...
InstallMethod(ReadDecompositionMatrix, "load matrix from library",
  [IsAlgebraObj,IsInt,IsBool],
  function(H, n, crystal) local d, file, SpechtDirectory;

    if not IsBound(SpechtDirectory) then SpechtDirectory:=""; fi;

    if crystal then
      if IsBound(H!.CrystalMatrices[n]) then
        d:=H!.CrystalMatrices[n];
        if d=fail and H!.info.SpechtDirectory=SpechtDirectory then
          return fail;
        elif d<>fail and ForAll([1..Length(d!.cols)],c->IsBound(d!.d[c]))
        then return d;
        fi;
      fi;
      file:=Concatenation("e",String(H!.e),"crys.",String(n));
    else
      file:=Concatenation(H!.HeckeRing,".",String(n));
    fi;
    d:=ReadDecompositionMatrix(H,file,crystal);
    if crystal then H!.CrystalMatrices[n]:=d; fi;
    return d;
  end
);

InstallMethod(ReadDecompositionMatrix, "load matrix from library",
  [IsAlgebraObj,IsString,IsBool],
  function(H, str, crystal)
    local msg, file, M, d, c, parts, coeffs, p, x, r, cm, rm;

    A_Specht_Decomposition_Matrix:=fail; ## just in case

    file:=Filename([Directory("."),SpechtDirectory,H!.info.Library],str);
    if crystal then
      msg:="ReadCrystalMatrix-";
    else
      msg:="ReadDecompositionMatrix-";
    fi;

    d:=fail;

    if file <> fail then Read(file); fi;

    if A_Specht_Decomposition_Matrix<>fail then   ## extract matrix from M
      M:=A_Specht_Decomposition_Matrix;
      A_Specht_Decomposition_Matrix:=fail;
      r:=Set(M.rows); c:=Set(M.cols);
      if not IsSchur(H) and r=c then
        d:=DecompositionMatrix(H,r,
              Filtered(c,x->IsERegular(H!.e,x)),not IsBound(M.crystal));
      elif not IsSchur(H) then
        d:=DecompositionMatrix(H,r,c,not IsBound(M.crystal));
      else
        d:=DecompositionMatrix(H,r,r,not IsBound(M.crystal));
      fi;
      if IsSet(M.rows) and IsSet(M.cols) then ## new format
        if IsBound(M.matname) then d!.matname:=M.matname; fi;
        for c in [1..Length(d!.cols)] do
          cm:=Position(M.cols, d!.cols[c]);
          if cm<>fail and IsBound(M.d[cm]) then
            x:=M.d[cm];
            parts:=[]; coeffs:=[];
            for rm in [1..Length(x)/2] do
              r:=Position(d!.rows,M.rows[x[rm+Length(x)/2]]);
              if r<>fail then
                Add(parts,r);
                if IsInt(x[rm]) then Add(coeffs,x[rm]);
                else
                  p:=LaurentPolynomialByCoefficients(
                    FamilyObj(1),
                    x[rm]{[2..Length(x[rm])]},x[rm]);
                  Add(coeffs,p);
                fi;
              fi;
            od;
            if parts<>[] then   ## paranoia
              SortParallel(parts,coeffs);
              d!.d[c]:=rec(parts:=parts,coeffs:=coeffs);
           fi;
          fi;
        od;
      else  ## old format
        d!.d:=List(c, r->rec(coeffs:=[], parts:=[]));
        ## next, we unravel the decomposition matrix
        for rm in [1..Length(M.rows)] do
          r:=Position(d!.rows,M.rows[rm]);
          if r<>fail then
            x:=1;
            while x<Length(M.d[rm]) do
              c:=Position(d!.cols,M.cols[M.d[rm][x]]);
              if c<>fail then
                Add(d!.d[c].coeffs, M.d[rm][x+1]);
                Add(d!.d[c].parts, r);
              fi;
              x:=x+2;
            od;
          fi;
        od;
        for c in [1..Length(d!.d)] do
          if d!.d[c].parts=[] then Unbind(d!.d[c]);
          else SortParallel(d!.d[c].parts, d!.d[c].coeffs);
          fi;
        od;
      fi;
    fi;
    return d;
  end
); # ReadDecompositionMatrix

## Look up the decomposition matrix in the library files and in the
## internal lists and return it if it is known.
## NOTE: this function does not use the crystal basis to calculate the
## decomposition matrix because it is called by the various conversion
## functions X->Y which will only need a small part of the matrix in
## general. The function FindDecompositionMatrix() also uses the crystal
## basis.
InstallMethod(KnownDecompositionMatrix,"looks for a known decomposition matrix",
  [IsAlgebraObj,IsInt],
  function(H,n)
    local d, x, r, c;

    if IsBound(H!.DecompositionMatrices[n]) then
      d:=H!.DecompositionMatrices[n];
      if ( d<>fail and ForAll([1..Length(d!.cols)],c->IsBound(d!.d[c])) )
      or H!.info.SpechtDirectory=SpechtDirectory then return d;
      elif H!.info.SpechtDirectory<>SpechtDirectory then
        for x in [1..Length(H!.DecompositionMatrices)] do
          if IsBound(H!.DecompositionMatrices[x]) and
          H!.DecompositionMatrices[x]=fail then
            Unbind(H!.DecompositionMatrices[x]);
          fi;
        od;
      fi;
    fi;
    d:=ReadDecompositionMatrix(H,n,false);

    ## next we look for crystal matrices
    if d=fail and not IsSchur(H) and IsZeroCharacteristic(H) then
      d:=ReadDecompositionMatrix(H,n,true);
      if d<>fail then d:=Specialized(d); fi;
    fi;

    if d=fail and n<2*H!.e then
      ## decomposition matrix can be calculated
      r:=Partitions(n);
      if not IsSchur(H) then c:=ERegularPartitions(H!.e,n);
      else c:=r;
      fi;
      d:=DecompositionMatrix(H, r, c, true);

      for x in [1..Length(d!.cols)] do
        if IsECore(H,d!.cols[x]) then    ## an e-core
          AddIndecomposable(d,
            Module(H,"S",1,d!.cols[x]),false);
        elif IsSimpleModule(H,d!.cols[x]) then ## first e-weight 1 partition
          c:=Module(H,"S",1,ECore(H,d!.cols[x]) )
                   *HeckeOmega(H,"S",H!.e);
          for r in [2..Length(c!.parts)] do
            AddIndecomposable(d,
                  Module(H,"S",[1,1],c!.parts{[r-1,r]}),false);
          od;
          if IsSchur(H) then
            AddIndecomposable(d,
              Module(H,"S",1,c!.parts[1]),false);
          fi;
        fi;
      od;
    elif IsBound(H!.DecompositionMatrices[n]) then ## partial answer only
      return H!.DecompositionMatrices[n];
    fi;
    H!.DecompositionMatrices[n]:=d;
    return d;
  end
);  # KnownDecompositionMatrix

## almost identical to KnownDecompositionMatrix except that if this function
## fails then the crystalized decomposition matrix is calculated.
InstallMethod(FindDecompositionMatrix,"find or calculate CDM",
  [IsAlgebraObj,IsInt],
  function(H,n) local d,c;
    d:=KnownDecompositionMatrix(H,n);
    if d=fail and not IsSchur(H) and IsZeroCharacteristic(H) then
      d:=DecompositionMatrix(H,
              Partitions(n),ERegularPartitions(H!.e,n),false);
      for c in [1..Length(d!.cols)] do
        if not IsBound(d!.d[c]) then
          AddIndecomposable(d,FindPq(H,d!.cols[c]),false);
        fi;
      od;
      H!.CrystalMatrices[n]:=d;
      d:=Specialized(d);
      H!.DecompositionMatrices[n]:=d;
    fi;
    return d;
  end
); # FindDecompositionMatrix

## Crystal basis elements ######################################################

#########################################################################
## Next, for fields of characteristic 0, we implement the LLT algorithm.
## Whenever a crystal basis element of the Fock space is calculated we
## store it in the relevant decomposition matrix n CrystalMatrices[].
## The actual LLT algorithm is contained in the function Pq (should
## really be called Pv...), but there are also functions Sq and Dq as
## well. These functions work as follows:
##   FindPq(mu) -> sum_nu d_{nu,mu} S(nu)
##   FindSq(mu) -> sum_nu d_{mu,nu} D(nu)
##   FindDq(mu) -> sum_nu d_{mu,nu}^{-1} S(nu)
## The later two functions will call Pq() as needed. The "modules" x
## returned by these functions have x.module equal to "Pq", "Sq" and
## "Dq" respectively to distinguish them from the specialized versions.
## Accordingly we need H.operations.Xq for X = S, P, and D; these are
## defined after Pq(), Sq(), and Dq() (which they make use of).

## Retrieves or calculates the crystal basis element Pq(mu)
InstallMethod(FindPq,"finds the crystal basis element Pq(mu)",
  [IsAlgebraObj,IsList],
  function(H,mu) local  n, c, CDM, i, r, s, x, v, val,coeffs,tmp;
    v:=H!.Indeterminate;

    if mu=[] then return Module(H,"Sq",v^0,[]); fi;
    n:=Sum(mu);

    ## first we see if we have already calculated Pq(mu)
    if not IsBound(H!.CrystalMatrices[n]) or H!.CrystalMatrices[n]=fail then
      x:=ReadDecompositionMatrix(H,n,true);
      if x=fail then
        H!.CrystalMatrices[n]:=DecompositionMatrix(H,
          Partitions(n), ERegularPartitions(H!.e,n), false);
      fi;
    fi;
    CDM:=H!.CrystalMatrices[n];
    c:=Position(CDM!.cols,mu);
    if IsBound(CDM!.d[c]) then
      return Module(H,"Sq",CDM!.d[c].coeffs,
                 List(CDM!.d[c].parts, s->CDM!.rows[s]) );
    fi;

    if IsECore(H!.e,mu) then
      x:=Module(H,"Sq",v^0,mu);
    elif EWeight(H!.e,mu)=1 then
      x:=Module(H,"Sq",v^0,ECore(H!.e,mu))
                 * HeckeOmega(H,"Sq",H!.e);
      r:=Position(x!.parts,mu);
      if r=1 then
        x!.parts:=x!.parts{[1]};
        x!.coeffs:=[v^0];
      else
        x!.parts:=x!.parts{[r-1,r]};
        x!.coeffs:=[v,v^0];
      fi;
    else  ## we calculate Pq(mu) recursively using LLT

      ## don't want to change the original mu
      mu:=StructuralCopy(mu);
      i:=1;
      while i<Length(mu) and mu[i]=mu[i+1] do
        i:=i+1;
      od;
      r:=(mu[i]-i) mod H!.e;
      mu[i]:=mu[i]-1;
      s:=1;
      i:=i+1;
      while i<=Length(mu) do
        while i<>Length(mu) and mu[i]=mu[i+1] do
          i:=i+1;
        od;
        if r=(mu[i]-i) mod H!.e then
          s:=s+1;
          mu[i]:=mu[i]-1;
          i:=i+1;
        else i:=Length(mu)+1;
        fi;
      od;
      if mu[Length(mu)]=0  then Unbind( mu[Length(mu)] ); fi;
      x:=qSInducedModule(H, FindPq(H,mu), s, r);
      n:=1;
      while n<Length(x!.parts) do
        tmp:=CoefficientsOfLaurentPolynomial(x!.coeffs[Length(x!.parts)-n]);
        if tmp[2]>0 then n:=n+1;
        else
          r := StructuralCopy( x!.coeffs[Length(x!.parts)-n] );
          tmp:=ShallowCopy(CoefficientsOfLaurentPolynomial(r));
          mu:=x!.parts[Length(x!.parts)-n];
          if Length(tmp[1]) < 1-tmp[1] then
            tmp[1]:=Concatenation(tmp[1],
              List([1..Length(tmp[1])-1-tmp[2]], i->0));
          fi;
          tmp[1]:=tmp[1]{[1..1-tmp[2]]};
          tmp[1]:=Concatenation(tmp[1], Reversed(tmp[1]{[1..-tmp[2]]}));
          r:=LaurentPolynomialByCoefficients(FamilyObj(1),tmp[1],tmp[2]);
          x := x-r*FindPq(H,mu);
          if mu in x!.parts then n:= n+1; fi;
        fi;
      od;
      r := List(x!.coeffs, s->s <> 0 * v ^ 0);
      if false in r then
        x!.coeffs := ListBlist( x!.coeffs, r );
        x!.parts := ListBlist( x!.parts, r );
      fi;
    fi;

    ## having found x we add it to CDM
    CDM!.d[c]:=rec(coeffs:=x!.coeffs,
                  parts:=List(x!.parts, r->Position(CDM!.rows,r)) );

    ## for good measure we also add the Mullineux image of Pq(mu) to CDM
    ## (see LLT Theorem 7.2)
    n:=Position(CDM!.cols,ConjugatePartition(x!.parts[1]));
    if c<>n then             ## not self-image under MullineuxMap
      r:=List(x!.coeffs*v^0, i->Value(i,v^-1));     ## v -> v^-1
      for i in [Length(r),Length(r)-1..1] do   ## multiply by r[1]
        tmp:=ShallowCopy(CoefficientsOfLaurentPolynomial(r[i]));
        tmp[2]:=tmp[2]-CoefficientsOfLaurentPolynomial(r[1])[2];
        r[i]:=LaurentPolynomialByCoefficients(FamilyObj(1),tmp[1],tmp[2]);
      od;
      s:=List(x!.parts, mu->Position(CDM!.rows,ConjugatePartition(mu)));
      SortParallel(s,r);
      CDM!.d[n]:=rec(coeffs:=r,parts:=s);
    fi;
    return x;
  end
); # FindPq

## Strictly speaking the functions Sq() and Dq() are superfluous as
## these functions could be carried out by the decomposition matrix
## operations; however the point is that the crystal matrices are
## allowed to have missing columns, and if any needed columns missing
## these functions calculate them on the fly via using Pq().

## writes S(mu), from the Fock space, as a sum of D(nu)
InstallMethod(FindSq,"write Sq(mu) as sum of Dq(mu)",
  [IsAlgebraObj,IsList],
  function(H,mu) local x, r, CDM, n;
    n:=Sum(mu);

    ## we need the list of e-regular partitions. as we will have to
    ## create CrystalMatrices[n] anyway, we get the columns from there.
    if not IsBound(H!.CrystalMatrices[n]) or H!.CrystalMatrices[n]=fail then
      r:=Partitions(n);
      H!.CrystalMatrices[n]:=DecompositionMatrix(H,
        r, ERegularPartitions(H!.e,n), false);
    fi;
    CDM:=H!.CrystalMatrices[n];

    x:=Module(H,"Dq",0,[]);
    r:=Length(CDM!.cols);
    mu:=Position(CDM!.rows,mu);
    while r>0 and CDM!.cols[r]>=CDM!.rows[mu] do
      if not IsBound(CDM!.d[r]) then FindPq(H,CDM!.cols[r]); fi;
      n:=Position(CDM!.d[r].parts,mu);
      if n<>fail then
        x:=x+Module(H,"Dq",CDM!.d[r].coeffs[n],CDM!.cols[r]);
      fi;
      r:=r-1;
    od;
    return x;
  end
);

## write D(mu), from the Fock space, as a sum of S(nu)
InstallMethod(FindDq, "write Dq(mu) as sum of Sq(mu)",
  [IsAlgebraObj,IsList],
  function(H,mu) local inv, x, c, CDM;
    inv:=Module(H,"Sq",1,mu);
    x:=Module(H,"Dq",1,mu)-FindSq(H,mu);
    CDM:=H!.CrystalMatrices[Sum(mu)];
    c:=Position(CDM!.cols,mu);

    if IsBound(CDM!.inverse[c]) then
      return Module(H,"Sq",CDM!.inverse[c].coeffs,
                                List(CDM!.inverse[c].parts,c->CDM!.cols[c]));
    fi;

    while x<>0*x do
      c:=Length(x!.parts);
      inv:=inv+Module(H,"Sq",x!.coeffs[c],x!.parts[c]);
      x:=x-x!.coeffs[c]*FindSq(H,x!.parts[c]);
    od;

    ## now place answer back in CDM
    CDM!.inverse[Position(CDM!.cols,mu)]:=rec(coeffs:=inv!.coeffs,
             parts:=List(inv!.parts,c->Position(CDM!.cols,c)) );
    return inv;
  end
); # FindDq

