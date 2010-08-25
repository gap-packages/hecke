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

## 3.0: June 2010:
##   - Translated to GAP4

## Change log
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
##     with correspondining operation sets H.operations.X. In particular,
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
##                   'H' modules repectively. These functions exist
##                   only when 'H' is a Hecke algebra record.
#### Use NewModule(H,"S",...) instead of H.S(...)
##   W(), P(), F() : these are the corresponding functions for Weyl
##                   modules, PIMs, and simple modules of Schur algbras.
#### Use NewModule(H,"W",...) instead of H.W(...)
##   info          : this is a record with components
##                     version: SPECHT version number,
##                     Library: path to SPECHT library files
##                     SpechtDirectory: addition directory searched by
##                            SPECHT (defaults to current directory)
##                     Indeterminate: the indedeterminate used by SPECHT
##                            in the LLT algorithm when H.p=0.
#### TODO
##   operations    : apart from the obvious things like the Print()
##                   function for 'H' this record also contains the
##                   operation records for the modules S(), P() etc.
##                   as well as functions for manipulating these records
##                   and accessing decomposition matrix files. The most
##                   most important of these are:
##                     S, P, D, Pq, Sq, Dq : operations records for modules
#### TODO
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
#### Use New(H,...) instead of H.operations.New(...)
##                     Collect: like New() except that  coeffs and parts
##                       need not be sets (and may contain repeats).
#### Use Collect(H,...) instead of H.operations.Collect(...)
##                     NewDecompositionMatrix : creates a decomposition
##                       matrix.
#### TODO
##                     ReadDecompositionMatrix : reads, and returns, a 
##                       decomposition matrix file.
#### TODO
##                     KnownDecompositionMatrix : returns a decomposition
##                       matrix of a given size; will either extract this
##                       matrix from Specht's internal lists or call
##                       ReadDecompositionMatrix(), or try to calculate 
##                       the decomposition matrix (without using the
##                       crystalized decomposition matrices).
#### TODO
##                     FindDecompositionMatrix : like KnownDM except that
##                       it will calculate the crystalized DM if needed.
#### TODO
##   Ordering      : a function for ordering partitions; controls how
##                   decomposition matrices for H are printed.
#### Use SetOrdering(...) to control the output
##   e             : order of <q> in <R>
#### Use OrderOfQ(...) to extract the e from an algebra or a module 
#### corresponding to an algebra 
##   p             : characteristic of <R>
#### TODO
##   valuation     : the valuation map of [JM2]; used primarily by
##                   the q-Schaper theorem.D
#### TODO
##   HeckeRing     : bookkeeping string used primarily in searching for
##                   library files.
#### TODO
##   Pq(), Sq()    : Functions for computing elements of the Fock space
##                   when H.p=0 (used in LLT algorithm). Note that there is
##                   no Dq; also unlike their counter parts S(), P(), and
##                   D() they accept only partitions as arguments.
#### TODO
##
## 2. The module functions S(), P() and D() (and Schur equivalents)
## These functions return record 'x' which represents some 'H'--module.
## 'x' is a record with the following components:
##   H      : a pointer back to the corresponding algebra
##   module : one of "S", "P", "D", "Sq", "Pq", or "Dq", (not "W", or "F").
##   coeffs : a *set* of coefficients
##   parts  : the corresponding *set* of partitions
##   operations :
##       + - * / : for algebric manipulations
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
##   dimensions : a list of the dimenions of the simle modules; again
##          comuted only as needed.
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
##              to write an X-module as a linear comination of Y-modules.
##     Invert : calculates D(mu) using 'd'.
##     IsNewIndecomposable : the heart of the 'IsNewIndecomposable'
##              function.
##     Induce : for inducing decomposition matrices (non--crystallized).
##   P : a short-hand for d.H.P('d',<mu>).

InstallMethod(Specht,"generate a Hecke-Algebra object",
  [IsInt,IsInt,IsFunction,IsString],
  function(e,p,valuation,HeckeRing)
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
      info:=rec(version:=PackageInfo("specht")[1].Version, 
                Library:=Directory(
                  Concatenation(DirectoriesPackageLibrary("specht")[1]![1],"e",
                  String(e),"/")),
      ## We keep a copy of SpechtDirectory in H so that we have a
      ## chance of finding new decomposition matrices when it changes.
                SpechtDirectory:=Directory(".")),
                
      ## This record will hold any decomposition matrices which Specht()
      ## (or rather its derivatives) read in. This used to be a public 
      ## record; it is now private because q-Schur algebra matrices and
      ## Hecke algebra matrices might need to coexist. 
      DecompositionMatrices:=[],

      ## ordering used when printing decomposition matrices
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

    H.S:=function(arg) return NewModule(H,"S",arg); end;
    H.P:=function(arg) return NewModule(H,"P",arg); end;
    H.D:=function(arg) return NewModule(H,"D",arg); end;

    Objectify(HeckeType,H);

    return H;
  end
);

InstallMethod(Specht,"generate a Hecke-Algebra object",
  [IsInt,IsInt,IsFunction],
  function(e,p,val) local ring;
    if not IsPrime(p) 
    then Error("Specht(<e>,<p>,<val>), <p> must be a prime number"); 
    fi;
    if e=p then
      ring:=Concatenation("p",String(p),"sym");
      return Specht(e,p,val,ring);
    else
      return Specht(e,p,val,"unknown");
    fi; 
  end
);

InstallMethod(Specht,"generate a Hecke-Algebra object",
  [IsInt,IsInt],
  function(e,p) local val, ring;
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
    return Specht(e,p,val,ring);
  end
);

InstallMethod(Specht,"generate a Hecke-Algebra object",
  [IsInt],
  function(e) local val;
      if e=0 then val:=x->x;
      else
        val:=function(x)
          if x mod e = 0 then return 1;
          else return 0;
          fi;
        end;
      fi;
      return Specht(e,0,val,Concatenation("e",String(e), "p0")); 
  end
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

#F Returns true if S(mu)=D(mu) - note that this implies that mu is e-regular
## (if mu is not e-regular, fail is returned).     -- see [JM2]
## IsSimle(H,mu)
##   ** uses H.valuation #### FIXME check this?
InstallMethod(IsSimpleModuleOp,
  "test whether the given partition defines a simple module",
  [IsAlgebraObj,IsList],
  function(H,mu) local mud, simple, r, c, v;
    if not IsERegular(H!.e,mu) then return fail; 
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

#F Split an element up into compontents which have the same core.
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
      if cpos=false then 
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
  [IsAlgebraObjModule,IsHeckeSpecht], ## TODO Is this really only for specht modules?
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
## the Kleshcehev(-James) algorihm, or the supplied decomposition matrix.
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
    if x!.module{[1]}="S" then
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
    MullineuxMapOp(H!.e,mu);
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
        return false;
      else return ConjugatePartition(x!.parts[1]);
      fi;
  end
); #MullineuxMap

#F Calculates the Specht modules in sum_{i>0}S^lambda(i) using the
## q-analogue of Schaper's theorem. 
## Uses H.valuation. FIXME check this?
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
              if s<>false then
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
); #Omega

## MODULES #####################################################################
InstallMethod(Module,"create new module",[IsAlgebraObj,IsString,IsList,IsList],
  function(H,m,c,p)
    local module;
    ## TODO: Argument tests!?
    module := rec(H:=H,module:=m,coeffs:=c,parts:=p);

    if m = "S" and IsHecke(H) then Objectify(HeckeSpechtType,module);
    elif m = "P" and IsHecke(H) then Objectify(HeckePIMType,module);
    elif m = "D" and IsHecke(H) then Objectify(HeckeSimpleType,module);
    elif m = "Sq" and IsHecke(H) then Objectify(HeckeSpechtFockType,module);
    elif m = "Pq" and IsHecke(H) then Objectify(HeckePIMFockType,module);
    elif m = "Dq" and IsHecke(H) then Objectify(HeckeSimpleFockType,module);
    fi;

    return module;
  end
);

InstallMethod(Module,"create new module",[IsAlgebraObj,IsString,IsInt,IsList],
  function(H,m,c,p)
    return Module(H,m,[c],[p]);
  end
);

InstallMethod(Module,"create new module",[IsAlgebraObj,IsString,IsUnivariatePolynomial,IsList],
  function(H,m,c,p)
    return Module(H,m,[c],[p]);
  end
);

## Takes two lists, one containing coefficients and the other the
## corresponding partitions, and orders them lexicogrphcailly collecting
## like terms on the way. We use a variation on quicksort which is
## induced by the lexicographic order (if parts contains partitions of
## different integers this can lead to an error - which we don't catch).
InstallMethod(Collect,"TODO",[IsAlgebraObj,IsString,IsList,IsList],
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
InstallMethod(MakeSpecht,"S()->S()",[IsHeckeSpecht,IsBool],
  function(x,silent) return x; end
);

## Here I only allow for linear combinations of projectives which 
## have non-negative coefficients; the reason for this is that I
## can't see how to do it otherwise. The problem is that in the
## Grothendieck ring there are many ways to write a given linear
## combination of Specht modules (or PIMs).
InstallMethod(MakePIM,"S()->P()",[IsHeckeSpecht,IsBool],
  function(x,silent) local proj;
    if x=fail or x=0*x then return x; 
    elif x!.parts=[[]] then return Module(x!.H,"P",x!.coeffs[1],[]);
    fi;
  
    proj:=Module(x!.H,"P",0,[]);
    while x<>false and x<>0*x and 
    ( not IsHeckeSpecht(x!.H) or IsERegular(x!.H!.e,x!.parts[Length(x!.parts)]) ) do
      proj:=proj+Module(x!.H,"P",x!.coeffs[Length(x!.parts)], 
                                      x!.parts[Length(x!.parts)]);
      x:=x+MakeSpecht(
                Module(x!.H,"P",-x!.coeffs[Length(x!.parts)],
                                      x!.parts[Length(x!.parts)]),true);
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

InstallMethod(MakeSimple,"S()->D()",[IsHeckeSpecht,IsBool],
  function(x,silent) return x; end ## FIXME
);

## The P->S functions are quite involved.
   
#F Writes x, which is a sum of indecomposables, as a sum of S(nu)'s if 
## possible. We first check to see if the decomposition matrix for x is
## stored somewhere, and if not we try to calculate what we need. If we
## can't do this we return false.
InstallMethod(MakeSpecht,"P()->S()",[IsHeckePIM,IsBool],
  function(x,silent) return x; end ## FIXME
);

InstallMethod(MakePIM,"P()->P()",[IsHeckePIM,IsBool],
  function(x,silent) return x; end
);

InstallMethod(MakeSimple,"P()->D()",[IsHeckePIM,IsBool],
    function(x,silent)
      x:=MakeSpecht(x,silent);
      if x=fail then return x;
      else return MakeSimple(x,silent); 
      fi;
    end
);

#F Writes D(mu) as a sum of S(nu)'s if possible. We first check to see
## if the decomposition matrix for Sum(mu) is stored in the library, and
## then try to calculate it directly. If we are unable to do this either
## we return fail.
InstallMethod(MakeSpecht,"D()->S()",[IsHeckeSimple,IsBool],
  function(x,silent) return x; end ## FIXME
);

InstallMethod(MakePIM,"D()->P()",[IsHeckeSimple,IsBool],
  function(x,silent) 
      x:=MakeSpecht(x,silent);
      if x=fail then return x;
      else return MakePIM(x,silent);
      fi;
    end
);

InstallMethod(MakeSimple,"D()->D()",[IsHeckeSimple,IsBool],
  function(x,silent) return x; end
);

## ARITHMETICS #################################################################
InstallMethod(\=,"compare modules",[IsAlgebraObjModule,IsAlgebraObjModule],
  function(a,b) return a!.H=b!.H and a!.module=b!.module 
    and Set(Zip(a!.coeffs,a!.parts))=Set(Zip(b!.coeffs,b!.parts)); end
);

InstallMethod(\+,"add modules",[IsAlgebraObjModule,IsAlgebraObjModule],
  function(a,b)
    local i, j, ab, x;
    
    if a=fail or b=fail then return false;
    elif a=0*a then return b; 
    elif b=0*b then return a; 
    elif a!.H<>b!.H then 
      Error("modules belong to different Grothendieck rings");
    fi;

    if a!.module<>b!.module then # only convert to Specht modules if different
      if Length(a!.module) <> Length(b!.module) then
        Error("AddModule(<a>,<b>): can only add modules of same type.");
      fi;
      a:=MakeSpecht(a,false);
      b:=MakeSpecht(b,false);
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
    if a=fail or b=fail then return false;
    elif a!.H<>b!.H then 
      Error("modules belong to different Grothendieck rings");
    fi;
    # a:=MakeSpecht(a,false); # TODO Reconsider?
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
    x:=MakeSpecht(a,false) * MakeSpecht(b,false);
    nx:=MakePIM(x,true);
    if nx<>fail then return nx; else return x; fi;
  end
);

InstallMethod(\*,"multiply simple modules",[IsHeckeSimple,IsHeckeSimple],
  function(a,b) local x, nx;
    x:=MakeSpecht(a,false) * MakeSpecht(b,false);
    nx:=MakeSimple(x,true);
    if nx<>fail then return nx; else return x; fi;
  end
); 

InstallMethod(\*,"multiply modules",[IsAlgebraObjModule,IsAlgebraObjModule],
  function(a,b) 
    return MakeSpecht(a,false) * MakeSpecht(b,false);
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
      a:=MakeSpecht(a,true);
      b:=MakeSpecht(b,true);
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
InstallMethod(Coefficient, "extract coefficient of a partition from module",
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
);  # PositiveCoefficients

#F Returns true if all coefficients are integral
InstallMethod(IntegralCoefficients, "test if all coefficients are integral",
  [IsAlgebraObjModule],
  function(x) local c;
    return ForAll(x!.coeffs, c->IsInt(c));
  end
); # IntegralCoefficients

## INDUCTION AND RESTRICTION ###################################################

## The next functions are for restricting and inducing Specht
## modules. They all assume that their arguments are indeed Specht
## modules; conversations are done in H.operations.X.Y() as necessary.

## r-induction: on Specht modules:
InstallMethod(RInducedModule, "r-induction on specht modules",
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
## ordinary indcution or q-induction (in the Fock space). We could
## write H.operations.X.SInduce to as to make this choice for us, or
## do q-induction always, setting v=1 afterwards, but this seems the
## better choice.
InstallMethod(SInducedModule,"string-induction on specht modules",
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
InstallMethod(RRestrictedModule,"r-restriction on specht modules",
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
InstallMethod(SRestrictedModule,"string-restriction on specht modules",
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
InstallMethod(RInducedModule, "r-induction for specht modules",
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

InstallMethod(RInducedModule, "r-induction for projective indecomposable modules",
  [IsAlgebraObj, IsHeckePIM, IsList],
  function(H, x, list) local nx;
    x:=RInducedModule(H,MakeSpecht(x,false),list);
    if x=fail or x=0*x then return x; fi;
    nx:=MakePIM(x,false);
    if nx<>fail then return nx; else return x; fi;
  end
); 

InstallMethod(RInducedModule, "r-induction for simple modules",
  [IsAlgebraObj, IsHeckeSimple, IsList],
  function(H, x, list) local nx;
    x:=RInducedModule(H,MakeSpecht(x,false),list);
    if x=fail or x=0*x then return x; fi;
    nx:=MakeSimple(x,false);
    if nx<>fail then return nx; else return x; fi;
  end
); # RInducedModule

InstallMethod(RRestrictedModule, "r-restriction for specht modules",
  [IsAlgebraObj, IsHeckeSpecht, IsList],
  function(H, x, list) local r;
    if x=fail or x=0*x then return x;
    elif list=[] then return RRestrictedModule(x,1,0);
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

InstallMethod(RRestrictedModule, "r-restriction for projective indecomposable modules",
  [IsAlgebraObj, IsHeckePIM, IsList],
  function(H, x, list) local nx;
    x:=RRestrictedModule(H,MakeSpecht(x,false),list);
    if x=fail or x=0*x then return x; fi;
    nx:=MakePIM(x,false);
    if nx<>fail then return nx; else return x; fi;
  end
);   

InstallMethod(RRestrictedModule, "r-restriction for simple modules",
  [IsAlgebraObj, IsHeckeSimple, IsList],
  function(H, x, list) local nx;
    x:=RRestrictedModule(H,MakeSpecht(x,false),list);
    if x=fail or x=0*x then return x; fi;
    nx:=MakeSimple(x,false);
    if nx<>fail then return nx; else return x; fi;
  end
); # RRestrictedModule

InstallMethod(SInducedModule,"string induction for specht modules",
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
      Error("SInduce, r-induction is not defined when e=0.");
    elif list[2]>H!.e or list[2]<0 then
      Error("SInduce, r-induction is defined only when 0<=r<e.\n");
    else return SInducedModule(H, x, H!.e, list[1], list[2]);
    fi;
  end
);  

InstallMethod(SInducedModule, "string induction for projective indecomposable modules",
  [IsAlgebraObj, IsHeckePIM, IsList],
  function(H, x, list) local nx;
    x:=SInducedModule(H,MakeSpecht(x,false),list);
    if x=fail or x=0*x then return x; fi;
    nx:=MakePIM(x,false);
    if nx<>fail then return nx; else return x; fi;
  end
);   

InstallMethod(SInducedModule, "string induction for simple modules",
  [IsAlgebraObj, IsHeckeSimple, IsList],
  function(H, x, list) local nx;
    x:=SInducedModule(H,MakeSpecht(x,false),list);
    if x=fail or x=0*x then return x; fi;
    nx:=MakeSimple(x,false);
    if nx<>fail then return nx; else return x; fi;
  end
); # SInducedModule

InstallMethod(SRestrictedModule,"string restriction for specht modules",
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

InstallMethod(SRestrictedModule, "string restriction for projective indecomposable modules",
  [IsAlgebraObj, IsHeckePIM, IsList],
  function(H, x, list) local nx;
    x:=SRestrictedModule(H,MakeSpecht(x,false),list);
    if x=fail or x=0*x then return x; fi;
    nx:=MakePIM(x,false);
    if nx<>fail then return nx; else return x; fi;
  end
);   

InstallMethod(SRestrictedModule, "string restriction for simple modules",
  [IsAlgebraObj, IsHeckeSimple, IsList],
  function(H, x, list) local nx;
    x:=SRestrictedModule(H,MakeSpecht(x,false),list);
    if x=fail or x=0*x then return x; fi;
    nx:=MakeSimple(x,false);
    if nx<>fail then return nx; else return x; fi;
  end
); #SRestrictedModule

## DECOMPOSITION MATRICES ######################################################

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
#### FIXME Necessary?
####             P:=function(d,mu)     ## a lazy helper
####               return d.operations.P.S(d, d.H.operations.New("P",1,mu));
####             end
    );
    
    if decompmat then 
      return Objectify(DecompositionMatrixType,d);
    else
      return Objectify(CrystalDecompositionMatrixType,d);
    fi;
  end
);   # DecompositonMatrix

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
    return ReadDecompositionMatrix(H,file,crystal);
  end
);

InstallMethod(ReadDecompositionMatrix, "load matrix from library",
  [IsAlgebraObj,IsString,IsBool],
  function(H, n, crystal)
    local msg, file, M, d, c, parts, coeffs, p, x, r, cm, rm,
      SpechtDirectory;
    
    ## The following directory is searched by ReadDecompositionMatrix()
    ## when it is looking for decomposition matrices. By default, it points
    ## to the current directory (if set, the current directory is not
    ## searched).
    if not IsBound(SpechtDirectory) then SpechtDirectory:=Directory("."); fi;

    A_Specht_Decomposition_Matrix:=fail; ## just in case
    
    file:=n;
    if crystal then
      msg:="ReadCrystalMatrix-";
    else 
      msg:="ReadDecompositionMatrix-";
    fi;

    d:=fail;
  
    Read(Filename([Directory("."),SpechtDirectory,H!.info.Library],file));
  
    if A_Specht_Decomposition_Matrix<>fail then   ## extract matrix from M
      M:=A_Specht_Decomposition_Matrix;
      A_Specht_Decomposition_Matrix:=fail;
      r:=Set(M.rows); c:=Set(M.cols);
      if IsHecke(H) and r=c then
        d:=DecompositionMatrix(H,r,
              Filtered(c,x->IsERegular(H!.e,x)),not IsBound(M.crystal));
      elif IsHecke(H) then
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
                    FamilyObj(One(H!.info.Indeterminate)),
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
              if c<>false then
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
    if crystal then H!.CrystalMatrices[n]:=d; fi;
    return d;
  end
); # ReadDecompositionMatrix

