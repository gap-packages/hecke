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
##  - Changed ReadDecompositionmatrix so that it automatically reads the
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

## 1. Specht()
## Specht() is the main function in specht.g, it returns a record 'H'
## which represents the family of Hecke algebras, or Schur algebras, 
## corresponding to some fixed field <R> and parameter <q>. 'H' has 
## the components:
##   IsSpecht      : is either 'true' or 'false' depending on whether
##                   'H' is a Hecke algebra or Schur algebra resp.
##   S(), P(), D() : these three functions return records which 
##                   represent Specht modules, PIMs, and simple
##                   'H' modules repectively. These functions exist
##                   only when 'H' is a Hecke algebra record.
##   W(), P(), F() : these are the corresponding functions for Weyl
##                   modules, PIMs, and simple modules of Schur algbras.
##   info          : this is a record with components
##                     version: SPECHT version number,
##                     Library: path to SPECHT library files
##                     SpechtDirectory: addition directory searched by
##                            SPECHT (defaults to current directory)
##                     Indeterminate: the indedeterminate used by SPECHT
##                            in the LLT algorithm when H.p=0.
##   operations    : apart from the obvious things like the Print()
##                   function for 'H' this record also contains the
##                   operation records for the modules S(), P() etc.
##                   as well as functions for manipulating these records
##                   and accessing decomposition matrix files. The most
##                   most important of these are:
##                     S, P, D, Pq, Sq, Dq : operations records for modules
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
##                     Collect: like New() except that  coeffs and parts
##                       need not be sets (and may contain repeats).
##                     NewDecompositionMatrix : creates a decomposition
##                       matrix.
##                     ReadDecompositionMatrix : reads, and returns, a 
##                       decomposition matrix file.
##                     KnownDecompositionMatrix : returns a decomposition
##                       matrix of a given size; will either extract this
##                       matrix from Specht's internal lists or call
##                       ReadDecompositionMatrix(), or try to calculate 
##                       the decomposition matrix (without using the
##                       crystalized decomposition matrices).
##                     FindDecompositionMatrix : like KnownDM except that
##                       it will calculate the crystalized DM if needed.
##   Ordering      : a function for ordering partitions; controls how
##                   decomposition matrices for H are printed.
##   e             : order of <q> in <R>
##   p             : characteristic of <R>
##   valuation     : the valuation map of [JM2]; used primarily by
##                   the q-Schaper theorem.D
##   HeckeRing     : bookkeeping string used primarily in searching for
##                   library files.
##   Pq(), Sq()    : Functions for computing elements of the Fock space
##                   when H.p=0 (used in LLT algorithm). Note that there is
##                   no Dq; also unlike their counter parts S(), P(), and
##                   D() they accept only partitions as arguments.
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

InstallMethod(
  Specht,
  "generate a Hecke-Algebra object",
  [IsInt],
  function(e)
    local H;
    
    H := rec(
      e:=e,
      ## bits and pieces about H
      info:=rec(version:=PackageInfo("specht")[1].Version),

      ## ordering used when printing decomposition matrices
      Ordering:=Lexicographic,
    );

    H.S:=function(arg) return NewModule(H,"S",arg); end;
    H.P:=function(arg) return NewModule(H,"P",arg); end;
    H.D:=function(arg) return NewModule(H,"D",arg); end;

    Objectify(HeckeType,H);

    return H;
  end
);

InstallMethod(OrderOfQ,"reading access to H.e",[IsAlgebraObj],
  function(H) return H!.e; end
);

InstallMethod(OrderOfQ,"reading access to x.H.e",[IsAlgebraObjModule],
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
InstallMethod(MakeSpecht,"S()->P()",[IsHeckeSpecht,IsBool],
  function(x,silent) return x; end ## FIXME
);

InstallMethod(MakeSpecht,"S()->D()",[IsHeckeSpecht,IsBool],
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

InstallMethod(MakeSpecht,"P()->P()",[IsHeckePIM,IsBool],
  function(x,silent) return x; end
);

InstallMethod(MakeSpecht,"P()->D()",[IsHeckePIM,IsBool],
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
## we return false.
InstallMethod(MakeSpecht,"D()->S()",[IsHeckeSimple,IsBool],
  function(x,silent) return x; end ## FIXME
);

InstallMethod(MakeSpecht,"D()->P()",[IsHeckeSimple,IsBool],
  function(x,silent) 
      x:=MakeSpecht(x,silent);
      if x=fail then return x;
      else return MakePIM(x,silent);
      fi;
    end
);

InstallMethod(MakeSpecht,"D()->D()",[IsHeckeSimple,IsBool],
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

InstallMethod(\*,"multiply module by scalar",[IsScalar,IsHeckeSpecht],
  function(n,b)
    if n = 0 
    then return Module(b!.H, b!.module, 0, []);
    else return Module(b!.H, b!.module, n*b!.coeffs, b!.parts);
    fi;
  end
);

InstallMethod(\*,"multiply module by scalar",[IsHeckeSpecht,IsScalar],
  function(a,n)
    if n = 0 
    then return Module(a!.H, a!.module, 0, []);
    else return Module(a!.H, a!.module, n*a!.coeffs, a!.parts);
    fi;
  end
);

InstallMethod(\*,"multiply modules",[IsHeckeSpecht,IsHeckeSpecht],
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
    else return Collect(b!.H, b!.module, ab[1], ab[2]); #TODO
    fi;
  end  # MultiplyModules
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

