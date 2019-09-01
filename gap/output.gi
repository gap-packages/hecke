#######################################################################
##  Hecke - output.gi : printing functions                           ##
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

InstallMethod(PrintObj, "simple algebra output", [IsAlgebraObj],
	function(x) Print(AlgebraString(x)); end
);

InstallMethod(ViewString, "compact algebra output", [IsHecke],
	function(H)
		return Concatenation("<Hecke algebra with e = ",String(H!.e),">");
	end
);

InstallMethod(ViewString, "compact algebra output", [IsSchur],
	function(S)
		return Concatenation("<Schur algebra with e = ",String(S!.e),">");
	end
);

InstallMethod(ViewObj, "compact algebra output", [IsAlgebraObj],
	function(H) Print(ViewString(H)); end
);

InstallMethod(DisplayString, "pretty algebra output", [IsAlgebraObj],
  function(x) return AlgebraString(x); end
);

InstallMethod(Display, "pretty algebra output", [IsAlgebraObj],
  function(x) Print(DisplayString(x),"\n"); end
);

InstallMethod(AlgebraString, "generic algebra output", [IsHecke],
	function(H) local p, pq, ring;
		if H!.p<>0
		then p:=Concatenation("p=",String(H!.p),", ");
		else p:="";
		fi;
		if IsBound(H!.pq)
		then pq:=", Pq()";
		else pq:="";
		fi;
		if H!.p<>H!.e and H!.p<>0
		then ring:=Concatenation(", HeckeRing=\"", String(H!.HeckeRing), "\")");
		else ring:=")";
		fi;

		return
			Concatenation("Specht(e=",String(H!.e),", ",p,"S(), P(), D()",pq,ring);
	end
);

InstallMethod(AlgebraString, "generic algebra output", [IsSchur],
	function(S) local p, pq, ring;
		if S!.p<>0
		then p:=Concatenation("p=",String(S!.p),", ");
		else p:="";
		fi;
		if IsBound(S!.pq)
		then pq:=", Pq()";
		else pq:="";
		fi;
		if S!.p<>S!.e and S!.p<>0
		then ring:=Concatenation(", HeckeRing=\"", String(S!.HeckeRing), "\")");
		else ring:=")";
		fi;

		return
			Concatenation("Schur(e=",String(S!.e),", ",p,"W(), P(), F()",pq,ring);
	end
);

InstallMethod(PrintObj, "simple module output", [IsAlgebraObjModule],
  function(m) Print(ModuleString(m,false)); end
);

InstallMethod(ViewString, "compact module output", [IsAlgebraObjModule],
  function(m)
    return Concatenation("<direct sum of ",
      String(Length(m!.parts))," ",m!.module,"-modules>");
  end
);

InstallMethod(ViewObj, "compact module output",[IsAlgebraObjModule],
  function(m) Print(ViewString(m)); end
);

InstallMethod(DisplayString, "pretty module output", [IsAlgebraObjModule],
  function(m) return ModuleString(m,true); end
);

InstallMethod(Display, "pretty module output", [IsAlgebraObjModule],
  function(m) Print(DisplayString(m),"\n"); end
);

InstallMethod(ModuleString, "generic module output", [IsAlgebraObjModule,IsBool],
  function(a,pp)
    local x, len, n, star, tmp, coefficients, valuation, str;

    str := "";
    if Length(a!.parts) > 1 then
      n:=Sum(a!.parts[1]);
      if not ForAll(a!.parts, x->Sum(x)=n) then
        Error(a!.module,"(x), <x> contains mixed partitions\n\n");
      fi;
    fi;
    if pp then star:=""; ## pretty printing
    else star:="*";
    fi;

    if Length(a!.parts)=0 or a!.coeffs=[0] then
      a!.coeffs:=[ 0 ]; a!.parts:=[ [] ];
      str:=StringFold(str,["0",star,a!.module,"()"]);
    else
      for x in [Length(a!.parts),Length(a!.parts)-1..1] do
        if IsLaurentPolynomial(a!.coeffs[x]) then
          tmp := CoefficientsOfLaurentPolynomial(a!.coeffs[x]);
          coefficients := tmp[1];
          valuation := tmp[2];
          if Length(coefficients)=1 then
            if valuation=0 then
              if coefficients=[-1] then Append(str," - ");
              elif coefficients[1]<0 then
                str:=StringFold(str,[coefficients[1],star]);
              else
                if x<Length(a!.parts) then Append(str," + "); fi;
                if coefficients<>[1] then
                  str:=StringFold(str,[coefficients[1],star]);
                fi;
              fi;
            else
              if x<Length(a!.parts) and coefficients[1]>0 then
                Append(str," + ");
              fi;
              str:=StringFold(str,[a!.coeffs[x],star]);
            fi;
          elif coefficients[Length(coefficients)]<0
          then str:=StringFold(str,[" - (",-a!.coeffs[x],")", star]);
          else
            if x<Length(a!.parts) then Append(str, " + "); fi;
            str:=StringFold(str,["(",a!.coeffs[x],")", star]);
          fi;
        else
          if a!.coeffs[x]=-1 then Append(str, " - ");
          elif a!.coeffs[x]<0 then str:=StringFold(str,[a!.coeffs[x],star]);
          else
            if x<Length(a!.parts) then Append(str, " + "); fi;
            if a!.coeffs[x]<>1
              then str:=StringFold(str,[a!.coeffs[x],star]); fi;
          fi;
        fi;
        if pp then
          tmp := LabelPartition(a!.parts[x]);
        else
          tmp := StringPartition(a!.parts[x]);
        fi;
        str:=StringFold(str,[a!.module,"(",tmp,")"]);
      od;
    fi;
    return str;
  end
);

InstallMethod(PrintObj, "simple decomposition matrix output", [IsDecompositionMatrix],
	function(d) Print("DecompositionMatrix(",d!.H,
	  ",",d!.rows,",",d!.cols,",",true,")"); end
);

InstallMethod(PrintObj, "simple decomposition matrix output", [IsCrystalDecompositionMatrix],
	function(d) Print("CrystalDecompositionMatrix(",d!.H,
	  ",",d!.rows,",",d!.cols,",",false,")"); end
);

InstallMethod(ViewString, "compact decomposition matrix output", [IsDecompositionMatrix],
	function(d)
	    return Concatenation("<",String(Length(d!.rows)),"x",String(Length(d!.cols)),
	      " decomposition matrix>");
	end
);

InstallMethod(ViewObj, "compact decomposition matrix output", [IsDecompositionMatrix],
	function(d) Print(ViewString(d)); end
);

InstallMethod(DisplayString, "pretty decomposition matrix output", [IsDecompositionMatrix],
  function(d) return DecompositionMatrixString(d); end
);

InstallMethod(Display, "pretty decomposition matrix output", [IsDecompositionMatrix],
  function(d) Print(DisplayString(d),"\n"); end
);

## Pretty printing of a (decomposition) matrix d.
## Actually, d can be any record having .d=matrix, .labels=[strings],
## row, and cols components (this is also used by KappaMatrix for example).
InstallMethod(DecompositionMatrixString,"generic decomposition matrix output",
  [IsDecompositionMatrix],
  function(d)
    local rows, cols, r, c, col, len, M, rowlabel,
          spacestr, dotstr, i, str;

    str:="";

    ## if have to fix up the ordering before printing d
    rows:=StructuralCopy(d!.rows);
    cols:=StructuralCopy(d!.cols);
    if d!.H!.Ordering=Lexicographic then
      rows:=rows{[Length(rows),Length(rows)-1..1]};
      cols:=cols{[Length(cols),Length(cols)-1..1]};
    else
      Sort(rows, d!.H!.Ordering);
      Sort(cols, d!.H!.Ordering);
    fi;
    rows:=List(rows, r->Position(d!.rows,r));
    cols:=List(cols, c->Position(d!.cols,c));

    rowlabel:=List(d!.rows, LabelPartition);

    M:=-Maximum( List(rows, r->Length(rowlabel[r])) );

    ## Find out how wide the columns have to be (very expensive for
    ## crystallized matrices - also slightly incorrect as String(<poly>)
    ## returns such wonders as (2)*v rather than 2*v).
    len:=1;
    for i in d!.d do
      if i.coeffs<>[] then
        len:=Maximum(len,Maximum(List(i.coeffs,
                               j->Length(String(j)))));
      fi;
    od;
    spacestr:=String("",len);
    dotstr:=String(".",len);
    col:=0;
    for r in rows do
      Append(str,Concatenation(String(rowlabel[r],M),"| "));
      if d!.rows[r] in d!.cols then col:=col+1; fi;
      for c in [1..Length(cols)] do
        if IsBound(d!.d[cols[c]]) and r in d!.d[cols[c]].parts then
            Append(str,String(d!.d[cols[c]].coeffs[Position(d!.d[cols[c]].parts,r)], len));
          if c<>cols[1] then Append(str," "); fi;
        elif c<=col then Append(str, Concatenation(dotstr, " "));
        else Append(str, Concatenation(spacestr, " "));
        fi;
      od;
      if rows[r]<>Length(d!.rows) then Append(str,"\n"); fi;
    od;
    return str;
  end
); # DecompositionMatrixString

## adding this string if it does not already exist.
InstallMethod(LabelPartition, "pretty partition output", [IsList],
  function(mu) local n, m, label, p;
    n:=Sum(mu);
    if n<2 then return String(n); fi;
    label:="";
    for p in Collected(mu) do
      if p[2]=1 then label:=Concatenation(String(p[1]),",",label);
      else label:=Concatenation(String(p[1]),"^",String(p[2]),",",label);
      fi;
    od;
    return label{[1..Length(label)-1]};
  end
);

#F Returns a string for ModuleString() from SpechtParts.labels, adding this
## string if it does not already exist.
InstallMethod(StringPartition, "ergonomic partition output", [IsList],
  function(mu) local m, string, p;
    if mu=[] or mu=[0] then return "0";
    else return TightStringList(mu);
    fi;
  end
);

#F Returns a string of the form "l1,l2,...,lk" for the list [l1,l2,..,lk]"
## which contains no spaces, and has first element 1.
InstallMethod(TightStringList, "ergonomic list output", [IsList],
  function(list) local s, l;
    if list=[] then return ""; fi;
    s:=String(list[1]);
    for l in [2..Length(list)] do
      s:=Concatenation(s,",",String(list[l]));
    od;
    return s;
  end
);

## Keep ourselves honest when inducing decomposition matrices
InstallMethod(BUGOp,"hopefully, no-one will ever see this function ;-)",
  [IsString,IsInt],
  function(msg,pos) BUGOp(msg,pos,[]); end
);

InstallMethod(BUGOp,"hopefully, no-one will ever see this function ;-)",
  [IsString,IsInt,IsList],
  function(msg,pos,list) local a;
     PrintTo("*errout*",
             "\n\n *** You have uncovered a bug in Hecke's function ",
              msg,"() - #", pos, "\n *** Please e-mail all possible ",
              "details to ", PackageInfo("hecke")[1].Persons[1].Email,"\n");
     for a in list do
       PrintTo("*errout*", a);
     od;
     Error("\n");
  end
);

