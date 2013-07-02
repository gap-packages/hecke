#######################################################################
##                                                                   ##
##   Hecke - tableaux.gi : Young tableaux                            ##
##                                                                   ##
##     This file contains functions for generating all kinds of      ##
##     young tableaux                                                ##
##                                                                   ##
##     Dmitriy Traytel                                               ##
##     (heavily using the GAP3-package SPECHT 2.4 by Andrew Mathas)  ##
##                                                                   ##
#######################################################################

## Hecke 1.4: June 2013
##   - fixed a bug (reported by Rudolf Tango) in
##     SemiStandardTableau(): ignore trailing zeros
##     in the type of a tableau.
##   - Tableau constructor can now actually create
##     the empty tableau (again spotted by Rudolf Tango)

## Hecke 1.0: June 2010:
##   - Translated to GAP4

## SPECHT Change log
## July 1997
##   o fixed a bug (reported by Schmuel Zelikson) in
##     SemiStandardTableau(); added Type and Shape
##     functions for tableaux and allowed the type of
##     a tableau to be a composition which has parts
##     which are zero.

## March 1996

#F constructor of Young tableaux
InstallMethod(TableauOp, "tableau constructor",[IsList],
  function(tab)
    local isSemiStandardTab, isStandardTab;

    isSemiStandardTab := function(tab)
      local d, r, c;

      if tab=[] then
        return true;
      else
        d:=List([1..Length(tab[1])], r->[]);
        for r in [1..Length(tab)] do
	  for c in [1..Length(tab[r])] do
	    d[c][r]:=tab[r][c];
	  od;
        od;
        return ForAll(tab, r->r=SortedList(r))
	  and ForAll(d, r->IsSet(r));
      fi;
    end;

    ## assuming semistandardness already been checked
    isStandardTab := function(tab)
      return ForAll(tab, r->IsSet(r));
    end;

    ## validity checks ##
    if not tab=[] and ForAll(tab, x -> IsList(x) and
        not ForAll(x, y-> ##not IsBound(y) or ## TODO: desirable
            IsInt(y))) then
      Error("argument must be a list of lists of integers.");
    elif IsSortedList(Reversed(List(tab,Length)))=false then
      Error("row lengths must decrease ...");
    fi;
    ## ############### ##

    if isSemiStandardTab(tab) then
      if isStandardTab(tab) then
        return Objectify(StandardTableauType,[tab]);
      else
        return Objectify(SemiStandardTableauType,[tab]);
      fi;
    else
      return Objectify(TableauType,[tab]);
    fi;
  end
);

## OUTPUT ######################################################################
InstallMethod(PrintObj,"for tableaux",[IsTableau],
  function(tab) Print(Concatenation("Tableau( ",String(tab![1])," )")); end
);

InstallMethod(
  Specht_PrettyPrintTableau,
  [IsTableau],
  function(tab)
    local t, r, c, str;

    t := tab![1];
    str := "Tableau:\n";
    for r in [1..Length(t)] do
      for c in [1..Length(t[r])] do
        str := Concatenation(str,String(t[r][c]),"\t");
      od;
      str := Concatenation(str,"\n");
    od;
    return str;
  end
);

InstallMethod(
  DisplayString,
  "for tableaux",
  [IsStandardTableau],
  function(tab)
    return Concatenation("Standard ",Specht_PrettyPrintTableau(tab));
  end
);

InstallMethod(
  DisplayString,
  "for tableaux",
  [IsSemiStandardTableau],
  function(tab)
    return Concatenation("Semi-Standard ",Specht_PrettyPrintTableau(tab));
  end
);

InstallMethod(DisplayString,"for tableaux",[IsTableau],
  Specht_PrettyPrintTableau
);

InstallMethod(Display,"for tableaux",[IsTableau],
  function(tab) Print(DisplayString(tab)); end
);

InstallMethod(ViewString,"for tableaux",[IsTableau],
  function(tab)
    return Concatenation("<tableau of shape ",String(ShapeTableau(tab)),">");
  end
);

InstallMethod(ViewObj,"for tableaux",[IsTableau],
  function(tab) Print(ViewString(tab)); end
);
## #############################################################################

#F the dual tableaux of t = [ [t_1, t_2, ...], [t_k, ...], ... ]
InstallMethod(
  ConjugateTableau,
  [IsTableau],
  function(tab)
    local  t, d, r, c;

    t:=tab![1];
    d:=List([1..Length(t[1])], r->[]);
    for r in [1..Length(t)] do
      for c in [1..Length(t[r])] do
        d[c][r]:=t[r][c];
      od;
    od;
    return Tableau(d);
  end   ## ConjugateTableau
);

#F Returns the list of semi-standard nu-tableaux of content mu.
##   Usage: SemiStandardTableaux(nu, mu) or SemiStandardTableaux(nu).
## If mu is omitted the list of all semistandard nu-tableaux is return;
## otherwise only those semistandard nu-tableaux of type mu is returned.
##   Nodes are placed recursively via FillTableau in such a way so as
## so avoid repeats.
InstallMethod(SemiStandardTableauxOp, [IsList,IsList],
  function(nu,mu)
    local FillTableau, ss, i, maxnz;

    maxnz:=Length(mu);
    while maxnz>0 and mu[maxnz]=0 do
      maxnz:=maxnz - 1;
    od;
    mu:=mu{[1..maxnz]};

    #F FillTableau adds upto <n>nodes of weight <i> into the tableau <t> on
    ## or below its <row>-th row (similar to the Littlewood-Richardson rule).
    FillTableau:=function(t, i, n, row)
      local max, nn, nodes, nt, r;

      if row>Length(mu) then return;
      elif n=0 then          # nodes from i-th row of mu have all been placed
        if i=Length(mu) then # t is completely full
          Add(ss, t);
          return;
        else
          while i<Length(mu) and n=0 do
            i:=i + 1;          # start next row of mu
            n:=mu[i];          # mu[i] nodes to go into FillTableau
          od;
          row:=1;            # starting from the first row
        fi;
      fi;
      for r in [row..Length(t)] do
        if Length(t[r]) < nu[r] then
          if r = 1 then max:=nu[1]-Length(t[1]);
          elif Length(t[r-1]) > Length(t[r]) and t[r-1][Length(t[r]+1)]<i
               and Length(t[r-1])>=n then
            max:=Position(t[r-1], i);
            if max=fail then max:=Length(t[r-1])-Length(t[r]); #no i in t[r-1]
            else max:=max-1-Length(t[r]);
            fi;
          else max:=0;
          fi;
          max:=Minimum(max, n, nu[r]-Length(t[r]));
          nodes:=[];
          for nn in [1..max] do
            Add(nodes, i);
            nt:=StructuralCopy(t);
            Append(nt[r], nodes);
            FillTableau(nt, i, n - nn, r + 1);
          od;
        fi;
      od;
      r:=Length(t);
      if r < Length(nu)  and n <= nu[r+1] and n <= Length(t[r])
      and t[r][n] < i  then
        Add(t, List([1..n], nn->i));
        if i < Length(mu) then FillTableau(t, i+1, mu[i+1], 1);
        else Add(ss, t);
        fi;
      fi;
    end;

    if Sum(nu) <> Sum(mu) then
      Error("<nu> and <mu> must be partitions of the same integer.\n");
    fi;

    ## no semi-standard nu-tableau with content mu
    if Dominates(mu, nu) then
      if mu<>nu then return [];
      else return [Tableau(List([1..Length(mu)], i->List([1..mu[i]], ss->i)))];
      fi;
    fi;

    ss:=[]; ## will hold the semi-standard tableaux
    FillTableau([ List([1..mu[1]],i->1) ], 2, mu[2], 1);

    return List(Set(ss),Tableau);
  end   ## SemiStandardTableaux
);

InstallMethod(SemiStandardTableauxOp, [IsList],
  function(nu)
    local ss, i, mu;

    ss:=[];
    for mu in Partitions(Sum(nu)) do
      if Dominates(mu, nu) then
        if mu = nu then
          Append(ss,
              [ Tableau(List([1..Length(mu)], i->List([1..mu[i]], ss->i)) )]);
          return ss;
        fi;
      else
          Append(ss,SemiStandardTableaux(nu,mu));
      fi;
    od;
    return ss;
  end
);

#F Standard tableau of shape ([nu,] mu)
InstallMethod(StandardTableauxOp,[IsList],
  function(lam)
    local i;
    return SemiStandardTableauxOp(lam, List([1..Sum(lam)], i->1) );
  end
);

#F return the type of the tableau <tab>; note the slight complication
## because the composition which is the type of <tab> may contain zeros.
InstallMethod(TypeTableau,[IsTableau],
  function(tab)
    local t;

    t:=Flat(tab![1]);
    Append(t, [1..Maximum(t)]);
    t:=Collected(t);
    return t{[1..Length(t)]}[2]-1;
  end
);

#F return the shape of the tableau <tab>
InstallMethod(ShapeTableau,[IsTableau],tab->List(tab![1], Length));

