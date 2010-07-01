########################################################################
## TeX()  GAP routines for TeXing (polynomials) GAP objects           ##
##                                                                    ##
##        In the case of records, the general handler TeX() checks    ##
##        for an <X>.operations.TeX() command and calls this function ##
##        if possible. Otherwise special functions are provided here  ##
##        only for polynomials and (trivially) matrices (actually     ## 
##        LaTeX code).                                                ##
##                                                                    ##
##        These routines are included with the SPECHT package in      ##
##        order to allow "crystalized decomposition matrices" (which  ##
##        have polynomial entries), to be TeXed.                      ##
##                                                                    ##
##        (So not much in the way of documentation...)                ##
##                                                                    ##
##        Andrew Mathas                                               ##
##                                                                    ##
########################################################################

## July 1997
##   o fixed minor bug TeXing lists.

## We define TeX() below; for now we need to know that it exists.
TeX:=Ignore;

#P This function does the real work in TeXing polynomials. The functions 
## may not work well for polynomials with non-rational coefficients; I've 
## never checked. The variable term accounts for the different contexts in 
## which the polynomial can appear (as a coefficient etc.). 
##   term=0: no special treatment
##   term>0: polynomials are TeXed with surrounding brackets unless they 
##           are a power of the indeterminate; a + suffix is added as
##           nescessary when term>0 (in addition, if the leading term of 
##           the poynomial is negative then it is TeXed as -(-p) ).
## Traditionally (eg. with Print()), "term" are made by the user
## functions; I think this makes more sense (object oriented).
__TeXPolynomial:=function(poly, term)
  local i, val, coeffs, polyname, TeXCoeff;
 
  # printing the term a_n x^n where
  #   a_n = coeff
  #   n = val[uation]
  #   x = polyname
  # and len is the number of terms in the polynomials being
  # printed.
  TeXCoeff:=function(coeff, polyname, val, len)
    if not IsRat(coeff) then 
      TeX(coeff);
      if val = 1 then val:="";
      else val:=Concatenation("^{", String(val), "}");
      fi;
      Print(polyname, val);
    elif coeff <> 0 then
      if val = 0 then
        if (coeff <> 1 and coeff <> -1) or len > 1 or term = 0 then
          Print(coeff);
        elif coeff = -1 then Print("-");
        fi;
      else
        if val = 1 then val:="";
        else val:=Concatenation("^{", String(val), "}");
        fi;
        if coeff = -1 then Print("-", polyname, val);
        elif coeff = 1 then Print(polyname, val);
        else Print(coeff, polyname, val);
        fi;
      fi;
    fi;
  end; # TeXCoeff

  if not IsPolynomial(poly) then
    TeX(poly);
    if poly = 0 then return false;
    else return true;
    fi;
  elif Length(poly.coefficients) = 0 then  # poly = 0
    Print("0");
    return false;
  else
    if not IsBound(poly.baseRing.indeterminate.name) then polyname:="X";
    else polyname:=poly.baseRing.indeterminate.name;  
    fi;

    coeffs:=Copy(poly.coefficients);
    coeffs:=coeffs{[Length(coeffs),Length(coeffs)-1..1]};
    val:=poly.valuation+Length(coeffs)-1;
    if Length(coeffs) > 1 and term <> 0 then
      if coeffs[1] < 0 then
        coeffs:=List([1..Length(coeffs)], i -> -coeffs[i]);
        Print(" - (");
      elif term > 1 then Print(" + (");
      else Print("(");
      fi;
    elif term > 1 then
      if coeffs[1] > 0 then Print( " + ");
      else Print(" ");
      fi;
    fi;
    TeXCoeff(coeffs[1], polyname, val, Length(coeffs));
    if Length(coeffs) > 1 then
      for i in [2..Length(coeffs)] do
        if coeffs[i] > 0 then Print("+"); fi;
        val:=val - 1;
        TeXCoeff(coeffs[i], polyname, val, Length(coeffs));
      od;
      if term <> 0 then Print(")"); fi;
    fi;
    return true;
  fi;
end;      # __TeXPolynomial

#P TeXPolynomial is a front end to  __TeXPolynomial()
TeXPolynomial:=function(arg) local poly;

  if Length(arg)=1 then __TeXPolynomial(arg[1],0);
  else
   for poly in arg do
     __TeXPolynomial(poly,0);
     Print("\n");
   od;
  fi;
end;

#P LaTeX code for a matrix
TeXMatrix:=function(M)
  local i, j;

  Print("\\left(\\begin{array}{*{", Length(M[1]),"}{l}}\n");
  for i in [1..Length(M)] do
    for j in [1..(Length(M[i])-1)] do
      TeX(M[i][j]);
      Print("& ");
      if ( j mod 10 = 0 ) then Print("\n"); fi;
    od;
    TeX(M[i][Length(M[i])]);
    Print("\\\\\n");
  od;
  Print("\\end{array}\\right)\n");
end;

#P Usage: TeXWideMatrix(m1 [,m2, ...], C)
##   C is the number of columns per page. All of the matrices given to 
## TeXWideMatrix() are printed on the same page (if this is possible).
## *** undocumented
TeXWideMatrix:=function(arg)
  local mats, m, C, r, c, col, cmax, arraytop;

  mats:=arg{[1..Length(arg)-1]};  ## matrices
  if mats=[] or not IsInt(arg[Length(arg)])  then
    PrintTo("*errout*", "usage: TeXWideMatrix(m1 [,m2...], col)\n");
    return;
  fi;
  C:=arg[Length(arg)];            ## columns per page

  cmax:=Maximum(List(mats, m->Maximum(List(m, r->Length(r)))));
  arraytop:=Concatenation("\\begin{array}{*{",String(cmax),"}{l}}\n");
  col:=[1-C..0];
  while col[C]<cmax do
    col:=col+C;
    for m in mats do
      if col[1]=1 then Print("$\\left(");
      else Print("$\\left.");
      fi;
      Print(arraytop);
      for r in [1..Length(m)] do
        if IsBound(m[r][col[1]]) then 
          TeX(m[r][col[1]]); 
          for c in col{[2..C]} do
            if c<=Length(m[r]) then TeX("& ",m[r][c]); fi;
          od;
        fi;
        Print("\\\\\n");
      od;
      if col[C]>=cmax then Print("\\end{array}\\right)$\n\n\n");
      else Print("\\end{array}\\right.$\n\n\\medskip\n\n");
      fi;
    od;
    if col[C]<cmax then Print("\\newpage\n\n");fi;
  od;
end;


#P The TeX() handler; calls the above routines and checks operations record
TeX:=function(arg) local a, i;
  for a in arg do
    if IsRat(a) or IsString(a) then Print(a);
    elif IsMatrix(a) then TeXMatrix(a);
    elif IsList(a) then 
      if a=[] then Print("[]"); 
      else
        Print("[ "); TeX(a[1]);
        for i in [2..Length(a)] do 
          Print(", "); TeX(a[i]); 
        od;
        Print(" ]");
      fi;
    elif IsPerm(a) then Print(a);
    elif IsPolynomial(a) then TeXPolynomial(a);
    elif IsRec(a) and IsBound(a.operations) and IsBound(a.operations.TeX)
    then a.operations.TeX(a);
    else PrintTo("*errout*", " *** error: ",
                 "TeX(<a>), don't know how to TeX <a>.\n");
    fi;
  od;
end;

