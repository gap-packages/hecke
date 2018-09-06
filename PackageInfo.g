#############################################################################
##
##  PackageInfo.g for the package hecke
##                                                            Dmitriy Traytel
##

SetPackageInfo( rec(

PackageName := "hecke",
Subtitle := "Hecke - Specht 2.4 ported to GAP 4",
Version := "1.4",

##  Release date of the current version in dd/mm/yyyy format.
Date := "02/07/2013",

ArchiveURL := Concatenation(
  "http://home.in.tum.de/~traytel/hecke/",
  "hecke1.4"),

ArchiveFormats := ".tar.gz",
BinaryFiles := ["doc/manual.pdf"],

Persons := [
  rec(
    LastName      := "Traytel",
    FirstNames    := "Dmitriy",
    IsAuthor      := true,
    IsMaintainer  := true,
    Email         := "traytel@in.tum.de",
    WWWHome       := "http://home.in.tum.de/~traytel/hecke/",
    Place         := "Munich",
    Institution   := "Technische Universität München"
  ),
],

Status := "deposited",

README_URL := Concatenation(
  "http://home.in.tum.de/~traytel/hecke/",
  "README.hecke"),
PackageInfoURL := Concatenation(
  "http://home.in.tum.de/~traytel/hecke/",
  "PackageInfo.g"),
SourceRepository := rec( 
  Type := "hg", 
  URL := "https://bitbucket.org/gap-system/hecke"
),

AbstractHTML :=
"The <span class=\"pkgname\">Hecke</span> package provides functions for \
calculating decomposition matrices of Hecke algebras of the symmetric groups \
and q-Schur algebras. Hecke is a port of the \
<span class=\"pkgname\">GAP 3</span> package \
<span class=\"pkgname\">Specht 2.4</span> to \
<span class=\"pkgname\">GAP 4</span>.",

PackageWWWHome := Concatenation(
  "http://home.in.tum.de/~traytel/hecke/",
  "index.html"),

PackageDoc := rec(
  BookName  := "hecke",
  ArchiveURLSubset := ["doc"],
  HTMLStart := "doc/chap0.html",
  PDFFile   := "doc/manual.pdf",
  SixFile   := "doc/manual.six",
  LongTitle := "Hecke - Specht 2.4 ported to GAP 4",
  Autoload  := true
),

Dependencies := rec(
  GAP := ">=4.4",
  NeededOtherPackages := [["GAPDoc", ">= 0.99"]],
  SuggestedOtherPackages := [],
  ExternalConditions := []
),

AvailabilityTest := function()
    return true;
  end,

#TestFile := "tst/testall.g",
Keywords := ["Hecke", "decomposition matrix", "Specht module", "Schur"],

));

