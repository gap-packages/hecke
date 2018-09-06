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

Persons := [
  rec(
    LastName      := "Traytel",
    FirstNames    := "Dmitriy",
    IsAuthor      := true,
    IsMaintainer  := false,
    Email         := "traytel@in.tum.de",
    WWWHome       := "http://home.in.tum.de/~traytel/hecke/",
    Place         := "Munich",
    Institution   := "Technische Universität München"
  ),

  rec(
    LastName      := "GAP Team",
    FirstNames    := "The",
    IsAuthor      := false,
    IsMaintainer  := true,
    Email         := "support@gap-system.org",
  ),
],

Status := "deposited",

PackageWWWHome  := "https://gap-packages.github.io/hecke/",
README_URL      := Concatenation( ~.PackageWWWHome, "README.md" ),
PackageInfoURL  := Concatenation( ~.PackageWWWHome, "PackageInfo.g" ),
SourceRepository := rec(
    Type := "git",
    URL := "https://github.com/gap-packages/hecke",
),
IssueTrackerURL := Concatenation( ~.SourceRepository.URL, "/issues" ),
ArchiveURL      := Concatenation( ~.SourceRepository.URL,
                                 "/releases/download/v", ~.Version,
                                 "/hecke-", ~.Version ),
ArchiveFormats := ".tar.gz",

AbstractHTML :=
"The <span class=\"pkgname\">Hecke</span> package provides functions for \
calculating decomposition matrices of Hecke algebras of the symmetric groups \
and q-Schur algebras. Hecke is a port of the \
<span class=\"pkgname\">GAP 3</span> package \
<span class=\"pkgname\">Specht 2.4</span> to \
<span class=\"pkgname\">GAP 4</span>.",

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
  GAP := ">=4.8",
  NeededOtherPackages := [],
  SuggestedOtherPackages := [],
  ExternalConditions := []
),

AvailabilityTest := function()
    return true;
  end,

#TestFile := "tst/testall.g",
Keywords := ["Hecke", "decomposition matrix", "Specht module", "Schur"],

));

