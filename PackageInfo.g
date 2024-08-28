#############################################################################
##
##  PackageInfo.g for the package hecke
##                                                            Dmitriy Traytel
##

SetPackageInfo( rec(

PackageName := "hecke",
Subtitle := "Calculating decomposition matrices of Hecke algebras",
Version := "1.5.4",
Date := "27/08/2024", # dd/mm/yyyy format
License := "GPL-2.0-or-later",

Persons := [
  rec(
    LastName      := "Traytel",
    FirstNames    := "Dmitriy",
    IsAuthor      := true,
    IsMaintainer  := false,
    Email         := "traytel@in.tum.de",
    WWWHome       := "https://home.in.tum.de/~traytel/hecke/",
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
  HTMLStart := "doc/chap0_mj.html",
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

TestFile := "tst/testall.g",
Keywords := ["Hecke", "decomposition matrix", "Specht module", "Schur"],


  AutoDoc := rec(
      TitlePage := rec(
          Copyright := """
            &copyright; 2010&ndash;2013 by Dmitriy Traytel<P/>

            This package may be distributed under the terms and conditions of the
            GNU Public License Version 2 or higher.
            """,
          Acknowledgements := """
            &Specht; is a port of the &GAP; 3 package <Package>Specht</Package> 2.4 to &GAP; 4.
            <Package>Specht</Package> 2.4 was written by Andrew Mathas, who allowed
            Dmitriy Traytel to use his source code as the basis for &specht;.
            """,
      ),
  ),

));

