LoadPackage( "hecke" );

# Some doc examples call SaveDecompositionMatrix, which writes files such as
# "e4p0.10" into the current directory; run from a scratch directory that GAP
# removes on exit, so the working tree stays clean.
ChangeDirectoryCurrent( Filename( DirectoryTemporary(), "" ) );

TestDirectory( DirectoriesPackageLibrary("hecke", "tst"), rec(exitGAP := true ) );
FORCE_QUIT_GAP(1);
