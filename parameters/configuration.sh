# Available atlases from Draw-EM
AVAILABLE_ATLASES=""
for d in `find $DRAWEMDIR/parameters/* -maxdepth 1 -type d`;do
  AVAILABLE_ATLASES="$AVAILABLE_ATLASES "`basename $d`;
done
