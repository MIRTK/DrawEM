# Available atlases from Draw-EM
atlases=""
for d in `find $DRAWEMDIR/parameters/* -maxdepth 1 -type d`;do
  atlases="$atlases "`basename $d`;
done
export AVAILABLE_ATLASES=$atlases

export SUPER_GM_LABEL=10000
export SUPER_WM_LABEL=20000

# log function
run()
{
  echo "$@"
  "$@"
  if [ ! $? -eq 0 ]; then
    echo "$@ : failed"
    # exit 1
  fi
}

# make run function global
typeset -fx run
