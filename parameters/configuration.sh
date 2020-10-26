# Available atlases from Draw-EM
export AVAILABLE_ATLASES="ALBERT MCRIB"
export AVAILABLE_TISSUE_ATLASES="neonatal"

export SUPER_GM_LABEL=10000
export SUPER_WM_LABEL=20000

# log function
run()
{
  echo "$@"
  "$@"
  if [ ! $? -eq 0 ]; then
    echo "$@ : failed"
    exit 1
  fi
}

# make run function global
typeset -fx run
