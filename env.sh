# Will want a separate directory for each GW's background sample, but for now
# just testing for a single skymap with one.
outhistdir=/pnfs/nova/scratch/users/mstrait/ligobg/
outhadddir=/nova/ana/users/mstrait/ligobgresults/

if ! echo $PATH | grep -qE "$SRT_PRIVATE_CONTEXT/ligo([:/]|$)"; then
  PATH+=:$SRT_PRIVATE_CONTEXT/ligo
fi
