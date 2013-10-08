
echo " --- PREPARE --- "

echo "running on: $( hostname )"
echo "start time: $( date )    "

source $( $HOME/bashlib/myroot.sh )
source $HEPPY/mymeanalysis.git/setup.sh

set +f
cp -rv $HEPPY/mymeanalysis.git/megrid/share/* .
set -f
