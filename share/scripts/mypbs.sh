
CLUSTER="lhc01 lhc02 lhc03 lhc04 lhc05 lhc06 lhc07 lhc08 lhc09 lhc10"

#####################################################################
# Evaluate a floating point number expression.
function float_eval()
{
    local stat=0
    local result=0.0
    if [[ $# -gt 0 ]]; then
        result=$(echo "scale=$float_scale; $*" | bc -q 2>/dev/null)
        stat=$?
        if [[ $stat -eq 0  &&  -z "$result" ]]; then stat=1; fi
    fi
    echo $result
    return $stat
}


#####################################################################
# Evaluate a floating point number conditional expression.
function float_cond()
{
    local cond=0
    if [[ $# -gt 0 ]]; then
        cond=$(echo "$*" | bc -q 2>/dev/null)
        if [[ -z "$cond" ]]; then cond=0; fi
        if [[ "$cond" != 0  &&  "$cond" != 1 ]]; then cond=0; fi
    fi
    local stat=$((cond == 0))
    return $stat
}

#####################################################################
qsub() {    
    for node in $CLUSTER; do
	load=( $( ssh dschoute@${node}.phys.sfu.ca "cat /proc/loadavg" ) )
	ncpu=( $( ssh dschoute@${node}.phys.sfu.ca "cat /proc/cpuinfo | grep -c processor" ) )
	if float_cond "${load[1]} < ${ncpu} / 4"; then
	    echo "good node: ${node} ${load[1]} ${ncpu}"
	fi
    done

    # ... launch jobs on good nodes 
}
