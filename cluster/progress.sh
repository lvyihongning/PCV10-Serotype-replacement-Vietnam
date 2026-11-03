#!/bin/bash
# ---------------------------------------------------------------
# Progress of a fit array (approaches 1-4 or simulations).
# Run on the cluster from the project directory:
#   bash "cluster/progress.sh" <array_job_id> [log_prefix]
# log_prefix is the start of the log file names: M1to4 (primary,
# default) or Sim_fit
# (high-carriage simulations, 20 tasks).
# or from the Mac:
#   ssh s-sc-frontend1 'cd "$HOME/PCV10 SR/final_code_pneumoalloc" && bash "cluster/progress.sh" <array_job_id>'
# ---------------------------------------------------------------
JOB=${1:?usage: progress.sh <array_job_id> [log_prefix]}
PREFIX=${2:-M1to4}

echo "== Slurm state"
sacct -j "$JOB" -X -n -o JobID%14,State%11,Elapsed%10,Timelimit%11

echo
echo "== Sampler progress (iteration of the slowest chain / warm-up + sampling)"
for f in logs/"${PREFIX}"_"${JOB}"_*.out; do
  [ -e "$f" ] || continue
  task=$(head -1 "$f" | sed 's/, [0-9]* CPUs on.*//')
  # slowest chain = lowest last-reported iteration across chains;
  # total iterations are read from the log ("Iteration: x / total")
  slow=$(grep -E "^Chain [0-9]+ Iteration" "$f" |
         awk '{c=$2; it=$4; tot=$6} {last[c]=it} END{m=""; for (c in last) if (m=="" || last[c]+0<m+0) m=last[c]; print (m==""?"-":m) "/" (tot==""?"?":tot)}')
  fin=$(grep -c "finished in" "$f")
  written=$(grep -h "^Written:" "${f%.out}.err" 2>/dev/null | sed 's/Written: results\///')
  printf "%-32s slowest chain at %10s | chains done %s/4 | %s\n" \
         "$task" "$slow" "$fin" "${written:-not written yet}"
done

echo
echo "== Errors"
grep -l -E "Execution halted|^Error" logs/"${PREFIX}"_"${JOB}"_*.err 2>/dev/null || echo "none"
