#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DEFAULT="Rscripts/species/run.R"   # <-- change this if you want a default
TAIL_LINES="${TAIL_LINES:-200}"           # how many lines to show before follow

SCRIPT="${1:-$SCRIPT_DEFAULT}"            # allow: ./run_job.sh path/to/script.R
LOGDIR="logs"
mkdir -p "$LOGDIR"                        # ensure logs folder exists

if [[ ! -f "$SCRIPT" ]]; then
  echo "ERROR: script not found: $SCRIPT"
  exit 1
fi

ts="$(date +%F_%H%M%S)"                   # timestamp
base="$(basename "$SCRIPT" .R)"           # script name without .R
log="$LOGDIR/${base}_${ts}.log"           # log file path
pidfile="$LOGDIR/${base}.pid"             # pid file path

nohup Rscript "$SCRIPT" > "$log" 2>&1 &   # run in background + log output
pid=$!                                    # capture PID
echo "$pid" > "$pidfile"                  # save PID
disown                                     # detach from this shell (logout-safe)

echo "$pid" > "$pidfile"                  # save PID
echo "Started: $SCRIPT"                   # short info
echo "PID: $pid  (saved in $pidfile)"     # short info
echo "Log: $log"                          # short info

tail -n "$TAIL_LINES" -f "$log"           # follow log (Ctrl+C stops tail only)