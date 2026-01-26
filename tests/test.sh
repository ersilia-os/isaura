#!/usr/bin/env bash
set -euo pipefail

MODEL_ID="eos3b5e"
MODEL_VERSION="v1"
TEST_BUCKET="test-bucket"
READ_BUCKET="isaura-public"

MINIO_URL="http://localhost:9000"
MINIO_HEALTH="$MINIO_URL/minio/health/ready"

MILVUS_INFO_URL="http://localhost:8080/info?collection=${MODEL_ID}_${MODEL_VERSION}"
MILVUS_PORTS=(19530 9091 2379)

READ_RUNS=3

need(){ command -v "$1" >/dev/null 2>&1 || { echo "Missing: $1" >&2; exit 127; }; }
need isaura; need curl; need jq; need python

INPUT_CSV="${1:-}"
[[ -n "$INPUT_CSV" && -f "$INPUT_CSV" ]] || { echo "Usage: $0 /path/to/input.csv" >&2; exit 2; }

TMPDIR="$(mktemp -d -t isaura_smoketest_XXXXXX)"
trap 'rm -rf "$TMPDIR"' EXIT

STATS_DIR="$TMPDIR/stats"; mkdir -p "$STATS_DIR"
BASELINE_OUT="$TMPDIR/read_baseline.csv"
RUN_OUT_PREFIX="$TMPDIR/read_run"
READ_METRICS_TSV="$TMPDIR/read_metrics.tsv"; : >"$READ_METRICS_TSV"
INSPECT_OUT="$TMPDIR/inspect_available.csv"
MILVUS_INFO_JSON="$TMPDIR/milvus_info.json"
TESTLOG="$TMPDIR/testlog.tsv"; : >"$TESTLOG"

add_row(){ printf "%s\t%s\t%s\n" "$1" "$2" "$3" >>"$TESTLOG"; }

ms_now(){ python - <<'PY'
import time
print(int(time.time()*1000))
PY
}

csv_rows(){ python - "$1" <<'PY'
import csv,sys
p=sys.argv[1]
with open(p,newline='',encoding='utf-8') as f:
  r=csv.reader(f); next(r,None)
  print(sum(1 for _ in r))
PY
}

timed_read(){
  local run_idx="$1" out_csv="$2"
  local t0 t1 ms rows rps
  t0="$(ms_now)"
  isaura read -m "$MODEL_ID" -v "$MODEL_VERSION" -pn "$READ_BUCKET" -i "$INPUT_CSV" -o "$out_csv" >/dev/null 2>&1
  t1="$(ms_now)"
  ms=$((t1 - t0))
  rows="$(csv_rows "$out_csv")"
  rps="$(python - "$ms" "$rows" <<'PY'
import sys
ms=int(sys.argv[1]); rows=int(sys.argv[2])
print(0.0 if ms<=0 else rows/(ms/1000.0))
PY
)"
  printf "%s\t%s\t%s\t%s\n" "$run_idx" "$ms" "$rows" "$rps" >>"$READ_METRICS_TSV"
}

compare_csvs_strict(){
  python - "$1" "$2" <<'PY'
import csv,sys
a,b=sys.argv[1],sys.argv[2]
def read(path):
  with open(path,newline='',encoding='utf-8') as f:
    r=csv.reader(f); hdr=next(r,None); rows=list(r)
  return hdr, rows
ha,ra=read(a); hb,rb=read(b)
if ha!=hb: raise SystemExit(2)
if ra!=rb: raise SystemExit(3)
PY
}

check_input_present_in_output(){
  python - "$1" "$2" <<'PY'
import csv,sys
inc,outc=sys.argv[1],sys.argv[2]
def load_inputs(p):
  with open(p,newline='',encoding='utf-8') as f:
    r=csv.DictReader(f); vals=[]
    for row in r:
      v=(row.get("input") or row.get("smiles") or "").strip()
      if v: vals.append(v)
  return vals
def load_out(p):
  with open(p,newline='',encoding='utf-8') as f:
    r=csv.DictReader(f); fn=r.fieldnames or []
    key="input" if "input" in fn else ("smiles" if "smiles" in fn else None)
    if not key: return None,set()
    s=set()
    for row in r:
      v=(row.get(key) or "").strip()
      if v: s.add(v)
    return key,s
ins=load_inputs(inc)
key,outs=load_out(outc)
if key is None: raise SystemExit(4)
missing=[v for v in ins if v not in outs]
if missing: raise SystemExit(5)
PY
}

compare_csv_cols_from_third(){
  python - "$1" "$2" <<'PY'
import csv,sys
a,b=sys.argv[1],sys.argv[2]
def load(path):
  with open(path,newline='',encoding='utf-8') as f:
    r=csv.reader(f); hdr=next(r,None) or []; rows=list(r)
  return hdr, rows
ha,ra=load(a); hb,rb=load(b)
if ha!=hb: raise SystemExit(2)
if len(ha) < 3: raise SystemExit(0)
idx=list(range(2,len(ha)))
def proj(rows):
  return [[(row[i] if i < len(row) else "") for i in idx] for row in rows]
if proj(ra) != proj(rb): raise SystemExit(3)
PY
}

step(){
  local name="$1"; shift
  "$@" >/dev/null 2>&1 && add_row "$name" "PASS" "" || { add_row "$name" "FAIL" ""; return 1; }
}

check_http(){
  local url="$1" name="$2"
  curl -fsS "$url" >/dev/null 2>&1 && add_row "$name" "PASS" "$url" || add_row "$name" "FAIL" "$url"
}

check_tcp(){
  local port="$1"
  (echo >/dev/tcp/localhost/"$port") >/dev/null 2>&1 && add_row "Milvus TCP ${port}" "PASS" "localhost:${port}" || add_row "Milvus TCP ${port}" "FAIL" "localhost:${port}"
}

check_http "$MINIO_HEALTH" "MinIO health"
for p in "${MILVUS_PORTS[@]}"; do check_tcp "$p"; done

python - "$MILVUS_INFO_URL" "$MILVUS_INFO_JSON" <<'PY'
import sys, json, subprocess, time
from rich.console import Console
from rich.spinner import Spinner
from rich.live import Live

url, outp = sys.argv[1], sys.argv[2]
console = Console()
spinner = Spinner("dots", text=f"curl {url}")
with Live(spinner, refresh_per_second=12, console=console):
  p = subprocess.run(["curl","-fsS","-X","POST",url], capture_output=True, text=True)
  time.sleep(0.12)

if p.returncode != 0 or not p.stdout.strip():
  sys.exit(1)

try:
  data = json.loads(p.stdout)
except Exception:
  sys.exit(2)

with open(outp, "w", encoding="utf-8") as f:
  json.dump(data, f, indent=2)

sys.exit(0)
PY
if [[ -s "$MILVUS_INFO_JSON" ]]; then add_row "Milvus info" "PASS" "$MILVUS_INFO_URL"; else add_row "Milvus info" "WARN" "$MILVUS_INFO_URL"; fi

step "isaura write" isaura write -m "$MODEL_ID" -v "$MODEL_VERSION" -pn "$TEST_BUCKET" --access public -i "$INPUT_CSV"
step "isaura copy" isaura copy -m "$MODEL_ID" -v "$MODEL_VERSION" -pn "$TEST_BUCKET"

for i in $(seq 1 "$READ_RUNS"); do
  out="${RUN_OUT_PREFIX}_${i}.csv"
  timed_read "$i" "$out"
  [[ "$i" -eq 1 ]] && cp "$out" "$BASELINE_OUT"
done

if check_input_present_in_output "$INPUT_CSV" "$BASELINE_OUT"; then add_row "read#1 inputs present" "PASS" ""; else add_row "read#1 inputs present" "FAIL" ""; fi
for i in $(seq 2 "$READ_RUNS"); do
  out="${RUN_OUT_PREFIX}_${i}.csv"
  if compare_csvs_strict "$BASELINE_OUT" "$out"; then add_row "read#${i} equals baseline" "PASS" ""; else add_row "read#${i} equals baseline" "FAIL" ""; fi
  if compare_csv_cols_from_third "$BASELINE_OUT" "$out"; then add_row "read#${i} cols[3..] match" "PASS" ""; else add_row "read#${i} cols[3..] match" "FAIL" ""; fi
done

isaura inspect -m "$MODEL_ID" -v "$MODEL_VERSION" -pn "$READ_BUCKET" --access public -o "$INSPECT_OUT" >/dev/null 2>&1 \
  && add_row "isaura inspect" "PASS" "" || add_row "isaura inspect" "WARN" ""

isaura stats -pn "$READ_BUCKET" --access public -o "$STATS_DIR" >/dev/null 2>&1 \
  && add_row "isaura stats" "PASS" "" || add_row "isaura stats" "WARN" ""

step "remove isaura-public" isaura remove -m "$MODEL_ID" -v "$MODEL_VERSION" -pn "$READ_BUCKET" --yes
step "move test-bucket" isaura move -m "$MODEL_ID" -v "$MODEL_VERSION" -pn "$TEST_BUCKET"
step "remove isaura-public (post-move)" isaura remove -m "$MODEL_ID" -v "$MODEL_VERSION" -pn "$READ_BUCKET" --yes

python - "$TESTLOG" "$READ_METRICS_TSV" "$MILVUS_INFO_JSON" "$MODEL_ID" "$MODEL_VERSION" "$TEST_BUCKET" "$READ_BUCKET" "$INPUT_CSV" "$TMPDIR" <<'PY'
import sys, csv, json, statistics
from rich.console import Console
from rich.table import Table
from rich import box
from rich.text import Text

testlog, metrics, milvus_json, mid, mv, tb, rb, incsv, tmpdir = sys.argv[1:]
console = Console()

colors = {"PASS":"green", "FAIL":"red", "WARN":"yellow"}
def cstat(s):
  s=s.upper()
  return f"[{colors.get(s,'white')}]{s}[/{colors.get(s,'white')}]"

def pctl(vals, p):
  if not vals: return None
  s=sorted(vals)
  i=max(0, min(len(s)-1, int(round((p/100)*(len(s)-1)))))
  return s[i]

perf=[]
with open(metrics, newline="", encoding="utf-8") as f:
  r=csv.reader(f, delimiter="\t")
  for row in r:
    if not row: continue
    run, ms, rows, rps = row
    perf.append((int(run), int(ms), int(rows), float(rps)))

ms_vals=[p[1] for p in perf] if perf else []
rows_val=perf[0][2] if perf else 0

perf_tbl = Table(box=box.SIMPLE_HEAVY, show_header=True, pad_edge=False)
perf_tbl.add_column("run", justify="right", style="bold")
perf_tbl.add_column("ms", justify="right", style="cyan")
perf_tbl.add_column("rows", justify="right", style="magenta")
perf_tbl.add_column("rows/s", justify="right", style="green")
for r,ms,rows,rps in perf:
  perf_tbl.add_row(str(r), str(ms), str(rows), f"{rps:.2f}")

perf_meta = Table(box=box.SIMPLE, show_header=False, pad_edge=False)
perf_meta.add_column("k", style="dim", no_wrap=True)
perf_meta.add_column("v")
if ms_vals:
  perf_meta.add_row("avg ms", f"[cyan]{statistics.mean(ms_vals):.1f}[/]")
  perf_meta.add_row("p95 ms", f"[cyan]{pctl(ms_vals,95)}[/]")
  perf_meta.add_row("rows", f"[magenta]{rows_val}[/]")
else:
  perf_meta.add_row("status", "[yellow]no data[/]")

perf_wrap = Table(box=box.MINIMAL, show_header=False, pad_edge=False)
perf_wrap.add_column("x")
perf_wrap.add_row(perf_tbl)
perf_wrap.add_row(perf_meta)

milvus_tbl = Table(box=box.SIMPLE_HEAVY, show_header=False, pad_edge=False)
milvus_tbl.add_column("k", style="dim", no_wrap=True)
milvus_tbl.add_column("v")

milvus_status="WARN"
try:
  with open(milvus_json, "r", encoding="utf-8") as f:
    d=json.load(f)
  idx=d.get("index") or {}
  milvus_status="PASS"
  milvus_tbl.add_row("collection", f"[bold]{d.get('collection')}[/]")
  milvus_tbl.add_row("row_count", f"[magenta]{d.get('row_count')}[/]")
  milvus_tbl.add_row("dim_bits", f"[cyan]{d.get('dim_bits')}[/]")
  milvus_tbl.add_row("estimated_bytes", f"[cyan]{d.get('estimated_bytes')}[/]")
  milvus_tbl.add_row("load_state", f"[cyan]{d.get('load_state')}[/]")
  milvus_tbl.add_row("loading_progress", f"[green]{d.get('loading_progress')}%[/]")
  milvus_tbl.add_row("index.exists", f"[green]{idx.get('exists')}[/]")
  milvus_tbl.add_row("index.state", f"[green]{idx.get('state')}[/]")
  milvus_tbl.add_row("index.progress", f"[green]{idx.get('progress_pct')}%[/]")
  milvus_tbl.add_row("indexed_rows", f"[magenta]{idx.get('indexed_rows')}[/] / [magenta]{idx.get('total_rows')}[/]")
  milvus_tbl.add_row("server_time", f"[dim]{d.get('server_time_rfc3339')}[/]")
except Exception:
  milvus_tbl.add_row("milvus", "[yellow]info not available[/]")

main = Table(title="ISAURA CLI SMOKE TEST REPORT", box=box.MINIMAL_HEAVY_HEAD, show_lines=True)
main.add_column("Step", style="bold", no_wrap=True)
main.add_column("Status", justify="center", no_wrap=True)
main.add_column("Details", overflow="fold")

context = (
  f"model: [bold]{mid}/{mv}[/]\n"
  f"write/copy/move bucket: [bold]{tb}[/]\n"
  f"read bucket: [bold]{rb}[/]\n"
  f"input csv: [dim]{incsv}[/]\n"
  f"work dir: [dim]{tmpdir}[/]"
)
main.add_row("Context", "[blue]INFO[/blue]", context)
main.add_row("Milvus info (summary)", cstat(milvus_status), milvus_tbl)
main.add_row("Read Performance", "[blue]INFO[/blue]", perf_wrap)

stats={"PASS":0,"FAIL":0,"WARN":0}
with open(testlog, newline="", encoding="utf-8") as f:
  r=csv.reader(f, delimiter="\t")
  for step, status, details in r:
    s=status.upper()
    stats[s]=stats.get(s,0)+1
    if step == "Milvus info":
      continue
    main.add_row(step, cstat(s), details or "")

console.print(main)

msg = Text()
msg.append("Summary: ", style="bold")
msg.append(f"{stats.get('PASS',0)} PASS", style="green")
msg.append("  ")
msg.append(f"{stats.get('WARN',0)} WARN", style="yellow")
msg.append("  ")
msg.append(f"{stats.get('FAIL',0)} FAIL", style="red")
console.print(msg)

sys.exit(1 if stats.get("FAIL",0) else 0)
PY
