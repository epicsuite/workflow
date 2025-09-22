#!/usr/bin/env nextflow
nextflow.enable.dsl=2
import groovy.yaml.YamlSlurper

/*====================== PARAMS ======================*/
params.yaml   = params.yaml ?: null
params.inp    = params.inp  ?: null
params.outdir = params.outdir ?: 'results'

/*================== YAML: parse once =================*/
def yaml_path = (params.yaml ?: params.inp ?: 'input.yaml') as String
def yaml_file = new File(yaml_path).getCanonicalFile()
assert yaml_file.exists() && yaml_file.isFile() : "YAML not found: ${yaml_file}"

def cfg = new YamlSlurper().parse(yaml_file)
assert cfg instanceof Map : "Top-level YAML must be a mapping"

def ens         = cfg.ensemble ?: [:]
def experiments = (ens.experiments as List) ?: []
assert experiments, "YAML missing or empty: ensemble.experiments"

/*=========== Reference ===========*/
def ref_block      = ens.reference ?: [:]
def ref_sequence   = (ref_block.sequence     ?: '').toString()
def ref_mito       = (ref_block.mitochondria ?: '').toString()
def ref_resolution = (ref_block.resolution   ?: '').toString()
assert ref_sequence,   "reference.sequence missing in YAML"
//assert ref_mito,       "reference.mitochondria missing in YAML"
assert ref_resolution, "reference.resolution missing in YAML"

/*----------- Chromosomes: TWO separate lists ----------*/
def chroms_block     = (ref_block.chromosomes ?: [:]) as Map
def genomelist1_path = chroms_block.genomelist1 ? chroms_block.genomelist1.toString() : null
def genomelist2_raw  = chroms_block.genomelist2 ? chroms_block.genomelist2.toString() : null
assert genomelist1_path : "reference.chromosomes.genomelist1 must be provided"
assert genomelist2_raw  : "reference.chromosomes.genomelist2 must be provided"

/* genomelist1: pass PATH/URL string through (do NOT read contents) */
def chroms_p1_path = genomelist1_path
println "[CHROMS] P1 path=${chroms_p1_path}"

/* genomelist2: file OR comma-separated list → expand to IDs */
def readChromFile = { String path ->
  def f = new File(path)
  assert f.exists() && f.isFile() : "Chromosome list not found: ${path}"
  f.readLines()
   .collect { it.trim() }
   .findAll { it && !it.startsWith('#') }
   .collect { line -> line.split(/\s+/)[0] }
   .findAll { it }
}

List<String> chroms_p2
def gl2 = new File(genomelist2_raw)
if (gl2.exists() && gl2.isFile()) {
  chroms_p2 = readChromFile(genomelist2_raw)
  println "[CHROMS] P2 count=${chroms_p2.size()} from file ${genomelist2_raw}"
}
else {
  chroms_p2 = genomelist2_raw
                .split(/\s*,\s*/)
                .collect { it.trim() }
                .findAll { it }
  println "[CHROMS] P2 count=${chroms_p2.size()} from comma-separated list"
}
assert chroms_p2, "genomelist2 produced no chromosome IDs"

/*====== helpers ======*/
def firstNonNull = { Map m, List<String> keys ->
  for (k in keys) { def v = m[k]; if (v != null && v.toString().trim()) return v.toString() }
  return null
}
def joinPath = { String base, String rel ->
  if (!base) return rel
  if (!rel)  return base
  new File(base, rel).getPath()
}
def structure_root = (ens.structure_root ?: ens.root ?: ens.base ?: null)?.toString()

/*====== Build EXACT 6-field tuples: (exp, ts, structure, ref, mito, res) ======*/
List<List> rows = []
experiments.each { e ->
  def exp_name = firstNonNull(e as Map, ['name','id','exp'])
  assert exp_name : "Experiment missing a name/id/exp"

  def exp_struct = firstNonNull(e as Map, ['structure','dir','path','url'])
  if (exp_struct && structure_root) exp_struct = joinPath(structure_root, exp_struct)

  def tss_any = (e.timesteps ?: e.steps ?: e.timepoints)
  assert tss_any != null : "Experiment '${exp_name}' has no timesteps"

  if (tss_any instanceof Map) {
    tss_any.each { k, v ->
      def ts_name = k?.toString()
      String ts_struct = (v instanceof Map) ? firstNonNull(v as Map, ['structure','dir','path','url']) : null
      if (!ts_struct) ts_struct = exp_struct
      if (ts_struct && structure_root) ts_struct = joinPath(structure_root, ts_struct)
      assert ts_struct : "Experiment '${exp_name}' timestep '${ts_name}' missing structure"
      rows << [exp_name, ts_name, ts_struct, ref_sequence, ref_mito, ref_resolution]
    }
  }
  else if (tss_any instanceof List) {
    tss_any.eachWithIndex { tsObj, i ->
      String ts_name, ts_struct
      if (tsObj instanceof Map) {
        ts_name   = firstNonNull(tsObj as Map, ['name','id','timestep','timepoint','ts'])
        ts_struct = firstNonNull(tsObj as Map, ['structure','dir','path','url'])
      } else {
        ts_name   = tsObj?.toString()
      }
      if (!ts_struct) ts_struct = exp_struct
      if (ts_struct && structure_root) ts_struct = joinPath(structure_root, ts_struct)
      assert ts_name   : "Experiment '${exp_name}' timestep #${i+1} missing name"
      assert ts_struct : "Experiment '${exp_name}' timestep '${ts_name}' missing structure"
      rows << [exp_name, ts_name, ts_struct, ref_sequence, ref_mito, ref_resolution]
    }
  }
  else {
    assert false : "Experiment '${exp_name}' timesteps must be a list or a map; got: ${tss_any?.getClass()?.name}"
  }
}

if (!rows || !rows.every { it instanceof List && it.size()==6 }) {
  def bad = rows.find { !(it instanceof List) || it.size()!=6 }
  throw new IllegalStateException("rows contains non-6 item: ${bad} (type=${bad?.getClass()?.name}); rows.size=${rows?.size()}")
}
println "[YAML] Built ${rows.size()} tuples (exp, ts, structure, ref, mito, res)"
rows.take(5).each { println "[YAML] row: ${it}" }

/* ---- Channel of true 6-field Nextflow tuples ---- */
def TSTEPS6 = Channel
  .from(rows)
  .map { lst -> tuple(lst[0], lst[1], lst[2], lst[3], lst[4], lst[5]) }

/*========================= WORKFLOW =========================*/
workflow {

  /* 0) BWA index once */
  Channel.fromPath(ref_sequence, checkIfExists: true).set { REF_CH }
  def idx_done = bwa_index(REF_CH).bwa_index_done

  /* Gate: combine adds token as 7th; keep only the first 6 fields */
  def gated_tsteps = TSTEPS6
    .combine(idx_done)                         // -> (exp, ts, structure, ref, mito, res, token)
    .map { t -> tuple(t[0], t[1], t[2], t[3], t[4], t[5]) }

  // Debug: shows ts like 12hpi / 24hpi
  gated_tsteps.view { t -> "GATED_TSTEPS -> ts=${t[1]} full=${t}" }

  /* 1) Provide genomelist1 PATH (string) to slurpy_hic */
  def CHROMS_P1_PATH_CH = Channel.value(chroms_p1_path)
  def std_input = gated_tsteps
    .combine(CHROMS_P1_PATH_CH)                // -> (exp, ts, structure, ref, mito, res, chroms_p1_path)
    .map { t -> tuple(t[0], t[1], t[2], t[3], t[4], t[5], t[6]) }

  def std_ch = slurpy_hic(std_input)           // -> (exp, ts, res, file 'hicfile.hic')

  /* 2) Fan-out PER CHROMOSOME using ONLY P2 list (flatMap + collect) */
  def per_chr_ch = std_ch.flatMap { left ->
    def (exp, ts, res, stdfile) = left
    chroms_p2.collect { chr -> tuple(exp, ts, res, stdfile, chr as String) }
  }

  /* 3) One hic2structure per chromosome (parallel) */
  def per_chr_dirs = hic2struct_one(per_chr_ch) // -> (exp, ts, res, chr, path 'hic_<chr>')

//  /* 4) Group by (exp, ts, res) for concat */
//  def grouped = per_chr_dirs
//    .map { t -> tuple([t[0], t[1], t[2]] as List, t[4]) }   // ([exp,ts,res], hicdir)
//    .groupTuple()
//    .map { pair ->
//      def key  = pair[0] as List
//      def dirs = pair[1]
//      tuple(key[0] as String, key[1] as String, key[2] as String, dirs)
//    }
//
  /* 4) Group by (exp, ts) for concat */
  def grouped = per_chr_dirs
    .map { t -> tuple([t[0], t[1]] as List, [t[2], t[3], t[4]]) } 
    // ([exp,ts], [res, chr, hicdir])
    .groupTuple()
    .map { pair ->
      def key   = pair[0] as List    // [exp, ts]
      def items = pair[1] as List    // [[res, chr, hicdir], ...]
      def res   = items[0][0]        // resolution is same for all
      def dirs  = items.collect { it[2] } // gather hic_<chr> dirs
      tuple(key[0], key[1], res, dirs)
    }
  /* 5) Final concat and publish (uses P2 list) */
  hic_concat(grouped, Channel.value(chroms_p2))
  
  /* 6) Save input file for data provenance */
  // Make channels for the two inputs
  def YAML_CH = Channel.fromPath(params.yaml, checkIfExists: true)
  def GL1_CH  = Channel.fromPath(genomelist1_path, checkIfExists: true)

  // Invoke the process with those channels
  save_inp(YAML_CH, GL1_CH)
}

/*======================== PROCESSES ========================*/

process bwa_index {
  tag "bwa_index"
  input:
  path reference
  output:
  path "bwa_index.done", emit: bwa_index_done
  path "bwa_index.log",  emit: bwa_index_log
  script:
  """
  set -euo pipefail
  ref_abs=\$(readlink -f "${reference}")
  ref_dir=\$(dirname -- "\$ref_abs")
  ref_bn=\$(basename -- "\$ref_abs")

  have_index() {
    local dir="\$1"; local bn="\$2"
    local -a exts=(amb ann bwt pac sa)
    for ext in "\${exts[@]}"; do
      if ! find -L "\$dir" -maxdepth 1 -type f -name "\$bn.\$ext" -size +0c -print -quit | grep -q . ; then
        return 1
      fi
    done
    return 0
  }

  : > bwa_index.log
  {
    echo "[INFO] Index target : \$ref_abs"
    if have_index "\$ref_dir" "\$ref_bn"; then
      echo "[INFO] BWA index already present — skipping build."
    else
      echo "[INFO] Running bwa index..."
      bwa index "\$ref_abs"
      have_index "\$ref_dir" "\$ref_bn" || { echo "[ERROR] Index missing after run" >&2; exit 1; }
    fi
  } | tee -a bwa_index.log

  echo "\$ref_abs" > bwa_index.done
  """
}

//process slurpy_hic {
//  tag "${exp}/${ts}"
//  input:
//  tuple val(exp), val(ts), val(structure), val(reference), val(mitochondria), val(ref_resolution), val(chroms_p1_path)
//  output:
//  tuple val(exp), val(ts), val(ref_resolution), file('hicfile.hic')
//  script:
//  """
//  set -euo pipefail
//  echo "[slurpy_hic] exp=${exp} ts=${ts} structure=${structure}" >&2
//  echo "[slurpy_hic] using genomelist1 PATH: ${chroms_p1_path}" >&2
//  ln -s "/panfs/biopan04/4DGENOMESEQ/EDGE_WORKFLOW/workflow/nextflow/SLURPY/" ./SLURPY
//  ln -s "${structure}" ./fastq
//  # run SLUR-(M)-PY
//  /panfs/biopan04/4DGENOMESEQ/EDGE_WORKFLOW/workflow/nextflow/SLURPY/slurm.py -r ${reference} -P fast,tb,gpu -M ${mitochondria} --fq ${structure} -G ${chroms_p1_path} -F 150000 15000000 -J /panfs/biopan04/4DGENOMESEQ/EDGE_WORKFLOW/workflow/nextflow/SLURPY/juicer_tools_1.22.01.jar
//  watch_pattern="aligned/*valid*hic"
//  echo "[slurpy_hic] waiting for pattern: \${watch_pattern}" >&2
//
//  get_mtime() { stat -c %Y "\$1" 2>/dev/null || stat -f %m "\$1" 2>/dev/null || echo 0; }
//  get_size()  { stat -c %s "\$1" 2>/dev/null || wc -c < "\$1"; }
//
//  target=""
//  while :; do
//    files=( \$(ls -1 "\$watch_pattern" 2>/dev/null || true) )
//    if [ "\${#files[@]}" -gt 0 ]; then
//      for f in "\${files[@]}"; do
//        s1=\$(get_size "\$f" 2>/dev/null || echo 0)
//        sleep 2
//        s2=\$(get_size "\$f" 2>/dev/null || echo 0)
//        if [ "\$s1" -gt 0 ] && [ "\$s1" -eq "\$s2" ]; then
//          target="\$f"
//          echo "[slurpy_hic] stable file: \$target (size=\$s2)" >&2
//          break 2
//        fi
//      done
//    fi
//    sleep 1
//  done
//
//  mv aligned/*valid*hic hicfile.hic
//  """
//}

process slurpy_hic {
  tag "${exp}/${ts}"

  input:
  tuple val(exp), val(ts), val(structure), val(reference), val(mitochondria), val(ref_resolution), val(chroms_p1_path)

  output:
  tuple val(exp), val(ts), val(ref_resolution), file('hicfile.hic')

  script:
  """
  set -euo pipefail

  # Symlink structure into work dir (by basename)
  ln -s "${structure}" "./fastq"
  ln -s "/panfs/biopan04/4DGENOMESEQ/EDGE_WORKFLOW/workflow/nextflow/SLURPY/" ./SLURPY
  #WATCH_DIR="aligned"
  #FIXED_REGEX=".*valid.hic"   # <--- your fixed regex; adjust as needed

  echo "[slurpy_hic] exp=${exp} ts=${ts} structure=${structure}" >&2
  echo "[slurpy_hic] genomelist1 PATH: ${chroms_p1_path}" >&2
  #mkdir -p "\$WATCH_DIR"
  export PATH=\$PATH:/panfs/biopan04/4DGENOMESEQ/EDGE_WORKFLOW/workflow/nextflow/SLURPY/
  # --- Run your Python producer here ---
  cmd="/panfs/biopan04/4DGENOMESEQ/EDGE_WORKFLOW/workflow/nextflow/SLURPY/slurm.py -r ${reference} -P fast,tb,gpu --fq ${structure} -G ${chroms_p1_path} -F 150000 15000000 -J /panfs/biopan04/4DGENOMESEQ/EDGE_WORKFLOW/workflow/nextflow/SLURPY/juicer_tools_1.22.01.jar"
    # add mitochondria flag only if present
  if [[ -n "${mitochondria}" ]]; then
    cmd="\$cmd --mtDNA ${mitochondria}"
  fi

  echo "[slurpy_hic] running: \$cmd" >&2
  eval \$cmd |& tee slurpy_hic.log

  # --- Watcher: wait forever for a regex match; ensure size stability (no timeout, no poll controls) ---
    # --- Watcher: wait indefinitely for FIXED_REGEX match, ensure size stability ---
  cat > watcher.sh <<'BASH'
  set -euo pipefail
  watch_dir="aligned"
  regex=".*valid.hic"

  get_size() { stat -c %s "\$1" 2>/dev/null || wc -c < "\$1"; }

  echo "[slurpy_hic] watching dir=\$watch_dir regex=\$regex" >&2
  target=""
  while :; do
    mapfile -t files < <(find "\$watch_dir" -type f -regextype posix-extended -regex "\$regex" -print 2>/dev/null || true)
    if (( \${#files[@]} > 0 )); then
      for f in "\${files[@]}"; do
        s1=\$(get_size "\$f" 2>/dev/null || echo 0)
        sleep 300
        s2=\$(get_size "\$f" 2>/dev/null || echo 0)
        if [[ "\$s1" -gt 0 && "\$s1" -eq "\$s2" ]]; then
          target="\$f"
          echo "[slurpy_hic] stable match: \$target (size=\$s2)" >&2
          break 2
        fi
      done
    fi
    sleep 60
  done

  mv "\$target" hicfile.hic
  BASH

  #export WATCH_DIR FIXED_REGEX
  bash watcher.sh
  """
}

// ---------------- Process 2: hic2struct_one (parallel per chromosome) ----------------
process hic2struct_one {
  tag "${exp}/${ts} [${chr}]"
  errorStrategy 'ignore'
  
  input:
  tuple val(exp), val(ts), val(res), path(stdfile), val(chr)

  output:
  tuple val(exp), val(ts), val(res), val(chr), path("hic_${chr}")

  script:
  """
  set -euo pipefail
  outdir="hic_${chr}"
  echo "Output will be put in \$outdir"
  mkdir -p "\$outdir"
  echo "Processing experiment ${exp}, time step ${ts}, chromosome ${chr}, at resolution of ${res} bp"
  module load openmpi
  # Run hic2structure once per chromosome, using resolution from YAML
  python -m hic2structure --verbose --resolution ${res} --chromosome ${chr} --bond-coeff 55 --count-threshold 10 -o "\$outdir" ${stdfile}
  module unload openmpi
  compgen -G "\$outdir/*" > /dev/null || { echo "[hic2struct_one] No output for chromosome ${chr}" >&2; exit 1; }
  """
}

// ---------------- Process 3: hic_concat (fan-in) ----------------
// Calls your concat_chroms.py with: --indir --chroms --resolution --outdir
process hic_concat {
  tag "${exp}/${ts}"

  publishDir params.outdir, mode: 'copy', saveAs: { produced ->
    "ensemble/experiments/${produced}"
  }

  input:
  tuple val(exp), val(ts), val(res), path(hicdirs)
  val chroms_list

  output:
  file("${exp}/${ts}/structure.csv")

  script:
  """
  set -euo pipefail
  
  out="${exp}/${ts}"
  mkdir -p "\$out"
  
  indir="gather_${ts}"
  mkdir -p "\$indir"
  
  # Build dynamic list only from per-chrom dirs that actually produced structure.csv
  : > chroms.dynamic
  for d in ${hicdirs}; do
    base=\$(basename "\$d")      # e.g. hic_NC_023642.1
    chr="\${base#hic_}"          # -> NC_023642.1
    if [ -s "\$d/structure.csv" ]; then
      mkdir -p "\$indir/\$chr"
      cp -a "\$d"/. "\$indir/\$chr/"
      echo "\$chr" >> chroms.dynamic   # <-- on its own line
    else
      echo "[hic_concat] Skipping chromosome \$chr (no structure.csv)" >&2
    fi
  done
  
  # If nothing succeeded for this timestep, exit cleanly
  if [ ! -s chroms.dynamic ]; then
    echo "[hic_concat] No valid chromosome data for ${exp}/${ts}" >&2
    exit 0
  fi


  # 4) Run concatenation
  python /panfs/biopan04/4DGENOMESEQ/EDGE_WORKFLOW/workflow/nextflow/concat_chroms.py \\
    --indir      "\$indir" \\
    --chroms     chroms.dynamic \\
    --resolution ${res} \\
    --outdir     "\$out"

  # 5) Verify expected output
  test -s "\$out/structure.csv" \
  || { echo "[hic_concat] Expected file not produced: \$out/structure.csv" >&2; exit 1; }
"""
}
//----------------- Process 4: copy input files for data provenance -----------------
/* ---------- process definition (no `from` here) ---------- */
process save_inp {
  tag "save_inp"

  input:
  path yaml_file
  path chroms_file

  output:
  file("input.yaml")
  file("results_autosomes.tsv")

  publishDir params.outdir, mode: 'copy', saveAs: { produced -> "ensemble/provenance/${produced}" }

  script:
  """
  set -euo pipefail
  cp "${yaml_file}"   input.yaml
  cp "${chroms_file}" results_autosomes.tsv
  """
}