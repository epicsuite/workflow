nextflow.enable.dsl=2

/************************************************************
 *                   PARAMETERS & IMPORTS
 ************************************************************/
params.yaml   = params.yaml   ?: 'input.yaml'
params.outdir = params.outdir ?: 'results'

import groovy.yaml.YamlSlurper
import java.util.zip.GZIPInputStream

/************************************************************
 *                   HELPERS
 ************************************************************/
def firstNonNull = { Map m, List<String> keys ->
  for (k in keys) { def v = m[k]; if (v != null && v.toString().trim()) return v.toString() }
  return null
}

/* Read first column from TSV/BED(.gz); skip comments & BED headers */
def readFirstColumn = { String path ->
  def f = new File(path)
  assert f.exists() && f.isFile() : "Contig list file not found: ${path}"
  InputStream is = new FileInputStream(f)
  if (f.name.toLowerCase().endsWith('.gz')) is = new GZIPInputStream(is)
  def br = new BufferedReader(new InputStreamReader(is))
  try {
    List<String> out = []
    String line
    while ((line = br.readLine()) != null) {
      line = line.trim()
      if (!line) continue
      if (line.startsWith('#') || line.startsWith('track') || line.startsWith('browser')) continue
      out << line.split('\t', -1)[0].trim()
    }
    return new LinkedHashSet<String>(out).toList() // dedupe, preserve order
  }
  finally { br?.close(); is?.close() }
}

/************************************************************
 *                   PARSE YAML (MAIN)
 ************************************************************/
def yamlFile = new File(params.yaml).getCanonicalFile()
assert yamlFile.exists() && yamlFile.isFile() : "YAML not found: ${yamlFile}"

def cfg = new YamlSlurper().parse(yamlFile)
assert cfg instanceof Map && cfg.ensemble instanceof Map : "Invalid YAML: missing 'ensemble'"

def refBlock    = (cfg.ensemble.reference ?: [:]) as Map
def refSeq      = (refBlock.sequence     ?: '').toString()
def refMito     = (refBlock.mitochondria ?: '')?.toString()   // optional
def refRes      = (refBlock.resolution   ?: '').toString()
def contigsFile = (refBlock.contigs      ?: '').toString()
assert refSeq && refRes && contigsFile : "YAML must provide reference.sequence, reference.resolution, reference.contigs"

List<String> CONTIGS = readFirstColumn(contigsFile)
println "[CONTIGS] Loaded ${CONTIGS.size()} contigs from ${contigsFile}"
assert CONTIGS && !CONTIGS.isEmpty() : "No contigs parsed from ${contigsFile}"

def experiments = cfg.ensemble.experiments
assert experiments instanceof List && !experiments.isEmpty() : "YAML: ensemble.experiments must be a non-empty list"

/* Build rows of: (exp, ts, structure, struct_stage, ref, mito, res) */
List<List> rows = []
experiments.each { e ->
  def expName = firstNonNull(e as Map, ['name','id','exp']) ?: 'UNKNOWN'
  def timesteps = e.timesteps ?: []
  assert timesteps instanceof List && !timesteps.isEmpty() : "Experiment '${expName}' needs a non-empty 'timesteps' list"
  timesteps.eachWithIndex { ts, i ->
    assert ts instanceof Map : "Experiment '${expName}' timestep #${i+1} must be a map"
    def tsName   = firstNonNull(ts as Map, ['name','id','timestep','timepoint','ts'])
    def tsStruct = firstNonNull(ts as Map, ['structure','dir','path','url'])
    def stage    = (ts.struct_stage ?: ts.proc_stage ?: 1).toString()
    assert tsName && tsStruct : "Experiment '${expName}' timestep missing name/structure"
    rows << [expName, tsName, tsStruct, stage, refSeq, refMito, refRes]
  }
}

println "[YAML] Built ${rows.size()} tuples (exp, ts, structure, stage, ref, mito, res)"

/* Global channel with tuple: (exp, ts, structure, struct_stage, reference, mitochondria, resolution) */
TSTEPS7 = Channel
  .from(rows)
  .map { r ->
    def exp  = r[0] as String
    def ts   = r[1] as String
    def dir  = r[2] as String
    def stg  = r[3]?.toString()
    def ref  = r[4] as String
    def mito = r[5]?.toString()
    def res  = r[6]?.toString()
    tuple(exp, ts, dir, stg, ref, mito, res)
  }

/************************************************************
 *                   PROCESSES
 ************************************************************/

/* ---------- BWA index ---------- */
process bwa_index {
  tag 'bwa_index'
  input:
    path ref_fa
  output:
    path 'bwa_index.done', emit: bwa_index_done
  script:
  """
  set -euo pipefail
  ref="${ref_fa}"
  log="bwa_index.log"
  {
    echo "[bwa_index] reference = \$ref"
    if [ -s "\${ref}.amb" ] && [ -s "\${ref}.ann" ] && [ -s "\${ref}.bwt" ] \
       && [ -s "\${ref}.pac" ] && [ -s "\${ref}.sa" ]; then
      echo "[bwa_index] Existing index detected; skipping."
    else
      echo "[bwa_index] Building BWA index..."
      bwa index "\$ref"
      echo "[bwa_index] Build complete."
    fi
  } | tee "\$log"
  echo "OK" > bwa_index.done
  """
}

/* ---------- slurpy_hic (original IO) WITH WATCHER ---------- *
 * input : tuple val(exp), val(ts), val(structure), val(reference), val(mitochondria), val(ref_resolution), val(contigs_file)
 * output: tuple val(exp), val(ts), val(ref_resolution), path('hicfile.hic')
 */
process slurpy_hic {
  tag "${exp}/${ts}"
  maxForks 2
  input:
    tuple val(exp), val(ts), val(structure), val(reference), val(mitochondria), val(ref_resolution), val(contigs_file)
  output:
    tuple val(exp), val(ts), val(ref_resolution), path('hicfile.hic')
  script:
  """
  set -euo pipefail
  echo "[slurpy_hic] waiting for BWA index of ${reference}" >&2
  while :; do
    if [ -s "${reference}.amb" ] && [ -s "${reference}.ann" ] && \
       [ -s "${reference}.bwt" ] && [ -s "${reference}.pac" ] && \
       [ -s "${reference}.sa" ]; then
      echo "[slurpy_hic] BWA index detected" >&2
      break
    fi
    sleep 10
  done
  # Symlink structure into work dir (by basename); -n tolerates existing links
  ln -sfn "${structure}" "./fastq"
  ln -sfn "/panfs/biopan04/4DGENOMESEQ/EDGE_WORKFLOW/workflow/nextflow/SLURPY/" "./SLURPY"

  echo "[slurpy_hic] exp=${exp} ts=${ts} structure=${structure}" >&2
  echo "[slurpy_hic] contigs (genomelist) PATH: ${contigs_file}" >&2

  # Ensure SLURPY scripts are reachable
  export PATH=\$PATH:/panfs/biopan04/4DGENOMESEQ/EDGE_WORKFLOW/workflow/nextflow/SLURPY/
  echo \$PATH
  # Build the command as an array (no eval, no backslash-continues)
  declare -a cmd
  cmd=( "/panfs/biopan04/4DGENOMESEQ/EDGE_WORKFLOW/workflow/nextflow/SLURPY/slurm.py"
        -r "${reference}"
        -P "fast,tb,gpu"
        --fq "${structure}"
        -G "${contigs_file}"
        -F 150000 15000000
        -J "/panfs/biopan04/4DGENOMESEQ/EDGE_WORKFLOW/workflow/nextflow/SLURPY/juicer_tools_1.22.01.jar"
      )

  # add mitochondria flag only if present
  if [[ -n "${mitochondria}" ]]; then
    cmd+=( --mtDNA "${mitochondria}" )
  fi

  echo "[slurpy_hic] running:" "\${cmd[@]}" >&2
  "\${cmd[@]}" |& tee slurpy_hic.log

  # Watcher: wait indefinitely; require size-stable file matching regex, then rename to hicfile.hic
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
  rm -r aligned bedpe merged splits
  BASH

  bash watcher.sh
  """
}

/* ---------- hic2struct_one (original IO) ---------- *
 * input : tuple val(exp), val(ts), val(ref_resolution), path(hicfile), val(contig)
 * output: tuple val(exp), val(ts), val(ref_resolution), val(contig), path("hic_<contig>") optional true
 */
process hic2struct_one {
  tag "${exp}/${ts} [${contig}]"
  errorStrategy { task.attempt < 20 ? 'retry' : 'ignore' }
  maxErrors 240
  maxForks 4
  input:
    tuple val(exp), val(ts), val(ref_resolution), path(hicfile), val(contig)
  output:
    tuple val(exp), val(ts), val(ref_resolution), val(contig), path("hic_${contig}") optional true
  script:
  """
  set -euo pipefail
  outdir="hic_${contig}"
  mkdir -p "\$outdir"
  echo "[hic2struct_one] EXP=${exp} TS=${ts} CONTIG=${contig} RES=${ref_resolution}" >&2

  module load openmpi
  # Run hic2structure once per chromosome, using resolution from YAML
  python -m hic2structure --verbose --resolution ${ref_resolution} --chromosome ${contig} --bond-coeff 55 --count-threshold 10 -o "\$outdir" ${hicfile}
  module unload openmpi
  compgen -G "\$outdir/*" > /dev/null || { echo "[hic2struct_one] No output for chromosome ${contig}" >&2; exit 1; }
  """
}

/* ---------- hic_concat ---------- */
process hic_concat {
  tag "${exp}/${ts}"
  publishDir params.outdir, mode: 'copy', saveAs: { produced -> "ensemble/experiments/${produced}" }
  input:
    tuple val(exp), val(ts), val(ref_resolution), path(hicdirs)
  output:
    path("${exp}/${ts}/structure.csv")
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
    --resolution ${ref_resolution} \\
    --outdir     "\$out"

  # 5) Verify expected output
  test -s "\$out/structure.csv" \
  || { echo "[hic_concat] Expected file not produced: \$out/structure.csv" >&2; exit 1; }
"""
}

/* ---------- provenance ---------- */
process save_provenance {
  tag 'save_provenance'
  publishDir params.outdir, mode: 'copy', saveAs: { produced -> "ensemble/provenance/${produced}" }
  input:
    path yaml_file
    path contig_file
  output:
    path 'input.yaml'
    path 'contigs.tsv'
  script:
  """
  set -euo pipefail
  cp "${yaml_file}" input.yaml
  cp "${contig_file}" contigs.tsv
  """
}

/************************************************************
 *                   SUB-WORKFLOWS
 ************************************************************/

/* ---------- Stage 1: FASTQ -> BWA index -> HIC -> per-contig -> concat ---------- */
workflow STAGE1_FULL {
  take:
    TSTEPS7_ch

  main:
    // 1) Build the BWA index once and gate Stage-1 work on its completion
    def bwa_token = bwa_index(Channel.fromPath(refSeq, checkIfExists:true)).bwa_index_done

    // 2) Select stage==1 timesteps and prepare tuples for slurpy_hic
    def S1 = TSTEPS7_ch
      .filter { exp, ts, structure, stage, ref, mito, res -> stage?.toString() == '1' }
      .map    { exp, ts, structure, stage, ref, mito, res ->
                 tuple(exp, ts, structure, ref, (mito ?: ''), res, contigsFile as String)
              }
      //.cross(bwa_token)        // gate on index completion
      //.map { left, _ -> left } // drop the token, keep only the left tuple

    // 3) Run slurpy_hic -> emits (exp, ts, res, hicfile)
    println("Running SLURPY")

    def std_ch = slurpy_hic(S1)
    // 4) Fan out per contig (defensive destructuring; avoid getAt on broadcasts)
    def per_contig = std_ch.flatMap { t ->
      if (!(t instanceof List) || t.size()!=4) {
        println "[DEBUG] std_ch bad item -> ${t} (${t?.getClass()?.name})"
        return []   // skip malformed items
      }
      def (exp, ts, res, hicfile) = t
      CONTIGS.collect { c -> tuple(exp, ts, res, hicfile, c as String) }
    }
    println("Calculating 3D structure")
    // 5) Per-contig structure -> (exp, ts, res, contig, dir)
    def pc_dirs = hic2struct_one(per_contig)

    // 6) Group per (exp, ts, res) and pass list of dirs to concat
    def grouped = pc_dirs
      .map { exp, ts, res, contig, dir -> tuple([exp, ts, res] as List, dir) }
      .groupTuple()
      .map { key, dirs ->
        def (exp, ts, res) = key
        tuple(exp as String, ts as String, res as String, dirs as List)
      }

    // 7) Concatenate per timestep
    println("Merging structure outputs")
    hic_concat(grouped)
}

/* ---------- Stage 2: start from .hic -> per-contig -> concat ---------- */
workflow STAGE2_HIC {
  take:
    TSTEPS7_ch

  main:
    // 1) Expand .hic files only for stage==2 timesteps — NO channels inside the stream
    def HIC_FILES = TSTEPS7_ch
      .filter { exp, ts, structure, stage, ref, mito, res -> stage?.toString() == '2' }
      .flatMap { exp, ts, structure, stage, ref, mito, res ->
        // List *.hic under 'structure' using plain File APIs to avoid leaking channels
        def sdir  = structure?.toString()
        def base  = sdir ? new File(sdir) : null
        def files = (base?.exists() && base.isDirectory())
                      ? (base.listFiles()?.findAll { it.name.endsWith('.hic') } ?: [])
                      : []
        // Return a List of 4-tuples; flatMap will flatten it
        files.collect { f -> tuple(exp as String, ts as String, (res ?: ''), file(f.absolutePath)) }
      }
    // HIC_FILES now emits: (exp, ts, res, hicfile)

    // 2) Fan out per contig (defensive destructuring)
    def per_contig = HIC_FILES.flatMap { t ->
      if (!(t instanceof List) || t.size()!=4) {
        println "[DEBUG] HIC_FILES bad item -> ${t} (${t?.getClass()?.name})"
        return []
      }
      def (exp, ts, res, hicfile) = t
      CONTIGS.collect { c -> tuple(exp, ts, res, hicfile, c as String) }
    }

    // 3) Per-contig structure -> (exp, ts, res, contig, dir)
    def pc_dirs = hic2struct_one(per_contig)

    // 4) Group per (exp, ts, res) and pass list of dirs to concat
    def grouped = pc_dirs
      .map { exp, ts, res, contig, dir -> tuple([exp, ts, res] as List, dir) }
      .groupTuple()
      .map { key, dirs ->
        def (exp, ts, res) = key
        tuple(exp as String, ts as String, res as String, dirs as List)
      }

    // 5) Concatenate per timestep
    hic_concat(grouped)
}

/************************************************************
 *                   MAIN (auto-routing)
 ************************************************************/
workflow MAIN {
  def hasS1 = rows.any { it[3]?.toString() == '1' }
  def hasS2 = rows.any { it[3]?.toString() == '2' }
  def T7_CH = TSTEPS7

  if (hasS1) STAGE1_FULL(T7_CH)
  if (hasS2) STAGE2_HIC(T7_CH)
  if (!hasS1 && !hasS2) log.warn "No timesteps with struct_stage 1 or 2; nothing to run."

  def Y = Channel.fromPath(params.yaml, checkIfExists:true)
  def C = Channel.fromPath(contigsFile,   checkIfExists:true)
  save_provenance(Y, C)
}

/* ---------- Default entry ---------- */
workflow { MAIN() }
