#!/usr/bin/env nextflow
nextflow.enable.dsl=2

import groovy.yaml.YamlSlurper

/* ===================== PARAMS (safe defaults) ===================== */
params.yaml         = params.yaml ?: null
params.inp          = params.inp  ?: null          // legacy knob; initialized to suppress WARN
params.outdir       = params.outdir ?: 'results'
params.chroms       = params.chroms ?: null
params.chroms_file  = params.chroms_file ?: null
workflow {
/* ===================== YAML: resolve & parse ===================== */
def yaml_path = (params.yaml ?: params.inp ?: 'input.yaml') as String
def yaml_file = new File(yaml_path).getCanonicalFile()
assert yaml_file.exists() && yaml_file.isFile() : "YAML not found: ${yaml_file}"

def cfg = new YamlSlurper().parse(yaml_file)
assert cfg instanceof Map : "Top-level YAML must be a mapping"
assert cfg?.ensemble?.experiments instanceof List && cfg.ensemble.experiments, "Missing ensemble.experiments (non-empty list)"

/* ===================== Global reference fields ===================== */
def ens            = cfg.ensemble
def ref_block      = ens.reference ?: [:]
def ref_sequence   = (ref_block.sequence     ?: '').toString()
def ref_mito       = (ref_block.mitochondria ?: '').toString()
def ref_resolution = (ref_block.resolution   ?: '').toString()

/* Chromosomes: YAML -> file -> param */
List<String> chroms_included = []
def chroms_raw = ref_block?.chromosomes?.included
if (chroms_raw instanceof List) {
  chroms_included = chroms_raw.collect { it?.toString()?.trim() }.findAll { it }
} else if (chroms_raw != null) {
  chroms_included = chroms_raw.toString().split(/\s*,\s*/).collect { it.trim() }.findAll { it }
}
if (chroms_included.isEmpty() && params.chroms_file) {
  def f = new File(params.chroms_file as String)
  assert f.exists() && f.isFile() : "Chromosome file not found: ${params.chroms_file}"
  chroms_included = f.readLines().collect { it.trim() }.findAll { it }
}
if (chroms_included.isEmpty() && params.chroms) {
  chroms_included = params.chroms.toString().split(/\s*,\s*/).collect { it.trim() }.findAll { it }
}
assert chroms_included, "No chromosomes provided via YAML/--chroms_file/--chroms"

/* ===================== Build EXACT 6-field tuples =====================
 * (exp, ts, structure, ref_sequence, ref_mito, ref_resolution)
 */
def TSTEPS6 = Channel
  .from( ens.experiments )
  .flatMap { exp ->
    def exp_name = exp?.name?.toString()
    assert exp_name : "Experiment missing 'name'"
    def tss = (exp?.timesteps ?: []) as List
    assert tss, "Experiment '${exp_name}' must have timesteps"

    tss.collect { ts ->
      String ts_name
      String structure
      if (ts instanceof Map) {
        ts_name  = ts.name?.toString()
        structure= ts.structure?.toString()
      } else {
        ts_name  = ts?.toString()
        structure= null
      }
      assert ts_name   : "Timestep missing 'name' in experiment '${exp_name}'"
      assert structure : "Timestep '${ts_name}' in '${exp_name}' missing 'structure'"

      tuple(exp_name, ts_name, structure, ref_sequence, ref_mito, ref_resolution)
    }
  }
  // 🔒 Validate shape early (helpful if something upstream is off)
  .map { t ->
    assert t instanceof List && t.size()==6 : "Expected 6 fields, got: ${t}"
    assert !t.any{ it==null } : "Tuple contains nulls: ${t}"
    t
  }

/* ===================== Group by experiment safely =====================
 * Pair the key with the ORIGINAL full tuple, group, then emit the stored tuples.
 * This avoids arity/shape surprises entirely.
 */
def EXP_GROUPED = TSTEPS6
  .map { t -> tuple(t[0], t) }   // (exp, original_6_tuple)
  .groupTuple()                  // -> (exp, [original_6_tuple, ...])

def TSTEPS_FOR_SLURPY = EXP_GROUPED.flatMap { pair ->
  (List) pair[1]                 // emit the list of original 6-tuples as-is
}

/* ===================== Chromosome channels ===================== */
Channel.value(chroms_included).set { CHR_LIST_CH }  // whole list (for concat)
Channel.from(chroms_included).set   { CHR_CH }      // each chromosome (for fan-out)

/* ============================ WORKFLOW ============================ */

  /* 1) Build standard.name; carry resolution forward as the 3rd value */
  def std_ch = slurpy_hic(TSTEPS_FOR_SLURPY)   // -> (exp, ts, ref_resolution, file 'standard.name')

  /* 2) Fan-out per chromosome: (exp, ts, res, stdfile) × chr → (exp, ts, res, stdfile, chr) */
  def per_chr_ch = std_ch.cross(CHR_CH).map { left, chr ->
    def (exp, ts, res, stdfile) = left
    tuple(exp, ts, res, stdfile, chr as String)
  }

  /* 3) One hic2structure per chromosome (parallel) */
  def per_chr_dirs = hic2struct_one(per_chr_ch)  // -> (exp, ts, res, chr, path 'hic_<chr>')

  /* 4) Group by composite key (exp, ts, res) for concat */
  def grouped = per_chr_dirs
    .map { t -> tuple([t[0], t[1], t[2]] as List, t[4]) }   // ([exp,ts,res], hicdir)
    .groupTuple()                                           // -> ([exp,ts,res], [dir1, dir2, ...])
    .map { pair ->                                          // destructure defensively
      def key  = pair[0] as List
      def dirs = pair[1]
      tuple(key[0] as String, key[1] as String, key[2] as String, dirs)   // (exp, ts, res, [dirs])
    }

  /* 5) Final concat and publish */
  hic_concat(grouped, CHR_LIST_CH)
}

// ---------------- Process 1 ----------------
// Creates 'abc/xyz.123', renames to 'standard.name' using values from YAML.
process slurpy_hic {
  tag "${exp}/${ts}"

  input:
  tuple val(exp), val(ts), val(structure), val(reference), val(mitochondria), val(res)

  output:
  // carry forward resolution so hic2structure can use it
  tuple val(exp), val(ts), val(res), file('hicfile.hic')

  script:
  """
  set -euo pipefail

  # Example payload that embeds the inputs into standard.name
  mkdir -p ${exp}/${ts}
  # generate the genome index for bwa aligner
  bwa index ${reference}

  # run SLUR-(M)-PY
  /panfs/biopan04/4DGENOMESEQ/EDGE_WORKFLOW/workflow/nextflow/SLURPY/slurm.py -r ${reference} -P fast,tb,gpu -M ${mitochondria} -q ${structure} -F 150000 15000000 -J /panfs/biopan04/4DGENOMESEQ/EDGE_WORKFLOW/workflow/nextflow/SLURPY/juicer_tools_1.22.01.jar

  mv aligned/*valid*hic hicfile.hic
  """
}

// ---------------- Process 2: hic2struct_one (parallel per chromosome) ----------------
process hic2struct_one {
  tag "${exp}/${ts} [${chr}]"

  input:
  tuple val(exp), val(ts), val(res), path(stdfile), val(chr)

  output:
  tuple val(exp), val(ts), val(res), val(chr), path("hic_${chr}")

  script:
  """
  set -euo pipefail
  outdir="hic_${chr}"
  mkdir -p ${outdir}

  # Run hic2structure once per chromosome, using resolution from YAML
  python -m hic2structure --verbose --resolution \${res} --chromosome "\${chr}" --bond-coeff 55 --count-threshold 10 \\
                        -o ${outdir} \\
                        ${stdfile}

  #compgen -G "\$outdir/*" > /dev/null || { echo "[hic2struct_one] No output for chromosome ${chr}" >&2; exit 1; }
  """
}

// ---------------- Process 3: hic_concat (fan-in) ----------------
// Calls your concat_chroms.py with: --indir --chroms --resolution --outdir
process hic_concat {
  tag "${exp}/${ts}"

  // publish only the final file to: <outdir>/ensemble/<exp>/<ts>/per_chroms.txt
  publishDir params.outdir, mode: 'copy', saveAs: { produced -> "ensemble/${produced}" }

  input:
  tuple val(exp), val(ts), val(res), path(hicdirs)   // grouped list of per-chrom dirs
  val chroms                                                // full chromosome list (single value)

  output:
  file("${exp}/${ts}/per_chroms.txt")

  script:
  """
  set -euo pipefail
  mkdir -p "${exp}/${ts}"

  # Persist chromosomes (order preserved)
  cat > chroms.list <<'EOF'
  ${chroms.join('\n')}
  EOF

  # Gather per-chrom outputs into a single dir expected by the concat script
  indir="gather_${ts}"
  mkdir -p "\$indir"
  for d in ${hicdirs}; do
    base=\$(basename "\$d")
    cp -a "\$d" "\$indir/\$base"
  done

  # Call your concatenation script with required arguments
  python concat_chroms.py \\
    --indir      "\$indir" \\
    --chroms     "chroms.list" \\
    --resolution ${res} \\
    --outdir     "\${exp}/\${ts}"

  test -s "${exp}/${ts}/per_chroms.txt"
  """
}