#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --time=00:15:00
#SBATCH -p hsph
#SBATCH --mem=4G
#SBATCH --job-name=diversity
#SBATCH --output=logs/diversity_%A_%a.out
#SBATCH --error=logs/diversity_%A_%a.err
#SBATCH --array=0-0   # fallback only; pass --array at submit ###** UPDATE **###

# Windowed nucleotide diversity (pi) and Shannon entropy (H) per sample,
# following Hill et al. (Science, aed6094) supplement. Their equation numbers
# are carried through as comments at each computation:
#   Eq. 1  pi_b       per-base nucleotide diversity
#   Eq. 2  pi_window  windowed mean of pi_b
#   Eq. 3  pi_ww      genome-wide mean across windows
#   Eq. 5  H_b        per-base Shannon entropy
#   Eq. 6  H_window   windowed mean of H_b
#   Eq. 7  H_ww       genome-wide mean across windows
# Eq. 4 and Eq. 8 (population-weighted spatial aggregation) are downstream R
# steps, not computed here.

# ids.txt lives in the batch working directory, so one script serves both
# batches without editing -- see informatics/REPROCESSING.md
SAMPLES="${SAMPLES:-ids.txt}"
sample=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$SAMPLES")

# guard: array range longer than the sample list, or --array not passed
if [ -z "$sample" ]; then
    echo "Error: no sample at index ${SLURM_ARRAY_TASK_ID} in $SAMPLES"
    exit 1
fi

# define input paths (output path is set below, once the thresholds are known)
variants_file="freyja_variants/${sample}_variants.tsv"
depth_file="freyja_variants/${sample}_depths.tsv"

# window geometry: Hill et al. use 1000 bp windows starting every 100 bp
WIN=${WIN:-1000}
STEP=${STEP:-100}

# allele-frequency floor, applied here rather than inherited from whatever
# threshold Freyja handed ivar, so the detection limit is identical across
# samples instead of drifting with depth
AF_FLOOR=${AF_FLOOR:-0.03}

# per-site depth defining a "callable" position; used only by the callable
# denominator below. Must be high enough that AF_FLOOR is actually resolvable
# (at 100x a 3% allele is 3 reads).
DEPTH_FLOOR=${DEPTH_FLOOR:-100}

# Output goes in a threshold-named subdirectory (e.g. diversity/AF03DF200/) so
# parameter sets sit side by side instead of overwriting each other -- rerunning
# at a second threshold is cheap and expected. AF tag drops the leading "0.":
# 0.03 -> 03, 0.05 -> 05, 0.1 -> 1. Override OUT_DIR to place it elsewhere.
# Note WIN/STEP are NOT in the tag; vary those and the runs collide. All four
# parameters are echoed into every output row, so a mismatch is detectable.
af_tag=$(printf '%s' "$AF_FLOOR" | sed 's/^0*//; s/^\.//; s/\.//g')
OUT_DIR=${OUT_DIR:-diversity/AF${af_tag}DF${DEPTH_FLOOR}}
outfile="$OUT_DIR/${sample}.diversity.tsv"

# TWO DENOMINATORS. Hill et al. divide each window by the fixed window length
# L, so an uncovered position contributes 0 exactly like a covered-but-
# invariant one. That makes pi scale roughly with coverage breadth, which
# matters here because breadth steps ~2x at the 2023-05-05 primer change.
# Both are emitted so the difference is measurable rather than assumed:
#   *_fixed    -- Hill's definition, denominator = window length
#   *_callable -- denominator = positions at >= DEPTH_FLOOR, with the
#                 numerator restricted to those same positions
# Whole-genome values are the primary output; regions come free once the
# per-base vector exists.

# regions from Hill et al. table S11 (name:start:end, 1-based inclusive)
REGIONS="nsp5:10063:10954,nsp6:10972:11842,methyltr2:20661:21549,s1_ntd:21598:22474,s1_rbd:22516:23185,spike:21563:25384"

# make output and log directories if they don't exist
mkdir -p "$OUT_DIR" logs

echo "Running diversity for $sample (AF_FLOOR=$AF_FLOOR DEPTH_FLOOR=$DEPTH_FLOOR) -> $OUT_DIR"

# check that both inputs exist (04 writes them as a pair)
for f in "$variants_file" "$depth_file"; do
    if [ ! -f "$f" ]; then
        echo "Error: input file not found: $f"
        exit 1
    fi
done

# FILENAME dispatch rather than the NR==FNR idiom: a zero-byte variants file
# would make NR==FNR silently true for the depths file and corrupt the run.
awk -v vfile="$variants_file" -v dfile="$depth_file" \
    -v af="$AF_FLOOR" -v dfloor="$DEPTH_FLOOR" -v win="$WIN" -v step="$STEP" \
    -v sample="$sample" -v regions="$REGIONS" '
function fmt(v) { return (v == "NA" ? "NA" : sprintf("%.8f", v)) }

# ---- pass 1: variants (ivar format, header row) ----------------------------
# columns: 1 REGION  2 POS  3 REF  4 ALT  ...  11 ALT_FREQ  12 TOTAL_DP
FILENAME == vfile {
    if ($2 !~ /^[0-9]+$/) next                   # header or malformed row
    if ($4 !~ /^[ACGTacgt]$/) next               # SNVs only; ivar indels are +A / -A
    f = $11 + 0
    if (f < af) next                             # below the shared AF floor
    p = $2 + 0
    altsum[p] += f                               # total ALT frequency at p
    altsq[p]  += f * f                           # sum of squared ALT frequencies
    alth[p]   += -f * log(f)                     # ALT terms of Shannon H
    if ($12 + 0 > dp[p]) dp[p] = $12 + 0         # reads spanning p
    next
}

# ---- pass 2: depths (CHROM POS REF DEPTH, no header) -----------------------
# Every genome position appears here, so this pass defines the coordinate
# space; the variants file only carries positions with a call.
FILENAME == dfile {
    p = $2 + 0
    if (p > glen) glen = p                    # genome length, from the depths file
    call[p] = ($4 + 0 >= dfloor)

    n = dp[p]
    if (n > 1 && altsum[p] > 0) {
        fref = 1 - altsum[p]
        if (fref < 0) fref = 0                   # multiallelic rounding guard
        # Eq. 1  pi_b = (n/(n-1)) * (1 - sum_i f_i^2), i over ALTs plus REF
        pib[p] = (n / (n - 1)) * (1 - (fref * fref + altsq[p]))
        # Eq. 5  H_b = -sum_i f_i * ln(f_i), same allele set, natural log
        hb[p] = alth[p] + (fref > 0 ? -fref * log(fref) : 0)
        nvar++
    }
    # positions with no call keep pib/hb unset, which awk reads as 0 --
    # "no variation" per Eq. 2, and also how an uncovered position enters
    # the fixed denominator
    if (call[p]) ncall++
}

END {
    # ---- whole genome: sliding windows -----------------------------------
    for (w = 1; w + win - 1 <= glen; w += step) {
        spi = 0; sh = 0; spic = 0; shc = 0; nc = 0
        for (b = w; b < w + win; b++) {
            spi += pib[b]; sh += hb[b]
            if (call[b]) { spic += pib[b]; shc += hb[b]; nc++ }
        }
        # Eq. 2  pi_window = (1/L) * sum_{b in window} pi_b
        # Eq. 6  H_window  = (1/L) * sum_{b in window} H_b
        # L in the paper is the denominator: win for fixed, nc for callable.
        # (The awk variable glen is genome length, a different quantity.)
        pifix += spi / win; hfix += sh / win; nwf++
        if (nc > 0) { pical += spic / nc; hcal += shc / nc; nwc++ }
    }
    # Eq. 3  pi_ww = (1/W) * sum_w pi_window,w  -- across windows, not bases
    # Eq. 7  H_ww  = (1/W) * sum_w H_window,w
    pi_ww_fixed = (nwf > 0 ? pifix / nwf : "NA")
    h_ww_fixed  = (nwf > 0 ? hfix  / nwf : "NA")
    pi_ww_call  = (nwc > 0 ? pical / nwc : "NA")
    h_ww_call   = (nwc > 0 ? hcal  / nwc : "NA")

    # ---- regions: direct means, not windowed ------------------------------
    # Eq. 2 / Eq. 6 with the region as the single window, and no Eq. 3 / Eq. 7
    # step: NSP5 (892 bp) and the RBD (670 bp) are shorter than one window, so
    # a sliding average is not defined for them.
    nr = split(regions, R, ",")
    for (i = 1; i <= nr; i++) {
        split(R[i], G, ":")
        rname[i] = G[1]; rs = G[2] + 0; re = G[3] + 0
        spi = 0; sh = 0; spic = 0; shc = 0; nc = 0
        for (b = rs; b <= re; b++) {
            spi += pib[b]; sh += hb[b]
            if (call[b]) { spic += pib[b]; shc += hb[b]; nc++ }
        }
        rlen = re - rs + 1
        rpif[i] = spi / rlen;  rhf[i] = sh / rlen
        rpic[i] = (nc > 0 ? spic / nc : "NA")
        rhc[i]  = (nc > 0 ? shc  / nc : "NA")
        rnc[i]  = nc
    }

    # ---- emit one header + one row, matching 04.1_coverage_qc.sh ----------
    h = "sample\tgenome_sites\tn_var_sites\tn_callable\taf_floor\tdepth_floor"
    h = h "\twin_size\twin_step\tn_windows\tn_windows_callable"
    h = h "\tpi_ww_fixed\tpi_ww_callable\th_ww_fixed\th_ww_callable"
    for (i = 1; i <= nr; i++)
        h = h "\tpi_" rname[i] "_fixed\tpi_" rname[i] "_callable" \
              "\th_" rname[i] "_fixed\th_" rname[i] "_callable" \
              "\tn_callable_" rname[i]
    print h

    r = sample "\t" glen+0 "\t" nvar+0 "\t" ncall+0 "\t" af "\t" dfloor
    r = r "\t" win "\t" step "\t" nwf+0 "\t" nwc+0
    r = r "\t" fmt(pi_ww_fixed) "\t" fmt(pi_ww_call)
    r = r "\t" fmt(h_ww_fixed)  "\t" fmt(h_ww_call)
    for (i = 1; i <= nr; i++)
        r = r "\t" fmt(rpif[i]) "\t" fmt(rpic[i]) "\t" fmt(rhf[i]) "\t" fmt(rhc[i]) "\t" rnc[i]+0
    print r
}' "$variants_file" "$depth_file" > "$outfile"

echo "Done with $sample: $outfile"

# ---- Running ----------------------------------------------------------------
#
# Once per batch, from any directory -- --chdir selects the batch:
#
#   STORE=/n/holylfs05/LABS/hanage_lab/Lab/hsphfs1/bschaeffer
#   REPO=$STORE/bphc_ww_repo
#   D=$STORE/bphc_ww/baseload_batch2_outputs
#
#   sbatch --chdir=$D --export=ALL,AF_FLOOR=0.03,DEPTH_FLOOR=200 \
#          --array=0-$(($(wc -l < $D/ids.txt)-1)) \
#          $REPO/informatics/scripts/04.2_diversity.sh