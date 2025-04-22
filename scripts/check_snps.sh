#!/usr/bin/env bash
# check_snps.sh  —  Query key SNPs and produce a color-coded report (green/yellow/red)
# Usage: ./check_snps.sh <annotated_vcf.rs.vcf.gz> [section]
#
# Available sections:
#   - metabolism
#   - inflammation
#   - cardiovascular
#   - neurotransmitters
#   - antioxidant
#   - vitamin
#   - drug
#   - athletic
#   - longevity
#   - sleep
#   - cognitive
#   - nutrition
#
# Example: ./check_snps.sh allSites.rs.vcf.gz cognitive

set -euo pipefail

# Path to VCF file - will be passed as first argument
if [[ $# -lt 1 ]]; then
  echo "Usage: $0 <vcf_file> [section_name] [--debug]"
  echo "Available sections: all, metabolism, inflammation, cardiovascular, neurotransmitters, antioxidant, vitamin, drug, athletic, longevity, sleep, cognitive, nutrition"
  exit 1
fi

VCF="$1"
REF_GENOME="data/GRCh38.p14.fa" # Path to reference genome

# Process arguments
SECTION="all"
DEBUG=false

shift 1  # Remove VCF file argument

# Process remaining arguments
while [[ $# -gt 0 ]]; do
  case "$1" in
    --debug)
      DEBUG=true
      shift
      ;;
    --*)
      echo "Unknown option: $1"
      exit 1
      ;;
    *)
      SECTION="$1"
      shift
      ;;
  esac
done

# Function to conditionally print debug info
debug_echo() {
  if $DEBUG; then
    echo "$@" >&2
  fi
}

if [[ ! -f $VCF ]]; then
  echo "Usage: $0 <annotated_vcf.vcf.gz> [section]"
  echo "Available sections: metabolism, inflammation, cardiovascular, neurotransmitters, antioxidant, vitamin, drug, athletic, longevity, sleep, cognitive, nutrition"
  exit 1
fi

# ANSI color codes
GREEN="\e[42m\e[30m"   # black text on green background
YELLOW="\e[43m\e[30m"  # black text on yellow background
RED="\e[41m\e[37m"     # white text on red background
NC="\e[0m"

# Function to print section header
print_section_header() {
  local section=$1
  echo ""
  echo "==== $section ===="
  echo ""
}

# Print header
printf "%s\t%s\t%s\t%s\n" "RSID" "Genotype" "Status" "Interpretation"

# Function to fetch genotype with diagnostics
fetch_gt() {
  local chrpos=$1
  local rsid=$2

  # For diagnostic purposes
  debug_echo "Checking $rsid at position $chrpos"

  # First, check if the VCF file contains this SNP by rsID (using exact word match)
  snp_count=$(bcftools view -H "$VCF" | grep -w "$rsid" | wc -l)
  # Clean up the count to ensure it's just a number
  snp_count=$(echo "$snp_count" | tr -d '\n' | tr -d ' ')
  debug_echo "snp_count: $snp_count"

  if [[ "$snp_count" -gt 0 ]]; then
    debug_echo "Found $rsid in VCF by ID"
    # Get the full line to extract both genotype and alleles
    vcf_line=$(bcftools view -H "$VCF" | grep -w "$rsid" | head -1)
    gt=$(echo "$vcf_line" | awk '{print $10}' | cut -d ':' -f 1)
    ref=$(echo "$vcf_line" | awk '{print $4}')
    alt=$(echo "$vcf_line" | awk '{print $5}')

    debug_echo "DEBUG: ref=$ref, alt=$alt, gt=$gt, vcf_line=$vcf_line"

    # If genotype is missing or empty, set it to 0/0
    if [[ -z "$gt" || "$gt" == "./." || "$gt" == "./" ]]; then
      gt="0/0"
      if [[ -n "$ref" ]]; then
        actual_gt="$ref/$ref"
      else
        actual_gt="Unknown/Unknown"
      fi
    else
      # Convert 0/0, 0/1, 1/1 to actual alleles
      if [[ "$gt" == "0/0" ]]; then
        if [[ -n "$ref" ]]; then
          # Convert to uppercase
          ref_upper=$(echo "$ref" | tr '[:lower:]' '[:upper:]')
          actual_gt="$ref_upper/$ref_upper"
        else
          actual_gt="Unknown/Unknown"
        fi
      elif [[ "$gt" == "1/1" ]]; then
        if [[ -n "$alt" ]]; then
          # Convert to uppercase
          alt_upper=$(echo "$alt" | tr '[:lower:]' '[:upper:]')
          actual_gt="$alt_upper/$alt_upper"
        else
          actual_gt="Unknown/Unknown"
        fi
      elif [[ "$gt" == "0/1" || "$gt" == "1/0" ]]; then
        if [[ -n "$ref" && -n "$alt" ]]; then
          # Convert to uppercase
          ref_upper=$(echo "$ref" | tr '[:lower:]' '[:upper:]')
          alt_upper=$(echo "$alt" | tr '[:lower:]' '[:upper:]')
          actual_gt="$ref_upper/$alt_upper"
        else
          actual_gt="Unknown/Unknown"
        fi
      else
        # For multi-allelic or other complex cases
        # Try to parse the genotype if possible
        IFS='/' read -r a1 a2 <<< "$gt"
        if [[ "$a1" == "0" && -n "$ref" ]]; then
          a1="$ref"
        elif [[ "$a1" == "1" && -n "$alt" ]]; then
          a1="$alt"
        fi

        if [[ "$a2" == "0" && -n "$ref" ]]; then
          a2="$ref"
        elif [[ "$a2" == "1" && -n "$alt" ]]; then
          a2="$alt"
        fi

        if [[ -n "$a1" && -n "$a2" ]]; then
          actual_gt="$a1/$a2"
        else
          actual_gt="Unknown/Unknown"
        fi
      fi
    fi

    debug_echo "Genotype: $gt ($actual_gt)"
    echo "$gt|$actual_gt"
    return
  fi

  # Next, try by position
  # Check if we need to convert chr notation
  if [[ "$chrpos" != chr* && "$chrpos" =~ ^[0-9]+: ]]; then
    # Convert notation from "22:123456" to "chr22:123456"
    chr_num=$(echo "$chrpos" | cut -d ':' -f1)
    pos=$(echo "$chrpos" | cut -d ':' -f2)
    chr_pos="chr$chr_num:$pos"
    debug_echo "Trying with chromosome notation: $chr_pos"
    position_count=$(bcftools view -H -r "$chr_pos" "$VCF" 2>/dev/null | wc -l)
  else
    position_count=$(bcftools view -H -r "$chrpos" "$VCF" 2>/dev/null | wc -l)
  fi

  position_count=$(echo "$position_count" | tr -d '\n' | tr -d ' ')
  debug_echo "position_count: $position_count"

  if [[ "$position_count" -gt 0 ]]; then
    # Get the full line to extract both genotype and alleles
    if [[ "$chrpos" != chr* && "$chrpos" =~ ^[0-9]+: ]]; then
      chr_num=$(echo "$chrpos" | cut -d ':' -f1)
      pos=$(echo "$chrpos" | cut -d ':' -f2)
      chr_pos="chr$chr_num:$pos"
      vcf_line=$(bcftools view -H -r "$chr_pos" "$VCF" | head -1)
      debug_echo "Found by position $chr_pos"
    else
      vcf_line=$(bcftools view -H -r "$chrpos" "$VCF" | head -1)
      debug_echo "Found by position $chrpos"
    fi

    gt=$(echo "$vcf_line" | awk '{print $10}' | cut -d ':' -f 1)
    ref=$(echo "$vcf_line" | awk '{print $4}')
    alt=$(echo "$vcf_line" | awk '{print $5}')

    debug_echo "DEBUG: ref=$ref, alt=$alt, gt=$gt, vcf_line=$vcf_line"

    # If genotype is missing or empty, set it to 0/0
    if [[ -z "$gt" || "$gt" == "./." || "$gt" == "./" ]]; then
      gt="0/0"
      if [[ -n "$ref" ]]; then
        actual_gt="$ref/$ref"
      else
        actual_gt="Unknown/Unknown"
      fi
    else
      # Convert 0/0, 0/1, 1/1 to actual alleles
      if [[ "$gt" == "0/0" ]]; then
        if [[ -n "$ref" ]]; then
          # Convert to uppercase
          ref_upper=$(echo "$ref" | tr '[:lower:]' '[:upper:]')
          actual_gt="$ref_upper/$ref_upper"
        else
          actual_gt="Unknown/Unknown"
        fi
      elif [[ "$gt" == "1/1" ]]; then
        if [[ -n "$alt" ]]; then
          # Convert to uppercase
          alt_upper=$(echo "$alt" | tr '[:lower:]' '[:upper:]')
          actual_gt="$alt_upper/$alt_upper"
        else
          actual_gt="Unknown/Unknown"
        fi
      elif [[ "$gt" == "0/1" || "$gt" == "1/0" ]]; then
        if [[ -n "$ref" && -n "$alt" ]]; then
          # Convert to uppercase
          ref_upper=$(echo "$ref" | tr '[:lower:]' '[:upper:]')
          alt_upper=$(echo "$alt" | tr '[:lower:]' '[:upper:]')
          actual_gt="$ref_upper/$alt_upper"
        else
          actual_gt="Unknown/Unknown"
        fi
      else
        # For multi-allelic or other complex cases
        # Try to parse the genotype if possible
        IFS='/' read -r a1 a2 <<< "$gt"
        if [[ "$a1" == "0" && -n "$ref" ]]; then
          a1="$ref"
        elif [[ "$a1" == "1" && -n "$alt" ]]; then
          a1="$alt"
        fi

        if [[ "$a2" == "0" && -n "$ref" ]]; then
          a2="$ref"
        elif [[ "$a2" == "1" && -n "$alt" ]]; then
          a2="$alt"
        fi

        if [[ -n "$a1" && -n "$a2" ]]; then
          actual_gt="$a1/$a2"
        else
          actual_gt="Unknown/Unknown"
        fi
      fi
    fi

    debug_echo "Genotype: $gt ($actual_gt)"
    echo "$gt|$actual_gt"
    return
  fi

  # If we get here, we couldn't find it in the VCF
  # Try to get the reference allele from the reference genome
  debug_echo "SNP $rsid at $chrpos not found in VCF"

  # Extract chromosome and position
  chr=$(echo "$chrpos" | cut -d ':' -f1)
  pos=$(echo "$chrpos" | cut -d ':' -f2)

  # Remove 'chr' prefix if present
  chr=${chr#chr}

  # Map chromosome number to NCBI accession
  case "$chr" in
    1) acc="NC_000001.11" ;;
    2) acc="NC_000002.12" ;;
    3) acc="NC_000003.12" ;;
    4) acc="NC_000004.12" ;;
    5) acc="NC_000005.10" ;;
    6) acc="NC_000006.12" ;;
    7) acc="NC_000007.14" ;;
    8) acc="NC_000008.11" ;;
    9) acc="NC_000009.12" ;;
    10) acc="NC_000010.11" ;;
    11) acc="NC_000011.10" ;;
    12) acc="NC_000012.12" ;;
    13) acc="NC_000013.11" ;;
    14) acc="NC_000014.9" ;;
    15) acc="NC_000015.10" ;;
    16) acc="NC_000016.10" ;;
    17) acc="NC_000017.11" ;;
    18) acc="NC_000018.10" ;;
    19) acc="NC_000019.10" ;;
    20) acc="NC_000020.11" ;;
    21) acc="NC_000021.9" ;;
    22) acc="NC_000022.11" ;;
    X) acc="NC_000023.11" ;;
    Y) acc="NC_000024.10" ;;
    M|MT) acc="NC_012920.1" ;;
    *) acc="" ;;
  esac

  # Query the reference genome using the accession
  if [[ -n "$acc" ]]; then
    ref_base=$(samtools faidx "$REF_GENOME" "$acc:$pos-$pos" 2>/dev/null | tail -n +2)
    # Convert to uppercase
    ref_base=$(echo "$ref_base" | tr '[:lower:]' '[:upper:]')
  else
    ref_base=""
  fi

  if [[ -n "$ref_base" ]]; then
    debug_echo "Reference allele at $chrpos: $ref_base"
    echo "0/0|$ref_base/$ref_base"
  else
    # Fall back to Unknown if we can't get the reference base
    debug_echo "Reference allele not found"
    echo "0/0|Unknown/Unknown"
  fi
}

# Function to categorize and display SNP
check_snp() {
  local snp_data=$1
  IFS='|' read -r rsid chrpos protective_gt risk_gt interpretation <<< "$snp_data"

  # Get actual genotype - pass both chrpos and rsid
  local gt_info=$(fetch_gt "$chrpos" "$rsid")

  # Split into encoded genotype and actual alleles
  IFS='|' read -r gt actual_gt <<< "$gt_info"

  # If genotype is empty, use 0/0
  if [[ -z "$gt" ]]; then
    gt="0/0"
    actual_gt="Unknown/Unknown"
  fi

  # Determine risk category
  local status
  local color
  if [[ "$gt" == "$protective_gt" ]]; then
    status="GOOD"
    color="$GREEN"
  elif [[ "$gt" == "$risk_gt" ]]; then
    status="RISK"
    color="$RED"
  else
    status="MODERATE"
    color="$YELLOW"
  fi

  # Print result
  printf "%s\t%s (%s)\t${color}%s${NC}\t%s\n" \
    "$rsid" "$gt" "$actual_gt" "$status" "$interpretation"
}

# Define separate arrays for each section
# Format: rsID|chr:pos|protective_gt|risk_gt|interpretation

# Metabolism SNPs
metabolism_snps=(
  "rs4680|22:19963748|0/0|1/1|COMT G/G: Better dopamine metabolism; A/A: Anxiety-prone"
  "rs1801133|1:11856378|0/0|1/1|MTHFR C677T: C/C normal folate metabolism; T/T reduced (~40% function)"
  "rs1801131|1:11854476|0/0|1/1|MTHFR A1298C: A/A normal; C/C reduced folate metabolism activity"
  "rs1050152|7:30645223|1/1|0/0|SLC22A4 1672C/T: T/T increased carnitine transport; C/C reduced"
  "rs1799930|8:18258103|0/0|1/1|NAT2 G590A: G/G rapid acetylator; A/A slow acetylator affecting drug metabolism"
  "rs1799931|8:18257854|0/0|1/1|NAT2 G857A: G/G rapid acetylator; A/A slow acetylator affecting drug metabolism"
  "rs1208|8:18258316|0/0|1/1|NAT2: A/A rapid acetylator; G/G slow acetylator affecting drug metabolism"
  "rs1801131|1:11854476|0/0|1/1|MTHFR A1298C: A/A normal; C/C reduced folate metabolism"
  "rs662|7:94937446|1/1|0/0|PON1 Q192R: G/G better ability to detoxify organophosphates; A/A reduced"

)

# Inflammation SNPs
inflammation_snps=(
  "rs20417|1:186722460|0/0|1/1|PTGS2 -765G>C: G/G associated with lower COX2 expression; C/C higher"
  "rs1800795|7:22766645|0/0|1/1|IL6 -174G>C: G/G associated with lower IL-6 inflammation; C/C higher"
  "rs1800629|6:31543059|0/0|1/1|TNF -308G>A: G/G associated with lower TNF-α production; A/A higher"
  "rs16944|2:113594387|0/0|1/1|IL1B -511C>T: C/C associated with lower IL-1β production; T/T higher"
  "rs1143634|2:113590390|0/0|1/1|IL1B +3954C>T: C/C associated with lower IL-1β production; T/T higher"
  "rs4553808|4:123377980|1/1|0/0|CTLA4 -1661A>G: A/A associated with better T-cell function; G/G reduced"
  "rs4580644|11:120154442|1/1|0/0|PAWR: G/G protective against inflammatory disease; A/A increased risk"

)

# Cardiovascular SNPs
cardiovascular_snps=(
  "rs1801282|3:12393125|1/1|0/0|PPARG Pro12Ala: G/G protective against Type 2 diabetes; C/C increased risk"
  "rs5186|3:148459988|0/0|1/1|AGTR1 A1166C: A/A lower risk of hypertension; C/C increased risk"
  "rs1799983|7:150999023|0/0|1/1|NOS3 G894T: G/G better nitric oxide production; T/T reduced production"
  "rs3025058|11:102842889|0/0|1/1|MMP3 -1171 5A/6A: 6A/6A lower MMP3 activity; 5A/5A higher activity"
  "rs1799752|17:64213675|0/0|1/1|ACE Insertion/Deletion: I/I better endurance; D/D better for power/strength"
  "rs9939609|16:53820527|0/0|1/1|FTO: T/T lower obesity risk; A/A higher risk (0.4kg/allele weight gain)"
  "rs7412|19:44908822|1/1|0/0|APOE rs7412: T/T (E2) lower cholesterol, reduced Alzheimer's risk"
  "rs429358|19:44908684|0/0|1/1|APOE rs429358: T/T lower Alzheimer's risk; C/C (E4) higher risk"

)

# Neurotransmitters SNPs
neurotransmitters_snps=(
  "rs6265|11:27679916|1/1|0/0|BDNF G196A: G/G better stress resilience and memory; A/A reduced"
  "rs6313|13:46895829|1/1|0/0|HTR2A T102C: T/T better response to antidepressants; C/C reduced"
  "rs6311|13:46897343|1/1|0/0|HTR2A -1438G>A: A/A higher serotonin receptor expression; G/G lower"
  "rs53576|3:8804371|1/1|0/0|OXTR: G/G better empathy, lower stress; A/A reduced social sensitivity"
  "rs1800497|11:113400106|0/0|1/1|ANKK1/DRD2 Taq1A: C/C higher dopamine receptor density; T/T lower"
  "rs4633|22:19963563|0/0|1/1|COMT: C/C higher COMT enzyme function; T/T lower function"
  "rs2868771|22:40852480|1/1|0/0|COMT: G/G higher frontal cognitive performance; T/T lower"

)

# Antioxidant SNPs
antioxidant_snps=(
  "rs1050450|3:49394834|0/0|1/1|GPX1 C599T: C/C better glutathione peroxidase function; T/T reduced"
  "rs1695|11:67585218|0/0|1/1|GSTP1 A313G: A/A better detoxification capacity; G/G reduced"
  "rs4880|6:159692840|0/0|1/1|SOD2 T175C: C/C better mitochondrial transport of SOD2; T/T reduced"
  "rs1001179|11:34460231|0/0|1/1|CAT -262C>T: C/C higher catalase activity; T/T lower activity"

)

# Vitamin Metabolism SNPs
vitamin_snps=(
  "rs10741657|11:14914878|1/1|0/0|CYP2R1: A/A higher Vitamin D levels; G/G lower levels"
  "rs7041|4:71750730|0/0|1/1|GC: G/G higher Vitamin D binding protein; T/T lower levels"
  "rs2228570|12:47879112|0/0|1/1|VDR FokI: C/C more active Vitamin D receptor; T/T less active"
  "rs1544410|12:47846052|1/1|0/0|VDR BsmI: A/A better bone mineral density; G/G reduced"
  "rs731236|12:47844610|1/1|0/0|VDR TaqI: T/T better calcium metabolism; C/C reduced"

)

# Drug Metabolism SNPs
drug_snps=(
  "rs12248560|10:94761900|1/1|0/0|CYP2C19*17: T/T ultra-rapid metabolizer of many drugs; C/C normal"
  "rs4244285|10:94781859|0/0|1/1|CYP2C19*2: G/G normal function; A/A reduced metabolism"
  "rs4986893|10:94780653|0/0|1/1|CYP2C19*3: G/G normal function; A/A reduced metabolism"
  "rs28399504|10:94781858|0/0|1/1|CYP2C19*4: A/A normal function; G/G reduced metabolism"
  "rs1799853|10:94942290|0/0|1/1|CYP2C9*2: C/C normal metabolizer of warfarin/NSAIDs; T/T reduced"
  "rs1057910|10:94981296|0/0|1/1|CYP2C9*3: A/A normal metabolizer of warfarin/NSAIDs; C/C reduced"
  "rs3892097|22:42128945|0/0|1/1|CYP2D6*4: G/G normal metabolizer of ~25% of drugs; A/A poor"
  "rs35742686|22:42130692|0/0|1/1|CYP2D6*3: T/T normal metabolizer; - (deletion) poor metabolizer"
  "rs5030655|22:42130655|0/0|1/1|CYP2D6*6: T/T normal metabolizer; - (deletion) poor metabolizer"
  "rs59421388|7:99652770|0/0|1/1|CYP3A4*3: T/T normal metabolizer of most drugs; C/C reduced"
  "rs776746|7:99652770|1/1|0/0|CYP3A5*3: A/A normal expression of CYP3A5; G/G reduced expression"
  "rs1045642|7:87138645|1/1|0/0|ABCB1 3435C>T: T/T lower P-glycoprotein activity; C/C higher"
  "rs2032582|7:87160618|0/0|1/1|ABCB1 2677G>T/A: G/G normal P-glycoprotein; T/T or A/A altered"
  "rs762551|15:74749576|1/1|0/0|CYP1A2*1F: A/A rapid caffeine metabolism; C/C slower metabolism"
  "rs1799931|8:18257854|0/0|1/1|NAT2*7: G/G faster drug acetylator; A/A slower acetylator"

)

# Athletic Performance SNPs
athletic_snps=(
  "rs1815739|11:66560624|1/1|0/0|ACTN3 R577X: C/C better for power/sprint performance; T/T endurance"
  "rs8192678|4:89015940|0/0|1/1|PPARGC1A G482S: G/G better aerobic capacity; A/A reduced"
  "rs9939609|16:53820527|0/0|1/1|FTO: T/T less obesity risk, better exercise response; A/A higher risk"
  "rs1042713|5:148826877|1/1|0/0|ADRB2 G46A: G/G better lung function response to albuterol; A/A reduced"
  "rs1042714|5:148826910|0/0|1/1|ADRB2 C79G: C/C protective against obesity; G/G increased risk"

)

# Longevity SNPs
longevity_snps=(
  "rs2802292|6:108587315|1/1|0/0|FOXO3: G/G associated with longevity; T/T normal aging"
  "rs1800795|7:22766645|0/0|1/1|IL6 -174G>C: G/G lower inflammation levels; C/C higher"
  "rs1061170|1:196659237|0/0|1/1|CFH Y402H: T/T lower macular degeneration risk; C/C higher risk"
  "rs1799945|6:26091179|0/0|1/1|HFE H63D: C/C normal iron metabolism; G/G increased iron absorption"
  "rs1800562|6:26093141|0/0|1/1|HFE C282Y: G/G normal iron metabolism; A/A hemochromatosis risk"
  "rs2070744|7:150992991|1/1|0/0|NOS3 -786T>C: T/T better nitric oxide production; C/C reduced"
  "rs1800629|6:31543031|0/0|1/1|TNF -308G>A: G/G lower inflammation markers; A/A higher"

)

# Sleep SNPs
sleep_snps=(
  "rs73598374|12:121348906|0/0|1/1|ADA G22A: G/G deeper sleep patterns; A/A lighter sleep"
  "rs12649507|4:183636156|0/0|1/1|CLOCK: A/A better sleep quality; G/G poorer sleep quality"
  "rs1801260|4:183635768|0/0|1/1|CLOCK 3111T>C: T/T better circadian rhythm; C/C delayed sleep phase"
  "rs10830963|11:92708710|0/0|1/1|MTNR1B: C/C better glycemic control; G/G increased fasting glucose"

)

# Cognitive Function SNPs
cognitive_snps=(
  "rs429358|19:44908684|0/0|1/1|APOE epsilon 4: T/T lower Alzheimer's risk; C/C higher risk"
  "rs7412|19:44908822|1/1|0/0|APOE epsilon 2: T/T lower Alzheimer's risk; C/C normal risk"
  "rs6265|11:27679916|1/1|0/0|BDNF G196A: G/G better cognitive function; A/A reduced"
  "rs1800497|11:113400106|0/0|1/1|ANKK1 Taq1A: C/C better dopamine signaling; T/T reduced"
  "rs4680|22:19963748|0/0|1/1|COMT G472A: G/G better working memory under stress; A/A reduced"

)

# Nutrition SNPs
nutrition_snps=(
  "rs9939609|16:53820527|0/0|1/1|FTO: T/T less obesity risk, better response to diet"
  "rs1801282|3:12393125|1/1|0/0|PPARG Pro12Ala: Ala/Ala better insulin sensitivity"
  "rs5082|1:161241930|1/1|0/0|APOA2: C/C sensitive to saturated fat intake"
  "rs4994|8:37823798|0/0|1/1|ADRB3 Trp64Arg: T/T better response to exercise"
  "rs1799883|4:120241902|0/0|1/1|FABP2 Ala54Thr: G/G better fat metabolism"
  "rs1051730|15:78601997|0/0|1/1|CHRNA3: G/G associated with lower nicotine dependence"
  "rs4988235|2:135851076|1/1|0/0|MCM6/LCT: T/T associated with lactose tolerance"
  "rs713598|7:141673345|1/1|0/0|TAS2R38 PAV: G/G bitter taste perception, vegetable intake"
  "rs1726866|7:141972804|0/0|1/1|TAS2R38 AVI: A/A better taste receptor function"
  "rs10246939|7:141972933|1/1|0/0|TAS2R38: C/C better taste receptor function"
)

# Process the requested section
if [[ "$SECTION" == "all" || "$SECTION" == "metabolism" ]]; then
  print_section_header "Metabolism and Detoxification"
  for snp in "${metabolism_snps[@]}"; do
    check_snp "$snp"
  done
fi

if [[ "$SECTION" == "all" || "$SECTION" == "inflammation" ]]; then
  print_section_header "Inflammation and Immunity"
  for snp in "${inflammation_snps[@]}"; do
    check_snp "$snp"
  done
fi

if [[ "$SECTION" == "all" || "$SECTION" == "cardiovascular" ]]; then
  print_section_header "Cardiovascular Health"
  for snp in "${cardiovascular_snps[@]}"; do
    check_snp "$snp"
  done
fi

if [[ "$SECTION" == "all" || "$SECTION" == "neurotransmitters" ]]; then
  print_section_header "Neurotransmitters and Mental Health"
  for snp in "${neurotransmitters_snps[@]}"; do
    check_snp "$snp"
  done
fi

if [[ "$SECTION" == "all" || "$SECTION" == "antioxidant" ]]; then
  print_section_header "Antioxidant Capacity"
  for snp in "${antioxidant_snps[@]}"; do
    check_snp "$snp"
  done
fi

if [[ "$SECTION" == "all" || "$SECTION" == "vitamin" ]]; then
  print_section_header "Vitamin Metabolism"
  for snp in "${vitamin_snps[@]}"; do
    check_snp "$snp"
  done
fi

if [[ "$SECTION" == "all" || "$SECTION" == "drug" ]]; then
  print_section_header "Drug Metabolism"
  for snp in "${drug_snps[@]}"; do
    check_snp "$snp"
  done
fi

if [[ "$SECTION" == "all" || "$SECTION" == "athletic" ]]; then
  print_section_header "Athletic Performance"
  for snp in "${athletic_snps[@]}"; do
    check_snp "$snp"
  done
fi

if [[ "$SECTION" == "all" || "$SECTION" == "longevity" ]]; then
  print_section_header "Longevity and Healthy Aging"
  for snp in "${longevity_snps[@]}"; do
    check_snp "$snp"
  done
fi

if [[ "$SECTION" == "all" || "$SECTION" == "sleep" ]]; then
  print_section_header "Sleep and Circadian Rhythm"
  for snp in "${sleep_snps[@]}"; do
    check_snp "$snp"
  done
fi

if [[ "$SECTION" == "all" || "$SECTION" == "cognitive" ]]; then
  print_section_header "Cognitive Function and Neurodegeneration"
  for snp in "${cognitive_snps[@]}"; do
    check_snp "$snp"
  done
fi

if [[ "$SECTION" == "all" || "$SECTION" == "nutrition" ]]; then
  print_section_header "Nutrition and Diet Response"
  for snp in "${nutrition_snps[@]}"; do
    check_snp "$snp"
  done
fi

echo -e "\nNote: This report is for informational purposes only and should not be used for medical diagnosis or treatment decisions."