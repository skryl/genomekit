#!/usr/bin/env bash
#
# lookup_snps.sh  —  colour‑coded SNP lookup using GNU parallel
# --------------------------------------------------------------------
# Usage:
#   ./lookup_snps.sh                 # all sections in parallel
#   ./lookup_snps.sh cardio          # single section
#   ./lookup_snps.sh -f my.txt sleep # different report + only “sleep”
#   ./lookup_snps.sh -l              # list section names
# --------------------------------------------------------------------

#########################  SNP CATALOGS  ####################################

# Define separate arrays for each section
# Format: rsID|chr:pos|protective_gt|risk_gt|interpretation

# ───────────────────── Metabolism ─────────────────────
# metabolism_snps=(
#   "rs1801133|1:11856378|C/C|T/T|MTHFR C677T: C/C normal folate; T/T ~40 % activity"
#   "rs1801131|1:11854476|A/A|C/C|MTHFR A1298C: A/A normal; C/C reduced activity"
#   "rs1050152|7:30645223|T/T|C/C|SLC22A4 1672C/T: T/T ↑carnitine transport; C/C ↓"
#   "rs1799930|8:18258103|G/G|A/A|NAT2 G590A: G/G rapid; A/A slow acetylator"
#   "rs1799931|8:18257854|G/G|A/A|NAT2 G857A: G/G rapid; A/A slow acetylator"
#   "rs662|7:94937446|G/G|A/A|PON1 Q192R: G/G better organophosphate detox; A/A reduced"

#   # Additional Methylation / Endocrine Variants
#   "rs1805087|1:236729215|A/A|G/G|MTR A2756G: A/A normal function; G/G higher homocysteine"
#   "rs1801394|5:7875012|A/A|G/G|MTRR A66G: A/A normal function; G/G reduced B12 regeneration"
#   "rs3733890|5:139541247|G/G|A/A|BHMT c.742G>A: G/G normal betaine remethylation; A/A lower"
#   "rs234706|21:44489385|C/C|T/T|CBS c.699C>T: C/C normal transsulfuration; T/T upregulated"
#   "rs743572|10:104513936|T/T|C/C|CYP17A1 −34T>C: T/T normal androgen synthesis; C/C ↑androgens"
#   "rs1799941|17:7632191|A/A|G/G|SHBG c.-68G>A: A/A higher SHBG; G/G lower"
# )

# # ─────────────────── Inflammation ───────────────────
# inflammation_snps=(
#   "rs20417|1:186722460|G/G|C/C|PTGS2 −765 G>C: G/G lower COX‑2; C/C higher"
#   "rs1800795|7:22766645|G/G|C/C|IL6 −174 G>C: G/G lower IL‑6; C/C higher"
#   "rs1800629|6:31543059|G/G|A/A|TNF −308 G>A: G/G lower TNF‑α; A/A higher"
#   "rs16944|2:113594387|C/C|T/T|IL1B −511 C>T: C/C lower IL‑1β; T/T higher"
#   "rs1143634|2:113590390|C/C|T/T|IL1B +3954 C>T: C/C lower IL‑1β; T/T higher"
#   "rs4553808|4:123377980|A/A|G/G|CTLA4 −1661 A>G: A/A better T‑cell control; G/G reduced"
#   "rs4580644|11:120154442|A/A|G/G|PAWR: A/A protective (↓NF‑κB); G/G higher inflam‑risk"

#   # Additional Inflammation Variants
#   "rs1800896|1:206945946|G/G|A/A|IL10 −1082 G>A: G/G higher IL‑10; A/A lower"
#   "rs1205|1:159682037|C/C|T/T|CRP 3' UTR variant: C/C lower CRP; T/T higher"
#   "rs4986790|9:117713086|C/C|T/T|TLR4 Asp299Gly: C/C normal LPS response; T/T altered"
# )

# # ─────────────────── Cardiovascular ───────────────────
# cardiovascular_snps=(
#   "rs1801282|3:12393125|G/G|C/C|PPARG Pro12Ala: G/G protective; C/C ↑T2D"
#   "rs429358|19:45411941|T/T|C/C|APOE ε4 marker; C allele ↑Alzheimer/lipids"
#   "rs7412|19:45412079|T/T|C/C|APOE ε2 marker; T allele ↓cholesterol"
#   "rs1799983|7:150696111|G/G|T/T|NOS3 Glu298Asp: G/G normal NO; T/T reduced"
#   "rs5186|3:148459988|A/A|C/C|AGTR1 A1166C: A/A benign; C/C ↑hypertension"
#   "rs1042713|5:148826877|G/G|A/A|ADRB2 Gly16Arg: G/G better β‑agonist response"
#   # "rs1799752|17:61566095|I/I|D/D|ACE I/D: I/I lower ACE; D/D higher – indel ‑ handle separately"

#   # Additional Cardiovascular Variants
#   "rs11591147|1:55505621|C/C|T/T|PCSK9 R46L: T allele ↓LDL; protective"
#   "rs2228671|19:11211158|C/C|T/T|LDLR c.1773C>T: T/T decreased LDL clearance"
#   "rs6025|1:169519049|G/G|A/A|Factor V Leiden R506Q: G/G normal clotting; A/A risk"
#   "rs1799963|11:4672404|G/G|A/A|Prothrombin G20210A: G/G normal clotting; A/A risk"
# )

# # ───────────────── Neurotransmitters ───────────────────
# neurotransmitters_snps=(
#   "rs6280|11:113283688|G/G|A/A|DRD3 Ser9Gly: G/G stronger binding; A/A weaker"
#   "rs4680|22:19963748|G/G|A/A|COMT Val158Met: G/G fast DA clearance; A/A slow"
#   "rs6265|11:27679916|G/G|A/A|BDNF Val66Met: G/G higher plasticity; A/A lower"
#   "rs1800497|11:113400106|G/G|A/A|ANKK1 Taq1A: G/G higher D2 density"
#   "rs1799971|6:154360797|A/A|G/G|OPRM1 A118G: A/A normal opioids; G/G altered"
#   "rs4633|22:19962712|C/C|T/T|COMT rs4633: C/C low DA; T/T higher"
# )

# # ───────────────────── Antioxidant ─────────────────────
# antioxidant_snps=(
#   "rs1050450|3:49394834|C/C|T/T|GPX1 Pro198Leu: C/C optimal GPX1; T/T reduced"
#   "rs1695|11:67352689|A/A|G/G|GSTP1 Ile105Val: A/A normal detox; G/G altered"
#   "rs4880|6:160113872|C/C|T/T|SOD2 Ala16Val: C/C efficient Mn‑SOD import; T/T reduced"
#   "rs1800668|3:49357401|A/A|G/G|GPX1 −602 A>G: A/A higher GPX1; G/G lower"
# )

# # ───────────────────── Vitamin D / Folate ──────────────
# vitamin_snps=(
#   "rs731236|12:48238757|T/T|C/C|VDR TaqI: T/T ↑VDR sensitivity; C/C lower"
#   "rs2228570|12:48272895|G/G|A/A|VDR FokI: G/G active VDR; A/A less active"
#   "rs7041|4:72618334|G/G|T/T|GC/DBP: G/G better D transport; T/T reduced"
#   "rs1544410|12:48239835|A/A|G/G|VDR BsmI: A/A higher receptor activity"
# )

# # ───────────────── Drug Metabolism ─────────────────────
# drug_snps=(
#   "rs12248560|10:94761900|T/T|C/C|CYP2C19*17: T/T ultra‑rapid; C/C normal"
#   "rs4244285|10:94781859|G/G|A/A|CYP2C19*2: G/G normal; A/A loss‑of‑fx"
#   "rs4986893|10:94780653|G/G|A/A|CYP2C19*3: G/G normal; A/A loss‑of‑fx"
#   "rs28399504|10:94781858|A/A|G/G|CYP2C19*4: A/A normal; G/G loss‑of‑fx"
#   "rs1057910|10:96741053|A/A|C/C|CYP2C9*3: A/A normal; C/C ↓activity"
#   "rs1799853|10:96702047|C/C|T/T|CYP2C9*2: C/C normal; T/T ↓activity"
#   "rs3892097|22:42524947|G/G|A/A|CYP2D6*4: G/G normal; A/A poor metaboliser"
#   "rs762551|15:75041917|A/A|C/C|CYP1A2*1F: A/A rapid caffeine clearance"
#   "rs776746|7:99361466|A/A|G/G|CYP3A5*1 expressor vs non‑expresser"
# )

# # ─────────────────── Athletic / Muscle ─────────────────
# athletic_snps=(
#   "rs1815739|11:66560624|C/C|T/T|ACTN3 R577X: C/C power; T/T endurance"
#   "rs8192678|4:89011240|G/G|A/A|PPARGC1A Gly482Ser: G/G endurance ↑"
#   "rs1042713|5:148826877|G/G|A/A|ADRB2 Gly16Arg: G/G better training response"
#   "rs4994|8:37823798|T/T|C/C|ADRB3 Trp64Arg: T/T normal thermogenesis"
#   "rs6265|11:27679916|G/G|A/A|BDNF Val66Met: G/G quicker recovery"
# )

# # ───────────────────── Longevity ───────────────────────
# longevity_snps=(
#   "rs2802292|6:108587315|G/G|T/T|FOXO3 G allele linked to longevity"
#   "rs1042522|17:7579472|G/G|C/C|TP53 Pro72Arg: G/G efficient apoptosis"
#   "rs929347|3:57126395|A/A|G/G|HESX1 A/A longer lifespan (pop‑specific)"
#   "rs10811661|9:22134094|C/C|T/T|CDKN2A/B T/T ↑T2D risk"
#   "rs4880|6:160113872|C/C|T/T|SOD2 Ala16Val duplicated for anti‑oxidant"
# )

# # ───────────────────── Sleep / Chrono ──────────────────
# sleep_snps=(
#   "rs73598374|6:90503222|G/G|A/A|ADA G22A: G/G normal adenosine; A deeper sleep"
#   "rs12413112|10:13124039|G/G|A/A|CTNNA3 G allele better sleep quality"
#   "rs324981|2:27494431|A/A|T/T|NPSR1 Asn107Ile: T/T longer sleep latency"
#   "rs10830963|11:92711137|C/C|G/G|MTNR1B C/C normal glucose chronotype"

#   # Additional Circadian Variants
#   "rs1801260|4:56158451|T/T|C/C|CLOCK 3111T>C: T/T morning preference; C/C delayed"
#   "rs2287161|12:26555282|C/C|T/T|CRY1 c.*682C>T: T/T delayed phase"
# )

# # ───────────────── Cognitive ───────────────────────────
# cognitive_snps=(
#   "rs6265|11:27679916|G/G|A/A|BDNF Val66Met duplicated for cognition"
#   "rs429358|19:45411941|T/T|C/C|APOE ε4 risk allele (see haplotype)"
#   "rs7412|19:45412079|T/T|C/C|APOE ε2 protective (haplotype context)"
#   "rs1800497|11:113400106|G/G|A/A|ANKK1 Taq1A D2 receptor density"
#   "rs4680|22:19963748|G/G|A/A|COMT Val158Met executive‑function edge"
#   "rs2760118|1:228960692|C/C|T/T|DISC1 neuronal signalling"
#   "rs1468363|12:82179232|T/T|C/C|KIBRA episodic memory"
# )

# # ───────────────── Nutrition / Diet ────────────────────
# nutrition_snps=(
#   "rs5082|1:230295691|T/T|C/C|APOA2: T/T benefit from low‑sat‑fat diet"
#   "rs4988235|2:136608646|T/T|C/C|LCT lactase persistence"
#   "rs713598|7:141673345|C/C|G/G|TAS2R38: C/C non‑taster (veg intake ↓)"
#   "rs1800562|6:26093141|G/G|A/A|HFE C282Y: A/A haemochromatosis risk"
#   "rs1801282|3:12393125|G/G|C/C|PPARG Pro12Ala duplicated for diet/insulin"
#   "rs9939609|16:53820527|T/T|A/A|FTO: A allele ↑BMI risk"

#   # Additional Nutrition Variants
#   "rs601338|19:49198677|G/G|A/A|FUT2 G428A: G/G secretor; A/A non‑secretor"
#   "rs1799945|6:26093088|C/C|G/G|HFE H63D: C/C normal iron; G/G overload risk"
#   "rs17782313|18:57884661|C/C|T/T|MC4R: T/T appetite dysregulation"
#   "rs58542926|19:19270377|C/C|T/T|TM6SF2 E167K: T/T ↑NAFLD risk"
# )

#########################  SNP CATALOGS  ####################################

# Define separate arrays for each section
# Format: rsID|chr:pos|protective_gt|risk_gt|interpretation

# ───────────────────── Metabolism ─────────────────────
metabolism_snps=(
  "rs1801133|1:11856378|C/C|T/T|MTHFR C677T: C/C normal folate; T/T ~40 % activity"
  "rs1801131|1:11854476|A/A|C/C|MTHFR A1298C: A/A normal; C/C reduced activity"
  "rs1050152|7:30645223|T/T|C/C|SLC22A4 1672C/T: T/T ↑carnitine transport; C/C ↓"
  "rs1799930|8:18258103|G/G|A/A|NAT2 G590A: G/G rapid; A/A slow acetylator"
  "rs1799931|8:18257854|G/G|A/A|NAT2 G857A: G/G rapid; A/A slow acetylator"
  "rs662|7:94937446|G/G|A/A|PON1 Q192R: G/G better organophosphate detox; A/A reduced"

  # Additional (Methylation / Hormone) Variants
  "rs1805087|1:236729215|A/A|G/G|MTR A2756G: A/A normal function; G/G higher homocysteine"
  "rs1801394|5:7875012|A/A|G/G|MTRR A66G: A/A normal function; G/G reduced B12 regeneration"
  "rs3733890|5:139541247|G/G|A/A|BHMT c.742G>A: G/G normal betaine remethylation; A/A lower"
  "rs234706|21:44489385|C/C|T/T|CBS c.699C>T: C/C normal transsulfuration; T/T upregulated"
  "rs743572|10:104513936|T/T|C/C|CYP17A1 −34T>C: T/T normal androgen synthesis; C/C ↑androgens"
  "rs1799941|17:7632191|A/A|G/G|SHBG c.-68G>A: A/A higher SHBG; G/G lower"
)

# ─────────────────── Inflammation ───────────────────
inflammation_snps=(
  "rs20417|1:186722460|G/G|C/C|PTGS2 −765 G>C: G/G lower COX‑2; C/C higher"
  "rs1800795|7:22766645|G/G|C/C|IL6 −174 G>C: G/G lower IL‑6; C/C higher"
  "rs1800629|6:31543059|G/G|A/A|TNF −308 G>A: G/G lower TNF‑α; A/A higher"
  "rs16944|2:113594387|C/C|T/T|IL1B −511 C>T: C/C lower IL‑1β; T/T higher"
  "rs1143634|2:113590390|C/C|T/T|IL1B +3954 C>T: C/C lower IL‑1β; T/T higher"
  "rs4553808|4:123377980|A/A|G/G|CTLA4 −1661 A>G: A/A better T‑cell control; G/G reduced"
  "rs4580644|11:120154442|A/A|G/G|PAWR: A/A protective (↓NF‑κB); G/G higher inflam‑risk"

  # Additional Inflammation Variants
  "rs1800896|1:206945946|G/G|A/A|IL10 −1082 G>A: G/G higher IL‑10; A/A lower"
  "rs1205|1:159682037|C/C|T/T|CRP 3' UTR variant: C/C lower CRP; T/T higher"
  "rs4986790|9:117713086|C/C|T/T|TLR4 Asp299Gly: C/C normal LPS response; T/T altered"
)

# ─────────────────── Cardiovascular ───────────────────
cardiovascular_snps=(
  "rs1801282|3:12393125|G/G|C/C|PPARG Pro12Ala: G/G protective; C/C ↑T2D"
  "rs429358|19:45411941|T/T|C/C|APOE ε4 marker; C allele ↑Alzheimer/lipids"
  "rs7412|19:45412079|T/T|C/C|APOE ε2 marker; T allele ↓cholesterol"
  "rs1799983|7:150696111|G/G|T/T|NOS3 Glu298Asp: G/G normal NO; T/T reduced"
  "rs5186|3:148459988|A/A|C/C|AGTR1 A1166C: A/A benign; C/C ↑hypertension"
  "rs1042713|5:148826877|G/G|A/A|ADRB2 Gly16Arg: G/G better β‑agonist response"
  # "rs1799752|17:61566095|I/I|D/D|ACE I/D: I/I lower ACE; D/D higher – indel ‑ handle separately"

  # Additional Cardiovascular Variants
  "rs11591147|1:55505621|T/T|C/C|PCSK9 R46L: T/T ↓LDL; protective"
  "rs2228671|19:11211158|C/C|T/T|LDLR c.1773C>T: C/C normal clearing; T/T decreased clearance"
  "rs6025|1:169519049|G/G|A/A|Factor V Leiden R506Q: G/G normal clotting; A/A risk"
  "rs1799963|11:4672404|G/G|A/A|Prothrombin G20210A: G/G normal clotting; A/A risk"
)

# ───────────────── Neurotransmitters ───────────────────
neurotransmitters_snps=(
  "rs6280|11:113283688|G/G|A/A|DRD3 Ser9Gly: G/G stronger binding; A/A weaker"
  "rs4680|22:19963748|G/G|A/A|COMT Val158Met: G/G fast DA clearance; A/A slow"
  "rs6265|11:27679916|G/G|A/A|BDNF Val66Met: G/G higher plasticity; A/A lower"
  "rs1800497|11:113400106|G/G|A/A|ANKK1 Taq1A: G/G higher D2 density"
  "rs1799971|6:154360797|A/A|G/G|OPRM1 A118G: A/A normal opioids; G/G altered"
  "rs4633|22:19962712|C/C|T/T|COMT rs4633: C/C low DA; T/T higher"
)

# ───────────────────── Antioxidant ─────────────────────
antioxidant_snps=(
  "rs1050450|3:49394834|C/C|T/T|GPX1 Pro198Leu: C/C optimal GPX1; T/T reduced"
  "rs1695|11:67352689|A/A|G/G|GSTP1 Ile105Val: A/A normal detox; G/G altered"
  "rs4880|6:160113872|C/C|T/T|SOD2 Ala16Val: C/C efficient Mn‑SOD import; T/T reduced"
  "rs1800668|3:49357401|A/A|G/G|GPX1 −602 A>G: A/A higher GPX1; G/G lower"
)

# ───────────────────── Vitamin D / Folate ──────────────
vitamin_snps=(
  "rs731236|12:48238757|T/T|C/C|VDR TaqI: T/T ↑VDR sensitivity; C/C lower"
  "rs2228570|12:48272895|G/G|A/A|VDR FokI: G/G active VDR; A/A less active"
  "rs7041|4:72618334|G/G|T/T|GC/DBP: G/G better D transport; T/T reduced"
  "rs1544410|12:48239835|A/A|G/G|VDR BsmI: A/A higher receptor activity"
)

# ───────────────── Drug Metabolism ─────────────────────
drug_snps=(
  "rs12248560|10:94761900|T/T|C/C|CYP2C19*17: T/T ultra‑rapid; C/C normal"
  "rs4244285|10:94781859|G/G|A/A|CYP2C19*2: G/G normal; A/A loss‑of‑fx"
  "rs4986893|10:94780653|G/G|A/A|CYP2C19*3: G/G normal; A/A loss‑of‑fx"
  "rs28399504|10:94781858|A/A|G/G|CYP2C19*4: A/A normal; G/G loss‑of‑fx"
  "rs1057910|10:96741053|A/A|C/C|CYP2C9*3: A/A normal; C/C ↓activity"
  "rs1799853|10:96702047|C/C|T/T|CYP2C9*2: C/C normal; T/T ↓activity"
  "rs3892097|22:42524947|G/G|A/A|CYP2D6*4: G/G normal; A/A poor metaboliser"
  "rs762551|15:75041917|A/A|C/C|CYP1A2*1F: A/A rapid caffeine clearance"
  "rs776746|7:99361466|A/A|G/G|CYP3A5*1 expressor vs non‑expresser"
)

# ─────────────────── Athletic / Muscle ─────────────────
athletic_snps=(
  "rs1815739|11:66560624|C/C|T/T|ACTN3 R577X: C/C power; T/T endurance"
  "rs8192678|4:89011240|G/G|A/A|PPARGC1A Gly482Ser: G/G endurance ↑"
  "rs1042713|5:148826877|G/G|A/A|ADRB2 Gly16Arg: G/G better training response"
  "rs4994|8:37823798|T/T|C/C|ADRB3 Trp64Arg: T/T normal thermogenesis"
  "rs6265|11:27679916|G/G|A/A|BDNF Val66Met: G/G quicker recovery"
)

# ───────────────────── Longevity ───────────────────────
longevity_snps=(
  "rs2802292|6:108587315|G/G|T/T|FOXO3 G allele linked to longevity"
  "rs1042522|17:7579472|G/G|C/C|TP53 Pro72Arg: G/G efficient apoptosis"
  "rs929347|3:57126395|A/A|G/G|HESX1 A/A longer lifespan (pop‑specific)"
  "rs10811661|9:22134094|C/C|T/T|CDKN2A/B T/T ↑T2D risk"
  "rs4880|6:160113872|C/C|T/T|SOD2 Ala16Val duplicated for anti‑oxidant"
)

# ───────────────────── Sleep / Chrono ──────────────────
sleep_snps=(
  "rs73598374|6:90503222|G/G|A/A|ADA G22A: G/G normal adenosine; A deeper sleep"
  "rs12413112|10:13124039|G/G|A/A|CTNNA3 G allele better sleep quality"
  "rs324981|2:27494431|A/A|T/T|NPSR1 Asn107Ile: A/A normal; T/T longer latency"
  "rs10830963|11:92711137|C/C|G/G|MTNR1B C/C normal glucose chronotype"

  # Additional Circadian Variants
  "rs1801260|4:56158451|T/T|C/C|CLOCK 3111T>C: T/T morning preference; C/C delayed"
  "rs2287161|12:26555282|C/C|T/T|CRY1 c.*682C>T: C/C normal; T/T delayed phase"
)

# ───────────────── Cognitive ───────────────────────────
cognitive_snps=(
  "rs6265|11:27679916|G/G|A/A|BDNF Val66Met duplicated for cognition"
  "rs429358|19:45411941|T/T|C/C|APOE ε4 risk allele (see haplotype)"
  "rs7412|19:45412079|T/T|C/C|APOE ε2 protective (haplotype context)"
  "rs1800497|11:113400106|G/G|A/A|ANKK1 Taq1A D2 receptor density"
  "rs4680|22:19963748|G/G|A/A|COMT Val158Met executive‑function edge"
  "rs2760118|1:228960692|C/C|T/T|DISC1 neuronal signalling"
  "rs1468363|12:82179232|T/T|C/C|KIBRA episodic memory"
)

# ───────────────── Nutrition / Diet ────────────────────
nutrition_snps=(
  "rs5082|1:230295691|T/T|C/C|APOA2: T/T benefit from low‑sat‑fat diet"
  "rs4988235|2:136608646|T/T|C/C|LCT lactase persistence"
  "rs713598|7:141673345|C/C|G/G|TAS2R38: C/C non‑taster (veg intake ↓)"
  "rs1800562|6:26093141|G/G|A/A|HFE C282Y: A/A haemochromatosis risk"
  "rs1801282|3:12393125|G/G|C/C|PPARG Pro12Ala duplicated for diet/insulin"
  "rs9939609|16:53820527|T/T|A/A|FTO: A allele ↑BMI risk"

  # Additional Nutrition Variants
  "rs601338|19:49198677|G/G|A/A|FUT2 G428A: G/G secretor; A/A non‑secretor"
  "rs1799945|6:26093088|C/C|G/G|HFE H63D: C/C normal iron; G/G overload risk"
  "rs17782313|18:57884661|C/C|T/T|MC4R near gene: T/T appetite dysregulation"
  "rs58542926|19:19270377|C/C|T/T|TM6SF2 E167K: T/T ↑NAFLD risk"
)

########################  config / colours  ##########################
DEFAULT_REPORT_FILE="./reports/microarray/all_snp_report.txt"
GREEN='\033[0;32m'; YELLOW='\033[1;33m'; RED='\033[0;31m'; GRAY='\033[0;37m'; NC='\033[0m'

########################  CLI parsing  ###############################
SECTION="all"; REPORT_FILE="$DEFAULT_REPORT_FILE"; LIST_ONLY=false
while getopts ":f:l" opt; do
  case $opt in
    f) REPORT_FILE="$OPTARG" ;;
    l) LIST_ONLY=true ;;
    *) echo "usage: $0 [-f file] [-l] [section]"; exit 1 ;;
  esac
done; shift $((OPTIND-1)); [[ $1 ]] && SECTION="$1"

SECTIONS=(metabolism inflammation cardiovascular neurotransmitters antioxidant \
          vitamin drug athletic longevity sleep cognitive nutrition all)

$LIST_ONLY && { printf "%s\n" "${SECTIONS[@]}"; exit 0; }
[[ -f $REPORT_FILE ]] || { echo "report not found: $REPORT_FILE" >&2; exit 1; }

########################  helper: DNA complement  ####################
get_complement(){ case $1 in A)echo T;;T)echo A;;C)echo G;;G)echo C;;*)echo "$1";; esac; }

lookup_snp(){
  IFS="|" read -r rsid _ protect risk interp <<<"$1"
  local rec; rec=$(awk -F'\t' -v id="$rsid" '$1==id{print;exit}' "$REPORT_FILE")
  [[ -z $rec ]] && {
    printf "%-10s ${GRAY}%-12s${NC} %s\n" "$rsid" "Not‑found" "$interp"
    return
  }

  read -r _ _ _ gt_raw <<<"$rec"
  # Validate genotype is e.g. 'AA', 'AG', etc.
  [[ ! $gt_raw =~ ^[ACGT]{2}$ ]] && {
    printf "%-10s ${GRAY}%-12s${NC} %s\n" "$rsid" "Invalid" "$interp"
    return
  }

  # Prepare all forms
  local g1=${gt_raw:0:1} g2=${gt_raw:1:1}
  local user="$g1/$g2"
  local rev="$g2/$g1"
  local cu="$(get_complement "$g1")/$(get_complement "$g2")"
  local cv="$(get_complement "$g2")/$(get_complement "$g1")"

  # Decide which representation to *display* depending on what matches protect/risk
  local display_gt="$user"
  local status="VARIANT"    # default
  local col="$GRAY"

  # 1) Check if user or rev already matches the protect/risk
  if [[ $user == $protect || $user == $risk ]]; then
    display_gt="$user"
  elif [[ $rev == $protect || $rev == $risk ]]; then
    display_gt="$rev"
  # 2) Else check complement forms
  elif [[ $cu == $protect || $cu == $risk ]]; then
    display_gt="$cu"
  elif [[ $cv == $protect || $cv == $risk ]]; then
    display_gt="$cv"
  fi

  # Now set status (GOOD, RISK, etc.) using the final display genotype
  if [[ $display_gt == $protect ]]; then
    status="GOOD"
    col=$GREEN
  elif [[ $display_gt == $risk ]]; then
    status="RISK"
    col=$RED
  elif [[ $g1 != "$g2" ]]; then
    # Heterozygote => if not matched protect/risk exactly, call it "CARRIER"
    status="CARRIER"
    col=$YELLOW
  fi

  printf "%-10s ${col}%-12s${NC} %s\n" "$rsid" "$display_gt" "$interp"
}

print_section(){ printf "\n===== %s =====\n\n" "$1"; }

################  SERIALISE arrays so child shells can use them  ######
SNP_TABLES=$(declare -p \
  metabolism_snps inflammation_snps cardiovascular_snps neurotransmitters_snps \
  antioxidant_snps vitamin_snps drug_snps athletic_snps longevity_snps \
  sleep_snps cognitive_snps nutrition_snps)
export SNP_TABLES

################  process_section (runs in parallel)  #################
process_section(){
  eval "$SNP_TABLES"                      # recreate arrays in child
  local sect="$1" list
  case $sect in
    metabolism)        list=("${metabolism_snps[@]}");;
    inflammation)      list=("${inflammation_snps[@]}");;
    cardiovascular)    list=("${cardiovascular_snps[@]}");;
    neurotransmitters) list=("${neurotransmitters_snps[@]}");;
    antioxidant)       list=("${antioxidant_snps[@]}");;
    vitamin)           list=("${vitamin_snps[@]}");;
    drug)              list=("${drug_snps[@]}");;
    athletic)          list=("${athletic_snps[@]}");;
    longevity)         list=("${longevity_snps[@]}");;
    sleep)             list=("${sleep_snps[@]}");;
    cognitive)         list=("${cognitive_snps[@]}");;
    nutrition)         list=("${nutrition_snps[@]}");;
    *) return;;
  esac
  [[ ${#list[@]} -eq 0 ]] && return
  print_section "${sect^}"
  for snp in "${list[@]}"; do lookup_snp "$snp"; done
}

export -f lookup_snp print_section process_section get_complement
export REPORT_FILE GREEN YELLOW RED GRAY NC

########################  MAIN  #######################################
printf "%-10s %-12s %s\n" "RSID" "Genotype" "Interpretation"

selected=()
for s in "${SECTIONS[@]}"; do
  [[ $SECTION == "all" || $SECTION == "$s" ]] && [[ $s != "all" ]] && selected+=("$s")
done

printf "%s\n" "${selected[@]}" | parallel --keep-order -j8 process_section {}

echo -e "\nLegend: ${GREEN}GOOD${NC}, ${RED}RISK${NC}, ${YELLOW}CARRIER${NC}, ${GRAY}VARIANT${NC}"