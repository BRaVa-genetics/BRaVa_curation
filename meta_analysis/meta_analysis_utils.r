# Renaming errors
renaming_group_list <- list(
	`damaging_missense_or_protein_altering` = c("damaging_missense", "missenseLC"),
	`other_missense_or_protein_altering` = c("other_missense"),
	`pLoF;damaging_missense_or_protein_altering` = c("pLoF;damaging_missense", "pLoF;missenseLC", "damaging_missense;pLoF"),
	`pLoF;damaging_missense_or_protein_altering;other_missense_or_protein_altering;synonymous` = c(
		"pLoF;damaging_missense;other_missense;synonymous",
		"pLoF;damaging_missense;other_missense;synonymous;pLoF"),
	Cauchy = c(NA)
)

renaming_header_list <- list(
	`Region` = c("gene"),
	`Group` = c("annot"),
	`max_MAF` = c("max_maf"),
	`Pvalue` = c("p_value"),
	`Pvalue_Burden` = c("p_value_burden"),
	`Pvalue_SKAT` = c("p_value_skat"),
	`BETA_Burden` = c("beta_burden"),
	`SE_Burden` = c("se_burden"),
	`MAC` = c("mac"),
	`MAC_case` = c("mac_case"),
	`MAC_control` = c("mac_control"), 
	`Number_rare` = c ("rare_var_count"),
	`Number_ultra_rare` = c("ultrarare_var_count")
)

renaming_variant_header_list <- list(
	`CHR` = c("chr", "chromosome"),
	`POS` = c("pos", "BP", "bp", "base_pair_location"),
	`MarkerID` = c("varID", "variant_id"),
	`Allele1` = c("A1", "other_allele"),
	`Allele2` = c("A2", "effect_allele"),
	`AC_Allele2` = c("AC_A2", "effect_allele_count"),
	`AF_Allele2` = c("AF_A2", "effect_allele_frequency"),
	`MissingRate` = c("missing", "missing_rate"),
	`BETA` = c("beta"),
	`SE` = c("se", "standard_error"),
	`Tstat` = c("t_stat", "t_statistic", "tstat"),
	`var` = c("Var", "variance"),
	`p.value` = c("P_value", "p_value"),
	`p.value.NA` = c("P_value.NA", "p_value_na"),
	`Is.SPA` = c("SPA", "is_spa_test"),
	`AF_case` = c("case_AF", "allele_freq_case"),
	`AF_ctrl` = c("ctrl_AF", "allele_freq_ctrl"),
	`N_case` = c("n_case"),
	`N_ctrl` = c("n_control", "n_ctrl"),
	`N_case_hom` = c("n_case_hom"),
	`N_case_het` = c("n_case_het"),
	`N_ctrl_hom` = c("n_control_hom", "n_ctrl_hom"),
	`N_ctrl_het` = c("n_control_het", "n_ctrl_het"),
	`N` = c("n")
)

# Change the above to restrict to the set that best represents the trait of interest
aou_curated <- list(
	AMD = "SO_374.51",
	Asth = "RE_475",
	AFib = "CV_416.21",
	BenCervUterNeo = "218",
	BenIntNeo = "208",
	BreastCanc = "174",
	CervCanc = "180.1",
	COPD = "RE_474",
	CRF = "585.3",
	ColonRectCanc = "CA_101.41",
	CAD = "CV_404.11",
	EFRMB = "GU_626.1",
	FemInf = "GU_629",
	Gout = "MS_703.1",
	HF = "CV_424",
	HTN = "CV_401",
	IBD = "GI_522.1",
	IFHern = "GI_520.11",
	ILDSarc = "RE_481",
	MatHem = "PP_904.1",
	NonRheuValv = "747.12",
	Pancreat = "GI_554.1",
	PeptUlcer = "GI_513",
	PAD = "CV_436",
	Psori = "DE_664.4",
	RheumHeaDis = "CV_400.2",
	RheumArth = "MS_705.1",
	Stroke = "CV_431.1",
	T2Diab = "EM_202.2",
	Urolith = "GU_585",
	VaricVeins = "CV_444.1",
	VTE = "CV_440.1"
)

# Expected annotation names
correct_names <- c("pLoF",
	"damaging_missense_or_protein_altering",
	"other_missense_or_protein_altering",
	"synonymous",
	"pLoF;damaging_missense_or_protein_altering",
	"pLoF;damaging_missense_or_protein_altering;other_missense_or_protein_altering;synonymous",
	"Cauchy")

requested_max_MAFs <- c(1e-4, 1e-3, 1e-2)

# Expected columns of gene results files
default_gene_result_columns <- list(
	binary = c(
		"Region", "Group", "max_MAF", "Pvalue",
		"Pvalue_Burden", "Pvalue_SKAT","BETA_Burden", "SE_Burden",
		"MAC", "MAC_case", "MAC_control", "Number_rare", "Number_ultra_rare"),
	continuous = c(
		"Region", "Group", "max_MAF", "Pvalue",
		"Pvalue_Burden", "Pvalue_SKAT","BETA_Burden", "SE_Burden",
		"MAC", "Number_rare", "Number_ultra_rare"),
	minimal = c(
		"Region", "Group", "max_MAF", "Pvalue",
		"Pvalue_Burden", "Pvalue_SKAT", "BETA_Burden", "SE_Burden")
)

# Expected columns of variant results files
default_variant_result_columns <- list(
	binary = c(
		"CHR", "POS", "MarkerID", "Allele1", "Allele2", "AC_Allele2",
		"AF_Allele2", "MissingRate", "BETA", "SE", "Tstat", "var",
		"p.value", "p.value.NA", "Is.SPA", "AF_case", "AF_ctrl",
		"N_case", "N_ctrl", "N_case_hom", "N_case_het",
		"N_ctrl_hom", "N_ctrl_het"),
	continuous = c(
		"CHR", "POS", "MarkerID", "Allele1", "Allele2", "AC_Allele2",
		"AF_Allele2", "MissingRate", "BETA", "SE", "Tstat", "var",
		"p.value", "N"),
	minimal = c(
		"CHR", "POS", "MarkerID", "Allele1", "Allele2", "BETA", "SE",
		"p.value")
)

renaming_plot_group_list <- list(
    damaging_missense_or_protein_altering = "Damaging missense or protein altering",
    other_missense_or_protein_altering = "Other missense or protein altering",
    synonymous = "Synonymous",
    pLoF = "pLoF",
    `pLoF;damaging_missense_or_protein_altering` = "pLoF; damaging missense or protein altering",
	`pLoF;damaging_missense_or_protein_altering;other_missense_or_protein_altering;synonymous` = "pLoF; damaging missense or protein altering; other missense or protein altering; synonymous"
)

renaming_plot_biobank_list <- list(
    `all-of-us` = "All of Us",
	`alspac` = "ALSPAC",
	`biome` = "BioMe",
	`bbj` = "Biobank Japan",
	`ckb` = "CKB",
	`ccpm` = "CCPM",
	`decode` = "DECODE",
	`egcut` = "EGCUT",
	`dan-rav` = "Dan-RaV",
	`genes-and-health` = "Genes & Health",
	`gel` = "Genomics England",
	`mgbb` = "MGBB",
	`pmbb` = "PMBB",
	`qatar-genomes` = "Qatar Genome",
	`uk-biobank` = "UK Biobank",
	`viking-genes` = "Viking Genes")

renaming_phenotype_list <- list(
	`AAA` = "Abdominal aortic aneurysm",
	`AcApp` = "Acute appendicitis",
	`AcuLymLeuk` = "Acute lymphoid leukemia",
	`Adenomy` = "Adenomyosis",
	`AMD` = "Age-related macular degeneration",
	`ALamy` = "AL amyloidosis",
	`AUD` = "Alcohol use disorder",
	`AloAre` = "Alopecia areata",
	`AnoNer` = "Anorexia nervosa",
	`AoSten` = "Aortic stenosis",
	`Asth` = "Asthma",
	`AtopDis` = "Atopic diseases",
	`AFib` = "Atrial fibrillation",
	`ADHD` = "ADHD",
	`ASD` = "Autism spectrum disorder",
	`BCLL` = "B-cell chronic lymphocytic leukemia",
	`BenCervUterNeo` = "Benign and in situ cervical and uterine neoplasms",
	`BenIntNeo` = "Benign and in situ intestinal neoplasms",
	`BenNodGoit` = "Benign nodular goiter",
	`BladCanc` = "Bladder cancer",
	`BrainCNSCanc` = "Brain and central nervous system cancer",
	`BreastCanc` = "Breast cancer",
	`BrugSynd` = "Brugada syndrome",
	`BuliNer` = "Bulimia nervosa",
	`BullPemph` = "Bullous pemphigoid",
	`CarShock` = "Cardiogenic shock",
	`HCM` = "Cardiomyopathy (hypertrophic, obstructive)",
	`CRVO` = "Central retinal vein occlusion",
	`CervCanc` = "Cervical cancer",
	`CML` = "Chronic myeloid leukemia",
	`COPD` = "Chronic obstructive pulmonary disease",
	`CRF` = "Chronic renal failure",
	`CoffSirSynd` = "Coffin-Siris syndrome",
	`ColonRectCanc` = "Colon and rectum cancer",
	`CAD` = "Coronary artery disease",
	`CCANS` = "Corpus callosum agenesis-neuronopathy syndrome",
	`EatDis` = "Eating disorders",
	`Endocar` = "Endocarditis",
	`Endometr` = "Endometriosis",
	`EsophCanc` = "Esophageal cancer",
	`EssThrom` = "Essential thrombocythemia",
	`EFRMB` = "Excess, frequent and irregular menstrual bleeding",
	`FSP` = "Familial spontaneous pneumothorax",
	`FemInf` = "Female infertility",
	`FemInfAC` = "Female infertility (anatomic causes only)",
	`FolLymph` = "Follicular lymphoma",
	`Gout` = "Gout",
	`GravesDis` = "Graves disease",
	`HemoChromo` = "Haemochromatosis",
	`HF` = "Heart failure",
	`HepCarcin` = "Hepatocellular carcinoma",
	`HTN` = "Hypertension",
	`HHD` = "Hypertensive heart disease",
	`HypoThyr` = "Hypothyroidism",
	`HypoThyrSec` = "Hypothyroidism, secondary",
	`IPF` ="Idiopathic pulmonary fibrosis",
	`ITP` = "Immune thrombocytopenia",
	`IBD` = "Inflammatory bowel disease",
	`IFHern` = "Inguinal, femoral, and abdominal hernia",
	`ILDSarc` = "Interstitial lung disease and pulmonary sarcoidosis",
	`IodDef` = "Iodine deficiency",
	`KabSynd` = "Kabuki Syndrome",
	`KidCanc` = "Kidney cancer",
	`KleefSynd` = "Kleefstra Syndrome",
	`LaryxCanc` = "Larynx cancer",
	`Leuk` = "Leukemia",
	`LiverCanc` = "Liver cancer",
	`LiverFibCirr` = "Liver fibrosis/cirrhosis",
	`LongQTSynd` = "Long QT syndrome",
	`LymphThyrit` = "Lymphocytic thyroiditis",
	`MalInf` = "Male infertility",
	`MatHem` = "Maternal hemorrhage",
	`MatHypDis` = "Maternal hypertensive disorders",
	`MODYDiab` = "MODY diabetes",
	`MultiMyel` = "Multiple myeloma",
	`MS` = "Multiple sclerosis",
	`MECS` = "Chronic fatigue syndrome",
	`Myocard` = "Myocarditis",
	`Narco1` = "Narcolepsy type 1",
	`NonFuncPitAd` = "Non-functioning pituitary adenoma",
	`NHL` = "Non-Hodgkin lymphoma",
	`NonPapTCCBlad` = "Non-papillary transitional cell carcinoma of the bladder",
	`NonRheuValv` = "Non-rheumatic valvular heart disease",
	`OUD` = "Opioid use disorder",
	`OCD` = "Osteochondritis dissecans",
	`OvCanc` = "Ovarian cancer",
	`Pancreat` = "Pancreatitis",
	`ParkDis` = "Parkinson's disease",
	`PeptUlcer` = "Peptic ulcer disease",
	`PAD` = "Peripheral artery disease",
	`PlacInsuf` = "Placental insufficiency",
	`PCOS` = "Polycystic ovarian syndrome",
	`PolycythVera` = "Polycythemia vera",
	`Preeclamps` = "Preeclampsia",
	`PregLoss` =  "Pregnancy loss",
	`POAG` = "Primary open angle glaucoma",
	`PrimSjoSynd` = "Primary Sjögren syndrome",
	`Prolactinom` = "Prolactinoma",
	`Psori` = "Psoriasis",
	`RheumHeaDis` = "Rheumatic heart disease",
	`RheumArth` = "Rheumatoid arthritis",
	`RomWardSynd` = "Romano-Ward syndrome",
	`Sarcoid` = "Sarcoidosis",
	`SebDerm` = "Seborrhoeic dermatitis",
	`SpinaBifAp` = "Spina bifida aperta",
	`StomCanc` =  "Stomach cancer",
	`Stroke` = "Stroke",
	`SLE` = "Systemic lupus erythematosus",
	`TAAD` = "Thoracic aortic aneurysm and dissection",
	`ThyroCanc` = "Thyroid cancer",
	`T2Diab` = "Type 2 diabetes",
	`Urolith` = "Urolithiasis",
	`UterCanc` = "Uterine cancer",
	`VaricVeins` = "Varicose veins",
	`VTE` = "Venous thromboembolism",
	`ALT` = "Alanine transaminase",
	`AlcCons` = "Alcohol consumption (drinks per week)",
	`AST` = "Aspartate aminotransferase",
	`BMI` = "Body mass index",
	`CRP` = "C-reactive protein",
	`CACS` = "Coronary artery calcium score",
	`CK` = "Creatine kinase",
	`HDLC` = "HDL cholesterol",
	`Height` = "Height",
	`LDLC` = "LDL cholesterol",
	`TChol` = "Total cholesterol",
	`TG` = "Triglycerides",
	`WHRBMI` = "Waist to hip ratio adjusted for BMI",
	`LVH` = "Left ventricular hypertrophy",
	`Append` = "Appendectomy",
	`HipRep` = "Hip replacement",
	`CogAbil` = "Cognitive ability at different ages",
	`EduAtt` = "Educational attainment",
	`PsySymp` = "Psychotic symptoms",
	`SchGrades` = "School grades",
	`SCDCAT` = "Social and communication disorders checklist/autistic traits"
)

phenotype_broad_categories <- list(
	AMD	= "Sense organs",
	Asth = "Respiratory",
	AFib = "Cardiovascular",
	BenCervUterNeo = "Neoplasms",
	BenIntNeo = "Neoplasms",
	BreastCanc = "Neoplasms",
	CervCanc = "Neoplasms",
	COPD = "Respiratory",
	CRF = "Genitourinary",
	ColonRectCanc = "Neoplasms",
	CAD = "Cardiovascular",
	EFRMB = "Genitourinary",
	FemInf = "Genitourinary",
	Gout = "Musculoskeletal",
	HF = "Cardiovascular",
	HTN = "Cardiovascular",
	IBD = "Gastrointestinal",
	IFHern = "Gastrointestinal",
	ILDSarc = "Respiratory",
	MatHem = "Pregnancy",
	NonRheuValv = "Cardiovascular",
	Pancreat = "Gastrointestinal",
	PeptUlcer = "Gastrointestinal",
	PAD = "Cardiovascular",
	Psori = "Dermatological",
	RheumHeaDis = "Cardiovascular",
	RheumArth = "Musculoskeletal",
	Stroke = "Cardiovascular",
	T2Diab = "Endocrine/Metabolic",
	Urolith = "Genitourinary",
	VaricVeins = "Cardiovascular",
	VTE = "Cardiovascular",
	ALT = "Gastrointestinal",
	AlcCons = "Endocrine/Metabolic",
	AST = "Gastrointestinal",
	BMI = "Endocrine/Metabolic",
	CRP = "Cardiovascular",
	HDLC = "Cardiovascular",
	Height = "Musculoskeletal",
	LDLC = "Cardiovascular",
	TChol = "Cardiovascular",
	TG = "Cardiovascular",
	WHRBMI = "Endocrine/Metabolic",
	HipRep = "Musculoskeletal"
)

phenotype_class <- list(
	binary = c("AAA","AcApp","AcuLymLeuk","Adenomy","AMD","ALamy","AUD","AloAre",
		"AnoNer","AoSten","Asth","AtopDis","AFib","ADHD","ASD","BCLL","BenCervUterNeo",
		"BenIntNeo","BenNodGoit","BladCanc","BrainCNSCanc","BreastCanc","BrugSynd",
		"BuliNer","BullPemph","CarShock","HCM","CRVO","CervCanc","CML","COPD","CRF",
		"CoffSirSynd","ColonRectCanc","CAD","CCANS","EatDis","Endocar","Endometr",
		"EsophCanc","EssThrom","EFRMB","FSP","FemInf","FemInfAC","FolLymph","Gout",
		"GravesDis","HemoChromo","HF","HepCarcin","HTN","HHD","HypoThyr","HypoThyrSec",
		"IPF","ITP","IBD","IFHern","ILDSarc","IodDef","KabSynd","KidCanc","KleefSynd",
		"LaryxCanc","Leuk","LiverCanc","LiverFibCirr","LongQTSynd","LymphThyrit",
		"MalInf","MatHem","MatHypDis","MODYDiab","MultiMyel","MS","MECS","Myocard",
		"Narco1","NonFuncPitAd","NHL","NonPapTCCBlad","NonRheuValv","OUD","OCD",
		"OvCanc","Pancreat","ParkDis","PeptUlcer","PAD","PlacInsuf","PCOS",
		"PolycythVera","Preeclamps","PregLoss","POAG","PrimSjoSynd","Prolactinom",
		"Psori","RheumHeaDis","RheumArth","RomWardSynd","Sarcoid","SebDerm",
		"SpinaBifAp","StomCanc","Stroke","SLE","TAAD","ThyroCanc","T2Diab","Urolith",
		"UterCanc","VaricVeins","VTE"),
	continuous = c("ALT", "AlcCons", "AST","BMI","CRP","CACS","CK","HDLC","Height",
		"LDLC","TChol","TG","WHRBMI", "LVH","Append","HipRep","CogAbil","EduAtt",
		"PsySymp","SchGrades","SCDCAT")
)

file_check_information <- list(
	dataset = list(
		`all-of-us` = c("all-of-us", "All-of-Us", "All-of-us"),
		`alspac` = c("alspac"),
		`biome` = c("biome", "BioMe"),
		`bbj` = c("bbj", "BBJ"),
		`ckb` = c("ckb", "CKB"),
		`ccpm` = c("ccpm", "CCPM"),
		`decode` = c("decode", "DECODE"),
		`egcut` = c("egcut", "EGCUT", "EstBB", "estbb", "ESTBB"),
		`dan-rav` = c("dan-rav", "Dan-RaV"),
		`genes-and-health` = c("genes-and-health", "GnH", "Genes-and-Health", "GNH"),
		`gel` = c("gel", "GEL"),
		`mgbb` = c("mgbb", "MGBB"),
		`pmbb` = c("pmbb", "PMBB"),
		`qatar-genomes` = c("qatar-genomes", "Qatar-genomes"),
		`uk-biobank` = c("uk-biobank", "ukbb"),
		`viking-genes` = c("viking-genes", "Viking-Genes") 
	),
	phenotype = list(
		`AAA` = c("AAA", "Abdominal aortic aneurysm (AAA)", "Abdominal_aortic_aneurysm_(AAA)"),
		`AcApp` = c("AcApp", "Acute appendicitis (AcApp)", "Acute_appendicitis_(AcApp)", "Acute_appendicitis_AcApp"),
		`AcuLymLeuk` = c("AcuLymLeuk", "Acute lymphoid leukemia", "Acute_lymphoid_leukemia"),
		`Adenomy` = c("Adenomy", "Adenomyosis"),
		`AMD` = c("AMD", "Age-related macular degeneration", "Age-related_macular_degeneration", "Age_related_macular_degeneration"),
		`ALamy` = c("ALamy", "AL amyloidosis", "AL_amyloidosis"),
		`AUD` = c("AUD", "Alcohol use disorder", "Alcohol_use_disorder"),
		`AloAre` = c("AloAre", "Alopecia Areata", "Alopecia_Areata"),
		`AnoNer` = c("AnoNer", "Anorexia nervosa", "Anorexia_nervosa"),
		`AoSten` = c("AoSten", "Aortic stenosis", "Aortic_stenosis"),
		`Asth` = c("Asth", "Asthma", "Asthma_Asthma"),
		`AtopDis` = c("AtopDis", "Atopic diseases", "Atopic_diseases"),
		`AFib` = c("AFib", "Atrial Fibrillation", "Atrial_Fibrillation"),
		`ADHD` = c("ADHD", "Attention-deficit/hyperactivity disorder (ADHD)", "Attention-deficit/hyperactivity_disorder_(ADHD)"),
		`ASD` = c("ASD", "Autism spectrum disorder (ASD)", "Autism_spectrum_disorder_(ASD)"),
		`BCLL` = c("BCLL", "B-cell chronic lymphocytic leukemia", "B-cell_chronic_lymphocytic_leukemia"),
		`BenCervUterNeo` = c("BenCervUterNeo", "Benign and in situ cervical and uterine neoplasms", "Benign_and_in_situ_cervical_and_uterine_neoplasms"),
		`BenIntNeo` = c("BenIntNeo", "Benign and in situ intestinal neoplasms", "Benign_and_in_situ_intestinal_neoplasms"),
		`BenNodGoit` = c("BenNodGoit", "Benign nodular goiter", "Benign_nodular_goiter"),
		`BladCanc` = c("BladCanc", "Bladder cancer", "Bladder_cancer"),
		`BrainCNSCanc` = c("BrainCNSCanc", "Brain and central nervous system cancer", "Brain_and_central_nervous_system_cancer"),
		`BreastCanc` = c("BreastCanc", "Breast cancer", "Breast_cancer"),
		`BrugSynd` = c("BrugSynd", "Brugada syndrome", "Brugada_syndrome"),
		`BuliNer` = c("BuliNer", "Bulimia nervosa", "Bulimia_nervosa"),
		`BullPemph` = c("BullPemph", "Bullous pemphigoid", "Bullous_pemphigoid"),
		`CarShock` = c("CarShock", "Cardiogenic shock", "Cardiogenic_shock"),
		`HCM` = c("HCM", "Cardiomyopathy (hypertrophic, obstructive) (HCM)", "Cardiomyopathy_(hypertrophic,_obstructive)_(HCM)"),
		`CRVO` = c("CRVO", "Central retinal vein occlusion", "Central_retinal_vein_occlusion"),
		`CervCanc` = c("CervCanc", "Cervical cancer", "Cervical_cancer"),
		`CML` = c("CML", "Chronic myeloid leukemia", "Chronic_myeloid_leukemia"),
		`COPD` = c("COPD", "Chronic obstructive pulmonary disease (COPD)", "Chronic_obstructive_pulmonary_disease_(COPD)", "Chronic_obstructive_pulmonary_disease_COPD", "Chronic_obstructive_pulmonary_disease"),
		`CRF` = c("CRF", "Chronic Renal Failure", "Chronic_Renal_Failure"),
		`CoffSirSynd` = c("CoffSirSynd", "Coffin-Sirus Syndrome", "Coffin-Sirus_Syndrome"),
		`ColonRectCanc` = c("ColonRectCanc", "Colon and rectum cancer", "Colon_and_rectum_cancer"),
		`CAD` = c("CAD", "Coronary artery disease", "Coronary_artery_disease"),
		`CCANS` = c("CCANS", "Corpus callosum agenesis-neuronopathy syndrome", "Corpus_callosum_agenesis-neuronopathy_syndrome"),
		`EatDis` = c("EatDis", "Eating disorders", "Eating_disorders"),
		`Endocar` = c("Endocar", "Endocarditis"),
		`Endometr` = c("Endometr", "Endometriosis"),
		`EsophCanc` = c("EsophCanc", "Esophageal cancer", "Esophageal_cancer"),
		`EssThrom` = c("EssThrom", "Essential thrombocythemia", "Essential_thrombocythemia"),
		`EFRMB` = c("EFRMB", "Excess, frequent and irregular menstrual bleeding", "Excess,_frequent_and_irregular_menstrual_bleeding", "Excess_frequent_and_irregular_menstrual_bleeding"),
		`FSP` = c("FSP", "Familial spontaneous pneumothorax", "Familial_spontaneous_pneumothorax"),
		`FemInf` = c("FemInf", "Female infertility", "Female_infertility"),
		`FemInfAC` = c("FemInfAC", "Female infertility (anatomic causes only)", "Female_infertility_(anatomic_causes_only)"),
		`FolLymph` = c("FolLymph", "Follicular lymphoma", "Follicular_lymphoma"),
		`Gout` = c("Gout"),
		`GravesDis` = c("GravesDis", "Graves disease ", "Graves_disease_"),
		`HemoChromo` = c("HemoChromo", "Haemochromatosis"),
		`HF` = c("HF", "Heart Failure (HF)", "Heart_Failure_(HF)", "Heart_Failure_HF", "Heart_Failure"),
		`HepCarcin` = c("HepCarcin", "Hepatocellular Carcinoma", "Hepatocellular_Carcinoma"),
		`HTN` = c("HTN", "Hypertension"),
		`HHD` = c("HHD", "Hypertensive heart disease", "Hypertensive_heart_disease"),
		`HypoThyr` = c("HypoThyr", "Hypothyroidism"),
		`HypoThyrSec` = c("HypoThyrSec", "Hypothyroidism, secondary (for defining exclusions in the hypothyroidism analysis)", "Hypothyroidism,_secondary_(for_defining_exclusions_in_the_hypothyroidism_analysis)"),
		`IPF` = c("IPF", "Idiopathic pulmonary fibrosis (IPF)", "Idiopathic_pulmonary_fibrosis_(IPF)"),
		`ITP` = c("ITP", "Immune thrombocytopenia", "Immune_thrombocytopenia"),
		`IBD` = c("IBD", "Inflammatory bowel disease", "Inflammatory_bowel_disease"),
		`IFHern` = c("IFHern", "Inguinal, femoral, and abdominal hernia", "Inguinal,_femoral,_and_abdominal_hernia", "Inguinal_femoral_and_abdominal_hernia"),
		`ILDSarc` = c("ILDSarc", "Interstitial lung disease and pulmonary sarcoidosis", "Interstitial_lung_disease_and_pulmonary_sarcoidosis"),
		`IodDef` = c("IodDef", "Iodine deficiency", "Iodine_deficiency"),
		`KabSynd` = c("KabSynd", "Kabuki Syndrome", "Kabuki_Syndrome"),
		`KidCanc` = c("KidCanc", "Kidney cancer", "Kidney_cancer"),
		`KleefSynd` = c("KleefSynd", "Kleefstra Syndrome", "Kleefstra_Syndrome"),
		`LaryxCanc` = c("LaryxCanc", "Larynx cancer", "Larynx_cancer"),
		`Leuk` = c("Leuk", "Leukemia"),
		`LiverCanc` = c("LiverCanc", "Liver cancer", "Liver_cancer"),
		`LiverFibCirr` = c("LiverFibCirr", "Liver fibrosis / cirrhosis", "Liver_fibrosis_/_cirrhosis"),
		`LongQTSynd` = c("LongQTSynd", "Long QT syndrome", "Long_QT_syndrome"),
		`LymphThyrit` = c("LymphThyrit", "Lymphocytic thyroiditis", "Lymphocytic_thyroiditis"),
		`MalInf` = c("MalInf", "Male infertility", "Male_infertility"),
		`MatHem` = c("MatHem", "Maternal hemorrhage", "Maternal_hemorrhage"),
		`MatHypDis` = c("MatHypDis", "Maternal hypertensive disorders", "Maternal_hypertensive_disorders"),
		`MODYDiab` = c("MODYDiab", "MODY diabetes", "MODY_diabetes"),
		`MultiMyel` = c("MultiMyel", "Multiple myeloma", "Multiple_myeloma"),
		`MS` = c("MS", "Multiple Sclerosis", "Multiple_Sclerosis"),
		`MECS` = c("MECS", "Myalgic Encephalomyelitis/Chronic Fatigue Syndrome (ME/CFS)", "Myalgic_Encephalomyelitis/Chronic_Fatigue_Syndrome_(ME/CFS)"),
		`Myocard` = c("Myocard", "Myocarditis"),
		`Narco1` = c("Narco1", "Narcolepsy type 1", "Narcolepsy_type_1"),
		`NonFuncPitAd` = c("NonFuncPitAd", "Non-functioning pituitary adenoma", "Non-functioning_pituitary_adenoma"),
		`NHL` = c("NHL", "Non-Hodgkin lymphoma", "Non-Hodgkin_lymphoma"),
		`NonPapTCCBlad` = c("NonPapTCCBlad", "Non-papillary transitional cell carcinoma of the bladder", "Non-papillary_transitional_cell_carcinoma_of_the_bladder"),
		`NonRheuValv` = c("NonRheuValv", "Non-rheumatic valvular heart disease", "Non-rheumatic_valvular_heart_disease", "Non_rheumatic_valvular_heart_disease"),
		`OUD` = c("OUD", "Opioid use disorder", "Opioid_use_disorder"),
		`OCD` = c("OCD", "Osteochondritis dissecans", "Osteochondritis_dissecans"),
		`OvCanc` = c("OvCanc", "Ovarian Cancer", "Ovarian_Cancer"),
		`Pancreat` = c("Pancreat", "Pancreatitis"),
		`ParkDis` = c("ParkDis", "Parkinson's disease", "Parkinson's_disease"),
		`PeptUlcer` = c("PeptUlcer", "Peptic ulcer disease", "Peptic_ulcer_disease"),
		`PAD` = c("PAD", "Peripheral artery disease", "Peripheral_artery_disease"),
		`PlacInsuf` = c("PlacInsuf", "Placental insufficiency", "Placental_insufficiency"),
		`PCOS` = c("PCOS", "Polycystic ovarian syndrome", "Polycystic_ovarian_syndrome"),
		`PolycythVera` = c("PolycythVera", "Polycythemia vera", "Polycythemia_vera"),
		`Preeclamps` = c("Preeclamps", "Preeclampsia"),
		`PregLoss` = c("PregLoss", "Pregnancy Loss", "Pregnancy_Loss"),
		`POAG` = c("POAG", "Primary open angle glaucoma (POAG)", "Primary_open_angle_glaucoma_(POAG)"),
		`PrimSjoSynd` = c("PrimSjoSynd", "Primary Sjögren syndrome", "Primary_Sjögren_syndrome"),
		`Prolactinom` = c("Prolactinom", "Prolactinoma"),
		`Psori` = c("Psori", "Psoriasis"),
		`RheumHeaDis` = c("RheumHeaDis", "Rheumatic heart disease", "Rheumatic_heart_disease"),
		`RheumArth` = c("RheumArth", "Rheumatoid arthritis", "Rheumatoid_arthritis"),
		`RomWardSynd` = c("RomWardSynd", "Romano-Ward syndrome", "Romano-Ward_syndrome"),
		`Sarcoid` = c("Sarcoid", "Sarcoidosis"),
		`SebDerm` = c("SebDerm", "Seborrhoeic dermatitis", "Seborrhoeic_dermatitis"),
		`SpinaBifAp` = c("SpinaBifAp", "Spina bifida aperta", "Spina_bifida_aperta"),
		`StomCanc` = c("StomCanc", "Stomach cancer", "Stomach_cancer"),
		`Stroke` = c("Stroke", "Stroke_Stroke"),
		`SLE` = c("SLE", "Systemic lupus erythematosus", "Systemic_lupus_erythematosus"),
		`TAAD` = c("TAAD", "Thoracic aortic aneurysm and dissection (TAAD)", "Thoracic_aortic_aneurysm_and_dissection_(TAAD)"),
		`ThyroCanc` = c("ThyroCanc", "Thyroid cancer (ThC)", "Thyroid_cancer_(ThC)"),
		`T2Diab` = c("T2Diab", "Type 2 diabetes", "Type_2_diabetes"),
		`Urolith` = c("Urolith", "Urolithiasis"),
		`UterCanc` = c("UterCanc", "Uterine cancer (UtC)", "Uterine_cancer_(UtC)"),
		`VaricVeins` = c("VaricVeins", "Varicose Veins", "Varicose_Veins"),
		`VTE` = c("VTE", "Venous Thromboembolism (VTE)", "Venous_Thromboembolism_(VTE)", "Venous_Thromboembolism_VTE", "Venous_Thromboembolism"),
		`ALT` = c("ALT", "Alanine transaminase", "Alanine_transaminase"),
		`AlcCons` = c("AlcCons", "Alcohol consumption (drinks per week)", "Alcohol_consumption_(drinks_per_week)", "Alcohol", "Alcohol_consumption_drinks_per_week"),
		`AST` = c("AST", "Aspartate aminotransferase (AST)", "Aspartate_aminotransferase_(AST)", "Aspartate_aminotransferase"),
		`BMI` = c("BMI"),
		`CRP` = c("CRP", "C-reactive protein (CRP)", "C-reactive_protein_(CRP)", "C_reactive_protein"),
		`CACS` = c("CACS", "Coronary artery calcium score", "Coronary_artery_calcium_score"),
		`CK` = c("CK", "Creatine kinase (CK)", "Creatine_kinase_(CK)"),
		`HDLC` = c("HDLC", "HDL"),
		`Height` = c("Height", "height"),
		`LDLC` = c("LDLC", "LDL"),
		`TChol` = c("TChol", "Total cholesterol", "Total_cholesterol"),
		`TG` = c("TG", "Triglycerides"),
		`WHRBMI` = c("WHRBMI", "WHR adjusted for BMI", "WHR_adjusted_for_BMI"),
		`LVH` = c("LVH", "Left Ventricular Hypertrophy", "Left_Ventricular_Hypertrophy"),
		`Append` = c("Append", "Appendectomy"),
		`HipRep` = c("HipRep", "Hip replacement (operation)", "Hip_replacement_(operation)", "Hip_replacement"),
		`CogAbil` = c("CogAbil", "Cognitive ability at different ages", "Cognitive_ability_at_different_ages"),
		`EduAtt` = c("EduAtt", "Educational attainment", "Educational_attainment"),
		`PsySymp` = c("PsySymp", "Psychotic symptoms", "Psychotic_symptoms"),
		`SchGrades` = c("SchGrades", "School grades", "School_grades"),
		`SCDCAT` = c("SCDCAT", "Social and Communication Disorders Checklist/Autistic traits", "Social_and_Communication_Disorders_Checklist/Autistic_traits")
	),
	sex = list(
		`ALL` = c("ALL", "all", "All"),
		`M` = c("M", "male", "Male", "MALE"),
		`F` = c("F", "female", "Female", "FEMALE")
	),
	ancestry = list(
		`AFR` = c("AFR", "afr"),
		`AMR` = c("AMR", "amr"),
		`EAS` = c("EAS", "eas"),
		`EUR` = c("EUR", "eur"),
		`MID` = c("MID", "mid"),
		`SAS` = c("SAS", "sas", "CSA", "csa")
	),
	type = list(
		`gene` = c("Gene", "gene"),
		`variant` = c("Variant", "var", "variant")
	)
)

color_amr <- '#ED1E24'
color_eur <- '#6AA5CD'
color_afr <- '#941494'
color_sas <- '#FF9912'
color_eas <- '#108C44'
color_oth <- '#ABB9B9'
color_mde <- '#33CC33'
color_asj <- 'coral'
color_nfe <- color_eur
color_fin <- '#002F6C'

pop_colors <- c('AFR' = color_afr,
               'AMR' = color_amr,
               'EAS' = color_eas,
               'FIN' = color_fin,
               'EUR' = color_nfe,
               'NEF' = color_nfe,
               'OTH' = color_oth,
               'SAS' = color_sas,
               'MDE' = color_mde,
               'ASJ' = color_asj,
               'uniform' = 'pink',
               'consanguineous' = 'pink',
               'SAS_non_consang' = 'orange')

pop_names <- c('OTH' = 'Other',
              'AFR' = 'African/African-American',
              'AMR' = 'Latino',
              'EAS' = 'East Asian',
              'FIN' = 'Finnish',
              'EUR' = 'European',
              'NFE' = 'European',
              'SAS' = 'South Asian',
              'MDE' = 'Middle Eastern',
              'ASJ' = 'Ashkenazi Jewish',
              'uniform' = 'Uniform',
              'SAS_non_consang' = 'South Asian (F < 0.05)',
              'consanguineous' = 'South Asian (F > 0.05)')

color_syn = '#AAAAAA'
color_mis = '#FF6103'
color_os = '#74099E'
color_lof = '#9D1309'
color_lc_lof = '#EE799F'
var_type_aliases = c('syn' = 'Synonymous', 'mis' = 'Missense', 'lof' = 'pLoF', 'os' = 'Other splice')
colors = c('synonymous_variant' = color_syn,
           'missense_variant' = color_mis,
           'stop_gained' = color_lof,
           'Synonymous' = color_syn,
           'Missense' = color_mis,
           'synonymous' = color_syn,
           'missense' = color_mis,
           'nonsense' = color_lof,
           'Other splice' = color_os,
           'LoF' = color_lof,
           'pLoF' = color_lof,
           'damaging_missense_or_protein_altering' = color_mis,
           'Damaging missense or protein altering' = color_mis)

case_ctrl <- c(
	"AAA",
	"AcApp",
	"AcuLymLeuk",
	"Adenomy",
	"AMD",
	"ALamy",
	"AUD",
	"AloAre",
	"AnoNer",
	"AoSten",
	"Asth",
	"AtopDis",
	"AFib",
	"ADHD",
	"ASD",
	"BCLL",
	"BenCervUterNeo",
	"BenIntNeo",
	"BenNodGoit",
	"BladCanc",
	"BrainCNSCanc",
	"BreastCanc",
	"BrugSynd",
	"BuliNer",
	"BullPemph",
	"CarShock",
	"HCM",
	"CRVO",
	"CervCanc",
	"CML",
	"COPD",
	"CRF",
	"CoffSirSynd",
	"ColonRectCanc",
	"CAD",
	"CCANS",
	"EatDis",
	"Endocar",
	"Endometr",
	"EsophCanc",
	"EssThrom",
	"EFRMB",
	"FSP",
	"FemInf",
	"FemInfAC",
	"FolLymph",
	"Gout",
	"GravesDis",
	"HemoChromo",
	"HF",
	"HepCarcin",
	"HTN",
	"HHD",
	"HypoThyr",
	"HypoThyrSec",
	"IPF",
	"ITP",
	"IBD",
	"IFHern",
	"ILDSarc",
	"IodDef",
	"KabSynd",
	"KidCanc",
	"KleefSynd",
	"LaryxCanc",
	"Leuk",
	"LiverCanc",
	"LiverFibCirr",
	"LongQTSynd",
	"LymphThyrit",
	"MalInf",
	"MatHem",
	"MatHypDis",
	"MODYDiab",
	"MultiMyel",
	"MS",
	"MECS",
	"Myocard",
	"Narco1",
	"NonFuncPitAd",
	"NHL",
	"NonPapTCCBlad",
	"NonRheuValv",
	"OUD",
	"OCD",
	"OvCanc",
	"Pancreat",
	"ParkDis",
	"PeptUlcer",
	"PAD",
	"PlacInsuf",
	"PCOS",
	"PolycythVera",
	"Preeclamps",
	"PregLoss",
	"POAG",
	"PrimSjoSynd",
	"Prolactinom",
	"Psori",
	"RheumHeaDis",
	"RheumArth",
	"RomWardSynd",
	"Sarcoid",
	"SebDerm",
	"SpinaBifAp",
	"StomCanc",
	"Stroke",
	"SLE",
	"TAAD",
	"ThyroCanc",
	"T2Diab",
	"Urolith",
	"UterCanc",
	"VaricVeins",
	"VTE"
)

cts <- c(
	"ALT",
	"AlcCons",
	"AST",
	"BMI",
	"CRP",
	"CACS",
	"CK",
	"HDLC",
	"Height",
	"LDLC",
	"TChol",
	"TG",
	"WHRBMI",
	"LVH",
	"Append",
	"HipRep",
	"CogAbil",
	"EduAtt",
	"PsySymp",
	"SchGrades",
	"SCDCAT"
	)

determine_binary_filename <- function(dataset, last_name, analysis_name, phenotype, sex,
    ancestry, n_cases, n_controls, type, date,
    method="SAIGE", freeze_number="JULY23Freeze") {
    # Format
    # [dataset].[last_name].[analysis_name].[phenotype].[freeze_number].[sex].[ancestry].[n_cases].[n_controls].[SAIGE].{gene,variant}.[YYYYMMDD].txt.gz
    return(paste(
        dataset,
        last_name,
        analysis_name,
        gsub("\\.", "_", phenotype),
        freeze_number,
        sex,
        ancestry,
        n_cases,
        n_controls,
        method,
        type,
        date,
        "txt.gz",
        sep="."))
}

determine_cts_filename <- function(dataset, last_name, analysis_name, phenotype, sex,
    ancestry, n, type, date,
    method="SAIGE", freeze_number="JULY23Freeze") {
    # Format
    # [dataset].[last_name].[analysis_name].[phenotype].[freeze_number].[sex].[ancestry].[n].[SAIGE].{gene,variant}.[YYYYMMDD].txt.gz
    return(paste(
        dataset,
        last_name,
        analysis_name,
        gsub("\\.", "_", phenotype),
        freeze_number,
        sex,
        ancestry,
        n,
        method,
        type,
        date,
        "txt.gz",
        sep="."))
}

search_for_files <- function(file)
{
    if (file.exists(paste0("gene/", file, ".txt.gz")) &
        file.exists(paste0("variant/", file, ".txt.singleAssoc.txt.gz"))) {
        gene_file <- paste0("gene/", file, ".txt.gz")
        variant_file <- paste0("variant/", file, ".txt.singleAssoc.txt.gz")
        gz <- TRUE
    } else if (file.exists(paste0("gene/", file, ".txt")) &
               file.exists(paste0("variant/", file, ".txt.singleAssoc.txt"))) {
        gene_file <- paste0("gene/", file, ".txt")
        variant_file <- paste0("variant/", file, ".txt.singleAssoc.txt")
        gz <- FALSE
    } else if (file.exists(paste0("gene/", file, "_F.txt.gz")) &
               file.exists(paste0("variant/", file, "_F.txt.singleAssoc.txt.gz"))) {
        gene_file <- paste0("gene/", file, "_F.txt.gz")
        variant_file <- paste0("variant/", file, "_F.txt.singleAssoc.txt.gz")
        gz <- TRUE
    } else if (file.exists(paste0("gene/", file, "_F.txt")) &
               file.exists(paste0("variant/", file, "_F.txt.singleAssoc.txt"))) {
        gene_file <- paste0("gene/", file, "_F.txt")
        variant_file <- paste0("variant/", file, "_F.txt.singleAssoc.txt")
        gz <- FALSE
    } else {
        return(NULL)
    }
    return(list(gene_file=gene_file, variant_file=variant_file, gz=gz))
}

checks <- function(file_info, file_info_template, no_sex_check=FALSE)
{
	sex_out <- ifelse(no_sex_check, warning, stop)
	
    if (file_info$phenotype != file_info_template$phenotype) {
    	cat(paste(paste(names(file_info), file_info, sep=": "), collapse="\n"), "\n")
        stop(paste0("phenotype ", file_info$phenotype, " does not match '",
        	file_info_template$phenotype, "' - check files or rename"))
    }
    if (file_info$sex != file_info_template$sex) {
        cat(paste(paste(names(file_info), file_info, sep=": "), collapse="\n"), "\n")
        sex_out(paste0("sex ", file_info$sex, " does not match '",
        	file_info_template$sex, "' - check files or rename"))
    }
    if (file_info$type != file_info_template$type) {
        cat(paste(paste(names(file_info), file_info, sep=": "), collapse="\n"), "\n")
        stop("attempting to meta-analyse gene based results with variant results - check files")
    }
    if (file_info$binary != file_info_template$binary) {
        cat(paste(paste(names(file_info), file_info, sep=": "), collapse="\n"), "\n")
        stop("attempting to analyse binary phenotype with a continuous phenotype")
    }
}

extract_file_info <- function(filename)
{
	gz <- ifelse(grepl(".gz$", filename), TRUE, FALSE)
	filename <- gsub(".gz$", "", filename)
	filename <- gsub("cleaned.", "", filename)
	file_info <- as.list(strsplit(filename, split="\\.")[[1]])

	if (!(length(file_info) %in% c(12, 13))) {
		print(file_info)
		stop("Incorrect file naming convention, please check filenames")
	}

	if (length(file_info) == 13) {
		binary <- TRUE
		names(file_info) <- c(
			"dataset",
			"last_name",
			"analysis_name",
			"phenotype",
			"freeze_number",
			"sex",
			"ancestry",
			"n_cases",
			"n_controls",
			"software",
			"type",
			"date",
			"split")
	} else {
		binary <- FALSE
		names(file_info) <- c(
			"dataset",
			"last_name",
			"analysis_name",
			"phenotype",
			"freeze_number",
			"sex",
			"ancestry",
			"n",
			"software",
			"type",
			"date",
			"split")
	}

	file_info$gz <- gz
	file_info$binary <- binary
	return(file_info)
}

add_N_using_filename <- function(file_info, dt)
{
    if (!file_info$binary) {
        dt$N_eff <- file_info$n
    } else {
        N_case <- file_info$n_cases
        N_control <- file_info$n_controls
        N_eff <- (4 * N_case * N_control) / (N_case + N_control)
        dt$N_eff <- N_eff
        dt$N_case <- N_case
        dt$N_control <- N_control
    }
    return(dt)
}

add_N_using_Neff_weights_file <- function(file_info,
	dt,
	Neff_weights_file="/well/lindgren/dpalmer/BRaVa_meta-analysis_inputs/Neff/Neff_weights_may25.tsv.gz") {
	if (!file.exists(Neff_weights_file)) {
		warning("Neff weights files not found")
		return(dt)
	}
	dt_weights <- fread(Neff_weights_file)
	# Replace with Neff if it is present, otherwise throw a warning.
	Neff_replace <- (dt_weights %>% filter(
		pheno == file_info$phenotype,
		ancestry == file_info$ancestry,
		dataset == file_info$dataset
		))$nglmm
	if (length(Neff_replace) == 1) {
		dt$N_eff <- Neff_replace
	} else if (length(Neff_replace) == 0) {
		warning("Nglmm is not present, reverting to assuming all samples are unrelated")
	} else {
		warning("Multiple matches to this (pheno, ancestry, dataset) tuple")
	}
	return(dt)
}

add_N <- function(file_info, dt_gene)
{
    # Read in the variant file to extract the sample size
    dt_variant <- fread(cmd = ifelse(file_info$gz,
        paste0("gzcat ", file_info$variant_file),
        paste0("cat ", file_info$variant_file)))

    if ("N" %in% colnames(dt_variant)) {
        N_eff <- mean(dt_variant$N)
        dt_gene$N_eff <- N_eff
        binary <- FALSE
    } else {
        N_case <- mean(dt_variant$N_case)
        N_control <- mean(dt_variant$N_ctrl)
        N_eff <- (4 * N_case * N_control) / (N_case + N_control)
        dt_gene$N_eff <- N_eff
        dt_gene$N_case <- N_case
        dt_gene$N_control <- N_control
        binary <- TRUE
    }
    return(list(dt_gene = dt_gene, binary=binary))
}

determine_null_correlation <- function(dt, pval_T=0.05, binary=TRUE)
{
	# Assume that correlation across ancestry labels is 0
	cor_count_pval <- list()
	# Restrict to what should be the null set of tests
	if (binary) {
		dt <- dt %>% filter(Group == "synonymous", max_MAF == 1e-3, MAC_case > 30, MAC_case < 1000)
	} else {
		dt <- dt %>% filter(Group == "synonymous", max_MAF == 1e-3, MAC > 50)
	}

	i <- 1
	for (anc in c("AFR", "AMR", "EAS", "EUR", "SAS")) {
		dt_tmp <- dt %>% filter(ancestry == anc)
		datasets <- unique(dt_tmp$dataset)
		for (d1 in datasets) {
			dt_i <- dt_tmp %>% filter(dataset == d1)
			setkeyv(dt_i, c("max_MAF", "Region"))
			for (d2 in datasets) {
				dt_j <- dt_tmp %>% filter(dataset == d2)
				setkeyv(dt_j, c("max_MAF", "Region"))
				dt_ij <- merge(dt_i, dt_j)
				if (nrow(dt_ij) < 100) {
					next
				}
				Z_i <- qnorm(dt_ij$`Pvalue.x`/2) * sign(dt_ij$`BETA_Burden.x`)
				Z_j <- qnorm(dt_ij$`Pvalue.y`/2) * sign(dt_ij$`BETA_Burden.y`)
				test <- cor.test(Z_i, Z_j, alternative = "greater")
				cor_count_pval[[i]] <- data.table(
					dataset=d1,
					dataset2=d2,
					ancestry=anc,
					cor=test$estimate,
					pval=test$`p.value`,
					count=nrow(dt_ij),
					Neff1=dt_i$N_eff[1],
					Neff2=dt_j$N_eff[2]
				) %>% mutate(summation = sign(cor)*sqrt(abs(cor)*Neff1)*sqrt(abs(cor)*Neff2)*(pval<pval_T))
				i <- i+1
			}
		}
	}
	cor_count_pval <- rbindlist(cor_count_pval)
	setkeyv(cor_count_pval, c("dataset", "dataset2", "ancestry"))
	return(cor_count_pval)
}


cauchy_combination <- function(p_values, weights=NULL)
{
	is.zero <- sum(p_values == 0) >= 1
	is.one <- sum(p_values > (1 - 1e-14)) >= 1

	if (is.zero) {
		return(0)
	}

	if (is.one) {
		p <- min(p_values) * length(p_values)
		if (p > 1) {
			return(-Inf)
		} else {
			return(qcauchy(p, lower.tail=FALSE))
		}
		return(min(1, (min(p_values)) * (length(p_values))))
	}

	if (is.null(weights)) {
		weights <- rep(1 / length(p_values), length(p_values))
	} else if (length(weights) != length(p_values)) {
		stop("The length of weights should be the same as that of the p-values!")
	} else if (sum(weights < 0) > 0) {
		stop("All the weights must be positive!")
	} else {
		weights <- weights / sum(weights)
	}

	is_small <- (p_values < 1e-16)
	if (sum(is_small) == 0) {
		cct_stat <- sum(weights * tan((0.5 - p_values) * pi))
	} else {
		cct_stat <- sum((weights[is_small] / p_values[is_small]) / pi)
		cct_stat <- cct_stat + 
			sum(weights[!is_small] * tan((0.5 - p_values[!is_small]) * pi))
	}

	return(cct_stat)
}

weighted_fisher <- function(p_values, weights=NULL, two_tail=FALSE, input_beta=NULL)
{
	if (is.null(weights)) {
		weights <- rep(1, length(p_values))
	}

	idx.na <- which(is.na(p_values))
	
	if (length(idx.na) > 0) {
		p_values <- p_values[-idx.na]
		weights <- weights[-idx.na]
		if (two_tail) {
			input_beta <- input_beta[-idx.na]
		}
	}

	NP <- length(p_values)
	NS <- length(weights)
	if (NP != NS) { stop("The length of p and weights vector must be identical.") }

	N <- NS
	Ntotal <- sum(weights)
	ratio <- weights / Ntotal
	Ns <- N * ratio
	G <- c()

	if (!two_tail) {
		for (i in 1:length(p_values)) {
		  G <- append(G, qgamma(p = p_values[i], shape = Ns[i], scale = 2, lower.tail=FALSE))
		}
		Gsum <- sum(G)
		resultP <- pgamma(q = Gsum, shape = N, scale = 2, lower.tail = FALSE)
	} else {
		p1 <- p2 <- p_values
		idx_pos <- which(input_beta > 0)
		idx_neg <- which(input_beta < 0)
		
		# Positive direction
		G <- c()
		p1[idx_pos] <- p_values[idx_pos] / 2
		p1[idx_neg] <- 1 - p_values[idx_neg] / 2

		for (i in 1:length(p1)) {
		  G <- append(G, qgamma(p = p1[i], shape = Ns[i], scale = 2, lower.tail = FALSE))
		}
		Gsum <- sum(G)
		resultP1 <- pgamma(q = Gsum, shape = N, scale = 2, lower.tail = FALSE)
		
		# Negative direction
		G <- c()
		p2[idx_pos] <- 1 - p_values[idx_pos] / 2
		p2[idx_neg] <- p_values[idx_neg] / 2

		for (i in 1:length(p2)) {
		  G <- append(G, qgamma(p = p2[i], shape = Ns[i], scale = 2, lower.tail = FALSE))
		}
		Gsum <- sum(G)
		resultP2 <- pgamma(q = Gsum, shape = N, scale = 2, lower.tail = FALSE)
		resultP <- 2 * min(resultP1, resultP2)
		if (resultP > 1.0) {
			resultP <- 1.0
		}
	}
	return(min(1, resultP))
}

run_weighted_fisher <- function(
	grouped_dt, n_eff_name, input_pvalues, output_meta_pvalue,
	two_tail = FALSE, input_beta = NULL)
{
	if (two_tail) {
		result <- grouped_dt %>% 
		summarise(
			"{output_meta_pvalue}" := weighted_fisher(
				.data[[input_pvalues]],
				weights = if (is.null(n_eff_name)) { NULL } else { .data[[n_eff_name]] },
				two_tail = two_tail, input_beta=.data[[input_beta]]
			)
		)
	} else {
		result <- grouped_dt %>% 
		summarise(
			"{output_meta_pvalue}" := weighted_fisher(
				.data[[input_pvalues]],
				weights = if (is.null(n_eff_name)) { NULL } else { .data[[n_eff_name]] }
			)
		)
	}
	return(result)
}

run_inv_var <- function(
	grouped_dt, input_beta_name, input_se_name, 
	output_beta_meta, output_se_meta,
	output_meta_pvalue
) {
	dt <- grouped_dt %>% 
		mutate(weight = 1/(.data[[input_se_name]]**2)) %>%
		mutate(effs_inv_var = .data[[input_beta_name]] * weight) %>%
		summarise(
			"{output_beta_meta}" := sum(effs_inv_var) / sum(weight),
			"{output_se_meta}" := sqrt(1/sum(weight))) %>% 
		mutate("{output_meta_pvalue}" := 2 * pnorm(
				abs(.data[[output_beta_meta]] / .data[[output_se_meta]]), lower=FALSE))
	
	# Deal with the edge cases
	# Send the pvalues of the data with infinite standard errors to 1, and the 
	# pvalues of the data with standard errors of 0 to 0.
	dt$Pvalue[which(is.na(dt$Pvalue) & (dt$SE_Burden == Inf))] <- 1
	dt$Pvalue[which(is.na(dt$Pvalue) & (dt$SE_Burden == 0))] <- 0

	return(dt)
}

create_gene_data_table <- function(group_files_lines)
{
	group_files_variant_line <- group_files_lines[1]
	group_files_group_line <- group_files_lines[2]
	to_munge_variant <- strsplit(group_files_variant_line, split=" ")[[1]]
	n <- length(to_munge_variant)
	to_munge_group <- strsplit(group_files_group_line, split=" ")[[1]]
	if (n != length(to_munge_group)) {
		stop("Error: number of variants in the region does not match the number of groups in the region")
	}
	return(data.table(Region=to_munge_variant[1],
		MarkerID=to_munge_variant[3:n],
		Group=to_munge_group[3:n]))
}

extract_group_file_information <- function(regexp_group_files) {
	# Extract the directory portion and the regular expression within the folder
	regexp_group_files <- strsplit(regexp_group_files, "/")[[1]]
	if (length(regexp_group_files) > 1) {
		directory <- paste(regexp_group_files[1:(length(regexp_group_files)-1)], collapse="/")
		regexp_group_files <- regexp_group_files[length(regexp_group_files)]
	} else {
		directory <- "."
	}
	group_files <- lapply(dir(directory, pattern=regexp_group_files), fread, sep="@", header=FALSE)
	group_files <- rbindlist(group_files)
	print(group_files)
	gene_data_table_list <- list()
	total_lines <- nrow(group_files)
	i <- 0
	for (i in 1:(nrow(group_files)/2)) {
		print(i)
		gene_data_table_list[[i]] <- create_gene_data_table(group_files$V1[(2*i-1):(2*i)])
	}
	gene_data_table <- rbindlist(gene_data_table_list)
	return(gene_data_table)
}

combine_gene_and_variant_information <- function(variant_file, gene_data_table)
{
	# Merge with the variant files to obtain the variants that have a MAC>10 in the analysed dataset
	dt_variant <- fread(variant_file)
	# Deal with the ultra-rare variants
	setkey(dt_variant, "MarkerID")
	setkey(gene_data_table, "MarkerID")

	dt <- merge(dt_variant, gene_data_table, all=TRUE)
	dt[, Region:=ifelse(CHR=="UR", gsub(":.*", "", MarkerID), Region)]
	dt[, Group:=ifelse(CHR=="UR", gsub(".*:(.*):.*", "\\1", MarkerID), Group)]
	dt[, Tstat := ifelse(AF_Allele2 > 0.5, -Tstat, Tstat)]
	dt[, MAF := pmin(AF_Allele2, 1-AF_Allele2)]
	return(dt)
}

extract_combined_annotations <- function(gene_file) {
	dt_gene <- fread("data/meta_analysis/gcloud/uk-biobank.palmer.PRELIMINARY.Coronary_artery_disease.JULY23Freeze.ALL.EUR.19915.382460.SAIGE.gene.20240110.cleaned.txt.gz")
	combined_annotations <- grep(";", unique(dt_gene$Group), value=TRUE)
	return(combined_annotations)
}

extract_beta_burden <- function(dt, combined_annotations, max_MAF)
{
	dt <- dt %>% filter(!is.na(Region))
	m <- max_MAF
	dt_common <- dt %>% 
		filter(CHR != "UR", MAF < max_MAF) %>% 
		group_by(Region, Group) %>% 
		summarise(
			numerator_common = sum(Tstat),
			denominator_common = sum(var))
	dt_rare <- dt %>% filter(CHR == "UR") %>% 
		mutate(max_MAF = as.numeric(gsub(".*:", "", MarkerID))) %>% 
		filter(max_MAF == m) %>% group_by(Region, Group) %>% 
		summarise(
			numerator_rare = sum(Tstat),
			denominator_rare = sum(var))
	dt_common <- data.table(dt_common)
	
	# Finally, if there are any combined annotations, then
	# include those

	# The below just needs to be done for dt_common, not for dt_rare.
	dt_common_combined_list <- list()
	for (str_comb_annot in combined_annotations) {
		combined_annotation_vector <- strsplit(str_comb_annot, split=";")[[1]]
		dt_common_combined_list[[str_comb_annot]] <- dt_common %>% group_by(Region) %>% 
			filter(Group %in% combined_annotation_vector) %>%
			summarise(
				numerator_common = sum(numerator_common),
				denominator_common = sum(denominator_common)
				) %>% mutate(Group = str_comb_annot)
	}
	rbind(rbindlist(dt_common_combined_list), dt_common)

	setkeyv(dt_common, c("Region", "Group"))
	dt_rare <- data.table(dt_rare)
	setkeyv(dt_rare, c("Region", "Group"))
	dt <- merge(dt_common, dt_rare, all=TRUE)
	
	# Numerator
	dt[, numerator_common := ifelse(is.na(numerator_common), 0, numerator_common)]
	dt[, numerator_rare := ifelse(is.na(numerator_rare), 0, numerator_rare)]

	# Denominator
	dt[, denominator_common := ifelse(is.na(denominator_common), 0, denominator_common)]
	dt[, denominator_rare := ifelse(is.na(denominator_rare), 0, denominator_rare)]
	return(dt)
}

extract_beta_burden_across_MAFs <- function(dt, max_MAFs, combined_annotations) {
	# Loop over max MAFs
	dt_list <- list()
	# unique(dt_gene$max_MAF[!is.na(dt_gene$max_MAF)])
	for (m in max_MAFs) {
		dt_list[[as.character(m)]] <- extract_beta_burden(dt, combined_annotations, m)
		dt_list[[as.character(m)]]$max_MAF <- m
	}
	dt_out <- rbindlist(dt_list) %>% group_by(Region, Group, max_MAF)
}

run_heterogeneity <- function(
	grouped_dt, n_eff_name, input_beta, output_meta_beta)
{
	grouped_dt <- grouped_dt %>%
	mutate(
		weights = sqrt(.data[[n_eff_name]]),
		beta = .data[[input_beta]]
	)
	summary_dt <- grouped_dt %>% 
	summarise(
		sum_weights = sum(weights),
		sum_betas = sum(weights * beta)
		) %>% mutate("{output_meta_beta}" := sum_betas/sum_weights)
	return(
		merge(grouped_dt, summary_dt) %>% 
		mutate(
			deviation = weights * (beta - .data[[output_meta_beta]])^2
		) %>% group_by(across(as.character(groups(grouped_dt)))) %>%
		summarise(sum_deviation = sum(deviation), n_studies=n()) %>% 
		filter(n_studies > 1) %>% mutate(p_het = pchisq(sum_deviation, n_studies-1, lower.tail=FALSE))
	)
}

weights <- function(dt, is_inv_var=TRUE,
	n_eff_name=NULL, se_name=NULL) {
	if (is_inv_var) {
		return(dt %>% mutate(weights = 1/(.data[[se_name]]^2)))
	} else {
		return(dt %>% mutate(weights = sqrt(.data[[n_eff_name]])))
	}
}

run_heterogeneity_test <- function(
	grouped_dt, input_beta, output_meta_beta) {
	grouped_dt %>% 
	mutate(weighted_effects = weights * .data[[input_beta]]) %>% 
	summarise(
		df = n()-1,
		sum_weights = sum(weights),
		"{output_meta_beta}" := sum(weighted_effects) / sum_weights,
		chisq_het = sum(
			weights * (.data[[input_beta]] - .data[[output_meta_beta]])^2),
		Pvalue_het = pchisq(chisq_het, df, lower.tail=FALSE)
	)
}

run_fisher <- function(
	grouped_dt, chi2_stat_name, input_pvalues, output_meta_pvalue)
{
	return(
		grouped_dt %>% 
		summarise(df = 2*n(),
			"{chi2_stat_name}" := -2*sum(log(.data[[input_pvalues]]))) %>%
		mutate("{output_meta_pvalue}" := pchisq(.data[[chi2_stat_name]],
			df=df, lower.tail=FALSE)) %>% select(-df)
	)
}

run_cauchy <- function(
	grouped_dt, n_eff_name, Cauchy_stat_name, input_pvalues, output_meta_pvalue) 
{
	return(
		grouped_dt %>% 
		summarise(
			"{Cauchy_stat_name}" := cauchy_combination(
				.data[[input_pvalues]],
				weights=sqrt(.data[[n_eff_name]])
			),
			number_of_pvals := n()
		) %>% mutate(
		  	"{output_meta_pvalue}" := ifelse(
				.data[[Cauchy_stat_name]] > 1e+15,
				(1 / .data[[Cauchy_stat_name]]) / pi,
				pcauchy(.data[[Cauchy_stat_name]], lower.tail=FALSE)
			)
		) %>% mutate(
			"{output_meta_pvalue}" := ifelse(
				.data[[output_meta_pvalue]] > (1 - 1e-10),
				(1 - 1/number_of_pvals), .data[[output_meta_pvalue]]
			)
		)
	)
}

run_stouffer <- function(
	grouped_dt, n_eff_name, weighted_Z_name,
	input_pvalues, output_meta_pvalue,
	two_tail = FALSE, input_beta = NULL
) {
	if (two_tail) {
		grouped_dt <- grouped_dt %>% 
			mutate("{input_pvalues}" := .data[[input_pvalues]]/2)
	} else {
		input_beta <- "beta_dummy"
		grouped_dt <- grouped_dt %>% mutate("{input_beta}" := 1)
	}
	
	result <- grouped_dt %>%
		mutate(
			weighted_Z_numerator = (
				sqrt(.data[[n_eff_name]]) * 
				(-qnorm(.data[[input_pvalues]])) * 
				sign(.data[[input_beta]])
			)
		)

	if (two_tail) {
		result <- result %>%
		summarise(
			"{weighted_Z_name}" := sum(weighted_Z_numerator) / 
				sqrt(sum(.data[[n_eff_name]])),
			"{output_meta_pvalue}" := 2 * pnorm(abs(.data[[weighted_Z_name]]), lower.tail=FALSE)
		)
	} else {
		result <- result %>%
		summarise(
			"{weighted_Z_name}" := sum(weighted_Z_numerator) / 
				sqrt(sum(.data[[n_eff_name]])),
			"{output_meta_pvalue}" := pnorm(.data[[weighted_Z_name]], lower.tail=FALSE)
		)
	}
	return(result)
}


run_stouffer_overlap <- function(
	grouped_dt, n_eff_name, weighted_Z_name, cross_terms,
	input_pvalues, output_meta_pvalue,
	two_tail = FALSE, input_beta = NULL
) {
	if (two_tail) {
		grouped_dt <- grouped_dt %>% 
			mutate("{input_pvalues}" := .data[[input_pvalues]]/2)
	} else {
		input_beta <- "beta_dummy"
		grouped_dt <- grouped_dt %>% mutate("{input_beta}" := 1)
	}
	
	result <- grouped_dt %>%
		mutate(
			weighted_Z_numerator = (
				sqrt(.data[[n_eff_name]]) * 
				(-qnorm(.data[[input_pvalues]])) * 
				sign(.data[[input_beta]])
			)
		)

	if (two_tail) {
		result <- result %>%
		summarise(
			"{weighted_Z_name}" := sum(weighted_Z_numerator) / 
				sqrt(sum(.data[[cross_terms]])),
			"{output_meta_pvalue}" := 2 * pnorm(abs(.data[[weighted_Z_name]]), lower.tail=FALSE)
		)
	} else {
		result <- result %>%
		summarise(
			"{weighted_Z_name}" := sum(weighted_Z_numerator) / 
				sqrt(sum(.data[[cross_terms]])),
			"{output_meta_pvalue}" := pnorm(.data[[weighted_Z_name]], lower.tail=FALSE)
		)
	}
	return(result)
}