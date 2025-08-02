db_config_array = [
    {"db_fasta": "IPAH/ipaH_150-mers.fasta", "db_marker": "ipaH", "identity": 98, "coverage": 95, "rfb": None},
    {"db_fasta": "RFB/RFB_serotypes_AtoC_150-mers_v2.fasta", "db_marker": "rfb", "identity": 98, "coverage": 95, "rfb": None},
    {"db_fasta": "RFB/galF_SB1.fasta", "db_marker": "additionalrfb", "identity": 98, "coverage": 95, "rfb": "C1"},
    {"db_fasta": "RFB/RFB_serogroup_D_150-mers.fasta", "db_marker": "additionalrfb", "identity": 100, "coverage": 100, "rfb": None},
    {"db_fasta": "RFB/POAC-genes_150-mers.fasta", "db_marker": "POAC", "identity": 98, "coverage": 95, "rfb": None},
    {"db_fasta": "FLIC/fliC_Shigella_v1.fasta", "db_marker": "flic", "identity": 98, "coverage": 95, "rfb": None},
    {"db_fasta": "CRISPR/CRISPR_spacers.fasta", "db_marker": "crispr", "identity": 100, "coverage": 100, "rfb": None}
]

db_alt_config_array = [
    {"db_fasta": "RFB/RFB_AprovBEDP02-5104_150-mers.fasta", "db_marker": "additionalrfb", "identity": 98, "coverage": 95, "rfb": "A2", "new_rfb": "AprovBEDP02-5104"},
    {"db_fasta": "RFB/RFB_A16_150-mers_v2.fasta", "db_marker": "additionalrfb", "identity": 98, "coverage": 95, "rfb": "A3a", "new_rfb": "A16"},
    {"db_fasta": "RFB/RFB_A16_150-mers_v2.fasta", "db_marker": "additionalrfb", "identity": 98, "coverage": 95, "rfb": "A3b", "new_rfb": "A16"},
    {"db_fasta": "RFB/taurine_SB6.fasta", "db_marker": "additionalrfb", "identity": 98, "coverage": 95, "rfb": "C10", "new_rfb": "C6"},
]

db_mlst_config_array = [
    {"db_fasta": "MLST/adk_len.fasta", "db_marker": "adk", "identity": 100, "coverage": 100},
    {"db_fasta": "MLST/fumC_len.fasta", "db_marker": "fumC", "identity": 100, "coverage": 100},
    {"db_fasta": "MLST/gyrB_len.fasta", "db_marker": "gyrB", "identity": 100, "coverage": 100},
    {"db_fasta": "MLST/icd_len.fasta", "db_marker": "icd", "identity": 100, "coverage": 100},
    {"db_fasta": "MLST/mdh_len.fasta", "db_marker": "mdh", "identity": 100, "coverage": 100},
    {"db_fasta": "MLST/purA_len.fasta", "db_marker": "purA", "identity": 100, "coverage": 100},
    {"db_fasta": "MLST/recA_len.fasta", "db_marker": "recA", "identity": 100, "coverage": 100}
]

flex_score_array = [
    {"flex": "1a", "score": 1},
    {"flex": "1b", "score": 129},
    {"flex": "1c(7a)", "score": 3},
    {"flex": "7b", "score": 131},
    {"flex": "2a", "score": 4},
    {"flex": "2b", "score": 36},
    {"flex": "3a", "score": 96},
    {"flex": "3b", "score": 64},
    {"flex": "3b atypical (oac1b)", "score": 128},
    {"flex": "4a", "score": 8},
    {"flex": "4av", "score": 264},
    {"flex": "4b", "score": 72},
    {"flex": "4bv", "score": 328},
    {"flex": "5a", "score": 16},
    {"flex": "5a", "score": 80},
    {"flex": "5b", "score": 48},
    {"flex": "5b", "score": 112},
    {"flex": "X", "score": 32},
    {"flex": "Xv", "score": 288},
    {"flex": "Y", "score": 0},
    {"flex": "Yv", "score": 256}
]
