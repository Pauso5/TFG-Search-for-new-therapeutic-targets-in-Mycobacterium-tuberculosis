# ============================================================
#  PIPELINE Dianes1 (R) — Resistència predita per BV-BCR (Antibiòtics)
#  FIGL-only + Resistència (xlsx) + Sinònims
#  + Associació (Fisher) ortho ↔ resistència per antibiòtic i/o MDR
#
#  Estructura:
#  projecte/
#   ├─ Dades/
#   │   ├─ Resistencia_soques.xlsx
#   │   ├─ sinonims_soques.csv
#   │   └─ FIGL/    (02/, 03/, 04/, ... amb 1.figl,2.figl,...)
#   └─ results/
#
#  Notes:
#   - IGs eliminat: aquest pipeline NO depèn d'IGs ni valida IGs vs FIGL.
#   - ortho_id és SEMPRE uniforme i NO depèn d'H37Rv.
#   - Si H37Rv hi és, guardem el seu gen a 'ref_gene' (opcional, per buscar).
# ============================================================

# -------------------------
# 0) Packages
# -------------------------
suppressPackageStartupMessages({
  library(readxl)
  library(readr)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(purrr)
  library(tibble)
})

# -------------------------
# 1) RUTES
# -------------------------
path_projecte <- "C:/Users/pauss/Desktop/Apunts/TFG/Dianes"

# IMPORTANT: el teu cas és "Dades" amb majúscula
path_dades <- file.path(path_projecte, "Dades")

# subcarpetes
path_figl  <- file.path(path_dades, "FIGL")

# resultats
path_out   <- file.path(path_projecte, "results")
if (!dir.exists(path_out)) dir.create(path_out, recursive = TRUE)

# fitxers
fitxer_resistencia <- file.path(path_dades, "Resistencia_soques.xlsx")
fitxer_sinonims    <- file.path(path_dades, "sinonims_soques.csv")

# Comprovar camins (opcional)
cat("path_dades:", file.exists(path_dades), "\n")
cat("path_figl:", file.exists(path_figl), "\n")
cat("fitxer_resistencia:", file.exists(fitxer_resistencia), "\n")
cat("fitxer_sinonims:", file.exists(fitxer_sinonims), "\n")

# -------------------------
# 2) FUNCIONS AUXILIARS
# -------------------------

clean_names_base <- function(nms) {
  # Normalitza noms de columnes
  nms <- iconv(nms, from = "UTF-8", to = "ASCII//TRANSLIT")
  nms %>%
    tolower() %>%
    gsub("[^a-z0-9]+", "_", .) %>%
    gsub("^_|_$", "", .)
}

clean_strain <- function(x) {
  x %>%
    as.character() %>%
    toupper() %>%
    str_replace_all("^MTU", "") %>%
    str_replace_all("MTU", "") %>%
    str_replace_all("[^A-Z0-9]", "") %>%
    na_if("")
}

read_synonyms_robust <- function(file) {
  if (!file.exists(file)) {
    stop("No trobo el fitxer de sinònims: ", file,
         "\nCrea Dades/sinonims_soques.csv amb columnes: from,to (nota opcional).")
  }
  
  # Prova CSV amb coma; si no, prova amb ';'
  syn <- tryCatch(readr::read_csv(file, show_col_types = FALSE), error = function(e) NULL)
  if (is.null(syn)) {
    syn <- readr::read_delim(file, delim = ";", show_col_types = FALSE)
  }
  
  names(syn) <- clean_names_base(names(syn))
  
  # 'nota' és opcional; from/to obligatori
  if (!all(c("from", "to") %in% names(syn))) {
    stop("El CSV de sinònims ha de tenir columnes: from,to (nota opcional).")
  }
  
  syn %>%
    mutate(
      from = clean_strain(from),
      to   = clean_strain(to)
    ) %>%
    filter(!is.na(from), !is.na(to))
}

apply_synonyms <- function(x, syn) {
  idx <- match(x, syn$from)
  out <- x
  out[!is.na(idx)] <- syn$to[idx[!is.na(idx)]]
  out
}

# ---------- FIGL ----------
get_size_from_path <- function(f) suppressWarnings(as.integer(basename(dirname(f))))
get_index_from_filename <- function(f) suppressWarnings(as.integer(tools::file_path_sans_ext(basename(f))))

# Parser robust per línies .figl:
# - Score és l'últim número a la línia, separat per espais (o tab).
# - Part esquerra conté tokens separats per '+', cada token és 'Strain:GeneId'.
# - Pot aparèixer '#' enganxat al gen (s'elimina).
# - ortho_id: SEMPRE uniforme i NO depèn d'H37Rv.
# - ref_gene: si hi ha H37Rv, guardem el seu gen aquí (opcional).
parse_figl_file <- function(f, syn) {
  size <- get_size_from_path(f)
  idx  <- get_index_from_filename(f)
  
  lines <- readLines(f, warn = FALSE, encoding = "UTF-8")
  lines <- lines[nzchar(trimws(lines))]
  
  purrr::map_df(seq_along(lines), function(i) {
    line <- lines[i]
    
    # 1) separar score (últim bloc numèric)
    m <- str_match(line, "^(.*?)[\\t ]+([0-9]+\\.?[0-9]*)\\s*$")
    left  <- m[,2]
    score <- suppressWarnings(as.numeric(m[,3]))
    if (is.na(left)) return(NULL)
    
    # 2) elimina prefix del tipus '12:'
    left <- str_replace(left, "^\\s*[0-9]+:\\s*", "")
    
    tokens <- str_split(left, "\\+")[[1]]
    tokens <- tokens[nzchar(tokens)]
    
    df_tok <- purrr::map_df(tokens, function(tk) {
      pieces <- str_split(tk, ":", n = 2)[[1]]
      if (length(pieces) < 2) return(NULL)
      tibble(
        strain_raw = pieces[1],
        gene_raw   = str_replace_all(pieces[2], "#", "")
      )
    })
    
    if (nrow(df_tok) == 0) return(NULL)
    
    df_tok <- df_tok %>%
      mutate(
        strain_clean_raw = clean_strain(strain_raw),
        strain_clean     = apply_synonyms(strain_clean_raw, syn)
      )
    
    # ortho_id UNIFORME (mai depèn d'H37Rv)
    ortho_id <- paste0("ORTHO__size", size, "__", basename(f), "__line", i)
    
    # ref_gene opcional (si hi ha H37Rv a la línia)
    h37 <- df_tok %>% filter(str_detect(toupper(strain_raw), "H37RV$"))
    ref_gene <- if (nrow(h37) > 0 && !is.na(h37$gene_raw[1]) && h37$gene_raw[1] != "") {
      h37$gene_raw[1]
    } else {
      NA_character_
    }
    
    df_tok %>%
      transmute(
        figl_file = basename(f),
        ig_size = size,
        ig_index = idx,
        line_no = i,
        score = score,
        ortho_id = ortho_id,
        ref_gene = ref_gene,
        strain_clean = strain_clean,
        gene_in_strain = gene_raw
      )
  })
}

# -------------------------
# 3) LLEGIR SINÒNIMS
# -------------------------
syn <- read_synonyms_robust(fitxer_sinonims)
write.csv(syn, file.path(path_out, "sinonims_nets.csv"), row.names = FALSE)

# -------------------------
# 4) LLEGIR RESISTÈNCIA (Excel)
# -------------------------
if (!file.exists(fitxer_resistencia)) stop("No trobo: ", fitxer_resistencia)

res_raw <- readxl::read_excel(fitxer_resistencia)
names(res_raw) <- clean_names_base(names(res_raw))
res_raw <- res_raw %>% select(-matches("^unnamed"))

# Comprovacions mínimes
if (!("organism_infraspecific_names_strain" %in% names(res_raw))) {
  stop("No trobo la columna 'Organism Infraspecific Names Strain' a l'Excel.")
}

# Alguns Excels porten 'linage' o 'lineage' o 'llinatge'... (farem robust)
lineage_col <- dplyr::case_when(
  "linage" %in% names(res_raw)  ~ "linage",
  "lineage" %in% names(res_raw) ~ "lineage",
  "llinatge" %in% names(res_raw) ~ "llinatge",
  TRUE ~ NA_character_
)

# Resistència score global: 'resistencia' segons el teu script
if (!("resistencia" %in% names(res_raw))) {
  warning("No trobo la columna 'resistencia' (score global). Continuo sense MDR_high.")
}

res <- res_raw %>%
  mutate(
    strain = organism_infraspecific_names_strain,
    strain_clean_raw = clean_strain(strain),
    strain_clean     = apply_synonyms(strain_clean_raw, syn)
  )

if (!is.na(lineage_col)) {
  res <- res %>% mutate(lineage = .data[[lineage_col]])
} else {
  res <- res %>% mutate(lineage = NA)
}

if ("resistencia" %in% names(res_raw)) {
  res <- res %>% mutate(resistencia_score = .data[["resistencia"]])
} else {
  res <- res %>% mutate(resistencia_score = NA)
}

# Columnes meta i columnes d'antibiòtic
cols_meta <- c("organism_name","organism_infraspecific_names_strain","strain","lineage",
               "resistencia","resistencia_score","strain_clean_raw","strain_clean")
ab_cols <- setdiff(names(res), cols_meta)

# Normalitza valors R/S (per si hi ha minúscules o espais)
res <- res %>%
  mutate(across(all_of(ab_cols), ~ toupper(trimws(as.character(.x)))))

# Comptes globals R/S per fila (opcional)
res <- res %>%
  mutate(
    n_R = rowSums(across(all_of(ab_cols), ~ .x == "R"), na.rm = TRUE),
    n_S = rowSums(across(all_of(ab_cols), ~ .x == "S"), na.rm = TRUE)
  )

write.csv(res, file.path(path_out, "resistencia_neta_amb_sinonims.csv"), row.names = FALSE)

# -------------------------
# 5) FIGL-only: LLEGIR I PARSEJAR FIGL
# -------------------------
figl_files <- list.files(path_figl, pattern = "\\.figl$", full.names = TRUE, recursive = TRUE)
if (length(figl_files) == 0) stop("No he trobat .figl a: ", path_figl)

# DEBUG duplicats reals: mateix ig_size i mateix filename repetit
dup_figl_names <- tibble(path = figl_files) %>%
  mutate(
    ig_size = suppressWarnings(as.integer(basename(dirname(path)))),
    name    = basename(path)
  ) %>%
  count(ig_size, name, sort = TRUE) %>%
  filter(n > 1)

write.csv(dup_figl_names,
          file.path(path_out, "DEBUG_duplicate_figl_filenames_by_size.csv"),
          row.names = FALSE)

if (nrow(dup_figl_names) > 0) {
  warning("⚠️ Duplicats reals detectats (mateix ig_size + mateix .figl). Mira DEBUG_duplicate_figl_filenames_by_size.csv")
}

figl_long <- map_df(figl_files, parse_figl_file, syn = syn)
write.csv(figl_long, file.path(path_out, "figl_long_parsed.csv"), row.names = FALSE)

# DEBUG soques FIGL sense fenotip (R/S) en cap antibiòtic? (opcional informatiu)
debug_figl_not_in_res <- figl_long %>%
  distinct(strain_clean) %>%
  anti_join(res %>% distinct(strain_clean), by = "strain_clean") %>%
  arrange(strain_clean)

write.csv(debug_figl_not_in_res, file.path(path_out, "DEBUG_figl_strains_not_in_excel.csv"), row.names = FALSE)

# -------------------------
# 6) Presència/absència ortho_id per soca (només soques amb dades a Excel)
# -------------------------
ortho_strain <- figl_long %>%
  filter(!is.na(ortho_id), !is.na(strain_clean), strain_clean != "") %>%
  distinct(strain_clean, ortho_id) %>%
  semi_join(res %>% select(strain_clean) %>% distinct(), by = "strain_clean")

write.csv(ortho_strain, file.path(path_out, "ortho_strain_presence.csv"), row.names = FALSE)

strains_by_ortho <- ortho_strain %>%
  group_by(ortho_id) %>%
  summarise(
    strains_present = list(unique(strain_clean)),
    n_present = n_distinct(strain_clean),
    .groups = "drop"
  )

# -------------------------
# 7) Associació (Fisher) + SCORING per antibiòtic
# -------------------------

min_present <- 3
min_R <- 3

strains_by_ortho2 <- strains_by_ortho %>%
  filter(n_present >= min_present)

run_fisher_for_ab <- function(ab) {
  
  ph <- res %>%
    select(strain_clean, all_of(ab)) %>%
    filter(!is.na(strain_clean), .data[[ab]] %in% c("R", "S")) %>%
    mutate(isR = (.data[[ab]] == "R"))
  
  if (nrow(ph) == 0) return(NULL)
  if (sum(ph$isR) < min_R) return(NULL)
  
  universe <- ph$strain_clean
  ph_vec <- ph$isR
  names(ph_vec) <- ph$strain_clean
  
  purrr::map_df(seq_len(nrow(strains_by_ortho2)), function(i) {
    
    o <- strains_by_ortho2$ortho_id[i]
    present_strains <- strains_by_ortho2$strains_present[[i]]
    
    present <- universe %in% present_strains
    resistant <- ph_vec[universe]
    
    tab <- table(present, resistant)
    
    ft <- tryCatch(fisher.test(tab), error = function(e) NULL)
    if (is.null(ft)) return(NULL)
    
    tibble(
      antibiotic = ab,
      ortho_id = o,
      odds_ratio = unname(ft$estimate),
      p_value = ft$p.value,
      n_strains_tested = length(universe),
      n_present = sum(present),
      n_R = sum(resistant)
    )
  })
}

# -------- Fisher per tots antibiòtics --------
assoc_all <- purrr::map_df(ab_cols, run_fisher_for_ab)

if (nrow(assoc_all) == 0) {
  warning("No hi ha resultats d'associació")
} else {
  
  # -------- Correcció FDR --------
  assoc_all <- assoc_all %>%
    mutate(q_value = p.adjust(p_value, method = "BH"))
  
  # -------- FUNCIO OR CORREGIDA (pseudocount 0.5) --------
  compute_or_corr <- function(ortho_id, antibiotic) {
    
    ph <- res %>%
      select(strain_clean, all_of(antibiotic)) %>%
      filter(!is.na(strain_clean), .data[[antibiotic]] %in% c("R", "S")) %>%
      mutate(isR = (.data[[antibiotic]] == "R"))
    
    universe <- ph$strain_clean
    ph_vec <- ph$isR
    names(ph_vec) <- ph$strain_clean
    
    present_strains <- strains_by_ortho$strains_present[
      match(ortho_id, strains_by_ortho$ortho_id)
    ][[1]]
    
    present <- universe %in% present_strains
    resistant <- ph_vec[universe]
    
    tab <- table(present, resistant)
    
    # forcem 2x2 sempre
    tab2 <- matrix(0, nrow = 2, ncol = 2)
    rownames(tab2) <- c("FALSE", "TRUE")
    colnames(tab2) <- c("FALSE", "TRUE")
    
    tab2[rownames(tab), colnames(tab)] <- tab
    
    a <- tab2["TRUE","TRUE"]
    b <- tab2["TRUE","FALSE"]
    c <- tab2["FALSE","TRUE"]
    d <- tab2["FALSE","FALSE"]
    
    ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))
  }
  
  # -------- SCORING COMPLET --------
  assoc_scored <- assoc_all %>%
    rowwise() %>%
    mutate(
      odds_ratio_corr = compute_or_corr(ortho_id, antibiotic),
      log2_or_corr    = log2(odds_ratio_corr),
      neglog10_p      = -log10(p_value + 1e-300),
      score_sym       = neglog10_p * abs(log2_or_corr)
    ) %>%
    ungroup() %>%
    arrange(q_value, desc(score_sym))
  
  # -------- Guardar resultats --------
  write.csv(
    assoc_all,
    file.path(path_out, "associacions_ortho_antibiotic_fisher.csv"),
    row.names = FALSE
  )
  
  write.csv(
    assoc_scored,
    file.path(path_out, "associacions_antibiotics_scored.csv"),
    row.names = FALSE
  )
}


# -------------------------
# 8) Rànquing global (multi-antibiòtic)
# -------------------------
global_hits <- assoc_all %>%
  mutate(sig = !is.na(q_value) & q_value <= 0.05 & odds_ratio > 1) %>%
  group_by(ortho_id) %>%
  summarise(
    n_ab_tested = n(),
    n_ab_sig = sum(sig),
    best_q = suppressWarnings(min(q_value, na.rm = TRUE)),
    best_or = odds_ratio[which.min(q_value)],
    antibiotics_sig = paste(unique(antibiotic[sig]), collapse = ";"),
    .groups = "drop"
  ) %>%
  arrange(desc(n_ab_sig), best_q)

write.csv(global_hits, file.path(path_out, "ranking_global_orthos.csv"), row.names = FALSE)

candidats <- global_hits %>% filter(n_ab_sig >= 2, best_q <= 0.05)
write.csv(candidats, file.path(path_out, "candidats_dianes.csv"), row.names = FALSE)

# -------------------------
# 9) Mapa locus tags per soca per als candidats
# -------------------------
orthomap <- figl_long %>%
  filter(!is.na(ortho_id), !is.na(strain_clean), !is.na(gene_in_strain), gene_in_strain != "") %>%
  distinct(ortho_id, strain_clean, gene_in_strain)

if (nrow(candidats) > 0) {
  cand_map <- orthomap %>%
    semi_join(candidats %>% select(ortho_id), by = "ortho_id") %>%
    arrange(ortho_id, strain_clean)
  write.csv(cand_map, file.path(path_out, "candidats_locus_tags_per_soca.csv"), row.names = FALSE)
}

# -------------------------
# 10) MDR binari a partir de resistencia_score (si existeix)
# -------------------------
mdr_threshold <- 10

if ("resistencia_score" %in% names(res)) {
  res_mdr <- res %>%
    filter(!is.na(resistencia_score), !is.na(strain_clean)) %>%
    mutate(MDR_high = as.numeric(resistencia_score) >= mdr_threshold) %>%
    select(strain_clean, MDR_high)
  
  universe <- res_mdr$strain_clean
  ph_vec <- res_mdr$MDR_high
  names(ph_vec) <- res_mdr$strain_clean
  
  out_mdr <- purrr::map_df(seq_len(nrow(strains_by_ortho2)), function(i) {
    o <- strains_by_ortho2$ortho_id[i]
    present_strains <- strains_by_ortho2$strains_present[[i]]
    
    present <- universe %in% present_strains
    high <- ph_vec[universe]
    
    tab <- table(present, high)
    ft <- tryCatch(fisher.test(tab), error = function(e) NULL)
    if (is.null(ft)) return(NULL)
    
    tibble(
      phenotype = paste0("MDR_high_ge_", mdr_threshold),
      ortho_id = o,
      odds_ratio = unname(ft$estimate),
      p_value = ft$p.value,
      n_strains_tested = length(universe),
      n_present = sum(present),
      n_high = sum(high)
    )
  })
  
  if (nrow(out_mdr) > 0) {
    out_mdr <- out_mdr %>%
      mutate(q_value = p.adjust(p_value, method = "BH")) %>%
      arrange(q_value, p_value)
    write.csv(out_mdr, file.path(path_out, "associacions_ortho_MDR_fisher.csv"), row.names = FALSE)
  }
}

cat("\n✅ OK (Dianes1 FIGL-only). Resultats a: ", normalizePath(path_out), "\n")
cat("Fitxer principal: associacions_ortho_antibiotic_fisher.csv\n")

# Resums ràpids
if (nrow(assoc_all) > 0) {
  print(summary(assoc_all$q_value))
  print(table(assoc_all$q_value < 0.1))
  print(summary(assoc_all$odds_ratio))
}

# -------------------------
# 11) BÚSQUEM LES MILLORS OPCIONS (exploratori)
# -------------------------
top_hits1 <- assoc_all %>%
  filter(q_value <= 0.1, odds_ratio > 1) %>%
  arrange(q_value, p_value) %>%
  slice_head(n = 20)

top_hits2 <- assoc_all %>%
  filter(odds_ratio >= 2, q_value <= 0.2) %>%
  arrange(q_value, p_value) %>%
  slice_head(n = 20)

top_inf <- assoc_all %>%
  filter(is.infinite(odds_ratio)) %>%
  arrange(q_value, p_value) %>%
  slice_head(n = 20)

write.csv(top_hits1, file.path(path_out, "TOP20_q010_or_gt1.csv"), row.names = FALSE)
write.csv(top_hits2, file.path(path_out, "TOP20_or_ge2_q020.csv"), row.names = FALSE)
write.csv(top_inf,  file.path(path_out, "TOP20_or_infinite.csv"), row.names = FALSE)

# -------------------------
# 12) GRÀFICS:
# -------------------------

#VOLCANO:

library(ggplot2)

ggplot(assoc_scored, aes(x = log2_or_corr, y = neglog10_p)) +
  geom_point(
    aes(size = n_present, color = score_sym),
    alpha = 0.75
  ) +
  
  # línies de referència
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  geom_vline(xintercept = c(-1, 1), linetype = "dotted", color = "grey40") +
  
  theme_minimal() +
  labs(
    title = "Associació IGs — resistència predita a antibiòtics",
    x = "log2(OR corregida)",
    y = "-log10(p-value)",
    color = "score_sym",
    size = "n_present"
  )

#VOLCANO NET
library(ggrepel)

top_labels <- assoc_scored %>%
  arrange(desc(score_sym)) %>%
  slice_head(n = 15)

ggplot(assoc_scored, aes(x = log2_or_corr, y = neglog10_p)) +
  geom_point(aes(size = n_present), alpha = 0.5) +
  geom_text_repel(data = top_labels, aes(label = ortho_id), size = 3) +
  theme_minimal()


#HEATMAP

library(Matrix)

# 1) Carrega figl_long (si no el tens en memòria)
figl <- read.csv("C:/Users/pauss/Desktop/Apunts/TFG/Dianes/results/figl_long_parsed.csv", stringsAsFactors = FALSE)

# 2) Agafa només ortho_id i strains vàlides
figl <- figl %>%
  filter(!is.na(ortho_id), !is.na(strain_clean), strain_clean != "") %>%
  distinct(ortho_id, strain_clean)

# 3) Definir nivells
ortho_levels  <- sort(unique(figl$ortho_id))
strain_levels <- sort(unique(figl$strain_clean))

# 4) Índexs
i <- match(figl$ortho_id, ortho_levels)
j <- match(figl$strain_clean, strain_levels)

# 5) Matriu presència/absència
M <- sparseMatrix(
  i = i,
  j = j,
  x = 1,
  dims = c(length(ortho_levels), length(strain_levels)),
  dimnames = list(ortho_levels, strain_levels)
)

TOP_N <- 30

top_hits <- assoc_scored %>%
  arrange(desc(score_sym)) %>%
  distinct(ortho_id, .keep_all = TRUE) %>%
  slice_head(n = TOP_N)

M_zoom <- M[rownames(M) %in% top_hits$ortho_id, , drop = FALSE]

library(pheatmap)


#Crear anotació de columnes (fenotip)
anno_col <- data.frame(
  Fenotip = res$isoniazid_inh[match(colnames(M_zoom), res$strain_clean)]
)
rownames(anno_col) <- colnames(M_zoom)

# Heatmap
pheatmap(
  as.matrix(M_zoom),
  color = c("white", "black"),
  breaks = c(-0.5, 0.5, 1.5),
  legend_breaks = c(0, 1),
  legend_labels = c("Absent", "Present"),
  annotation_col = anno_col,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  main = "Top ortho_id (score_sym)"
)


