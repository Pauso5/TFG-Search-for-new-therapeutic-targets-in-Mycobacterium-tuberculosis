# ============================================================
#  PIPELINE Dianes2 (R) — Resistència confirmada en literatura (global) — FIGL-only
#  Entrada:
#    - Dianes2/Dades/Resistencia_soques2.xlsx        (taula resistències confirmades + columna RESISTENCIA)
#    - Dianes2/Dades/sinonims_soques.csv             (mapping manual de soques)
#    - Dianes2/Dades/FIGL/                           (02/03/... amb N.figl)
#  Sortida:
#    - Dianes2/results2/ (CSVs amb associacions gen ↔ resistència global)
#
#  Objectiu:
#    1) Parsejar FIGL -> ortho_id × soca (presència/absència)
#    2) Llegir resistència confirmada i construir fenotip global (resistent vs sensible)
#    3) Test Fisher per cada ortho_id contra fenotip global
# ============================================================

# -------------------------
# 0) Packages (paquets utilitzats per l'anàlisi)
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
# 1) RUTES (EDITA NOMÉS path_projecte)
# -------------------------
path_projecte <- "C:/Users/pauss/Desktop/Apunts/TFG/Dianes2"

# Marcar els paths
path_dades <- file.path(path_projecte, "Dades")
path_figl  <- file.path(path_dades, "FIGL")

# Resultats en una carpeta nova per no barrejar amb els anteriors
path_out   <- file.path(path_projecte, "results2")
if (!dir.exists(path_out)) dir.create(path_out, recursive = TRUE)

fitxer_resistencia <- file.path(path_dades, "Resistencia_soques2.xlsx")
fitxer_sinonims    <- file.path(path_dades, "sinonims_soques.csv")


# -------------------------
# 2) FUNCIONS
# -------------------------
clean_names_base <- function(nms) {
  # Treu accents (RESISTÈNCIA -> RESISTENCIA) i normalitza
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
         "\nCrea o copia 'sinonims_soques.csv' dins de Dades/.")
  }
  
  # Prova CSV amb coma; si no, prova amb ';'
  syn <- tryCatch(readr::read_csv(file, show_col_types = FALSE), error = function(e) NULL)
  if (is.null(syn)) {
    syn <- readr::read_delim(file, delim = ";", show_col_types = FALSE)
  }
  
  names(syn) <- clean_names_base(names(syn))
  
  # 'nota' és opcional; només exigim from/to
  if (!all(c("from", "to") %in% names(syn))) {
    stop("El CSV de sinònims ha de tenir columnes: from,to (nota és opcional).")
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

# -------- FIGL --------
get_size_from_path <- function(f) suppressWarnings(as.integer(basename(dirname(f))))
get_index_from_filename <- function(f) suppressWarnings(as.integer(tools::file_path_sans_ext(basename(f))))

parse_figl_file <- function(f, syn) {
  size <- get_size_from_path(f)
  idx  <- get_index_from_filename(f)
  
  lines <- readLines(f, warn = FALSE, encoding = "UTF-8")
  lines <- lines[nzchar(trimws(lines))]
  
  purrr::map_df(seq_along(lines), function(i) {
    line <- lines[i]
    
    # separa score final (espai/tab)
    m <- str_match(line, "^(.*?)[\\t ]+([0-9]+\\.?[0-9]*)\\s*$")
    left  <- m[,2]
    score <- suppressWarnings(as.numeric(m[,3]))
    if (is.na(left)) return(NULL)
    
    # elimina prefix tipus '12:'
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
    
    # -----------------------------
    # IMPORTANT: ortho_id UNIFORME (no depèn d'H37Rv)
    # -----------------------------
    ortho_id <- paste0("ORTHO__size", size, "__", basename(f), "__line", i)
    
    # Opcional: guardar el gen de H37Rv en una columna separada (per facilitar interpretació/cerca)
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
        ref_gene = ref_gene,     # <-- NOVA COLUMNA (opcional però recomanada)
        strain_clean = strain_clean,
        gene_in_strain = gene_raw
      )
  })
}

# -------------------------
# 3) CARREGAR SINÒNIMS
# -------------------------
syn <- read_synonyms_robust(fitxer_sinonims)

# -------------------------
# 4) LLEGIR RESISTÈNCIA CONFIRMADA (EXCEL)
# -------------------------
if (!file.exists(fitxer_resistencia)) stop("No trobo: ", fitxer_resistencia)

res2_raw <- readxl::read_excel(fitxer_resistencia)
names(res2_raw) <- clean_names_base(names(res2_raw))

# Treu columnes buides "unnamed" o numèriques
res2_raw <- res2_raw %>% select(-matches("^unnamed"), -matches("^[0-9]+$"))

# Esperem una columna de soca i una columna global RESISTENCIA
# A la teva taula: 'organism_infraspecific_names_strain' i 'resistencia'
if (!("organism_infraspecific_names_strain" %in% names(res2_raw))) {
  stop("No trobo la columna 'Organism Infraspecific Names Strain' a l'Excel.")
}
if (!("resistencia" %in% names(res2_raw))) {
  stop("No trobo la columna global 'RESISTENCIA' a l'Excel. Assegura't que existeix i s'escriu igual.")
}

res2 <- res2_raw %>%
  rename(
    strain = organism_infraspecific_names_strain,
    resistencia_global = resistencia
  ) %>%
  mutate(
    strain_clean_raw = clean_strain(strain),
    strain_clean     = apply_synonyms(strain_clean_raw, syn),
    resistencia_global = toupper(trimws(as.character(resistencia_global)))
  )

# Normalitza valors globals: SENSIBLE vs RESISTENT
# Accepta: SENSIBLE / MDR / XDR / R
res2 <- res2 %>%
  mutate(
    resistant = case_when(
      resistencia_global %in% c("MDR","XDR","R","RESISTENT","RESISTENTE","RESISTANT") ~ TRUE,
      resistencia_global %in% c("SENSIBLE","SENSITIVE","SUSCEPTIBLE","S") ~ FALSE,
      TRUE ~ NA
    )
  )

# Informe ràpid de valors desconeguts
unknown_vals <- res2 %>% filter(is.na(resistant) & !is.na(resistencia_global)) %>% distinct(resistencia_global)
if (nrow(unknown_vals) > 0) {
  warning("Valors de RESISTENCIA no reconeguts (s'ignoren com NA): ", paste(unknown_vals$resistencia_global, collapse=", "))
}

# Filtra només soques amb fenotip global definit
res2_f <- res2 %>% filter(!is.na(strain_clean), !is.na(resistant))
write.csv(res2_f, file.path(path_out, "resistencia_confirmada_global_neta.csv"), row.names = FALSE)

# -------------------------
# 5) (FIGL-only) — NO llegim IGs ni fem validació IGs vs FIGL
# -------------------------
igs_long <- tibble()
message("Mode FIGL-only: s'ignoren fitxers IGs i no es fa validació IGs vs FIGL.")

# -------------------------
# 6) LLEGIR I PARSEJAR FIGL
# -------------------------
figl_files <- list.files(path_figl, pattern = "\\.figl$", full.names = TRUE, recursive = TRUE)
if (length(figl_files) == 0) stop("No he trobat .figl a: ", path_figl)

figl_long <- map_df(figl_files, parse_figl_file, syn = syn)
write.csv(figl_long, file.path(path_out, "figl_long_parsed.csv"), row.names = FALSE)

# -------------------------
# 6b) DEBUG (FIGL-only): coherència i possibles barreges de runs
# -------------------------

# (1) Detecta noms duplicats de .figl (pot indicar runs barrejats si recursive=TRUE)
dup_figl_names <- tibble(path = figl_files) %>%
  mutate(name = basename(path)) %>%
  count(name, sort = TRUE) %>%
  filter(n > 1)

write.csv(dup_figl_names, file.path(path_out, "DEBUG_duplicate_figl_filenames.csv"), row.names = FALSE)

if (nrow(dup_figl_names) > 0) {
  warning("⚠️ Hi ha fitxers .figl amb el mateix nom en rutes diferents. Mira DEBUG_duplicate_figl_filenames.csv (possible barreja de runs).")
}

# (2) Soques que apareixen en FIGL però NO tenen fenotip usable
debug_figl_no_pheno <- figl_long %>%
  distinct(strain_clean) %>%
  anti_join(res2_f %>% distinct(strain_clean), by = "strain_clean") %>%
  arrange(strain_clean)

write.csv(debug_figl_no_pheno, file.path(path_out, "DEBUG_figl_strains_without_pheno.csv"), row.names = FALSE)

# (3) Soques amb fenotip però que NO surten a cap FIGL
debug_pheno_no_figl <- res2_f %>%
  distinct(strain_clean) %>%
  anti_join(figl_long %>% distinct(strain_clean), by = "strain_clean") %>%
  arrange(strain_clean)

write.csv(debug_pheno_no_figl, file.path(path_out, "DEBUG_pheno_strains_not_in_figl.csv"), row.names = FALSE)

# -------------------------
# 7) PRESÈNCIA/ABSÈNCIA ORTHO × SOCA
# -------------------------
ortho_strain <- figl_long %>%
  filter(!is.na(ortho_id), !is.na(strain_clean), strain_clean != "") %>%
  distinct(strain_clean, ortho_id)

# Limita a soques amb fenotip global
ortho_strain <- ortho_strain %>% semi_join(res2_f %>% select(strain_clean) %>% distinct(), by = "strain_clean")
write.csv(ortho_strain, file.path(path_out, "ortho_strain_presence.csv"), row.names = FALSE)

strains_by_ortho <- ortho_strain %>%
  group_by(ortho_id) %>%
  summarise(
    strains_present = list(unique(strain_clean)),
    n_present = n_distinct(strain_clean),
    .groups = "drop"
  )

# -------------------------
# 8) TEST FISHER: ortho_id ↔ RESISTÈNCIA GLOBAL
# -------------------------
min_present <- 3
min_resistant <- 3

strains_by_ortho2 <- strains_by_ortho %>% filter(n_present >= min_present)

# Fenotip
ph <- res2_f %>% select(strain_clean, resistant)
universe <- ph$strain_clean
ph_vec <- ph$resistant
names(ph_vec) <- ph$strain_clean

if (sum(ph_vec) < min_resistant) {
  warning("Hi ha pocs resistents (", sum(ph_vec), ") per fer tests robustos.")
}

assoc_global <- purrr::map_df(seq_len(nrow(strains_by_ortho2)), function(i) {
  o <- strains_by_ortho2$ortho_id[i]
  present_strains <- strains_by_ortho2$strains_present[[i]]
  
  present <- universe %in% present_strains
  resistant <- ph_vec[universe]
  
  tab <- table(present, resistant)
  ft <- tryCatch(fisher.test(tab), error = function(e) NULL)
  if (is.null(ft)) return(NULL)
  
  tibble(
    phenotype = "GLOBAL_RESISTANCE",
    ortho_id = o,
    odds_ratio = unname(ft$estimate),
    p_value = ft$p.value,
    n_strains_tested = length(universe),
    n_present = sum(present),
    n_resistant = sum(resistant)
  )
})

assoc_global <- assoc_global %>%
  mutate(q_value = p.adjust(p_value, method = "BH")) %>%
  arrange(q_value, p_value)

write.csv(assoc_global, file.path(path_out, "associacions_ortho_global_fisher.csv"), row.names = FALSE)

# TOP candidats (exploratori)
top_global <- assoc_global %>%
  filter(q_value <= 0.1, (is.infinite(odds_ratio) | odds_ratio >= 2)) %>%
  arrange(q_value)
write.csv(top_global, file.path(path_out, "TOP_candidats_global.csv"), row.names = FALSE)

# -------------------------
# 9) MAPA DE LOCUS TAGS PER SOCA (interpretació dels candidats)
# -------------------------
orthomap <- figl_long %>%
  filter(!is.na(ortho_id), !is.na(strain_clean), !is.na(gene_in_strain), gene_in_strain != "") %>%
  distinct(ortho_id, strain_clean, gene_in_strain)

if (nrow(top_global) > 0) {
  cand_map <- orthomap %>%
    semi_join(top_global %>% select(ortho_id) %>% distinct(), by = "ortho_id") %>%
    arrange(ortho_id, strain_clean)
  write.csv(cand_map, file.path(path_out, "TOP_candidats_global_locus_tags.csv"), row.names = FALSE)
}

cat("\n✅ OK (Dianes2 — FIGL-only). Resultats a: ", normalizePath(path_out), "\n")
cat("Fitxer principal: associacions_ortho_global_fisher.csv\n")

View(top_global)

# NO DONA CAP CANDIDAT LA CONFIANÇA MARCADA
# BUSQUEM ELS MILLORS RESULTATS ENCARA QUE NO COMPLEIXIN ELS PARÀMETRES

# -------------------------
# 10) BÚSQUEM LES MILLORS OPCIONS
# -------------------------
assoc_global <- assoc_global %>%
  mutate(q_value = p.adjust(p_value, method = "BH")) %>%
  arrange(q_value, p_value)

# Guarda tots els resultats (complet)
write.csv(assoc_global, file.path(path_out, "associacions_ortho_global_fisher.csv"), row.names = FALSE)

# 1) TOP per q_value (FDR) i per p_value (sense correcció)
top_by_q <- assoc_global %>% arrange(q_value, p_value) %>% slice_head(n = 50)
top_by_p <- assoc_global %>% arrange(p_value, q_value) %>% slice_head(n = 50)

write.csv(top_by_q, file.path(path_out, "TOP50_by_qvalue.csv"), row.names = FALSE)
write.csv(top_by_p, file.path(path_out, "TOP50_by_pvalue.csv"), row.names = FALSE)

# 2) TOP per OR (tenint en compte que Inf podria dominar)
top_by_or <- assoc_global %>%
  arrange(desc(is.infinite(odds_ratio)), desc(odds_ratio), q_value) %>%
  slice_head(n = 50)

write.csv(top_by_or, file.path(path_out, "TOP50_by_oddsratio.csv"), row.names = FALSE)

# 3) Llistes per llindars (estricte / bo / exploratori)
strict  <- assoc_global %>% filter(q_value <= 0.05)
good    <- assoc_global %>% filter(q_value <= 0.10)
explore <- assoc_global %>% filter(q_value <= 0.20)

write.csv(strict,  file.path(path_out, "CANDIDATS_q005.csv"), row.names = FALSE)
write.csv(good,    file.path(path_out, "CANDIDATS_q010.csv"), row.names = FALSE)
write.csv(explore, file.path(path_out, "CANDIDATS_q020.csv"), row.names = FALSE)

# 4) Score combinat (p + OR) per rankejar millor
calc_or_corrected <- function(present_strains, universe, resistant_vec) {
  present <- universe %in% present_strains
  resistant <- resistant_vec[universe]
  
  a <- sum(present & resistant, na.rm = TRUE)
  b <- sum(present & !resistant, na.rm = TRUE)
  c <- sum(!present & resistant, na.rm = TRUE)
  d <- sum(!present & !resistant, na.rm = TRUE)
  
  ( (a + 0.5) * (d + 0.5) ) / ( (b + 0.5) * (c + 0.5) )
}

ortho_map <- strains_by_ortho2 %>% select(ortho_id, strains_present)

assoc_scored <- assoc_global %>%
  left_join(ortho_map, by = "ortho_id") %>%
  mutate(
    odds_ratio_corr = map_dbl(strains_present, calc_or_corrected, universe = universe, resistant_vec = ph_vec),
    log2_or_corr = log2(odds_ratio_corr),
    neglog10_p = -log10(p_value + 1e-300),
    score = neglog10_p * pmax(log2_or_corr, 0)
  ) %>%
  arrange(desc(score), q_value, p_value)

assoc_scored_fixed <- assoc_scored %>%
  mutate(
    strains_present = map_chr(strains_present, ~ paste(.x, collapse = ";"))
  )

write.csv(assoc_scored_fixed, file.path(path_out, "associacions_scored.csv"), row.names = FALSE)

top_score <- assoc_scored %>% slice_head(n = 50)
top_score_clean <- top_score %>% select(-strains_present)

write.csv(top_score_clean, file.path(path_out, "TOP50_by_score.csv"), row.names = FALSE)

# 5) “TOP candidates” relaxat però amb OR mínim
top_like_before <- assoc_global %>%
  filter(q_value <= 0.20, (is.infinite(odds_ratio) | odds_ratio >= 2)) %>%
  arrange(q_value, p_value)

write.csv(top_like_before, file.path(path_out, "TOP_candidats_global_relaxats.csv"), row.names = FALSE)

View(top_score_clean)


# -------------------------
# 10) HEATMAP
# -------------------------

library(tidyverse)
library(readxl)
library(pheatmap)
library(Matrix)

# 1) Carrega presència per ortho_id i soca (FIGL parsed)
figl <- read.csv("C:/Users/pauss/Desktop/Apunts/TFG/Dianes2/results2/figl_long_parsed.csv", stringsAsFactors = FALSE)

# IMPORTANT: si tu filtres per fenotip al pipeline, i vols exactament el mateix univers:
# carrega resistència neta que has guardat
res2_f <- read.csv("C:/Users/pauss/Desktop/Apunts/TFG/Dianes2/results2/resistencia_confirmada_global_neta.csv", stringsAsFactors = FALSE)
phenotyped_strains <- unique(res2_f$strain_clean)

# Filtra figl només a soques amb fenotip (mateix univers que Fisher)
figl <- figl %>%
  filter(strain_clean %in% phenotyped_strains) %>%
  distinct(ortho_id, strain_clean)

# 2) Construeix matriu esparsa (estalvia memòria)
ortho_levels  <- sort(unique(figl$ortho_id))
strain_levels <- sort(unique(figl$strain_clean))

i <- match(figl$ortho_id, ortho_levels)
j <- match(figl$strain_clean, strain_levels)

M <- sparseMatrix(i = i, j = j, x = 1,
                  dims = c(length(ortho_levels), length(strain_levels)),
                  dimnames = list(ortho_levels, strain_levels))

# 3) Anotació de columnes (fenotip)
anno_col <- data.frame(
  Fenotip = res2_f$resistant[match(colnames(M), res2_f$strain_clean)],
  row.names = colnames(M)
)
anno_col$Fenotip <- ifelse(anno_col$Fenotip, "Resistent", "Sensible")

# 4) Dibuix (recomanat: sense clustering de files si hi ha milers)
pheatmap(
  as.matrix(M),                      # si peta memòria, mira Opció B (ComplexHeatmap)
  color = c("white", "black"),
  cluster_rows = FALSE,              # TRUE pot ser molt lent amb milers
  cluster_cols = TRUE,
  show_rownames = FALSE,             # amb molts ortho_id, millor amagar
  show_colnames = TRUE,
  annotation_col = anno_col,
  fontsize_col = 7,
  main = "GLOBAL_RESISTANCE — Presència/absència (tots els ortho_id)"
)


#ALTRES TIPUS DE GRÀFICS:

#VOLCANO GRÀFIC:

library(tidyverse)
library(ggplot2)

sc <- read.csv("C:/Users/pauss/Desktop/Apunts/TFG/Dianes2/results2/associacions_scored.csv") %>%
  filter(phenotype == "GLOBAL_RESISTANCE") %>%
  mutate(neglog10p = -log10(p_value + 1e-300))

ggplot(sc, aes(x = log2_or_corr, y = neglog10p)) +
  geom_point(aes(size = n_present, color = score), alpha = 0.7) +
  
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  geom_vline(xintercept = c(-1, 1), linetype = "dotted", color = "grey40") +
  
  theme_minimal() +
  labs(title = "Associació IGs — resistència soques",
       x = "log2(OR corregida)",
       y = "-log10(p-value)",
       color = "Score",
       size = "n_present")

#BARPLOT:
topN <- sc %>%
  arrange(p_value) %>%
  distinct(ortho_id, .keep_all = TRUE) %>%
  slice_head(n = 30)

ggplot(topN, aes(x = reorder(ortho_id, neglog10p), y = neglog10p, fill = log2_or_corr > 0)) +
  geom_col() +
  coord_flip() +
  theme_minimal() +
  labs(title="Top 30 ortho_id per p-value (GLOBAL_RESISTANCE)",
       x="ortho_id", y="-log10(p-value)", fill="log2OR>0")
