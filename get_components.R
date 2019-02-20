#!/usr/bin/env Rscript

# Check if basic packages are installed -----------------------------------

is.installed <- function(pkg){
  is.element(pkg, installed.packages()[,1])
}

if (!is.installed("crayon") || !is.installed("optparse") || !is.installed("config")){
  cat("We will try to install the packages crayon, optparse and config... (this will be only be done once)\n")
  Sys.sleep(5)
  if (!is.installed("crayon")){
    suppressMessages(install.packages("crayon", repos = "http://cran.us.r-project.org"))
  }
  if (!is.installed("optparse")){
    suppressMessages(install.packages("optparse", repos = "http://cran.us.r-project.org"))
  }
  if (!is.installed("config")){
    suppressMessages(install.packages("config", repos = "http://cran.us.r-project.org"))
  }
}

suppressMessages(library(crayon))
suppressMessages(library(optparse))


# Check if packaged are installed -----------------------------------------

cat("\nChecking if all packages are installed...\n\n")


needed <- c("igraph", "tidyverse", "pbmcapply", "maditr", "data.table",
            "tidygraph", "bigreadr", "unixtools")

missing_package <- FALSE
# For loop to run through each of the packages
for (p in 1:length(needed)){
  if(is.installed(needed[p])){
    cat(sprintf("%-10s: %s", needed[p], green("Installed\n")))
  }else{
    cat(sprintf("%-10s: %s", needed[p], red("Not installed\n")))
    missing_package <- TRUE
  }
}

quit_not_installed <- function(){
  cat("\nMissing packages, please install them.\n")
  quit(save = "no", status = 1)
}

if (missing_package) {
  quit_not_installed()
}else{
  cat("\nAll packages installed.\n")
}

Sys.sleep(2)
system("clear")



# Load libraries ----------------------------------------------------------
cat("Loading libraries...")
silent <- suppressMessages(lapply(needed, function(X) {require(X, character.only = TRUE)}))
rm(silent)
cat(" done\n")

# Get config file

# Script command line options ---------------------------------------------

option_list <- list(
  make_option(c("-c", "--config"), type="character", default=NULL,
              help="YAML config file", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list);
opt <- parse_args(opt_parser);


if (is.null(opt$config)){
  print_help(opt_parser)
  stop("At least one arguments must be supplied (YAML config file).n", call.=FALSE)
}

cat("Getting config parameters...")
cfg <- config::get(file = opt$config, use_parent = FALSE)
cat(" done\n")

tmp <- file.path(cfg$wd, "tmp")
cat(paste0("Creating TMP folder at ", tmp, "..."))
if (dir.exists(tmp)) unlink(tmp, recursive = TRUE, force = TRUE)
dir.create(tmp, recursive = TRUE)
# Set a new teppdir
set.tempdir(tmp)
cat("done\n")

# Create a new environment
lo_env <- new.env()

# Define data.table parameters
dt_nthreads <- cfg$dt_cores
options(datatable.optimize=Inf)
options(datatable.verbose = cfg$verbose)
# setDTthreads(dt_nthreads)
# getDTthreads()


# library(foreach)
# library(doSNOW)
# library(itertools)
# nproc <- 64
# cl <- makeCluster(nproc, rscript_args = c("--no-init-file", "--no-site-file", "--no-environ"), outfile = "/scratch/antonio/unk_C_SC/cl_out.log")
# registerDoSNOW(cl)
# getDoParWorkers()
# getDoParRegistered()
# getDoParName()

# Read data ---------------------------------------------------------------

#load("/bioinf/projects/megx/UNKNOWNS/chiara/architect/components/knowns_DA.rda",
#    verbose = TRUE)


# Number of k clusters
# 1,122,105
cat("Reading cluster categories...")
cl_cat <- fread(file = cfg$cl_cat, header = FALSE, showProgress = TRUE)
setnames(cl_cat, names(cl_cat), c("cl_name", "category"))
cat(" done\n")

cat("Reading cluster completeness...")
cl_compl <- fread(file = cfg$cl_compl, header = TRUE, showProgress = TRUE)
cat(" done\n")

cat("Reading PFAM31 clan information...")
p_clan <- fread(file = cfg$p_clan, header = FALSE, showProgress = TRUE) %>%
  dt_select(V2, V4) %>%
  setnames(c("V2", "V4"), c("clan", "pfam"))
cat(" done\n")

cat("Reading simplified domain architecture information...")
p_doms <- fread(input = cfg$p_doms, header = TRUE, showProgress = TRUE) %>%
  mutate(cl_name = as.character(cl_name))
cat(" done\n")

k_comp <- p_doms %>% dt_filter(cl_name %in% (cl_cat %>% dt_filter(category == "K") %>% .$cl_name))


lo_env$k_hhblits_all <- big_fread1("/scratch/antonio/unk_C_SC/data/k_hhblits.tsv",
                                   every_nlines = 500e6,
                                   data.table = TRUE,
                                   .combine = list,
                                   nThread = 78,
                                   verbose = TRUE,
                                   header = FALSE
)

k_hhblits_all_nrow <- map(unlist(lo_env$k_hhblits_all, recursive = FALSE), nrow) %>% unlist() %>% sum()
cat(paste("Read", scales::comma(k_hhblits_all_nrow), "lines\n"))

# Let's filter for prob >= 50, and cov >= 0.6
lo_env$k_hhblits <-  data.table::rbindlist(map(unlist(lo_env$k_hhblits_all, recursive = FALSE), function(X){
  X <- X %>% dt_filter(V1 != V2, V3 >= 50, V13 > 0.6, V14 > 0.6)
  setnames(X, names(X),
           c("cl_name1", "cl_name2", "probability", "e-value", "Score",
             "Cols", "q_start", "q_stop", "t_start", "t_stop", "q_len",
             "t_len", "q_cov", "t_cov"))
  X %>% maditr::dt_mutate(cl_name1 = as.character(cl_name1),
                          cl_name2 = as.character(cl_name2),
                          score_col = Score/Cols)
}))

# How many clusters do we have (removed self-hits)
# 1,122,105
# 1,100,297 p50;c0.6
unique_hhblits_k_cl <- c(lo_env$k_hhblits$cl_name1, lo_env$k_hhblits$cl_name2) %>% unique()
length(unique_hhblits_k_cl)

# We are missing 21,808
# Which ones
setdiff(cl_cat %>% dt_filter(category == "K") %>% .$cl_name, unique_hhblits_k_cl) %>% length()

# Find missing clusters ---------------------------------------------------
missing_ids <- setdiff(cl_cat %>% dt_filter(category == "K") %>% .$cl_name, unique_hhblits_k_cl)

lo_env$k_hhblits_missing <-  data.table::rbindlist(map(unlist(lo_env$k_hhblits_all, recursive = FALSE), function(X){
  X <- X %>% dt_filter(V1 %in% missing_ids | V2 %in% missing_ids) %>%
    dt_filter(V1 != V2)
  setnames(X, names(X),
           c("cl_name1", "cl_name2", "probability", "e-value", "Score",
             "Cols", "q_start", "q_stop", "t_start", "t_stop", "q_len",
             "t_len", "q_cov", "t_cov"))
  X %>% maditr::dt_mutate(cl_name1 = as.character(cl_name1),
                          cl_name2 = as.character(cl_name2),
                          score_col = Score/Cols)
}))

# Prepare data for network analysis ---------------------------------------

lo_env$k_hhb_bh_score <- lo_env$k_hhblits %>% dt_select(cl_name1, cl_name2, score_col) %>% as_tibble()

k_hh_g <- lo_env$k_hhblits %>%
  dt_select(cl_name1, cl_name2, score_col) %>%
  as_tibble() %>%
  rename(weight = score_col) %>%
  igraph::graph_from_data_frame(directed = FALSE) %>%
  igraph::simplify(remove.multiple = TRUE, remove.loops = TRUE, edge.attr.comb = list("max")) %>%
  as_tbl_graph()

vcount(k_hh_g)
vcount(k_hh_g %>% filter(node_is_isolated()))

k_hh_g <- k_hh_g %>%
  activate(nodes) %>%
  #select(-archit) %>%
  inner_join(k_comp %>% as_tibble() %>% rename(name = cl_name) %>% mutate(name = as.character(name))) %>%
  filter(!node_is_isolated())



# Identify communities using MCL ------------------------------------------
source("/scratch/antonio/unk_C_SC/graph_lib.R")

gx <- k_hh_g

lo_env$w <- E(gx)$weight
lo_env$w  <- lo_env$w - min(lo_env$w ) + 0.001

E(gx)$weight <- lo_env$w

inflation_list <- seq(1.2, 3, 0.1)
mcl_bin <- "mcl"
k_g_mcl_list <- pbmcapply::pbmclapply(inflation_list, optimal_mcl, G = k_hh_g, mcl_cores = 32,
                                      Gx = gx, max.vector.size = 1e+07, mc.cores = 2, mc.cleanup = TRUE, mc.silent = TRUE)

names(g_cml_list) <- inflation_list
#
#save(g_cml_list, file = '/scratch/antonio/unk_C_SC/g_cml_list_c0.6_p50_all.Rda')
inflation_list <- seq(1.2, 3, 0.1)
load(file = '/scratch/antonio/unk_C_SC/g_cml_list_c0.6_p50_all.Rda')
names(g_cml_list) <- inflation_list

# Contract identified communities -----------------------------------------
k_gc <- pbmcapply::pbmclapply(g_cml_list, contract_graphs, G = k_hh_g, max.vector.size = 3e+09,  mc.cores = 19)
names(k_gc) <- inflation_list
# save(k_gc, file = "/scratch/antonio/unk_C_SC/k_gc.Rda")
load(file = "/scratch/antonio/unk_C_SC/k_gc.Rda", verbose = TRUE)
# Get some stats of the different communities -----------------------------
# Modularity
k_hh_g_modularity <- map_df(g_cml_list, function(X){
  vnames <- V(k_hh_g)$name
  tibble(modularity = modularity(k_hh_g, X$coms[match(vnames,X$coms$vertex),]$com))
}, .id = "inflation")

# Number of different DAs in each MCL cluster
k_hh_gc_da_sg <- lapply(k_gc, function(X) {
  ppbmmclapply(X$da$com, evaluate_da_components, da = X$da, mc.cores = 78) %>%
    bind_rows() %>% inner_join(X$da, by = "com")
})
names(k_hh_gc_da_sg) <- inflation_list
#save(k_hh_gc_da_sg, file = "/scratch/antonio/unk_C_SC/k_hh_gc_da_sg.Rda", compress = FALSE)
load(file = "/scratch/antonio/unk_C_SC/k_hh_gc_da_sg.Rda", verbose = TRUE)

# Entropy
k_hh_g_da <- k_hh_g %>%
  activate(nodes) %>%
  as_tibble() %>%
  separate_rows(archit, sep = "__") %>%
  as.data.table()

gc()
k_hh_gc_com_entropy <- lapply(k_hh_gc_da_sg, function(X) {
  pbmcapply::pbmclapply(X$com, FUN = evaluate_component_entropy, g = k_hh_g_da, df = X,  mc.cores = 32, max.vector.size = 3e+08) %>% bind_rows()
})

#save(k_hh_gc_com_entropy, file = "/scratch/antonio/unk_C_SC/k_hh_gc_com_entropy.Rda")
load(file = "/scratch/antonio/unk_C_SC/k_hh_gc_com_entropy.Rda")

k_hh_gc_com_entropy_summary <- map_df(k_hh_gc_com_entropy, function(X){
  X %>% mutate(e_0 = ifelse(entropy == 0, "eq", "gt")) %>%
    group_by(e_0) %>%
    count() %>%
    ungroup() %>%
    mutate(p = n/sum(n))
}, .id = "inflation")


# Inter and intra score modes
k_hh_g_dt <- k_hh_g %>% as_data_frame(what="edges") %>% rename(cl_name1 = from, cl_name2 = to, score_col = weight) %>% as.data.table()
k_orig_hh_mode <- pbmcapply::pbmclapply(g_cml_list, function(X){
  k_hh_g_dt %>%
    dt_left_join(X$coms %>% dt_select(vertex, com) %>% dt_mutate(com = as.character(com), vertex = as.character(vertex)) %>% rename(cl_name1 = vertex) %>% unique() %>% as.data.table, by = "cl_name1") %>%
    dt_left_join(X$coms %>% dt_select(vertex, com) %>% dt_mutate(com = as.character(com), vertex = as.character(vertex)) %>% rename(cl_name2 = vertex) %>% unique() %>% as.data.table, by = "cl_name2") %>%
    dt_filter(com.x == com.y) %>% mutate(com = com.x) %>% group_by(com) %>% summarise(mode = estimate_mode(score_col))
}, mc.cores = 19)

com_orig_intra_score <- map_df(k_orig_hh_mode, function(X){tibble(mode = estimate_mode(X$mode))}, .id = "inflation")
com_orig_inter_score <- map_df(k_gc, function(X){ tibble(mode = estimate_mode(E(X$graph)$weight))}, .id = "inflation")

k_orig_hh_mode_summary <- com_orig_intra_score %>%
  mutate(class = "intra") %>%
  bind_rows(com_orig_inter_score %>% mutate(class = "inter"))


# Get non-redundant set of DAs
# library(zoo)
# k_da_c <- map(p_doms %>% filter(class != "mono") %>% .$exp_rep %>% unique %>% strsplit("\\|"),
#               function(x) {
#                 g <- rollapply(data = x, 2, by=1, c) %>%
#                   graph_from_data_frame(directed = TRUE) %>% as_tbl_graph()
#               })


# Prefilter combinations that are too distant
# Get list of domains used during the MCL clustering
p_doms %>% filter(cl_name %in% unique_hhblits_k_cl)

d <- p_doms %>% select(class, exp_rep) %>% unique() %>% as.data.table()
lo_env$da_dist <- stringdist::stringdistmatrix(d$exp_rep, useNames = TRUE, method = "cosine", q = 3, nthread = 64) %>%
  broom::tidy() %>%
  as.data.table()
#fwrite(lo_env$da_dist, file = "/scratch/antonio/unk_C_SC/da_dist.tsv", col.names = TRUE)
lo_env$da_dist <- fread(input = "/scratch/antonio/unk_C_SC/da_dist.tsv", header = TRUE)

da_dist_0.9 <- lo_env$da_dist %>% dt_filter(distance < 0.9) %>% dt_arrange(-distance)

d[, n := seq_len(nrow(d))]

d.2 <- da_dist_0.9 %>% as_tibble() %>% mutate(idx1 = plyr::mapvalues(item1, from = d$exp_rep, to = d$n),
                                              idx2 = plyr::mapvalues(item2, from = d$exp_rep, to = d$n))

d.2 %>% filter(is.na(item1))
# Remove partial DAs that are contained in larger architectures


#iterations <- 1000
pb <- txtProgressBar(max = nproc, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

d.2 <- d.2 %>% as.data.table()

pb <- txtProgressBar(min = 0, max = nrow(d.2), style = 3)
k_das_refinement <- d.2[, c("n1", "n2", "V1_t", "V2_t", "DET") := {setTxtProgressBar(pb, .GRP);
  n1 = str_count(item1, pattern = "\\|");
  n2 = str_count(item2, pattern = "\\|");
  V1_t = ifelse(str_count(item1, pattern = "\\|") <= str_count(item2, pattern = "\\|"), as.character(item1), as.character(item2));
  V2_t = ifelse(str_count(item1, pattern = "\\|") > str_count(item2, pattern = "\\|"), as.character(item1), as.character(item2));
  DET = ifelse(grepl(V1_t, V2_t, fixed = TRUE), TRUE, FALSE); list(n1, n2, V1_t, V2_t, DET)}, seq_len(nrow(d.2))]
close(pb)

d.2

k_das_refinement <- foreach(i=isplitRows(d.2, chunks=4), .combine = bind_rows, .packages = c("tidyverse", "igraph", "tidygraph", "future", "furrr", "stringr"), .options.snow = opts) %do% {
  plan(multicore, workers = 64)
  future_pmap_dfr(i %>% as_tibble(),
                  ~tibble(idx1 = ..4,
                          idx2 = ..5,
                          item1 = ..1,
                          item2 = ..2,
                          n1 = str_count(item1, pattern = "\\|"),
                          n2 = str_count(item2, pattern = "\\|"),
                          V1_t = ifelse(str_count(item1, pattern = "\\|") <= str_count(item2, pattern = "\\|"), as.character(item1), as.character(item2)),
                          V2_t = ifelse(str_count(item1, pattern = "\\|") > str_count(item2, pattern = "\\|"), as.character(item1), as.character(item2)),
                          DET = ifelse(grepl(V1_t, V2_t, fixed = TRUE), TRUE, FALSE)), .progress = TRUE)
  #strn = ecount(graph.intersection(c[[..4]] %>% igraph::simplify(), c[[..5]] %>% igraph::simplify()))))

}

close(pb)
stopCluster(cl)

save(k_das_refinement, file = "/scratch/antonio/unk_C_SC/k_das_refinement.Rda", compress = FALSE)
#load(file = "/scratch/antonio/unk_C_SC/k_das_refinement.Rda", verbose = TRUE)
# Domains that are part of a larger domain
sub_da <- k_das_refinement  %>% dt_filter(DET == TRUE) %>% dt_select(V1_t) %>% unique()

# Large domains
large_da <- k_das_refinement %>% dt_filter(!(V2_t %in% sub_da$V1_t)) %>% dt_select(V2_t) %>% unique()

# get completeness
completness <- p_doms %>%
  group_by(archit) %>%
  summarise(complete = sum(complete),
            partial = sum(partial),
            N = complete + partial) %>%
  mutate(complete = complete/N,
         partial = partial/N)

# DAs that are part of a larger DA but at least they are seen in more than 75% of
# complete ORFs
sub_da <- k_das_refinement %>%
  dt_filter(DET == TRUE) %>%
  dt_select(V1_t) %>% dplyr::rename(archit = V1_t) %>% unique() %>%
  inner_join(completness) %>% arrange(desc(complete)) %>% filter(complete < 0.75) %>% as_tibble()

# create list with minimun set of DAs
# all DAs
final_da <- p_doms %>%
  select(exp_rep) %>%
  unique %>%
  filter(!(exp_rep %in% sub_da$archit))


# Summarize results -------------------------------------------------------
k_partition_stats <- map_df(g_cml_list, function(X) {
  tibble(com_orig_intra_score = (estimate_mode(X$intra_scores$mode)),
         com_orig_n = length(unique(X$coms$com)),
         com_orig_1mem = X$coms %>% group_by(com) %>% count() %>% filter(n == 1) %>% nrow())
}, .id = "inflation")


k_partition_stats <- k_partition_stats %>%
  inner_join(map_df(k_hh_gc_da_sg, function(X) {
    tibble(com_orig_1comp = X %>% filter(n_comp == 1) %>% nrow())
  }, .id = "inflation")) %>%
  mutate(com_orig_1com_prop = com_orig_1comp/com_orig_n,
         com_orig_1mem_prop = com_orig_1mem/com_orig_n)


# Plot results
# ggpubr::ggarrange(k_orig_hh_mode_summary %>%
#                     select(-mode) %>%
#                     mutate(inflation = fct_rev(inflation)) %>%
#                     ggplot(aes(x=inter, xend=intra, y=inflation, group=inflation)) +
#                     geom_dumbbell(size=1, color="#e3e2e1",
#                                   colour_x = "#5b8124", colour_xend = "#bad744",
#                                   dot_guide=TRUE, dot_guide_size=0.25) +
#                     theme_light(),
#                   k_partition_stats %>%
#                     mutate(inflation = fct_rev(inflation)) %>%
#                     ggplot(aes(inflation, com_orig_n)) +
#                     geom_lollipop() +
#                     ggpubr::rotate() +
#                     theme_light(),
#                   k_partition_stats %>%
#                     mutate(inflation = fct_rev(inflation)) %>%
#                     ggplot(aes(inflation, com_orig_1mem_prop)) +
#                     geom_lollipop() +
#                     ggpubr::rotate() +
#                     theme_light(),
#                   k_partition_stats %>%
#                     mutate(inflation = fct_rev(inflation)) %>%
#                     ggplot(aes(inflation, com_orig_1com_prop)) +
#                     geom_lollipop() +
#                     ggpubr::rotate() +
#                     theme_light(),
#                   k_hh_gc_com_entropy_summary %>%
#                     filter(e_0 == "eq") %>%
#                     mutate(inflation = fct_rev(inflation)) %>%
#                     ggplot(aes(inflation, p)) +
#                     geom_lollipop() +
#                     ggpubr::rotate() +
#                     theme_light(),
#                   nrow = 1, ncol = 5, align = "hv"
# )


# Make a decision for the best inflation value ----------------------------

narchs <- nrow(final_da)
k_partition_stats_eval <- k_partition_stats %>%
  inner_join(k_hh_gc_com_entropy_summary %>% filter(e_0 == "eq")) %>%
  mutate(com_orig_1mem_prop = 1 - com_orig_1mem_prop,
         com_orig_n = 1/abs(com_orig_n - narchs),
         com_orig_n = com_orig_n/max(com_orig_n),
         com_orig_inter_score = com_orig_intra_score/max(com_orig_intra_score)) %>%
  select(inflation, com_orig_1com_prop, com_orig_1mem_prop, p, com_orig_inter_score, com_orig_n) %>%
  rename(grp = inflation) %>%
  gather(var, value, -grp) %>%
  group_by(grp) %>%
  mutate(nextval = lead(value, default = value[1]),
         angle = (1/5)*(2*pi),
         area = value*nextval*sin(angle)/2) %>%
  mutate(total = sum(area)) %>%
  ungroup()


coord_radar <- function (theta = "x", start = 0, direction = 1)
{
  theta <- match.arg(theta, c("x", "y"))
  r <- if (theta == "x")
    "y"
  else "x"
  ggproto("CoordRadar", CoordPolar, theta = theta, r = r, start = start,
          direction = sign(direction),
          is_linear = function(coord) TRUE)
}

k_partition_stats_eval %>%
  mutate(var1 = plyr::mapvalues(var, from = c("com_orig_1com_prop", "com_orig_1mem_prop", "p", "com_orig_inter_score", "com_orig_n"),
                                to = 1:5), var1 = as.numeric(var1)) %>% ungroup() %>%
  ggplot(aes(var1, value)) +
  #geom_polygon() +
  geom_point() +
  geom_polygon(, fill = "grey", size = 0.2, color = "black", alpha = 0.5) +
  geom_text(aes(0,0, label = round(total, 2)), color = "red") +
  facet_wrap(~grp ) +
  scale_y_continuous("", limits = c(0, 1), expand = c(0,0)) +
  scale_x_continuous("", breaks = 1:5, expand = c(0,0)) +
  theme_bw() +
  coord_radar()

best_inflation <- k_partition_stats_eval %>%
  ungroup() %>%
  select(grp, total) %>%
  unique() %>%
  top_n(1) %>%
  .$grp

# Add missing clusters to MCL components ----------------------------------

# Add the missing nodes to existing communities
# Try to find the best hits for the missing nodes
# Based on coverage/probability/score per position



missing_dt <- lo_env$k_hhblits_missing %>% ungroup()
mcl_coms <- g_cml_list[[best_inflation]]
max_com <- max(mcl_coms$coms$com)

missing_ids %>% length()

# missing_ids -> 7047

# Try to identify these cluster that can have some homology to existing
# We only take the queries for missing ids
# We check that the target is in the MCL components
missing_d <- missing_dt %>%
  as.data.table() %>%
  dt_filter(cl_name1 %in% missing_ids) %>%
  dt_filter(cl_name2 %in% mcl_coms$coms$vertex) %>%
  dt_left_join(mcl_coms$coms %>% mutate(vertex = as.character(vertex)) %>% rename(cl_name2 = vertex))

# MCL clusters with more relaxed filters (p50 and cov >= 40)
# we just keep the best hit
missing_d_pass <- missing_d %>%
  dt_filter(probability >= 50, (q_cov >= 0.4 | t_cov >= 0.4)) %>%
  as_tibble() %>%
  group_by(cl_name1) %>%
  arrange(-probability, -q_cov, -t_cov) %>%
  #mutate(com1 = majority_vote(com)$majority,
  mutate(nhit = row_number()) %>%
  arrange(cl_name1, nhit) %>%
  ungroup() %>%
  #filter(probability >= 30, nhit <= 3) %>%
  filter(nhit == 1) %>%
  select(cl_name1, cl_name2, score_col, com) %>%
  rename(weight = score_col)

missing_d_pass$cl_name1 %>% unique() %>% length()

# Find the ones we couldn't classify
# 299
to_assign <- setdiff(missing_ids, missing_d_pass$cl_name1 %>% unique())

assigned <- setdiff(missing_ids, to_assign)

# We need to check if any of the remaining no assigned clusters have any homology
# to the just classified using more relaxed parameters

missing_d_pass_1 <- bind_rows(lo_env$k_hhblits_missing %>%
                                dt_inner_join(missing_d_pass %>% select(cl_name1, com)) %>%
                                dt_filter(cl_name2 %in% to_assign) %>% dt_mutate(cl_name = cl_name2, com_f = com),
                              lo_env$k_hhblits_missing %>%
                                dt_inner_join(missing_d_pass %>% select(cl_name2, com)) %>%
                                dt_filter(cl_name1 %in% to_assign) %>% dt_mutate(cl_name = cl_name1, com_f = com)) %>%
  dt_filter(probability > 50, (q_cov > 0.4 | t_cov > 0.4)) %>%
  group_by(cl_name) %>%
  arrange(-probability, -q_cov, -t_cov) %>%
  #mutate(com1 = majority_vote(com)$majority,
  mutate(nhit = row_number()) %>%
  arrange(cl_name, nhit) %>%
  ungroup() %>%
  #filter(probability >= 30, nhit <= 3) %>%
  filter(nhit == 1) %>%
  dt_select(cl_name, com_f) %>%
  dplyr::rename(com = com_f) %>%
  as_tibble()



# The ones that cannot be assigned, we will try to find any existing connections between no classified
# and run MCL with the best inflation value and create new clusters

missing_d_pass_2 <- bind_rows(lo_env$k_hhblits_missing %>%
                                dt_filter(cl_name1 %in% to_assign, cl_name2 %in% to_assign) %>%
                                dt_filter(probability > 50, (q_cov > 0.4 | t_cov > 0.4)) %>%
                                dt_inner_join(missing_d_pass_1 %>% select(cl_name, com) %>% rename(cl_name1 = cl_name)) %>%
                                dt_mutate(cl_name = cl_name1, com_f = com),
                              lo_env$k_hhblits_missing %>%
                                dt_filter(cl_name1 %in% to_assign, cl_name2 %in% to_assign) %>%
                                dt_filter(probability > 50, (q_cov > 0.4 | t_cov > 0.4)) %>%
                                dt_inner_join(missing_d_pass_1 %>% select(cl_name, com) %>% rename(cl_name2 = cl_name)) %>%
                                dt_mutate(cl_name = cl_name2, com_f = com)) %>%
  group_by(cl_name) %>%
  arrange(-probability, -q_cov, -t_cov) %>%
  #mutate(com1 = majority_vote(com)$majority,
  mutate(nhit = row_number()) %>%
  arrange(cl_name, nhit) %>%
  ungroup() %>%
  #filter(probability >= 30, nhit <= 3) %>%
  filter(nhit == 1) %>%
  dt_select(cl_name, com_f) %>%
  dplyr::rename(com = com_f) %>%
  as_tibble()

missing_ids_1 <- setdiff(missing_ids, c(missing_d_pass$cl_name1 %>% unique, missing_d_pass_1$cl_name, missing_d_pass_2$cl_name))

missing_d_pass_3 <- tibble(cl_name = missing_ids_1) %>% mutate(com = max_com + row_number())

k_components <-  bind_rows(missing_d_pass %>% select(cl_name1, com) %>% rename(cl_name = cl_name1),
                           missing_d_pass_1,
                           missing_d_pass_2,
                           missing_d_pass_3 %>% mutate(cl_name = as.character(cl_name)),
                           mcl_coms$coms %>% rename(cl_name = vertex) %>% mutate(cl_name = as.character(cl_name)) %>% select(cl_name, com)) %>% unique()

# Do we have all components
nrow(k_comp) == nrow(k_components)

readr::write_tsv(k_components, path = "/scratch/antonio/PRs/WF/EPA-NG/NEW/k_components.tsv", col_names = FALSE)

###########################################################################
# Get components for the KWP ----------------------------------------------
###########################################################################
lo_env$kwp_hhblits_all <- fread(input = "/scratch/antonio/unk_C_SC/data/kwp_hhblits.tsv.gz",
                                header = FALSE, verbose = TRUE, nThread = 80)
names(lo_env$kwp_hhblits_all) <- c("cl_name1", "cl_name2", "probability", "e-value", "Score",
                                   "Cols", "q_start", "q_stop", "t_start", "t_stop", "q_len", "t_len", "q_cov",
                                   "t_cov")

lo_env$kwp_hhblits <- lo_env$kwp_hhblits_all %>% dt_filter(probability > 50, q_cov >= 0.6, t_cov >= 0.6)
lo_env$kwp_hhblits <- lo_env$kwp_hhblits %>% dt_filter(cl_name1 != cl_name2)
# How many clusters do we have (removed self-hits)
# 522,765 p50;c0.6
kwp_unique_hhblits_cl <- c(lo_env$kwp_hhblits$cl_name1, lo_env$kwp_hhblits$cl_name2) %>% unique()
length(kwp_unique_hhblits_cl)

# Missing clusters
# Num: 27,244
kwp_missing_ids <- setdiff(cl_cat %>% filter(category == "KWP") %>% .$cl_name, kwp_unique_hhblits_cl)
length(kwp_missing_ids)

# Let's create a graph
lo_env$kwp_hhb_bh <- lo_env$kwp_hhblits %>%
  as_tibble() %>%
  mutate(cl_name1 = as.character(cl_name1),
         cl_name2 = as.character(cl_name2)) %>%
  mutate(score_col = Score/Cols)

kwp_hh_g <- lo_env$kwp_hhb_bh %>%
  select(cl_name1, cl_name2, score_col) %>%
  rename(weight = score_col) %>%
  igraph::graph_from_data_frame(directed = FALSE) %>%
  igraph::simplify(remove.multiple = TRUE, remove.loops = TRUE, edge.attr.comb = list("max")) %>%
  as_tbl_graph()

gx <- kwp_hh_g

w <- E(gx)$weight
w <- w - min(w) + 0.001

E(gx)$weight <- w


# We run mcl with different inflation parameters:
mcl_bin <- "mcl"


inflation_list <- seq(1.2, 3, 0.1)
mcl_bin <- "mcl"
kwp_g_cml_list <- pbmcapply::pbmclapply(inflation_list, optimal_mcl, G = kwp_hh_g,
                                        Gx = gx, max.vector.size = 1e+07, mc.cores = 2, mc.cleanup = TRUE, mc.silent = TRUE)

kwp_g_cml_list <- g_cml_list_kwp
names(kwp_g_cml_list) <- inflation_list
#save(kwp_g_cml_list, file = '/scratch/antonio/unk_C_SC/kwp_g_cml_list_c0.6_p50_all.Rda')

load(file = '/scratch/antonio/unk_C_SC/kwp_g_cml_list_c0.6_p50_all.Rda', verbose = TRUE)


kwp_partition_stats <- map_df(kwp_g_cml_list, function(X) {
  tibble(intra = (estimate_mode(X$intra_scores$mode)), ncomps = length(X$intra_scores$com))
}, .id = "inflation")

# Contract identified communities -----------------------------------------
kwp_gc <- pbmcapply::pbmclapply(kwp_g_cml_list, contract_graphs, G = kwp_hh_g, max.vector.size = 2e+09,  mc.cores = 19)
names(kwp_gc) <- inflation_list


# Get some stats of the different communities -----------------------------
# Modularity
kwp_hh_g_modularity <- map_df(kwp_g_cml_list, function(X){
  vnames <- V(kwp_hh_g)$name
  tibble(modularity = modularity(kwp_hh_g, X$coms[match(vnames,X$coms$vertex),]$com))
}, .id = "inflation")

# Inter and intra score modes

kwp_hh_g_dt <- kwp_hh_g %>% as_data_frame(what="edges") %>% rename(cl_name1 = from, cl_name2 = to, score_col = weight) %>% as.data.table()
kwp_orig_hh_mode <- pbmcapply::pbmclapply(kwp_g_cml_list, function(X){
  kwp_hh_g_dt %>%
    dt_left_join(X$coms %>% select(vertex, com) %>% mutate(com = as.character(com), vertex = as.character(vertex)) %>% rename(cl_name1 = vertex) %>% unique() %>% as.data.table, by = "cl_name1") %>%
    dt_left_join(X$coms %>% select(vertex, com) %>% mutate(com = as.character(com), vertex = as.character(vertex)) %>% rename(cl_name2 = vertex) %>% unique() %>% as.data.table, by = "cl_name2") %>%
    dt_filter(com.x == com.y) %>% mutate(com = com.x) %>% group_by(com) %>% summarise(mode = estimate_mode(score_col))
}, mc.cores = 19)

com_orig_intra_score_kwp <- map_df(kwp_orig_hh_mode, function(X){tibble(mode = estimate_mode(X$mode))}, .id = "inflation")
com_orig_inter_score_kwp <- map_df(kwp_gc, function(X){tibble(mode = estimate_mode(E(X$graph)$weight))}, .id = "inflation")

kwp_orig_hh_mode_summary <- com_orig_intra_score_kwp %>%
  mutate(class = "intra") %>%
  bind_rows(com_orig_inter_score_kwp %>% mutate(class = "inter"))

kwp_orig_hh_mode_summary %>%
  ggplot(aes(inflation, mode, group = class)) +
  geom_line() + facet_wrap(~class, scales = "free") +
  geom_point(shape = 21, fill = "grey", color = "black", alpha = 0.8)

kwp_partition_stats <- map_df(kwp_g_cml_list, function(X) {
  tibble(com_orig_intra_score = (estimate_mode(X$intra_scores$mode)),
         com_orig_n = length(unique(X$coms$com)),
         com_orig_1mem = X$coms %>% group_by(com) %>% count() %>% filter(n == 1) %>% nrow())
}, .id = "inflation")


# Add missing clusters to MCL components ----------------------------------

# Add the missing nodes to existing communities
# Try to find the best hits for the missing nodes
# Based on coverage/probability/score per position


lo_env$kwp_hhblits_missing <- lo_env$kwp_hhblits_all %>%
  dt_mutate(cl_name1 = as.character(cl_name1), cl_name2 = as.character(cl_name2)) %>%
  dt_filter(cl_name1 %in% kwp_missing_ids | cl_name2 %in% kwp_missing_ids) %>%
  as_tibble() %>% mutate(score_col = Score/Cols)


kwp_missing_dt <- lo_env$kwp_hhblits_missing %>% ungroup()
kwp_mcl_coms <- kwp_g_cml_list[[best_inflation]]
kwp_max_com <- max(kwp_mcl_coms$coms$com)

kwp_missing_ids %>% length()

# missing_ids -> 27244

# Try to identify these cluster that can have some homology to existing
# We only take the queries for missing ids
# We check that the target is in the MCL components
kwp_missing_d <- kwp_missing_dt %>%
  as.data.table() %>%
  dt_filter(cl_name1 %in% kwp_missing_ids) %>%
  dt_filter(cl_name2 %in% kwp_mcl_coms$coms$vertex) %>%
  dt_left_join(kwp_mcl_coms$coms %>% mutate(vertex = as.character(vertex)) %>% rename(cl_name2 = vertex))

# MCL clusters with more relaxed filters (p50 and cov >= 40)
# we just keep the best hit
kwp_missing_d_pass <- kwp_missing_d %>%
  dt_filter(probability >= 50, (q_cov >= 0.4 | t_cov >= 0.4)) %>%
  as_tibble() %>%
  group_by(cl_name1) %>%
  arrange(-probability, -q_cov, -t_cov) %>%
  #mutate(com1 = majority_vote(com)$majority,
  mutate(nhit = row_number()) %>%
  arrange(cl_name1, nhit) %>%
  ungroup() %>%
  #filter(probability >= 30, nhit <= 3) %>%
  filter(nhit == 1) %>%
  select(cl_name1, cl_name2, score_col, com) %>%
  rename(weight = score_col)

kwp_missing_d_pass$cl_name1 %>% unique() %>% length()

# Find the ones we couldn't classify
# 299
kwp_to_assign <- setdiff(kwp_missing_ids, kwp_missing_d_pass$cl_name1 %>% unique())

kwp_assigned <- setdiff(kwp_missing_ids, kwp_to_assign)

# We need to check if any of the remaining no assigned clusters have any homology
# to the just classified using more relaxed parameters

kwp_missing_d_pass_1 <- bind_rows(lo_env$kwp_hhblits_missing %>%
                                    dt_inner_join(kwp_missing_d_pass %>% select(cl_name1, com)) %>%
                                    dt_filter(cl_name2 %in% to_assign) %>% dt_mutate(cl_name = cl_name2, com_f = com),
                                  lo_env$kwp_hhblits_missing %>%
                                    dt_inner_join(kwp_missing_d_pass %>% select(cl_name2, com)) %>%
                                    dt_filter(cl_name1 %in% kwp_to_assign) %>% dt_mutate(cl_name = cl_name1, com_f = com)) %>%
  dt_filter(probability > 50, (q_cov > 0.4 | t_cov > 0.4)) %>%
  group_by(cl_name) %>%
  arrange(-probability, -q_cov, -t_cov) %>%
  #mutate(com1 = majority_vote(com)$majority,
  mutate(nhit = row_number()) %>%
  arrange(cl_name, nhit) %>%
  ungroup() %>%
  #filter(probability >= 30, nhit <= 3) %>%
  filter(nhit == 1) %>%
  dt_select(cl_name, com_f) %>%
  dplyr::rename(com = com_f) %>%
  as_tibble()



# The ones that cannot be assigned, we will try to find any existing connections between no classified
# and run MCL with the best inflation value and create new clusters

kwp_missing_d_pass_2 <- bind_rows(lo_env$kwp_hhblits_missing %>%
                                    dt_filter(cl_name1 %in% kwp_to_assign, cl_name2 %in% kwp_to_assign) %>%
                                    dt_filter(probability > 50, (q_cov > 0.4 | t_cov > 0.4)) %>%
                                    dt_inner_join(kwp_missing_d_pass_1 %>% select(cl_name, com) %>% rename(cl_name1 = cl_name)) %>%
                                    dt_mutate(cl_name = cl_name1, com_f = com),
                                  lo_env$kwp_hhblits_missing %>%
                                    dt_filter(cl_name1 %in% kwp_to_assign, cl_name2 %in% kwp_to_assign) %>%
                                    dt_filter(probability > 50, (q_cov > 0.4 | t_cov > 0.4)) %>%
                                    dt_inner_join(kwp_missing_d_pass_1 %>% select(cl_name, com) %>% rename(cl_name2 = cl_name)) %>%
                                    dt_mutate(cl_name = cl_name2, com_f = com)) %>%
  group_by(cl_name) %>%
  arrange(-probability, -q_cov, -t_cov) %>%
  #mutate(com1 = majority_vote(com)$majority,
  mutate(nhit = row_number()) %>%
  arrange(cl_name, nhit) %>%
  ungroup() %>%
  #filter(probability >= 30, nhit <= 3) %>%
  filter(nhit == 1) %>%
  dt_select(cl_name, com_f) %>%
  dplyr::rename(com = com_f) %>%
  as_tibble()

kwp_missing_ids_1 <- setdiff(kwp_missing_ids, c(kwp_missing_d_pass$cl_name1 %>% unique, kwp_missing_d_pass_1$cl_name, kwp_missing_d_pass_2$cl_name))

kwp_missing_d_pass_3 <- tibble(cl_name = kwp_missing_ids_1) %>% mutate(com = kwp_max_com + row_number())

kwp_components <-  bind_rows(kwp_missing_d_pass %>% select(cl_name1, com) %>% rename(cl_name = cl_name1),
                             kwp_missing_d_pass_1,
                             kwp_missing_d_pass_2,
                             kwp_missing_d_pass_3 %>% mutate(cl_name = as.character(cl_name)),
                             kwp_mcl_coms$coms %>% rename(cl_name = vertex) %>% mutate(cl_name = as.character(cl_name)) %>% select(cl_name, com)) %>% unique()

# Do we have all components
nrow(cl_cat %>% filter(category == "KWP")) == nrow(kwp_components)

readr::write_tsv(kwp_components, path = "/scratch/antonio/PRs/WF/EPA-NG/NEW/kwp_components.tsv", col_names = FALSE)


###########################################################################
# Get components for the GU  ----------------------------------------------
###########################################################################
lo_env$gu_hhblits_all <- fread(input = "zcat /bioinf/projects/megx/UNKNOWNS/chiara/unkn_hhblits/gu/gu_hhblits.tsv.gz",
                               header = FALSE, verbose = TRUE, nThread = 80)
names(lo_env$gu_hhblits_all) <- c("cl_name1", "cl_name2", "probability", "e-value", "Score",
                                  "Cols", "q_start", "q_stop", "t_start", "t_stop", "q_len", "t_len", "q_cov",
                                  "t_cov")

lo_env$gu_hhblits <- lo_env$gu_hhblits_all %>% dt_filter(probability > 50, q_cov >= 0.6, t_cov >= 0.6)
lo_env$gu_hhblits <- lo_env$gu_hhblits %>% dt_filter(cl_name1 != cl_name2)
# How many clusters do we have (removed self-hits)
# 765,848 p50;c0.6
gu_unique_hhblits_cl <- c(lo_env$gu_hhblits$cl_name1, lo_env$gu_hhblits$cl_name2) %>% unique()
length(gu_unique_hhblits_cl)

# Missing clusters
# Num: 26,203
gu_missing_ids <- setdiff(cl_cat %>% filter(category == "GU") %>% .$cl_name, gu_unique_hhblits_cl)


# Let's create a graph
lo_env$gu_hhb_bh <- lo_env$gu_hhblits %>%
  as_tibble() %>%
  mutate(cl_name1 = as.character(cl_name1),
         cl_name2 = as.character(cl_name2)) %>%
  mutate(score_col = Score/Cols)

gu_hh_g <- lo_env$gu_hhb_bh %>%
  select(cl_name1, cl_name2, score_col) %>%
  rename(weight = score_col) %>%
  igraph::graph_from_data_frame(directed = FALSE) %>%
  igraph::simplify(remove.multiple = TRUE, remove.loops = TRUE, edge.attr.comb = list("max")) %>%
  as_tbl_graph()

gx <- gu_hh_g

w <- E(gx)$weight
w <- w - min(w) + 0.001

E(gx)$weight <- w


# We run mcl with different inflation parameters:
mcl_bin <- "mcl"


inflation_list <- seq(1.2, 3, 0.1)
mcl_bin <- "mcl"
gu_g_cml_list <- pbmcapply::pbmclapply(inflation_list, optimal_mcl, G = gu_hh_g,
                                       Gx = gx, max.vector.size = 1e+07, mc.cores = 2, mc.cleanup = TRUE, mc.silent = TRUE)

names(gu_g_cml_list) <- inflation_list
save(gu_g_cml_list, file = '/scratch/antonio/unk_C_SC/gu_g_cml_list_c0.6_p50_all.Rda')

#load(file = '/scratch/antonio/unk_C_SC/g_cml_list_gu_c0.6_p50_all.Rda', verbose = TRUE)


gu_partition_stats <- map_df(gu_g_cml_list, function(X) {
  tibble(intra = (estimate_mode(X$intra_scores$mode)), ncomps = length(X$intra_scores$com))
}, .id = "inflation")

# Contract identified communities -----------------------------------------
gu_gc <- pbmcapply::pbmclapply(gu_g_cml_list, contract_graphs, G = gu_hh_g, max.vector.size = 2e+09,  mc.cores = 19)
names(gu_gc) <- inflation_list


# Get some stats of the different communities -----------------------------
# Modularity
gu_hh_g_modularity <- map_df(gu_g_cml_list, function(X){
  vnames <- V(gu_hh_g)$name
  tibble(modularity = modularity(gu_hh_g, X$coms[match(vnames,X$coms$vertex),]$com))
}, .id = "inflation")

# Inter and intra score modes

gu_hh_g_dt <- gu_hh_g %>% as_data_frame(what="edges") %>% rename(cl_name1 = from, cl_name2 = to, score_col = weight) %>% as.data.table()
gu_orig_hh_mode <- pbmcapply::pbmclapply(gu_g_cml_list, function(X){
  gu_hh_g_dt %>%
    dt_left_join(X$coms %>% select(vertex, com) %>% mutate(com = as.character(com), vertex = as.character(vertex)) %>% rename(cl_name1 = vertex) %>% unique() %>% as.data.table, by = "cl_name1") %>%
    dt_left_join(X$coms %>% select(vertex, com) %>% mutate(com = as.character(com), vertex = as.character(vertex)) %>% rename(cl_name2 = vertex) %>% unique() %>% as.data.table, by = "cl_name2") %>%
    dt_filter(com.x == com.y) %>% mutate(com = com.x) %>% group_by(com) %>% summarise(mode = estimate_mode(score_col))
}, mc.cores = 19)

com_orig_intra_score_gu <- map_df(gu_orig_hh_mode, function(X){tibble(mode = estimate_mode(X$mode))}, .id = "inflation")
com_orig_inter_score_gu <- map_df(gu_gc, function(X){tibble(mode = estimate_mode(E(X$graph)$weight))}, .id = "inflation")

gu_orig_hh_mode_summary <- com_orig_intra_score_gu %>%
  mutate(class = "intra") %>%
  bind_rows(com_orig_inter_score_gu %>% mutate(class = "inter"))

gu_orig_hh_mode_summary %>%
  ggplot(aes(inflation, mode, group = class)) +
  geom_line() + facet_wrap(~class, scales = "free") +
  geom_point(shape = 21, fill = "grey", color = "black", alpha = 0.8)

gu_partition_stats <- map_df(gu_g_cml_list, function(X) {
  tibble(com_orig_intra_score = (estimate_mode(X$intra_scores$mode)),
         com_orig_n = length(unique(X$coms$com)),
         com_orig_1mem = X$coms %>% group_by(com) %>% count() %>% filter(n == 1) %>% nrow())
}, .id = "inflation")


# Add missing clusters to MCL components ----------------------------------

# Add the missing nodes to existing communities
# Try to find the best hits for the missing nodes
# Based on coverage/probability/score per position


lo_env$gu_hhblits_missing <- lo_env$gu_hhblits_all %>%
  dt_mutate(cl_name1 = as.character(cl_name1), cl_name2 = as.character(cl_name2)) %>%
  dt_filter(cl_name1 %in% gu_missing_ids| cl_name2 %in% gu_missing_ids) %>%
  as_tibble() %>% mutate(score_col = Score/Cols)


gu_missing_dt <- lo_env$gu_hhblits_missing %>% ungroup()
gu_mcl_coms <- gu_g_cml_list[[best_inflation]]
gu_max_com <- max(gu_mcl_coms$coms$com)

gu_missing_ids %>% length()

# missing_ids -> 7047

# Try to identify these cluster that can have some homology to existing
# We only take the queries for missing ids
# We check that the target is in the MCL components
gu_missing_d <- gu_missing_dt %>%
  as.data.table() %>%
  dt_filter(cl_name1 %in% gu_missing_ids) %>%
  dt_filter(cl_name2 %in% gu_mcl_coms$coms$vertex) %>%
  dt_left_join(gu_mcl_coms$coms %>% mutate(vertex = as.character(vertex)) %>% rename(cl_name2 = vertex))

# MCL clusters with more relaxed filters (p50 and cov >= 40)
# we just keep the best hit
gu_missing_d_pass <- gu_missing_d %>%
  dt_filter(probability >= 50, (q_cov >= 0.4 | t_cov >= 0.4)) %>%
  as_tibble() %>%
  group_by(cl_name1) %>%
  arrange(-probability, -q_cov, -t_cov) %>%
  #mutate(com1 = majority_vote(com)$majority,
  mutate(nhit = row_number()) %>%
  arrange(cl_name1, nhit) %>%
  ungroup() %>%
  #filter(probability >= 30, nhit <= 3) %>%
  filter(nhit == 1) %>%
  select(cl_name1, cl_name2, score_col, com) %>%
  rename(weight = score_col)

gu_missing_d_pass$cl_name1 %>% unique() %>% length()

# Find the ones we couldn't classify
# 299
gu_to_assign <- setdiff(gu_missing_ids, gu_missing_d_pass$cl_name1 %>% unique())

gu_assigned <- setdiff(gu_missing_ids, gu_to_assign)

# We need to check if any of the remaining no assigned clusters have any homology
# to the just classified using more relaxed parameters

gu_missing_d_pass_1 <- bind_rows(lo_env$gu_hhblits_missing %>%
                                   dt_inner_join(gu_missing_d_pass %>% select(cl_name1, com)) %>%
                                   dt_filter(cl_name2 %in% to_assign) %>% dt_mutate(cl_name = cl_name2, com_f = com),
                                 lo_env$gu_hhblits_missing %>%
                                   dt_inner_join(gu_missing_d_pass %>% select(cl_name2, com)) %>%
                                   dt_filter(cl_name1 %in% gu_to_assign) %>% dt_mutate(cl_name = cl_name1, com_f = com)) %>%
  dt_filter(probability > 50, (q_cov > 0.4 | t_cov > 0.4)) %>%
  group_by(cl_name) %>%
  arrange(-probability, -q_cov, -t_cov) %>%
  #mutate(com1 = majority_vote(com)$majority,
  mutate(nhit = row_number()) %>%
  arrange(cl_name, nhit) %>%
  ungroup() %>%
  #filter(probability >= 30, nhit <= 3) %>%
  filter(nhit == 1) %>%
  dt_select(cl_name, com_f) %>%
  dplyr::rename(com = com_f) %>%
  as_tibble()



# The ones that cannot be assigned, we will try to find any existing connections between no classified
# and run MCL with the best inflation value and create new clusters

gu_missing_d_pass_2 <- bind_rows(lo_env$gu_hhblits_missing %>%
                                   dt_filter(cl_name1 %in% gu_to_assign, cl_name2 %in% gu_to_assign) %>%
                                   dt_filter(probability > 50, (q_cov > 0.4 | t_cov > 0.4)) %>%
                                   dt_inner_join(gu_missing_d_pass_1 %>% select(cl_name, com) %>% rename(cl_name1 = cl_name)) %>%
                                   dt_mutate(cl_name = cl_name1, com_f = com),
                                 lo_env$gu_hhblits_missing %>%
                                   dt_filter(cl_name1 %in% gu_to_assign, cl_name2 %in% gu_to_assign) %>%
                                   dt_filter(probability > 50, (q_cov > 0.4 | t_cov > 0.4)) %>%
                                   dt_inner_join(gu_missing_d_pass_1 %>% select(cl_name, com) %>% rename(cl_name2 = cl_name)) %>%
                                   dt_mutate(cl_name = cl_name2, com_f = com)) %>%
  group_by(cl_name) %>%
  arrange(-probability, -q_cov, -t_cov) %>%
  #mutate(com1 = majority_vote(com)$majority,
  mutate(nhit = row_number()) %>%
  arrange(cl_name, nhit) %>%
  ungroup() %>%
  #filter(probability >= 30, nhit <= 3) %>%
  filter(nhit == 1) %>%
  dt_select(cl_name, com_f) %>%
  dplyr::rename(com = com_f) %>%
  as_tibble()

gu_missing_ids_1 <- setdiff(gu_missing_ids, c(gu_missing_d_pass$cl_name1 %>% unique, gu_missing_d_pass_1$cl_name, gu_missing_d_pass_2$cl_name))

gu_missing_d_pass_3 <- tibble(cl_name = gu_missing_ids_1) %>% mutate(com = gu_max_com + row_number())

gu_components <-  bind_rows(gu_missing_d_pass %>% select(cl_name1, com) %>% rename(cl_name = cl_name1),
                            gu_missing_d_pass_1,
                            gu_missing_d_pass_2,
                            gu_missing_d_pass_3 %>% mutate(cl_name = as.character(cl_name)),
                            gu_mcl_coms$coms %>% rename(cl_name = vertex) %>% mutate(cl_name = as.character(cl_name)) %>% select(cl_name, com)) %>% unique()

# Do we have all components
nrow(cl_cat %>% filter(category == "GU")) == nrow(gu_components)

readr::write_tsv(gu_components, path = "/scratch/antonio/PRs/WF/EPA-NG/NEW/gu_components.tsv", col_names = FALSE)


###########################################################################
# Get components for the EU  ----------------------------------------------
###########################################################################
lo_env$eu_hhblits_all <- fread(input = "zcat /bioinf/projects/megx/UNKNOWNS/chiara/unkn_hhblits/eu/eu_hhblits.tsv.gz",
                               header = FALSE, verbose = TRUE, nThread = 80)
names(lo_env$eu_hhblits_all) <- c("cl_name1", "cl_name2", "probability", "e-value", "Score",
                                  "Cols", "q_start", "q_stop", "t_start", "t_stop", "q_len", "t_len", "q_cov",
                                  "t_cov")

lo_env$eu_hhblits <- lo_env$eu_hhblits_all %>% dt_filter(probability > 50, q_cov >= 0.6, t_cov >= 0.6)
lo_env$eu_hhblits <- lo_env$eu_hhblits %>% dt_filter(cl_name1 != cl_name2)
# How many clusters do we have (removed self-hits)
# 765,848 p50;c0.6
eu_unique_hhblits_cl <- c(lo_env$eu_hhblits$cl_name1, lo_env$eu_hhblits$cl_name2) %>% unique()
length(eu_unique_hhblits_cl)

# Missing clusters
# Num: 26,203
eu_missing_ids <- setdiff(cl_cat %>% filter(category == "EU") %>% .$cl_name, eu_unique_hhblits_cl)


# Let's create a graph
lo_env$eu_hhb_bh <- lo_env$eu_hhblits %>%
  as_tibble() %>%
  mutate(cl_name1 = as.character(cl_name1),
         cl_name2 = as.character(cl_name2)) %>%
  mutate(score_col = Score/Cols)

eu_hh_g <- lo_env$eu_hhb_bh %>%
  select(cl_name1, cl_name2, score_col) %>%
  rename(weight = score_col) %>%
  igraph::graph_from_data_frame(directed = FALSE) %>%
  igraph::simplify(remove.multiple = TRUE, remove.loops = TRUE, edge.attr.comb = list("max")) %>%
  as_tbl_graph()

gx <- eu_hh_g

w <- E(gx)$weight
w <- w - min(w) + 0.001

E(gx)$weight <- w


# We run mcl with different inflation parameters:
mcl_bin <- "mcl"


inflation_list <- seq(1.2, 3, 0.1)
mcl_bin <- "mcl"
eu_g_cml_list <- pbmcapply::pbmclapply(inflation_list, optimal_mcl, G = eu_hh_g,
                                       Gx = gx, max.vector.size = 1e+07, mc.cores = 2, mc.cleanup = TRUE, mc.silent = TRUE)

names(eu_g_cml_list) <- inflation_list
# save(eu_g_cml_list, file = '/scratch/antonio/unk_C_SC/eu_g_cml_list_c0.6_p50_all.Rda')

#load(file = '/scratch/antonio/unk_C_SC/g_cml_list_eu_c0.6_p50_all.Rda', verbose = TRUE)


eu_partition_stats <- map_df(eu_g_cml_list, function(X) {
  tibble(intra = (estimate_mode(X$intra_scores$mode)), ncomps = length(X$intra_scores$com))
}, .id = "inflation")

# Contract identified communities -----------------------------------------
eu_gc <- pbmcapply::pbmclapply(eu_g_cml_list, contract_graphs, G = eu_hh_g, max.vector.size = 2e+09,  mc.cores = 19)
names(eu_gc) <- inflation_list


# Get some stats of the different communities -----------------------------
# Modularity
eu_hh_g_modularity <- map_df(eu_g_cml_list, function(X){
  vnames <- V(eu_hh_g)$name
  tibble(modularity = modularity(eu_hh_g, X$coms[match(vnames,X$coms$vertex),]$com))
}, .id = "inflation")

# Inter and intra score modes

eu_hh_g_dt <- eu_hh_g %>% as_data_frame(what="edges") %>% rename(cl_name1 = from, cl_name2 = to, score_col = weight) %>% as.data.table()
eu_orig_hh_mode <- pbmcapply::pbmclapply(eu_g_cml_list, function(X){
  eu_hh_g_dt %>%
    dt_left_join(X$coms %>% select(vertex, com) %>% mutate(com = as.character(com), vertex = as.character(vertex)) %>% rename(cl_name1 = vertex) %>% unique() %>% as.data.table, by = "cl_name1") %>%
    dt_left_join(X$coms %>% select(vertex, com) %>% mutate(com = as.character(com), vertex = as.character(vertex)) %>% rename(cl_name2 = vertex) %>% unique() %>% as.data.table, by = "cl_name2") %>%
    dt_filter(com.x == com.y) %>% mutate(com = com.x) %>% group_by(com) %>% summarise(mode = estimate_mode(score_col))
}, mc.cores = 19)

com_orig_intra_score_eu <- map_df(eu_orig_hh_mode, function(X){tibble(mode = estimate_mode(X$mode))}, .id = "inflation")
com_orig_inter_score_eu <- map_df(eu_gc, function(X){tibble(mode = estimate_mode(E(X$graph)$weight))}, .id = "inflation")

eu_orig_hh_mode_summary <- com_orig_intra_score_eu %>%
  mutate(class = "intra") %>%
  bind_rows(com_orig_inter_score_eu %>% mutate(class = "inter"))

eu_orig_hh_mode_summary %>%
  ggplot(aes(inflation, mode, group = class)) +
  geom_line() + facet_wrap(~class, scales = "free") +
  geom_point(shape = 21, fill = "grey", color = "black", alpha = 0.8)

eu_partition_stats <- map_df(eu_g_cml_list, function(X) {
  tibble(com_orig_intra_score = (estimate_mode(X$intra_scores$mode)),
         com_orig_n = length(unique(X$coms$com)),
         com_orig_1mem = X$coms %>% group_by(com) %>% count() %>% filter(n == 1) %>% nrow())
}, .id = "inflation")


# Add missing clusters to MCL components ----------------------------------

# Add the missing nodes to existing communities
# Try to find the best hits for the missing nodes
# Based on coverage/probability/score per position


lo_env$eu_hhblits_missing <- lo_env$eu_hhblits_all %>%
  dt_mutate(cl_name1 = as.character(cl_name1), cl_name2 = as.character(cl_name2)) %>%
  dt_filter(cl_name1 %in% eu_missing_ids| cl_name2 %in% eu_missing_ids) %>%
  as_tibble() %>% mutate(score_col = Score/Cols)


eu_missing_dt <- lo_env$eu_hhblits_missing %>% ungroup()
eu_mcl_coms <- eu_g_cml_list[[best_inflation]]
eu_max_com <- max(eu_mcl_coms$coms$com)

eu_missing_ids %>% length()

# missing_ids -> 7047

# Try to identify these cluster that can have some homology to existing
# We only take the queries for missing ids
# We check that the target is in the MCL components
eu_missing_d <- eu_missing_dt %>%
  as.data.table() %>%
  dt_filter(cl_name1 %in% eu_missing_ids) %>%
  dt_filter(cl_name2 %in% eu_mcl_coms$coms$vertex) %>%
  dt_left_join(eu_mcl_coms$coms %>% mutate(vertex = as.character(vertex)) %>% rename(cl_name2 = vertex))

# MCL clusters with more relaxed filters (p50 and cov >= 40)
# we just keep the best hit
eu_missing_d_pass <- eu_missing_d %>%
  dt_filter(probability >= 50, (q_cov >= 0.4 | t_cov >= 0.4)) %>%
  as_tibble() %>%
  group_by(cl_name1) %>%
  arrange(-probability, -q_cov, -t_cov) %>%
  #mutate(com1 = majority_vote(com)$majority,
  mutate(nhit = row_number()) %>%
  arrange(cl_name1, nhit) %>%
  ungroup() %>%
  #filter(probability >= 30, nhit <= 3) %>%
  filter(nhit == 1) %>%
  select(cl_name1, cl_name2, score_col, com) %>%
  rename(weight = score_col)

eu_missing_d_pass$cl_name1 %>% unique() %>% length()

# Find the ones we couldn't classify
# 299
eu_to_assign <- setdiff(eu_missing_ids, eu_missing_d_pass$cl_name1 %>% unique())

eu_assigned <- setdiff(eu_missing_ids, eu_to_assign)

# We need to check if any of the remaining no assigned clusters have any homology
# to the just classified using more relaxed parameters

eu_missing_d_pass_1 <- bind_rows(lo_env$eu_hhblits_missing %>%
                                   dt_inner_join(eu_missing_d_pass %>% select(cl_name1, com)) %>%
                                   dt_filter(cl_name2 %in% to_assign) %>% dt_mutate(cl_name = cl_name2, com_f = com),
                                 lo_env$eu_hhblits_missing %>%
                                   dt_inner_join(eu_missing_d_pass %>% select(cl_name2, com)) %>%
                                   dt_filter(cl_name1 %in% eu_to_assign) %>% dt_mutate(cl_name = cl_name1, com_f = com)) %>%
  dt_filter(probability > 50, (q_cov > 0.4 | t_cov > 0.4)) %>%
  group_by(cl_name) %>%
  arrange(-probability, -q_cov, -t_cov) %>%
  #mutate(com1 = majority_vote(com)$majority,
  mutate(nhit = row_number()) %>%
  arrange(cl_name, nhit) %>%
  ungroup() %>%
  #filter(probability >= 30, nhit <= 3) %>%
  filter(nhit == 1) %>%
  dt_select(cl_name, com_f) %>%
  dplyr::rename(com = com_f) %>%
  as_tibble()



# The ones that cannot be assigned, we will try to find any existing connections between no classified
# and run MCL with the best inflation value and create new clusters

eu_missing_d_pass_2 <- bind_rows(lo_env$eu_hhblits_missing %>%
                                   dt_filter(cl_name1 %in% eu_to_assign, cl_name2 %in% eu_to_assign) %>%
                                   dt_filter(probability > 50, (q_cov > 0.4 | t_cov > 0.4)) %>%
                                   dt_inner_join(eu_missing_d_pass_1 %>% select(cl_name, com) %>% rename(cl_name1 = cl_name)) %>%
                                   dt_mutate(cl_name = cl_name1, com_f = com),
                                 lo_env$eu_hhblits_missing %>%
                                   dt_filter(cl_name1 %in% eu_to_assign, cl_name2 %in% eu_to_assign) %>%
                                   dt_filter(probability > 50, (q_cov > 0.4 | t_cov > 0.4)) %>%
                                   dt_inner_join(eu_missing_d_pass_1 %>% select(cl_name, com) %>% rename(cl_name2 = cl_name)) %>%
                                   dt_mutate(cl_name = cl_name2, com_f = com)) %>%
  group_by(cl_name) %>%
  arrange(-probability, -q_cov, -t_cov) %>%
  #mutate(com1 = majority_vote(com)$majority,
  mutate(nhit = row_number()) %>%
  arrange(cl_name, nhit) %>%
  ungroup() %>%
  #filter(probability >= 30, nhit <= 3) %>%
  filter(nhit == 1) %>%
  dt_select(cl_name, com_f) %>%
  dplyr::rename(com = com_f) %>%
  as_tibble()

eu_missing_ids_1 <- setdiff(eu_missing_ids, c(eu_missing_d_pass$cl_name1 %>% unique, eu_missing_d_pass_1$cl_name, eu_missing_d_pass_2$cl_name))

eu_missing_d_pass_3 <- tibble(cl_name = eu_missing_ids_1) %>% mutate(com = eu_max_com + row_number())

eu_components <-  bind_rows(eu_missing_d_pass %>% select(cl_name1, com) %>% rename(cl_name = cl_name1),
                            eu_missing_d_pass_1,
                            eu_missing_d_pass_2,
                            eu_missing_d_pass_3 %>% mutate(cl_name = as.character(cl_name)),
                            eu_mcl_coms$coms %>% rename(cl_name = vertex) %>% mutate(cl_name = as.character(cl_name)) %>% select(cl_name, com)) %>% unique()

# Do we have all components
nrow(cl_cat %>% filter(category == "EU")) == nrow(eu_components)

readr::write_tsv(eu_components, path = "/scratch/antonio/PRs/WF/EPA-NG/NEW/eu_components.tsv", col_names = FALSE)

bind_rows(k_components %>% mutate(com = paste0("k_c_", com), category = "k"),
          kwp_components %>% mutate(com = paste0("kwp_c_", com), category = "kwp"),
          gu_components %>% mutate(com = paste0("gu_c_", com), category = "gu"),
          eu_components %>% mutate(com = paste0("eu_c_", com), category = "eu")
) %>%
  write_tsv(path = "/scratch/antonio/GTDB/components_20180209.tsv")

