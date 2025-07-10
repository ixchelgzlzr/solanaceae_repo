##############################################################
# Contrast ages for time homogeneous and time heterogeneous
##############################################################

library(tidyr)

# Read time homo and time hetero EXTANT traces
  
  # specify files paths
path_THom = "output/final_runs/time_homogeneous/extant_trees_combined.trees"
path_THet = "output/final_runs/time_heterogeneous/extant_trees_combined.trees"

   # read trees
tHom <- read.tree(path_THom)
tHet <- read.tree(path_THet)

#---------------
# define clades
#---------------

Solanaceae <- tHom[[1]]$tip.label


Solanum <- c( "Solanum_annuum", "Solanum_barbisetum", "Solanum_caripense", "Solanum_corymbosum", "Solanum_crinitum", "Solanum_crispum", "Solanum_evolvulifolium", "Solanum_evolvuloides", "Solanum_herculeum", "Solanum_hieronymi", "Solanum_jamaicense", "Solanum_laciniatum", "Solanum_mahoriense", "Solanum_microdontum", "Solanum_montanum", "Solanum_oedipus", "Solanum_palustre", "Solanum_peruvianum", "Solanum_polygamum", "Solanum_salamancae", "Solanum_sejunctum", "Solanum_sinuatirecurvum", "Solanum_sisymbriifolium", "Solanum_terminale", "Solanum_thelopodium", "Solanum_torvum", "Solanum_trizygum", "Solanum_uncinellum", "Solanum_valdiviense", "Solanum_weddellii", "Solanum_pimpinellifolium", "Solanum_pinnatum")


Datureae <- c("Trompettia_cardenasiana", "Datura_ferox", "Brugmansia_suaveolens")

Hyoscyaminae<- c("Anisodus_tanguticus", "Atropanthe_sinensis", "Hyoscyamus_niger", "Physochlaina_orientalis", "Przewalskia_tangutica", "Scoolia_carniolica")

Cestroideae<- c("Sessea_stipulata", "Cestrum_bracteatum", "Vestia_foetida", "Reyesia_chilensis", "Salpiglossis_sinuata", "Streptosolen_jamesonii", "Browallia_americana", "Protoschwenkia_mandonii")

Capsiceae<- c("Lycianthes_acapulcensis", "Lycianthes_rantonnetii", "Capsicum_chacoense")

Physalideae<- c("Physalis_angulata", "Physalis_peruviana", "Quincula_lobata", "Chamaesaracha_arida", "Calliphysalis_carpenteri", "Alkekengi_officinarum", "Capsicophysalis_potosina", "Leucophysalis_grandiflora", "Tzeltalia_amphitricha", "Schraderanthus_viscosus", "Witheringia_solanaceae", "Brachistus_nelsonii", "Oryctes_nevadensis", "Physaliastrum_heterophyllum", "Tubocapsicum_anomalum", "Archyphysalis_chamaesarachoides", "Nothocestrum_latifolium", "Discopodium_penninervium", "Withania_somnifera", "Mellissia_begoniifolia", "Vassobia_breviflora", "Eriolarynx_australis", "Iochroma_cyaneum", "Saracha_quitensis", "Dunalia_brachyacantha", "Trozellia_umbellata", "Deprea_pauciflora", "Deprea_abra_patriciae", "Darcyanthus_spruceanus", "Athenaea_cuspidata", "Athenaea_brasiliana", "Cuatresia_harlingiana")

Solanoideae<-  c("Quincula_lobata", "Chamaesaracha_arida", "Physalis_angulata", "Physalis_peruviana", "Alkekengi_officinarum", "Calliphysalis_carpenteri", "Capsicophysalis_potosina", "Leucophysalis_grandiflora", "Schraderanthus_viscosus", "Tzeltalia_amphitricha", "Brachistus_nelsonii", "Witheringia_solanacea", "Oryctes_nevadensis", "Physaliastrum_heterophyllum", "Tubocapsicum_anomalum", "Nothocestrum_latifolium", "Archiphysalis_chamaesarachoides", "Discopodium_penninervium", "Withania_somnifera", "Mellissia_begoniifolia", "Eriolarynx_australis", "Iochroma_cyaneum", "Vassobia_breviflora", "Saracha_quitensis", "Dunalia_brachyacantha", "Trozelia_umbellata", "Deprea_abra_patriciae", "Darcyanthus_spruceanus", "Deprea_pauciflora", "Athenaea_cuspidata", "Athenaea_brasiliana", "Cuatresia_harlingiana", "Lycianthes_acapulcensis", "Lycianthes_rantonnetii", "Capsicum_chacoense", "Nectouxia_formosa", "Salpichroa_tristis", "Solanum_annuum", "Solanum_salamancae", "Solanum_weddellii", "Solanum_sinuatirecurvum", "Solanum_corymbosum", "Solanum_crispum", "Solanum_uncinellum", "Solanum_valdiviense", "Solanum_terminale", "Solanum_herculeum", "Solanum_laciniatum", "Solanum_peruvianum", "Solanum_pimpinellifolium", "Solanum_palustre", "Solanum_microdontum", "Solanum_caripense", "Solanum_evolvulifolium", "Solanum_trizygum", "Solanum_pinnatum", "Solanum_montanum", "Solanum_mahoriense", "Solanum_sejunctum", "Solanum_barbisetum", "Solanum_oedipus", "Solanum_hieronymi", "Solanum_jamaicense", "Solanum_torvum", "Solanum_sisymbriifolium", "Solanum_crinitum", "Solanum_polygamum", "Solanum_evolvuloides", "Solanum_thelopodium", "Jaltomata_bicolor", "Jaltomata_grandiflora", "Datura_ferox", "Brugmansia_suaveolens", "Trompettia_cardenasiana", "Hawkesiophyton_panamense", "Juanulloa_mexicana", "Merinthopodium_neuranthum", "Markea_coccinea", "Dyssochroma_viridiflorum", "Solandra_grandiflora", "Schultesianthus_leucanthus", "Exodeconus_miersii", "Nicandra_physalodes", "Mandragora_officinarum", "Przewalskia_tangutica", "Physochlaina_orientalis", "Scopolia_carniolica", "Hyoscyamus_niger", "Atropanthe_sinensis", "Anisodus_tanguticus", "Atropa_belladonna", "Jaborosa_laciniata", "Sclerophylax_spinescens", "Nolana_divaricata", "Lycium_americanum", "Latua_pubiflora")

Dodecachroma<- c("Quincula_lobata", "Chamaesaracha_arida", "Physalis_angulata", "Physalis_peruviana", "Alkekengi_officinarum", "Calliphysalis_carpenteri", "Capsicophysalis_potosina", "Leucophysalis_grandiflora", "Schraderanthus_viscosus", "Tzeltalia_amphitricha", "Brachistus_nelsonii", "Witheringia_solanacea", "Oryctes_nevadensis", "Physaliastrum_heterophyllum", "Tubocapsicum_anomalum", "Nothocestrum_latifolium", "Archiphysalis_chamaesarachoides", "Discopodium_penninervium", "Withania_somnifera", "Mellissia_begoniifolia", "Eriolarynx_australis", "Iochroma_cyaneum", "Vassobia_breviflora", "Saracha_quitensis", "Dunalia_brachyacantha", "Trozelia_umbellata", "Deprea_abra_patriciae", "Darcyanthus_spruceanus", "Deprea_pauciflora", "Athenaea_cuspidata", "Athenaea_brasiliana", "Cuatresia_harlingiana", "Lycianthes_acapulcensis", "Lycianthes_rantonnetii", "Capsicum_chacoense", "Nectouxia_formosa", "Salpichroa_tristis", "Solanum_annuum", "Solanum_salamancae", "Solanum_weddellii", "Solanum_sinuatirecurvum", "Solanum_corymbosum", "Solanum_crispum", "Solanum_uncinellum", "Solanum_valdiviense", "Solanum_terminale", "Solanum_herculeum", "Solanum_laciniatum", "Solanum_peruvianum", "Solanum_pimpinellifolium", "Solanum_palustre", "Solanum_microdontum", "Solanum_caripense", "Solanum_evolvulifolium", "Solanum_trizygum", "Solanum_pinnatum", "Solanum_montanum", "Solanum_mahoriense", "Solanum_sejunctum", "Solanum_barbisetum", "Solanum_oedipus", "Solanum_hieronymi", "Solanum_jamaicense", "Solanum_torvum", "Solanum_sisymbriifolium", "Solanum_crinitum", "Solanum_polygamum", "Solanum_evolvuloides", "Solanum_thelopodium", "Jaltomata_bicolor", "Jaltomata_grandiflora", "Datura_ferox", "Brugmansia_suaveolens", "Trompettia_cardenasiana", "Hawkesiophyton_panamense", "Juanulloa_mexicana", "Merinthopodium_neuranthum", "Markea_coccinea", "Dyssochroma_viridiflorum", "Solandra_grandiflora", "Schultesianthus_leucanthus", "Exodeconus_miersii", "Nicandra_physalodes", "Mandragora_officinarum", "Przewalskia_tangutica", "Physochlaina_orientalis", "Scopolia_carniolica", "Hyoscyamus_niger", "Atropanthe_sinensis", "Anisodus_tanguticus", "Atropa_belladonna", "Jaborosa_laciniata", "Sclerophylax_spinescens", "Nolana_divaricata", "Lycium_americanum", "Latua_pubiflora", "Duboisia_myoporoides", "Cyphanthera_anthocercidea", "Crenidium_spinescens", "Grammosolen_dixonii", "Anthotroche_pannosa", "Anthocercis_intricata", "Nicotiana_wigandioides", "Symonanthus_aromaticus")


# make a vector with all the clades

clades <- list(Solanaceae, Solanoideae, Solanum, Dodecachroma, Datureae, Cestroideae, Capsiceae, Physalideae)
names(clades) <- c("Solanaceae", "Solanoideae", "Solanum", "Dodecachroma", "Datureae", "Cestroideae", "Capsiceae", "Physalideae")


#---------------------------------------
# Get the samples for time homogeneous
#----------------------------------------

# get a sample of trees
trees_tHom = sample(tHom, 100)

df_tHom = data.frame(row.names = c(1:length(trees_tHom)))

# for each clade
for (this_clade in 1:length(clades)){
  
  # make an empty vector
  clade_ages_tHom = vector(length = length(trees_tHom))
  
  # for each tree in the sample
  for (t in 1:length(trees_tHom)){
    
    # find the common ancestor of this clade
    this_mrca     = getMRCA(trees_tHom[[t]], clades[[this_clade]])
    this_subtree  = extract.clade(trees_tHom[[t]], node = this_mrca)
    
    # and get its age
    ages     = dispRity::tree.age(this_subtree)
    this_age = max(ages$ages)
    
    # and save it 
    clade_ages_tHom[t] = this_age 
  }
  
  # save in df
  df_tHom[ , this_clade ]        = clade_ages_tHom
  colnames(df_tHom)[ this_clade] = names(clades[this_clade])
}

#---------------------------------------
# Get the samples for time heterogeneous
#----------------------------------------

# get a sample of trees
trees_tHet = sample(tHet, 100)

# make empty df
df_tHet = data.frame(row.names = c(1:length(trees_tHet)))

# for each clade
for (this_clade in 1:length(clades)){
  
  # make an empty vector
  clade_ages_tHet = vector(length = length(trees_tHet))
  
  # for each tree in the sample
  for (t in 1:length(trees_tHet)){
    
    # find the common ancestor of this clade
    this_mrca     = getMRCA(trees_tHet[[t]], clades[[this_clade]])
    this_subtree  = extract.clade(trees_tHet[[t]], node = this_mrca)
    
    # and get its age
    ages     = dispRity::tree.age(this_subtree)
    this_age = max(ages$ages)
    
    # and save it 
    clade_ages_tHet[t] = this_age 
  }
  
  # save in df
  df_tHet[ , this_clade ]        = clade_ages_tHet
  colnames(df_tHet)[this_clade]  = names(clades[this_clade])
}



# combine into a single table and expand for plot formatting

# add model column
df_tHet$model = "Time heterogeneous"
df_tHom$model = "Time homogeneous"

# combine tables
df = rbind(df_tHom, df_tHet)

# get it in long format
df_long = pivot_longer(df,
                       cols = c(1:(ncol(df)-1)),
                       names_to = "clade",
                       values_to = "age")

df_long$clade = as.factor(df_long$clade)


# plot
plot <- ggplot(df_long, aes(x = age, y = clade, fill = model)) +
  geom_violin(position = position_dodge(width = 1)) +
  stat_summary(
    fun = mean, geom = "point", size = 3,
    position = position_dodge(width = 1),
    aes(color = model)
  ) +
  theme_bw() + 
  xlab("Age (Mya)") +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  ) +
  scale_fill_viridis(discrete = TRUE) +
  scale_color_manual(values = c("white", "black")) +  
  facet_grid(rows = vars(clade), scales = "free", switch = "y")

plot




