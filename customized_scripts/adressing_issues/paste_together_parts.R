
library(beepr)

if (exists('LOCAL')) {
    base_dir = '/Users/m.wehrens/Data/_2019_02_HCM_SCS/2021_HPC_analysis/'
    script_dir = '/Users/m.wehrens/Documents/git_repos/SCS_More_analyses/'
    data_dir2 = '/Volumes/fastq_m.wehrens/Mapping/WANG2/counttables/'
} else {
    # HPC
    script_dir = '/hpc/hub_oudenaarden/mwehrens/scripts/SCS_HCM_analysis/'
    base_dir = '/hpc/hub_oudenaarden/mwehrens/data/Wang/'
    data_dir2='/hpc/hub_oudenaarden/mwehrens/fastq/WANG2/mapping_may25/counttables/'
}


source(paste0(script_dir,'Functions/Load-Pool-Scale_Simple_MW.R'))


# small script to paste together data that was mapped separately

partlist_string = paste0('part.',0:9,'.')
file_list_parts = paste0(data_dir2, partlist_string, 'GSM2970358_N2_LV_cat_nc_uniaggGenes_spliced.UFICounts.tsv')
names(file_list_parts) = partlist_string

SCS_df_list_data_raw = loadData_MW_parallel(file_list_parts, mc.cores=MYMCCORES, prefix=F)
pryr::object_size(SCS_df_list_data_raw)

SCS_df_list_data_raw_pooled = pool_df_mw(SCS_df_list_data_raw)
beepr::beep()

    # this takes ages, ±15 mins/ dataset (each ±100mb) 
    # so expected time ±150 mins or 2.5 hr

############################################################
# This can also be done using Seurat

if (!exists('SCS_df_list_data_raw')) {
    SCS_df_list_data_raw = loadData_MW_parallel(file_list_parts, mc.cores=MYMCCORES, prefix=F)
}

# simple name replace (underscores not allowed)
SCS_df_list_data_raw_renamed = 
    mclapply(SCS_df_list_data_raw, function(x) {rownames(x) = gsub(pattern = '_',replacement = '==',x = rownames(x)); return(x)})

Seurat_object_list_parts = 
    mclapply(1:length(SCS_df_list_data_raw_renamed), 
                function(idx) {
                    object=CreateSeuratObject(counts = SCS_df_list_data_raw_renamed[[idx]], project = names(SCS_df_list_data_raw_renamed)[idx])
                    print(paste0(names(SCS_df_list_data_raw_renamed)[idx],' done .'))
                    return(object)
                    }, mc.cores = MYMCCORES)

Seurat_object_merged <- merge(Seurat_object_list_parts[[1]], y = unlist(Seurat_object_list_parts)[2:length(Seurat_object_list_parts)], add.cell.ids = names(Seurat_object_list_parts), project = "W_N2")

# extract table
merged_data_table = Seurat_object_merged@assays$RNA@counts

# fix names back again
rownames(merged_data_table) = gsub(pattern = '==',replacement = '_',x = rownames(merged_data_table))
    # rownames(merged_data_table)[1:10]

# now export the table
write.table(x=merged_data_table, file = paste0(data_dir2,'merged.','GSM2970358_N2_LV_cat_nc_uniaggGenes_spliced.UFICounts.tsv'))
    # This led to an error, since the table was too large
    # In fact, the table had 1 mln rows, so it indeed was large

##########
# Let's try to clean up the rows first then ..

SCS_df_list_data = lapply(SCS_df_list_data_raw, preprocess_convertAAnames_toSymbol, revert_to_hgnc=T, script_dir=script_dir)
    # dimensions seem more sensible now

# Let's try again
Seurat_object_list_parts = 
    mclapply(1:length(SCS_df_list_data), 
                function(idx) {
                    object=CreateSeuratObject(counts = SCS_df_list_data[[idx]], project = names(SCS_df_list_data)[idx])
                    print(paste0(names(SCS_df_list_data)[idx],' done .'))
                    return(object)
                    }, mc.cores = MYMCCORES)

Seurat_object_merged <- merge(Seurat_object_list_parts[[1]], y = unlist(Seurat_object_list_parts)[2:length(Seurat_object_list_parts)], add.cell.ids = names(Seurat_object_list_parts), project = "W_N2")
    # dim(Seurat_object_merged)
    # [1] 49993  1240

merged_data_table = Seurat_object_merged@assays$RNA@counts

write.table(x=merged_data_table, file = paste0(data_dir2,'merged.','GSM2970358_N2_LV_cat_nc_uniaggGenes_spliced.UFICounts.tsv'))

