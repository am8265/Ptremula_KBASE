library(DESeq2)
library(gage)
library(pathview)

setwd("/Users/ehu/UCKalluri/KBASE/Ptremula/Ptrem_pathway_Analysis")
#Ctrl_Axillary_bud_IAA752_Axillary_bud
#Ctrl_shoot vs Ctrl_Shoot_IAA7_52
#CtrlA_AxillaryBud_717c_1, CtrlA_AxillaryBud_717c_2, CtrlA_AxillaryBud_717c_3
#IAA752_AxillaryBud_2, IAA752_AxillaryBud_3

#CtrlA_ShootTip_717c_1", "CtrlA_ShootTip_717c_2", "CtrlA_ShootTip_717c_3",
#IAA7_54_ShootTip_1, IAA7_54_ShootTip_2, IAA7_54_ShootTip_3

cnts = read.csv("Ctrl_Axillary_bud_IAA754_Axillary_bud/gene_count_matrix.csv", header = T)

cnts = read.csv("Ctrl_Shoot_tip_IAA754_Shoot_tip/gene_count_matrix.csv", header = T)
rownames(cnts) = cnts$gene_id
cnts$gene_id <- NULL
colnames(cnts) <- c("CtrlA_ShootTip_717c_1", "CtrlA_ShootTip_717c_2", 
                    "CtrlA_ShootTip_717c_3",
                    "IAA7_54_ShootTip_1", "IAA7_54_ShootTip_2", "IAA7_54_ShootTip_3")

colnames(cnts) <- c()

colnames(cnts) <- c("CtrlA_AxillaryBud_717c_1", "CtrlA_AxillaryBud_717c_2", "CtrlA_AxillaryBud_717c_3",
                    "IAA754_AxillaryBud_1", "IAA752_AxillaryBud_2")
colnames(cnts) <- c("CtrlA_ShootTip_717c_1", "CtrlA_ShootTip_717c_2", "CtrlA_ShootTip_717c_3",
                   "IAA752_ShootTip_2", "IAA752_ShootTip_3")
dim(cnts)

#removing 'All 0' counts from the dataset
cnts = cnts[rowSums(cnts) != 0,]
dim(cnts)
cnts <- as.matrix(cnts)

coldat = data.frame(condition = factor(c(rep("CtrlA_ShootTip",3), rep("IAA754_ShootTip",3))))

coldat = data.frame(condition = factor(c(rep("CtrlA_AxillaryBud",3), rep("IAA754_AxillaryBud",2))))

coldat = data.frame(condition = factor(c(rep("CtrlA_ShootTip",3), rep("IAA752_ShootTip",2))))

dds = DESeqDataSetFromMatrix(cnts, colData = coldat, design = ~ condition) 
dds <- DESeq(dds)
deseq2.res <- results(dds)
deseq2.fc = deseq2.res$log2FoldChange

#read the mapping table 'map_gids.txt' for Poplar geneIDs conversion
map_gids = read.table("map_gids.txt")
rownames(map_gids) = map_gids$V1
map_gids$V1 <- NULL
saveRDS(map_gids, "map_gids.rds")


gids_mapped = as.character(map_gids[rownames(deseq2.res),1])
names(deseq2.fc) <- gids_mapped
#remove 'NA' names in deseq2.fc ..NA names appear where JGI ids could not be mapped!
deseq2.fc <- deseq2.fc[!is.na(names(deseq2.fc))]
names(deseq2.fc) = paste0(names(deseq2.fc),"g")
exp.fc <- deseq2.fc
out.suffix <- "deseq2"

kegg_ptrem = kegg.gsets(species = "pop",check.new = T)
kegg_ptrem_sigmet = kegg_ptrem$kg.sets[kegg_ptrem$sigmet.idx]
kegg.gs = kegg_ptrem_sigmet

kegg.gs = readRDS("kegg.pop.sigmet.gsets.rds")
fc.kegg.p <- gage(exp.fc, gsets = kegg.gs, ref = NULL, samp = NULL)


##using DESEQ2 variance stablised data-VSD
vsd <- vst(dds, blind=FALSE)
ref_idx = c(1:3)
samp_idx = c(4:5)
gids_mapped = as.character(map_gids[rownames(vsd),1])
vsd_mat = assay(vsd)
rownames(vsd_mat) = gids_mapped
vsd_mat = vsd_mat[!is.na(rownames(vsd_mat)),]


vsd.cnts.kegg.p <- gage(vsd_mat, gsets = kegg.gs, ref = ref_idx, samp = samp_idx,
                        compare = "unpaired", )



##Let's process using raw counts using GAGE native workflow

cnts = assay(dds)
libsizes = colSums(cnts)
size.factor = libsizes/exp(mean(log(libsizes)))
cnts.norm=t(t(cnts)/size.factor)
range(cnts.norm)
cnts.norm=log2(cnts.norm+8)
range(cnts.norm)

ref_idx = c(1:3)
samp_idx = c(4:5)
gids_mapped = as.character(map_gids[rownames(vsd),1])
rownames(cnts.norm) = gids_mapped
cnts.norm = cnts.norm[!is.na(rownames(cnts.norm)),]
kegg.ptrem = kegg.gsets(species = "pop", check.new = T)
kegg.ptrem.sigmet = kegg.ptrem$kg.sets[kegg.ptrem$sigmet.idx]
kegg.gs = kegg.ptrem.sigmet
cnts.norm.kegg.p <- gage(cnts.norm, gsets = kegg.gs, ref = ref_idx, samp = samp_idx,
                        compare = "unpaired")



##PathView Analysis

sel.g = fc.kegg.p$greater[,"q.val"] < 0.1 & !is.na(fc.kegg.p$greater[,"q.val"])
sum(sel.g)
path.ids = rownames(fc.kegg.p$greater)[sel.g]



sel.l = fc.kegg.p$less[,"q.val"] < 0.1 & !is.na(fc.kegg.p$less[,"q.val"])
sum(sel.l)
path.ids.l = rownames(fc.kegg.p$less)[sel.l]

path.ids2 = substr(c(path.ids,path.ids.l),1, 8)


kegg_dir = "Shoot_54"
pv.out.list <- sapply(path.ids2, 
                      function(pid) 
                        pathview(gene.data = exp.fc, pathway.id = pid, 
                                 species = "pop", gene.idtype = "kegg", 
                                 out.suffix = out.suffix, kegg.dir = kegg_dir))

