allen_prep = function(recompute = F) {
    if (recompute) {
        load("./0_DATA/SCRNA/Seurat.ss.rda")
        meta = read.csv("./0_DATA/SCRNA/annotation 20200913-Table 1.csv",
            header = T)
        my_ref = subset(ss.seurat, region_label %in% c("ALM", "SSp", "VISp",
            "SSs", "VIS"))

        ins = subset(my_ref, class_label == "GABAergic")
        exc = subset(my_ref, class_label == "Glutamatergic")
        nn = subset(my_ref, class_label == "Non-Neuronal")

        ins@meta.data$supertype_label = recode(ins@meta.data$cluster_label,
            `114_Pvalb` = "Pvalb Lpl", `43_Vip` = "Vip Pcdh11x", `12_Lamp5` = "Lamp5 Pdlim5",
            `69_Sst` = "Sst Syndig1l", `45_Vip` = "Vip Pcdh11x", `51_Vip` = "Vip Cp Rspo1",
            `97_Sst` = "Sst Hpse", `64_Sst Chodl` = "Sst Chodl", `42_Vip` = "Vip Pcdh11x",
            `16_Lamp5` = "Lamp5 Egln3", `60_Vip Igfbp6` = "Vip Igfbp6",
            `44_Vip` = "Vip Mybpc1", `18_Lamp5` = "Lamp5 Egln3", `59_Vip Igfbp6` = "Vip Igfbp6",
            `46_Vip` = "Vip Pcdh11x", `115_Pvalb` = "Pvalb Lpl", `50_Vip` = "Vip Cp Rspo1",
            `57_Vip Igfbp6` = "Vip Igfbp6", `37_Sncg` = "Sncg Jam2 Npy2r",
            `61_Vip Igfbp6` = "Vip Igfbp6", `11_Lamp5` = "Lamp5 Pdlim5",
            `15_Lamp5` = "Lamp5 Egln3", `109_Pvalb` = "Pvalb Th", `14_Lamp5` = "Lamp5 Egln3",
            `68_Sst` = "Sst Syndig1l", `116_Pvalb` = "Pvalb Lpl", `119_Pvalb` = "Pvalb Lpl",
            `73_Sst` = "Sst Crh", `111_Pvalb` = "Pvalb Lpl", `118_Pvalb` = "Pvalb Lpl",
            `112_Pvalb` = "Pvalb Lpl", `86_Sst` = "Sst Etv1", `84_Sst` = "Sst Myh8",
            `53_Vip` = "Vip Lmo1", `71_Sst` = "Sst Crh", `82_Sst` = "Sst Myh8",
            `35_Sncg` = "Sncg Jam2 Npy2r", `74_Sst` = "Sst Crh", `88_Sst` = "Sst Nmbr",
            `108_Pvalb` = "Pvalb Th", `7_Lamp5 Lhx6` = "Lamp5 Lhx6", `85_Sst` = "Sst Etv1",
            `33_Sncg` = "Sncg Krt73", `13_Lamp5` = "Lamp5 Egln3", `123_Pvalb Vipr2` = "Pvalb Vipr2",
            `19_Pax6` = "Pax6", `62_Vip Igfbp6` = "Vip Igfbp6", `40_Vip` = "Vip Serpinf1",
            `113_Pvalb` = "Pvalb Lpl", `81_Sst` = "Sst Myh8", `22_Sncg` = "Sncg Serpinf1",
            `98_Sst` = "Sst Calb2", `76_Sst` = "Sst Nts", `66_Sst` = "Sst Syndig1l",
            `24_Sncg` = "Sncg Serpinf1", `99_Sst` = "Sst Mme", `100_Sst` = "Sst Mme",
            `79_Sst` = "Sst Myh8", `96_Sst` = "Sst Hpse", `36_Sncg` = "Sncg Jam2 Npy2r",
            `52_Vip` = "Vip Lmo1", `75_Sst` = "Sst Crh", `48_Vip` = "Vip Cp Rspo1",
            `49_Vip` = "Vip Cp Rspo1", `90_Sst` = "Sst Nmbr", `6_Lamp5 Lhx6` = "Lamp5 Lhx6",
            `47_Vip` = "Vip Pcdh11x", `89_Sst` = "Sst Nmbr", `91_Sst` = "Sst Nmbr",
            `92_Sst` = "Sst Nmbr", `83_Sst` = "Sst Myh8", `70_Sst` = "Sst Crh",
            `93_Sst` = "Sst Nmbr", `67_Sst` = "Sst Syndig1l", `94_Sst` = "Sst Nmbr",
            `80_Sst` = "Sst Myh8", `65_Sst Chodl` = "Sst Chodl", `41_Vip` = "Vip Serpinf1",
            `23_Sncg` = "Sncg Serpinf1", `63_Sst Chodl` = "Sst Chodl",
            `117_Pvalb` = "Pvalb Lpl", `58_Vip Igfbp6` = "Vip Igfbp6",
            `31_Sncg` = "Sncg Krt73", `3_Meis2` = "Meis2", `38_Sncg` = "Sncg Jam2 Npy2r",
            `95_Sst` = "Sst Nmbr", `110_Pvalb` = "Pvalb Th", `10_Lamp5` = "Lamp5 Pdlim5",
            `21_Sncg` = "Sncg Serpinf1", `20_Pax6` = "Pax6", `56_Vip HPF` = "Vip Cbln4 HPF",
            `28_Ntng1 HPF` = "Ntng1 HPF", `105_Sst HPF` = "Sst Ctsc HPF",
            `27_Ntng1 HPF` = "Ntng1 HPF", `17_Lamp5` = "Lamp5 Egln3", `2_Meis2` = "Meis2",
            `25_Sncg` = "Sncg Krt73", `121_Pvalb` = "Pvalb Lpl", `72_Sst` = "Sst Crh")

        ins = subset(ins, supertype_label %nin% c("Meis2", "Ntng1 HPF",
            "Sst Ctsc HPF", "Vip Cbln4 HPF"))
        exc = subset(exc, subclass_label %nin% c("CR", "DG", "L2/3 IT PPP",
            "L5/6 IT TPE-ENT", "L6b/CT ENT"))
        exc = subset(exc, region_label %in% c("SSp", "SSs"))

        cmb = subset(my_ref, cells = c(colnames(nn), colnames(exc), colnames(ins)))
        rep = ""
        for (x in unique(cmb@meta.data$cluster_label)) {
            sl = subset(meta, cluster_label == x, supertype_label)
            rep = paste0(rep, x, "\"=\"", sl, "\",\"")
        }
        rep = strtrim(rep, nchar(rep) - 2)

        eval(parse(text = paste0("cmb@meta.data$supertype_label = recode(cmb@meta.data$cluster_label,\"",
            rep, ")")))
        cmb@meta.data$supertype_label = recode(cmb@meta.data$supertype_label,
            `Sst Nts` = "PV/SST Th", `Pvalb Th` = "PV/SST Th", Pax6 = "Lamp5 Pax6")
        cmb = subset(cmb, region_label != "VIS")

        cmb = RunUMAP(FindNeighbors(RunPCA(ScaleData(FindVariableFeatures(NormalizeData(cmb,
            verbose = F), verbose = F), verbose = F), npcs = 50, verbose = F),
            verbose = F), reduction = "pca", dims = 1:50, verbose = F)
    } else {
        cmb = readRDS("./0_DATA/SCRNA/cmb_jun28.RDS")
    }
    return(cmb)
}

fishell_prep = function(recompute = F) {
    if (recompute) {
        P28 = readRDS("./0_DATA/SCRNA/CX_P28_filtered2_05312022.rds")
        P28 = FindClusters(P28, resolution = 2.1, verbose = F)
    } else {
        P28 = readRDS("./0_DATA/SCRNA/P28_jun28.RDS")
    }
    return(P28)
}

data_alignement = function(cmb, P28, recompute = F) {
    if (recompute) {
        smrt_in = subset(cmb, class_label == "GABAergic")
        smrt_in = RunUMAP(FindNeighbors(RunPCA(ScaleData(FindVariableFeatures(smrt_in,
            verbose = F), verbose = F), npcs = 50, verbose = F), verbose = F),
            reduction = "pca", dims = 1:50, verbose = F)
        anchors = FindIntegrationAnchors(list(P28 = P28, smrt_in = smrt_in),
            dims = 1:100, verbose = F)
        integrated = IntegrateData(anchorset = anchors, dims = 1:100, verbose = F)
        DefaultAssay(integrated) = "integrated"
        integrated = ScaleData(integrated, verbose = FALSE)
        integrated = RunPCA(integrated, npcs = 50, verbose = FALSE)
        integrated = ProjectDim(integrated, verbose = F)
        integrated = FindNeighbors(integrated, verbose = F)
        integrated = FindClusters(integrated, resolution = 1.7, random.seed = 42,
            verbose = F)
        integrated = RunUMAP(integrated, reduction = "pca", dims = 1:50,
            verbose = F)

    } else {
        integrated = readRDS("./0_DATA/SCRNA/integrated_jun28.RDS")
    }
    return(integrated)
}


load_pucks = function(recompute = F) {

    if (recompute) {
        wd = "./0_DATA/PUCKS/A_RAW/"
        pucks = list(F03 = list(dge = paste0(wd, "F03/Puck_200306_03.digital_expression.txt.bz2"),
            bead = paste0(wd, "F03/Puck_200306_03_coords.csv.bz2")), SW6 = list(dge = paste0(wd,
            "SW6/MappedDGEForR.csv.bz2"), bead = paste0(wd, "SW6/BeadLocationsForR.csv.bz2")),
            SW56 = list(dge = paste0(wd, "SW56/Puck_210824_02.matched.digital_expression.txt.gz"),
                bead = paste0(wd, "SW56/Puck_210824_02_barcode_matching.txt.gz")),
            SW55 = list(dge = paste0(wd, "SW55/Puck_210824_01.matched.digital_expression.txt.gz"),
                bead = paste0(wd, "SW55/Puck_210824_01_barcode_matching.txt.gz")),
            SW33 = list(dge = paste0(wd, "SW33/Puck_210605_38.matched.digital_expression.txt.bz2"),
                bead = paste0(wd, "SW33/Puck_210605_38_barcode_matching.txt.bz2")),
            SW30 = list(dge = paste0(wd, "SW30/Puck_210605_35.matched.digital_expression.txt.bz2"),
                bead = paste0(wd, "SW30/Puck_210605_35_barcode_matching.txt.bz2")),
            SW31 = list(dge = paste0(wd, "SW31/Puck_210605_36.matched.digital_expression.txt.bz2"),
                bead = paste0(wd, "SW31/Puck_210605_36_barcode_matching.txt.bz2")))

        spatial_data = lapply(pucks, function(p) {
            return(list(fread(p$dge), fread(p$bead)))
        })
        names(spatial_data) = names(pucks)

        mats = lapply(spatial_data, function(p) {
            setkey(p[[1]], "GENE")
            return(as.matrix(p[[1]], rownames = TRUE))
        })
        names(mats) = names(spatial_data)

        for (p in names(spatial_data)) {
            spatial_data[[p]][[1]] = NULL
        }
        gc()

        for (p in names(spatial_data)) {
            spatial_data[[p]] = data.frame(spatial_data[[p]])
            if ("x" %in% colnames(spatial_data[[p]])) {

            } else {
                if ("V2" %in% colnames(spatial_data[[p]])) {
                  spatial_data[[p]] = spatial_data[[p]][which(!duplicated(spatial_data[[p]]$V2)),
                    ]
                  rownames(spatial_data[[p]]) = spatial_data[[p]]$V2
                  spatial_data[[p]] = spatial_data[[p]][, c("V3", "V4")]
                  colnames(spatial_data[[p]]) = c("x", "y")
                } else {
                  if ("barcodes" %in% colnames(spatial_data[[p]])) {
                    rownames(spatial_data[[p]]) = spatial_data[[p]]$barcodes
                  } else {
                    rownames(spatial_data[[p]]) = spatial_data[[p]]$barcode
                  }
                  spatial_data[[p]] = spatial_data[[p]][, c("xcoord", "ycoord")]
                  colnames(spatial_data[[p]]) = c("x", "y")

                }
            }
        }
        puck_rctd = lapply(names(pucks), function(x) {
            SpatialRNA(spatial_data[[x]], mats[[x]])
        })
        names(puck_rctd) = names(spatial_data)
        
    } else {
        puck_rctd = list()
        pucks = c('F03','SW6','SW56','SW55','SW33','SW30','SW31')
        for(p in pucks){
            puck_rctd[[p]] = readRDS(paste0('./0_DATA/PUCKS/B_PREPROCESSED/rctd_ref_',p,'.RDS'))
            
        }
    }
    return(puck_rctd)    
}

build_pl_df = function(rdf=rdf, spatial=spatial)
{
    a_res = subset(rdf, spot_class %in% c('doublet_certain','singlet'))[,c(1:3)]
    a_res = cbind(spatial[rownames(a_res),],a_res)
    sgs = subset(a_res, spot_class == 'singlet')
    dbs = subset(a_res, spot_class != 'singlet')
    dbs_d = dbs
    dbs_d$first_type = dbs_d$second_type
    rownames(dbs_d) = paste(rownames(dbs),'1',sep = '_')
    sgs = rbind(sgs,rbind(dbs, dbs_d))
    tpl = sgs[,c('x','y','first_type')]
    rm(list=c('dbs','sgs'));gc();

    colnames(tpl) = c('x','y','id')
    tpl[,'my_t'] = ''
    tpl[grep('Sst|SST',tpl$id),'my_t'] = 'SST'
    tpl[grep('PV|Pvalb',tpl$id),'my_t'] = 'PV'
    tpl[grep('Lamp5',tpl$id),'my_t'] = 'Lamp5'
    tpl[grep('Vip',tpl$id),'my_t'] = 'Vip'
    tpl[rownames(subset(tpl, my_t=='')[grep('L',subset(tpl, my_t=='')$id),]),'my_t']='EXC'
    tpl[rownames(subset(tpl, my_t=='')),'my_t'] = 'OTH'
    
    tpl[,'L'] = ''
    tpl[,'L'] = ''
    tpl[rownames(subset(tpl, my_t == 'EXC' & grepl("^L2", id))),'L'] = 'L2'
    tpl[rownames(subset(tpl, my_t == 'EXC' & grepl("^L2_3|^L23", id))),'L'] = 'L2_3'
    tpl[rownames(subset(tpl, my_t == 'EXC' & grepl("^L3", id))),'L'] = 'L3'
    tpl[rownames(subset(tpl, my_t == 'EXC' & grepl("^L4_5|L45", id))),'L'] = 'L4_5'
    tpl[rownames(subset(tpl, my_t == 'EXC' & grepl("^L4", id))),'L'] = 'L4'
    tpl[rownames(subset(tpl, my_t == 'EXC' & grepl("^L5_6|L56", id))),'L'] = 'L5_6'
    tpl[rownames(subset(tpl, my_t == 'EXC' & grepl("^L5", id))),'L'] = 'L5'
    tpl[rownames(subset(tpl, my_t == 'EXC' & grepl("^L6b", id))),'L'] = 'L6b'
    tpl[rownames(subset(tpl, my_t == 'EXC' & grepl("^L6", id))),'L'] = 'L6'

    return(tpl)
}


