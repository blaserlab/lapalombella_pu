blaseRtemplates::project_data("~/network/X/Labs/Blaser/share/collaborators/lapalombella_pu_network/datapkg")
datadir <- fs::path("../lapalombella.pu.datapkg/data")
extdatadir <- fs::path("../lapalombella.pu.datapkg/inst/extdata")
blaseRtemplates::save_monocle_disk(cds_disk = cds_main,
                                   data_directory = datadir,
                                   extdata_directory = extdatadir)

blaseRtemplates::save_monocle_disk(cds_disk = cds_WT_AML_bALL,
                                   data_directory = datadir,
                                   extdata_directory = extdatadir)

blaseRtemplates::save_monocle_disk(cds_disk = cds_wt_aml_marrow,
                                   data_directory = datadir,
                                   extdata_directory = extdatadir)

blaseRtemplates::save_monocle_disk(cds_disk = cds_p568,
                                   data_directory = datadir,
                                   extdata_directory = extdatadir)

blaseRtemplates::save_monocle_disk(cds_disk = muench_cds,
                                   data_directory = datadir,
                                   extdata_directory = extdatadir)
