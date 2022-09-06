#!/bin/bash

# Run all local jobs

# directory with scripts
currentpath=/Users/m.wehrens/Documents/git_repos/mapping_aa_private_vasa/customized_scripts/

# Run the main script for each sample
${currentpath}L_TS_submit_vasaplate_map_LOCAL.sh MW-TS-S1-Hybr-NB_HTMH2BGXF_S23_cat HUMAN 74 XXX
${currentpath}L_TS_submit_vasaplate_map_LOCAL.sh MW-TS-S1A-Hybr-Bst_HTMH2BGXF_S25_cat HUMAN 74 XXX
${currentpath}L_TS_submit_vasaplate_map_LOCAL.sh MW-TS-S2A-TargOnly-Bst_HTMH2BGXF_S26_cat HUMAN 74 XXX
${currentpath}L_TS_submit_vasaplate_map_LOCAL.sh MW-TS-S2-TargOnly-NB_HTMH2BGXF_S24_cat HUMAN 74 XXX
${currentpath}L_TS_submit_vasaplate_map_LOCAL.sh MW-TS-S4A-Vasa-Bst_HTMH2BGXF_S27_cat HUMAN 74 XXX
${currentpath}L_TS_submit_vasaplate_map_LOCAL.sh MW-TS-S4-Vasa-NB_HTMH2BGXF_S28_cat HUMAN 74 XXX

