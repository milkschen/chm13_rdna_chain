# Chain files for the CHM13 rDNA flanking sequences

1. Prepare CHM13-v1.0 and CHM13-v1.1 FASTAs and the rDNA flanking region sequences wrt to each reference version (see `rdna_to_chain.sh` for details).
2. Run `rdna_to_chain.sh` to build the chain file, including rDNA flanking regions, from CHM13-v1.0 to CHM13-v1.1 and CHM13-v1.1 to CHM13-v1.0.

Outputs:
- v1.0_to_v1.1_rdna_merged.chain
- v1.1_to_v1.0_rdna_merged.chain
