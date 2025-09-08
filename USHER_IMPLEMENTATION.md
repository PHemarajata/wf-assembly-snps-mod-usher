# UShER-Based Tree Grafting Implementation

## Overview

This implementation replaces the previous `gotree`-based tree grafting approach with **UShER (Ultrafast Sample placement on Existing tRees)** for phylogenetic placement. This provides several advantages:

- **No label matching required**: Placement is based on genotypes, not tip labels
- **Backbone topology preserved**: UShER strictly maintains the backbone tree structure
- **Scalable**: Designed for very large datasets
- **Deterministic**: Reproducible results based on parsimony placement
- **Robust**: Handles sample naming inconsistencies gracefully

## Implementation Details

### New Modules

1. **`INTEGRATE_CORE_SNPS`** (`modules/local/integrate_core_snps/main.nf`)
   - Integrates core SNPs from all cluster alignments into a single alignment
   - Creates sample-to-cluster mapping
   - Handles missing data appropriately

2. **`PLACE_WITH_USHER`** (`modules/local/place_with_usher/main.nf`)
   - Converts integrated alignment to VCF format using `snp-sites`
   - Splits samples into backbone vs non-backbone subsets
   - Builds mutation-annotated tree (MAT) preserving backbone topology
   - Places non-backbone samples using maximum parsimony
   - Exports final tree in Newick format

### Workflow Changes

The recombination-aware workflow (`workflows/recombination_aware_snps.nf`) has been updated:

**Previous approach (Steps 7):**
```
7. Graft subtrees onto backbone (gotree-based)
```

**New approach (Steps 7-8):**
```
7. Integrate core SNPs from all clusters
8. UShER-based phylogenetic placement
```

### Configuration Parameters

New parameters in `conf/params.config`:

```groovy
// UShER-based phylogenetic placement parameters
use_usher_placement               = true     // Use UShER instead of gotree
usher_retain_branch_lengths       = true     // Retain input branch lengths
usher_optimize_final_tree         = true     // Optimize final tree after placement
usher_container                   = null     // Override UShER container if needed
```

## Algorithm Flow

1. **Clustering**: Genomes are clustered using Mash distances
2. **Per-cluster analysis**: Each cluster undergoes:
   - Whole genome alignment (Snippy/Parsnp)
   - Recombination detection (Gubbins)
   - Phylogenetic reconstruction (IQ-TREE)
3. **Representative selection**: One representative per cluster
4. **Backbone tree**: Built from cluster representatives
5. **Core SNPs integration**: Variable positions from all clusters combined
6. **VCF generation**: Integrated alignment converted to VCF format
7. **MAT construction**: Mutation-annotated tree built with backbone topology
8. **Sample placement**: Non-backbone samples placed using parsimony
9. **Tree export**: Final tree exported in Newick format

## Advantages Over Previous Approach

### Genotype-Based Placement
- No need for complex label matching and normalization
- Robust to sample naming inconsistencies
- Placement based on actual genetic variation

### Topology Preservation
- Backbone tree structure is strictly maintained
- Branch lengths can be preserved if desired
- No risk of topology distortion during grafting

### Scalability
- UShER is designed for very large phylogenies
- Efficient parsimony placement algorithm
- Memory-efficient mutation-annotated tree format

### Reproducibility
- Deterministic placement based on parsimony
- No dependency on label matching heuristics
- Consistent results across runs

## Requirements

### Software Dependencies
- **UShER** (≥0.6.7): Phylogenetic placement
- **matUtils** (≥0.6.7): MAT manipulation
- **bcftools** (≥1.18): VCF processing
- **snp-sites** (≥2.5.1): VCF generation from alignment
- **gotree** (≥0.4.4): Tree manipulation utilities

### Input Requirements
- **Reference genome**: Required for VCF generation and MAT annotation
- **Integrated alignment**: Core SNPs from all clusters
- **Backbone tree**: Tree of cluster representatives
- **Sample mapping**: Cluster assignments for all samples

## Usage

The UShER-based approach is automatically used in the recombination-aware mode:

```bash
nextflow run . --recombination_aware_mode \
    --input assemblies/ \
    --ref reference.fasta \
    --outdir results \
    -profile docker
```

## Output Files

### Final Results (`Final_Results/`)
- `global_grafted.treefile`: Final phylogenetic tree with all samples
- `grafting_report.txt`: Detailed placement report
- `grafting_log.txt`: Processing log

### Integrated Results (`Integrated_Results/`)
- `integrated_core_snps.fa`: Combined core SNPs alignment
- `sample_cluster_mapping.tsv`: Sample-to-cluster assignments
- `integration_report.txt`: Integration statistics

## Troubleshooting

### Common Issues

1. **No backbone samples in VCF**: Ensure sample names match between backbone tree and alignment
2. **UShER placement fails**: Check VCF format and reference genome compatibility
3. **Empty integrated alignment**: Verify cluster alignments contain variable positions

### Debug Information

The `grafting_log.txt` file contains detailed information about:
- Sample counts and assignments
- VCF generation and processing
- MAT construction progress
- Placement statistics
- Error messages and warnings

## Performance Considerations

- **Memory usage**: Scales with number of variable positions
- **Runtime**: Linear with number of samples to place
- **Disk space**: VCF files can be large for many samples
- **CPU usage**: Parallelizable across available cores

## References

- **UShER**: Turakhia et al. "Ultrafast Sample placement on Existing tRees (UShER) enables real-time phylogenetics for the SARS-CoV-2 pandemic." Nature Genetics (2021)
- **Mutation-annotated trees**: McBroome et al. "A daily-updated database and tools for comprehensive SARS-CoV-2 mutation-annotated trees." Molecular Biology and Evolution (2021)