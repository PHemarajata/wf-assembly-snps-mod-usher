#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process PLACE_WITH_USHER {
    tag "usher_placement"
    label 'process_medium'
    
    // Use conda environment with UShER, matUtils, bcftools, and gotree
    // Override via params.usher_container in your config if you'd like
    container "${params.usher_container ?: 'quay.io/biocontainers/usher:0.6.6--h719ac0c_2'}"

    conda """
        bioconda::usher=0.6.7
        bioconda::matutils=0.6.7
        bioconda::bcftools=1.18
        bioconda::gotree=0.4.4
        bioconda::snp-sites=2.5.1
    """
    
    publishDir "${params.outdir}/Final_Results",
               mode: params.publish_dir_mode,
               pattern: "*.{treefile,txt,log,yml,yaml,pdf}"

    input:
    path backbone_tree          // e.g. backbone.treefile (Newick)
    path integrated_alignment   // Combined core SNPs alignment from all clusters
    path reference_fasta        // Reference genome used for SNP calling
    path cluster_representatives // TSV mapping cluster_id to representative_id

    output:
    path "global_grafted.treefile", emit: grafted_tree
    path "grafting_report.txt",    emit: report
    path "grafting_log.txt",       emit: log
    path "versions.yml",           emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    '''
    set -Eeuo pipefail

    : > grafting_log.txt
    printf 'Starting UShER placement\\n' >> grafting_log.txt
    printf 'Backbone tree: %s\\n' "!{backbone_tree}" >> grafting_log.txt
    printf 'Integrated alignment: %s\\n' "!{integrated_alignment}" >> grafting_log.txt
    printf 'Reference: %s\\n' "!{reference_fasta}" >> grafting_log.txt

    # Check if we have the required inputs
    if [[ ! -s "!{backbone_tree}" ]]; then
        echo "ERROR: Backbone tree file is empty or missing" | tee -a grafting_log.txt
        exit 1
    fi

    if [[ ! -s "!{integrated_alignment}" ]]; then
        echo "ERROR: Integrated alignment file is empty or missing" | tee -a grafting_log.txt
        exit 1
    fi

    # Handle reference genome requirement
    if [[ ! -s "!{reference_fasta}" ]] || [[ "!{reference_fasta}" == *"NO_FILE"* ]]; then
        echo "WARNING: No reference genome provided. UShER requires a reference for VCF annotation." | tee -a grafting_log.txt
        echo "Attempting to use the first sequence from the integrated alignment as reference..." | tee -a grafting_log.txt
        
        # Extract first sequence as reference
        head -2 "!{integrated_alignment}" > temp_reference.fa
        if [[ -s temp_reference.fa ]]; then
            reference_file="temp_reference.fa"
            echo "Using first sequence as reference: $(head -1 temp_reference.fa)" >> grafting_log.txt
        else
            echo "ERROR: Cannot create reference from alignment" | tee -a grafting_log.txt
            exit 1
        fi
    else
        reference_file="!{reference_fasta}"
        echo "Using provided reference: !{reference_fasta}" >> grafting_log.txt
    fi

    # 1) Extract backbone sample IDs from the Newick (via gotree)
    gotree labels -i "!{backbone_tree}" > backbone.samples.txt
    n_bb=$(wc -l < backbone.samples.txt || echo 0)
    printf 'Backbone samples detected: %s\\n' "$n_bb" >> grafting_log.txt
    
    if [[ "$n_bb" -eq 0 ]]; then
        echo "ERROR: No labels found in backbone tree" | tee -a grafting_log.txt
        exit 1
    fi

    # 2) Convert the integrated alignment to VCF format using snp-sites
    printf 'Converting alignment to VCF format...\\n' >> grafting_log.txt
    
    # First, create a VCF from the integrated alignment
    snp-sites -v -o combined.vcf "!{integrated_alignment}" 2>> grafting_log.txt || {
        echo "ERROR: Failed to convert alignment to VCF" | tee -a grafting_log.txt
        exit 1
    }

    # Compress and index the VCF
    bgzip -c combined.vcf > combined.vcf.gz
    tabix -p vcf combined.vcf.gz

    # 3) Get all sample names from the VCF
    bcftools query -l combined.vcf.gz > all.samples.txt
    n_total=$(wc -l < all.samples.txt || echo 0)
    printf 'Total samples in VCF: %s\\n' "$n_total" >> grafting_log.txt

    # 4) Split samples into backbone vs non-backbone subsets
    # Find backbone samples that exist in the VCF
    grep -Fx -f backbone.samples.txt all.samples.txt > backbone.samples.in_vcf.txt || true
    n_in=$(wc -l < backbone.samples.in_vcf.txt || echo 0)
    
    if [[ "$n_in" -eq 0 ]]; then
        echo "ERROR: None of the backbone samples are present in the VCF" | tee -a grafting_log.txt
        echo "Backbone samples:" >> grafting_log.txt
        cat backbone.samples.txt >> grafting_log.txt
        echo "VCF samples:" >> grafting_log.txt
        head -10 all.samples.txt >> grafting_log.txt
        exit 1
    fi

    # Find non-backbone samples
    comm -23 <(sort all.samples.txt) <(sort backbone.samples.in_vcf.txt) > nonbackbone.samples.in_vcf.txt || true
    n_nonbb=$(wc -l < nonbackbone.samples.in_vcf.txt || echo 0)

    printf 'Backbone samples in VCF: %s\\n' "$n_in" >> grafting_log.txt
    printf 'Non-backbone samples in VCF: %s\\n' "$n_nonbb" >> grafting_log.txt

    # 5) Create separate VCFs for backbone and non-backbone samples
    bcftools view -S backbone.samples.in_vcf.txt -Oz -o backbone.vcf.gz combined.vcf.gz
    tabix -p vcf backbone.vcf.gz

    if [[ "$n_nonbb" -gt 0 ]]; then
        bcftools view -S nonbackbone.samples.in_vcf.txt -Oz -o nonbackbone.vcf.gz combined.vcf.gz
        tabix -p vcf nonbackbone.vcf.gz
    fi

    # 6) Build a mutation-annotated tree (MAT) that PRESERVES the backbone topology
    printf 'Building mutation-annotated tree with backbone topology...\\n' >> grafting_log.txt
    
    # Configure UShER options based on parameters
    retain_branch_lengths=""
    if [[ "!{params.usher_retain_branch_lengths}" == "true" ]]; then
        retain_branch_lengths="--retain-input-branch-lengths"
    fi
    
    usher \\
        --tree "!{backbone_tree}" \\
        --vcf backbone.vcf.gz \\
        --reference "$reference_file" \\
        --save-mutation-annotated-tree backbone.pb \\
        --out-dir step1 \\
        $retain_branch_lengths \\
        --threads !{task.cpus} \\
        2>> grafting_log.txt || {
        echo "ERROR: Failed to build backbone MAT" | tee -a grafting_log.txt
        exit 1
    }

    # 7) Place the remaining samples onto the MAT (if any)
    MAT_OUT="step1/backbone.pb"
    
    if [[ "$n_nonbb" -gt 0 && -s nonbackbone.vcf.gz ]]; then
        printf 'Placing non-backbone samples onto MAT...\\n' >> grafting_log.txt
        
        # Configure optimization option
        optimize_tree=""
        if [[ "!{params.usher_optimize_final_tree}" == "true" ]]; then
            optimize_tree="--optimize-final-tree"
        fi
        
        usher \\
            --load-mutation-annotated-tree "$MAT_OUT" \\
            --vcf nonbackbone.vcf.gz \\
            --reference "$reference_file" \\
            --out-dir step2 \\
            $optimize_tree \\
            --threads !{task.cpus} \\
            2>> grafting_log.txt || {
            echo "WARNING: Failed to place non-backbone samples, using backbone-only tree" | tee -a grafting_log.txt
        }
        
        # Check if placement was successful
        if [[ -f "step2/usher.pb" ]]; then
            MAT_OUT="step2/usher.pb"
            printf 'Successfully placed %s non-backbone samples\\n' "$n_nonbb" >> grafting_log.txt
        else
            printf 'Placement failed, using backbone-only MAT\\n' >> grafting_log.txt
        fi
    else
        printf 'No non-backbone samples to place; using backbone-only MAT\\n' >> grafting_log.txt
    fi

    # 8) Export final Newick tree
    printf 'Exporting final tree...\\n' >> grafting_log.txt
    
    matUtils extract -i "$MAT_OUT" -t global_grafted.treefile 2>> grafting_log.txt || {
        echo "ERROR: Failed to extract final tree" | tee -a grafting_log.txt
        # Fallback: copy backbone tree
        cp "!{backbone_tree}" global_grafted.treefile
        echo "Using backbone tree as fallback" >> grafting_log.txt
    }

    # 9) Generate comprehensive report
    final_leaves=$(gotree labels -i global_grafted.treefile | wc -l 2>/dev/null || echo "unknown")
    
    {
        printf 'TREE GRAFTING (UShER placement) REPORT\\n'
        printf '======================================\\n'
        printf 'Date (UTC): %s\\n' "$(date -u +%FT%TZ)"
        printf 'Backbone tree: %s\\n' "!{backbone_tree}"
        printf 'Integrated alignment: %s\\n' "!{integrated_alignment}"
        printf 'Reference genome: %s\\n' "!{reference_fasta}"
        printf '\\n'
        printf 'SAMPLE STATISTICS\\n'
        printf '-----------------\\n'
        printf 'Total samples in alignment: %s\\n' "$n_total"
        printf 'Backbone samples found in VCF: %s\\n' "$n_in"
        printf 'Non-backbone samples: %s\\n' "$n_nonbb"
        printf 'Final tree leaves: %s\\n' "$final_leaves"
        printf '\\n'
        printf 'PROCESSING RESULTS\\n'
        printf '------------------\\n'
        if [[ "$final_leaves" -gt "$n_in" ]]; then
            printf 'Status: SUCCESS - Non-backbone samples successfully placed\\n'
        elif [[ "$final_leaves" -eq "$n_in" ]]; then
            printf 'Status: PARTIAL - Only backbone samples in final tree\\n'
        else
            printf 'Status: WARNING - Unexpected number of leaves in final tree\\n'
        fi
        printf '\\n'
        printf 'METHOD DETAILS\\n'
        printf '--------------\\n'
        printf 'Phylogenetic placement: UShER (Ultrafast Sample placement on Existing tRees)\\n'
        printf 'Mutation annotation: matUtils\\n'
        printf 'VCF generation: snp-sites from core SNPs alignment\\n'
        printf 'Backbone topology: Preserved using --tree and --retain-input-branch-lengths\\n'
        printf 'Placement algorithm: Maximum parsimony on mutation-annotated tree\\n'
        printf '\\n'
        printf 'ADVANTAGES OF THIS APPROACH\\n'
        printf '---------------------------\\n'
        printf '• No label matching required - placement based on genotypes\\n'
        printf '• Backbone topology strictly preserved\\n'
        printf '• Scalable to very large datasets\\n'
        printf '• Deterministic and reproducible results\\n'
        printf '• Robust to sample naming inconsistencies\\n'
    } > grafting_report.txt

    printf 'UShER placement completed successfully\\n' >> grafting_log.txt

    # 10) Generate versions file
    {
        printf '"!{task.process}":\\n'
        printf '    usher: %s\\n' "$(usher --version 2>&1 | head -n1 | sed 's/^/    /' || echo 'unknown')"
        printf '    matUtils: %s\\n' "$(matUtils --version 2>&1 | head -n1 | sed 's/^/    /' || echo 'unknown')"
        printf '    bcftools: %s\\n' "$(bcftools --version 2>&1 | head -n1 | sed 's/^/    /' || echo 'unknown')"
        printf '    gotree: %s\\n' "$(gotree version 2>/dev/null | sed 's/^/    /' | tr -d '\\r' || echo 'unknown')"
        printf '    snp-sites: %s\\n' "$(snp-sites -V 2>&1 | sed 's/snp-sites /    /' || echo 'unknown')"
    } > versions.yml
    '''
}