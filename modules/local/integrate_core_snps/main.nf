#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process INTEGRATE_CORE_SNPS {
    tag "integrate_snps"
    label 'process_low'
    container "quay.io/biocontainers/python:3.9--1"
    
    publishDir "${params.outdir}/Integrated_Results", mode: params.publish_dir_mode, pattern: "*.{fa,tsv,txt}"

    input:
    path cluster_alignments
    path cluster_info
    path reference_fasta

    output:
    path "integrated_core_snps.fa", emit: integrated_alignment
    path "sample_cluster_mapping.tsv", emit: sample_mapping
    path "integration_report.txt", emit: report
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # Install compatible versions of numpy and pandas
    pip install --upgrade numpy>=1.15.4
    pip install pandas biopython

    python3 << 'EOF'
import os
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict
import glob

def integrate_core_snps():
    # Integrate core SNPs from all cluster alignments into a single alignment
    
    print("Starting core SNPs integration...")
    
    # Find all alignment files
    alignment_files = []
    for pattern in ['*.aln', '*.fa', '*.fasta', '*.core.full.aln']:
        alignment_files.extend(glob.glob(pattern))
    
    print(f"Found {len(alignment_files)} alignment files")
    
    if not alignment_files:
        print("WARNING: No alignment files found")
        # Create empty outputs
        with open("integrated_core_snps.fa", 'w') as f:
            pass
        with open("sample_cluster_mapping.tsv", 'w') as f:
            f.write("sample_id\\tcluster_id\\ttotal_snps\\n")
        with open("integration_report.txt", 'w') as f:
            f.write("No alignment files found for integration\\n")
        return
    
    # Read cluster information
    try:
        clusters_df = pd.read_csv("${cluster_info}", sep='\\t')
        sample_to_cluster = dict(zip(clusters_df['sample_id'], clusters_df['cluster_id']))
        print(f"Loaded cluster assignments for {len(sample_to_cluster)} samples")
    except Exception as e:
        print(f"Warning: Could not read cluster info: {e}")
        sample_to_cluster = {}
    
    # Process each alignment file
    all_sequences = {}
    all_positions = set()
    cluster_positions = defaultdict(set)
    sample_cluster_map = {}
    
    for alignment_file in alignment_files:
        print(f"Processing alignment: {alignment_file}")
        
        # Try to determine cluster ID from filename
        cluster_id = None
        for possible_cluster in set(sample_to_cluster.values()):
            if str(possible_cluster) in alignment_file:
                cluster_id = possible_cluster
                break
        
        if cluster_id is None:
            # Extract from filename pattern
            base_name = os.path.basename(alignment_file)
            if base_name.startswith('cluster_'):
                cluster_id = base_name.split('_')[1].split('.')[0]
            else:
                cluster_id = base_name.split('.')[0]
        
        print(f"Assigned cluster ID: {cluster_id}")
        
        try:
            # Read sequences from alignment
            sequences = {}
            for record in SeqIO.parse(alignment_file, "fasta"):
                sequences[record.id] = str(record.seq)
                sample_cluster_map[record.id] = cluster_id
            
            if not sequences:
                print(f"WARNING: No sequences found in {alignment_file}")
                continue
            
            print(f"Read {len(sequences)} sequences from cluster {cluster_id}")
            
            # Get alignment length
            seq_names = list(sequences.keys())
            if seq_names:
                alignment_length = len(sequences[seq_names[0]])
                print(f"Alignment length: {alignment_length} bp")
                
                # Find variable positions (SNPs) in this cluster
                cluster_snps = []
                for pos in range(alignment_length):
                    bases_at_pos = set()
                    for name in seq_names:
                        base = sequences[name][pos]
                        if base not in ['-', 'N', 'n']:
                            bases_at_pos.add(base.upper())
                    
                    # If more than one base type, it's a SNP
                    if len(bases_at_pos) > 1:
                        global_pos = f"{cluster_id}_{pos}"
                        all_positions.add(global_pos)
                        cluster_positions[cluster_id].add(global_pos)
                        cluster_snps.append(pos)
                        
                        # Store sequences for this position
                        for name in seq_names:
                            if name not in all_sequences:
                                all_sequences[name] = {}
                            base = sequences[name][pos]
                            all_sequences[name][global_pos] = base.upper() if base not in ['-', 'N', 'n'] else 'N'
                
                print(f"Found {len(cluster_snps)} variable positions in cluster {cluster_id}")
            
        except Exception as e:
            print(f"ERROR processing {alignment_file}: {e}")
            continue
    
    print(f"Total samples: {len(all_sequences)}")
    print(f"Total variable positions: {len(all_positions)}")
    
    if not all_sequences or not all_positions:
        print("WARNING: No sequences or variable positions found")
        # Create minimal outputs
        with open("integrated_core_snps.fa", 'w') as f:
            f.write(">dummy\\nN\\n")
        
        sample_data = []
        for sample_id, cluster_id in sample_cluster_map.items():
            sample_data.append({
                'sample_id': sample_id,
                'cluster_id': cluster_id,
                'total_snps': 0
            })
        
        if sample_data:
            sample_df = pd.DataFrame(sample_data)
            sample_df.to_csv("sample_cluster_mapping.tsv", sep='\\t', index=False)
        else:
            with open("sample_cluster_mapping.tsv", 'w') as f:
                f.write("sample_id\\tcluster_id\\ttotal_snps\\n")
        
        with open("integration_report.txt", 'w') as f:
            f.write("No valid sequences or variable positions found\\n")
        return
    
    # Sort positions for consistent output
    sorted_positions = sorted(all_positions)
    
    # Create integrated alignment
    print("Creating integrated core SNPs alignment...")
    
    with open("integrated_core_snps.fa", 'w') as f:
        for sample_id in sorted(all_sequences.keys()):
            # Build sequence from all variable positions
            integrated_seq = []
            for pos in sorted_positions:
                if pos in all_sequences[sample_id]:
                    integrated_seq.append(all_sequences[sample_id][pos])
                else:
                    integrated_seq.append('N')  # Missing data
            
            sequence = ''.join(integrated_seq)
            f.write(f">{sample_id}\\n{sequence}\\n")
    
    # Create sample mapping
    sample_data = []
    for sample_id in all_sequences.keys():
        cluster_id = sample_cluster_map.get(sample_id, 'unknown')
        # Count non-N positions for this sample
        total_snps = sum(1 for pos in sorted_positions 
                        if pos in all_sequences[sample_id] and all_sequences[sample_id][pos] != 'N')
        
        sample_data.append({
            'sample_id': sample_id,
            'cluster_id': cluster_id,
            'total_snps': total_snps
        })
    
    sample_df = pd.DataFrame(sample_data)
    sample_df.to_csv("sample_cluster_mapping.tsv", sep='\\t', index=False)
    
    # Create integration report
    with open("integration_report.txt", 'w') as f:
        f.write("CORE SNPs INTEGRATION REPORT\\n")
        f.write("=" * 50 + "\\n\\n")
        
        f.write(f"Total samples integrated: {len(all_sequences)}\\n")
        f.write(f"Total variable positions: {len(all_positions)}\\n")
        f.write(f"Clusters processed: {len(set(sample_cluster_map.values()))}\\n\\n")
        
        f.write("Per-cluster statistics:\\n")
        f.write("-" * 30 + "\\n")
        
        cluster_stats = defaultdict(int)
        for sample_id, cluster_id in sample_cluster_map.items():
            cluster_stats[cluster_id] += 1
        
        for cluster_id in sorted(cluster_stats.keys()):
            sample_count = cluster_stats[cluster_id]
            snp_count = len(cluster_positions[cluster_id])
            f.write(f"Cluster {cluster_id}: {sample_count} samples, {snp_count} variable positions\\n")
        
        f.write("\\nIntegration method: Concatenation of cluster-specific variable positions\\n")
        f.write("Missing data handling: Filled with 'N' for positions not present in sample's cluster\\n")
        f.write("Output format: FASTA alignment with all samples and integrated SNP positions\\n")
    
    print("Core SNPs integration completed successfully!")

# Run integration
integrate_core_snps()
EOF

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
        biopython: \$(python -c "import Bio; print(Bio.__version__)")
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
        numpy: \$(python -c "import numpy; print(numpy.__version__)")
    END_VERSIONS
    """
}