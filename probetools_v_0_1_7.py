#!/usr/bin/env python3


import sys
import subprocess
import os


def main():
    version = '0.1.7'
    # Parse command line arguments
    module, args = parse_args(sys.argv, version)
    # Set path to output directory and name to append to output files
    out_path, name = os.path.split(args['-o'])
    out_path = '.' if out_path == '' else out_path
    # Make sure output directory exists
    if os.path.exists(out_path) == False and os.path.isdir(out_path) == False:
        print(f'\nERROR: Output path {out_path} does not exist.\n')
        exit(1)
    # Run top-level function for selected module
    if module == 'clusterkmers':
        print(f'\nProbeTools ClusterKmers v{version}')
        print('https://github.com/KevinKuchinski/ProbeTools\n')
        cluster_kmers(out_path, name, args['-t'], args['-k'], args['-i'], args['-s'], args['-p'], args['-n'], args['-T'])
    elif module == 'capture':
        print(f'\nProbeTools Capture v{version}')
        print('https://github.com/KevinKuchinski/ProbeTools\n')
        capture(out_path, name, args['-t'], args['-p'], args['-i'], args['-l'], args['-T'])
    elif module == 'getlowcov':
        print(f'\nProbeTools GetLowCov v{version}')
        print('https://github.com/KevinKuchinski/ProbeTools\n')
        get_low_cov(out_path, name, args['-i'], args['-k'], args['-D'], args['-L'])
    elif module == 'stats':
        print(f'\nProbeTools Stats v{version}')
        print('https://github.com/KevinKuchinski/ProbeTools\n')
        stats(out_path, name, args['-i'])
    elif module == 'makeprobes':
        print(f'\nProbeTools MakeProbes v{version}')
        print('https://github.com/KevinKuchinski/ProbeTools\n')
        make_probes(out_path, name, args['-t'], args['-b'], args['-m'], args['-c'], args['-k'], args['-i'], args['-s'], 
                    args['-D'], args['-L'], args['-i'], args['-l'], args['-T'])
    elif module == 'merge':
        print(f'\nProbeTools Merge v{version}')
        print('https://github.com/KevinKuchinski/ProbeTools\n')
        merge(args['-i'], args['-I'], args['-o'])
    print('\nDone.\n')
    exit(0)


########## Command line interface functions ##########
def parse_args(args, version):
    if len(args) < 2:
        print('\nERROR: A ProbeTools module must be selected.\n')
        print_usage(None, version)
        exit(1)
    else:
        module = args[1]
    arg_values = {}
    for arg_1, arg_2 in zip(args[1:-1], args[2:]):
        if arg_1[0] == '-':
            if arg_2[0] != '-':
                arg_values[arg_1] = arg_2
            else:
                arg_values[arg_1] = ''
    if args[-1][0] == '-':
        arg_values[args[-1]] = ''
    # Set defaults, mins, and maxs depending on selected module
    if module == 'clusterkmers':
        required_args = {'-t', '-o'}
        arg_value_types = {'-t': str, '-o': str, '-k': int, '-s': int, '-d': int, '-i': float, '-p': str, '-n': int, '-T': int}
        min_arg_values = {'-k': 32, '-s': 1, '-d': 0, '-i': 50, '-n': 1, '-T': 1}
        max_arg_values = {'-i': 100}
        default_arg_values = {'-k': 120, '-s': 1, '-d': 0, '-i': 90, '-p': None, '-n': 'MAX', '-T': 0}
    elif module == 'capture':
        required_args = {'-t', '-p', '-o'}
        arg_value_types = {'-t': str, '-p': str, '-o':str, '-i': float, '-l': int, '-T': int}
        min_arg_values = {'-i': 50, '-l': 1, '-T': 1}
        max_arg_values = {'-i': 100}
        default_arg_values = {'-i': 90, '-l': 60, '-T': 1}
    elif module == 'getlowcov':
        required_args = {'-i', '-o'}
        arg_value_types = {'-i': str, '-o': str, '-k': int, '-D': int, '-L': int}
        min_arg_values = {'-k': 32, '-D': 0, '-L': 1}
        max_arg_values = {}
        default_arg_values = {'-k': 120, '-D': 0, '-L': 40}
    elif module == 'stats':
        required_args = {'-i', '-o'}
        arg_value_types = {'-i': str, '-o': str}
        min_arg_values = {}
        max_arg_values = {}
        default_arg_values = {}
    elif module == 'makeprobes':
        required_args = {'-t', '-b', '-o'}
        arg_value_types = {'-t': str, '-b': int, '-o': str, '-m': int, '-c': float, '-k': int, 
                           '-s': int, '-d': int, '-D': int, '-L': int, '-i': float, '-l': int, '-T': int}
        min_arg_values = {'-m': 1, '-c': 0, '-k': 32, '-s': 1, '-d': 0, '-D': 0, '-L': 1, '-i': 50, '-l': 1, '-T': 1}
        max_arg_values = {'-c': 100, '-i': 100}
        default_arg_values = {'-m': 'MAX', '-c': 90, '-k': 120, '-s': 1, '-d': 0, '-D': 0, '-L': 40, '-i': 90, '-l': 60, '-T': 0}
    elif module == 'merge':
        required_args = {'-i', '-I', '-o'}
        arg_value_types = {'-i': str, '-I': str, '-o': str}
        min_arg_values = {}
        max_arg_values = {}
        default_arg_values = {}
    else:
        print('\nERROR: Module not recognized.')
        print_usage(None, version)
        exit(1)
    # Check if all required arguments were provided
    missing_args = set()
    for required_arg in required_args:
        if required_arg not in arg_values.keys() or arg_values[required_arg] == '':
            missing_args = missing_args | {required_arg}
    if missing_args != set():
        print(f'\nERROR: Values must be provided for the argument following arguments: {", ".join(sorted(missing_args))}')
        print_usage(module, version)
        exit(1)
    # Check if unrecognized arguments were provided
    recognized_args = required_args | set(arg_value_types.keys())
    unrecognized_args = set()
    for provided_arg in arg_values.keys():
        if provided_arg not in recognized_args:
            unrecognized_args = unrecognized_args | {provided_arg}
    if unrecognized_args != set():
        print(f'\nERROR: The following arguments are not recognized: {", ".join(sorted(unrecognized_args))}')
        print_usage(module, version)
        exit(1)
    # Check if arguments were provided without values
    empty_args = set()
    for arg, value in arg_values.items():
        if value == '':
            empty_args = empty_args | {arg}
    if empty_args != set():
        print(f'\nERROR: The following arguments were provided without values: {", ".join(sorted(empty_args))}')
        print_usage(module, version)
        exit(1)
    # Check if provided values are of the correct type
    for arg, value in arg_values.items():
        try:
            arg_values[arg] = arg_value_types[arg](value)
        except ValueError:
            print(f'\nERROR: Value for argument {arg} must be of type {str(arg_value_types[arg].__name__)}')
            print_usage(module, version)
            exit(1)
    # Check if provided values are within the correct range
    for arg, value in arg_values.items():
        if arg in min_arg_values.keys() and value < min_arg_values[arg]:
            print(f'\nERROR: Value for argument {arg} must be at least {min_arg_values[arg]}')
            print_usage(module, version)
            exit(1)
        if arg in max_arg_values.keys() and value > max_arg_values[arg]:
            print(f'\nERROR: Value for argument {arg} must not exceed {max_arg_values[arg]}')
            print_usage(module, version)
            exit(1)        
    # Assign default values to unspecified arguments
    for arg, value in default_arg_values.items():
        if arg not in arg_values.keys():
            arg_values[arg] = value
    # Return keyword args and their values
    return module, arg_values


def print_usage(module, version):
    if module == None:
        print(f'\nProbeTools v{version}')
        print('https://github.com/KevinKuchinski/ProbeTools\n')
        print('Available modules:')
        print('makeprobes - probe panel design using a general purpose incremental strategy')
        print('clusterkmers - enumerate and cluster kmers from target sequences')
        print('capture - assess probe panel coverage of target sequences')
        print('getlowcov - extract low coverage sequences from target space')
        print('stats - calculate probe coverage and depth stats')
        print('merge - merge two sets of capture results\n')
    elif module == 'clusterkmers':
        print(f'\nProbeTools ClusterKmers v{version}')
        print('https://github.com/KevinKuchinski/ProbeTools\n')
        print('Usage: probetools clusterkmers -t <target seqs> -o <path to output directory>/<output name> [<optional args>]\n')
        print('Required arguments:')
        print(' -t : path to target seqs in FASTA file')
        print(' -o : path to output directory and name to append to output files')
        print('Optional arguments:')
        print(' -k : length of kmers to enumerate (default=120, min=32)')
        print(' -s : number of bases separating each kmer (default=1, min=1)')
        print(' -d : number of degenerate bases to permit in probes (default=0, min=0)')
        print(' -i : nucleotide sequence identity (%) threshold used for kmer clustering (default=90, min=50, max=100)')
        print(' -p : path to FASTA file containing previously-generated probe sequences to filter from new probes')
        print(' -n : number of probe candidates to return (default=MAX, min=1)')
        print(' -T : number of threads used by VSEARCH for clustering kmers (default=MAX, min=1)\n')
    elif module == 'capture':
        print(f'\nProbeTools Capture v{version}')
        print('https://github.com/KevinKuchinski/ProbeTools\n')
        print('Usage: probetools capture -t <target seqs> -p <probe seqs> -o <path to output directory>/<output name> [<optional args>]\n')
        print('Required arguments:')
        print(' -t : path to target seqs in FASTA file')
        print(' -p : path to probe sequences in FASTA file')
        print(' -o : path to output directory and name to append to output files')
        print('Optional arguments:')
        print(' -i : nucleotide sequence identity (%) threshold used for probe-target alignments (default=90, min=50, max=100)')
        print(' -l : minimum length for probe-target alignments (default=60, min=1)')
        print(' -T : number of threads used by BLASTn for aligning probes to targets (default=1, min=1)\n')
    elif module == 'getlowcov':
        print(f'\nProbeTools GetLowCov v{version}')
        print('https://github.com/KevinKuchinski/ProbeTools\n')
        print('Usage: probetools getlowcov -i <input file> -o <path to output directory>/<output name> [<optional args>]')
        print('Required arguments:')
        print(' -i : path to capture results in PT file')
        print(' -o : path to output directory and name to append to output files')
        print('Optional arguments:')
        print(' -k : minimum sub-sequence length extracted, should be same as kmer length used for making probes (default=120, min=32)')
        print(' -D : minimum probe depth threshold used to define low coverage sub-sequences (default=0, min=0)')
        print(' -L : minimum number of consecutive bases below probe depth threshold to define a low coverage sub-sequence (default=40, min=1)')
    elif module == 'stats':
        print(f'\nProbeTools Stats v{version}')
        print('https://github.com/KevinKuchinski/ProbeTools\n')
        print('Usage: probetools stats -i <input file> -o <path to output directory>/<output name>')
        print('Required arguments:')
        print(' -i : path to capture results in PT file')
        print(' -o : path to output directory and name to append to output files\n')
    elif module == 'makeprobes':
        print(f'\nProbeTools MakeProbes v{version}')
        print('https://github.com/KevinKuchinski/ProbeTools\n')
        print('Usage: probetools incrementalprobes -t <target seqs> -b <batch size> -o <path to output directory>/<output name> [<optional args>]')
        print('Required arguments:')
        print(' -t : path to target sequences in FASTA file')
        print(' -b : number of probes in each batch (min=1)')
        print(' -o : path to output directory and name to append to output files')
        print('Optional arguments:')
        print(' -m : max number of probes to add to panel (default=MAX, min=1)')
        print(' -c : target for 10th percentile of probe coverage (default=90, min=1, max=100)')
        print(' -k : length of probes to generate (default=120, min=32)')
        print(' -s : number of bases separating each kmer (default=1, min=1)')
        print(' -d : number of degenerate bases to permit in probes (default=0, min=0)')
        print(' -i : nucleotide sequence identity (%) threshold used for kmer clustering and probe-target alignments (default=90, min=50, min=100)')
        print(' -l : minimum length for probe-target alignments (default=60, min=1)')
        print(' -D : minimum probe depth threshold used to define low coverage sub-sequences (default=0, min=0)')
        print(' -L : minimum number of consecutive bases below probe depth threshold to define a low coverage sub-sequence (default=40, min=1)')
        print(' -T : number of threads used by VSEARCH and BLASTn for clustering kmers and aligning'
              f' probes to targets (default=MAX for VSEARCH, default=1 for BLASTn, min=1)\n')
    elif module == 'merge':
        print(f'\nProbeTools Merge v{version}')
        print('https://github.com/KevinKuchinski/ProbeTools\n')
        print('Usage: probetools merge -i <input file 1> -I <input file 2> -o <path to output file>')
        print('Required arguments:')
        print(' -i : path to capture results in PT file')
        print(' -I : path to other capture results in PT file')
        print(' -o : path to merged capture results PT output file\n')


########## Helper functions ##########
def check_input(files):
    '''Checks list of input file paths to make sure they exist.'''
    for file in files:
        if os.path.exists(file) == False:
            print(f'\nERROR: Input file {file} does not exist.\n')
            exit(1)


def load_fasta(fasta_path):
    '''Loads contents of a FASTA file into two lists (one for headers, and another for their
    corresponding seqs).'''
    # Check if targets FASTA is valid                                                                                                                                           
    check_input([fasta_path])
    # Open targets FASTA and load headers and seqs into lists                                                                                                                   
    with open(fasta_path, 'r') as input_file:
        entry_counter = 0
        header = ''
        headers, seqs = [], []
        for line in input_file:
            if line[0] == '>':
                entry_counter += 1
                header = line.strip().lstrip('>')
                headers.append(header)
                seqs.append('')
            elif header != '':
                seqs[-1] += line.strip()
    for header, seq in zip(headers, seqs):
        if header == '':
            print('\nWARNING: Seq is present without header in {fasta_path}.')
        if seq == '':
            print('\nWARNING: Header {header} in {fasta_path} has no accompanying seq.')
    return headers, seqs


def append_fasta(fasta_path, new_fasta_path):
    '''Appends contents of new FASTA file to existing FASTA file then deletes new FASTA file.'''
    with open(new_fasta_path, 'r') as input_file, open(fasta_path, 'a') as output_file:
        for line in input_file:
            output_file.write(line)
    # GARBAGE COLLECTION - Delete new FASTA file after contents have been appended
    os.remove(new_fasta_path)


######### Top-level functions for modules ##########
def cluster_kmers(out_path, name, targets_path, k, cluster_id, step, prev_probes_path, num_probes, threads):
    check_input([targets_path])
    kmers_path = os.path.join(out_path, name + '_kmers.fa')
    enum_kmers(targets_path, kmers_path, k, step)
    centroids_path = os.path.join(out_path, name + '_centroids.fa')
    cluster_kmers_with_VSEARCH(kmers_path, centroids_path, cluster_id, threads)
    remove_prev_probes(centroids_path, prev_probes_path)
    probe_names, probe_seqs = rank_probe_candidates(centroids_path)
    probes_path = os.path.join(out_path, name + '_probes.fa')
    potential_probes, probes_writen = write_top_probes(probe_names, probe_seqs, num_probes, probes_path, name)
    return potential_probes, probes_writen


def capture(out_path, name, targets_path, probes_path, min_id, min_length, threads):
    check_input([targets_path, probes_path])
    blast_path = os.path.join(out_path, name + '_blast_results.tsv')
    align_probes_to_targets_with_BLAST(targets_path, probes_path, blast_path, threads)
    capture_data = create_empty_capture_data(targets_path)
    capture_data = add_BLAST_results_to_capture_data(blast_path, capture_data, min_id, min_length)
    capture_path = os.path.join(out_path, name + '_capture.pt')
    write_capture_data(capture_path, capture_data)


def get_low_cov(out_path, name, capture_path, k, min_depth, min_length):
    check_input([capture_path])
    capture_data = load_capture_data(capture_path)
    low_cov_path = os.path.join(out_path, name + '_low_cov_seqs.fa')
    seqs_writen = write_low_cov_seqs(capture_data, low_cov_path, k, min_depth, min_length)


def stats(out_path, name, capture_path):
    '''Top-level function for stats module.'''
    check_input([capture_path])
    capture_data = load_capture_data(capture_path)
    report_path = os.path.join(out_path, name + '_summary_stats_report.tsv')
    write_summary_report(capture_data, report_path, name)
    report_path = os.path.join(out_path, name + '_long_stats_report.tsv')
    write_long_report(capture_data, report_path, name)


def make_probes(out_path, name, targets_path, batch_size, max_probes, cov_target, k, cluster_id, step,
                min_depth, min_low_cov_length, min_id, min_capture_length, threads):
    check_input([targets_path])
    # Set variable for counting rounds of incremental probe design
    round_counter = 0
    # Initialize variables for panel size (# probes) and 10th percentile of panel coverage
    panel_size, panel_cov = 0, 0
    # If panel size is not specified, set max_panel_size as 1 probe larger than current panel size so loop doesn't break
    max_panel_size = panel_size + batch_size if max_probes == 'MAX' else max_probes
    # Create empty FASTA file for final probe panel
    final_probes_path = os.path.join(out_path, name + '_probes.fa')
    with open(final_probes_path, 'w') as output_file:
        pass
    # Create empty capture dict from target space
    capture_data = create_empty_capture_data(targets_path)
    print()
    # Enter main incremental design loop
    while panel_size < max_panel_size and panel_cov < cov_target:
        round_counter += 1
        print('*' * 20, f'ROUND {round_counter}', '*' * 20 + '\n')
        # Write low cov seqs to FASTA as target space for this round
        low_cov_path = os.path.join(out_path, name + '_low_cov_seqs.fa')
        seqs_writen = write_low_cov_seqs(capture_data, low_cov_path, k, min_depth, min_low_cov_length)
        print()
        # Break design loop if no low seqs were writen
        if seqs_writen == 0:
            break
        # Make probes from target space
        num_probes = min(batch_size, max_panel_size - panel_size)
        round_name = name + f'_round_{round_counter}'
        potential_probes, probes_writen = cluster_kmers(out_path, round_name, low_cov_path, k, cluster_id, step, final_probes_path, num_probes, threads)
        print()
        # GARBAGE COLLECTION - Delete low cov seqs
        os.remove(low_cov_path)
        # Update panel size and max panel size
        panel_size += probes_writen
        max_panel_size = panel_size + batch_size if max_probes == 'MAX' else max_probes
        # Break design loop if no potential probes or no probes writen
        if potential_probes == 0 or probes_writen == 0:
            os.remove(low_cov_path)
            break
        # Capture probes against original targets
        blast_path = os.path.join(out_path, name + '_blast_results.tsv')
        new_probes_path = os.path.join(out_path, round_name + '_probes.fa')
        align_probes_to_targets_with_BLAST(targets_path, new_probes_path, blast_path, threads)
        capture_data = add_BLAST_results_to_capture_data(blast_path, capture_data, min_id, min_capture_length)
        print()
        # Append new probes to final panel
        append_fasta(final_probes_path, new_probes_path)
        # Check if panel cov exceeds cov target
        panel_cov = calc_cov_percentiles(capture_data, percentiles=(0.1,))[0]
        print(f'10th percentile of target coverage: {panel_cov}%\n')
    print('Incremental probe design finished.')
    if panel_size >= max_panel_size:
        print(' Maximum panel size reached.')
    if panel_cov >= cov_target:
        print(' Coverage target for 10th percentile of targets reached.')
    if seqs_writen == 0 or potential_probes == 0 or probes_writen == 0:
        print(' No additional probes could be designed.')
    capture_path = os.path.join(out_path, name + '_capture.pt')
    print(f'\nWriting capture results to {capture_path}...')
    write_capture_data(capture_path, capture_data)
    report_path = os.path.join(out_path, name + '_summary_stats_report.tsv')
    write_summary_report(capture_data, report_path, name)
    report_path = os.path.join(out_path, name + '_long_stats_report.tsv')
    write_long_report(capture_data, report_path, name)


def merge(capture_path, other_capture_path, merged_capture_path):
    capture_data = load_capture_data(capture_path)
    other_capture_data = load_capture_data(other_capture_path)
    merged_capture_data = merge_capture_results(capture_data, other_capture_data)
    write_capture_data(merged_capture_path, merged_capture_data)


########## clusterkmers functions ##########
def enum_kmers(targets_path, kmers_path, k, step):
    print(f'Enumerating all {k}-mers in {targets_path}...')
    # Load contents of targets FASTA
    headers, seqs = load_fasta(targets_path)
    print(f' Loaded {"{:,}".format(len(seqs))} target seqs in {targets_path}...')
    # Create FASTA file for kmers output and enumerate kmers from target seqs
    with open(kmers_path, 'w') as output_file:
        target_counter = 0
        kmer_counter = 0
        for header, seq in zip(headers, seqs):
            if seq == '':
                print(f' WARNING --- Target {header} has no sequence!')
            elif len(seq) < k:
                print(f' WARNING --- Target {header} is shorter than the desired probe length!')
            elif len(seq) >= k and seq != '':
                target_counter += 1
                for i in range(0, len(seq) - k + 1, step):
                    kmer_counter += 1
                    output_file.write(f'>kmer_{kmer_counter}\n')
                    output_file.write(seq[i:i+k] + '\n')
    print(f' Enumerated {"{:,}".format(kmer_counter)} k-mers from {targets_path}.')


def cluster_kmers_with_VSEARCH(kmers_path, centroids_path, cluster_id, threads):
    print(f' Clustering k-mers at {cluster_id}% identity...')
    # Create and run terminal command for VSEARCH
    terminal_command = (f'vsearch --cluster_fast {kmers_path} --id {cluster_id / 100} --centroids {centroids_path}'
                        f' --fasta_width 0 --sizeout --qmask none --threads {threads}')
    finished_process = subprocess.call(terminal_command, stdout=subprocess.DEVNULL, 
                                       stderr=subprocess.DEVNULL, shell=True)
    if finished_process != 0:
        print(f'\nERROR: vsearch cluster_fast terminated with errors while clustering k-mers'
              f' (Error code: {finished_process}).\n')
        exit(1)
    # GARBAGE COLLECTION - Delete kmers FASTA file
    os.remove(kmers_path)
    print(f' K-mer clustering finished.')


def remove_prev_probes(centroids_path, prev_probes_path):
    if prev_probes_path != None:
        print(f'Removing any previously designed probes...')
        prev_probe_headers, prev_probe_seqs = load_fasta(prev_probes_path)
        print(f' Loaded {"{:,}".format(len(prev_probe_seqs))} previously designed probes in {prev_probes_path}')
        centroid_headers, centroid_seqs = load_fasta(centroids_path)
        with open(centroids_path, 'w') as output_file:
            for header, seq in zip(centroid_headers, centroid_seqs):
                if seq not in prev_probe_seqs:
                    output_file.write('>' + header + '\n')
                    output_file.write(seq + '\n')


def rank_probe_candidates(centroids_path):
    print(f'Ranking probe candidates...')
    # Load centroids FASTA file
    headers, seqs = load_fasta(centroids_path)
    print(f' Loaded {"{:,}".format(len(seqs))} probe candidates.')
    # Create dicts for looking up each candidate probe seq's cluster size and probe name
    cluster_sizes = {seq: int(header.split('size=')[1].split(';')[0]) for header, seq in zip(headers, seqs)}
    probe_names = {seq: header.split(';')[0] for header, seq in zip(headers, seqs)}
    # Rank candidate probe seqs based on cluster size
    probe_seqs = sorted(seqs, key=lambda seq: cluster_sizes[seq], reverse=True)
    # GARBAGE COLLCTION - Delete centroids FASTA file
    os.remove(centroids_path)
    return probe_names, probe_seqs


def write_top_probes(probe_names, probe_seqs, num_probes, probes_path, name):
    print(f'Writing top probe candidates to {probes_path}...')
    potential_probes = len(probe_seqs)
    # Check if there are any probes to write
    if potential_probes == 0:
        print(' No probe candidates to write to FASTA.')
        return 0, 0
    # Determine how many probes to write to FASTA file
    if num_probes == 'MAX':
        max_probes = len(probe_seqs)
    else:
        max_probes = num_probes
    # Write probes to FASTA file
    probe_counter = 0
    with open(probes_path, 'w') as output_file:
        while probe_counter < len(probe_seqs) and probe_counter < max_probes:
            probe_seq = probe_seqs[probe_counter]
            probe_counter += 1
            output_file.write(f'>{name}_probe_{probe_counter}\n')
            output_file.write(probe_seq + '\n')
    print(f' Wrote {"{:,}".format(probe_counter)} probes to {probes_path}.')
    return potential_probes, probe_counter


########## capture functions ##########
def align_probes_to_targets_with_BLAST(targets_path, probes_path, blast_path, threads):
    print(f'Aligning probes in {probes_path} to targets in {targets_path}...')
    # Check that BLAST db files exist for target seqs
    if any([os.path.exists(targets_path + '.' + suffix) == False for suffix in ['nhr', 'nin' , 'nsq']]):
        print(' WARNING: blastn db files do not exist for target seqs. Creating blastn db...')
        # Create terminal command for making BLASTn db
        terminal_command = f'makeblastdb -in {targets_path} -dbtype nucl'
        # Run terminal command and redirect stdout and stderr to DEVNULL
        finished_process = subprocess.run(terminal_command, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, shell=True)
        # Exit with errors if BLASTn terminates with errors
        if finished_process.returncode != 0:
            print(f'\nERROR: makeblastdb terminated with errors while making db for target seqs (Error code: {finished_process.returncode}).\n')
            exit(1)
    # Count target seqs
    headers, seqs = load_fasta(targets_path)
    num_targets = len(seqs)
    # Create and run terminal command for BLASTn
    threads = 1 if threads == 0 else threads
    terminal_command = (f"blastn -db {targets_path} -query {probes_path} -max_target_seqs {num_targets}"
                        f" -num_threads {threads} -outfmt 6 > {blast_path}")
    finished_process = subprocess.run(terminal_command, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, shell=True)
    if finished_process.returncode != 0:
        print(f'\nERROR: blastn terminated with errors while aligning probes to target seqs (Error code: {finished_process.returncode}).\n')
        exit(1)


def create_empty_capture_data(targets_path):
    print(f'Extracting target names and seqs from {targets_path}...')
    with open(targets_path, 'r') as input_file:
        headers, seqs = [], []
        for line in input_file:
            if line[0] == '>':
                header = line.strip().lstrip('>')
                headers.append(header)
                seqs.append('')
            else:
                seqs[-1] += line.strip()
    for header in headers:
        if headers.count(header) > 1:
            print(f'\nERROR: Header {header} appears more than once in {targets_path}.\n')
            exit(1)
        if ' ' in header:
            print(f'\nERROR: Header {header} contains spaces.\n')
            exit(1)
    capture_data = {header: (seq, [0] * len(seq)) for header, seq in zip(headers, seqs)}
    print(f' Total targets loaded: {"{:,}".format(len(capture_data))}')
    return capture_data


def load_capture_data(capture_path):
    print(f'Loading capture data from {capture_path}...')
    with open(capture_path, 'r') as input_file:
        headers, seqs, depths = [], [], []
        for line in input_file:
            if line[0] == '>':
                header = line.strip().lstrip('>')
                headers.append(header)
                seqs.append('')
                depths.append('')
            elif line[0] == '$':
                seqs[-1] += line.strip()
            elif line[0] == '#':
                depths[-1] += line.strip()
    seqs = [seq.lstrip('$') for seq in seqs]
    depths = [[int(d) for d in depth.lstrip('#').split(',')] for depth in depths]
    if len(set(len(c) for c in [headers, seqs, depths])) != 1:
        print(f'\nERROR: The number of headers, seqs, and probe depth lists do not match in {capture_path}.')
        exit(1)
    for header in headers:
        if headers.count(header) > 1:
            print(f'\nERROR: Header {header} appears more than once in {capture_path}.\n')
            exit(1)
    capture_data = {header: (seq, depth) for header, seq, depth in zip(headers, seqs, depths)}
    for header in capture_data.keys():
        if len(capture_data[header][0]) != len(capture_data[header][1]):
            print(f'\nERROR: Seq length and probe depth list length do not match for entry {header}.')
            exit(1)
    print(f' Total targets loaded: {"{:,}".format(len(capture_data))}')
    return capture_data


def merge_capture_results(capture_data, other_capture_data):
    headers, other_headers = set(capture_data.keys()), set(other_capture_data.keys())
    shared_headers = headers & other_headers
    merged_capture_data = {}
    for header in shared_headers:
        if capture_data[header][0] != other_capture_data[header][0]:
            print(f'\nERROR: Header {header} appears in both capture results but they do not have the same seq.')
            exit(1)
        if len(capture_data[header][1]) != len(other_capture_data[header][1]):
            print(f'\nERROR: Header {header} appears in both capture results but their probe depth lists are different lengths.')
            exit(1)
        merged_capture_data[header] = (capture_data[header][0], [a + b for a, b in zip(capture_data[header][1], other_capture_data[header][1])])
    for header in headers - shared_headers:
        merged_capture_data[header] = capture_data[header]
    for header in other_headers - shared_headers:
        merged_capture_data[header] = other_capture_data[header]
    return merged_capture_data


def add_BLAST_results_to_capture_data(blast_path, capture_data, min_id, min_length):
    print(f'Extracting probe coverage and probe depth from alignment of probes against targets...')
    with open(blast_path, 'r') as input_file:
        for line in input_file:
            line = line.split('\t')
            if float(line[2]) >= min_id and float(line[3]) >= min_length:
                target = line[1].strip()
                aln_start = int(line[8].strip()) - 1
                aln_end = int(line[9].strip()) - 1
                aln_start, aln_end = sorted([aln_start, aln_end])
                for position in range(aln_start, aln_end + 1):
                    capture_data[target][1][position] += 1
    # GARBAGE COLLECTION - Delete BLAST results
    os.remove(blast_path)
    return capture_data


def write_capture_data(capture_path, capture_data):
    print(f'Writing capture results to {capture_path}...')
    with open(capture_path, 'w') as output_file:
        for header, (seq, depth) in capture_data.items():
            output_file.write('>' + header + '\n')
            output_file.write('$' + seq + '\n')
            output_file.write('#' + ','.join([str(d) for d in depth]) + '\n')
    print(f' Wrote capture results for {"{:,}".format(len(capture_data))} targets.')


########## getlowcov functions ##########
def write_low_cov_seqs(capture_data, low_cov_path, k, min_depth, min_length):
    print(f'Extracting low coverage areas...')
    # Write low cov seqs to FASTA file
    with open(low_cov_path, 'w') as output_file: 
        total_low_cov_counter = 0
        for header, (seq, depth) in capture_data.items():
            if len(seq) < k:
                print(f' WARNING: Seq {header} is shorter than the window length so low coverage seqs cannot be extracted.')
            else:
                low_cov_counter = 0
                low_cov_start = 0
                while low_cov_start < len(seq):
                    while low_cov_start < len(seq) and depth[low_cov_start] > min_depth:
                        low_cov_start += 1
                    low_cov_end = low_cov_start
                    while low_cov_end < len(seq) and depth[low_cov_end] <= min_depth:
                        low_cov_end += 1
                    low_cov_seq = seq[low_cov_start:low_cov_end]
                    if len(low_cov_seq) > min_length and len(low_cov_seq) >= k:
                        total_low_cov_counter += 1
                        low_cov_counter += 1
                        output_file.write(f'>{header}_position_{low_cov_start + 1}_to_{low_cov_end}\n')
                        output_file.write(low_cov_seq + '\n')
                    elif len(low_cov_seq) > min_length and len(low_cov_seq) < k:
                        window_start = low_cov_start - int((k - len(low_cov_seq)) / 2)
                        window_end = window_start + k
                        if window_start < 0:
                            window_end = window_end - window_start
                            window_start = 0
                        if window_end > len(seq):
                            window_start = window_start - window_end + len(seq) 
                            window_end = len(seq)
                        extended_low_cov_seq = seq[window_start:window_end]
                        total_low_cov_counter += 1
                        low_cov_counter += 1
                        output_file.write(f'>{header}_position_{window_start + 1}_to_{window_end}\n')
                        output_file.write(extended_low_cov_seq + '\n')
                    low_cov_start = low_cov_end
    print(f' Wrote {total_low_cov_counter} low coverage seqs from {"{:,}".format(len(capture_data))} targets.')
    seqs_writen = total_low_cov_counter
    return seqs_writen


########## stats functions ##########
def calc_cov_percentiles(capture_data, percentiles=(0, 0.05, 0.1, 0.25, 0.5, 0.75, 1)):
    """Takes a capture dict and tuple of perentiles, and returns a list of those percentile
    values."""
    cov_values = []
    for header,(seq, depth) in capture_data.items():
        cov_values.append((len(depth) - depth.count(0)) * 100 / len(depth))
    cov_values = sorted(cov_values)
    percentile_values = []
    for percentile in percentiles:
        values_under_percentile = (len(cov_values) - 1) * percentile
        if values_under_percentile == int(values_under_percentile):
            percentile_values.append(cov_values[int(values_under_percentile)])
        else:
            percentile_value = cov_values[int(values_under_percentile)]
            percentile_value += cov_values[int(values_under_percentile) + 1]
            percentile_value = percentile_value / 2
            percentile_values.append(percentile_value)
    return [round(percentile_value,2) for percentile_value in percentile_values] 


def calc_total_probe_depth(capture_data):
    """Takes a capture dict and returns a tuple containing the percentage of nucleotide positions
    in the target space covered by 0, 1, 2, 3, 4, and 5+ probes."""
    total = 0
    total_0 = 0
    total_1 = 0
    total_2 = 0
    total_3 = 0
    total_4 = 0
    total_5 = 0
    for header,(seq, depth) in capture_data.items():
        total += len(depth)
        total_0 += depth.count(0)
        total_1 += depth.count(1)
        total_2 += depth.count(2)
        total_3 += depth.count(3)
        total_4 += depth.count(4)
        total_5 += len([d for d in depth if d >= 5])
    total_0 = round(total_0 * 100 / total, 2)
    total_1 = round(total_1 * 100 / total, 2)
    total_2 = round(total_2 * 100 / total, 2)
    total_3 = round(total_3 * 100 / total, 2)
    total_4 = round(total_4 * 100 / total, 2)
    total_5 = round(total_5 * 100 / total, 2)
    return (total_0, total_1, total_2, total_3, total_4, total_5)


def write_summary_report(capture_data, report_path, name):
    """Takes a capture dict and generates the summary report."""
    print(f'Writing probe coverage and probe depth summary stats to {report_path}...')
    with open(report_path, 'w') as output_file:
        header = ['name', 'total_targets']
        header += ['cov_' + s for s in ['min', '5%tile', '10%tile', 'Q1', 'med', 'Q3', 'max']]
        header += ['depth_' + s for s in ['0', '1', '2', '3', '4', '5+']]
        output_file.write('\t'.join(header) + '\n')
        line = [name, str(len(capture_data))]
        line += [str(f) for f in calc_cov_percentiles(capture_data)]
        line += [str(f) for f in calc_total_probe_depth(capture_data)]
        output_file.write('\t'.join(line) + '\n')
        

def write_long_report(capture_data, report_path, name):
    """Takes a capture dict and generates the long-form report."""
    print(f'Writing probe coverage and probe depth stats for each target to {report_path}...')
    with open(report_path, 'w') as output_file:
        header = ['name', 'target', 'length']
        header += ['bases_' + str(d) + 'X' for d in ['0', '1', '2', '3', '4', '5+']]
        header += ['%_cov_' + str(d) + 'X' for d in ['1', '2', '3', '4', '5']]
        output_file.write('\t'.join(header) + '\n')
        for header,(seq, depth) in capture_data.items():
            line = [name, header, str(len(depth))]
            line += [str(depth.count(d)) for d in [0, 1, 2, 3, 4]]
            line += [str(len([d for d in depth if d >= 5]))]
            line += [str(round(len([d for d in depth if d >= dd]) * 100 / len(depth), 1)) for dd in [1, 2, 3, 4, 5]]
            output_file.write('\t'.join(line) + '\n')


######### call main function ##########
if __name__ == '__main__':
    main()
