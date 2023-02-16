

import random
import argparse
import os
import sys
import yaml
import attrmap as ap


def sectorise_genome(input_string, min_length=1, max_length=3):
    length = len(input_string)
    random_splits = []
    start = 0
    while start + min_length < length:
        split_length = random.randint(min_length, min(max_length, length - start))
        split = input_string[start:start + split_length]
        random_splits.append(split)
        start += split_length
    return random_splits


def add_snps(input_string, substitution_percent=10):
    ACGT = ["A", "C", "G", "T"]
    length = len(input_string)
    output_string = list(input_string)
    for i in range(length):
        if random.randint(1, 100) <= substitution_percent:
            output_string[i] = random.choice(ACGT)
    return "".join(output_string)


def add_indels(input_string, min_insert_size=1, max_insert_size=10, insertion_frequency=10, DNA_bases=["A", "C", "G", "T"]):
    length = len(input_string)
    output_string = list(input_string)
    for i in range(length):
        if random.randint(1, 100) <= insertion_frequency:
            insert_length = random.randint(min_insert_size, max_insert_size)
            if random.randint(0,1)==0:
                output_string = output_string[:i] + output_string[i+insert_length:]
                length -= insert_length
            else:
                insert = "".join([random.choice(DNA_bases) for j in range(insert_length)])
                output_string[i:i] = list(insert)
                length += insert_length
    return "".join(output_string)


def invert_sequence(input_string):
    complement_dict = {
        "A": "T",
        "T": "A",
        "C": "G",
        "G": "C"
    }
    reverse_complement = [complement_dict[base] for base in input_string[::-1]]
    return "".join(reverse_complement)


def parse_fasta(file):
    with open(file, 'r') as f:
        sequences = {}
        header = ''
        sequence = ''
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if header:
                    sequences[header] = sequence
                header = line[1:]
                sequence = ''
            else:
                sequence += line
        sequences[header] = sequence
    return sequences


def insert_sector(lst, value):
    index = random.randint(0, len(lst))
    lst.insert(index, value)
    return lst


def read_random_sequence(folder, min_length, max_length):
    files = [f for f in os.listdir(folder) if os.path.isfile(os.path.join(folder, f))]
    if not files:
        return None
    file = random.choice(files)
    seq = str()
    with open(os.path.join(folder, file), 'r') as f:
        for line in f:
            if line.startswith('>'):
                continue
            seq = seq + line.strip()
    split_length = random.randint(min_length, min(max_length, len(seq)))
    start = random.randint(0, len(seq) - split_length)
    return seq[start:start + split_length]


def mutate_seq(seq, config):
    mutated = []
    sectors = sectorise_genome(
        seq,
        min_length=config.sectorise.min_len,
        max_length=config.sectorise.max_len
    )
    for sector in sectors:
        if random.randint(1, 100) <= config.model.delete:
            continue
        if random.randint(1, 100) <= config.model.snpindel:
            sector = add_snps(
                sector,
                substitution_percent=random.randint(config.ani.min_snp_perc, config.ani.max_snp_perc)
            )
            sector = add_indels(
                sector,
                min_insert_size=config.ani.min_indel_len,
                max_insert_size=config.ani.max_indel_len,
                insertion_frequency=random.randint(config.ani.min_indel_freq, config.ani.max_indel_freq)
            )
        if random.randint(1, 100) <= config.model.insert:
            mutated.append(
                read_random_sequence(config.donors, config.sectorise.min_len, config.sectorise.max_len)
            )
        if random.randint(1, 100) <= config.model.invert:
            sector = invert_sequence(sector)
        if random.randint(1, 100) <= config.model.shuffle:
            mutated = insert_sector(mutated, sector)
        else:
            mutated.append(sector)
    return ''.join(mutated)


def wrap_dna_sequence(sequence):
    wrapped_sequence = ''
    for i in range(0, len(sequence), 80):
        wrapped_sequence += sequence[i:i+80] + '\n'
    return wrapped_sequence


def circularise_seq(s):
    if len(s) <= 600:
        return s
    else:
        first_300 = s[:300]
        last_300 = s[-300:]
        return s + 'N'*300 + first_300 + last_300


def main():
    parser = argparse.ArgumentParser(description='Read in a YAML config file and a text file.')
    parser.add_argument('-c', '--config_file', nargs='?', default='config.yaml', help='The YAML config file (default: config.yaml)')
    parser.add_argument('-f', '--fasta_file', help='Fasta file')
    parser.add_argument('-d', '--genomes_dir', help='directory of genomes')
    parser.add_argument('--circ', action='store_true', help='Circular genome support for simulating reads')
    args = parser.parse_args()

    if args.config_file == 'config.yaml':
        script_dir = os.path.dirname(os.path.realpath(__file__))
        config_file_path = os.path.join(script_dir, 'config.yaml')
        args.config_file = open(config_file_path, 'r')

    config = yaml.safe_load(args.config_file)
    config = ap.AttrMap(config)
    config.donors = args.genomes_dir

    seqs = parse_fasta(args.fasta_file)

    for id, seq in seqs.items():
        for i in range(0,random.randint(config.mutate.min, config.mutate.max)):
            mutated_seq = mutate_seq(seq, config)
            if args.circ:
                mutated_seq = circularise_seq(mutated_seq)
            mutated_seq = wrap_dna_sequence(mutated_seq)
            sys.stdout.write(f'>{id}.{i}\n{mutated_seq}')
            seq = mutated_seq.replace('\n','')


if __name__ == '__main__':
    main()

