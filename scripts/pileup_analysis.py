import sys
import csv

def parse_pileup_line(line):
    fields = line.strip().split()
    if len(fields) < 5:
        return None  # Skip malformed lines

    pos = fields[1]
    ref = fields[2]
    try:
        depth = int(fields[3])
    except ValueError:
        depth = 0
    read_bases = fields[4]

    # Initialize nucleotide counts
    base_counts = {'A': 0, 'T': 0, 'G': 0, 'C': 0, '*': 0}

    # Initialize insertion/deletion counts
    ins_types = {}  # inserted sequence -> count
    del_types = {}  # deleted sequence -> count

    i = 0
    while i < len(read_bases):
        char = read_bases[i]

        # Skip start-of-read marker
        if char == '^':
            i += 2
            continue

        # Skip end-of-read marker
        if char == '$':
            i += 1
            continue

        # Handle insertions and deletions
        if char in ('+', '-'):
            sign = char
            i += 1
            num_str = ""
            while i < len(read_bases) and read_bases[i].isdigit():
                num_str += read_bases[i]
                i += 1
            if num_str == "":
                continue
            try:
                indel_length = int(num_str)
            except ValueError:
                continue

            indel_seq = read_bases[i:i+indel_length]
            i += indel_length

            if sign == '+':
                ins_types[indel_seq] = ins_types.get(indel_seq, 0) + 1
            else:
                del_types[indel_seq] = del_types.get(indel_seq, 0) + 1
            continue

        # Process nucleotide calls
        if char in ('.', ','):
            base = ref.upper()
        else:
            base = char.upper()

        if base in base_counts:
            base_counts[base] += 1

        i += 1

    # Determine the best nucleotide call
    best_call = max(base_counts, key=base_counts.get)
    best_count = base_counts[best_call]

    # Sort top 3 insertions and deletions
    top_insertions = sorted(ins_types.items(), key=lambda x: x[1], reverse=True)[:3]
    top_deletions = sorted(del_types.items(), key=lambda x: x[1], reverse=True)[:3]

    # **NEW FIX:** Only count the top 3 insertions/deletions
    total_ins = sum(count for _, count in top_insertions)
    total_del = sum(count for _, count in top_deletions)

    # **Corrected heteroplasmy formula:**
    # total = A + T + G + C + ins + del (excluding unmapped '*' reads)
    # heteroplasmy = (total - ref_count) / total
    total_bases = (base_counts['A'] + base_counts['T'] +
                   base_counts['G'] + base_counts['C'] +
                   total_ins + total_del)
    ref_count = base_counts.get(ref.upper(), 0)
    heteroplasmy = 0.0
    if total_bases > 0:
        heteroplasmy = (total_bases - ref_count) / total_bases

    # Prepare output dictionary
    result = {
        "Pos": pos,
        "RefCall": ref,
        "Depth": depth,
        "A": base_counts['A'],
        "T": base_counts['T'],
        "G": base_counts['G'],
        "C": base_counts['C'],
        "*": base_counts['*'],
        "ins": total_ins,
        "del": total_del,
        "BestCall": best_call,
        "Heteroplasmy": heteroplasmy
    }

    # Add top insertions and deletions
    for i in range(3):
        if i < len(top_insertions):
            result[f"ins_type{i+1}"] = top_insertions[i][0]
            result[f"ins_count{i+1}"] = top_insertions[i][1]
        else:
            result[f"ins_type{i+1}"] = ""
            result[f"ins_count{i+1}"] = 0

        if i < len(top_deletions):
            result[f"del_type{i+1}"] = top_deletions[i][0]
            result[f"del_count{i+1}"] = top_deletions[i][1]
        else:
            result[f"del_type{i+1}"] = ""
            result[f"del_count{i+1}"] = 0

    return result

def process_pileup_file(input_file, output_file):
    header = ["Pos", "RefCall", "Depth", "A", "T", "G", "C", "*",
              "ins", "del", "BestCall", "Heteroplasmy",
              "ins_type1", "ins_count1", "ins_type2", "ins_count2", "ins_type3", "ins_count3",
              "del_type1", "del_count1", "del_type2", "del_count2", "del_type3", "del_count3"]

    with open(input_file, 'r') as infile, open(output_file, 'w', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=header)
        writer.writeheader()

        for line in infile:
            if not line.strip():
                continue
            record = parse_pileup_line(line)
            if record is not None:
                writer.writerow(record)

def main():
    if len(sys.argv) != 3:
        sys.exit("Usage: python parse_pileup.py input.pileup output.csv")
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    process_pileup_file(input_file, output_file)
    print(f"CSV file '{output_file}' has been created from '{input_file}'.")

if __name__ == '__main__':
    main()

