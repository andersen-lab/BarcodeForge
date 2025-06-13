import pandas as pd
from Bio import SeqIO
from collections import OrderedDict
import re
from rich.console import Console
from .utils import STYLES  # Assuming STYLES is in utils.py

console = Console()


def _load_sample_mutations(path):
    df = pd.read_csv(path, sep="\t", names=["sample", "mutations"])
    return df


MUT_PATTERN = re.compile(r"([^:]+):([A-Za-z0-9,]+)")


def _extract_mutations(sample):
    muts_str = sample.get("mutations", "")
    if pd.isnull(muts_str) or not muts_str:
        return OrderedDict()
    # parse and build an OrderedDict in reversed order
    return OrderedDict(
        (k.strip(), v.split(",")) for k, v in reversed(MUT_PATTERN.findall(muts_str))
    )


def _reverse_mutations_to_root(muts):
    """Reverse the mutations to the root node."""
    root_muts = {}
    if muts == {}:
        root_muts[0] = {
            "base": "",
            "mut": "",
        }
        return root_muts
    for value in muts.values():
        for i in value:
            nuc_loc = int(i[1:-1])
            # Here the tip mutations are reversed to the root node.
            # The base is the nucleotide of the tip node.
            if nuc_loc in root_muts.keys():
                root_muts[nuc_loc]["mut"] = i[0]
            else:
                root_muts[nuc_loc] = {
                    "base": i[-1],
                    "mut": i[0],
                }
    return root_muts


def _construct_root_sequence(root_muts, seq):
    """Generate the root sequence."""
    for key, value in root_muts.items():
        seq.seq = seq.seq[: int(key) - 1] + value["mut"] + seq.seq[int(key) :]
    # SeqIO.write(seq, "root.fasta", "fasta")
    return seq


def _compare_sequences(ref, root):
    """Compare the root sequence to the mutated sequence."""
    additional_muts = {}
    for i, nuc in enumerate(
        zip(
            ref.seq,
            root,
        )
    ):
        ref_nuc = nuc[0]
        root_nuc = nuc[1]
        if (
            ref_nuc.upper() != root_nuc.upper()
            and ref_nuc.upper() != "N"
            and root_nuc.upper() != "N"
        ):
            additional_muts[i + 1] = {"ref": ref_nuc.upper(), "root": root_nuc.upper()}
    return additional_muts


# generate a consensus root sequence
def _derive_root_sequence(root_seqs):
    """Generate a consensus root sequence."""
    # for each position in the sequence get the most common nucleotide and add it to the consensus sequence
    consensus_root = ""
    for i in range(len(root_seqs[0].seq)):
        nucs = []
        for seq in root_seqs:
            nucs.append(seq.seq[i])
        # eliminate nucleotides that are not A, T, C, or G or if all nucleotides are the same
        nucs = [
            x
            for x in nucs
            if x.upper() in ["A", "T", "C", "G", "N"] or len(set(nucs)) == 1
        ]
        consensus_root += max(set(nucs), key=nucs.count)
    return consensus_root


def _parse_tree_paths(df):
    df = df.set_index("clade")
    df["from_tree_root"] = df["from_tree_root"].apply(lambda x: x.split(" "))
    return df


def _sanitize_mutation_data(mutations):
    tmp_mutations = {}
    for i in mutations.keys():
        # Note: Insertions and deletions are not included in the additional mutations list
        if "-" in mutations[i]["root"] or "-" in mutations[i]["ref"]:
            continue
        # NOTE: Reversions are not included in the additional mutations list
        if mutations[i]["ref"] == mutations[i]["root"]:
            continue
        tmp_mutations[i] = mutations[i]
    return tmp_mutations


def process_and_reroot_lineages(
    debug: bool,
    sample_muts_path: str,
    reference_fasta_path: str,
    sequences_fasta_path: str,
    input_lineage_paths_path: str,
    output_additional_muts_path: str,
    output_rerooted_lineage_paths_path: str,
):
    """
    Processes mutation data, identifies additional mutations relative to a reference or
    a generated root, and updates lineage path files.
    """
    sample_muts_df = _load_sample_mutations(sample_muts_path)
    seqs = SeqIO.to_dict(SeqIO.parse(sequences_fasta_path, "fasta"))
    ref = SeqIO.read(reference_fasta_path, "fasta")

    # if reference in the sample mutations file, use that as the root
    if sample_muts_df[sample_muts_df["sample"] == ref.id].shape[0] > 0:
        console.print(
            f"[{STYLES['success']}]Reference {ref.id} is present in sample mutations file.[/{STYLES['success']}]"
        )
        additional_muts = _reverse_mutations_to_root(
            _extract_mutations(
                sample_muts_df[sample_muts_df["sample"] == ref.id].iloc[0]
            )
        )
        # change base key to ref and mut key to root
        for i in additional_muts.keys():
            additional_muts[i]["ref"] = additional_muts[i].pop("base")
            additional_muts[i]["root"] = additional_muts[i].pop("mut")

        if debug:
            console.print(
                f"[{STYLES['debug']}]Additional mutations derived from reference {ref.id}: {additional_muts}[/{STYLES['debug']}]"
            )
    # else generate the root sequence
    else:
        console.print(
            f"[{STYLES['warning']}]Reference {ref.id} not present in sample mutations file. Inferring root sequence.[/{STYLES['warning']}]"
        )
        # Pre‑filter samples with non‑null mutations
        valid = sample_muts_df.loc[
            sample_muts_df["mutations"].notnull(), ["sample", "mutations"]
        ]
        root_seqs = []

        for sample_id, muts_str in valid.itertuples(index=False, name=None):
            # build root mutations and fetch the sequence by direct dict lookup
            root_muts = _reverse_mutations_to_root(
                _extract_mutations({"mutations": muts_str})
            )
            seq = seqs.get(sample_id, None)
            if seq is None:
                # It's better to raise an error or handle this case explicitly
                console.print(
                    f"[{STYLES['warning']}]Warning: Sample {sample_id} not found in FASTA file. Skipping.[/{STYLES['warning']}]"
                )
                break
            root_seqs.append(_construct_root_sequence(root_muts, seq))

        if not root_seqs:
            raise ValueError(
                "No valid root sequences could be generated. Check input FASTA and sample mutations."
            )

        root = _derive_root_sequence(root_seqs)
        additional_muts = _compare_sequences(ref, root)

    additional_muts = _sanitize_mutation_data(additional_muts)
    # convert to dataframe and save as csv
    df = pd.DataFrame.from_dict(additional_muts, orient="index")
    df.to_csv(output_additional_muts_path, sep="\t", index_label="position")
    console.print(
        f"[{STYLES['success']}]Additional mutations saved to {output_additional_muts_path}[/{STYLES['success']}]"
    )

    lineage_paths_df = _parse_tree_paths(
        pd.read_csv(input_lineage_paths_path, sep="\t").fillna("")
    )
    # append additional mutations to the begining from_tree_root column list
    additional_muts_list = []
    for i in additional_muts.keys():
        additional_muts_list.append(
            str(additional_muts[i]["ref"] + str(i) + additional_muts[i]["root"])
        )

    if not additional_muts_list:
        console.print(
            f"[{STYLES['warning']}]No additional mutations found to add to lineage paths.[/{STYLES['warning']}]"
        )
        # If no additional mutations, write the original lineage paths content or handle as needed
        # For now, let's just save the parsed (and potentially slightly reformatted) df
        lineage_paths_df["from_tree_root"] = lineage_paths_df["from_tree_root"].apply(
            lambda x: " ".join(x)
        )
    else:
        console.print(
            f"[{STYLES['info']}]Found {len(additional_muts_list)} additional mutations to incorporate into lineage paths.[/{STYLES['info']}]"
        )

        # add the additional mutations to the lineage paths after the first item
        # Ensure the logic for constructing the path string is robust
        def update_path(path_list):
            if not path_list:  # Handle empty path_list
                return ",".join(additional_muts_list)
            first_node = path_list[0]
            remaining_nodes = " ".join(path_list[1:])
            # Construct the new path string carefully
            new_path_parts = [first_node]
            if additional_muts_list:  # Add mutations only if they exist
                new_path_parts.append(",".join(additional_muts_list))
            if remaining_nodes:  # Add remaining nodes only if they exist
                new_path_parts.append(remaining_nodes)
            return " > ".join(new_path_parts)

        lineage_paths_df["from_tree_root"] = lineage_paths_df["from_tree_root"].apply(
            update_path
        )

    lineage_paths_df.to_csv(output_rerooted_lineage_paths_path, sep="\t")
    console.print(
        f"[{STYLES['success']}]Rerooted lineage paths saved to {output_rerooted_lineage_paths_path}[/{STYLES['success']}]"
    )
