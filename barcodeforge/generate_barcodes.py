"""Generates barcodes from UShER lineage path files."""

# This script is a modified version of the https://github.com/andersen-lab/Freyja/blob/main/freyja/convert_paths2barcodes.py

import pandas as pd
from rich.console import Console
from .utils import sortFun, STYLES

console = Console()


def parse_tree_paths(df: pd.DataFrame) -> pd.DataFrame:
    """
    Parses lineage paths from a DataFrame and extracts clade information.
    This function processes a DataFrame containing lineage paths, extracting clade names
    and their corresponding mutations from the 'from_tree_root' column. It creates a new DataFrame
    with clade names as the index and a list of mutations for each clade.
    Args:
        df (pd.DataFrame): DataFrame containing lineage paths with columns 'clade' and 'from_tree_root'.
    Returns:
        pd.DataFrame: A DataFrame with clade names as the index and a list of mutations for each clade.
    """
    df = df.set_index("clade")
    # Make sure to check with new tree versions, lineages could get trimmed.
    df = df.drop_duplicates(keep="last")
    df["from_tree_root"] = df["from_tree_root"].fillna("")
    df["from_tree_root"] = df["from_tree_root"].apply(
        lambda x: x.replace(" ", "").strip(">").split(">")
    )
    return df


def convert_to_barcodes(df: pd.DataFrame) -> pd.DataFrame:
    """
    Converts lineage paths DataFrame to barcodes DataFrame.
    This function takes a DataFrame with lineage paths and converts it into a
    barcodes DataFrame, where each column represents a mutation and each row
    represents a clade. The mutations are encoded as binary values, indicating
    whether a mutation is present in a clade or not. The function also handles
    the separation of combined mutations and checks for reversion pairs.
    Args:
        df (pd.DataFrame): DataFrame containing lineage paths with columns 'clade' and 'from_tree_root'.
    Returns:
        pd.DataFrame: A DataFrame where each column represents a mutation and each row represents a clade.
                      The values are binary, indicating the presence of mutations.
    """
    # builds simple barcodes, not accounting for reversions
    df_barcodes = pd.DataFrame()
    for clade in df.index:
        # sparse,binary encoding
        cladeSeries = pd.Series(
            {
                c: df.loc[clade, "from_tree_root"].count(c)
                for c in df.loc[clade, "from_tree_root"]
            },
            name=clade,
        )
        df_barcodes = pd.concat((df_barcodes, cladeSeries), axis=1)

    # console.print('separating combined splits') # Original print, can be made conditional with a verbose flag if needed
    df_barcodes = df_barcodes.T
    # dropped since no '' column this time.
    # df_barcodes = df_barcodes.drop(columns='')
    df_barcodes = df_barcodes.fillna(0)
    temp = pd.DataFrame()
    dropList = []
    for c in df_barcodes.columns:
        # if column includes multiple mutations,
        # split into separate columns and concatenates
        if "," in c:
            for mt in c.split(","):
                if mt not in temp.columns:
                    temp = pd.concat((temp, df_barcodes[c].rename(mt)), axis=1)
                else:
                    # to handle multiple different groups with mut
                    temp[mt] += df_barcodes[c]
            dropList.append(c)
    df_barcodes = df_barcodes.drop(columns=dropList)
    df_barcodes = pd.concat((df_barcodes, temp), axis=1)
    df_barcodes = df_barcodes.T.groupby(level=0).sum().T

    # drop columns with empty strings
    # Warning: this is a hack to deal with empty strings in the O/P from matUtils extract.
    if "" in df_barcodes.columns:
        df_barcodes = df_barcodes.drop(columns="")
    return df_barcodes


def reversion_checking(df_barcodes: pd.DataFrame) -> pd.DataFrame:
    """
    Check for reversion pairs in the barcodes DataFrame.
    This function identifies pairs of mutations that are reversions of each other
    and adjusts the counts accordingly.
    Args:
        df_barcodes (pd.DataFrame): DataFrame containing barcodes with mutations as columns.
    Returns:
        pd.DataFrame: The DataFrame with reversion pairs adjusted.
    """
    flipPairs = [
        (d, d[-1] + d[1 : len(d) - 1] + d[0])
        for d in df_barcodes.columns
        if (d[-1] + d[1 : len(d) - 1] + d[0]) in df_barcodes.columns
    ]
    flipPairs = [list(fp) for fp in list(set(flipPairs))]
    # subtract lower of two pair counts to get the lineage defining mutations
    for fp in flipPairs:
        df_barcodes[fp] = df_barcodes[fp].subtract(df_barcodes[fp].min(axis=1), axis=0)
    # drop all unused mutations (i.e. paired mutations with reversions)
    df_barcodes = df_barcodes.drop(
        columns=df_barcodes.columns[df_barcodes.sum(axis=0) == 0]
    )
    return df_barcodes


def check_no_flip_pairs(barcode_file: str):
    """
    Test if there are any flip pairs in the generated barcode file.
    Args:
        barcode_file (str): Path to the barcode file to be tested.
    Raises:
        Exception: If flip pairs are found in the barcode file.
    """
    df_barcodes = pd.read_csv(barcode_file, index_col=0)
    flipPairs = [
        (d, d[-1] + d[1 : len(d) - 1] + d[0])
        for d in df_barcodes.columns
        if (d[-1] + d[1 : len(d) - 1] + d[0]) in df_barcodes.columns
    ]
    if len(flipPairs) == 0:
        console.print(
            f"[{STYLES['success']}]PASS: no flip pairs found in the generated barcode file.[/{STYLES['success']}]"
        )
    else:
        # This should ideally not happen if the logic is correct, consider logging or raising a more specific error.
        console.print(
            f"[{STYLES['error']}]FAIL: flip pairs found: {flipPairs}[/{STYLES['error']}]"
        )
        raise Exception(f"FAIL: flip pairs found: {flipPairs}")


def identify_chains(df_barcodes: pd.DataFrame) -> list:
    """
    Identify sequential mutations in the barcodes DataFrame.
    This function looks for mutations that can be sequentially combined
    based on the last character of the mutation strings.
    Args:
        df_barcodes (pd.DataFrame): DataFrame containing barcodes with mutations as columns.
    Returns:
        list: A list of lists, where each inner list contains three elements:
              [original mutation, next mutation, combined mutation].
    """

    sites = [d[0 : len(d) - 1] for d in df_barcodes.columns]
    flip_sites = [d[-1] + d[1 : len(d) - 1] for d in df_barcodes.columns]
    # for each mutation, find possible sequential mutations
    seq_muts = [
        [d, df_barcodes.columns[j], d[0 : len(d) - 1] + df_barcodes.columns[j][-1]]
        for i, d in enumerate(df_barcodes.columns)
        for j, d2 in enumerate(sites)
        if (
            (flip_sites[i] == sites[j])
            and (d[-1] + d[1 : len(d) - 1] + d[0]) != df_barcodes.columns[j]
        )
    ]

    # confirm that mutation sequence is actually observed
    seq_muts = [
        sm
        for sm in seq_muts
        if df_barcodes[(df_barcodes[sm[0]] > 0) & (df_barcodes[sm[1]] > 0)].shape[0] > 0
    ]

    mut_sites = [sortFun(sm[2]) for sm in seq_muts]
    # return only one mutation per site for each iteration
    seq_muts = [
        seq_muts[i] for i, ms in enumerate(mut_sites) if ms not in mut_sites[:i]
    ]
    return seq_muts


def check_mutation_chain(df_barcodes: pd.DataFrame) -> pd.DataFrame:
    # case when (non-reversion) mutation happens in site with existing mutation
    seq_muts = identify_chains(df_barcodes)
    while len(seq_muts) > 0:
        # combine mutations string into single mutation
        for i, sm in enumerate(seq_muts):
            lin_seq = df_barcodes[(df_barcodes[sm[0]] > 0) & (df_barcodes[sm[1]] > 0)]
            if sm[2] not in df_barcodes.columns:
                # combination leads to new mutation
                newCol = pd.Series(
                    [(1 if dfi in lin_seq.index else 0) for dfi in df_barcodes.index],
                    name=sm[2],
                    index=df_barcodes.index,
                )
                df_barcodes = pd.concat([df_barcodes, newCol], axis=1)
            else:
                # combining leads to already existing mutation
                # just add in that mutation
                df_barcodes.loc[lin_seq.index, sm[2]] = 1
            # remove constituent mutations
            df_barcodes.loc[lin_seq.index, sm[0:2]] -= 1
        # drop all unused mutations
        # print('before_trim\n',df_barcodes)
        df_barcodes = df_barcodes.drop(
            columns=df_barcodes.columns[df_barcodes.sum(axis=0) == 0]
        )
        # in case mutation path leads to a return to the reference.
        df_barcodes = reversion_checking(df_barcodes)
        seq_muts = identify_chains(df_barcodes)
    return df_barcodes


def replace_underscore_with_dash(df: pd.DataFrame) -> pd.DataFrame:
    """
    Replace underscores with dashes in the DataFrame index labels.
    Parameters:
        df (pandas.DataFrame): A DataFrame whose index labels may contain underscores.
    Returns:
        pandas.DataFrame: The same DataFrame with all underscores in its index labels replaced by dashes.
    """

    df.index = [i.replace("_", "-") for i in df.index]
    return df


def create_barcodes_from_lineage_paths(
    input_file_path: str, output_file_path: str, prefix: str = ""
):
    """
    Generates a barcode CSV file from a lineage paths TSV file.

    Args:
        input_file_path: Path to the input lineage paths file (TSV format from matUtils extract -C).
        output_file_path: Path to save the generated barcodes (CSV format).
        prefix: Optional prefix to add to lineage names in the barcode file.
    """
    console.print(
        f"[{STYLES['info']}]Reading lineage paths from: {input_file_path}[/{STYLES['info']}]"
    )
    df = pd.read_csv(input_file_path, sep="\t")

    console.print(f"[{STYLES['info']}]Parsing tree paths...[/{STYLES['info']}]")
    df = parse_tree_paths(df)

    console.print(f"[{STYLES['info']}]Converting to barcodes...[/{STYLES['info']}]")
    df_barcodes = convert_to_barcodes(df)

    if prefix and prefix.strip() != "":
        console.print(
            f"[{STYLES['info']}]Adding prefix '{prefix}' to lineage names...[/{STYLES['info']}]"
        )
        df_barcodes.index = [prefix + "-" + str(i) for i in df_barcodes.index]

    console.print(
        f"[{STYLES['info']}]Performing reversion checking...[/{STYLES['info']}]"
    )
    df_barcodes = reversion_checking(df_barcodes)

    console.print(f"[{STYLES['info']}]Checking mutation chains...[/{STYLES['info']}]")
    df_barcodes = check_mutation_chain(df_barcodes)

    console.print(
        f"[{STYLES['info']}]Replacing underscores with dashes in lineage names...[/{STYLES['info']}]"
    )
    df_barcodes = replace_underscore_with_dash(df_barcodes)

    console.print(f"[{STYLES['info']}]Sorting barcode columns...[/{STYLES['info']}]")
    df_barcodes = df_barcodes.reindex(sorted(df_barcodes.columns, key=sortFun), axis=1)

    # Drop unclassified lineage if it exists
    if "unclassified" in df_barcodes.index:
        df_barcodes = df_barcodes.drop(index="unclassified")

    df_barcodes.to_csv(output_file_path)
    console.print(
        f"[{STYLES['success']}]Barcode file saved to: {output_file_path}[/{STYLES['success']}]"
    )

    # Test for flip pairs in the final output file
    try:
        check_no_flip_pairs(output_file_path)
    except Exception as e:
        console.print(
            f"[{STYLES['error']}]Error during final flip pair test: {e}[/{STYLES['error']}]"
        )
        # Depending on desired behavior, you might re-raise the exception or just log it.
