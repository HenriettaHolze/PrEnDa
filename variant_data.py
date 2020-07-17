from copy import deepcopy
import pandas as pd
from sys import exit

# new exception class for checking input sequence
class SequenceInconsistency(Exception):
    def __init__(self, message, sequence):
        super().__init__(message)
        self.sequence = sequence


def apply_mutations(sequence, mutations, parent=True):
    """Get FASTA-mutationindices and apply to a sequence"""

    if type(mutations) != list:
        raise Exception("please input mutations as list")

    site_list = []
    sequence_list = list(sequence)

    # potentially multiple mutations, in Fireprot only 1 per variant
    for mutation in mutations:
        # -1 because of one-indexed FASTA-mutation index
        site = int(mutation[1:-1]) - 1
        # print(site)
        # check if multiple mutations in same site
        if str(site) in site_list:
            print("muliple mutations at site", site + 1)
            raise SequenceInconsistency(
                "muliple mutations at site {}".format(site + 1), sequence
            )
            # exit(1)
        else:
            site_list.append(str(site))

        # if given sequence is supposed to be parent, check if start aa fits
        if parent:
            # print('old sequence:', sequence[site].upper())
            if sequence[site].upper() != mutation[0].upper():
                print("not parent sequence at site", site + 1)
                raise SequenceInconsistency(
                    "not parent sequence at site {}".format(site + 1), sequence
                )
                # exit(1)

        # change sequence
        # sequence_list = list(sequence)
        sequence_list[site] = mutation[-1]
    sequence = "".join(sequence_list)

    return sequence


def graph_to_df(df, variantGraph):
    """
    Add information from a variant graph structure to a variant df

    parentAll = set of all parents of a variant (set of string)
    mutationAll = set of mutations of variant and all it's parents (set of string)
    sequence = wt sequence with all mutations from mutationAll (string)
    numMutations = number of mutationAll (integer)
    """

    df_extended = deepcopy(df)
    # initiate column and convert to object, so cells can contain sets
    df_extended["mutationAll"] = None
    df_extended["mutationAll"] = df_extended["mutationAll"].astype("object")
    df_extended["sequence"] = ""
    df_extended["numMutation"] = 0
    df_extended["parentAll"] = None
    df_extended["parentAll"] = df_extended["parentAll"].astype("object")

    for variant in variantGraph:
        # in the graph, parents are represented as variant objects -> convert to set of strings
        parentAllSet = set()
        for parent in variant.parentAll:
            parentAllSet.add(parent.variantId)
        # get index of variant in df
        i = df_extended[df_extended["variant"] == variant.variantId].index.tolist()[0]
        # add mutationAll, sequence, numMutations and parentAll to row
        # this line raises warning: A value is trying to be set on a copy of a slice from a DataFrame
        df_extended["mutationAll"].iloc[i] = list(variant.mutationAll)
        df_extended["sequence"].iloc[i] = variant.sequence
        df_extended["numMutation"].iloc[i] = variant.numMutation
        df_extended["parentAll"].loc[i] = list(parentAllSet)

    return df_extended


def df_to_graph(df, wtSequence):
    """
    Converts a df with variant data into a directed graph. 
    Each node is a variant object and points to it's parent(s).
    There can be multiple initial nodes, because as of now not the wt is the initial node but the first variants.
    """
    # initiate variant graph
    # could also be class
    variantGraph = set()

    # structure for a single variant
    class Variant:
        def __init__(self, variantId):
            self.variantId = variantId
            self.parent = set()
            self.parentAll = set()
            self.mutationNew = set()
            self.mutationAll = set()
            self.numMutation = 0
            self.sequence = ""

    # add variants to graph as Variant objects
    for i, row in df.iterrows():
        variant = Variant(row["variant"])
        # add own mutations to object
        # if not pd.isnull(row['mutation']):
        if row["mutation"] != None:
            variant.mutationNew.update(row["mutation"].split(" "))
            variant.mutationAll.update(row["mutation"].split(" "))
        variantGraph.add(variant)

    # add pointer to parents
    for variant in variantGraph:
        # check if variant has a parent
        # if not pd.isnull(df[df['variant'] == variant.variantId]['parent'].tolist()[0]):
        if df[df["variant"] == variant.variantId]["parent"].values[0] != None:
            parents = (
                df[df["variant"] == variant.variantId]["parent"].values[0].split(" ")
            )
            # iterate over parents
            for parentId in parents:
                # find the parent in the graph and add the link
                for potentialParent in variantGraph:
                    if potentialParent.variantId == parentId:
                        variant.parent.add(potentialParent)
                        variant.parentAll.add(potentialParent)

    # add all mutations by backtracking
    # basically a search algorithm
    for variant in variantGraph:
        parentsToExplore = []
        parentsExplored = []
        for parent in variant.parent:
            parentsToExplore.append(parent)
        # as long as there are parents to explore
        while parentsToExplore != []:
            # remove one item
            cursor = parentsToExplore.pop()
            # add item to explored list
            parentsExplored.append(cursor)
            # add mutations of item to variants mutationAll ant item to list of all parents
            variant.mutationAll.update(cursor.mutationNew)
            variant.parentAll.add(cursor)
            # add parents of cursor to list of parents to explore if they haven't already been explored
            for parent in cursor.parent:
                if parent not in parentsExplored:
                    parentsToExplore.append(parent)

    # apply mutations to wt sequence and add number of mutations
    for variant in variantGraph:
        # print(variant.variantId)
        variant.sequence = apply_mutations(wtSequence, list(variant.mutationAll))
        variant.numMutation = len(variant.mutationAll)

    return variantGraph
