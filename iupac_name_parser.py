import re
from molecule_builder import *
from config import *


# Tokenize an organic molecule's name so its parts can be easily distinguished.
def tokenize(name):
    print("Name:", name)
    # Pre-process to fix up -oates
    if name[-4:] == "oate":
        mult = 2 if name[-6:-4] == "di" else 1
        name = name[:-6] if mult == 2 else name[:-4]
        sub, root = name.split(" ")
        if mult == 1:
            name = "1-" + sub + "oate" + root + "e"
        elif mult == 2:
            pointer = len(root) - 1
            while root[pointer].isdigit() or root[pointer] in {"(", ",", ")", "-"}:
                pointer -= 1
            name = root[pointer + 2:len(root) - 1] + "-di" + sub + "oate" + root[:pointer + 1] + "e"
    # Grab the locations of the substituents
    formula = name.replace("-", "").replace(" ", "")

    # Generate name parts
    parts = []
    base = 0
    while base < len(formula):
        # If the base char is a digit, we have hit some substituent locations. Grab them.
        if formula[base].isdigit():
            pointer = base + 1
            while formula[pointer].isdigit() or formula[pointer] == ",":
                pointer += 1
            parts.append(tuple(int(c) for c in formula[base:pointer].split(",")))
            pointer += len(MULTIPLIERS[len(parts[-1]) - 1]) if len(parts[-1]) > 1 else 0
        # If the base char is a bracket, add it as a token.
        elif formula[base] in {"[", "]"}:
            parts.append(formula[base])
            pointer = base + 1
        # Otherwise, move forward until a token is found.
        else:
            pointer = base + 2
            while not (pointer >= len(formula)
                       or formula[base:pointer] in MULTIPLIERS or formula[base:pointer] in PREFIXES
                       or formula[base:pointer] in ROOTS or formula[base:pointer] in SUFFIXES
                       or formula[base:pointer] in MEDIXES or formula[base:pointer] in RADICALS):
                pointer += 1
            # Check for the -a difference between an alkyl and a multiplier
            if formula[base:pointer] in RADICALS:
                if formula[pointer] == "a" and \
                        (formula[pointer + 1] not in {"n", "l"} or
                         (formula[pointer + 1:pointer + 4] == "non" and formula[pointer + 4] != "e")):
                    pointer += 1
            # Check for 14+ multipliers & alkyls
            if formula[base:pointer] != "di" and formula[base:pointer] in MULTIPLIERS \
                    and formula[pointer:pointer + 3] == "dec":
                pointer += 3
            # Check for carboxy/carboxylicacid and eth/ether ambiguity
            elif formula[base:] in ["carboxylicacid", "ether"]:
                pointer = len(formula)
            # Append to list of parts
            parts.append(formula[base:pointer])
        base = pointer

    # # Regex out some valid parts
    # options = set(RADICALS) | set(MULTIPLIERS[1:]) | set(ROOTS) | \
    #           set(REDUNDANT_SUFFIXES.keys()) | set(SUFFIXES.keys()) | \
    #           set(PREFIXES.keys()) | set(MEDIXES.keys()) | {"e"}
    # print("Options:", options)
    # number_reg = "(?:[0-9]+(?:,[0-9]+)*)"
    # char_reg = r"[\[\]-]"
    # reg = r"\A(({}|{}|{})+)\Z".format("|".join(options), number_reg, char_reg)
    # # Regex out a valid set of parts
    # from regex import regex
    # regged_parts = []
    # i = 0
    # reg = r"\A(({}|{}|{})(?1)*)\Z".format("|".join(options), number_reg, char_reg)
    # while i < len(formula):
    #     next_res = regex.match(reg, formula[i:]).group(2)
    #     regged_parts.append(next_res)
    #     i += len(next_res)
    # # Kill useless tokens and intify substituent numberings
    # regged_parts = [tuple([int(d) for d in part.split(',')]) if part[0].isdigit() else part
    #                 for part in regged_parts]
    # # Fix tridec
    #
    #
    # # if not (
    # #         i > 0
    # #         and regged_parts[i - 1][0].isdigit()
    # #         and MULTIPLIERS[len(regged_parts[i - 1].split(",")) - 1] == part
    # # )
    #
    # print("PARTS:", parts)
    # print("NEWS:", regged_parts)
    # print("WORKS:", parts == regged_parts)

    # Replace all remaining multipliers with 1s
    if parts[-1] in {"al", "oicacid"} and parts[-2] == "di":
        k = -1
        while parts[k] not in RADICALS:
            k -= 1
        parts[-2] = (1, RADICALS.index(parts[k]) + 1)
    for i in range(len(parts)):
        if isinstance(parts[i], str) and parts[i] in MULTIPLIERS:
            parts[i] = tuple([1] * (MULTIPLIERS.index(parts[i]) + 1))
    # Check whether a root alkyl exists before the SUFFIX.
    root_alkyl_exists = True
    if len(parts) == 1:
        root_alkyl_exists = False
    elif parts[-1] in SUFFIXES or parts[-1] in REDUNDANT_SUFFIXES:
        i = len(parts) - 1
        while i >= 0:
            if parts[i] in RADICALS or parts[i] in PREFIXES:
                if parts[i] in PREFIXES or len([part for part in parts[i:] if part[-2:] == "yl"]) > 0:
                    root_alkyl_exists = False
                break
            i -= 1
    # We no longer need the "yl" and "an" to disprove end roots, so we can kill them off.
    parts = [part for part in parts if part not in {"an", "yl", "e"}]
    # Post-process to add tuple notation to those items that lack it and guarantee a root.
    new_parts = []
    i = len(parts) - 1
    rooted = False
    # Seed the parts with a root SUFFIX it it exists.
    if parts[i] in SUFFIXES or parts[i] in REDUNDANT_SUFFIXES:
        new_parts.append(parts[i])
        if not isinstance(parts[i - 1], tuple):
            if root_alkyl_exists:
                new_parts.append((1,))
                rooted = False
            else:
                new_parts.append((0,))
                rooted = True
        i -= 1
    # Iterate over substituents until complete.
    while i >= 0:
        if (parts[i] in RADICALS or parts[i] in PREFIXES or parts[i] in ROOTS) \
                and (i == 0 or (type(parts[i - 1]) != tuple)):
            # If it is a prefix, it would have a number if it were not on C1
            if parts[i] in PREFIXES:
                new_parts.append(parts[i])
                if i == 0 or parts[i - 1] != "]":
                    new_parts.append((1,))
            elif parts[i] in MEDIXES:
                new_parts.append(parts[i])
            elif parts[i] in RADICALS:
                new_parts.append(parts[i])
                # If there is an end bracket before it, its numbering occurs before the bracket
                if i != 0 and parts[i - 1] == "]":
                    pass
                # Account for triethylamine
                elif rooted:
                    # Account for cyclo prefix
                    if i != 0 and parts[i - 1] == "cyclo":
                        new_parts.append("cyclo")
                        i -= 1
                    # If there is a multiplier ahead, account for it
                    if i != 0 and parts[i - 1] in MULTIPLIERS:
                        new_parts.append(tuple([1 for i in range(MULTIPLIERS.index(parts[i - 1]) + 1)]))
                        i -= 1
                    # If we moved forward, we might now be at a bracket.
                    elif i != 0 and parts[i - 1] == "]" or isinstance(parts[i - 1], tuple):
                        pass
                    # Otherwise, assume we are at a root-level unlabeled substituent
                    else:
                        new_parts.append((1,))
                # If an alkyl substituent has not yet been seen, the first one must be the root.
                else:
                    # Account for cyclo prefix.
                    if i != 0 and parts[i - 1] == "cyclo":
                        new_parts.append("cyclo")
                        new_parts.append((0,))
                        rooted = True
                        i -= 1
                    else:
                        new_parts.append((0,))
                        rooted = True
            elif parts[i] in ROOTS and parts[i] != "an":
                new_parts.append(parts[i])
                new_parts.append((1,))
        else:
            new_parts.append(parts[i])
        i -= 1
    # Add a preceeding 1 if none exists
    if not isinstance(new_parts[-1], tuple) or parts[-1] == "[":
        new_parts.append((1,))
    # Reverse the constructed list
    parts = list(reversed(new_parts))
    # Kill erroneous end tuples left over from popping the "an".
    if isinstance(parts[-1], tuple):
        parts.pop()
    # Insert 1 between contiguous [[.
    i = 1
    new_parts = [parts[0]]
    while i < len(parts):
        if parts[i] == "[":
            if not isinstance(parts[i - 1], tuple):
                new_parts.append((1,))
            new_parts.append(parts[i])
            if i < len(parts) - 1 and not (isinstance(parts[i + 1], tuple) or parts[i + 1] == "["):
                new_parts.append((1,))
        else:
            new_parts.append(parts[i])
        i += 1
    parts = new_parts
    # Replace the SUFFIX with a PREFIX if possible to reduce reconstruction overhead.
    if parts[-1] in REDUNDANT_SUFFIXES:
        parts[-1] = REDUNDANT_SUFFIXES[parts[-1]]
    # Print final parts:
    print("Final Parts:", parts)
    # Return the parsed parts
    return parts


# Print the molecule tree.
def recursively_print(node, indent=0):
    print("  " * indent, node.title)
    for substituent in node.substituents:
        print("  " * indent, substituent[0])
        recursively_print(substituent[1], indent + 1)


# A class for the nodes of the ast
class ASTNode:
    def __init__(self, title="LOREM_IPSUM"):
        self.title = title
        self.substituents = []

    # Multiply substituent atom counts by the number of times the substituent appears.
    @staticmethod
    def atomic_multiplier(mult, atoms):
        return {atom: atoms[atom] * mult for atom in atoms}

    @staticmethod
    # Sum up the various atom counts of a substituent's substituents.
    def atomic_summer(atoms_list):
        total_dict = {}
        for atoms in atoms_list:
            for atom in atoms:
                if atom not in total_dict:
                    total_dict[atom] = 0
                total_dict[atom] += atoms[atom]
        return total_dict

    # Recursively descend through the AST to find the total atom count.
    def recursively_atomize(self, indent=0):
        multed_subs = [ASTNode.atomic_multiplier(len(substituent[0]), substituent[1].recursively_atomize(indent + 1))
                       for substituent in self.substituents]
        if self.title in RADICALS:
            carbon_count = RADICALS.index(self.title) + 1
            this_node_value = {"C": carbon_count, "H": carbon_count * 2}
        elif self.title in SUFFIXES:
            this_node_value = SUFFIXES[self.title]
        elif self.title in PREFIXES:
            this_node_value = PREFIXES[self.title]
        elif self.title in MEDIXES:
            this_node_value = MEDIXES[self.title]
        elif self.title in ROOTS:
            desaturation_count = ROOTS.index(self.title)
            this_node_value = {"H": -2 * desaturation_count}
        else:
            this_node_value = {}
        total_atoms = ASTNode.atomic_summer(multed_subs + [this_node_value])
        return total_atoms

    # Convert the Abstract Syntax Tree into an atomic formula by descending through the AST from its root.
    def atomize(self):
        # Call the recursive function and subtract two Hydrogens from its result to match the root.
        return ASTNode.atomic_summer([self.recursively_atomize()] + [{"H": 2}])


# Convert the tokenized name into an Abstract Syntax Tree.
def astify(parts):
    # Get parentheses pairs
    p_d = {}
    p_q = []
    for p_i in range(len(parts)):
        if parts[p_i] == "[":
            p_q.append(p_i)
        elif parts[p_i] == "]":
            p_d[p_q.pop()] = p_i

    root = ASTNode()
    q = [root]
    i = 0
    while i < len(parts):
        # Account for MEDIXES
        if parts[i] in MEDIXES:
            numbering, progress_root = q[-1].substituents.pop()
            while parts[i] in MEDIXES:
                if parts[i] == "oyl":
                    progress_root.substituents.append(((1,), ASTNode(parts[i])))
                else:
                    next_root = ASTNode(parts[i])
                    next_root.substituents.append(((1,), progress_root))
                    progress_root = next_root
                i += 1
            q[-1].substituents.append((numbering, progress_root))
        # Look for the next substituent.
        # If the next token is an open bracket, generate a new node.
        if parts[i + 1] == "[":
            q.append(ASTNode())
            q[-2].substituents.append((parts[i], q[-1]))
        # If the next token is a closed bracket, close the latest node.
        elif parts[i] == "]":
            if parts[i + 1] == "cyclo":
                q[-1].title = parts[i + 2]
                q[-1].substituents.append(((1,), ASTNode("cyclo")))
                q.pop(-1)
                i += 1
            else:
                q[-1].title = parts[i + 1]
                q.pop(-1)
        # If the next token is a root, append it to the appropriate substituent.
        elif parts[i + 1] in ROOTS:
            q[-1].substituents[-1][1].substituents.append((parts[i], ASTNode(parts[i + 1])))
        # If the next token is a cyclo, move up to the following alkyl and add the cyclo as a substituent.
        elif parts[i + 1] == "cyclo":
            next_node = ASTNode(parts[i + 2])
            next_node.substituents.append(((1,), ASTNode("cyclo")))
            q[-1].substituents.append((parts[i], next_node))
            i += 1
        # If none of the above special cases apply, simply add the next token as a substituent to the latest node.
        else:
            q[-1].substituents.append((parts[i], ASTNode(parts[i + 1])))
        i += 2
    # Account for root's name
    i = len(root.substituents) - 1
    while True:
        if root.substituents[i][0] == (0,):
            actual_root = root.substituents.pop(i)[1]
            root.title = actual_root.title
            root.substituents.extend(actual_root.substituents)
            break
        i -= 1
    # Return the completed root
    return root


# Declare a class that applies the tokenization & AST logic.
class ParseHer:

    def __init__(self, name):
        # Grab the parts of the name
        tokens = tokenize(name)
        # Convert into an Abstract Syntax Tree
        ast = astify(tokens)
        # Reduce AST to atoms
        self.atoms = Molecule().recursively_populate(ast).closer().formula_dict

    def parse(self):
        # Parse the name given as argument in the constructor and output the dict representing the raw formula
        return self.atoms


tokenize("hexahexahexhexhexahexandiol")
print()
print()
