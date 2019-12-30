
# FOR REFERENCE
RADICALS = ["meth", "eth", "prop", "but", "pent", "hex", "hept", "oct", "non", "dec", "undec", "dodec",
            "tridec", "tetradec", "pentadec", "hexadec", "heptadec", "octadec", "nonadec"]
MULTIPLIERS = ["", "di", "tri", "tetra", "penta", "hexa", "hepta", "octa", "nona", "deca", "undeca", "dodeca",
               "trideca", "tetradeca", "pentadeca", "hexadeca", "heptadeca", "octadeca", "nonadeca"]
ROOTS = [
    "an", "en", "yn"
]
REDUNDANT_SUFFIXES = {
    "ol": "hydroxy",
    "one": "oxo",
    "amide": "amido",
    "amine": "amino",
    "imine": "imino",
    "benzene": "phenyl",
    "thiol": "mercapto",
    "phosphine": "phosphino",
    "arsine": "arsino",
    "carboxylicacid": "carboxy"
}
SUFFIXES = {
    "al": {"O": 1, "H": -2},
    "oicacid": {"O": 2, "H": -2},
    "ether": {"O": 1},
}
PREFIXES = {
    "cyclo": {"H": -2},
    "hydroxy": {"O": 1},
    "oxo": {"O": 1, "H": -2},
    "carboxy": {"O": 2, "C": 1},
    "formyl": {"C": 1, "O": 1},
    "amido": {"N": 1, "O": 1, "H": -1},
    "amino": {"N": 1, "H": 1},
    "imino": {"N": 1, "H": -1},
    "phenyl": {"C": 6, "H": 4},
    "mercapto": {"S": 1},
    "phosphino": {"P": 1, "H": 1},
    "arsino": {"As": 1, "H": 1},
    "fluoro": {"F": 1, "H": -1},
    "chloro": {"Cl": 1, "H": -1},
    "bromo": {"Br": 1, "H": -1},
    "iodo": {"I": 1, "H": -1}
}
MEDIXES = {
    "oxy": {"O": 1},
    "yl": {},
    "oyl": {"O": 1, "H": -2},
    "carbonyl": {"C": 1, "O": 1},
    "oate": {"O": 2, "H": -2},
}


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
    # Return the parsed parts
    return parts


# Print the molecule tree.
def recursively_print(node, indent=0):
    print("  " * indent, node.title)
    for substituent in node.substituents:
        print("  " * indent, substituent[0])
        recursively_print(substituent[1], indent + 1)


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

    class ASTNode:
        def __init__(self, title="LOREM_IPSUM"):
            self.title = title
            self.substituents = []

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


# Convert the Abstract Syntax Tree into an atomic formula by descending through the AST from its root.
def atomize(ast):
    # Multiply substituent atom counts by the number of times the substituent appears.
    def atomic_multiplier(mult, atoms):
        return {atom: atoms[atom] * mult for atom in atoms}

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
    def recursively_atomize(node, indent=0):
        multed_subs = [atomic_multiplier(len(substituent[0]), recursively_atomize(substituent[1], indent + 1))
                       for substituent in node.substituents]
        if node.title in RADICALS:
            carbon_count = RADICALS.index(node.title) + 1
            this_node_value = {"C": carbon_count, "H": carbon_count * 2}
        elif node.title in SUFFIXES:
            this_node_value = SUFFIXES[node.title]
        elif node.title in PREFIXES:
            this_node_value = PREFIXES[node.title]
        elif node.title in MEDIXES:
            this_node_value = MEDIXES[node.title]
        elif node.title in ROOTS:
            desaturation_count = ROOTS.index(node.title)
            this_node_value = {"H": -2 * desaturation_count}
        else:
            this_node_value = {}
        total_atoms = atomic_summer(multed_subs + [this_node_value])
        return total_atoms

    # Call the recursive function and subtract two Hydrogens from its result to match the root.
    return atomic_summer([recursively_atomize(ast)] + [{"H": 2}])


import ctypes

PERIODIC_TABLE = {
    "H": {"valence": 1, "weight": 1},
    "B": {"valence": 3, "weight": 10.8},
    "C": {"valence": 4, "weight": 12},
    "N": {"valence": 3, "weight": 14},
    "As": {"valence": 3, "weight": 74.9},
    "O": {"valence": 2, "weight": 16},
    "F": {"valence": 1, "weight": 19},
    "Mg": {"valence": 2, "weight": 24.3},
    "P": {"valence": 3, "weight": 31},
    "S": {"valence": 2, "weight": 32.1},
    "Cl": {"valence": 1, "weight": 35.5},
    "Br": {"valence": 1, "weight": 80},
    "I": {"valence": 1, "weight": 126.9}
}


def obj(n):
    return ctypes.cast(n, ctypes.py_object).value


def order_elements(elements):
    elts_order = []
    for elt in ["C", "H", "O"]:
        if elt in elements:
            elts_order.append(elt)
            elements.remove(elt)
    elts_order.extend(sorted(elements))
    return elts_order


def stringify_atoms_dict(formula_parts):
    elts_order = order_elements(set(formula_parts.keys()))
    formula = "".join(
        [elt + str(formula_parts[elt]) if formula_parts[elt] != 1 else elt for elt in elts_order])
    return formula


class InvalidBond(Exception):
    pass


class LockedMolecule(Exception):
    pass


class UnlockedMolecule(Exception):
    pass


class EmptyMolecule(Exception):
    pass


class Atom:

    def __init__(self, elt, idn):
        self.element = elt
        self.id = idn
        self.connections = {}

    def __hash__(self):
        return self.id

    def __eq__(self, other):
        return self.id == other.id

    def __str__(self):
        elts_set = set([obj(n).element for n in self.connections]).difference({"H"})
        elts_order = order_elements(elts_set)
        conns_formatted = []
        for elt in elts_order:
            matches = sorted([n for n in self.connections if obj(n).element == elt], key=lambda x: obj(x).id)
            for match in matches:
                for i in range(self.connections[match]):
                    conns_formatted.append(obj(match).element + str(obj(match).id))
        conns_formatted.extend(["H" for n in self.connections if obj(n).element == "H"])
        conns = ",".join(conns_formatted)
        return "Atom({}.{}{}{})".format(self.element, self.id, ": " * (len(conns) > 0), conns)


def connect(a1, a2):
    if a1 == a2:
        raise InvalidBond
    for u, t in [(a1, a2), (a2, a1)]:
        if sum([u.connections[n] for n in u.connections]) == PERIODIC_TABLE[u.element]["valence"]:
            raise InvalidBond
    for u, t in [(a1, a2), (a2, a1)]:
        if id(t) not in u.connections:
            u.connections[id(t)] = 0
        u.connections[id(t)] += 1


class Molecule:
    def __init__(self, name=""):
        print("Initiating!")
        self.name = name
        self.atoms = []
        self.branches = []
        self.closed = False
        self._molecular_weight = 0
        self._formula = {}

    @property
    def formula(self):
        if not self.closed:
            raise UnlockedMolecule
        return self._formula

    @property
    def molecular_weight(self):
        if not self.closed:
            raise UnlockedMolecule
        return self._molecular_weight

    def create(self, element):
        if self.closed:
            raise LockedMolecule
        self.atoms.append(Atom(element, len(self.atoms) + 1))

    def brancher(self, *new_branches):
        self.log("Brancher", new_branches)
        if self.closed:
            raise LockedMolecule
        for new_branch in new_branches:
            self.create("C")
            self.branches.append([id(self.atoms[-1])])
            for i in range(new_branch - 1):
                self.create("C")
                connect(self.atoms[-1], self.atoms[-2])
                self.branches[-1].append(id(self.atoms[-1]))
        return self

    def bounder(self, *bounds):
        self.log("Bounder", bounds)
        if self.closed:
            raise LockedMolecule
        for bound in bounds:
            c1, b1, c2, b2 = bound
            connect(obj(self.branches[b1 - 1][c1 - 1]), obj(self.branches[b2 - 1][c2 - 1]))
        return self

    def mutate(self, *mutations):
        self.log("Mutating", mutations)
        if self.closed:
            raise LockedMolecule
        for mutation in mutations:
            nc, nb, elt = mutation
            target = obj(self.branches[nb - 1][nc - 1])
            if PERIODIC_TABLE[elt]["valence"] < sum([target.connections[neighbor] for neighbor in target.connections]):
                raise InvalidBond
            target.element = elt
        return self

    def add(self, *additions):
        self.log("Adding", additions)
        if self.closed:
            raise LockedMolecule
        for addition in additions:
            nc, nb, elt = addition
            target = obj(self.branches[nb - 1][nc - 1])
            if PERIODIC_TABLE[target.element]["valence"] <= sum(
                    [target.connections[neighbor] for neighbor in target.connections]):
                raise InvalidBond
            self.create(elt)
            connect(target, self.atoms[-1])
        return self

    def add_chaining(self, *args):
        self.log("Adding chains", args)
        if self.closed:
            raise LockedMolecule
        nc, nb, newbs = args[0], args[1], args[2:]
        for newb in newbs[:-1]:
            if PERIODIC_TABLE[newb]["valence"] == 1:
                raise InvalidBond
        root = obj(self.branches[nb - 1][nc - 1])
        if PERIODIC_TABLE[root.element]["valence"] <= sum(
                [root.connections[neighbor] for neighbor in root.connections]):
            raise InvalidBond
        for i, elt in enumerate(newbs):
            self.create(elt)
            connect(root, self.atoms[-1])
            root = self.atoms[-1]
        return self

    def closer(self):
        self.log("Closing")
        if self.closed:
            raise LockedMolecule
        formula_parts = {}
        for atom in self.atoms.copy():
            for i in range(
                    PERIODIC_TABLE[atom.element]["valence"] - sum([atom.connections[a] for a in atom.connections])):
                self.create("H")
                connect(atom, self.atoms[-1])
        for atom in self.atoms:
            if atom.element not in formula_parts:
                formula_parts[atom.element] = 0
            formula_parts[atom.element] += 1
        self._formula = stringify_atoms_dict(formula_parts)
        self.formula_dict = formula_parts
        self._molecular_weight = sum(formula_parts[elt] * PERIODIC_TABLE[elt]["weight"] for elt in formula_parts)
        self.closed = True
        return self

    def unlock(self):
        self.log("Opening")
        if not self.closed:
            raise UnlockedMolecule
        # Sift out empty branches
        i = len(self.branches) - 1
        while i >= 0:
            j = len(self.branches[i]) - 1
            while j >= 0:
                if obj(self.branches[i][j]).element == "H":
                    self.branches[i].pop(j)
                j -= 1
            if len(self.branches[i]) == 0:
                self.branches.pop(i)
            i -= 1
        if len(self.branches) == 0:
            raise EmptyMolecule
        # Sift out Hydrogens
        i = len(self.atoms) - 1
        while i >= 0:
            if self.atoms[i].element == "H":
                for n in self.atoms[i].connections:
                    obj(n).connections = {t: obj(n).connections[t] for t in obj(n).connections if obj(t).element != "H"}
                self.atoms.pop(i)
            i -= 1
        # Renumber
        for i, atom in enumerate(self.atoms):
            atom.id = i + 1
        # Set open
        self.closed = False
        return self

    def recursively_populate(self, node, root_branch_id=0, root_carbon_number=0):
        # Add the current node's title to the molecule
        this_branch_id = len(self.branches) + 1
        if node.title in RADICALS:
            self.brancher(RADICALS.index(node.title) + 1)
        elif node.title == "mercapto":
            self.brancher(1)
            self.mutate((1, this_branch_id, "S"))
        elif node.title in {"oxo", "al", "oyl", "oxy", "hydroxy", "ether"}:
            self.brancher(1)
            self.mutate((1, this_branch_id, "O"))
        elif node.title == "fluoro":
            self.brancher(1)
            self.mutate((1, this_branch_id, "F"))
        elif node.title == "chloro":
            self.brancher(1)
            self.mutate((1, this_branch_id, "Cl"))
        elif node.title == "bromo":
            self.brancher(1)
            self.mutate((1, this_branch_id, "Br"))
        elif node.title == "iodo":
            self.brancher(1)
            self.mutate((1, this_branch_id, "I"))
        elif node.title == "carboxy":
            self.brancher(2)
            self.mutate((2, this_branch_id, "O"))
            self.bounder((1, this_branch_id, 2, this_branch_id))
            self.add((1, this_branch_id, "O"))
        elif node.title == "formyl":
            self.brancher(2)
            self.mutate((2, this_branch_id, "O"))
            self.bounder((1, this_branch_id, 2, this_branch_id))
        elif node.title == "amido":
            self.brancher(1)
            self.mutate((1, len(self.branches), "N"))
            self.brancher(1)
            self.mutate((1, len(self.branches), "O"))
            self.bounder((root_carbon_number, root_branch_id, 1, this_branch_id + 1))
            self.bounder((root_carbon_number, root_branch_id, 1, this_branch_id + 1))
        elif node.title == "amino" or node.title == "imino":
            self.brancher(1)
            self.mutate((1, this_branch_id, "N"))
        elif node.title == "phenyl":
            self.brancher(6)
            self.bounder((1, this_branch_id, 6, this_branch_id))
            self.bounder((1, this_branch_id, 2, this_branch_id))
            self.bounder((3, this_branch_id, 4, this_branch_id))
            self.bounder((5, this_branch_id, 6, this_branch_id))
        elif node.title == "phosphino":
            self.brancher(1)
            self.mutate((1, this_branch_id, "P"))
        elif node.title == "arsino":
            self.brancher(1)
            self.mutate((1, this_branch_id, "As"))
        elif node.title == "carbonyl":
            self.brancher(2)
            self.mutate((2, this_branch_id, "O"))
            self.bounder((1, this_branch_id, 2, this_branch_id))
        elif node.title == "oate" or node.title == "oicacid":
            self.brancher(1)
            self.mutate((1, len(self.branches), "O"))
            self.brancher(1)
            self.mutate((1, len(self.branches), "O"))
            self.bounder((root_carbon_number, root_branch_id, 1, this_branch_id + 1))
            self.bounder((root_carbon_number, root_branch_id, 1, this_branch_id + 1))
        # Link the new branch to the root at the correct location
        if root_branch_id and root_carbon_number:
            self.bounder((root_carbon_number, root_branch_id, 1, this_branch_id))
            if node.title in {"oxo", "al", "imino", "oyl"}:
                self.bounder((root_carbon_number, root_branch_id, 1, this_branch_id))
        # Loop over substituents
        for locations, substituent in node.substituents:
            for location in locations:
                if substituent.title == "cyclo":
                    self.bounder((1, this_branch_id, len(self.branches[this_branch_id - 1]), this_branch_id))
                elif substituent.title in ROOTS:
                    for i in range(ROOTS.index(substituent.title)):
                        self.bounder((location, this_branch_id, location + 1, this_branch_id))
                else:
                    self.recursively_populate(substituent, this_branch_id, location)
        return self

    def log(self, *args):
        #         print("After last command:", [str(atom) for atom in self.atoms])
        #         name, content = args[0], args[1:]
        #         print(name, content)
        #         if name == "Adding":
        #             for trio in content[0]:
        #                 nc, nb, elt = trio
        #                 target = obj(self.branches[nb-1][nc-1])
        #                 print(target.element, target.id)
        pass


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